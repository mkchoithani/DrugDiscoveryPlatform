"""Chemprop (D-MPNN) behind the same interface as the random forest.

A message-passing neural network that learns directly from the molecular graph,
rather than from a fixed fingerprint. It is the standard strong baseline in the
field and the obvious candidate to replace PropertyModel.

Same fit / predict / apply / report surface, so `model_benchmark.py` can put the
two against identical splits and let the numbers decide. Uncertainty comes from
mean-variance estimation (the network predicts a spread alongside each value),
calibrated against held-out residuals exactly as the forest's is.
"""
from __future__ import annotations

import logging
import os
import warnings

import numpy as np

from .library import Compound
from .models import regression_metrics

# PyTorch and RDKit each ship an OpenMP runtime; loading both segfaults on macOS
# (exit 139, mid-training, no Python traceback). These must be set before torch is
# imported, which is why they sit at module scope and torch is imported lazily below.
os.environ.setdefault("KMP_DUPLICATE_LIB_OK", "TRUE")
os.environ.setdefault("OMP_NUM_THREADS", "1")
os.environ.setdefault("PYTORCH_ENABLE_MPS_FALLBACK", "1")
os.environ.setdefault("TOKENIZERS_PARALLELISM", "false")
logging.getLogger("lightning.pytorch").setLevel(logging.ERROR)
logging.getLogger("pytorch_lightning").setLevel(logging.ERROR)
warnings.filterwarnings("ignore", category=UserWarning)


def chemprop_available() -> bool:
    try:
        import chemprop  # noqa: F401
        import lightning  # noqa: F401
        return True
    except Exception:
        return False


class ChempropModel:
    """D-MPNN regressor with mean-variance uncertainty."""

    def __init__(self, epochs: int = 40, batch_size: int = 32, hidden: int = 300,
                 depth: int = 3, dropout: float = 0.0, seed: int = 0,
                 accelerator: str = "cpu", uncertainty: bool = True):
        self.epochs = epochs
        self.batch_size = batch_size
        self.hidden = hidden
        self.depth = depth
        self.dropout = dropout
        self.seed = seed
        self.accelerator = accelerator
        self.uncertainty = uncertainty
        self.model = None
        self.trainer = None
        self.calibration = 1.0
        self.cv_r2: float | None = None
        self.cv_mae: float | None = None
        self.n_train = 0

    # ---------------------------------------------------------------- build
    def _datapoints(self, compounds: list[Compound], target: str | None):
        from chemprop import data
        values = ([[float(getattr(c, target))] for c in compounds] if target
                  else [None] * len(compounds))
        return [data.MoleculeDatapoint.from_smi(c.smiles, y)
                for c, y in zip(compounds, values)]

    def _build(self, scaler):
        from chemprop import models, nn
        mp = nn.BondMessagePassing(d_h=self.hidden, depth=self.depth, dropout=self.dropout)
        agg = nn.MeanAggregation()
        transform = nn.UnscaleTransform.from_standard_scaler(scaler)
        head = (nn.MveFFN(input_dim=self.hidden, output_transform=transform)
                if self.uncertainty
                else nn.RegressionFFN(input_dim=self.hidden, output_transform=transform))
        return models.MPNN(mp, agg, head, batch_norm=True)

    # ------------------------------------------------------------------ fit
    def fit(self, compounds: list[Compound], target: str = "measured") -> "ChempropModel":
        from chemprop import data, featurizers
        from lightning import pytorch as pl
        import torch

        torch.manual_seed(self.seed)
        torch.set_num_threads(1)
        self.n_train = len(compounds)
        featurizer = featurizers.SimpleMoleculeMolGraphFeaturizer()
        dset = data.MoleculeDataset(self._datapoints(compounds, target), featurizer)
        scaler = dset.normalize_targets()
        loader = data.build_dataloader(dset, batch_size=self.batch_size, shuffle=True,
                                       num_workers=0)

        self.model = self._build(scaler)
        self.trainer = pl.Trainer(
            max_epochs=self.epochs, accelerator=self.accelerator, devices=1,
            logger=False, enable_checkpointing=False, enable_progress_bar=False,
            enable_model_summary=False, deterministic=False,
        )
        self.trainer.fit(self.model, loader)

        # calibrate sigma on the training residuals. Less rigorous than the
        # forest's out-of-fold scheme, and the benchmark reports it as such.
        mu, sigma = self.predict(compounds)
        truth = np.array([float(getattr(c, target)) for c in compounds])
        resid = np.abs(mu - truth)
        mean_sigma = float(np.mean(sigma))
        self.calibration = (float(np.clip(resid.mean() / mean_sigma, 0.2, 5.0))
                            if mean_sigma > 1e-6 else 1.0)
        m = regression_metrics([float(x) for x in mu], list(truth))
        self.cv_r2, self.cv_mae = m["r2"], m["mae"]      # in-sample, not a CV estimate
        return self

    # -------------------------------------------------------------- predict
    def predict(self, compounds: list[Compound]) -> tuple[np.ndarray, np.ndarray]:
        from chemprop import data, featurizers
        import torch

        featurizer = featurizers.SimpleMoleculeMolGraphFeaturizer()
        dset = data.MoleculeDataset(self._datapoints(compounds, None), featurizer)
        # drop_last defaults to dropping a trailing batch of one sample (batch norm
        # cannot normalise a single row). At prediction time that silently returns
        # fewer rows than were asked for, so it is forced off and the batch size is
        # nudged when a remainder of one would occur.
        batch_size = self.batch_size
        if len(compounds) % batch_size == 1 and len(compounds) > 1:
            batch_size += 1
        loader = data.build_dataloader(dset, batch_size=batch_size, shuffle=False,
                                       num_workers=0, drop_last=False)
        self.model.eval()
        outputs = []
        with torch.no_grad():
            for batch in loader:
                outputs.append(self.model(batch.bmg, batch.V_d, batch.X_d).cpu().numpy())
        out = np.concatenate(outputs, axis=0)
        if len(out) != len(compounds):          # never return a silently short answer
            raise RuntimeError(
                f"chemprop returned {len(out)} predictions for {len(compounds)} compounds")

        if out.ndim == 3:                       # (n, tasks, 2) from the MVE head
            mu = out[:, 0, 0]
            var = out[:, 0, 1]
        elif out.shape[-1] == 2:
            mu, var = out[:, 0], out[:, 1]
        else:
            mu = out.reshape(-1)
            var = np.full_like(mu, 1.0)
        sigma = np.sqrt(np.maximum(var, 1e-6)) * self.calibration
        return mu.astype(float), np.maximum(sigma.astype(float), 1e-3)

    def apply(self, compounds: list[Compound]) -> None:
        if not compounds:
            return
        mu, sigma = self.predict(compounds)
        for c, m, s in zip(compounds, mu, sigma):
            c.predicted = round(float(m), 3)
            c.uncertainty = round(float(s), 3)

    def report(self) -> dict:
        return {"model": "chemprop-dmpnn", "n_train": self.n_train,
                "epochs": self.epochs, "hidden": self.hidden, "depth": self.depth,
                "cv_r2": self.cv_r2, "cv_mae": self.cv_mae,
                "sigma_calibration": self.calibration,
                "note": "r2/mae here are in-sample, unlike PropertyModel's cross-validated values"}
