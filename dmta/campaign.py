"""The orchestrator: one object that runs Design → Make → Test and remembers everything."""
from __future__ import annotations

from dataclasses import dataclass, field
from functools import lru_cache

from .assay import AssayProvider, RetrospectiveAssay, SimulatedAssay
from .library import Compound, enumerate_library, group_by_intermediate
from .models import PropertyModel, regression_metrics
from .scoring import DEFAULT_PROFILE, Objective, profile_as_dicts, score_all
from .select import diversity_report, select_batch
from .store import Store


@lru_cache(maxsize=1)
def _base_library() -> tuple[Compound, ...]:
    return tuple(enumerate_library())


def fresh_library() -> list[Compound]:
    """A detached copy of the enumerated library — campaigns must not share state."""
    import copy
    return [copy.deepcopy(c) for c in _base_library()]


@dataclass
class CycleSettings:
    batch_size: int = 16
    acquisition: str = "ucb"
    beta: float = 0.25   # tuned: see README, exploration hurts at small budgets
    diversity: float = 0.15
    intermediate_bonus: float = 0.03
    budget: float | None = None
    max_per_scaffold: int | None = None
    assay_noise: float = 0.2
    oracle_seed: int = 42
    hit_threshold: float = 7.5

    def as_dict(self) -> dict:
        return self.__dict__.copy()


@dataclass
class Campaign:
    name: str
    target: str = "simulated kinase"
    settings: CycleSettings = field(default_factory=CycleSettings)
    profile: list[Objective] = field(default_factory=lambda: list(DEFAULT_PROFILE))
    store: Store | None = None
    assay: AssayProvider | None = None
    library: list[Compound] | None = None

    def __post_init__(self):
        if self.library is None:
            self.library = fresh_library()
        self.by_id = {c.id: c for c in self.library}
        self.model: PropertyModel | None = None
        self.history: list[dict] = []
        if self.assay is None:
            # a library carrying real measurements replays them; otherwise simulate
            if all(getattr(c, "truth", None) is not None for c in self.library):
                self.assay = RetrospectiveAssay(self.library)
            else:
                self.assay = SimulatedAssay(self.library, seed=self.settings.oracle_seed,
                                            noise=self.settings.assay_noise)
        if self.store is not None:
            existing = self.store.get_campaign(self.name)
            self.campaign_id = existing["id"] if existing else self.store.create_campaign(
                self.name, self.target, profile_as_dicts(self.profile), self.settings.as_dict())
            if existing:
                self._replay()
        else:
            self.campaign_id = None
        score_all(self.library, self.profile)

    # ------------------------------------------------------------ replay
    def _replay(self):
        for m in self.store.measurements(self.campaign_id):
            c = self.by_id.get(m["compound_id"])
            if c:
                c.measured = m["value"]
                c.cycle_made = m["cycle"]
        # Rebuild full cycle summaries, not just the batch rows: a replayed campaign
        # must present exactly the same history as a live one, or every consumer of
        # it (charts, reports, exports) sees missing fields.
        runs = {r["cycle"]: r for r in self.store.model_runs(self.campaign_id)}
        measured_by_cycle: dict[int, list] = {}
        for c in self.measured:
            measured_by_cycle.setdefault(c.cycle_made, []).append(c)

        self.history = []
        cumulative = 0
        for b in self.store.batches(self.campaign_id):
            batch = measured_by_cycle.get(b["cycle"], [])
            cumulative += len(batch)
            values = [c.measured for c in batch]
            best_so_far = max(
                (c.measured for c in self.measured
                 if c.cycle_made is not None and c.cycle_made <= b["cycle"]),
                default=None)
            run = runs.get(b["cycle"], {})
            self.history.append({
                "cycle": b["cycle"],
                "acquisition": b["acquisition"],
                "batch": [c.id for c in batch],
                "made": len(batch),
                "tested_total": cumulative,
                "library_size": len(self.library),
                "best_measured": round(best_so_far, 3) if best_so_far is not None else None,
                "batch_best": round(max(values), 3) if values else None,
                "hits_in_batch": sum(1 for v in values if v >= self.settings.hit_threshold),
                "plan": b["plan"],
                "diversity": b["diversity"],
                "prospective": run.get("prospective"),
                "model": run.get("report"),
                "batch_id": b["id"],
            })
        if self.measured:
            self._retrain()

    # -------------------------------------------------------- properties
    @property
    def measured(self) -> list[Compound]:
        return [c for c in self.library if c.measured is not None]

    @property
    def untested(self) -> list[Compound]:
        return [c for c in self.library if c.measured is None]

    @property
    def cycle(self) -> int:
        return len(self.history)

    # ------------------------------------------------------------ design
    def propose(self, size: int | None = None) -> list[Compound]:
        s = self.settings
        if not self.measured:
            # cold start: diverse, model-free — you cannot learn from nothing
            return select_batch(self.library, size or s.batch_size, acquisition="random",
                                beta=0, diversity=0.5, intermediate_bonus=0.0,
                                budget=s.budget, max_per_scaffold=s.max_per_scaffold,
                                seed=s.oracle_seed)
        return select_batch(self.library, size or s.batch_size, acquisition=s.acquisition,
                            beta=s.beta, diversity=s.diversity,
                            intermediate_bonus=s.intermediate_bonus, budget=s.budget,
                            max_per_scaffold=s.max_per_scaffold, seed=s.oracle_seed + self.cycle)

    # -------------------------------------------------------------- make
    def plan(self, batch: list[Compound]) -> dict:
        return group_by_intermediate(batch)

    # -------------------------------------------------------------- test
    def run_cycle(self, size: int | None = None, assay: AssayProvider | None = None) -> dict:
        batch = self.propose(size)
        if not batch:
            raise RuntimeError("library exhausted — nothing left to make")
        plan = self.plan(batch)
        div = diversity_report(batch)

        # predictions made BEFORE results: the only honest accuracy check
        prospective_pairs = [(c.predicted, c.id) for c in batch if c.predicted is not None]

        results = (assay or self.assay).measure(batch)
        cycle = self.cycle
        for c in batch:
            c.measured = results[c.id]
            c.cycle_made = cycle

        prospective = None
        if prospective_pairs:
            prospective = regression_metrics(
                [p for p, _ in prospective_pairs],
                [self.by_id[i].measured for _, i in prospective_pairs])

        batch_id = None
        if self.store:
            batch_id = self.store.add_batch(self.campaign_id, cycle, self.settings.acquisition, plan, div)
            self.store.add_measurements(self.campaign_id, batch_id, [
                {"compound_id": c.id, "smiles": c.smiles, "value": c.measured,
                 "predicted": c.predicted, "uncertainty": c.uncertainty} for c in batch
            ], source=(assay or self.assay).name)

        model_report = self._retrain()
        if self.store:
            self.store.add_model_run(self.campaign_id, cycle, model_report, prospective)

        summary = {
            "cycle": cycle,
            "acquisition": self.settings.acquisition if self.measured else "cold start",
            "batch": [c.id for c in batch],
            "made": len(batch),
            "tested_total": len(self.measured),
            "library_size": len(self.library),
            "best_measured": round(max(c.measured for c in self.measured), 3),
            "batch_best": round(max(c.measured for c in batch), 3),
            "hits_in_batch": sum(1 for c in batch if c.measured >= self.settings.hit_threshold),
            "plan": plan,
            "diversity": div,
            "prospective": prospective,
            "model": model_report,
            "batch_id": batch_id,
        }
        self.history.append(summary)
        return summary

    # ------------------------------------------------------------- learn
    def _retrain(self) -> dict:
        self.model = PropertyModel(seed=0).fit(self.measured)
        self.model.apply(self.library)
        score_all(self.library, self.profile)
        return self.model.report()

    # ------------------------------------------------------------ report
    def leaderboard(self, n: int = 20, tested_only: bool = True) -> list[Compound]:
        pool = self.measured if tested_only else self.library
        key = (lambda c: c.measured) if tested_only else (lambda c: c.score or 0)
        return sorted(pool, key=key, reverse=True)[:n]

    def state(self) -> dict:
        m = self.measured
        return {
            "name": self.name, "target": self.target, "cycle": self.cycle,
            "library_size": len(self.library), "tested": len(m),
            "best_measured": round(max((c.measured for c in m), default=0), 3) if m else None,
            "hits": sum(1 for c in m if c.measured >= self.settings.hit_threshold),
            "settings": self.settings.as_dict(),
            "profile": profile_as_dicts(self.profile),
            "model": self.model.report() if self.model else None,
            "history": [{k: v for k, v in h.items() if k != "batch"} for h in self.history],
        }
