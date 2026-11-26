#!/usr/bin/env python3
"""Does the loop run with the network cut off, and is the claim worth making?

"Your molecules never leave your machine" is the one competitive property this platform
already has and has never claimed. Palo Alto Networks paid $4.5B for a company whose
pitch is substantially that sentence, so the demand is real; what is missing here is any
way for a customer to check it.

This runs the complete Design-Make-Test cycle inside dmta.egress.no_egress(), which
raises on any connection to a non-loopback address. Three things have to hold, and the
third is the one that makes it evidence rather than a slogan:

  1. the cycle completes, so the claim is true today;
  2. the guard actually fires when something does try to leave, so a passing run means
     "nothing was contacted" and not "the guard is broken";
  3. the guard is not quietly wider than the claim -- it must still permit the AF_UNIX
     and loopback traffic that scikit-learn's worker processes use, or it would be
     proving something about joblib rather than about us.

THE TRAP THIS IS BUILT AROUND
-----------------------------
A test that cuts the network and passes proves nothing on its own. It passes identically
if the guard is installed, if the guard is a no-op, and if the cycle silently skipped the
work. So this script runs a positive control: a deliberate connection to a remote address,
which MUST raise. If the control does not raise, the whole run is void and says so.

The second trap is the import boundary. The guard is installed at call time, so anything
a dependency does at import time -- a licence check, a telemetry ping -- happens before
it and is invisible here. The cycle is therefore exercised with every module already
imported, which is the honest arrangement: it measures the cycle, and the summary says
so rather than implying more.

Run:  python3 egress_test.py
"""
from __future__ import annotations

import ast
import socket
import sys
import time
from pathlib import Path

from dmta.campaign import Campaign, CycleSettings
from dmta.egress import EgressAttempt, no_egress

# The two modules allowed to reach the network, both explicit one-off fetch steps the
# user invokes by name. Every other module must be clean, and the scan below fails if
# that stops being true -- which is the part that keeps the claim true tomorrow.
FETCHERS = {"chembl.py", "structure.py"}
NETWORK_CALLS = {"urlopen", "urlretrieve", "Request", "create_connection"}
NETWORK_MODULES = {"requests", "httpx", "urllib", "aiohttp", "socket", "ftplib", "http"}


def scan_for_network_use(package: Path) -> dict[str, list[str]]:
    """Static scan: which modules import or call something that can reach a network.

    Deliberately crude and deliberately over-sensitive. A false positive costs someone
    thirty seconds of reading; a false negative costs the claim.
    """
    found: dict[str, list[str]] = {}
    for path in sorted(package.glob("*.py")):
        if path.name == "egress.py":
            continue                      # this module imports socket on purpose
        hits: list[str] = []
        tree = ast.parse(path.read_text())
        for node in ast.walk(tree):
            if isinstance(node, ast.Import):
                for a in node.names:
                    if a.name.split(".")[0] in NETWORK_MODULES:
                        hits.append(f"import {a.name}")
            elif isinstance(node, ast.ImportFrom):
                if node.module and node.module.split(".")[0] in NETWORK_MODULES:
                    hits.append(f"from {node.module} import ...")
            elif isinstance(node, ast.Call):
                f = node.func
                name = getattr(f, "attr", None) or getattr(f, "id", None)
                if name in NETWORK_CALLS:
                    hits.append(f"{name}() at line {node.lineno}")
        if hits:
            found[path.name] = hits
    return found


def control_the_guard_actually_fires() -> bool:
    """The positive control. Without this the main result is unfalsifiable."""
    with no_egress():
        try:
            socket.create_connection(("93.184.216.34", 80), timeout=1)
        except EgressAttempt:
            return True
        except Exception as exc:                      # pragma: no cover - diagnostic
            print(f"  control raised {type(exc).__name__}, not EgressAttempt: {exc}")
            return False
    return False


def control_local_traffic_still_works() -> bool:
    """The guard must not be wider than the claim: loopback and AF_UNIX stay open.

    This opens real sockets rather than asserting on is_local_address, because the
    first version of this control did the latter and was wrong in a way that still
    reported a failure: for an AF_UNIX socket the guard never calls is_local_address
    at all -- it checks the socket family first and returns -- so the control was
    exercising a code path that does not run. A control has to use the thing the
    customer uses, or it is testing a helper nobody called.
    """
    with no_egress():
        try:
            listener = socket.socket(socket.AF_INET, socket.SOCK_STREAM)
            listener.bind(("127.0.0.1", 0))
            listener.listen(1)
            client = socket.create_connection(listener.getsockname(), timeout=2)
            client.close()
            listener.close()
        except EgressAttempt:
            print("  loopback was blocked — the UI would not run under this guard")
            return False
        try:
            a, b = socket.socketpair(socket.AF_UNIX, socket.SOCK_STREAM)
            a.close(); b.close()
        except EgressAttempt:
            print("  AF_UNIX was blocked — scikit-learn's workers would not run")
            return False
        except (AttributeError, OSError):
            pass                     # no AF_UNIX on this platform; not our failure
    return True


def main() -> int:
    print("Can the Design-Make-Test cycle run with the network cut off?\n")

    print("Static scan of dmta/ for anything that can reach a network")
    found = scan_for_network_use(Path("dmta"))
    for name, hits in found.items():
        mark = "opt-in fetcher" if name in FETCHERS else "*** UNEXPECTED ***"
        print(f"  {name:<22} {mark:<18} {'; '.join(hits[:3])}")
    unexpected = set(found) - FETCHERS
    missing = FETCHERS - set(found)
    if missing:
        print(f"\n  NOTE: {sorted(missing)} no longer reach the network. "
              "Update FETCHERS; the claim gets stronger, not weaker.")
    print(f"\n  {len(found)} of {len(list(Path('dmta').glob('*.py')))} modules touch the "
          f"network, {len(unexpected)} of them unexpectedly.\n")

    print("Controls")
    fires = control_the_guard_actually_fires()
    local_ok = control_local_traffic_still_works()
    print(f"  guard raises on a remote address:        {'yes' if fires else 'NO — RESULT IS VOID'}")
    print(f"  loopback and AF_UNIX still permitted:    {'yes' if local_ok else 'NO — too wide'}")
    if not fires:
        print("\nThe guard did not fire on a deliberate remote connection. Nothing below "
              "means anything, because a no-op guard passes every test.")
        return 2
    print()

    print("The cycle itself, with every non-loopback address blocked")
    t0 = time.time()
    try:
        with no_egress() as contacted:
            c = Campaign("f13-egress", settings=CycleSettings(batch_size=12, oracle_seed=7))
            summaries = [c.run_cycle(), c.run_cycle(), c.run_cycle()]
    except EgressAttempt as exc:
        print(f"  FAILED — the cycle tried to leave the machine:\n    {exc}")
        return 1
    secs = time.time() - t0

    last = summaries[-1]
    print(f"  3 cycles completed in {secs:.1f}s")
    print(f"  compounds made and tested: {last['tested_total']} of {last['library_size']}")
    print(f"  model retrained each cycle, prospective check: "
          f"{'present' if last.get('prospective') else 'cold start only'}")
    print(f"  remote hosts contacted: {len(contacted)}")
    print()

    ok = fires and local_ok and not unexpected
    print("=" * 78)
    if ok:
        print("CLAIM SUPPORTED, within the limits below.")
        print("  The full cycle -- enumerate, score, select, plan, measure, retrain --")
        print("  ran to completion with every non-loopback connection blocked, and the")
        print("  guard is known to fire because a deliberate remote connection raised.")
    else:
        print("CLAIM NOT SUPPORTED.")
        if unexpected:
            print(f"  Unexpected network use in: {sorted(unexpected)}")
        if not local_ok:
            print("  The guard blocks local traffic too, so it proves more than we claim")
            print("  and would break a customer running the UI.")
    print()
    print("WHAT THIS DOES NOT SHOW, stated so nobody has to find it")
    print("  - One process, run-time only. Import-time behaviour of dependencies happens")
    print("    before the guard is installed and is not covered.")
    print("  - Python's socket module only. A C extension opening a raw socket would not")
    print("    be seen.")
    print("  - Nothing about disk. A molecule written to a file that some other process")
    print("    uploads has still left, and no test here would notice.")
    print("  - The two fetchers above DO reach the network when you invoke them. The")
    print("    claim is about the cycle, not about the whole program.")
    return 0 if ok else 1


if __name__ == "__main__":
    sys.exit(main())
