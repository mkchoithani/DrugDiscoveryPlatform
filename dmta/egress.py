"""Block every network connection that could leave this machine, and prove it.

A customer's compound structures are the asset. "Your molecules never leave your
machine" is the kind of claim that is trivial to write on a slide and impossible for
the person reading the slide to check, so this module exists to make it checkable:
run the Design-Make-Test cycle inside ``no_egress()`` and if any code tries to reach a
remote host, the call raises instead of succeeding.

WHAT COUNTS AS EGRESS, AND WHY THE LINE IS WHERE IT IS
------------------------------------------------------
Blocking *all* sockets would be a stronger-sounding guarantee and a worse test, because
it would fail for reasons that have nothing to do with data leaving:

  - ``AF_UNIX`` sockets cannot leave the host. They are a filesystem object. scikit-learn
    with ``n_jobs=-1`` spawns worker processes that talk over exactly these, so banning
    them would ban our own model fitting and the test would be measuring joblib.
  - Loopback (127.0.0.0/8, ::1) cannot leave the host either. The API server binds there;
    a customer running the UI locally is not exfiltrating anything.

So the rule is: **any connection to an address that is not loopback, over AF_INET or
AF_INET6, is egress.** That is exactly the set of connections capable of carrying a
molecule off the machine, and nothing else is blocked.

DNS is included. A hostname lookup does not carry a structure, but it is the step that
precedes one, and blocking it means a leak is caught at the first call rather than the
second -- with a message naming the host, which is more useful in a failure than
``connection refused``.

WHAT THIS DOES NOT PROVE
------------------------
This is a run-time check on one process. It does not prove the absence of a leak by
some channel it does not cover -- writing a file that something else uploads, a C
extension opening a raw socket without going through Python's ``socket`` module, or a
dependency that phones home at import time, before the guard is installed. It proves
that the cycle, as this process runs it, completes without contacting anything. That is
a real claim and a limited one, and ``egress_test.py`` states both.
"""
from __future__ import annotations

import contextlib
import ipaddress
import socket
from typing import Iterator

__all__ = ["EgressAttempt", "no_egress", "is_local_address"]


class EgressAttempt(RuntimeError):
    """Raised when code under ``no_egress()`` tries to reach a remote host."""


def is_local_address(host: object) -> bool:
    """True if *host* cannot carry data off this machine.

    Accepts what the socket module accepts: an IP string, a hostname, or a path for
    AF_UNIX. Anything it cannot parse as an address is treated as remote, because the
    safe default for an unrecognised destination is to block it and be told.
    """
    if not isinstance(host, (str, bytes)):
        return False
    if isinstance(host, bytes):
        try:
            host = host.decode()
        except UnicodeDecodeError:
            return False
    if host in ("localhost", "localhost.localdomain", "", "0.0.0.0", "::"):
        # "" and the unspecified addresses appear when binding a listener, which sends
        # nothing anywhere; an actual connection to them resolves to loopback.
        return True
    try:
        return ipaddress.ip_address(host.split("%")[0]).is_loopback
    except ValueError:
        return False


def _host_of(address: object) -> object:
    """Pull the host out of the several shapes a socket address takes."""
    if isinstance(address, tuple) and address:
        return address[0]
    return address          # AF_UNIX: a path, which is local by construction


@contextlib.contextmanager
def no_egress(allow: tuple[str, ...] = ()) -> Iterator[list[str]]:
    """Forbid connections to anything but loopback for the duration of the block.

    *allow* names hosts that may be contacted anyway -- for a caller deliberately
    exercising an opt-in fetch step. It is empty by default, because a default that
    allows something is a default that will be copied without being read.

    Yields the list of destinations that were allowed through, so a caller can assert
    on it rather than trusting that nothing happened.
    """
    contacted: list[str] = []
    real_connect = socket.socket.connect
    real_connect_ex = socket.socket.connect_ex
    real_create = socket.create_connection
    real_getaddrinfo = socket.getaddrinfo

    def check(host: object, where: str) -> None:
        if is_local_address(host) or (isinstance(host, str) and host in allow):
            if not is_local_address(host):
                contacted.append(str(host))
            return
        raise EgressAttempt(
            f"{where} tried to reach {host!r}, which is not loopback. "
            "Nothing in the Design-Make-Test cycle should contact a remote host; "
            "if this is a deliberate opt-in fetch, it belongs outside no_egress()."
        )

    def guarded_connect(self, address):                      # type: ignore[no-untyped-def]
        if self.family in (socket.AF_INET, socket.AF_INET6):
            check(_host_of(address), "socket.connect")
        return real_connect(self, address)

    def guarded_connect_ex(self, address):                   # type: ignore[no-untyped-def]
        if self.family in (socket.AF_INET, socket.AF_INET6):
            check(_host_of(address), "socket.connect_ex")
        return real_connect_ex(self, address)

    def guarded_create(address, *args, **kwargs):            # type: ignore[no-untyped-def]
        check(_host_of(address), "socket.create_connection")
        return real_create(address, *args, **kwargs)

    def guarded_getaddrinfo(host, *args, **kwargs):          # type: ignore[no-untyped-def]
        check(host, "socket.getaddrinfo")
        return real_getaddrinfo(host, *args, **kwargs)

    socket.socket.connect = guarded_connect                  # type: ignore[method-assign]
    socket.socket.connect_ex = guarded_connect_ex            # type: ignore[method-assign]
    socket.create_connection = guarded_create                # type: ignore[assignment]
    socket.getaddrinfo = guarded_getaddrinfo                 # type: ignore[assignment]
    try:
        yield contacted
    finally:
        socket.socket.connect = real_connect                 # type: ignore[method-assign]
        socket.socket.connect_ex = real_connect_ex           # type: ignore[method-assign]
        socket.create_connection = real_create               # type: ignore[assignment]
        socket.getaddrinfo = real_getaddrinfo                # type: ignore[assignment]
