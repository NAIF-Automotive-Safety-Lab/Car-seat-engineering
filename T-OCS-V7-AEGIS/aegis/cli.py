from __future__ import annotations
import sys
from .gate_status import main as gate_status_main

def main() -> int:
    if sys.argv[1:] == ["gate-status"]:
        return gate_status_main()
    print("usage: aegis gate-status", file=sys.stderr)
    return 2

if __name__ == "__main__":
    raise SystemExit(main())
