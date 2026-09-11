# Remote Runner API

Remote runners are an abstraction only until a trusted runner registry and signed transport are configured. A future runner must register an identity, version, capabilities, status, and heartbeat before dispatch.

A result must bind `RUN_ID`, request hash, input hashes, output hash, engine version, runner identity, timestamp, and signature. A runner cannot promote a result to PASS; AEGIS-X validates evidence, audit, and gate policy first. Unknown runners and unsigned results are `BLOCKED`.
