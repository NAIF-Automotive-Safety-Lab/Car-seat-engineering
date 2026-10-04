# Runbook
1. Verify R4.1 SHA-256.
2. Build native Project Chrono.
3. Run the generic engine probe twice and compare hashes.
4. Keep G2 BLOCKED until authoritative T-OCS joint/constraint closure is closed.
5. Keep G3 BLOCKED until mandatory physical inputs and OpenRadioss runtime are closed.
6. Run G4 hashing/ledger self-check.
7. Execute `aegis gate-status`.
Never manually promote a gate.