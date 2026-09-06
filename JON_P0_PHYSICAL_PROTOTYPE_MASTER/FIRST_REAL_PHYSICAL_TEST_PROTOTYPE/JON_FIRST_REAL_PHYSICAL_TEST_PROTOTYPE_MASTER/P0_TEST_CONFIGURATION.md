# P0 Test Configuration

Every run gets a frozen configuration manifest. Minimum identity:

`P0 article -> design revision -> serial -> configuration ID -> module serials -> fixture ID -> sensor/calibration IDs -> test command -> operator -> date/time`

## Configuration rules
1. P0-A/B/C/D isolate one acquisition function or subsystem.
2. P0-E is the integrated mechanism.
3. P0-A/B/C/D results are established before integrated dynamic screening.
4. P0-B/C/D modules are interchangeable so A/B comparisons change one controlled variable at a time.
5. P0-E results marked MULTIVARIABLE are evidence of combined behavior only; they are not used for causal attribution.

## Dynamic command
Use a laboratory-controlled actuator/sled/launcher command whose actual motion/force is measured. Do not label a laboratory command as a vehicle crash pulse unless it is sourced from an authoritative external test and validated for the intended case.
