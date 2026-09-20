# Test strategy

Unit tests cover AEOK transition denial and fail-closed release gates. Integrity/security tests cover SHA-256 upload records, traversal rejection, HMAC verification, and nonce replay rejection. Existing `engineering_recovery` tests remain regression coverage. Integration against real CAD/mesh/solver, physical evidence, GPT/DeepSeek providers, federation peers, production hosting, and live smoke tests require externally provisioned runtimes and credentials and are not simulated.

AEGIS tests now cover deterministic impact traversal, protocol compilation, explicit external worker blocking, physical-intake separation, unknown CAE/test correlation, independent provider-unavailable records, and contradiction blocking. The full repository regression command continues to require the OCP runtime for pre-existing B-Rep tests.
