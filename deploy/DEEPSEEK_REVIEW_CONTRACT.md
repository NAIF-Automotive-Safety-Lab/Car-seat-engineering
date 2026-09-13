# DeepSeek integration contract

DeepSeek is used only as the advisory adversarial reviewer. It must never be treated as validation evidence, truth-state authority, merge authority, or release authority.

## API
- OpenAI-compatible base URL: `https://api.deepseek.com`
- Key variable: `DEEPSEEK_API_KEY`
- Default development model in this bootstrap: `deepseek-v4-flash`

## Review payload binding
Every review request must include:
- canonical Git SHA
- R4.1 SHA-256
- V7 package SHA-256
- review purpose
- explicit evidence boundary

The server must reject a packet whose hashes do not match the canonical values configured on the host.

## Required output classes
`VERIFIED_FACTS`, `INFERENCES`, `GAPS`, `CONTRADICTIONS_RISKS`, `REQUIRED_EVIDENCE`, `NEXT_GATE`.

The review result remains `UNVERIFIED` and advisory-only.
