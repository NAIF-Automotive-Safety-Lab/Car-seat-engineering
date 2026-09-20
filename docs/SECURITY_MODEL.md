# Security model

Mutating and ledger-read APIs require `Authorization: Bearer <AEGIS_AUTH_TOKEN>`; when this server-side variable is absent the APIs fail closed with `AUTH_CONFIGURATION_REQUIRED`. The public health endpoint discloses configuration state only. Uploads use basename-only paths, extension allowlisting, size limits, PDF signature checks, byte hashes, and do not execute uploaded content. The static UI contains no keys. Native processing must occur in a separately sandboxed worker; this repository does not claim such a worker is provisioned.

GPT and DeepSeek URLs, keys, and model IDs are read only from server environment variables (`AEGIS_GPT_*` and `AEGIS_DEEPSEEK_*`). Provider failures are persisted as external-runtime-required rather than exposed or retried as fabricated content.
