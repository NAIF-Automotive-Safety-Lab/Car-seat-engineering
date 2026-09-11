# Git Recovery

- LOCAL COMMIT: `aaa037a65169592cf975c8df66015c038e93e127`
- BRANCH: `main`
- REMOTE: `https://7f0548f6b28a1758d9bf7be630d36cb4.artifacts.cloudflare.net/git/prod/iMguyThUG2dP8xXixCH9da.git`
- EXPECTED REMOTE HEAD: `aaa037a65169592cf975c8df66015c038e93e127`
- LOCAL HEAD: verified locally as the commit above
- REMOTE HEAD: NOT_PROVEN; direct `git push` failed because shell credentials are unavailable
- WORKTREE: clean at the preceding commit before adding this recovery package
- MANIFEST: `V7_AEGIS_CANONICAL_MANIFEST.json`
- TEST STATE: `V7_AEGIS_TEST_STATE.json`
- BLOCKERS: `V7_AEGIS_BLOCKERS.json`

Recovery sequence when authorization is available: validate, run `pnpm test`, confirm `git diff --check`, push `main`, fetch `origin/main`, then require exact equality between local HEAD and remote HEAD. Do not call GitHub synchronization PASS before that equality is observed.
