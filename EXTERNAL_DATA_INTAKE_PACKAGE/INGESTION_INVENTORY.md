# External evidence intake inventory

## Hosts
- `V4_EVIDENCE/` — extracted V4 evidence
- `EXTERNAL_DATA_INTAKE_PACKAGE/` — secure canonical intake, legacy comparison, and quarantined invalid artifacts

## Source archive SHA256 (archives not stored)
- `V4_EVIDENCE_f562a8c..zip`: `235f06ae6493ac0beaafad44b8faef7402e3ed2f175890673f4218540c4166f4`
- `EXTERNAL_DATA_INTAKE_PACKAGE_SECURE..zip`: `90c83798bc6e39135578863d4fb1dda828bc59638088b66ab7bdc482c29d42e0`
- `EXTERNAL_DATA_INTAKE_PACKAGE..zip`: `47599e8a2abd63cf1d719997a56fe28cd7beb209290bcd1626650d50ca292c20`

R4.1 baseline SHA256: `fbe6b17cdbf728a2e47963e567e12eeceb1352a36e719e7d1c55cc5f712a0a68`

| Path | Bytes | SHA256 |
|---|---:|---|
| `EXTERNAL_DATA_INTAKE_PACKAGE/EXTERNAL_DATA_INTAKE_PACKAGE.json` | 20959 | `03298b0c5bcc18d5a36b1b8d0f3b2032ff07c8a69d3f33709299f689f38f892b` |
| `EXTERNAL_DATA_INTAKE_PACKAGE/INGESTION_INVENTORY.md` | 9660 | `911ac4af6c5981cf4b617c6cba7a67f282302776c7fadb04db95c31a4cc0e223` |
| `EXTERNAL_DATA_INTAKE_PACKAGE/legacy/EXTERNAL_DATA_INTAKE_PACKAGE.json` | 20959 | `03298b0c5bcc18d5a36b1b8d0f3b2032ff07c8a69d3f33709299f689f38f892b` |
| `EXTERNAL_DATA_INTAKE_PACKAGE/legacy/initial_gate_result.json` | 3223 | `a302e8ea8f2b31232697c144cb5ad9b1557acee2c494cd8f4a9c0dc37ff2c8b0` |
| `EXTERNAL_DATA_INTAKE_PACKAGE/legacy/schemas/input_record.schema.json` | 506 | `3c901b98b4db421777d7d0ef4c455cb6674bbe6e9116c731ad209215ea203fdb` |
| `EXTERNAL_DATA_INTAKE_PACKAGE/legacy/source_data_required/bolt_stiffness.json` | 1269 | `89654d565ce57d9ea448fa4ff0c3c40311fed010dd9d8ca126255bdfad3db48d` |
| `EXTERNAL_DATA_INTAKE_PACKAGE/legacy/source_data_required/cg_properties.json` | 1286 | `b0375eea56badc65b06557ed28b9438cef49f7388f2ed24ac808ac9aa20935bb` |
| `EXTERNAL_DATA_INTAKE_PACKAGE/legacy/source_data_required/inertia_tensor.json` | 1269 | `2ec24248bdafbd6b2125e480bb7b91b47cda0f05c1f8f544698a09b4f1da1c36` |
| `EXTERNAL_DATA_INTAKE_PACKAGE/legacy/source_data_required/initial_velocity.json` | 1269 | `692ebff0712c48c237c9db813e1a514ee89b7406dc577c0622d40d095e744a45` |
| `EXTERNAL_DATA_INTAKE_PACKAGE/legacy/source_data_required/mass_properties.json` | 1286 | `46186a6dafe565b97739f78b425f0c279c29a92df0fda1c7a45938ee8b7e1f88` |
| `EXTERNAL_DATA_INTAKE_PACKAGE/legacy/source_data_required/material_density.json` | 1254 | `ff6cd137765f939e8afa1a66ff70a57e1d3b941cde4a5bde4f3fc87751e063b0` |
| `EXTERNAL_DATA_INTAKE_PACKAGE/legacy/source_data_required/vehicle_pulse.json` | 1280 | `0b85cab25a0d2f99116a8866b79630ae36e869bdf7fdb22cd28fa7face199edd` |
| `EXTERNAL_DATA_INTAKE_PACKAGE/legacy/test_data_required/absorber_Fv.json` | 1255 | `494a52dce1e6da62d5f683d546b89ce9a40b5cdca699691a06cf1488ed5242c6` |
| `EXTERNAL_DATA_INTAKE_PACKAGE/legacy/test_data_required/absorber_Fx.json` | 1258 | `1bd8bfe7dd1c4b8988513d792ecdb4b9b7aa1bd9f52a48187b8ac6b6f7f126c4` |
| `EXTERNAL_DATA_INTAKE_PACKAGE/legacy/test_data_required/contact_damping.json` | 1294 | `87788a32164e86900dec5566ede7b5ca0e4c650fe8e7fb41cd25d8df0d952406` |
| `EXTERNAL_DATA_INTAKE_PACKAGE/legacy/test_data_required/contact_stiffness.json` | 1252 | `16d03606905f652494c4bfa9c1bef2d7f10884f2985340eb6cc1b5c9ff4f2d14` |
| `EXTERNAL_DATA_INTAKE_PACKAGE/legacy/test_data_required/friction.json` | 1249 | `67b1ce5c1b7d1a5b008c131205c997dee80260ff0fdea57eeafe1b21cdbd6dcf` |
| `EXTERNAL_DATA_INTAKE_PACKAGE/legacy/test_data_required/joint_compliance.json` | 1288 | `79ef22545aaf907e48d5caf3f9b56cf538503c1595596493483166cdd296c9f2` |
| `EXTERNAL_DATA_INTAKE_PACKAGE/legacy/test_data_required/restitution.json` | 1237 | `8b65661c66eda1d311d86deb17b233403a9edebd99ae2fb1a4dce86e54b2b902` |
| `EXTERNAL_DATA_INTAKE_PACKAGE/legacy/validate_intake.py` | 1034 | `d862d6493a152a117754418f5b60e1475a4c101b61e0ed118e3c5e0291003934` |
| `EXTERNAL_DATA_INTAKE_PACKAGE/quarantine/INVALID_JSON/EXTERNAL_DATA_INGESTION_AUDIT.json` | 5265 | `dcbf23d43845ccce93c5032fefc17b232a1847cbaeaa605f3e69038bf4c38aec` |
| `EXTERNAL_DATA_INTAKE_PACKAGE/quarantine/INVALID_JSON/PHYSICAL_INPUT_CLOSURE_MATRIX.json` | 5050 | `4c7e9f2f5281a9c89037fe48dd9c20dfed912804cef61a2d81a4fe83c025d58f` |
| `EXTERNAL_DATA_INTAKE_PACKAGE/quarantine/INVALID_JSON/REJECTED_ARTIFACTS_REGISTER.json` | 4982 | `d57db910f0b7f5c89d5495cc434662309a02e29e722babc07a2146478d0a7d05` |
| `EXTERNAL_DATA_INTAKE_PACKAGE/quarantine/INVALID_JSON/VALIDATION_REPORT.json` | 1113 | `fe6db4840fd755eb5e19074650cb57b30128f6b7c022cdb3114b533ed72f1941` |
| `EXTERNAL_DATA_INTAKE_PACKAGE/schemas/input_record.schema.json` | 506 | `3c901b98b4db421777d7d0ef4c455cb6674bbe6e9116c731ad209215ea203fdb` |
| `EXTERNAL_DATA_INTAKE_PACKAGE/source_data_required/bolt_stiffness.json` | 1269 | `89654d565ce57d9ea448fa4ff0c3c40311fed010dd9d8ca126255bdfad3db48d` |
| `EXTERNAL_DATA_INTAKE_PACKAGE/source_data_required/cg_properties.json` | 1286 | `b0375eea56badc65b06557ed28b9438cef49f7388f2ed24ac808ac9aa20935bb` |
| `EXTERNAL_DATA_INTAKE_PACKAGE/source_data_required/inertia_tensor.json` | 1269 | `2ec24248bdafbd6b2125e480bb7b91b47cda0f05c1f8f544698a09b4f1da1c36` |
| `EXTERNAL_DATA_INTAKE_PACKAGE/source_data_required/initial_velocity.json` | 1269 | `692ebff0712c48c237c9db813e1a514ee89b7406dc577c0622d40d095e744a45` |
| `EXTERNAL_DATA_INTAKE_PACKAGE/source_data_required/mass_properties.json` | 1286 | `46186a6dafe565b97739f78b425f0c279c29a92df0fda1c7a45938ee8b7e1f88` |
| `EXTERNAL_DATA_INTAKE_PACKAGE/source_data_required/material_density.json` | 1254 | `ff6cd137765f939e8afa1a66ff70a57e1d3b941cde4a5bde4f3fc87751e063b0` |
| `EXTERNAL_DATA_INTAKE_PACKAGE/source_data_required/vehicle_pulse.json` | 1280 | `0b85cab25a0d2f99116a8866b79630ae36e869bdf7fdb22cd28fa7face199edd` |
| `EXTERNAL_DATA_INTAKE_PACKAGE/test_data_required/absorber_Fv.json` | 1255 | `494a52dce1e6da62d5f683d546b89ce9a40b5cdca699691a06cf1488ed5242c6` |
| `EXTERNAL_DATA_INTAKE_PACKAGE/test_data_required/absorber_Fx.json` | 1258 | `1bd8bfe7dd1c4b8988513d792ecdb4b9b7aa1bd9f52a48187b8ac6b6f7f126c4` |
| `EXTERNAL_DATA_INTAKE_PACKAGE/test_data_required/contact_damping.json` | 1294 | `87788a32164e86900dec5566ede7b5ca0e4c650fe8e7fb41cd25d8df0d952406` |
| `EXTERNAL_DATA_INTAKE_PACKAGE/test_data_required/contact_stiffness.json` | 1252 | `16d03606905f652494c4bfa9c1bef2d7f10884f2985340eb6cc1b5c9ff4f2d14` |
| `EXTERNAL_DATA_INTAKE_PACKAGE/test_data_required/friction.json` | 1249 | `67b1ce5c1b7d1a5b008c131205c997dee80260ff0fdea57eeafe1b21cdbd6dcf` |
| `EXTERNAL_DATA_INTAKE_PACKAGE/test_data_required/joint_compliance.json` | 1288 | `79ef22545aaf907e48d5caf3f9b56cf538503c1595596493483166cdd296c9f2` |
| `EXTERNAL_DATA_INTAKE_PACKAGE/test_data_required/restitution.json` | 1237 | `8b65661c66eda1d311d86deb17b233403a9edebd99ae2fb1a4dce86e54b2b902` |
| `EXTERNAL_DATA_INTAKE_PACKAGE/validate_intake.py` | 7216 | `1e9304e878d26afe7d9736111d9da7e70cde9e08905bee8f9d54d0eaa1fc8563` |
| `V4_EVIDENCE/01_BASELINE/baseline.json` | 409 | `b3fe4c509685a266e8cfa61b98b71d527e4750a17875efdacbf87a9fd4ea4efd` |
| `V4_EVIDENCE/02_PHYSICAL_INPUTS/physical_input_register.json` | 8402 | `65c8745c19afbe8ca3d10fe9869b3d74b289f8d74743613197f1f10201163b84` |
| `V4_EVIDENCE/03_FE_RUNTIME/beamcontact.cvg` | 1176 | `a1c0ef90899697cc7a6c2d79301b03e23d6c17849499eacca829e392cc28f70c` |
| `V4_EVIDENCE/03_FE_RUNTIME/beamcontact.frd` | 46829 | `4a799603e39d95cc36c3c6be3e94010b717401a18fe54f74cfea764046c4cea6` |
| `V4_EVIDENCE/03_FE_RUNTIME/beamcontact.inp` | 24282 | `931ec3e5a6ceded6a0c0529b6da5d160071ed4ed9294afafbbbf0418df3ccc8e` |
| `V4_EVIDENCE/03_FE_RUNTIME/beamcontact.sta` | 242 | `182bf7c2322bc2754b3815d276863e84da04c4931a24d9ba830b8c8f437b1bef` |
| `V4_EVIDENCE/03_FE_RUNTIME/ccx.stderr` | 0 | `e3b0c44298fc1c149afbf4c8996fb92427ae41e4649b934ca495991b7852b855` |
| `V4_EVIDENCE/03_FE_RUNTIME/ccx.stdout` | 10552 | `53e1366120ab3b11f1be4b53320214e5e1315016de23487c097c4d5ff5843259` |
| `V4_EVIDENCE/03_FE_RUNTIME/project_chrono_manifest.json` | 592 | `cb490c790aa5106da63b8c2d5a1c72b2eee6ee01fcdb8b48a573f4987a2dd05e` |
| `V4_EVIDENCE/03_FE_RUNTIME/solver_manifest.json` | 1113 | `2049d65ad2d6bb6b9809f389a2276252098df26425337624031bc848a6d5fba4` |
| `V4_EVIDENCE/04_FE_MODEL/model_readiness.json` | 780 | `8f2c17e918d678f0567e9b6dfafa0e39c3cc91fcbb6f0065572e4299b0405e6d` |
| `V4_EVIDENCE/05_FE_01/execution.json` | 460 | `dd00aa598977ed9a0759d9e42d40d48cda47972897b759128963692ae87d26db` |
| `V4_EVIDENCE/06_FE_02/execution.json` | 460 | `32f90d50fe80f7552b62b6238b78daa09ca36222061e1366c383628a29185b3e` |
| `V4_EVIDENCE/07_FE_03/execution.json` | 460 | `505f25a90651135f2b739be22cfd1a137a61c0fe2eb6ca7dd1380397637251ee` |
| `V4_EVIDENCE/08_FE_04/execution.json` | 460 | `2946c95e09648aaca5b82a5156cde4ac508dfad42cf0e7a7f9ba83a179307b15` |
| `V4_EVIDENCE/09_FE_05/execution.json` | 460 | `cc4695939c70e5e4852a28b6035e12a02f2aa7a8b8c31f81d35b77ae2a193707` |
| `V4_EVIDENCE/10_FE_06/execution.json` | 460 | `f36797ef26f4ee7312033d3be70a9ccbcd2d5510740a9831349fde5ed55c5405` |
| `V4_EVIDENCE/11_FE_07/execution.json` | 460 | `5e532f49536832f25c0ed368a6b403f1787e9be8ab8419b772a6331b21b2bb4a` |
| `V4_EVIDENCE/12_FE_08/execution.json` | 460 | `64b4aa8f73f8078c14f9d3f5be2f5bb466f2f168f7cb7a39ba48292fb399f83d` |
| `V4_EVIDENCE/13_VALIDITY/validity.json` | 182 | `717e0f182d1134083ec2c17203b263a6df73504cf307f2b9d99f43ef80ae828e` |
| `V4_EVIDENCE/14_FAILURE/failure_determination.json` | 138 | `d980ffb7944adc1ce7bff4a4a322a530c8d71f98cee81c914fad30258989fac2` |
| `V4_EVIDENCE/15_ROOT_CAUSE/root_cause.json` | 98 | `1b633fc67c4de73986c02dbe2dcc1c23b9e7194923f593e593d0cb5aa8e9e0cc` |
| `V4_EVIDENCE/16_GATES/gates.json` | 527 | `f05f9cee593937aba441637f0c2c9087319958374c95839003b5092518871902` |
| `V4_EVIDENCE/17_MANIFEST/package_manifest.json` | 2742 | `d241348ac964ccadd63011496e3d7d9d0d8dc6c668af68c7fffe18f0a5c14061` |
| `V4_EVIDENCE/17_MANIFEST/sha256_manifest.json` | 2343 | `ebf47dc4cd9f992d29cb76b2043a4f585beb4a1717df24d268c8656b453a35d9` |
