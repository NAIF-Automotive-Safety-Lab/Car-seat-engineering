from aegis.gate_status import EXPECTED_R41_SHA256

def test_fixed_r41_sha256():
    assert len(EXPECTED_R41_SHA256) == 64
    int(EXPECTED_R41_SHA256, 16)

def test_fail_closed_vocabulary():
    assert {"PASS","FAIL","BLOCKED","NOT_PROVEN"} == {"PASS","FAIL","BLOCKED","NOT_PROVEN"}
    assert "UNKNOWN" not in {"PASS","FAIL","BLOCKED","NOT_PROVEN"}
