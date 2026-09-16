import pytest
from aegis.kernel import AEOK, Artifact, Decision, EvidenceGraph, Phase, ReviewFinding

H = "a" * 64

def artifact(i, kind): return Artifact(i, H, kind, f"store://{i}", "worker", "signer")

def test_graph_is_canonical_and_rejects_overwrite():
    graph = EvidenceGraph(); graph.add_artifact(artifact("a", "native-cad"))
    assert graph.digest() == graph.digest()
    with pytest.raises(ValueError): graph.add_artifact(Artifact("a", "b" * 64, "native-cad", "store://a", "worker"))

def test_release_fails_closed_without_evidence():
    kernel = AEOK(); event = kernel.transition(Phase.RELEASE, Decision.PASS, "request release")
    assert event.decision is Decision.BLOCK
    assert "missing required evidence" in event.reason

def test_release_blocks_unresolved_critical_review():
    kernel = AEOK()
    for index, kind in enumerate(kernel.policy.REQUIRED_KINDS): kernel.record_artifact(artifact(str(index), kind))
    kernel.add_finding(ReviewFinding("f1", "independent-auditor", ("0",), "CRITICAL", "contradiction"))
    assert kernel.transition(Phase.RELEASE, Decision.PASS, "request release").decision is Decision.BLOCK

def test_transition_rejects_unknown_artifact():
    with pytest.raises(ValueError): AEOK().transition(Phase.CAD, Decision.PASS, "x", ("missing",))
