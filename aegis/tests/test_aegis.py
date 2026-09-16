import tempfile
from aegis.core import AEOK, AegisStore, Lifecycle
from aegis.federation import sign, verify_handshake
from aegis.ingestion import validate_upload

def test_aeok_only_allows_ordered_transitions():
 s=AegisStore(":memory:"); k=AEOK(s); r=s.create("Requirement", {"text":"rail load"}, "architect", status="DRAFT")
 a=k.transition(r["id"], Lifecycle.ARCHITECTED, "architect")
 assert a["status"] == "ARCHITECTED"
 try: k.transition(r["id"], Lifecycle.RELEASED, "architect")
 except ValueError: pass
 else: assert False

def test_release_gate_is_fail_closed():
 s=AegisStore(":memory:"); result=AEOK(s).release_gate("candidate", {"build":"PASS"}, "auditor")
 assert result["status"] == "BLOCK" and result["decision"] == "BLOCK"

def test_upload_hash_and_path_protection():
 item=validate_upload("report.json", b'{"a":1}', "application/json")
 assert len(item.sha256) == 64
 try: validate_upload("../x.json", b"{}", "application/json")
 except ValueError: pass
 else: assert False

def test_signed_federation_replay_is_blocked():
 s=AegisStore(":memory:"); e={"source_id":"a","destination_endpoint":"https://example.invalid","protocol_version":"1","request_id":"r","timestamp":"2026-01-01T00:00:00+00:00","nonce":"n","payload_hash":"abc"}; e["signature"]=sign(e,"secret")
 assert verify_handshake(e,"secret",s)=="ESTABLISHED"
 assert verify_handshake(e,"secret",s)=="BLOCKED"
from aegis.workflows import change_impact, compile_protocol, intake_physical, correlate, review, contradiction, WorkerFabric

def test_change_impact_protocol_and_external_reviews_are_persisted():
 s=AegisStore(":memory:")
 req=s.create("Requirement", {"acceptance_criteria":["force <= released limit"]}, "architect", status="DRAFT")
 protocol=compile_protocol(s, req["id"], "test")
 cad=s.create("CADArtifact", {"name":"rail"}, "cad", input_refs=[req["content_hash"]])
 change=s.create("ChangeSet", {"changed_refs":[req["content_hash"]]}, "architect")
 impact=change_impact(s, change["id"], "architect")
 assert protocol["physical_test_status"] == "PHYSICAL_EVIDENCE_REQUIRED"
 assert cad["id"] in {x["id"] for x in impact["affected"]["cad"]}
 gpt=review(s,"GPT",[req["content_hash"]],"audit", "reviewer")
 deep=review(s,"DEEPSEEK",[req["content_hash"]],"audit", "reviewer")
 assert gpt["status"] == deep["status"] == "EXTERNAL_RUNTIME_REQUIRED"
 assert contradiction(s,gpt["id"],deep["id"],"auditor")["status"] == "BLOCK"

def test_worker_and_physical_correlation_do_not_promote_unknown_data():
 s=AegisStore(":memory:"); job=WorkerFabric(s).create("solver", [], "cae")
 assert WorkerFabric(s).external_runtime_required(job["id"],"CalculiX","cae")["status"] == "BLOCKED"
 measurement=intake_physical(s,{"test_id":"T-1","source":"lab","captured_at":"2026-09-16T00:00:00Z","measurement_refs":["sensor-1"],"filename":"run.json","content_hex":"7b7d"},"test")
 simulation=s.create("CAEResult", {"kind":"SIMULATION"}, "cae")
 result=correlate(s,simulation["id"],measurement["id"],"cae")
 assert result["correlation_status"] == "UNKNOWN"
