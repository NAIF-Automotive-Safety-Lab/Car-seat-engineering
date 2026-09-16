"""Executable, deterministic AEGIS workflow services; external execution is never simulated."""
from __future__ import annotations
import hashlib, json, os, urllib.request, uuid
from dataclasses import asdict, dataclass
from typing import Any
from .core import AegisStore, now, canonical_json
from .ingestion import validate_upload

JOB_STATES = frozenset({"QUEUED", "RUNNING", "SUCCEEDED", "FAILED", "BLOCKED", "CANCELLED"})

class WorkerFabric:
    def __init__(self, store: AegisStore): self.store = store
    def create(self, kind: str, input_refs: list[str], actor: str) -> dict[str, Any]:
        return self.store.create("WorkerJob", {"kind":kind,"worker_id":None,"input_refs":input_refs,"output_refs":[],"started_at":None,"ended_at":None,"failure_reason":None}, actor, status="QUEUED", input_refs=input_refs)
    def external_runtime_required(self, job_id: str, runtime: str, actor: str) -> dict[str, Any]:
        job=self.store.get(job_id)
        if not job or job["status"] != "QUEUED": raise ValueError("job is not queueable")
        return self.store.create("WorkerJob", {"kind":job["kind"],"worker_id":None,"input_refs":job["input_refs"],"output_refs":[],"started_at":None,"ended_at":now(),"failure_reason":f"EXTERNAL_RUNTIME_REQUIRED:{runtime}"}, actor, status="BLOCKED", parent_id=job_id, input_refs=[job["content_hash"]])

def change_impact(store: AegisStore, changeset_id: str, actor: str) -> dict[str, Any]:
    change=store.get(changeset_id)
    if not change or change["type"] != "ChangeSet": raise ValueError("ChangeSet not found")
    roots=set(change.get("changed_refs", [])); entities=store.list(); affected=[]; frontier=set(roots)
    while frontier:
        ref=frontier.pop()
        for entity in entities:
            if entity["id"] in roots or entity["id"] in {x["id"] for x in affected}: continue
            if ref in entity.get("provenance",{}).get("input_refs",[]) or ref == entity.get("parent_id"):
                affected.append(entity); frontier.add(entity["id"]); frontier.add(entity["content_hash"])
    groups={"requirements":[],"cad":[],"cae":[],"tests":[],"compliance":[],"release_gates":[]}
    mapping={"Requirement":"requirements","CADArtifact":"cad","STEPArtifact":"cad","Drawing":"cad","Mesh":"cae","FEModel":"cae","SolverRun":"cae","CAEResult":"cae","TestProtocol":"tests","TestRun":"tests","TestMeasurement":"tests","RegulatoryRequirement":"compliance","RegulatoryEvidence":"compliance","ReleaseGate":"release_gates"}
    for e in affected:
        if e["type"] in mapping: groups[mapping[e["type"]]].append({"id":e["id"],"hash":e["content_hash"],"type":e["type"]})
    status="UNKNOWN" if not roots else "PASS"
    return store.create("ChangeSet", {"changed_refs":sorted(roots),"affected":groups,"analysis_status":status}, actor, status=status, parent_id=changeset_id, input_refs=sorted(roots))

def compile_protocol(store: AegisStore, requirement_id: str, actor: str) -> dict[str, Any]:
    req=store.get(requirement_id)
    if not req or req["type"] != "Requirement": raise ValueError("Requirement not found")
    acceptance=req.get("acceptance_criteria")
    if not isinstance(acceptance, list) or not acceptance: raise ValueError("acceptance criteria required")
    return store.create("TestProtocol", {"requirement_id":requirement_id,"acceptance_criteria":acceptance,"machine_readable_protocol":{"steps":[{"id":f"AC-{i+1}","criterion":c,"evidence_required":"hashed raw measurement + source + timestamp"} for i,c in enumerate(acceptance)]},"physical_test_status":"PHYSICAL_EVIDENCE_REQUIRED"}, actor, status="UNKNOWN", input_refs=[req["content_hash"]])

def intake_physical(store: AegisStore, payload: dict[str, Any], actor: str) -> dict[str, Any]:
    required={"test_id","source","captured_at","measurement_refs","filename","content_hex"}
    if not required <= payload.keys() or not payload["measurement_refs"]: raise ValueError("PHYSICAL_EVIDENCE_REQUIRED: missing identity, source, timestamp, or measurements")
    item=validate_upload(payload["filename"], bytes.fromhex(payload["content_hex"]), payload.get("media_type","application/octet-stream"))
    return store.create("TestMeasurement", {"test_id":payload["test_id"],"source":payload["source"],"captured_at":payload["captured_at"],"measurement_refs":payload["measurement_refs"],"file":asdict(item),"evidence_class":"MEASURED_PHYSICALLY"}, actor, status="UNKNOWN", input_refs=payload["measurement_refs"])

def correlate(store: AegisStore, simulation_id: str, physical_id: str, actor: str) -> dict[str, Any]:
    simulation, physical=store.get(simulation_id),store.get(physical_id)
    if not simulation or simulation["type"] not in {"CAEResult","SolverRun"}: raise ValueError("simulation evidence required")
    if not physical or physical["type"] != "TestMeasurement" or physical.get("evidence_class") != "MEASURED_PHYSICALLY": raise ValueError("PHYSICAL_EVIDENCE_REQUIRED")
    # Actual numeric correlation requires matching released channel data; absence stays UNKNOWN.
    return store.create("CAEResult", {"kind":"CORRELATION","simulation_id":simulation_id,"physical_id":physical_id,"correlation_status":"UNKNOWN","reason":"matching calibrated numeric channels not supplied"}, actor, status="UNKNOWN", input_refs=[simulation["content_hash"],physical["content_hash"]])

def review(store: AegisStore, provider: str, artifact_refs: list[str], prompt: str, actor: str) -> dict[str, Any]:
    config={"GPT":("AEGIS_GPT_URL","AEGIS_GPT_API_KEY","AEGIS_GPT_MODEL"),"DEEPSEEK":("AEGIS_DEEPSEEK_URL","AEGIS_DEEPSEEK_API_KEY","AEGIS_DEEPSEEK_MODEL")}
    if provider not in config: raise ValueError("unsupported provider")
    url,key,model=(os.getenv(k,"") for k in config[provider]); request_id=str(uuid.uuid4()); context_hash=hashlib.sha256(canonical_json({"prompt":prompt,"artifact_refs":artifact_refs})).hexdigest()
    base={"provider":provider,"model":model or None,"request_id":request_id,"prompt_context_hash":context_hash,"artifact_refs":artifact_refs,"timestamp":now(),"evidence_refs":artifact_refs,"contradictions":[]}
    if not (url and key and model): return store.create("AdversarialReview", {**base,"output":None,"verdict":"EXTERNAL_RUNTIME_REQUIRED","confidence":None}, actor, status="EXTERNAL_RUNTIME_REQUIRED", input_refs=artifact_refs)
    request=urllib.request.Request(url, data=canonical_json({"model":model,"input":prompt,"artifact_refs":artifact_refs}), headers={"Authorization":f"Bearer {key}","Content-Type":"application/json","X-Aegis-Request-Id":request_id}, method="POST")
    try:
        with urllib.request.urlopen(request, timeout=30) as response: output=json.loads(response.read())
    except Exception as exc:
        return store.create("AdversarialReview", {**base,"output":None,"verdict":"EXTERNAL_RUNTIME_REQUIRED","confidence":None,"error":type(exc).__name__}, actor, status="EXTERNAL_RUNTIME_REQUIRED", input_refs=artifact_refs)
    verdict=output.get("verdict") if isinstance(output,dict) else None
    if verdict not in {"PASS","BLOCK","UNKNOWN"}: verdict="UNKNOWN"
    return store.create("AdversarialReview", {**base,"output":output,"verdict":verdict,"confidence":output.get("confidence") if isinstance(output,dict) else None}, actor, status=verdict, input_refs=artifact_refs)

def contradiction(store: AegisStore, first_id: str, second_id: str, actor: str) -> dict[str, Any]:
    a,b=store.get(first_id),store.get(second_id)
    if not a or not b or a["type"] != "AdversarialReview" or b["type"] != "AdversarialReview": raise ValueError("review IDs required")
    av,bv=a.get("verdict"),b.get("verdict")
    if "EXTERNAL_RUNTIME_REQUIRED" in {av,bv} or "UNKNOWN" in {av,bv,None}: classification="INSUFFICIENT_EVIDENCE"
    elif av==bv: classification="AGREE"
    elif {av,bv} <= {"PASS","BLOCK"}: classification="DISAGREE"
    else: classification="UNRESOLVED"
    return store.create("Contradiction", {"review_ids":[first_id,second_id],"artifact_refs":sorted(set(a.get("artifact_refs",[])+b.get("artifact_refs",[]))),"classification":classification,"evidence_refs":[a["content_hash"],b["content_hash"]],"unresolved":classification in {"DISAGREE","INSUFFICIENT_EVIDENCE","UNRESOLVED"}}, actor, status="BLOCK" if classification != "AGREE" else "PASS", input_refs=[a["content_hash"],b["content_hash"]])
