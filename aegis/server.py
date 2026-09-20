from __future__ import annotations
import json, os
from http import HTTPStatus
from http.server import BaseHTTPRequestHandler, ThreadingHTTPServer
from pathlib import Path
from urllib.parse import parse_qs, urlparse
from .core import AEOK, AegisStore, Lifecycle
from .federation import verify_handshake
from .ingestion import validate_upload
from .workflows import WorkerFabric, change_impact, compile_protocol, intake_physical, correlate, review, contradiction

class App:
    def __init__(self, db_path: str): self.store=AegisStore(db_path); self.aeok=AEOK(self.store)
    def authorize(self, headers) -> str:
        token=os.getenv("AEGIS_AUTH_TOKEN")
        if not token: raise PermissionError("AUTH_CONFIGURATION_REQUIRED")
        if headers.get("Authorization") != f"Bearer {token}": raise PermissionError("UNAUTHORIZED")
        return headers.get("X-Aegis-Actor", "unknown")

def handler(app: App):
 class Handler(BaseHTTPRequestHandler):
    def log_message(self, format, *args): pass
    def send(self, status, body, content_type="application/json"):
        raw = json.dumps(body, sort_keys=True).encode() if content_type == "application/json" else body
        self.send_response(status); self.send_header("Content-Type", content_type); self.send_header("Content-Length", str(len(raw))); self.end_headers(); self.wfile.write(raw)
    def json(self):
        length=int(self.headers.get("Content-Length", "0"));
        if length > 26*1024*1024: raise ValueError("request too large")
        return json.loads(self.rfile.read(length))
    def do_GET(self):
        path=urlparse(self.path).path
        if path == "/health": return self.send(200, {"status":"UNKNOWN", "service":"V7-AEGIS", "auth":"REQUIRES_CONFIGURATION" if not os.getenv("AEGIS_AUTH_TOKEN") else "CONFIGURED"})
        if path == "/":
            return self.send(200, (Path(__file__).parent/"static"/"index.html").read_bytes(), "text/html; charset=utf-8")
        try:
            self.actor=app.authorize(self.headers)
            if path == "/api/entities": return self.send(200, app.store.list(parse_qs(urlparse(self.path).query).get("type", [None])[0]))
            if path.startswith("/api/entities/"):
                item=app.store.get(path.rsplit("/",1)[1]); return self.send(200, item) if item else self.send(404,{"error":"NOT_FOUND"})
            self.send(404,{"error":"NOT_FOUND"})
        except PermissionError as e: self.send(503 if str(e)=="AUTH_CONFIGURATION_REQUIRED" else 401,{"error":str(e)})
    def do_POST(self):
        try:
            actor=app.authorize(self.headers); path=urlparse(self.path).path; body=self.json()
            if path == "/api/entities":
                typ=body["type"]; content=body.get("content", {}); status=body.get("status", "DRAFT" if typ=="Requirement" else "UNKNOWN")
                return self.send(201, app.store.create(typ, content, actor, status=status, input_refs=body.get("input_refs", [])))
            if path.startswith("/api/transition/"):
                return self.send(201, app.aeok.transition(path.rsplit("/",1)[1], Lifecycle(body["target"]), actor))
            if path == "/api/release-gates": return self.send(201, app.aeok.release_gate(body["candidate_id"], body["checks"], actor))
            if path == "/api/jobs": return self.send(201, WorkerFabric(app.store).create(body["kind"], body.get("input_refs", []), actor))
            if path.startswith("/api/jobs/") and path.endswith("/external-runtime-required"): return self.send(201, WorkerFabric(app.store).external_runtime_required(path.split("/")[3], body["runtime"], actor))
            if path == "/api/change-impact": return self.send(201, change_impact(app.store, body["changeset_id"], actor))
            if path == "/api/test-protocols": return self.send(201, compile_protocol(app.store, body["requirement_id"], actor))
            if path == "/api/physical-evidence": return self.send(201, intake_physical(app.store, body, actor))
            if path == "/api/correlation": return self.send(201, correlate(app.store, body["simulation_id"], body["physical_id"], actor))
            if path == "/api/adversarial-reviews": return self.send(201, review(app.store, body["provider"], body.get("artifact_refs", []), body["prompt"], actor))
            if path == "/api/contradictions": return self.send(201, contradiction(app.store, body["first_review_id"], body["second_review_id"], actor))
            if path == "/api/upload":
                data=bytes.fromhex(body["content_hex"]); record=validate_upload(body["filename"], data, body.get("media_type", "application/octet-stream"))
                return self.send(201, app.store.create(body.get("entity_type", "CADArtifact"), record.__dict__, actor, status="UNKNOWN"))
            if path == "/api/federation/handshake":
                result=verify_handshake(body, os.getenv("AEGIS_FEDERATION_SHARED_SECRET", ""), app.store)
                return self.send(201, app.store.create("FederationHandshake", {"envelope": {k:v for k,v in body.items() if k != "signature"}, "verification":result}, actor, status=result))
            self.send(404,{"error":"NOT_FOUND"})
        except PermissionError as e: self.send(503 if str(e)=="AUTH_CONFIGURATION_REQUIRED" else 401,{"error":str(e)})
        except (KeyError, ValueError, json.JSONDecodeError) as e: self.send(400,{"error":"INVALID_REQUEST", "detail":str(e)})
        except Exception: self.send(500,{"error":"INTERNAL_ERROR"})
 return Handler

def serve(host="127.0.0.1", port=8080, db_path="aegis.db"):
    server=ThreadingHTTPServer((host,port), handler(App(db_path))); server.serve_forever()
if __name__ == "__main__": serve(host=os.getenv("AEGIS_HOST", "127.0.0.1"), port=int(os.getenv("PORT", "8080")), db_path=os.getenv("AEGIS_DB", "aegis.db"))
