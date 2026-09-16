from __future__ import annotations
from enum import StrEnum
class Role(StrEnum):
 ARCHITECT="ARCHITECT"; CAD_ENGINEER="CAD_ENGINEER"; CAE_ENGINEER="CAE_ENGINEER"; TEST_ENGINEER="TEST_ENGINEER"; MANUFACTURING_ENGINEER="MANUFACTURING_ENGINEER"; COMPLIANCE_REVIEWER="COMPLIANCE_REVIEWER"; INDEPENDENT_AUDITOR="INDEPENDENT_AUDITOR"
PERMISSIONS={
 Role.ARCHITECT:{"Requirement","ArchitectureDecision","WorkerJob"}, Role.CAD_ENGINEER:{"CADArtifact","STEPArtifact","Drawing","BOM"}, Role.CAE_ENGINEER:{"Mesh","FEModel","SolverRun","CAEResult"}, Role.TEST_ENGINEER:{"TestProtocol","TestRun","TestMeasurement"}, Role.MANUFACTURING_ENGINEER:{"SupplierEvidence","OEMEvidence","BOM"}, Role.COMPLIANCE_REVIEWER:{"RegulatoryRequirement","RegulatoryEvidence","AdversarialReview"}, Role.INDEPENDENT_AUDITOR:{"Audit","Contradiction","ReleaseGate"}}
def allowed(role: Role, entity_type: str) -> bool: return entity_type in PERMISSIONS[role]
