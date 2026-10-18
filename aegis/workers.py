"""Typed worker fabric. Workers must emit artifacts; unavailable workers block."""
from __future__ import annotations
from dataclasses import dataclass
from .kernel import Artifact, Decision, Phase

@dataclass(frozen=True)
class WorkerResult:
    worker: str
    phase: Phase
    decision: Decision
    reason: str
    artifacts: tuple[Artifact, ...] = ()

class Worker:
    name: str
    phase: Phase
    def run(self) -> WorkerResult:
        return WorkerResult(self.name, self.phase, Decision.BLOCK, "worker adapter not configured")

class CADWorker(Worker): name, phase = "cad", Phase.CAD
class CAEWorker(Worker): name, phase = "cae", Phase.CAE
class DFMWorker(Worker): name, phase = "dfm", Phase.DFM
class TestWorker(Worker): name, phase = "test", Phase.TEST
class ComplianceWorker(Worker): name, phase = "compliance", Phase.COMPLIANCE
class ManufacturingWorker(Worker): name, phase = "manufacturing", Phase.MANUFACTURING
