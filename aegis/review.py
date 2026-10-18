"""Provider-neutral adversarial review contracts; LLMs never authorize release."""
from __future__ import annotations
from dataclasses import dataclass
from .kernel import ReviewFinding

@dataclass(frozen=True)
class ReviewRequest:
    claim: str
    evidence_ids: tuple[str, ...]
    prompt_digest: str

class AdversarialReviewer:
    provider: str
    def review(self, request: ReviewRequest) -> tuple[ReviewFinding, ...]:
        raise NotImplementedError("configure a GPT, DeepSeek, or independent audit provider")
