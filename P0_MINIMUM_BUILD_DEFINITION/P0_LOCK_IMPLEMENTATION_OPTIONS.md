[
  {
    "OPTION": "Mechanical pawl + capture tooth",
    "SOURCE_BASIS": "V6 lock board / V7 multi-state lock architecture",
    "FUNCTION": "Positive constraint engagement/disengagement",
    "ADVANTAGE": "Directly testable mechanical state",
    "RISK": "Engagement tolerance, wear, jamming",
    "REQUIRED_TEST": "engagement, incomplete engagement, delayed engagement, reverse-load hold",
    "AUTHORITY_REQUIRED": "Yes — final geometry/loads"
  },
  {
    "OPTION": "Lock pin + capture pocket",
    "SOURCE_BASIS": "V6 lock component concept includes lock pin",
    "FUNCTION": "Positive capture of defined mechanism state",
    "ADVANTAGE": "Simple discrete state verification",
    "RISK": "Shear/load concentration; alignment sensitivity",
    "REQUIRED_TEST": "alignment, engagement, reverse-load hold, release",
    "AUTHORITY_REQUIRED": "Yes"
  },
  {
    "OPTION": "Inertia-driven pawl/cam",
    "SOURCE_BASIS": "V6 inertial mass + drive cam + pawl concept",
    "FUNCTION": "Passive event-driven transition",
    "ADVANTAGE": "No external power required in concept",
    "RISK": "Threshold, latency, sensitivity to orientation",
    "REQUIRED_TEST": "trigger latency, false trigger, missed trigger, fail-safe",
    "AUTHORITY_REQUIRED": "Yes"
  }
]