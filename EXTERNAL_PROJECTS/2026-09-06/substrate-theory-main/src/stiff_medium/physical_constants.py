"""Shared physical constants used by the stiff-medium audit code."""

from __future__ import annotations

from typing import Final

#: NIST/CODATA 2022 inverse fine-structure constant at zero momentum transfer.
#: Source: https://physics.nist.gov/cgi-bin/cuu/Value?eqalphinv=
INV_ALPHA_THOMSON_CODATA_2022: Final[float] = 137.035999177

#: Standard uncertainty in the CODATA 2022 inverse fine-structure constant.
INV_ALPHA_THOMSON_CODATA_2022_UNCERTAINTY: Final[float] = 0.000000021

#: Fine-structure constant at zero momentum transfer.
ALPHA_THOMSON_CODATA_2022: Final[float] = 1.0 / INV_ALPHA_THOMSON_CODATA_2022
