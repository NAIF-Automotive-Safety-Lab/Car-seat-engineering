from __future__ import annotations

import sys
import unittest
from pathlib import Path

sys.path.insert(0, str(Path(__file__).resolve().parents[1]))

from bridge.r41_paths import R41_STEP  # noqa: E402
from model.r41_model import import_step_with_ocp  # noqa: E402


class R41StepBrepTests(unittest.TestCase):
    def test_step_read_and_62_solids(self) -> None:
        shape = import_step_with_ocp()
        from OCP.TopAbs import TopAbs_SOLID
        from OCP.TopExp import TopExp_Explorer

        explorer = TopExp_Explorer(shape, TopAbs_SOLID)
        count = 0
        while explorer.More():
            count += 1
            explorer.Next()
        self.assertEqual(count, 62)
        self.assertTrue(R41_STEP.is_file())


if __name__ == "__main__":
    unittest.main()
