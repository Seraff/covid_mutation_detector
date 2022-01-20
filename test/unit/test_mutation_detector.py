#!/usr/bin/env python3

import os
import sys
from pathlib import Path
import unittest

ROOT_PATH = str(Path(os.path.dirname(
    os.path.realpath(__file__))).parent.parent)
sys.path.insert(1, ROOT_PATH)

from lib.mutation_detector import MutationDetector

class TestMutationDetector(unittest.TestCase):
    def test_shit(self):
        """
        Test shit
        """
        ref_seq   = 'NANADAGAGADUPA'
        query_seq = 'NANA--VAGA--PA'

        # self.assertEqual(6, 6)


if __name__ == '__main__':
    unittest.main()
