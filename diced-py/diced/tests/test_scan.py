import unittest

from ..lib import scan


class TestScan(unittest.TestCase):
    def test_empty(self):
        crisprs = list(scan(""))
        self.assertEqual(len(crisprs), 0)

    def test_default(self):
        seq = "".join(
            [
                "TTTTACAATCTGCGTTTTAACTCCACACGGTACATTAGAAACCATCTGCAACATATT",
                "CAAGTTCAGCTTCAAAACCTTGTTTTAACTCCACACGGTACATTAGAAACTTCGTCA",
                "AGCTTTACCTCAAAAGTCCTCTCAAACCTGTTTTAACTCCACACGGTACATTAGAAA",
                "CAATAATCAACAACTCTTTGATTTTGTGAAATGGAAGAAGTTTTAACTCCACACGGT",
                "ACATTAGAAACAGAACTCTCAGAAGAACCGAGAGCTTTTTCTATTAACGTTTTAACT",
                "CCACACGGTACATTAGAAACCCTGCGTGCCTGTGTCTAAAAAATA",
            ]
        )
        crisprs = list(scan(seq))
        self.assertEqual(len(crisprs), 1)
        self.assertEqual(crisprs[0].start, 13)
        self.assertEqual(crisprs[0].end, 305)

        crispr = crisprs[0]
        repeats = crispr.repeats
        self.assertEqual(len(repeats), 5)
        self.assertEqual(repeats[0].start, 13)
        self.assertEqual(repeats[0].end, 42)
