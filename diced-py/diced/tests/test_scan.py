import unittest

from ..lib import scan


SEQUENCE = "".join(
    [
        "TTTTACAATCTGCGTTTTAACTCCACACGGTACATTAGAAACCATCTGCAACATATT",
        "CAAGTTCAGCTTCAAAACCTTGTTTTAACTCCACACGGTACATTAGAAACTTCGTCA",
        "AGCTTTACCTCAAAAGTCCTCTCAAACCTGTTTTAACTCCACACGGTACATTAGAAA",
        "CAATAATCAACAACTCTTTGATTTTGTGAAATGGAAGAAGTTTTAACTCCACACGGT",
        "ACATTAGAAACAGAACTCTCAGAAGAACCGAGAGCTTTTTCTATTAACGTTTTAACT",
        "CCACACGGTACATTAGAAACCCTGCGTGCCTGTGTCTAAAAAATA",
    ]
)


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

    def test_parameters(self):
        # 5 repeats of 29 bp separated by spacers of 32 to 40 bp
        self.assertEqual(len(list(scan(SEQUENCE, min_repeat_count=5))), 1)
        self.assertEqual(len(list(scan(SEQUENCE, min_repeat_count=6))), 0)
        self.assertEqual(len(list(scan(SEQUENCE, max_repeat_length=25))), 0)
        self.assertEqual(len(list(scan(SEQUENCE, min_spacer_length=45))), 0)
        crisprs = list(scan(SEQUENCE, min_repeat_length=19, max_repeat_length=38, min_spacer_length=19, max_spacer_length=48))
        self.assertEqual([(c.start, c.end) for c in crisprs], [(13, 305)])
