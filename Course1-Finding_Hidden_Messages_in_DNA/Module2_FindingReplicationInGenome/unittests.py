import unittest
import module2

class TestSkew(unittest.TestCase):
    def test_skew(self):
        self.assertEqual(module2.skew('CATGGGCATCGGCCATACGCC', 14), [0, -1, -1, -1, 0, 1, 2, 1, 1, 1, 0, 1, 2, 1, 0, 0])
    
    def test_minimum_skey(self):
        self.assertEqual(module2.minimum_skew("TAAAGACTGCCGAGAGGCCAACACGAGTGCTAGAACGAGGGGCGTAAACGCGGGTCCGAT"), [11, 24])

class TestHammingDistance(unittest.TestCase):
    def test_hamming_distance(self):
        self.assertEqual(module2.hamming_distance('GGGCCGTTGGT', 'GGACCGTTGAC'), 3)

class TestApproximatePatternMatching(unittest.TestCase):
    def test_approximate_pattern_matching(self):
        self.assertEqual(module2.approximate_pattern_matching('ATTCTGGA', 'CGCCCGAATCCAGAACGCATTCCCATATTTCGGGACCACTGGCCTCCACGGTACGGACGTCAATCAAAT', 3), [6, 7, 26, 27])
        self.assertEqual(module2.approximate_pattern_matching_count('AAAAA', 'AACAAGCTGATAAACATTTAAAGAG', 1), 4)
        self.assertEqual(module2.approximate_pattern_matching_count('GAGG', 'TTTAGAGCCTTCAGAGG', 2), 4)

class TestListPatternNeighbors(unittest.TestCase):
    def test_list_pattern_neighbors(self):
        self.assertEqual((module2.list_pattern_neighbors('ACG', 1)).sort(), ['CCG', 'TCG', 'GCG', 'AAG', 'ATG', 'AGG', 'ACA', 'ACC', 'ACT', 'ACG'].sort())

class TestFindKmersWithMismatches(unittest.TestCase):
    def test_find_kmers_with_mismatches(self):
        self.assertEqual(module2.find_kmers_with_mismatches('AACAAGCTGATAAACATTTAAAGAG', 5,1), ['AAAAA'])

class TestFindKmersWithMismatchesAndReverseComplements(unittest.TestCase):
    def test_find_kmers_with_mismatches_and_reverse_complements(self):
        self.assertEqual(module2.find_kmers_with_mismatches_and_reverse_complements('ACGTTGCATGTCGCATGATGCATGAGAGCT', 4, 1).sort(), ['ATGT', 'ACAT'].sort())

if __name__ == "__main__":
    unittest.main()
