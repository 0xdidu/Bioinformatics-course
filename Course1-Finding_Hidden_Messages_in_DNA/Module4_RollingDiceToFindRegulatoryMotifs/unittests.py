import module4
import unittest

class TestRandomizedMotifSearch(unittest.TestCase):
    def test_randomized_motif_search(self):
        self.assertEqual(
            module4.randomized_motif_search(
                ["CGCCCCTCTCGGGGGTGTTCAGTAAACGGCCA", 
                 "GGGCGAGGTATGTGTAAGTGCCAAGGTGCCAG", 
                 "TAGTACCGAGACCGAAAGAAGTATACAGGCGT", 
                 "TAGATCAAGTTTCAGGTGCACGTCGGTGAACC", 
                 "AATCCACCAGCTCCACGTGCAATGTTGGCCTA"], 
                8,
                5
            ), ["TCTCGGGG", "CCAAGGTG", "TACAGGCG", "TTCAGGTG", "TCCACGTG"]
        )

class TestGibbsSampler(unittest.TestCase):
    def test_gibbs_sampler(self):
        self.assertEqual(
            module4.gibbs_sampler(
                ["CGCCCCTCTCGGGGGTGTTCAGTAACCGGCCA",
                 "GGGCGAGGTATGTGTAAGTGCCAAGGTGCCAG",
                 "TAGTACCGAGACCGAAAGAAGTATACAGGCGT",
                 "TAGATCAAGTTTCAGGTGCACGTCGGTGAACC",
                 "AATCCACCAGCTCCACGTGCAATGTTGGCCTA"], 
                8, 5, 100
            ), ["TCTCGGGG", "CCAAGGTG", "TACAGGCG", "TTCAGGTG", "TCCACGTG"]
        )

if __name__ == "__main__":
    unittest.main()