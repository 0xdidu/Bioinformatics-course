import module3
import unittest

class TestMotifEnumeration(unittest.TestCase):
    def test_motif_enumeration_brute_force(self):
        self.assertEqual(module3.motif_enumeration_brute_force(["ATTTGGC", "TGCCTTA", "CGGTATC", "GAAAATT"], 3, 1).sort(), ["ATA", "ATT", "GTT", "TTT"].sort())

class TestMedianString(unittest.TestCase):
    def test_generate_all_patterns(self):
        self.assertEqual(module3.generate_all_patterns(2).sort(), ["AA", "AT", "AC", "AG", "TT", "TA", "TC", "TG", "GG", "GC", "GA", "GT", "CC", "CG", "CA", "CT"].sort())

    def test_distance_between_pattern_and_strings(self):
        self.assertEqual(module3.distance_between_pattern_and_strings(["TTACCTTAAC", "GATATCTGTC", "ACGGCGTTCG", "CCCTAAAGAG", "CGTCAGAGGT"], "AAA"), 5)
    
    def test_median_string_brute_force(self):
        # There are 2 solutions: GAC and ACG, equally optimal, keeping the last one
        self.assertEqual(module3.median_string_brute_force(["AAATTGACGCAT", "GACGACCACGTT", "CGTCAGCGCCTG", "GCTGAGCACCGG", "AGTTCGGGACAG"], 3), "GAC")
        self.assertEqual(module3.median_string_brute_force(["AAATTGACGCAT","GACGACCACGTT","CGTCAGCGCCTG","GCTGAGCACCGG","AGTACGGGACAG"], 3), "ACG")
        self.assertEqual(module3.median_string_brute_force(["ACGT","ACGT","ACGT"], 3), "ACG")
        self.assertEqual(module3.median_string_brute_force(["ATA","ACA","AGA","AAT","AAC"], 3), "AAA")
        self.assertEqual(module3.median_string_brute_force(["AAG","AAT"], 3), "AAG")

class TestGreedyMotifSearchHelpers(unittest.TestCase):
    def test_compute_profile_matrix(self):
        self.assertEqual(
            module3.compute_profile_matrix(
                ["ACCT", "AGGA"]
            ), [[1, 0, 0, 0.5], [0, 0.5, 0.5, 0], [0, 0.5, 0.5, 0], [0, 0, 0, 0.5]]
        )
    def test_score_motifs(self):
        self.assertEqual(
            module3.score_motifs([
                "TCGGGGGTTTTT",
                "CCGGTGACTTAC",
                "ACGGGGATTTTC",
                "TTGGGGACTTTT",
                "AAGGGGACTTCC",
                "TTGGGGACTTCC",
                "TCGGGGATTCAT", 
                "TCGGGGATTCCT", 
                "TAGGGGAACTAC", 
                "TCGGGTATAACC"
            ]), 30
        )

class TestGreedyMotifSearch(unittest.TestCase):
    def test_profile_most_probable_kmer(self):
        self.assertEqual(
            module3.profile_most_probable_kmer(
                "ACCTGTTTATTGCCTAAGTTCCGAACAAACCCAATATAGCCCGAGGGCCT",
                5,
                [
                    [0.2, 0.2, 0.3, 0.2, 0.3],
                    [0.4, 0.3, 0.1, 0.5, 0.1],
                    [0.3, 0.3, 0.5, 0.2, 0.4],
                    [0.1, 0.2, 0.1, 0.1, 0.2],
                ]
            ),
            "CCGAG"
        )
    
    def test_greedy_motif_search(self): 
        self.assertEqual(
            module3.greedy_motif_search(
                ["GGCGTTCAGGCA", "AAGAATCAGTCA", "CAAGGAGTTCGC", "CACGTCAATCAC", "CAATAATATTCG"], 3, 5),
            ["CAG", "CAG", "CAA", "CAA", "CAA"]
        )
        self.assertEqual(
            module3.greedy_motif_search(
            ["GCCCAA","GGCCTG","AACCTA","TTCCTT"], 3, 4),
            ["GCC","GCC", "AAC", "TTC"]
        )
        self.assertEqual(
            module3.greedy_motif_search(
                [
                    "GAGGCGCACATCATTATCGATAACGATTCGCCGCATTGCC", "TCATCGAATCCGATAACTGACACCTGCTCTGGCACCGCTC", "TCGGCGGTATAGCCAGAAAGCGTAGTGCCAATAATTTCCT", "GAGTCGTGGTGAAGTGTGGGTTATGGGGAAAGGCAGACTG", 
                    "GACGGCAACTACGGTTACAACGCAGCAACCGAAGAATATT", "TCTGTTGTTGCTAACACCGTTAAAGGCGGCGACGGCAACT", "AAGCGGCCAACGTAGGCGCGGCTTGGCATCTCGGTGTGTG", "AATTGAAAGGCGCATCTTACTCTTTTCGCTTTCAAAAAAA"
                ], 5, 8
            ), 
            ["GAGGC", "TCATC", "TCGGC", "GAGTC", "GCAGC", "GCGGC", "GCGGC", "GCATC"]
        )
        self.assertEqual(
            module3.greedy_motif_search(
                [
                    "GCAGGTTAATACCGCGGATCAGCTGAGAAACCGGAATGTGCGT",
                    "CCTGCATGCCCGGTTTGAGGAACATCAGCGAAGAACTGTGCGT",
                    "GCGCCAGTAACCCGTGCCAGTCAGGTTAATGGCAGTAACATTT",
                    "AACCCGTGCCAGTCAGGTTAATGGCAGTAACATTTATGCCTTC",
                    "ATGCCTTCCGCGCCAATTGTTCGTATCGTCGCCACTTCGAGTG"
                ], 
                6, 5
            ), 
            ["GTGCGT","GTGCGT", "GCGCCA","GTGCCA","GCGCCA"]
        )

        self.assertEqual(
            module3.greedy_motif_search(
                [
                    "GACCTACGGTTACAACGCAGCAACCGAAGAATATTGGCAA",
                    "TCATTATCGATAACGATTCGCCGGAGGCCATTGCCGCACA",
                    "GGAGTCTGGTGAAGTGTGGGTTATGGGGCAGACTGGGAAA",
                    "GAATCCGATAACTGACACCTGCTCTGGCACCGCTCTCATC",
                    "AAGCGCGTAGGCGCGGCTTGGCATCTCGGTGTGTGGCCAA",
                    "AATTGAAAGGCGCATCTTACTCTTTTCGCTTAAAATCAAA",
                    "GGTATAGCCAGAAAGCGTAGTTAATTTCGGCTCCTGCCAA",
                    "TCTGTTGTTGCTAACACCGTTAAAGGCGGCGACGGCAACT"
                ],
                5, 8
            ), 
            ["GCAGC", "TCATT", "GGAGT", "TCATC", "GCATC", "GCATC", "GGTAT", "GCAAC"]
        )

class TestGreedyMotifSearchWithPseudocounts(unittest.TestCase):
    def test_greedy_motif_search_with_pseudocounts(self):
        self.assertEqual(
            module3.greedy_motif_search_with_pseudocounts(
                [
                    "GGCGTTCAGGCA", "AAGAATCAGTCA", "CAAGGAGTTCGC", "CACGTCAATCAC", "CAATAATATTCG"
                ], 3, 5
            ),
            ["TTC", "ATC", "TTC", "ATC", "TTC"]
        )
if __name__ == "__main__":
    unittest.main()