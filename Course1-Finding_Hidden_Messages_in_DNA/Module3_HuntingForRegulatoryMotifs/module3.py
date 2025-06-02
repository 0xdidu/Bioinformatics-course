'''Module 3 of Finding hidden messages in DNA by University of San Diego.'''

nucleotides = {'A', 'T', 'C', 'G'}

# From module 2

def hamming_distance(str1: str, str2: str) -> int:
    # Optimization: return sum(x != y for x, y in zip(str1, str2, strict=True))
    if len(str1) != len(str2):
        raise ValueError("Mismatch in string lengths")
    result = 0
    for x, y in zip(str1, str2):
        if x != y:
            result += 1
    return result

def list_pattern_neighbors(pattern:str, distance:int) -> list:
    '''Find the d-neighborhood of a string 
    (patterns that differ from d characters from the original pattern).
    Input: A string Pattern and an integer d.
    Output: The collection of strings Neighbors(Pattern, d).
    '''
    if distance == 0:
        return [pattern]

    if len(pattern) == 1:
        return ['A', 'T', 'C', 'G']
    
    neighbors = []

    suffix_neighbors = list_pattern_neighbors(pattern[1:], distance)
    for sn in suffix_neighbors:
        if hamming_distance(sn, pattern[1:]) < distance:
            for nucl in nucleotides:
                neighbors.append(nucl + sn)
        else:
            neighbors.append(pattern[0] + sn)
    
    return neighbors

def motif_enumeration_brute_force(dna_strings:list[str], k:int, distance:int) -> list[str]:
    """ Naive implementation to find motives in a given set of DNA strings. 
        A (k,d)-motif appears in every string from Dna with at most d mismatches.
        Input: A collection of strings Dna, and integers k and d.
        Output: All (k, d)-motifs in Dna.
    """
    motives = []
    neighbors = []

    if len(dna_strings) < 2:
        raise ValueError("There should be at least 2 DNA strings")

    # Find all neighbors of potential k-mers on first line of DNA
    first_dna_string = dna_strings[0]
    for i in range(len(first_dna_string) - k + 1):
        word = first_dna_string[i : i + k]
        neighbors += list_pattern_neighbors(word, distance)
    neighbors = list(dict.fromkeys(neighbors))

    # Check if those candidates are at hamming distance "distance" of at least one word per dna string, for all the strings
    for n in neighbors:
        should_add = 1
        for d in dna_strings[1:]:
            for i in range(len(d) - k + 1):
                if hamming_distance(d[i:i+k], n) <= distance:
                    should_add += 1
                    break
        if should_add == len(dna_strings):
            motives.append(n)

    return motives

def generate_all_patterns(k:int) -> list[str]:
    return list_pattern_neighbors("A"*k, k)

def distance_between_pattern_and_strings(dna_strings:list[str], pattern:str) -> int:
    distance = 0
    k = len(pattern)
    for d in dna_strings:
        local_distance = k
        for j in range(len(d) - k + 1):
            # Find the local minimal distance for pattern p
            hd = hamming_distance(pattern, d[j : j + k])
            local_distance = min(hd, local_distance)
        distance += local_distance
    return distance

def median_string_brute_force(dna_strings:list[str], k:int) -> str:
    """ Median String Problem: Find a median string that minimizes the hamming distance for one motif in each dna string.
        Input: A collection of strings Dna and an integer k.
        Output: A k-mer Pattern that minimizes d(Pattern, Dna) among all possible choices of k-mers.
    """
    # Arbitrary initialization 
    minimal_distance = k * len(dna_strings)
    pattern = "A" * k

    # Generate all candidates
    possible_patterns = generate_all_patterns(k)

    # Loop through the candidates, compute the sum of the hamming distances with all 
    for p in possible_patterns:
        distance = distance_between_pattern_and_strings(dna_strings, p)

        if minimal_distance >= distance: # Could be only > but it looks like there is a bug in the online solver
            minimal_distance = distance
            print("distance: ", distance, "pattern: ", pattern)
            pattern = p
            # print("minimal distance: ", distance, "pattern: ", pattern)

    return pattern

def profile_most_probable_kmer(dna_string:str, k:int, probabilities_matrix:list[list[float]]) -> str:
    """
    Profile-most Probable k-mer Problem: Find a Profile-most probable k-mer in a string.

    Input: A string Text, an integer k, and a 4 × k matrix Profile.
    Output: A Profile-most probable k-mer in Text.
    If there are multiple Profile-most probable k-mers in Text, then we select the first such k-mer occurring in Text
    Arbitrary order: A // C // G // T
    """

    if len(dna_string) < k:
        raise ValueError("Dna string shorter than k.")
    
    # 2 options: pre-score the characters of the string or compute each possible k-mer in a sliding window
    result = dna_string[0:k]
    score = 0

    for i in range(len(dna_string) - k + 1):
        candidate = dna_string[i : i + k]
        local_score = 1
        for j in range(k):
            if candidate[j] == "A":
                local_score *= probabilities_matrix[0][j]
            if candidate[j] == "C":
                local_score *= probabilities_matrix[1][j]
            if candidate[j] == "G":
                local_score *= probabilities_matrix[2][j]
            if candidate[j] == "T":
                local_score *= probabilities_matrix[3][j]
        if local_score > score:
            result = candidate
            score = local_score
        
    return result

def compute_profile_matrix(motifs:list[str])->list[list[float]]:
    # Arbitrary order: A // C // G // T
    k = len(motifs[0])
    n = len(motifs)

    result = [[0.0] * k for _ in range(4)]

    for c in range(k):
        for m in motifs:
            if m[c] == "A":
                result[0][c] += 1/n
            if m[c] == "C":
                result[1][c] += 1/n
            if m[c] == "G":
                result[2][c] += 1/n
            if m[c] == "T":
                result[3][c] += 1/n
    return result

def score_motifs(motifs:list[str])->int:
    """
    Calculate the score of given motifs, which reflects the difference between the most common nucleotide 
    and the rest in each column of the motifs.

    Parameters:
    motifs (list): A list of strings, where each string represents a motif.

    Returns:
    int: The calculated score.
    """
    # Initialize the score variable to 0
    score = 0

    # Initialize variables
    k = len(motifs[0])
    number_of_motifs = len(motifs)
    
    # Iterate over the columns
    for i in range(k):   
        # Initialize a dictionary to count nucleotides
        nuc = {"A": 0, "C": 0, "G": 0, "T": 0}
        max_key = 'A'
        # Iterate over the lines
        for line in range(number_of_motifs):
            c = motifs[line][i]
            nuc[c] += 1
        
        # Find the nucleotide with the maximum count
        for key in nuc.keys():
            if nuc[key] == max(nuc.values()):
                max_key = key
        
        # Calculate the score by summing up the counts of nucleotides that are not the most frequent
        for key in nuc.keys():
            if key != max_key:
                score += nuc[key]
    
    # Return the calculated score
    return score


def greedy_motif_search(dna_strings:list[str], k:int, t:int)->list[str]:
    """
    Input: Integers k and t, followed by a space-separated collection of strings Dna.
    Output: A collection of strings BestMotifs resulting from applying GreedyMotifSearch(Dna, k, t). If at any step you find more than one Profile-most probable k-mer in a given string, use the one occurring first.
    """
    number_dna_strings = len(dna_strings)
    len_dna_string = len(dna_strings[0])
    assert t == number_dna_strings
    assert t > 1

    # Init of best motifs as the first k-mers of every dna string
    best_motifs = []
    for i in range(number_dna_strings):
        best_motifs.append(dna_strings[i][0:k])

    # Init profile matrix
    # profile_matrix = [[0.0] * k for _ in range(4)]

    for i in range(len_dna_string - k + 1):
        local_motifs = []
        local_motifs.append(dna_strings[0][i : i + k])
        # print("Iteration DNA String 0 #i: ", i, ", local_motifs: ", local_motifs)
        # Iteration on each DNA string
        for j in range(1, t):
            d = dna_strings[j]
            # Iteration on each k-mer and evaluation of the score
            profile_matrix = compute_profile_matrix(local_motifs)
            d_kmer = profile_most_probable_kmer(d, k, profile_matrix)
            local_motifs.append(d_kmer)
            # print("Iteration DNA Strings #j: ", j, ", local_motifs: ", local_motifs, ", profile matrix: ", profile_matrix)

        # print("Score: local ", score_motifs(local_motifs), " best ", score_motifs(best_motifs))
        if score_motifs(local_motifs) < score_motifs(best_motifs):
            best_motifs = local_motifs

    return best_motifs

def compute_profile_matrix_with_pseudocounts(motifs:list[str])->list[list[float]]:
    # Arbitrary order: A // C // G // T
    k = len(motifs[0])
    n = len(motifs)

    result = [[1/(2*n)] * k for _ in range(4)]

    for c in range(k):
        for m in motifs:
            if m[c] == "A":
                result[0][c] += 1/(2*n)
            if m[c] == "C":
                result[1][c] += 1/(2*n)
            if m[c] == "G":
                result[2][c] += 1/(2*n)
            if m[c] == "T":
                result[3][c] += 1/(2*n)
        
    return result

def greedy_motif_search_with_pseudocounts(dna_strings:list[str], k:int, t:int)->list[str]:
    """
    Input: Integers k and t, followed by a space-separated collection of strings Dna.
    Output: A collection of strings BestMotifs resulting from applying GreedyMotifSearch(Dna, k, t). If at any step you find more than one Profile-most probable k-mer in a given string, use the one occurring first.
    """
    number_dna_strings = len(dna_strings)
    len_dna_string = len(dna_strings[0])
    assert t == number_dna_strings
    assert t > 1

    # Init of best motifs as the first k-mers of every dna string
    best_motifs = []
    for i in range(number_dna_strings):
        best_motifs.append(dna_strings[i][0:k])

    # Init profile matrix
    # profile_matrix = [[0.0] * k for _ in range(4)]

    for i in range(len_dna_string - k + 1):
        local_motifs = []
        local_motifs.append(dna_strings[0][i : i + k])
        # print("Iteration DNA String 0 #i: ", i, ", local_motifs: ", local_motifs)
        # Iteration on each DNA string
        for j in range(1, t):
            d = dna_strings[j]
            # Iteration on each k-mer and evaluation of the score
            profile_matrix = compute_profile_matrix_with_pseudocounts(local_motifs)
            d_kmer = profile_most_probable_kmer(d, k, profile_matrix)
            local_motifs.append(d_kmer)
            # print("Iteration DNA Strings #j: ", j, ", local_motifs: ", local_motifs, ", profile matrix: ", profile_matrix)

        # print("Score: local ", score_motifs(local_motifs), " best ", score_motifs(best_motifs))
        if score_motifs(local_motifs) < score_motifs(best_motifs):
            best_motifs = local_motifs

    return best_motifs