'''Module 4 of Finding hidden messages in DNA by University of San Diego.'''
from random import randint, uniform

nucleotides = {'A', 'T', 'C', 'G'}

# From module 3
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


# New code

def randomized_motif_search(dna_strings:list[str], k: int, t: int) -> list[str]:
    """
    RandomizedMotifSearch
    Input: Integers k and t, followed by a space-separated collection of strings Dna.
    Output: A collection BestMotifs resulting from running RandomizedMotifSearch(Dna, k, t) 1,000 times. Remember to use pseudocounts!
    """
    number_dna_strings = len(dna_strings)
    len_dna_string = len(dna_strings[0])
    assert t == number_dna_strings
    assert t > 1

    best_motifs = []

    for _ in range(1000):
        # Init of motifs as randomly selected k-mers of every dna string
        local_motifs = []
        best_local_motifs = []

        for i in range(number_dna_strings):
            random_index = randint(0, len_dna_string - k)
            best_local_motifs.append(dna_strings[i][random_index: random_index + k])

        while(True):
            local_motifs.clear()
            profile_matrix = compute_profile_matrix_with_pseudocounts(best_local_motifs)
            for d in dna_strings:
                d_kmer = profile_most_probable_kmer(d, k, profile_matrix)
                local_motifs.append(d_kmer)
            if score_motifs(local_motifs) < score_motifs(best_local_motifs):
                best_local_motifs = local_motifs.copy()
            else:
                break
        
        if len(best_motifs) == 0 or score_motifs(best_local_motifs) < score_motifs(best_motifs):
            # if len(best_motifs) > 0:
                # print("Swap: former best motifs", best_motifs, ", new best motifs: ", best_local_motifs, " ", score_motifs(best_local_motifs), " < ", score_motifs(best_motifs))
            best_motifs = best_local_motifs.copy()
            print("Update: ", best_motifs, "new score: ", score_motifs(best_local_motifs))
    
    return best_motifs

def pick_random_candidate(dna_string:str, k:int, probabilities_matrix:list[list[float]]) -> str:
    """
    Compute probabilities for k-mers in a string and pick randomly with weighed probabilities
    
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
    probabilities = []

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
        probabilities.append(local_score)
        score += local_score
    
    chosen_value = uniform(0.0, score)
    # print("dna_string: ", dna_string, ", probabilities: ", probabilities, ", total_score: ", score)

    previous_cumulated_value = 0
    current_value = 0
    for i in range(len(probabilities)):
        current_value += probabilities[i]
        if chosen_value > previous_cumulated_value and chosen_value < current_value:
            result = dna_string[i: i + k]
            # print("Value: ", chosen_value, "chosen kmer: ", result, "at index: ", i)
            break
        previous_cumulated_value = current_value
    
    return result

def gibbs_sampler(dna_strings:list[str], k: int, t:int, N:int) -> list[str]:
    """
    Input: Integers k, t, and N, followed by a space-separated collection of strings Dna.
    Output: The strings BestMotifs resulting from running GibbsSampler(Dna, k, t, N) with 20 random starts. Remember to use pseudocounts!
    """

    number_dna_strings = len(dna_strings)
    len_dna_string = len(dna_strings[0])
    assert t == number_dna_strings
    assert t > 1

    best_motifs = []

    for _ in range(20):
        # Init of motifs as randomly selected k-mers of every dna string
        local_motifs = []
        best_local_motifs = []

        for i in range(number_dna_strings):
            random_index = randint(0, len_dna_string - k)
            best_local_motifs.append(dna_strings[i][random_index: random_index + k])

        local_motifs = best_local_motifs.copy()

        for i in range(N):
            # print("local: ", local_motifs)
            removed_index = randint(0, t - 1)
            # print("removed index: ", removed_index)
            local_motifs.pop(removed_index)
            # Compute the profile for the best local motifs with one line removed
            profile_matrix = compute_profile_matrix_with_pseudocounts(local_motifs)
            # Choose a new candidate k-mer for line "removed_index"
            kmer = pick_random_candidate(dna_strings[removed_index], k, profile_matrix)
            local_motifs.insert(removed_index, kmer)
            # print("local: ", local_motifs)

            if score_motifs(local_motifs) < score_motifs(best_local_motifs):
                best_local_motifs = local_motifs.copy()
        
        if len(best_motifs) == 0 or score_motifs(best_local_motifs) < score_motifs(best_motifs):
            best_motifs = best_local_motifs.copy()
            print("Update: ", best_motifs, "new score: ", score_motifs(best_local_motifs))

    return best_motifs 