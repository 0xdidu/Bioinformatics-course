from random import randint, uniform

nucleotides = {'A', 'T', 'C', 'G'}

# Module 1

def PatternCount(Text, Pattern):
    count = 0
    lt = len(Text)
    lp = len(Pattern)
    for i in range(lt-lp):
        if Text[i : i + len(Pattern)] == Pattern:
            count = count + 1
    print(count)
    return count

def frequencyTable (Text:str, k:int):
    frequency = dict()
    lenText = len(Text)
    for i in range(len(Text)-k +1):
        pattern = Text[i:i+k]
        if not pattern in frequency:
            frequency[pattern] = 1
        else:
            frequency[pattern] += 1
    # print("Text: ", Text, ", k: ", k, ", result: ", frequency)
    return frequency


def maxFrequencyTuple (frequency: dict):
    max = 0
    result = []
    for k, v in frequency.items():
        if v > max:
            max = v
            result.clear()
            result.append(k)
        elif v == max:
            result.append(k)
        else:
            pass
    return (max, result)

def reverseDNA(DNA:str) -> str:
    DNA = DNA.upper()
    reverse = ""
    for c in DNA:
        if c == 'A': 
            reverse += 'T'
        elif c == 'T':
            reverse += 'A'
        elif c == 'G':
            reverse += 'C'
        elif c == 'C':
            reverse += 'G'
        else:
            print("Unknown nucleotide ", c, ", aborting")
            return ""
    reverse = reverse[::-1]
    print ("Reverse DNA result: ", reverse)
    return reverse

def patternMatcher(Pattern:str, Genome:str) -> list:
    if Pattern == "" or Genome == "":
        return []
    if len(Pattern) > len(Genome):
        return []
    result = []
    for index in range(len(Genome) - len(Pattern) + 1):
        if Genome[index: index + len(Pattern)] == Pattern:
            result.append(index)
    return result

def FindClumps(Text:str, k:int, L:int, t:int) -> list:
    Patterns = []
    n = len(Text)
    for index in range(n-L):
        Window = Text[index:index+L]
        freqMap = frequencyTable(Window, k)
        for key in freqMap:
            if freqMap[key] >= t: 
                if not key in Patterns:
                    Patterns.append(key)
    # remove duplicates from Patterns
    return Patterns

# MODULE 2

def skew(sequence:str, i:int) -> list[int]:
    skew_array = []
    skew_array.append(0)
    # Note: to optimize memory access: 
    # 1/ store skew_array[j] in a local variable, 
    # 2/ preallocate the array with the right length and populate after to avoid implicit allocations

    for j in range(i + 1):
        if sequence[j] == 'A' or sequence[j] == 'T':
            skew_array.append(skew_array [j])
        elif sequence[j] == 'C':
            skew_array.append(skew_array [j] - 1)
        elif sequence[j] == 'G':
            skew_array.append(skew_array [j] + 1)
    print(skew_array)
    return skew_array

def minimum_skew(genome: str) -> list[int]:
    # Optimization: calculate min while in the loop: avoids 2n complexity, only n
    skew_array = skew(genome, len(genome)-1)
    genome_min = min(skew_array)
    indexes = []
    for i, val in enumerate(skew_array):
        if val == genome_min:
            indexes.append(i)
    print(indexes)
    return indexes

def hamming_distance(str1: str, str2: str) -> int:
    # Optimization: return sum(x != y for x, y in zip(str1, str2, strict=True))
    if len(str1) != len(str2):
        raise ValueError("Mismatch in string lengths")
    result = 0
    for x, y in zip(str1, str2):
        if x != y:
            result += 1
    return result

def approximate_pattern_matching(pattern:str, genome:str, distance:int) -> list:
    result = []
    if len(genome) < len(pattern):
        print("Genome shorter than pattern.")
        return None
    strl = len(pattern)
    for i in range(len(genome) - strl + 1):
        if hamming_distance(genome[i:i+strl], pattern) <= distance:
            result.append(i)
    return result

def approximate_pattern_matching_count(pattern:str, genome:str, distance:int) -> int:
    return len(approximate_pattern_matching(pattern, genome, distance))

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

def max_frequency (frequency: dict):
    maximum = 0
    result = []
    for k, v in frequency.items():
        if v > maximum:
            maximum = v
            result.clear()
            result.append(k)
        elif v == maximum:
            result.append(k)
        else:
            pass
    return result

def find_kmers_with_mismatches(genome:str, k:int, distance:int) -> list:
    ''' Find the most frequent k-mers with mismatches in a string.
    Input: A string Text as well as integers k and d.
    Output: All most frequent k-mers with up to d mismatches in Text.'''
    most_frequent_kmers = []
    frequency_map = dict()

    if k > len(genome):
        print ("k longer than genome. Returning.")
        return None
    
    for i in range(len(genome) - k + 1):
        pattern = genome[i:i+k]
        neighbors = list_pattern_neighbors(pattern, distance)
        for neighbor in neighbors:
            if not neighbor in frequency_map:
                frequency_map[neighbor] = 1
            else:
                frequency_map[neighbor] += 1

    most_frequent_kmers = max_frequency(frequency_map)

    return most_frequent_kmers

def find_kmers_with_mismatches_and_reverse_complements(genome:str, k:int, distance:int) -> list:
    '''
    Find the most frequent k-mers (with mismatches and reverse complements) in a string.
    Input: A DNA string Text as well as integers k and d.
    Output: All k-mers Pattern maximizing the sum Countd(Text, Pattern)+ Countd(Text, Patternrc) over all possible k-mers.
    '''
    most_frequent_kmers = []
    frequency_map = dict()

    if k > len(genome):
        print ("k longer than genome. Returning.")
        return None
    
    for i in range(len(genome) - k + 1):
        pattern = genome[i:i+k]
        rcpattern = reverseDNA(pattern)

        neighbors = list_pattern_neighbors(pattern, distance)
        neighbors += list_pattern_neighbors(rcpattern, distance)

        for neighbor in neighbors:
            if not neighbor in frequency_map:
                frequency_map[neighbor] = 1
            else:
                frequency_map[neighbor] += 1

    most_frequent_kmers = max_frequency(frequency_map)

    return most_frequent_kmers

# MODULE 3

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

# MODULE 4

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