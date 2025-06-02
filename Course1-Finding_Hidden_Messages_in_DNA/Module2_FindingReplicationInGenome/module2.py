'''Module 2 of Finding hidden messages in DNA by University of San Diego.'''

# from ..Module1_FindingOriC import module1
# import sys
# from pathlib import Path
# sys.path.insert(0, str(Path(__file__).resolve().parent.parent / "Module1_FindingOriC"))
# import module1

nucleotides = {'A', 'T', 'C', 'G'}

# These recursive functions should not be used: quickly reaching max number of recursion
def skew_recursive(sequence:str, i:int) -> int:
    '''Skew computes the number of G-C at a given index i in a DNA sequence. Solved recursively.'''
    if i == -1:
        result = 0
    elif sequence[i] == 'A' or sequence[i] == 'T':
        result = skew_recursive(sequence, i-1)
    elif sequence[i] == 'C':
        result = skew_recursive(sequence, i-1) - 1
    elif sequence[i] == 'G':
        result = skew_recursive(sequence, i-1) + 1
    else:
        print("Wrong input at index i = ", i)
        result = -1
    print(result)
    return result

def minimum_skew_recursive(genome: str) -> list:
    genome_min = 0
    indexes = []
    for j in range(len(genome)):
        res = skew_recursive(genome[0:j], j -1)
        if res < genome_min:
            genome_min = res
            indexes.clear()
            indexes.append(j)
        elif res == genome_min:
            indexes.append(j)
    print("Minimum skew indexes: ", indexes)
    return indexes

# Used code starting here

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
    return reverse

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

if __name__ == "__main__":
    # skew_recursive('CATGGGCATCGGCCATACGCC', 14)
    skew('CATGGGCATCGGCCATACGCC', 14)
    minimum_skew("TAAAGACTGCCGAGAGGCCAACACGAGTGCTAGAACGAGGGGCGTAAACGCGGGTCCGAT")
    print(approximate_pattern_matching('ATTCTGGA', 'CGCCCGAATCCAGAACGCATTCCCATATTTCGGGACCACTGGCCTCCACGGTACGGACGTCAATCAAAT', 3))
