print ("####################### PATTERN COUNT ###########################")

Dataset = "CATCGTCTCCGTAAGCAAAAAGCAAATAAGCAAAAAGCAAATCCATTCACTTGAGGGCGTCAAAGCAAAAAGCAAAAAGCAAATAAGCAAAAAGCAAAAAGCAAATAAGCAAAGAAAAGCAAAAAGCAAATTAAGCAAAAAGCAAACGTAAGCAAACAAGCAAAAAGCAAAGGAAGCAAAATTAAGCAAAAAGCAAAAAGCAAATAAGCAAAATGAAGCAAAGAAGCAAAAAGCAAAGAAGCAAAGAAGCAAAACTTAAGCAAACAAGCAAAAAACGTAAGCAAAGTAACTAAGCAAATAGAAAGCAAAGGGAAGCAAATCAAAAAGCAAATAGCCAAGCAAAAGAAGCAAACGTCATGGAAGCAAAAAGCAAATTAAGCAAAGTGAAGCAAAAAGCAAAGTCTCCCATGAAGCAAATAAGCAAAAAAGCAAAAAGCAAAAAAGCAAAGAACATTCAAAGCAAACGAAGCAAAAAGCAAAAAGCAAAAGCCAAAGCAAATAGCTAATAAGCAAAAAGCAAAAAAGCAAAGAGGAAGCAAAAAGCAAAGGGAAAAGCAAACCAATTGAAGCAAAGCCGGCAAGCAAAAAGCAAAGAGCTAAGCAAAAAGCAAATAGGAAAAAAGCAAAGAAGCAAATTGGTAAAGCAAAAAGCAAAAAGCAAAAGTCCCACTCAAGCAAAGATGGGAAGCAAAAGAGTCGAAGCAAATCAAGCAAAAAGCAAAAAGCAAACCATTAAGCAAAAAGCAAATAAGCAAAAGGCGAAAAGCAAAACAAGCAAAGAAGCAAAGATAAGCAAAGTAAGCAAAATAAGCAAACCTAAGCAAAAAGCAAATCTAAGCAAAGAAGCAAAAAGCAAATATATGCGCAATCATAGCAAGCAAAAAGCAAAAGAAGCAAAGAAAAAGCAAAGTAAGCAAATTCAATCTT"
Pattern = "AAGCAAAAA"

def PatternCount(Text, Pattern):
  count = 0
  lt = len(Text)
  lp = len(Pattern)
  for i in range(lt-lp):
    if Text[i : i + len(Pattern)] == Pattern:
      count = count + 1
  print(count)
  return count

PatternCount(Dataset, Pattern)

print ("####################### FIND K-MERS ###########################")

def findKmer(s:str, k: int) -> tuple[str, int] | None:
    result = ("", 0)
    if k > len(s):
        return result
    for i in range(len(s)-k):
        count = 0 
        slice = s[i: i+k]
        print(slice)
        for j in range(i+1, len(s)-k+1):
            if s[j:j+k] == slice:
                count+=1
            if count > result[1]:
                result = (slice, count)
    print("Done, result: ", result)
    return result

findKmer("abcd", 2)
findKmer("abcabcde", 2)

print ("############### Better algo #############")

def frequencyTable (Text:str, k:int):
    frequency = dict()
    lenText = len(Text)
    for i in range(len(Text)-k +1):
        pattern = Text[i:i+k]
        if not pattern in frequency:
            frequency[pattern] = 1
        else:
            frequency[pattern] += 1
    print("Text: ", Text, ", k: ", k, ", result: ", frequency)
    return frequency

def maxFrequency (frequency: dict):
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

def solve(Text:str, k:int):
    return maxFrequency(frequencyTable(Text, k))

# (solve("ABCTUYABCTUTIOTCTU", 3))

#input = "GAAGTCCAGAAGTCCAACAGATCGCTCATACTGCGAGTAAGAGTGAGTAAGAGTGAGTAAGAGTTACCACGGTACCACGGCATACTGCTACCACGGACAGATCGCTTACCACGGACAGATCGCTGAGTAAGAGTGAGTAAGAGTTACCACGGTACCACGGACAGATCGCTACAGATCGCTACAGATCGCTGAAGTCCACATACTGCTACCACGGCATACTGCGAGTAAGAGTTACCACGGGAAGTCCAGAGTAAGAGTGAGTAAGAGTGAAGTCCACATACTGCGAGTAAGAGTTACCACGGACAGATCGCTTACCACGGCATACTGCACAGATCGCTTACCACGGACAGATCGCTGAGTAAGAGTCATACTGCCATACTGCCATACTGCGAGTAAGAGTACAGATCGCTCATACTGCTACCACGGGAAGTCCATACCACGGACAGATCGCTCATACTGCCATACTGCACAGATCGCTGAGTAAGAGTGAAGTCCATACCACGGGAGTAAGAGTGAAGTCCAGAAGTCCAGAGTAAGAGTCATACTGCTACCACGGTACCACGGACAGATCGCTGAAGTCCAGAAGTCCACATACTGCCATACTGCCATACTGCACAGATCGCTACAGATCGCTACAGATCGCTTACCACGGACAGATCGCTTACCACGGGAAGTCCAGAAGTCCACATACTGCTACCACGGTACCACGGTACCACGGCATACTGCGAAGTCCACATACTGCCATACTGCACAGATCGCTCATACTGCCATACTGCTACCACGGTACCACGGCATACTGCCATACTGCGAGTAAGAGTCATACTGCACAGATCGCTACAGATCGCTACAGATCGCTCATACTGCTACCACGGGAGTAAGAGTGAAGTCCAGAGTAAGAGTACAGATCGCTACAGATCGCTACAGATCGCTCATACTGCCATACTGCTACCACGGGAAGTCCAGAGTAAGAGT"
#print(solve(input, 14))

# print(solve("CGCCTAAATAGCCTCGCGGAGCCTTATGTCATACTCGTCCT", 3))

print ("############### Reverse DNA #############")

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

reverseDNA("abcd")

reverseDNA("CGTA")

print ("############### Find a given pattern in a string #############")

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

print ("pattern Matcher: ", patternMatcher("ab", "dabdabuab"))
pm = patternMatcher("GTGCCAAGT","GTGCCAAGTGTGCCAAGTGTGCCAAGTGTGCCAAGT")
print(" ".join(str(x) for x in pm))
