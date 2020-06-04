f = open('20a.fa', 'r')

# process fasta file into string of DNA
def process_fasta(file):
    s = file.read()
    
    # split by new line
    s = s.split('\n')

    # get rid of first line (identifier)
    s = s[1:]

    # return the rest joined as one long string
    return ''.join(s)

mapping = {
    'G': 'C',
    'C': 'G',
    'A': 'T',
    'T': 'A'
}

# process reverse
def rev(s):
    res = ''
    s = reversed(s)
    for c in s:
        res += mapping[c]

    return res

# actual function that churns and finds the repeats of length 13 and their reverses
# takes in a long string of base pairs
def process(s):
    
    # this dictionary maps from a 13-char substring of the DNA to the number of times it occurs in the string
    memo = dict()

    # results
    res = [] 

    for i in range(len(s) - 13 + 1):
        curr_s = s[i : i+13]
        memo[curr_s] = memo.get(curr_s, 0) + 1

    # we don't reprocess those already found
    found = set()

    # now we iterate through the dictionary 
    for (sub, cnt) in memo.items():
        if cnt >= 2 and not (rev(sub) in found):
            sub_rev = rev(sub)

            a = [sub, cnt]
            b = [sub_rev, memo.get(sub_rev, 0)]
            found.add(sub)
            res.append([a, b])


    return res

DNA = process_fasta(f)
results = process(DNA)

for res in results:
    print(res)