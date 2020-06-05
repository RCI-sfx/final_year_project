'''
This programme uses a brute force method to predict 
many inversions/shuffling combinations until the 
DNA binding arrays match with a continuous section of the target sequence. 

To account for splicing please put your target sequence intersperced by Ns.

The script can only invert the string starting at an odd position and 
ending at an even position thus the length of inversion events 
possible were: 2,6,10,14, 18 and 22.  

define an 'action' as a reversal of a substring, starting from index 2*k,
of length 2+4n, of a string. n and k are integers >= 0

define 'fit', as in a 'fits' into b, len(a) <= len(b), as b having a substring
that can 'bind' to a

find the minimum number of actions required by a source string to 'fit' into a
destination string

'''
#gives use the IUPAC nucleotides found at https://www.bioinformatics.org/sms/iupac.html to enter modules. 

fitDict = {
    'G': {'G'},
    'C': {'C'},
    'A': {'A'},
    'T': {'T'},
    'M': {'A', 'C'},
    'N': {'A', 'T', 'G', 'C'},
    'W': {'A', 'T'},
    'R': {'G', 'A'},
    'Y': {'C', 'T'},
    'k': {'G', 'T'}
}

# returns all the strings possible and current flips made in form [i, j] where flip on index i,j
# given the current string and the previous flips made
def helper(s_pair):
    res = []
    s, curr_pairs = s_pair
    for i in range(0, len(s)+1, 2):
        for j in range(i+2, len(s)+1, 4):
            curr_s = s[:i] + s[i:j][::-1] + s[j:]
            res.append((curr_s, curr_pairs + [[i,j]]))

    return res

# returns (True, start_index) if a fits into b
# returns (False, None) else
# len(a) <= len(b)
def fits(a, b):
    if (len(b) < len(a)):
        raise Exception('b shorter than a')
    
    # return True if even-numbered index characters are the same
    # len(c) == len(d) assumed
    def help(c, d):
        if (len(c) != len(d)):
            raise Exception('c and d not equal in length')

        for i in range(0, len(c), 2):
            if (not (d[i] in fitDict[c[i]])): 
                return False

        return True

    for i in range(0, len(b) - len(a) + 1, 2):
        if (help(a, b[i:i+len(a)])):
            return (True, i)
    return (False, None)


# BFS brute force
# returns all moves (pairs of indexes of flips) made and the start index
# in which src fits in dst
def func(src, dst):
    cnt = 0
    srcs = [(src, [])]
    next_srcs = []

    tried = set()
    
    while True:
        for s_pair in srcs:
            s, moves = s_pair
            if (s in tried):
                continue
            tried.add(s)
            fit, idx = fits(s, dst)
            if (fit):
                return (moves, idx)
        
            next_srcs += helper(s_pair)

        cnt += 1
        print(f'Processing...{cnt}')
        srcs = next_srcs
        next_srcs = []

# aids visualising the flips
def vis(src, dst, moves, start_idx):
    cnt = 1
    for move in moves:
        i, j = move
        print(f'Move {cnt}, {move}')
        a = (f'{src[0:i]} >{src[i:j]}< {src[j:]}')
        b = (f'{src[0:i]} >{src[i:j][::-1]}< {src[j:]}')
        print(f'{a} ====> {b}')
        print()

        src = src[0:i] + src[i:j][::-1] + src[j:]
        cnt += 1


    print(f'Binding at index {start_idx}')
    print(f'src: {src}')
    print(f'dst: {dst[:start_idx]} >{dst[start_idx:start_idx+len(src)]}< {dst[start_idx+len(src):]}')
    


src = 'GCATGACTACGTTGCATCAGTACG'
dst = 'CNGNGNANGNTNANCNTNGNTNCNCNTNCNCNGN'
print(f'src: {src}')
print(f'dst: {dst}')
moves, idx = func(src, dst)

print(vis(src, dst, moves, idx))

print(f'Minimum: {len(moves)} moves')
