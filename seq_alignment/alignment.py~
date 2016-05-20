import numpy as np
import operator as op
MATCH = 2
MISMATCH = -1
GAP = -1

class Cell:
    # This class is for use in our grid of sequence alignment, 
    # each piece of information is later utilized in other functions
    def __init__(self):
        self.val = 0
        self.row = None 
        self.col = None
        self.prev = None

    def __str__(self):
        return(str(self.val))

    def __repr__(self):
        return(str(self.val))

    def __lt__(self, other):
        return self.val < other.val

    def __eq__(self, other):
        if type(other) != Cell:
            return False
        return self.val == other.val
    # These are used for the purpose of sorting in max functions

class SNP:

    def __init__(self, ref_allele, var_allele, ref_pos, coverage, prob):
        self.ref_allele = ref_allele
        self.var_allele = var_allele
        self.ref_pos = ref_pos
        self.coverage = coverage
        self.prob = prob

    def __str__(self):
        return('{} -> {}'.format(self.ref_allele, self.var_allele))

    def __repr__(self):
        return('{} -> {}'.format(self.ref_allele, self.var_allele))

# This function takes a reference sequence and a 'read' sequence to compare to.
# If looking for a local sequence alignment the matrix is not allowed to be negative
# Otherwise the grid created is for global match only with negatives allowed
def place_read(ref_seq, read, local = False):

    ref_len = len(ref_seq)
    read_len = len(read)
    #grid = np.array([[Cell()] * (ref_len + 1)] * (read_len + 1))
    grid = np.zeros((read_len + 1, ref_len + 1), dtype = object)
        # the grid needs an extra column and row in the case that the alignment starts
        # with a gap or gaps

    # double for loop sets up the initialization of the alignment matrix
    for i in range(grid.shape[0]):
        grid[i, 0] = Cell()
        if local:
            grid[i, 0].val = max(0, GAP * i)
        else:
            grid[i, 0].val = GAP * i
        grid[i, 0].row = i
        grid[i, 0].col = 0
        #print(grid)
        #print(grid[i, 0])
        if  i != 0:
            grid[i, 0].prev = grid[i-1, 0]
        for j in range(1, grid.shape[1]):
            grid[i, j] = Cell()
            #print('j:', j)
            grid[i, j].row = i
            grid[i, j].col = j
            if i == 0:
                if local:
                    grid[i,j].val = max(0, GAP * j)
                    grid[i,j].prev = None
                    # this prev init sets up terminating condition for traceback (local)
                else:
                    grid[i,j].val = GAP * j
                    grid[i,j].prev = grid[i,j-1]
            else:
                # This takes into account the three paths from which a value can be determined
                # and selects the one with the highest value
                ref_nucl = ref_seq[j-1]
                read_nucl = read[i-1]
                if ref_nucl == read_nucl:
                    diagonal = MATCH
                else:
                    diagonal = MISMATCH
                max_val = diagonal + grid[i-1, j-1].val
                grid[i,j].prev = grid[i-1,j-1]

                if GAP + grid[i-1,j].val > max_val:
                    max_val = GAP + grid[i-1, j].val
                    grid[i, j].prev = grid[i-1, j]

                if GAP + grid[i, j-1].val > max_val:
                    max_val = GAP + grid[i, j-1].val
                    grid[i, j].prev = grid[i, j-1]

                if local:
                    grid[i,j].val = max(0, max_val)
                else:
                    grid[i,j].val = max_val
                
    return grid

# This function is designed to 'walk' through the matrix and give back a list of cells 
# corresponding to the best alignment possible. In the case of local alignment we search the
# grid for the highest value and use the 'prev' values to traceback until we reach the top row.
# In a global alignment we start at the bottom right corner and continue until we reach the top left.


def find_start_cell(grid):
    max_val = 0
    start_cell = None
    for i in range(grid.shape[0]):
        for j in range(grid.shape[1]):
            if grid[i,j].val >= max_val:
                max_val = grid[i,j].val
                start_cell = grid[i,j]
    return start_cell


def traceback(grid, ref_seq, read, local = False, mr = True, align_dict = None):
    if local:
        start = find_start_cell(grid[-1:,]) #we need to go over this again to make sure
    else:
        start = grid[-1,-1]

    curr = start
    while curr.prev != None:
    # the grid creation makes sure that this terminating condition is correct
    # it is incorrect to call this function with local = True when place_read had
    # local = False
        ref_seq_pos = curr.col - 1
        if align_dict != None and ref_seq_pos not in align_dict:
            align_dict[ref_seq_pos] = {}
        prev_cell = curr.prev
        if curr.col == prev_cell.col:
            # in this case there is a gap in the reference
            insertion = read[curr.row - 1] 
            curr = prev_cell
            prev_cell = curr.prev
            while curr.col == prev_cell.col:
                insertion = read[curr.row - 1] + insertion
                curr = curr.prev
                prev_cell = curr.prev
            insertion = read[curr.row - 1] + insertion
            if align_dict != None and (insertion not in align_dict[ref_seq_pos]):
                align_dict[ref_seq_pos][insertion] = 0
            if mr:
                yield ref_seq_pos, insertion
            if align_dict != None:
                align_dict[ref_seq_pos][insertion] += 1
            
            curr = curr.prev

        elif curr.row == prev_cell.row:
            # in this case there is a gap in the read
            nucl = '-'
            if align_dict != None and (nucl not in align_dict[ref_seq_pos]):
                align_dict[ref_seq_pos][nucl] = 0
            if mr:
                yield ref_seq_pos, nucl
            if align_dict != None:
                align_dict[ref_seq_pos][nucl] += 1
            
            curr = prev_cell 
        else:
            # in this case there are no gaps
            nucl = read[curr.row - 1]
            if not mr and (nucl not in align_dict[ref_seq_pos]):
                align_dict[ref_seq_pos][nucl] = 0
            if mr:
                yield ref_seq_pos, nucl
            if align_dict != None:
                align_dict[ref_seq_pos][nucl] += 1
             
            curr = prev_cell

    if mr:
        return
    if align_dict != None:
        return align_dict

def genotype(align_dict, ref_seq):
    #rv =[]
    #for i in range(len(ref_dic)):
    #    rv[i] = max(align_dict[i].iteritems(),key=op.itemgetter(1))[0]
    #return rv

    snps = []
    for ref_pos in align_dict:
        max_val = 0
        most_probable = None
        coverage = 0
        for allele in align_dict[ref_pos]:
            ## Right now we're determining the genotype based on the 
            ## highest count. 

            # The coverage identifies how many reads were aligned
            # at this position. If the coverage is too low, the
            # alignment may not be reliable.
            coverage += align_dict[ref_pos][allele]    

            if align_dict[ref_pos][allele] >= max_val:
                max_val = align_dict[ref_pos][allele]
                most_probable = allele
        # The proportion of the most-probable allele out of all
        # the alleles that were aligned to this position of ref_seq.
        # If prob is low, the alignment may not be reliable
        prob = align_dict[ref_pos][most_probable] / coverage
        ref_allele = ref_seq[ref_pos]

        if most_probable == '-':


        if most_probable != ref_allele:
            snp = SNP(ref_allele, allele, ref_pos, coverage, prob)
            snps.append(snp)
    return snps







