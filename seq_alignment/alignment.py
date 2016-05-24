import numpy as np
import operator as op
MATCH = 2
MISMATCH = -1
GAP = -1

class Cell:
    # This class is for use in our grid of sequence alignment, 
    # each piece of information is later utilized in other functions
    def __init__(self, val = 0, row = None, col = None, ref_allele = None, read_allele = None, prev = None):
        self.val = val
        self.row = row
        self.col = col
        self.ref_allele = ref_allele
        self.read_allele = read_allele
        self.prev = prev

    def get_prev_cell(self, grid):
        if self.prev == None:
            return None
        prev_row = self.prev[0]
        prev_col = self.prev[1]
        prev_cell = grid[prev_row, prev_col]
        return prev_cell

    def __lt__(self, other):
        return self.val < other.val

    def __eq__(self, other):
        if type(other) != Cell:
            return False
        return self.val == other.val
    # These are used for the purpose of sorting in max functions

    def __str__(self):
        val_str = "val: {}\n".format(self.val)
        coor_str = "coor: {}\n".format((self.row, self.col))
        ref_allele_str = "ref allele: {}\n".format(self.ref_allele)
        read_allele_str = "read allele: {}\n".format(self.read_allele) 
        prev_str = "prev: {}\n".format(self.prev)
        string = val_str + coor_str + ref_allele_str + read_allele_str + prev_str
        return string 

    def __repr__(self):
        return self.__str__()

class SNP:

    def __init__(self, ref_allele, var_allele, ref_pos, coverage, prob):
        self.ref_allele = ref_allele
        self.var_allele = var_allele
        self.ref_pos = ref_pos
        self.coverage = coverage
        self.prob = prob

    def __str__(self):
        string = "{},{},{},{},{}\n".format(self.ref_pos, self.ref_allele,
                                               self.var_allele, self.coverage,
                                               self.prob)
        return string
        #return('{} -> {}'.format(self.ref_allele, self.var_allele))

    def __repr__(self):
        return self.__str__()

# This function takes a reference sequence and a 'read' sequence to compare to.
# If looking for a local sequence alignment the matrix is not allowed to be negative
# Otherwise the grid created is for global match only with negatives allowed
def place_read(ref_seq, read, local = False):

    ref_len = len(ref_seq)
    read_len = len(read)
    grid = np.zeros((read_len + 1, ref_len + 1), dtype = object)
        # the grid needs an extra column and row in the case that the alignment starts
        # with a gap or gaps

    # double for loop sets up the initialization of the alignment matrix
    for i in range(grid.shape[0]):
        grid[i, 0] = Cell(ref_allele = ref_seq[0], read_allele = read[i-1])
        if local:
            grid[i, 0].val = max(0, GAP * i)
        else:
            grid[i, 0].val = GAP * i
        grid[i, 0].row = i
        grid[i, 0].col = 0
        if  i != 0:
            grid[i, 0].prev = (grid[i-1, 0].row, grid[i-1, 0].col)
        for j in range(1, grid.shape[1]):
            grid[i, j] = Cell(ref_allele = ref_seq[j - 1], read_allele = read[i - 1])
            grid[i, j].row = i
            grid[i, j].col = j
            if i == 0:
                if local:
                    grid[i,j].val = max(0, GAP * j)
                    grid[i,j].prev = None
                    # this prev init sets up terminating condition for traceback (local)
                else:
                    grid[i, j].val = GAP * j
                    grid[i, j].prev = (grid[i,j-1].row, grid[i,j-1].col)
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
                grid[i, j].prev = (grid[i-1, j-1].row, grid[i-1, j-1].col)

                if GAP + grid[i-1, j].val > max_val:
                    max_val = GAP + grid[i-1, j].val
                    grid[i, j].prev = (grid[i-1, j].row, grid[i-1, j].col)

                if GAP + grid[i, j-1].val > max_val:
                    max_val = GAP + grid[i, j-1].val
                    grid[i, j].prev = (grid[i, j-1].row, grid[i, j-1].col)

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

def find_alignment_score(grid):
    start_cell = find_start_cell(grid)
    return start_cell.val

# need to rework traceback to store ref_seq nucl and 
# read nucl so that ref_seq and read do not need to be passed in
def traceback(grid, align_dict, local = False, chrom = None):
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
        if chrom != None:
            ref_seq_pos = "chr{}:{}".format(chrom, ref_seq_pos)
        ref_seq_allele = curr.ref_allele
        if (ref_seq_allele, ref_seq_pos) not in align_dict:
            align_dict[(ref_seq_allele, ref_seq_pos)] = {}
        prev_cell = curr.get_prev_cell(grid)
        
        if curr.col == prev_cell.col:
            # in this case there is a gap in the reference
            insertion = curr.read_allele 
            curr = prev_cell
            prev_cell = curr.get_prev_cell(grid)
            while curr.col == prev_cell.col:
                insertion = curr.read_allele + insertion
                curr = curr.prev
                prev_cell = curr.get_prev_cell(grid)
            insertion = curr.read_allele + insertion
            if insertion not in align_dict[ref_seq_pos]:
                align_dict[(ref_seq_allele, ref_seq_pos)][insertion] = 0
            align_dict[(ref_seq_allele, ref_seq_pos)][insertion] += 1            
            
            curr = prev_cell

        elif curr.row == prev_cell.row:
            # in this case there is a gap in the read
            nucl = "-"
            if nucl not in align_dict[(ref_seq_allele, ref_seq_pos)]:
                align_dict[(ref_seq_allele, ref_seq_pos)][nucl] = 0
            align_dict[(ref_seq_allele, ref_seq_pos)][nucl] += 1
            
            curr = prev_cell 
        
        else:
            # in this case there are no gaps
            nucl = curr.read_allele
            if nucl not in align_dict[(ref_seq_allele, ref_seq_pos)]:
                align_dict[(ref_seq_allele, ref_seq_pos)][nucl] = 0
            align_dict[(ref_seq_allele, ref_seq_pos)][nucl] += 1
             
            curr = prev_cell
    #print("inside function align_dict:", align_dict)
    return align_dict

# rework genotype so it doesn't need to take in ref_seq
def genotype(align_dict, save_file = None):
    #rv =[]
    #for i in range(len(ref_dic)):
    #    rv[i] = max(align_dict[i].iteritems(),key=op.itemgetter(1))[0]
    #return rv
    if save_file != None:
        f = open(save_file, "w")
        f.write("Ref Position,Ref Allele,Var Allele,Read Coverage,Proportion\n")
    else:
        snps = []

    for (ref_allele, ref_pos) in align_dict:
        max_val = 0
        most_probable = None
        coverage = 0
        for allele in align_dict[(ref_allele, ref_pos)]:
            ## Right now we're determining the genotype based on the 
            ## highest count. 

            # The coverage identifies how many reads were aligned
            # at this position. If the coverage is too low, the
            # alignment may not be reliable.
            coverage += align_dict[(ref_allele, ref_pos)][allele]    

            if align_dict[(ref_allele, ref_pos)][allele] >= max_val:
                max_val = align_dict[(ref_allele, ref_pos)][allele]
                most_probable = allele
        # The proportion of the most-probable allele out of all
        # the alleles that were aligned to this position of ref_seq.
        # If prob is low, the alignment may not be reliable
        prob = align_dict[(ref_allele, ref_pos)][most_probable] / coverage

        if most_probable == '-':
            snp = SNP(ref_allele, allele, ref_pos, coverage, prob)
            gap_alleles.append(snp)

        elif most_probable != ref_allele:
            snp = SNP(ref_allele, allele, ref_pos, coverage, prob)
            if save_file != None:
                f.write(str(snp))
            else:
                snps.append(snp)

    # write a gap-handling function that takes in gap_alleles, a list of snps
    # and combines the snps based on the reference position.
    # ie: chr1:2, chr1:3, chr1:4 -> chr1:2-4

    if save_file != None:
        f.close()
        return
    else:
        return snps


# THINGS WE NEED TO DO

# Be able to collapse/expand grid to/from a single list
# Handle gaps/deletions



