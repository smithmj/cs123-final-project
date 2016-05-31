import numpy as np

# The scoring system for the alignment 
# algorithm.

# match between base pairs is +2
MATCH = 2
# a mismatch between base pairs
# is penalize -1 
MISMATCH = -1
# a gap between base pairs is
# penalized -1
GAP = -1

class Cell:
    # This class is for use in our grid of sequence alignment, 
    # each piece of information is later utilized in other functions
    def __init__(self, val = 0, row = None, col = None, ref_allele = None, read_allele = None, prev = None):
        
        # store the sequence aligment score for the cell
        self.val = val

        # stores the grid row position of the cell
        self.row = row

        # stores the grid column position of the cell
        self.col = col

        # the allele of the reference sequence corresponding to
        # the column of the grid
        self.ref_allele = ref_allele
        
        # the allele of the read corresponding to the row
        # of the grid 
        self.read_allele = read_allele

        # gives the grid coordinates as a tuple for the cell that points
        # to the current cell
        self.prev = prev

    # Returns the cell corresponding to the coordinate self.prev.
    # Used in the traceback function to construct the alignment
    # from the grid.
    def get_prev_cell(self, grid):
        if self.prev == None:
            return None
        prev_row = self.prev[0]
        prev_col = self.prev[1]
        prev_cell = grid[prev_row, prev_col]
        return prev_cell

    # The __lt__ and __eq__ methods sort the cell based on the 
    # value attribute and are used for the purpose of sorting 
    # in max functions
    def __lt__(self, other):
        return self.val < other.val

    def __eq__(self, other):
        if type(other) != Cell:
            return False
        return self.val == other.val
    

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
    # This class stores information about Single Nucleotide Polymorphisms (SNPs)
    def __init__(self, ref_allele, var_allele, ref_pos, coverage, prop):
        # The reference allele
        self.ref_allele = ref_allele

        # The variant allele
        self.var_allele = var_allele

        # The position of the SNP. Stored as a tuple of (chromosome, base pair)
        self.ref_pos = ref_pos

        # This stores the number of reads that got mapped to the reference 
        # position. This can be used as a quality control measure: if the 
        # coverage is very low, then the alignment was not very reliable
        self.coverage = coverage

        # This stores the proportion of all the reads that were mapped
        # to the reference position that were var_allele. This can be
        # used as a quality control measure: if the proportion is low,
        # the alignment may not be reliable.
        self.prop = prop

    def __str__(self):
        if type(self.ref_pos[1]) == tuple:
            # in this case the SNP is a gap
            string = "chr{} {} {} ".format(self.ref_pos[0], self.ref_pos[1][0], self.ref_pos[1][1])
        else:
            string = "chr{} {} {} ".format(self.ref_pos[0], self.ref_pos[1], self.ref_pos[1])
        string += "{} {} {} {}\n".format(self.ref_allele,self.var_allele, 
                                         self.coverage,self.prob)
        return string

    def __lt__(self, other):
        if (self.ref_pos[0] > other.ref_pos[0]):
        #chrom check
            return False
        elif (self.ref_pos[0] == other.ref_pos[0]):
            return self.ref_pos[1] < other.ref_pos[1] #chrom_pos check if equal
        else:
            return True

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
        # with a gap or multiple gaps

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

# This function is designed to 'walk' through the matrix and give back a list of cells 
# corresponding to the best alignment possible. In the case of local alignment we search the
# grid for the highest value and use the 'prev' values to traceback until we reach the top row.
# In a global alignment we start at the bottom right corner and continue until we reach the top left.
def traceback(grid, align_dict, local = False, chrom = None):
    if local:
        start = find_start_cell(grid[-1:,]) #we need to go over this again to make sure
    else:
        start = grid[-1,-1]

    curr = start
    
    while curr != None and curr.get_prev_cell(grid) != None:
    # the grid creation makes sure that this terminating condition is correct
    # it is incorrect to call this function with local = True when place_read had
    # local = False
        ref_seq_pos = curr.col
        if chrom != None:
            ref_seq_pos = (chrom, ref_seq_pos)
        ref_seq_allele = curr.ref_allele
        if (ref_seq_allele, ref_seq_pos) not in align_dict:
            align_dict[(ref_seq_allele, ref_seq_pos)] = {}
        prev_cell = curr.get_prev_cell(grid)

        if curr.col == prev_cell.col:
            # in this case there is a gap in the reference
            insertion = curr.read_allele 
            curr = prev_cell
            prev_cell = curr.get_prev_cell(grid)
            while prev_cell != None and curr.col == prev_cell.col:
                insertion = curr.read_allele + insertion
                curr = prev_cell
                prev_cell = curr.get_prev_cell(grid)
            insertion = curr.read_allele + insertion
            if insertion not in align_dict[(ref_seq_allele, ref_seq_pos)]:
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
    return align_dict

# Givien the alignment dictionary, this function finds the genotype for 
# each position in the reference genome.
def genotype(align_dict, save_file = None):
    if save_file != None:
        f = open(save_file, "w")
        f.write("Ref-Position Ref-Allele Var-Allele Read-Coverage Proportion\n")
    else:
        snps = []

    gap_snps = []
    for (ref_allele, ref_pos) in align_dict:
        max_val = 0
        most_probable = None
        coverage = 0
        for allele in align_dict[(ref_allele, ref_pos)]:
            # The genotype is determined by the allele with 
            # the highest count. 

            # The coverage identifies how many reads were aligned
            # at this position. If the coverage is too low, the
            # alignment may not be reliable.
            coverage += align_dict[(ref_allele, ref_pos)][allele]    

            if align_dict[(ref_allele, ref_pos)][allele] >= max_val:
                max_val = align_dict[(ref_allele, ref_pos)][allele]
                most_probable = allele
        # The proportion of the most-probable allele out of all
        # the alleles that were aligned to this position of ref_seq.
        # If prop is low, the alignment may not be reliable
        prop = align_dict[(ref_allele, ref_pos)][most_probable] / coverage

        if most_probable == '-':
            snp = SNP(ref_allele, allele, ref_pos, coverage, prop)
            gap_snps.append(snp)

        # "N" is the symbol assigned asigned to alleles during the actual
        # DNA sequencing when there is not enough information to determine
        # the nucleotide. We are only interested in alleles that are
        # "A", "G", "T", or "C" and positions where ref_allele is not the
        # same as the ref_allele.
        elif most_probable != "N" and ref_allele != "N" and most_probable != ref_allele:
            snp = SNP(ref_allele, allele, ref_pos, coverage, prop)
            if save_file != None:
                f.write(str(snp))
            else:
                snps.append(snp)

    # This command condenses all of the adjacent gap SNPs into a single SNP.
    # For example, if there is a gap in the reads at chromosome 1 base at 
    # base pairs 1, 2, 3, only one SNP will be returned instead of three 
    # individual SNPs.
    combined_gap_snps = gap_handler(gap_snps)
    for snp in combined_gap_snps:
        if save_file != None:
            f.write(str(snp))
        else:
            snps.append(snp) 

    if save_file != None:
        f.close()
        return
    else:
        return snps

def avg(l):
    return sum(l) / len(l)

# This function handles the special case of SNP gaps 
def gap_handler(gaps):
    gaps.sort()    
    glen = len(gaps)
    i = 0
    rv = []
    while(i < glen - 1):
        start = i
        ra = ""
        covs = []
        props = []
        while (gaps[i].ref_pos[0] == gaps[i+1].ref_pos[0] and gaps[i].ref_pos[1] == (gaps[i+1].ref_pos[1] - 1)): #on same chrom and adjacent
            ra += gaps[i].ref_allele
            covs.append(gaps[i].coverage)
            props.append(gaps[i].prop)
            i += 1
        if (i == start):
            rv.append(gaps[i])
        else:
            ra += gaps[i].ref_allele
            covs.append(gaps[i].coverage)
            props.append(gaps[i].prop)
            print(ra)
            print(len(ra))
            var_allele = "-" * len(ra)
            r = SNP(ra, var_allele,(gaps[i].ref_pos[0],(gaps[start].ref_pos[1],gaps[i].ref_pos[1])),avg(covs),avg(props))
            rv.append(r)

        i += 1
    if (i == glen - 1):
        # this checks if the last one was a singleton
        rv.append(gaps[i])
    return rv

# this function was used for testing purposes so we could visualize
# what an alignment looked like.
def print_alignment(align_dict, ref_seq):
    align_list = []
    for ref_allele, ref_pos in align_dict:
        max_val = 0
        most_probable = None
        for allele in align_dict[(ref_allele, ref_pos)]:    
            if align_dict[(ref_allele, ref_pos)][allele] >= max_val:
                max_val = align_dict[(ref_allele, ref_pos)][allele]
                most_probable = allele
        align_list.append((ref_pos, ref_allele, most_probable))
    align_list.sort()

    # this is the first position of the reference sequence that the 
    # read has been aligned to. If it is not 1, then there are initial
    # gaps in the read sequence
    ref_str = ""
    read_str = ""
    first_ref_pos = align_list[0][0][1]
    for i in range(1, first_ref_pos):
        ref_str += ref_seq[i - 1]
        read_str += "-"

    for i in range(len(align_list)):
        ref_str += align_list[i][1]
    print(ref_str)
    
    for i in range(len(align_list)):
        read_str += align_list[i][2]
    print(read_str)
    
    return
