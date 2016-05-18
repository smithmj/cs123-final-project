import numpy as np

MATCH = 2
MISMATCH = -1
GAP = -1

class Cell:
    # This class is for use in our grid of sequence alignment, 
    # each piece of information is later utilized in other functions
    def __init__(self):
        self.val = 0
        self.xcoor = None 
        self.ycoor = None
        self.prev = None

    def __str__(self):
        return(str(self.val))

    def __repr__(self):
        return(str(self.val))

    def __lt__(self, other):
        return self.val < other.val

    def __eq__(self, other):
        return self.val == other.val
    # These are used for the purpose of sorting in max functions

# This function takes a reference sequence and a 'read' sequence to compare to.
# If looking for a local sequence alignment the matrix is not allowed to be negative
# Otherwise the grid created is for global match only with negatives allowed
def place_read(ref_Seq, read, local = False):

    ref_len = len(ref_Seq)
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
        grid[i, 0].xcoor = i
        grid[i, 0].ycoor = 0
        #print(grid)
        #print(grid[i, 0])
        if  i != 0:
            grid[i, 0].prev = grid[i-1, 0]
        for j in range(1, grid.shape[1]):
            grid[i, j] = Cell()
            #print('j:', j)
            grid[i, j].xcoor = i
            grid[i, j].ycoor = j
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
                nucl1 = ref_Seq[j-1]
                nucl2 = read[i-1]
                if nucl1 == nucl2:
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
            #print(grid[i,j])
            #print(grid)
                
    return grid

# This function is designed to 'walk' through the matrix and give back a list of cells 
# corresponding to the best alignment possible. In the case of local alignment we search the
# grid for the highest value and use the 'prev' values to traceback until we reach the top row.
# In a global alignment we start at the bottom right corner and continue until we reach the top left.
def traceback(grid, local = False):
    if local:
        start = max(grid[-1:,]) #we need to go over this again to make sure
    else:
        start = grid[-1,-1]

    rv = []
    curr = start
    while curr != None:
    # the grid creation makes sure that this terminating condition is correct
    # it is incorrect to call this function with local = True when place_read had
    # local = False
        rv.append(curr)
        curr = curr.prev

    return rv


