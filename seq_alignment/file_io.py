import re
import numpy as np
import alignment

# This function takes in the reference sequence
# filename and returns files for the sequence
# for each chromosome
def split_ref_seq(seq_filename):
    chrom_header_re = "^\d+|X|Y"

    f = open(seq_filename, 'r')
    complete = False
    while not complete:
        line = f.readline()
        if re.findall(chrom_header_re, line) != []:
            # in this case we have reached a new chromosome
            # and we want to start a new file
            chrom = re.findall[0]
            new_file_name = "chrom_{}.txt".format(chrom)
            new_file = open("chrom_{}.txt", "w")
        new_file.write(line + "\n")
        line = f.readline()
        if line == '':
            complete = True
    return

def convert_seq_file_to_str(filename):
    f = open(filename, "r")
    # don't need the header line
    f.readline()
    seq = ""
    for line in f:
        seq += line.strip()
    return seq

def prepair_reads(filename):
    read_re = "[A,G,C,T,N]+"
    chrom_list = list(range(1, 23)) + ['X', 'Y']
    f = open(filename, "r")
    new_file = open("mrjob_short_reads.txt", "w")
    read_num = 1
    for line in f:
        if re.match(read_re, line) != None:
            for chrom in chrom_list:
                new_file.write("{}, {}, {}".format(read_num, chrom, line))
        read_num += 1
    return

def convert_grid_to_str(chrom, grid):
    grid_rows = []
    for i in range(grid.shape[0]):
        row = []
        for j in range(grid.shape[1]):
            cell = grid[i,j]
            row.append(str(cell))
        grid_rows.append("+".join(row))
    string = "{}|".format(chrom) +  "|".join(grid_rows)
    return string 

def convert_str_to_grid(string):
    string_info = string.split("|")
    chrom = string_info[0]
    grid_str_rows = string_info[1:]
    grid_rows = []
    for str_row in grid_str_rows:
        str_cells = str_row.split("+")
        cells = []
        for str_cell in str_cells:
            cell_info = str_cell.split("\n")
            val_re = "\d+"
            val = int(re.findall(val_re, cell_info[0])[0])
            coor_re = "\d+"
            coor = re.findall(coor_re, cell_info[1])
            row = int(coor[0])
            col = int(coor[1])
            prev_re = "\d+|None"
            prev_coor = re.findall(prev_re, cell_info[2])
            if len(prev_coor) == 2:
                # in this case the cell has a traceback 
                prev = (int(prev_coor[0]), int(prev_coor[1]))
            else:
                # in this case the cell has no traceback
                prev = None
            cell = alignment.Cell(val, row, col, prev)
            cells.append(cell)
        grid_rows.append(cells)
    grid = np.array(grid_rows, dtype = object)
    return chrom, grid

