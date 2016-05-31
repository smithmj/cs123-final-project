import re
import numpy as np
import alignment

# This function takes in the reference sequence
# filename and returns files for the sequence
# for each chromosome
def split_ref_seq(seq_filename):
    chrom_header_re = "^\d+|X|Y"

    f = open(seq_filename, 'r')
    new_file = None
    line = f.readline()
    if line == '':
        complete = True
    else:
        complete = False
    while not complete:
        #line = f.readline()
        line_re = re.findall(chrom_header_re, line)
        if line_re != []:
            # in this case we have reached a new chromosome
            # and we want to start a new file
            if new_file != None:
                new_file.close()
            chrom = line_re[0]
            new_file_name = "chrom_{}.txt".format(chrom)
            new_file = open(new_file_name, "w")
        else:
            new_file.write(line)
        line = f.readline()
        if line == '':
            complete = True
    f.close()
    new_file.close()
    return

# This function was used for our own testing
# purposes. It converts a .txt file which we 
# assume to be entirely sequence characters
# to a string.
def convert_txt_to_str(filename):
    f = open(filename, "r")
    seq = ""
    for line in f:
        seq += line.strip()
    return seq

# A .fasta file has one header line and then 
# multiple lines of sequence. This function
# converts a .fasta format file to a
# string sequence.
def convert_fasta_to_str(filename):
    f = open(filename, "r")
    f.readline()
    seq = ""
    for line in f:
        seq += line.strip()
    return seq

# This function prepares the file for MRJob
# by creating a file that creates one line for
# every read-chromosome pair. The first field 
# is the read number, the second field is the
# chromosome name, and the third field is the 
# actual read nucleotide sequence.
def prepair_reads(filename, new_filename):
    chrom_list = list(range(1, 23)) + ['X', 'Y']
    f = open(filename, "r")
    new_file = open(new_filename, "w")
    counter = 0
    read_num = 1
    for line in f:
        if (counter - 1) % 4 == 0:
            for chrom in chrom_list:
                new_file.write("{}, {}, {}".format(read_num, chrom, line))
            read_num += 1
        counter += 1
    new_file.close()
    f.close()
    return

# This function is used in MRJob to convert an alignment grid to a string.
# It is used to pass the grid from the mapper to the combiner/reducer so
# that the alignment does not need to be repeated. 
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

# This function is used in MRJob to convert the string back to 
# the grid in the combiner and reducer steps. 
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
            
            allele_re = "[A,G,T,C,N,-]+"
            ref_allele = re.findall(allele_re, cell_info[2])[0]
            read_allele = re.findall(allele_re, cell_info[3])[0]

            prev_re = "\d+|None"
            prev_coor = re.findall(prev_re, cell_info[4])
            if len(prev_coor) == 2:
                # in this case the cell has a traceback 
                prev = (int(prev_coor[0]), int(prev_coor[1]))
            else:
                # in this case the cell has no traceback
                prev = None
            cell = alignment.Cell(val, row, col, ref_allele, read_allele, prev)
            cells.append(cell)
        grid_rows.append(cells)
    grid = np.array(grid_rows, dtype = object)
    return chrom, grid

