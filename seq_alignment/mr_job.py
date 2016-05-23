from mrjob.job import MRJob
import os
import alignment
import file_io

class MRAlignment(MRJob):

    def mapper(self, _, line):
        read_num_index = 0
        chrom_index = 1
        read_index = 2

        pair = line.split(',')
        read_num = pair[read_num_index].strip()
        chrom = pair[chrom_index].strip()
        #print("chrom:", chrom)
        read = pair[read_index].strip()

        # right now this only works on my laptop. MRJob is working in 
        # a temporary directory and we need to figure out how to get
        # the relative or absolute path to the chromosome file
        ref_seq_file = "/home/student/cs123-final-project/seq_alignment/chrom_{}.txt".format(chrom)
        #ref_seq_path = os.path.abspath(ref_seq_file)
        ref_seq = file_io.convert_seq_file_to_str(ref_seq_file)

        grid = alignment.place_read(ref_seq, read, local = True)
        #score = alignment.find_alignment_score(grid)
        #print("score:", score)
        #print(grid)
        alignment_str = file_io.convert_grid_to_str(chrom, grid)
        yield read_num, alignment_str

    def combiner(self, read_num, grid_strs):
        max_score = -float("inf")
        best_grid_str = None
        for string in grid_strs:
            chrom, grid = file_io.convert_str_to_grid(string)
            score = alignment.find_alignment_score(grid)
            if score >= max_score:
                max_score = score
                best_grid_str = string
        yield read_num, best_grid_str

    def reducer_init(self):
        self.align_dict = {}

    def reducer(self, read_num, grid_strs):
        max_score = -float("inf")
        best_grid = None
        best_chrom = None
        for string in grid_strs:
            chrom, grid = file_io.convert_str_to_grid(string)
            score = alignment.find_alignment_score(grid)
            if score >= max_score:
                max_score = score
                best_grid = grid
                best_chrom = chrom
        #print("read {} is aligned to chromosome {}".format(read_num, best_chrom))
        #print("align_dict:", self.align_dict)
        #print("best_grid:", best_grid)
        self.align_dict = alignment.traceback(best_grid, self.align_dict, local = True, chrom = best_chrom)
        
    def reducer_final(self):
        # again, right now this will only work for my laptop 
        alignment.genotype(self.align_dict, save_file = "/home/student/cs123-final-project/seq_alignment/snp.txt")

if __name__ == "__main__":
    MRAlignment.run()
# THINGS WE NEED TO DO

# pre-pair 24 copies of Short reads with each chromosome
# adapt for single job that filters out based on copy of short read
