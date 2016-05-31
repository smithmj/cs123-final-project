from mrjob.job import MRJob
import numpy as np
import re
import alignment
import file_io


class MRAlignment(MRJob):

    def mapper_init(self):
        self.chromosomes = {}
        
        # The following is the full list of chromosomes.
        # We had to scale it down so we only used the first
        # 4 chromosomes.
        # *** #
        # chrom_list = list(range(1,23)) + ["X", "Y"]
        # *** #

        chrom_list = [1, 2, 3, 4]

        for chrom in chrom_list:
            
            chrom_file = "GRCh38.p2_chr{}.fasta".format(chrom)
            chrom_seq = convert_fasta_to_str(chrom_file)

            # for scaling down purposes, we only use the first 
            # 100,000 base pairs of the chromosome sequence
            # *** #
            # self.chromosomes[str(chrom)] = chrom_seq
            # *** #

            self.chromosomes[str(chrom)] = chrom_seq[:10000]

    def mapper(self, _, line):
        read_num_index = 0
        chrom_index = 1
        read_index = 2

        pair = line.split(',')
        read_num = pair[read_num_index].strip()
        chrom = pair[chrom_index].strip()
        read = pair[read_index].strip()
        
        ref_seq = self.chromosomes[chrom]

        grid = place_read(ref_seq, read, local = True)
        
        alignment_str = convert_grid_to_str(chrom, grid)
        yield read_num, alignment_str

    def combiner(self, read_num, grid_strs):
        max_score = -float("inf")
        best_grid_str = None
        for string in grid_strs:
            chrom, grid = convert_str_to_grid(string)
            score = find_alignment_score(grid)
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
            chrom, grid = convert_str_to_grid(string)
            score = find_alignment_score(grid)
            if score >= max_score:
                max_score = score
                best_grid = grid
                best_chrom = chrom
        self.align_dict = traceback(best_grid, self.align_dict, local = True, chrom = best_chrom)
        
    def reducer_final(self):        
        genotype(self.align_dict, save_file = "snps.csv")


if __name__ == "__main__":
    MRAlignment.run()


