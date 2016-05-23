from mrjob.job import MRJob
import alignment
import file_io

class MRJobAlignment(MRJob):


    def mapper(self, _, line):
        read_num_index = 0
        chrom_index = 1
        read_index = 2

        pair = line.split(',')
        read_num = pair[read_num_index].strip()
        chrom = pair[chrom_index].strip()
        read = pair[read_index].strip()

        ref_seq_file = "chrom_{}.txt".format(chrom)
        ref_seq = file_io.convert_seq_file_to_str(ref_seq_file)

        grid = aligment.place_read(ref_seq, read, local = True)
        #alignment.traceback(grid, REF_SEQ, read, local = True, mr = True, align_dict = None)
        alignment_str = convert_grid_to_string(chrom, grid)
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

    def reducer(self, ref_seq_pos, nucls):
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
        # call traceback and genotyping
        # yeild SNPs

    def steps(self):
        return [
            MRStep(mapper_pre_filter = "sed -n '2~4p'",
                   mapper = self.mapper,
                   combiner = self.combiner,
                   reducer = self.reducer
                ]

# THINGS WE NEED TO DO

# pre-pair 24 copies of Short reads with each chromosome
# adapt for single job that filters out based on copy of short read
