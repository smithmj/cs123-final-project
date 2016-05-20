from mrjob.job import MRJob
import alignment

## How do we handle the reference sequence?????
## Can we pass a file into MRJob somehow?
##REF_SEQ = ...

class MRJobAlignment(MRJob):

    def mapper(self, _, line):
        read = line.strip()
        grid = aligment.place_read(REF_SEQ, read, local = True)
        alignment.traceback(grid, REF_SEQ, read, local = True, mr = True, align_dict = None)

    #def combiner(self, ref_seq_pos, nucl):

    def reducer(self, ref_seq_pos, nucls):
        coverage = 0
        counts = {}
        for nucl in nucls:
            coverage += 1
            if nucl not in counts:
                counts[nucl] = 0
            counts[nucl] += 1
        max_val = 0
        most_prob = None
        for nucl in counts:
            if counts[nucl] >= max_val:
                max_val = counts[nucl]
                most_prob = nucl

        if most_prob != REF_SEQ[ref_seq_pos]:
            ## can we change a setting so that we can output an 
            ## SNP object?? Wachs said that the key-value pairs had
            ## to be 'simple' outputs, like strings or numbers...
            yield ref_seq_pos, _____


    def steps(self):
        return [
            MRStep(mapper_pre_filter = "sed -n '2~4p'",
                   mapper = self.mapper,
                   combiner = self.combiner,
                   reducer = self.reducer
                ]