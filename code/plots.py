import matplotlib.pyplot as plt
import numpy as np
    
def plot_histograms(filename):
  f = open(filename,'r')
  data ={}
  chroms = ['1','2','3','4','5','6','7','8','9','10','11','12',
            '13','14','15','16','17','18','19','20','21','22','X','Y']
  for line in f:
    if line[0] != "#":
      line_array = line.split('\t')
      key = line_array[0]
      if key in chroms:
        val = int(line_array[1])
        if key in data:
          data[key].append(val)
        else:
          data[key] = [val]

  #make the plots for the distributions of SNPs over a chromosome:
  for keys in data:
    plt.hist(data[keys])
    plt.title("Chromosome {} SNP distribution".format(keys))
    plt.xlabel("Positions")
    plt.ylabel("Frequency")
    plt.savefig("chrom{}.png".format(keys))
    plt.close()
  
  # make the plot for the number of SNPs on each chromosome 
  data2 = [len(data[chrom]) for chrom in chroms]
  y_pos = np.arange(len(chroms))
  plt.bar(y_pos, data2)
  plt.xticks(y_pos,chroms)
  plt.title("SNPs Across Chromosomes")
  plt.xlabel("Chromosomes")
  plt.ylabel("Number of SNPs")
  plt.savefig("SNP across chroms.png")
  plt.close()

if __name__ == "__main__":
  plot_histograms("SRR710115.flt.hq.vcf")
