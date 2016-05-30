import matplotlib.pyplot as plt
import numpy as np

#def plot_chrom(chrom,f):
#  data = []
#  for line in f:
#    line_array = line.split()
#    if line_array[0][0] != chrom:
#      continue
    


if __name__ == "__main__":
  f = open("SRR710115.flt.hq.vcf",'r')
  data ={}
  data2 =[]
  chroms = ['1','2','3','4','5','6','7','8','9','10','11','12',
            '13','14','15','16','17','18','19','20','21','22','X','Y']
  for line in f:
    line_array = line.split()
    key = line_array[0]
    if key in chroms:
      val = int(line_array[1])
      if key in data:
        data[key].append(val)
      else:
        data[key] = [val]
    else:
      continue
  for keys in data:
    plt.hist(data[keys])
    plt.title("Chromosome {} SNP distribution".format(keys))
    plt.xlabel("Positions")
    plt.ylabel("Frequency")
    plt.savefig("chrom{}.png".format(keys))
    plt.close()
    data2.append(len(data[keys]))
  y_pos = np.arange(len(chroms))
  plt.bar(y_pos,data2)
  plt.xticks(y_pos,chroms)
  plt.title("SNPs Across Chromosomes")
  plt.xlabel("Chromosomes")
  plt.ylabel("Number of SNPs")
  plt.savefig("SNP across chroms.png")
  plt.close()
