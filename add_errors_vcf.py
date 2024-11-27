# adds errors to the vcf file. 0 to 1 flip rate is 0.000379 and 1 to 0 flip rate is 0.063983
# return csv file with three rows: pos, hap, allele
# hap uses a 0 based system. hap 0 and hap 1= individual 1
# allele column is used to represent the original (non-error version) of an allele we flipped

import numpy as np
import sys
import pandas as pd

def flip(allele):
    if allele == 0:
        if np.random.uniform(0, 1) <= 0.000379: #this is per site rate
            return(True)
        else:
            return(False)
    else:
        if np.random.uniform(0, 1) <= 0.063983: 
            return(True)
        else:
            return(False)

def main():
    vcf_file = sys.argv[1]
    noisy_file = "with-error"+sys.argv[1]
    save=[]
    with open(vcf_file, 'r') as f:
        with open(noisy_file, 'w') as nf:
            for line in f:
                if line.startswith('#'):
                    nf.write(line)
                else:
                    toks = line.split("\t")
                    gt_data = toks[9:]
                    new_toks = list(toks[:9])
                    for i,t in enumerate(gt_data):
                        allele1 = int(t[0])
                        allele2 = int(t[2])
                        if flip(allele1):
                            save.append({'pos':toks[1],'hap':2*(i),'allele':allele1})
                            new_allele=1 if allele1==0 else 0
                            repl_str=str(new_allele) + "|"
                        else:
                            repl_str=str(allele1) + "|"

                        if flip(allele2):
                            save.append({'pos':toks[1],'hap':2*(i)+1,'allele':allele2})
                            new_allele2=1 if allele2==0 else 0
                            repl_str+=str(new_allele2)
                        else:
                            repl_str+=str(allele2)

                        new_toks.append(repl_str)
                    new_line = "\t".join(new_toks) + "\n"
                    nf.write(new_line)
    
    df = pd.DataFrame(save)
    df.to_csv("edits"+sys.argv[1]+".csv", index=False)
    return None

if __name__ == '__main__':
    main()