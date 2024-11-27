from pyigd import IGDReader, IGDTransformer
import argparse
import random
import numpy as np
from tqdm import tqdm
import sys

if __name__ == "__main__":
    parser = argparse.ArgumentParser()
    parser.add_argument("infile", help="The input IGD file")
    parser.add_argument("outfile", help="The output IGD file")
    parser.add_argument("-p", "--zero-to-one", type=float, default=0.001,
        help="The probability of flipping a 0 to a 1")
    parser.add_argument("-n", "--one-to-zero", type=float, default=0.01,
        help="The probability of flipping a 1 to a 0")
    # parser.add_argument("--precise", action="store_true",
    #     help="Use the slower, more precise method of adding noise.")
    args = parser.parse_args()

    with open(args.infile, "rb") as f:
        r = IGDReader(f)
        total_variants = r.num_variants
        total_samples = r.num_samples
    progress_bar = tqdm(total=total_variants)

    # The precise way to add noise, to ensure that the rates are exactly as specified
    class AddNoiseToBV(IGDTransformer):
        def modify_samples(self, position, is_missing, samples):
            with open("edits.csv",'a') as edits:
                global progress_bar 
                progress_bar.update()
                num_samples = len(samples)
                rolled_dice = np.random.uniform(0, 1, size=num_samples)
                for i in range(num_samples):
                    if samples[i] == 1 and rolled_dice[i] <= args.one_to_zero:
                        samples[i] = 0
                        edits.write(f"{position},i,1")
                    elif samples[i] == 0 and rolled_dice[i] <= args.zero_to_one:
                        samples[i] = 1
                        edits.write(f"{position},i,0")
            return samples

    # The approximate way to add noise, much faster but the rates may vary a little more!
    class AddNoiseApprox(IGDTransformer):
        def modify_samples(self, position, is_missing, samples):
            global progress_bar 
            progress_bar.update()
            ones = len(samples)
            zeros = total_samples - ones
            flip01 = int(zeros * args.zero_to_one)
            flip10 = int(ones * args.one_to_zero)
            random.shuffle(samples)
            del samples[ones - flip10:]  # Flip 1's to 0's
            samples.extend(np.random.randint(0, total_samples-1, flip01))
            return sorted(set(samples))


    with open(sys.argv[1], "rb") as fin, open(sys.argv[2], "wb") as fout:
        #if args.precise:
        xformer = AddNoiseToBV(fin, fout, use_bitvectors=True)
        #else:
        #    xformer = AddNoiseApprox(fin, fout, use_bitvectors=False)
        xformer.transform()