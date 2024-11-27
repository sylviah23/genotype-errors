import pandas as pd
import BKTree
import argparse
import pyigd
import multiprocessing
import time
import collections
import pickle

def find_diffs(row):
    """
    Finds the differences between the consensus sequence and the query and returns their indices in a list
    """
    diff_markers=[row['start_index']+i for i in range(len(row['query'])) if \
                    row['query'][i] != row['consensus_seq'][i]]
   
    return pd.Series([diff_markers])

def get_consensus_seq(matches):
    """
    Returns the consensus sequence obtained from nearest neighbor matches. ie the string that concats the most
    common character at each position; ties are broken by using the alternate allele ("1")
    """ 
    strings=[''.join([str(e) for e in x]) for x in matches]
    consensus = ''
    for i in range(len(strings[0])):
        column = [s[i] for s in strings]
        counter = collections.Counter(column)
        most_common = counter.most_common(2)
        if len(most_common) == 1:
          consensus += most_common[0][0]
        elif most_common[1][1] == most_common[0][1]:
          consensus += "1"
        else:
          consensus += most_common[0][0]
        
    return consensus

def hamming_distance(vector1, vector2):
    """
    Hamming distance between two lists
    """
    return sum(el1 != el2 for el1, el2 in zip(vector1, vector2))

def write(results,index):
    """
    Write the results into a pickle file
    """
    with open(f"collated_{index}.pkl","wb") as f:
        pickle.dump(results,f)

def collate(matches,index):
    """
    Create consensus sequence and find predicted error positions
    """
    results = pd.DataFrame.from_dict(matches)
    results['consensus_seq'] = results['matches'].apply(lambda x: pd.Series(get_consensus_seq(x)))
    results['diff_markers'] = results.apply(lambda x: find_diffs(x),axis=1)
    results['neighborhood_size'] = results['matches'].apply(len)
    write(results,index)

def match(haplotypes,start_index,end_index,samples,relatives,trios,child_find,index):
    """
    Find nearest neighbors matches for each haplotype
    """
    with open(samples,"r") as s: #samples = ukbiobank_header.txt (all sample names of individuals in biobank data)
        samp_names = s.readline().split()

    with open(relatives,"rb") as r: #relatives = pickle file of dictionary {sample name: close relative sample names to skip when querying the key}
        skipping_dict = pickle.load(r)

    with open(trios,"rb") as t: #trios = pickle file with a set of all the trios indices in the igd
        trio_set = pickle.load(t)
    
    with open(child_find,"rb") as c: #child_find = pickle file of dictionary {trio sample names: associated child sample name(s) in a list} if trio sample is a child then key = value
        child_finder = pickle.load(c)

    # Initialize an empty BKTreeNode
    root_node = BKTree.BKTreeNode.make_empty()

    #Insert reference vectors into the BK tree
    for i,genotype in enumerate(haplotypes):
        BKTree.bk_tree_insert(root_node, [samp_names[i // 2]], genotype, hamming_distance)

    #Look up nearest neighbors for each ADMIXgenotype
    matches=[]
    for i,query_vector in enumerate(haplotypes):
        if (i // 2) in trio_set:
            results, dist_best, exact_matches=BKTree.bk_tree_lookup(root_node, query_vector, hamming_distance, skipping_dict[samp_names[i//2]])
            matches.append({\
                    'start_index':start_index,\
                    'end_index':end_index,\
                    'query':''.join([str(x) for x in query_vector]),\
                    'matches':[x.vector for x in results],\
                    'edit_distance':dist_best, \
                    'exact_matches':[x.vector for x in exact_matches],\
                    'hap_index': i,\
                    'sample_name':samp_names[i//2],\
                    'child_name(s)': child_finder[samp_names[i//2]]\
                    })
            
    return collate(matches,index)

def window(igd_file,start_index,end_index,samples, relatives,trios,child_find,index):
    """
    Create a list of haplotypes for a particular window starting at variant start_index and ending at variant
    end_index
    """
    start_time = time.time()
    print(f"starting window {index}",flush=True)

    with pyigd.IGDFile(igd_file) as igd:
        haplotypes=[[0 for _ in range(end_index-start_index)] for _ in range(igd.num_samples)]

        for i in range(end_index-start_index): #it is important not to read the end marker itself to avoid overlaps
            for variant in igd.get_samples(i+start_index)[2]:
                haplotypes[variant][i] = 1

    match(haplotypes,start_index,end_index,samples,relatives,trios,child_find,index)
    end_time=time.time()
    print(f"time to parse igd and match a window {index}: {(end_time)-(start_time)} seconds",flush=True)

def get_input(igd_file,window_size,samples,relatives,trios,child_find):
    """
    Returns iterable used for starmap()
    """

    with pyigd.IGDFile(igd_file) as f:
      df_positions = [0]

      counter = window_size+f.get_position_and_flags(0)[0]
      for i in range(1,f.num_variants-1):
          if (f.get_position_and_flags(i+1)[0] < counter):
              continue
          else:
              df_positions.append(i+1)
              df_positions.append(i+1)
              counter = f.get_position_and_flags(i+1)[0]+window_size
            

    num_windows = len(df_positions)
    positions_iter = iter(df_positions)
    igd_files = [igd_file for _ in range(num_windows)]
    samples_files = [samples for _ in range(num_windows)]
    relatives_files = [relatives for _ in range(num_windows)]
    trios_files = [trios for _ in range(num_windows)]
    child_find_files = [child_find for _ in range(num_windows)]
    index = list(range(num_windows))

    return zip(igd_files,positions_iter,positions_iter,samples_files,relatives_files,trios_files,child_find_files,index)

def multipool(input):
    number_of_cores = 64
    with multiprocessing.Pool(number_of_cores) as pool:
        # distribute computations and collect results:
        pool.starmap(window, input)


def main():
    parser = argparse.ArgumentParser(description="match admixed haplotypes with reference panel in windows using BK Trees")
    parser.add_argument('-w', '--window_size', type=int, required=True, help="window size to run matching with")
    parser.add_argument('-i', '--igd_file', required=True, help="igd file of data")
    parser.add_argument('-s','--samples',required=True,help="txt file with a list of all sample names in the igd file")
    parser.add_argument('-r','--relatives',required=True,help="pickle file of dictionary with keys as samples, values as the close relatives to skip when querying the key")
    parser.add_argument('-t','--trios',required=True,help="pickle file with a set of all the trios indices in the igd")
    parser.add_argument('-d','--child_find',required=True,help="pickle file of dictionary where keys are trio samples and values are the associated child")
    args = parser.parse_args()

    multipool(get_input(args.igd_file,args.window_size,args.samples,args.relatives,args.trios,args.child_find))

if __name__ == "__main__":
    main()
