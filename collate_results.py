import pandas as pd
import argparse
import multiprocessing
import pickle

def calculate_metrics(df):
    """
    Returns the number of true errors, predicted errors, and overlap
    """
    num_true_errors = 0
    num_predicted_errors = 0
    overlap = 0

    for i,row in df.iterrows():
        true_err = row.error_positions
        predicted = set(row.diff_markers)
        predicted_errors = [int(x) for x in predicted]
        true_errors = [int(x) for x in true_err]

        num_true_errors += len(true_errors)
        num_predicted_errors += len(predicted_errors)

        for err in true_errors:
            if err in predicted_errors:
                overlap += 1
    return num_true_errors,num_predicted_errors,overlap

def collate(errors_csv,index):
    """
    Compare predicted errors with ground truth errors for each haplotype  
    """
    with open(f"collated_{index}.pkl","rb") as f:
        results = pickle.load(f)
    with open(errors_csv) as c:
        edits = pd.read_csv(c)

    results = results.explode('child_name(s)')
    results.rename(columns={'child_name(s)': 'child_sample_name'}, inplace=True)
    grouped_results = results.groupby('child_sample_name')['diff_markers'].agg(lambda x: sum(x, []))
    grouped_results_df = grouped_results.reset_index()
    grouped_results_df['child_sample_name'] = grouped_results_df['child_sample_name'].astype(int)

    start_index = int(results['start_index'].iloc[0])
    end_index = int(results['end_index'].iloc[0])

    within = edits[(edits.igd_index>= start_index) & (edits.igd_index < end_index)].copy()
    within['error_positions'] = within['igd_index']
    grouped_hap = within.groupby('child_sample_name')['error_positions'].apply(list)
    grouped_hap_df = grouped_hap.reset_index()

    final_results = pd.merge(grouped_results_df, grouped_hap_df, on='child_sample_name', how='outer')
    final_results['error_positions'] = final_results['error_positions'].apply(lambda x: x if isinstance(x, list) else [])
    final_results['diff_markers'] = final_results['diff_markers'].apply(lambda x: x if isinstance(x, list) else [])
    
    #with open(f"final_results_{index}.pkl","wb") as f:
        #pickle.dump(final_results,f)

    return calculate_metrics(final_results)
    
def get_input(csv_file):
    """
    Returns iterable for starmap()
    """
    num_windows =  1509

    csv_files = [csv_file for _ in range(num_windows)]
    index = list(range(num_windows))

    return zip(csv_files,index)

def multipool(input):
    number_of_cores = 100

    num_true_errors = 0
    num_predicted_errors = 0
    overlap = 0

    with multiprocessing.Pool(number_of_cores) as pool:
        # distribute computations and collect results:

        results = pool.starmap(collate, input)
        with open("results_per_window_new.txt","w") as r:

            for result in results:
                r.write(f"{result[0]},{result[1]},{result[2]}\n")
                num_true_errors += result[0]

                num_predicted_errors += result[1]
                overlap += result[2]
    
    with open("results_new.txt","w") as f:
        f.write(f"Number of true errors: {num_true_errors}")
        f.write(f"Number of predicted errors: {num_predicted_errors}")
        f.write(f"Number of overlap: {overlap}")


def main():
    parser = argparse.ArgumentParser(description="match admixed haplotypes with reference panel in windows using BK Trees")
    parser.add_argument('-c', '--csv_file', type=str, required=True, help="csv file with errors")
    args = parser.parse_args()

    multipool(get_input(args.csv_file))

if __name__ == "__main__":
    main()
