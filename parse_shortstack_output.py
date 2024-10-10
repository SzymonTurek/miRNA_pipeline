import argparse
import pandas as pd


parser = argparse.ArgumentParser(description="Read a CSV file")
parser.add_argument('results_file', type=str, help='Path to the Results.txt file')
parser.add_argument('counts_file', type=str, help='Path to the counts.txt file')
parser.add_argument('sample_info_file', type=str, help='Path to the sample_info.csv file')
args = parser.parse_args()
    


def load_csv_file(shortstack_output_file, separator):  
    
    # Read the CSV file using pandas
    try:
        df = pd.read_csv(shortstack_output_file, sep=separator)
        print(f"Successfully loaded {shortstack_output_file}")
      #  print(df.head())  # Display the first few rows
    except FileNotFoundError:
        print(f"Error: File '{shortstack_output_file}' not found.")
    except pd.errors.EmptyDataError:
        print(f"Error: File '{shortstack_output_file}' is empty.")
    except Exception as e:
        print(f"Error: {e}")
    return df


def main():
    results = load_csv_file(args.results_file, "\t")
    counts = load_csv_file(args.counts_file, "\t")
    sample_info = load_csv_file(args.sample_info_file, ",")
    full_counts = counts

    counts = counts.iloc[:, 3:]
    counts.insert(0, "known_miRNAs" ,results["known_miRNAs"])
    counts.insert(0, "MIRNA" ,results["MIRNA"])
    Y_counts = counts[counts['MIRNA'] == "Y"]
    Y_counts = Y_counts.drop('MIRNA', axis=1)

    sample_info['read1'] = sample_info['read1'].str.replace(r'input_fastq\/', '', regex=True)
    sample_info['read1'] = sample_info['read1'].str.replace(r'.fastq.gz', '', regex=True)
    sample_info['read1'] = sample_info['read1'].str.replace(r'.fastq', '', regex=True)
    sample_info['read1'] = sample_info['read1'].str.replace(r'.fq', '', regex=True)
    sample_info['read1'] = sample_info['read1'].str.replace(r'.fq.gz', '', regex=True)
    sample_info['read1'] = sample_info['read1'].str.replace(r'trimmed', 'trimmed_condensed', regex=True)

    my_dict = dict(zip(sample_info['read1'], sample_info['sample_id']))
    counts.rename(columns=my_dict, inplace=True)
    Y_counts.rename(columns=my_dict, inplace=True)
    full_counts.rename(columns=my_dict, inplace=True)
   
    Y_counts.to_csv('named_Y_counts.csv', index=False)
    #counts.to_csv('named_counts.csv', index=False)
    full_counts.to_csv('Counts_with_names.csv', index=False)


if __name__ == "__main__":
    main()