import os
import pandas as pd
import re
import argparse

# Define the column names explicitly for the data lines
columns = [
    "Trial", "status", "num_train_epochs", "learning_rate", "weight_decay",
    "lr_scheduler_type", "warmup_steps", "seed", "train_batch_size",
    "iter", "total_time_s", "eval_loss", "eval_accuracy", "eval_macro_f1", "eval_runtime"
]

# Regular expressions
line_regex = re.compile(r'^(\S+)\s+(.*?)\s+(\S+)$')
split_regex = re.compile(r'\t|\s{2,}')
test_size_regex = re.compile(r'test_size\s*=\s*([\d.]+)')

# Numeric column indices for type conversion in final summary
numeric_col_indices = {
    "num_train_epochs": 2,
    "learning_rate": 3,
    "weight_decay": 4,
    "warmup_steps": 6,
    "seed": 7,
    "train_batch_size": 8,
    "iter": 9,
    "total_time_s": 10,
    "eval_loss": 11,
    "eval_accuracy": 12,
    "eval_macro_f1": 13,
    "eval_runtime": 14
}


def process_file(filepath, max_col="eval_accuracy"):
    # Prepare variables
    max_value = float('-inf')
    max_row_data = None
    test_size_val = None

    # List of encodings to try in order
    encodings_to_try = ["utf-8", "cp437"]

    # We'll store any exception we hit in case both fail
    last_exception = None

    for enc in encodings_to_try:
        try:
            with open(filepath, 'r', encoding=enc) as f:
                # Read the first line to extract test_size
                first_line = f.readline().strip()
                match_test_size = test_size_regex.search(first_line)
                if match_test_size:
                    test_size_val = match_test_size.group(1)

                # Now process the rest of the file for the data rows
                for line in f:
                    match = line_regex.match(line)
                    if match:
                        content = match.group(2)
                        row_data = split_regex.split(content.strip())
                        if len(row_data) == len(columns):
                            try:
                                val = float(row_data[numeric_col_indices[max_col]])
                                if val > max_value:
                                    max_value = val
                                    max_row_data = row_data
                            except ValueError:
                                continue
            # If we got here, reading + parsing succeeded with this encoding
            break

        except UnicodeDecodeError as e:
            # Save the exception, then try the next encoding
            last_exception = e
            continue
        except Exception as e:
            # Some other error (file not found, permission, etc.)
            print(f"Could not read file {filepath}, error: {e}")
            return None
    else:
        # If we exit the for-loop normally (no break),
        # it means both encodings failed with Unicode errors.
        print(f"Could not decode file {filepath} with UTF-8 or CP437. Last error: {last_exception}")
        return None

    # Build a dictionary for the summary row if we found one
    if max_row_data:
        row_dict = {
            "filename": os.path.basename(filepath),
            "test_size": test_size_val
        }
        # Attach all columns from the row
        for col_name, col_value in zip(columns, max_row_data):
            row_dict[col_name] = col_value
        return row_dict
    else:
        return None


def main():
    parser = argparse.ArgumentParser(description='Process all files in a directory and find max rows, compiling a summary.')
    parser.add_argument('--input_dir', type=str, required=True, help='Path to the slurm.out directory')
    parser.add_argument('--max_column', type=str, default='eval_accuracy',
                        help='Column name to select maximum value from')
    parser.add_argument('--output_csv', type=str, default='hyperopt_results.csv',
                        help='Optional path to a single summary CSV for all files')
    args = parser.parse_args()

    # Verify input directory exists
    if not os.path.isdir(args.input_dir):
        print(f"Error: {args.input_dir} is not a directory or does not exist.")
        return

    # Collect all files in the directory (depth=1)
    files = [
        f for f in os.listdir(args.input_dir)
        if os.path.isfile(os.path.join(args.input_dir, f))
    ]

    if not files:
        print(f"No files found in directory: {args.input_dir}")
        return

    summary_rows = []
    for filename in files:
        filepath = os.path.join(args.input_dir, filename)
        # If you only want certain extensions, uncomment:
        if not filename.endswith(".out"):
            continue

        row_dict = process_file(filepath, max_col=args.max_column)
        if row_dict is not None:
            summary_rows.append(row_dict)
            print("="*70)
            print(f"File: {filename}")
            print(f"test_size: {row_dict['test_size']}")
            print(f"Max row for '{args.max_column}':")
            row_series = pd.Series(row_dict)
            print(row_series.to_string())
        else:
            print("="*70)
            print(f"File: {filename} -> No valid rows found or could not be read.")

    # If we have summary rows, optionally compile and output them to a single CSV
    if summary_rows and args.output_csv:
        df_summary = pd.DataFrame(summary_rows)

        # Convert numeric columns
        for col in numeric_col_indices:
            if col in df_summary.columns:
                df_summary[col] = pd.to_numeric(df_summary[col], errors='coerce')

        # Convert test_size if present
        if 'test_size' in df_summary.columns:
            df_summary['test_size'] = pd.to_numeric(df_summary['test_size'], errors='coerce')
            # Sort the rows based on test_size (ascending order by default)
            df_summary.sort_values(by='test_size', inplace=True)

        df_summary.to_csv(args.output_csv, index=False)
        print(f"Summary of all files saved to '{args.output_csv}'")


if __name__ == '__main__':
    main()
