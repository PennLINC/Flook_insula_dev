import argparse
import os
import os.path as op
import sys
from collections import defaultdict

import h5py
import numpy as np
import pandas as pd
from tqdm import tqdm

def get_row(row_num, h5_file):
    with h5py.File(h5_file, "r") as h5f:
        
        if "insula_pearsonr" in h5f and isinstance(h5f["insula_pearsonr"], h5py.Dataset):
            dataset = h5f["insula_pearsonr"]
            # Access as a simple array
            data = dataset[row_num]
        else:
            raise ValueError("The key 'insula_pearsonr' is not a dataset or does not exist.")
    
    return data

def write_hdf5(
    insula_voxel_row_num, cohort_file, output_h5="voxeldb.h5", relative_root="/"
):
    """
    Load all volume data.
    Parameters
    -----------
    insula_voxel_row_num: int
        The row from which to grab the whole-brain connectivity data.
    cohort_file: str
        path to a csv with demographic info and paths to data
    output_h5: str
        path to a new .h5 file to be written
    relative_root: str
        path to which group_mask_file and cohort_file (and its contents) are relative
    """
    # gather cohort data
    cohort_df = pd.read_csv(op.join(relative_root, cohort_file))

    # upload each cohort's data
    scalars = defaultdict(list)
    sources_lists = defaultdict(list)
    print(f"Extracting connectivity data for row {insula_voxel_row_num}...")
    for _, row in tqdm(cohort_df.iterrows(), total=cohort_df.shape[0]):
        scalar_file = op.join(relative_root, row["source_file"])
        scalar_data = get_row(insula_voxel_row_num, scalar_file)  # Use row_num directly
        # append to specific scalar_name
        scalars[row["scalar_name"]].append(scalar_data)
        # append source mif filename to specific scalar_name
        sources_lists[row["scalar_name"]].append(row["source_file"])

    # Write the output:
    output_file = op.join(relative_root, output_h5)
    # make dir if not exist:
    output_dir = op.dirname(output_file)
    if not op.exists(output_dir):
        os.makedirs(output_dir, exist_ok=True)
    # initialize HDF5 file:
    with h5py.File(output_file, "w") as f:
        for scalar_name, data in scalars.items():
            one_scalar_h5 = f.create_dataset(
                f"scalars/{scalar_name}/values",
                data=np.vstack(data),
            )
            one_scalar_h5.attrs["column_names"] = sources_lists[scalar_name]  # column names: list of source .mif filenames

    return int(not op.exists(output_file))

def get_parser():
    parser = argparse.ArgumentParser(description="Create multiple HDF5 files of volume data")
    parser.add_argument(
        "--cohort_file",
        help="Path to a csv with demographic info and paths to data.",
        required=True,
    )
    parser.add_argument(
        "--relative_root",
        help="Root to which all paths are relative, i.e. defining the "
        "(absolute) path to root directory of group_mask_file, cohort_file, and output_hdf5.",
        type=op.abspath,
        default="/inputs/",
    )
    parser.add_argument(
        "--output_prefix",
        help="Prefix for HDF5 (.h5) files where outputs will be saved.",
        default="voxel",
    )
    parser.add_argument(
        "--start_row",
        help="Starting integer for row index (0-based).",
        type=int,
        default=0,
    )
    parser.add_argument(
        "--end_row",
        help="Ending integer for row index (0-based).",
        type=int,
        required=True,
    )
    return parser

def main():
    parser = get_parser()
    args = parser.parse_args()
    
    # Loop through the range of rows and create HDF5 files
    for row_num in range(args.start_row, args.end_row + 1):
        output_hdf5 = f"{args.output_prefix}{row_num + 1}.h5"  # Adjust file naming for 1-based output
        status = write_hdf5(
            insula_voxel_row_num=row_num,
            cohort_file=args.cohort_file,
            output_h5=output_hdf5,
            relative_root=args.relative_root,
        )
        if status:
            print(f"Failed to create {output_hdf5}")
        else:
            print(f"Successfully created {output_hdf5}")

    return 0

if __name__ == "__main__":
    sys.exit(main())
