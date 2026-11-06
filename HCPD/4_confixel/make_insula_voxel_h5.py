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
    h5f = h5py.File(h5_file, "r")
    data = h5f["insula_pearsonr"][row_num][:]
    h5f.close()
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
    print("Extracting connectivity data...")
    for _, row in tqdm(cohort_df.iterrows(), total=cohort_df.shape[0]):
        scalar_file = op.join(relative_root, row["source_file"])
        scalar_data = get_row(insula_voxel_row_num, scalar_file)
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
    f = h5py.File(output_file, "w")

    for scalar_name in scalars.keys():
        one_scalar_h5 = f.create_dataset(
            "scalars/{}/values".format(scalar_name),
            data=np.vstack(scalars[scalar_name]),
        )
        one_scalar_h5.attrs["column_names"] = list(
            sources_lists[scalar_name]
        )  # column names: list of source .mif filenames
    f.close()
    return int(not op.exists(output_file))


def get_parser():

    parser = argparse.ArgumentParser(description="Create a hdf5 file of volume data")
    parser.add_argument(
        "--cohort-file",
        "--cohort_file",
        help="Path to a csv with demographic info and paths to data.",
        required=True,
    )
    parser.add_argument(
        "--relative-root",
        "--relative_root",
        help="Root to which all paths are relative, i.e. defining the "
        "(absolute) path to root directory of group_mask_file, cohort_file, and output_hdf5.",
        type=op.abspath,
        default="/inputs/",
    )
    parser.add_argument(
        "--output-hdf5",
        "--output_hdf5",
        help="Name of HDF5 (.h5) file where outputs will be saved.",
        default="fixelarray.h5",
    )
    parser.add_argument(
        "--insula_voxel_row",
        "--insula_voxel_row",
        help="Integer representing the row to create the h5 file for.",
    )
    return parser


def main():
    parser = get_parser()
    args = parser.parse_args()
    status = write_hdf5(
        insula_voxel_row_num=args.insula_voxel_row,
        cohort_file=args.cohort_file,
        output_h5=args.output_hdf5,
        relative_root=args.relative_root,
    )

    return status


if __name__ == "__main__":
    sys.exit(main())

