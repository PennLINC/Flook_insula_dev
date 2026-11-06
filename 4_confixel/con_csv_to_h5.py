import sys
import h5py
import numpy as np


in_file = sys.argv[1]
out_file = in_file.replace(".csv", ".h5")
print(f"converting {in_file} to {out_file}")
xx = np.genfromtxt(in_file, delimiter=",", skip_header=1)
h5out = h5py.File(out_file, "w")
ds = h5out.create_dataset("insula_pearsonr", data=xx.T)
h5out.close()

