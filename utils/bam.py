#!/usr/bin/env python

## Convert BAM Rh files into SXS format
import argparse
import os
import glob

import h5py
import numpy as np


parser = argparse.ArgumentParser('bam')
parser.add_argument("--run-directory", type=str, help="BAM run directory")
parser.add_argument("--output-directory", type=str, help="Output directory", default=".")
parser.add_argument("--extraction-radius", type=float, help="Extraction radius", default=1100.0)
parser.add_argument("--level", type=str, help="AMR level", default=1)
args = parser.parse_args()

name = os.path.basename(args.run_directory)

f = h5py.File(f"{args.output_directory}/rh_{name}.h5", "w")
f.create_group("Extrapolated_N2.dir")

pattern = f"{args.run_directory}/postprocessed/Rh_*_r{args.extraction_radius}_l{args.level}.txt"
for filename in sorted(glob.glob(pattern)):
    basename = os.path.basename(filename)
    sp =  basename.split("_")
    # print(sp)
    l = int(sp[1].removeprefix("l"))
    m_sign = 1
    m_str = sp[2].removeprefix("m")
    if m_str.startswith("m"):
        m_sign = -1
        m_str = m_str.removeprefix("m")
    m = m_sign * int(m_str)

    rh_data = np.loadtxt(filename)
    ds = f.create_dataset(f"Extrapolated_N2.dir/Y_l{l}_m{m}.dat", (len(rh_data),3), dtype='<f8')
    ds[:,:] = rh_data[:,:3]

f.close()
