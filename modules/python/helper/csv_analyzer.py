import pandas as pd
import os
import sys
import h5py
from tqdm import tqdm

csv_path = sys.argv[1]

tmp_df = pd.read_csv(csv_path, header=None)
img_files = list(tmp_df[0])
allele_path = list(tmp_df[2])
for i in tqdm(range(len(img_files))):
    hdf5_file_path = img_files[i]
    allele_dict_path = allele_path[i]
    if os.path.isfile(hdf5_file_path) is False:
        print("INVALID FILE PATH: ", hdf5_file_path)
        exit()
    elif os.path.isfile(allele_dict_path) is False:
        print("INVALID FILE PATH: ", allele_dict_path)
        exit()
    else:
        hdf5_file = h5py.File(hdf5_file_path, 'r')
        if 'images' not in hdf5_file.keys():
            print("NO IMAGES IN HDF5", hdf5_file_path)
            exit()
