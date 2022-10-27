import umap
import pynndescent
import numpy as np
import pandas as pd
import sys
import os

in_file = sys.argv[1]
out_file = sys.argv[2]
data = pd.read_csv(filepath_or_buffer = in_file)

try:
    os.remove(in_file)
except OSError:
    pass

map = umap.UMAP(n_neighbors = 30, init = 'spectral', min_dist = 0.1, low_memory = True, random_state = 123, transform_seed = 123, verbose = True)
map_output = map.fit_transform(data)
map_df = pd.DataFrame(map_output)
map_df.columns = ["UMAP1","UMAP2"]
map_df.to_csv(out_file + "/__tmp_umap__.csv", index = False)