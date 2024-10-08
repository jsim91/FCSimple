from openTSNE import TSNE
import numpy as np
import pandas as pd
import sys
import os

in_file = sys.argv[1]
out_file = sys.argv[2]
n_threads = sys.argv[3]
perpl = sys.argv[4]

data = pd.read_csv(filepath_or_buffer = in_file)

try:
    os.remove(in_file)
except OSError:
    pass

tsne = TSNE(perplexity = int(perpl), metric = "euclidean", n_jobs = int(n_threads), random_state = 123, verbose = True)

map_output = tsne.fit(data.to_numpy())
map_df = pd.DataFrame(map_output)
map_df.columns = ["tSNE1","tSNE2"]
map_df.to_csv(out_file + "/__tmp_tsne__.csv", index = False)
