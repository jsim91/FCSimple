# python E:\sample_FCS\run_tsne.py
from openTSNE import TSNE
import numpy as np
import pandas as pd
# from matplotlib import pyplot as plt
import sys
import os

in_file = sys.argv[1]
out_file = sys.argv[2]
data = pd.read_csv(filepath_or_buffer = in_file)

try:
    os.remove(in_file)
except OSError:
    pass

#tsne = TSNE(
#    perplexity=30,
#    metric="euclidean",
#    # n_jobs=8,
#    random_state=123,
#    verbose=True,
#)

#map_output = tsne.fit(data)
#map_df = pd.DataFrame(map)
#map_df.columns = ["tSNE1","tSNE2"]
#print(embedding[:10])

#plt.plot(embedding[:,0],embedding[:,1],'ro')
#plt.show()

#map_df.to_csv(out_file + "/__tmp_tsne__.csv", index = False)