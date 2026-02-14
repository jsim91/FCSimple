import sys
import os
try:
    import pandas as pd
    import numpy as np
    import umap
except ImportError as e:
    sys.stderr.write(f"Missing Python dependency: {e}.\n")
    sys.stderr.write("Please run in R: FCSimple::fcs_install_python_dependencies(install=TRUE, precompile=TRUE)\n")
    sys.exit(2)

in_file = sys.argv[1]
out_file = sys.argv[2]
umap_nn = int(sys.argv[3])
mdist = float(sys.argv[4])

data = pd.read_csv(filepath_or_buffer = in_file)
try:
    os.remove(in_file)
except OSError:
    pass

map = umap.UMAP(n_neighbors = umap_nn, init = 'spectral', min_dist = mdist, 
                low_memory = True, n_jobs = 1, verbose = False)
map_output = map.fit_transform(data)
map_df = pd.DataFrame(map_output)
map_df.columns = ["UMAP1","UMAP2"]
map_df.to_csv(out_file + "/__tmp_umap__.csv", index = False)
