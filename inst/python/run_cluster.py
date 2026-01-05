# check dependencies
import importlib
import subprocess
import sys

REQUIRED_PACKAGES = ['pandas','igraph','scipy','os','leidenalg']
for package in REQUIRED_PACKAGES:
    try:
        importlib.import_module(package)
        print(f'{package} is installed')
    except ImportError:
        print(f'{package} not installed. Installing now...')
        try:
            subprocess.check_call([sys.executable, "-m", "pip", "install", package])
        except (subprocess.CalledProcessError, FileNotFoundError):
            print("pip install failed; trying pip3...")
            try:
                subprocess.check_call(["pip3", "install", package])
            except (subprocess.CalledProcessError, FileNotFoundError) as e:
                print(f"Failed to install {package} with pip and pip3: {e}")
                continue
        print(f'{package} installed successfully')

# import dependencies
from scipy.io import mmread
import igraph as ig
import pandas as pd
import os

in_file = sys.argv[1]
out_file = sys.argv[2]
algo = sys.argv[3]
cl_res = float(sys.argv[4])

adjm = mmread(in_file)

try:
    os.remove(in_file)
except OSError:
    pass

# // create graph from adj matrix
sources, targets = adjm.nonzero()
G = ig.Graph(directed=False)
G.add_vertices(adjm.shape[0])
edges = list(zip(sources, targets))
G.add_edges(edges)

if algo=='leiden':
    import leidenalg as la
    # // cluster with leiden: Traag, V.A., Waltman. L., Van Eck, N.-J. (2018). From Louvain to Leiden: guaranteeing well-connected communities
    partition = la.find_partition(G, la.RBConfigurationVertexPartition, resolution_parameter = cl_res, seed = 123)
elif algo=='louvain':
    # // cluster with louvain: VD Blondel, J-L Guillaume, R Lambiotte and E Lefebvre: Fast unfolding of community hierarchies in large networks, J Stat Mech P10008 (2008)
    partition = G.community_multilevel()
cluster_membership = pd.DataFrame(partition.membership)
# // write results
cluster_membership.to_csv(path_or_buf = out_file + "/__tmp_cl__.csv", index = False)
