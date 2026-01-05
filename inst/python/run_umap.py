# check dependencies
import importlib
import subprocess
import sys

REQUIRED_PACKAGES = ['pandas','numpy','pynndescent','umap-learn','os']
for package in REQUIRED_PACKAGES:
    try:
        importlib.import_module(package)
        print(f'{package} is installed')
    except ImportError:
        print(f'{package} not installed. Installing now...')
        pkg_spec = package
        installed = False
        try:
            subprocess.check_call([sys.executable, "-m", "pip", "install", pkg_spec])
            installed = True
        except (subprocess.CalledProcessError, FileNotFoundError):
            print("pip install failed; trying pip3...")
            try:
                subprocess.check_call(["pip3", "install", pkg_spec])
                installed = True
            except (subprocess.CalledProcessError, FileNotFoundError) as e:
                print(f"Failed to install {package} with pip and pip3: {e}")

        if installed:
            print(f'{package} installed successfully')
        else:
            print(f'Could not install {package}; continuing without it.')

# finish imports
import umap
import pynndescent
import numpy as np
import pandas as pd
import os

in_file = sys.argv[1]
out_file = sys.argv[2]
umap_nn = sys.argv[3]
mdist = sys.argv[4]

data = pd.read_csv(filepath_or_buffer = in_file)

try:
    os.remove(in_file)
except OSError:
    pass

map = umap.UMAP(n_neighbors = int(umap_nn), init = 'spectral', min_dist = float(mdist), low_memory = True, random_state = 123, transform_seed = 123, verbose = True)
map_output = map.fit_transform(data)
map_df = pd.DataFrame(map_output)
map_df.columns = ["UMAP1","UMAP2"]
map_df.to_csv(out_file + "/__tmp_umap__.csv", index = False)
