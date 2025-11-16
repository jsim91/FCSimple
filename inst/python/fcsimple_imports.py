# checks if python dependencies for clustering, umap, tsne are available/installable given system env
import importlib
import subprocess
import sys

REQUIRED_PACKAGES = ['pandas','igraph','scipy','os','leidenalg','numpy','pynndescent','umap-learn','opentsne']
for package in REQUIRED_PACKAGES:
    try:
        importlib.import_module(package)
        print(f'{package} is installed')
    except ImportError:
        print(f'{package} not installed. Installing now...')
        subprocess.check_call([sys.executable, "-m", "pip", "install", package])
        print(f'{package} installed successfully')

# import dependencies
from scipy.io import mmread
import igraph as ig
import pandas as pd
import os
import leidenalg as la
import umap
import pynndescent
import numpy as np
from openTSNE import TSNE