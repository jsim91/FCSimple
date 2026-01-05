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
import leidenalg as la
import umap
import pynndescent
import numpy as np

from openTSNE import TSNE
