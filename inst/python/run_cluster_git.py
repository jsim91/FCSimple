# check dependencies
import importlib
import subprocess
import sys

REQUIRED_PACKAGES = ['pandas','numpy','git_cluster','os']
for package in REQUIRED_PACKAGES:
    try:
        importlib.import_module(package)
        print(f'{package} is installed')
    except ImportError:
        print(f'{package} not installed. Installing now...')
        if package=='git_cluster':
            subprocess.check_call([sys.executable, "-m", "pip", "install", "git+https://github.com/gaozhangyang/GIT"])
        else:
            subprocess.check_call([sys.executable, "-m", "pip", "install", package])
        print(f'{package} installed successfully')

# finish imports
from git_cluster import GIT
# https://github.com/gaozhangyang/GIT
import numpy as np
import pandas as pd
import os

in_file = sys.argv[1]
out_file = sys.argv[2]
k_value = int(sys.argv[3])

my_data = pd.read_csv(in_file)

try:
    os.remove(in_file)
except OSError:
    pass

git = GIT(k = k_value)
git_fit = git.fit_predict(my_data)
cluster_membership = pd.DataFrame(git_fit)
cluster_membership.to_csv(path_or_buf = out_file + "/__tmp_cl__.csv", index = False)
