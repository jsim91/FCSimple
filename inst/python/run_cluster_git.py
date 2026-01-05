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
        if package == 'git_cluster':
            pkg_spec = "git+https://github.com/gaozhangyang/GIT"
        else:
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
