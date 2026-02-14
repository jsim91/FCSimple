# check dependencies
import importlib
import subprocess
import sys

REQUIRED_PACKAGES = ['pandas','numpy','opentsne','os']
for package in REQUIRED_PACKAGES:
    try:
        importlib.import_module(package)
        print(f'{package} is installed')
    except ImportError:
        print(f'{package} not installed. Installing now...')
        pkg_spec = package  # install spec (plain package name)
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

try:
    # avoid importing heavy modules at install-time; if missing, instruct user to run installer
    from openTSNE import TSNE
    import numpy as np
    import pandas as pd
    import os
except ImportError as e:
    import sys
    sys.stderr.write(f"Missing Python dependency: {e}.\n")
    sys.stderr.write("Please run in R: FCSimple::fcs_install_python_dependencies(install=TRUE, precompile=TRUE)\n")
    sys.exit(2)

in_file = sys.argv[1]
out_file = sys.argv[2]
n_threads = sys.argv[3]
perpl = sys.argv[4]
seed_arg = sys.argv[5] if len(sys.argv) > 5 else "NA"

data = pd.read_csv(filepath_or_buffer = in_file)

try:
    os.remove(in_file)
except OSError:
    pass

# Set random state only if seed is provided
if seed_arg != "NA":
    random_state = int(seed_arg)
    tsne = TSNE(perplexity = int(perpl), metric = "euclidean", n_jobs = int(n_threads), 
                random_state = random_state, verbose = True)
else:
    # No seed - allows better multi-threading performance
    tsne = TSNE(perplexity = int(perpl), metric = "euclidean", n_jobs = int(n_threads), verbose = True)

map_output = tsne.fit(data.to_numpy())
map_df = pd.DataFrame(map_output)
map_df.columns = ["tSNE1","tSNE2"]
map_df.to_csv(out_file + "/__tmp_tsne__.csv", index = False)
