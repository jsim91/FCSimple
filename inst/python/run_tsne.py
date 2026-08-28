# openTSNE fallback backend for FCSimple::fcs_reduce_dimensions().
# Used automatically when the opt-SNE (MulticoreTSNE) backend is not available.

# check dependencies
import importlib
import subprocess
import sys

REQUIRED_PACKAGES = ['pandas', 'numpy', 'opentsne']
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

try:
    from openTSNE import TSNE
    import numpy as np
    import pandas as pd
    import os
except ImportError as e:
    sys.stderr.write(f"Missing Python dependency: {e}.\n")
    sys.stderr.write("Please run in R: FCSimple::fcs_install_python_dependencies(install=TRUE)\n")
    sys.exit(2)

in_file = sys.argv[1]
out_file = sys.argv[2]
n_threads = int(sys.argv[3])
perpl = float(sys.argv[4])
seed_arg = sys.argv[5] if len(sys.argv) > 5 else "NA"
n_components = int(sys.argv[6]) if len(sys.argv) > 6 else 2

# openTSNE's FFT interpolation only supports 2D; use Barnes-Hut for 3D+
if n_components > 2:
    negative_gradient_method = "bh"
else:
    negative_gradient_method = "fft"

data = pd.read_csv(filepath_or_buffer=in_file)

try:
    os.remove(in_file)
except OSError:
    pass

try:
    if seed_arg != "NA":
        random_state = int(seed_arg)
        tsne = TSNE(n_components=n_components,
                    perplexity=int(perpl), metric="euclidean",
                    n_jobs=n_threads,
                    negative_gradient_method=negative_gradient_method,
                    random_state=random_state, verbose=True)
    else:
        tsne = TSNE(n_components=n_components,
                    perplexity=int(perpl), metric="euclidean",
                    n_jobs=n_threads,
                    negative_gradient_method=negative_gradient_method,
                    verbose=True)

    map_output = tsne.fit(data.to_numpy())
    map_df = pd.DataFrame(map_output)
    map_df.columns = ["tSNE%d" % (i + 1) for i in range(n_components)]
    map_df.to_csv(out_file + "/__tmp_tsne__.csv", index=False)
except Exception as e:
    sys.stderr.write(f"t-SNE failed: {e}\n")
    sys.exit(2)