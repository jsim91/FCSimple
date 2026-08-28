# run opt-SNE (forked from DmitryUlyanov/Multicore-TSNE by omiq-ai/Multicore-opt-SNE)
# CSV -> opt-SNE -> CSV command line helper for FCSimple::fcs_reduce_dimensions().

import importlib
import subprocess
import sys

REQUIRED_PACKAGES = ['numpy', 'cffi']
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
    from MulticoreTSNE import MulticoreTSNE as TSNE
    import numpy as np
    import pandas as pd
    import os
except ImportError as e:
    sys.stderr.write(f"Missing Python dependency: {e}.\n")
    sys.stderr.write("Please run in R: FCSimple::fcs_install_python_dependencies(install=TRUE, build_optsne=TRUE)\n")
    sys.exit(2)

in_file = sys.argv[1]
out_file = sys.argv[2]
n_threads = int(sys.argv[3])
perpl = float(sys.argv[4])
seed_arg = sys.argv[5] if len(sys.argv) > 5 else "NA"
n_components = int(sys.argv[6]) if len(sys.argv) > 6 else 2

data = pd.read_csv(filepath_or_buffer=in_file)

try:
    os.remove(in_file)
except OSError:
    pass

try:
    # -1 signals "no fixed seed" (upstream MulticoreTSNE convention).
    random_state = -1 if seed_arg == "NA" else int(seed_arg)

    early_exaggeration = 12
    # opt-SNE mode requires the optimal initial learning rate to be
    # N / early_exaggeration. The upstream run_optsne.py computes this
    # explicitly (the C++ core only warns, it does not override).
    learning_rate = data.shape[0] / early_exaggeration

    # auto_iter=True enables opt-SNE mode: the algorithm auto-selects the
    # early-exaggeration iterations and total iterations by monitoring
    # KL-divergence rate of change. auto_iter_end and early_exaggeration
    # mirror the upstream run_optsne.py defaults.
    tsne = TSNE(n_components=n_components,
                perplexity=int(perpl),
                metric='euclidean',
                n_jobs=n_threads,
                angle=0.5,
                learning_rate=learning_rate,
                auto_iter=True,
                auto_iter_end=5000,
                early_exaggeration=early_exaggeration,
                random_state=random_state,
                verbose=True)

    map_output = tsne.fit_transform(data.to_numpy())
    map_df = pd.DataFrame(map_output)
    map_df.columns = ["tSNE%d" % (i + 1) for i in range(n_components)]
    map_df.to_csv(out_file + "/__tmp_tsne__.csv", index=False)
except Exception as e:
    sys.stderr.write(f"opt-SNE failed: {e}\n")
    sys.exit(2)