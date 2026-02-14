"""
Precompile numba/umap components single-threaded to avoid long JIT compile during
first parallel run. Run once per environment (or after package updates).
"""
import os
import sys

# Force single-threaded numba during compilation
os.environ['NUMBA_NUM_THREADS'] = '1'
print(f"NUMBA_NUM_THREADS set to {os.environ['NUMBA_NUM_THREADS']}")

print('Importing numba...')
import numba
print('Importing pynndescent...')
import pynndescent
print('Importing umap...')
import umap
print('Import complete. Numba cache should be populated.')

# Report actual numba threads
try:
    print('numba.get_num_threads() =', numba.get_num_threads())
except Exception as e:
    print('Could not query numba thread count:', e)

print('Done. You can now run the normal UMAP script with multiple threads.')
