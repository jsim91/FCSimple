# Lightweight readiness probe for the opt-SNE backend.
#
# Verifies that the MulticoreTSNE package (and its compiled C++ shared
# library) can be imported under the same Python interpreter that executes
# this script. Prints a human-readable status and exits:
#   0 = ready
#   2 = not importable (missing package or failed to build/load native lib)

import sys

try:
    from MulticoreTSNE import MulticoreTSNE as TSNE  # noqa: F401
except ImportError as e:
    sys.stderr.write(f"opt-SNE not importable: {e}\n")
    sys.stderr.write("Build/install it with: FCSimple::fcs_install_python_dependencies(install=TRUE, build_optsne=TRUE)\n")
    sys.exit(2)

print("opt-SNE (MulticoreTSNE) is importable.")
sys.exit(0)