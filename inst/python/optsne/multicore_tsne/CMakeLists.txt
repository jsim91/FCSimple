# Portable CMake build for opt-SNE's MulticoreTSNE extension.
#
# This replaces the upstream opt-SNE CMakeLists.txt, which assumed Linux/macOS
# (hardcoded g++/g++-8, Unix linker flags, and CMake < 3). This version works
# across Linux, macOS, and Windows (MSVC or MinGW).

cmake_minimum_required(VERSION 3.5 FATAL_ERROR)

project(tsne_multicore CXX)

# OpenMP is optional but strongly recommended for multicore performance.
find_package(OpenMP)

if(OpenMP_CXX_FOUND)
    set(TSNE_OPENMP_FLAGS ${OpenMP_CXX_FLAGS})
else()
    message(WARNING "OpenMP not found; performance will suffer.")
endif()

if(MSVC)
    set(TSNE_WARNING_FLAGS /W3)
    set(TSNE_RELEASE_FLAGS /fp:fast /O2)
else()
    set(TSNE_WARNING_FLAGS -Wall)
    set(TSNE_RELEASE_FLAGS -O3 -fPIC -ffast-math -funroll-loops)
endif()

file(GLOB TSNE_SOURCES *.cpp)

add_library(tsne_multicore MODULE ${TSNE_SOURCES})

# Match the Python wrapper loader (which globs libtsne*.so / tsne*.so /
# *tsne*.dll / *tsne*.dylib) by forcing a predictable "tsne_multicore" name
# without a platform-specific prefix.
set_target_properties(tsne_multicore PROPERTIES PREFIX "")

target_compile_options(tsne_multicore PRIVATE
    ${TSNE_WARNING_FLAGS}
    ${TSNE_OPENMP_FLAGS}
    $<$<CONFIG:Release>:${TSNE_RELEASE_FLAGS}>
)

if(OpenMP_CXX_FOUND)
    target_link_libraries(tsne_multicore PRIVATE OpenMP::OpenMP_CXX)
endif()

if(UNIX)
    target_link_libraries(tsne_multicore PRIVATE m)
endif()