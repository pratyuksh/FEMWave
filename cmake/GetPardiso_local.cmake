find_library(PARDISO NAMES libpardiso600-GNU800-X86-64.so
  PATHS
    ${CMAKE_CURRENT_SOURCE_DIR}/third_party/pardiso
)

list(APPEND PARDISO gfortran)
