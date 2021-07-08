# Space-time DG and FEM for the wave equation
This repository consists of numerical solvers for the wave equation based on space-time Discontinuous Galerkin and Finite Element discretisations. The code is written in C++ and it is serial.

## Dependencies
1. MFEM-4.1, you can read more [info](https://mfem.org) and [download](https://mfem.org/download).
Follow the [instructions](https://mfem.org/building/) for the serial-build of MFEM.
It is built with GLVis for visualisation.

1. The linear solver package [PARDISO](https://www.pardiso-project.org/) is used.

1. Eigen library, [info](https://eigen.tuxfamily.org), download version - 3.3.7.

1. For writing the config files [JSON](https://github.com/nlohmann/json), download version - 3.7.0
1. The formatting library [fmt](https://fmt.dev/6.0.0).

1. Unit testing with [GoogleTest](https://github.com/google/googletest) framework.


## Mesh generation
A few sample meshes are given in the folder "meshes". The meshes are generated using the script "run_bisection_mesh_generator_2d" from the repository [PolytopeMeshGen](https://github.com/pratyuksh/PolytopeMeshGen). The meshes generated as such are in the hdf5 format, they are converted to "MFEM mesh v1.0" format (which is readable by the code in the current repo) using the script "run_converter_hdf5_to_mfem_2d.py" and the files in geometry_files folder in the repo [PolytopeMeshGen](https://github.com/pratyuksh/PolytopeMeshGen).