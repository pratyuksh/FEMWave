#include "../include/includes.hpp"
#include "../include/core/config.hpp"
#include "../include/core/utilities.hpp"
#include "../include/waveO1/solver.hpp"
#include "../include/waveO1/sparse_grids_handler.hpp"
#include "../include/mymfem/local_mesh_refinement.hpp"

#include <iostream>
#include <Eigen/Core>


void run_waveFG (const nlohmann::json& config,
                 std::string base_mesh_dir,
                 bool load_init_mesh=false)
{
    const int lt = config["level_t"];
    const int lx = config["level_x"];
    std::string sub_mesh_dir = config["mesh_dir"];
    const std::string mesh_dir = base_mesh_dir+sub_mesh_dir;

    auto testCase = make_waveO1_test_case(config);
    WaveO1Solver wave_solverFG(config, testCase,
                               mesh_dir, lx, lt,
                               load_init_mesh);

    auto [ndofs, ht_max, hx_max, errSol] = wave_solverFG();
    std::cout << "\n\nError: "
              << errSol.transpose() << std::endl;
    std::cout << "tMesh size: " << ht_max << std::endl;
    std::cout << "xMesh size: " << hx_max << std::endl;
    std::cout << "#Dofs: " << ndofs << std::endl;
}

void run_waveSG (const nlohmann::json& config,
                 std::string base_mesh_dir,
                 bool load_init_mesh=false)
{
    const int Lx = config["sg_max_level_x"];
    const int L0x = config["sg_min_level_x"];
    const int L0t = config["sg_min_level_t"];

    std::string sub_mesh_dir = config["mesh_dir"];
    const std::string mesh_dir = base_mesh_dir+sub_mesh_dir;

    auto testCase = make_waveO1_test_case(config);
    SparseGridsHandler heat_solverSG(config, testCase,
                                     mesh_dir,
                                     Lx, L0x, L0t,
                                     load_init_mesh);

    auto [ndofs, ht_max, hx_max, errSol] = heat_solverSG();
    std::cout << "\n\nError: "
              << errSol.transpose() << std::endl;
    std::cout << "tMesh size: " << ht_max << std::endl;
    std::cout << "xMesh size: " << hx_max << std::endl;
    std::cout << "#Dofs: " << ndofs << std::endl;
}

void run_waveProjFG (const nlohmann::json& config,
                     std::string base_mesh_dir,
                     bool load_init_mesh=false)
{
    const int lt = config["level_t"];
    const int lx = config["level_x"];
    std::string sub_mesh_dir = config["mesh_dir"];
    const std::string mesh_dir = base_mesh_dir+sub_mesh_dir;

    auto testCase = make_waveO1_test_case(config);
    WaveO1Solver wave_solverFG(config, testCase,
                               mesh_dir, lx, lt,
                               load_init_mesh);

    auto [ndofs, ht_max, hx_max, errSol]
            = wave_solverFG.projection();
    std::cout << "\n\nError: "
              << errSol.transpose() << std::endl;
    std::cout << "tMesh size: " << ht_max << std::endl;
    std::cout << "xMesh size: " << hx_max << std::endl;
    std::cout << "#Dofs: " << ndofs << std::endl;
}

void run_waveProjSG
(const nlohmann::json& config, std::string base_mesh_dir,
 bool load_init_mesh=false)
{
    const int Lx = config["sg_max_level_x"];
    const int L0x = config["sg_min_level_x"];
    const int L0t = config["sg_min_level_t"];

    std::string sub_mesh_dir = config["mesh_dir"];
    const std::string mesh_dir = base_mesh_dir+sub_mesh_dir;

    auto testCase = make_waveO1_test_case(config);
    SparseGridsHandler heat_solverSG(config, testCase,
                                     mesh_dir,
                                     Lx, L0x, L0t,
                                     load_init_mesh);

    auto [ndofs, ht_max, hx_max, errSol]
            = heat_solverSG.compute_projection();
    std::cout << "\n\nError: "
              << errSol.transpose() << std::endl;
    std::cout << "tMesh size: " << ht_max << std::endl;
    std::cout << "xMesh size: " << hx_max << std::endl;
    std::cout << "#Dofs: " << ndofs << std::endl;
}

void run_convergence_waveFG (const nlohmann::json& config,
                             std::string base_mesh_dir,
                             bool load_init_mesh=false)
{
    const int Lt0 = config["min_level_t"];
    const int Lx0 = config["min_level_x"];
    const int Lx = config["max_level_x"];
    std::string sub_mesh_dir = config["mesh_dir"];
    const std::string mesh_dir = base_mesh_dir+sub_mesh_dir;

    int num_levels = Lx-Lx0+1;
    Eigen::VectorXi ndofs(num_levels);
    Eigen::VectorXd h_max(num_levels);
    Eigen::MatrixXd errSol(2,num_levels);
    errSol.setZero();

    auto testCase = make_waveO1_test_case(config);
    WaveO1Solver *solver = nullptr;

    double ht_max, hx_max;
    Eigen::VectorXd errSol_(2);
    for (int k=0; k<num_levels; k++)
    {
        int lt = Lt0+k;
        int lx = Lx0+k;
        solver = new WaveO1Solver(config, testCase,
                                  mesh_dir, lx, lt,
                                  load_init_mesh);

        std::tie(ndofs(k), ht_max, hx_max, errSol_)
                =  (*solver)();
        h_max(k) = std::max(ht_max, hx_max);
        errSol.col(k) = errSol_;
        std::cout << "Level: " << lx
                  << ", ht_max: " << ht_max
                  << ", hx_max: " << hx_max
                  << ", ndofs: " << ndofs(k)
                  << std::endl;
        std::cout << "Error: " << errSol_.transpose()
                  << "\n" << std::endl;
        delete solver;
    }
    std::cout << "\n\nError:\n" << errSol << std::endl;

    write_json_file("waveFG", config, h_max, ndofs, errSol);
}

void run_convergence_waveSG (const nlohmann::json& config,
                             std::string base_mesh_dir,
                             bool load_init_mesh=false)
{
    const int Lx0 = config["min_level_x"];
    const int Lx = config["max_level_x"];
    const int sgL0t = config["sg_min_level_t"];
    const int sgL0x = config["sg_min_level_x"];
    std::string sub_mesh_dir = config["mesh_dir"];
    const std::string mesh_dir = base_mesh_dir+sub_mesh_dir;

    int num_levels = Lx-Lx0+1;
    Eigen::VectorXi ndofs(num_levels);
    Eigen::VectorXd h_max(num_levels);
    Eigen::MatrixXd errSol(2,num_levels);
    errSol.setZero();

    auto testCase = make_waveO1_test_case(config);
    SparseGridsHandler *solver = nullptr;

    double ht_max, hx_max;
    Eigen::VectorXd errSol_(2);
    for (int k=0; k<num_levels; k++)
    {
        int sgLx = Lx0+k;
        solver = new SparseGridsHandler (config, testCase,
                                         mesh_dir,
                                         sgLx, sgL0x, sgL0t,
                                         load_init_mesh);

        std::tie(ndofs(k), ht_max, hx_max, errSol_)
                =  (*solver)();
        h_max(k) = std::max(ht_max, hx_max);
        errSol.col(k) = errSol_;
        std::cout << "Level: " << sgLx
                  << ", ht_max: " << ht_max
                  << ", hx_max: " << hx_max
                  << ", ndofs: " << ndofs(k)
                  << std::endl;
        std::cout << "Error: " << errSol_.transpose()
                  << "\n" << std::endl;
        delete solver;
    }
    std::cout << "\n\nError:\n" << errSol << std::endl;

    write_json_file("waveSG", config, h_max, ndofs, errSol);
}

void run_convergence_waveProjFG
(const nlohmann::json& config, std::string base_mesh_dir,
 bool load_init_mesh=false)
{
    const int Lt0 = config["min_level_t"];
    const int Lx0 = config["min_level_x"];
    const int Lx = config["max_level_x"];
    std::string sub_mesh_dir = config["mesh_dir"];
    const std::string mesh_dir = base_mesh_dir+sub_mesh_dir;

    int num_levels = Lx-Lx0+1;
    Eigen::VectorXi ndofs(num_levels);
    Eigen::VectorXd h_max(num_levels);
    Eigen::MatrixXd errSol(2,num_levels);
    errSol.setZero();

    auto testCase = make_waveO1_test_case(config);
    WaveO1Solver *solver = nullptr;

    double ht_max, hx_max;
    Eigen::VectorXd errSol_(2);
    for (int k=0; k<num_levels; k++)
    {
        int lt = Lt0+k;
        int lx = Lx0+k;
        solver = new WaveO1Solver(config, testCase,
                                  mesh_dir, lx, lt,
                                  load_init_mesh);

        std::tie(ndofs(k), ht_max, hx_max, errSol_)
                =  solver->projection();
        h_max(k) = std::max(ht_max, hx_max);
        errSol.col(k) = errSol_;
        std::cout << "Level: " << lx
                  << ", ht_max: " << ht_max
                  << ", hx_max: " << hx_max
                  << ", ndofs: " << ndofs(k)
                  << std::endl;
        std::cout << "Error: " << errSol_.transpose()
                  << "\n" << std::endl;
        delete solver;
    }
    std::cout << "\n\nError:\n" << errSol << std::endl;

    write_json_file("waveFG", config, h_max, ndofs, errSol);
}

void run_convergence_waveProjSG
(const nlohmann::json& config, std::string base_mesh_dir,
 bool load_init_mesh=false)
{
    const int Lx0 = config["min_level_x"];
    const int Lx = config["max_level_x"];
    const int sgL0t = config["sg_min_level_t"];
    const int sgL0x = config["sg_min_level_x"];
    std::string sub_mesh_dir = config["mesh_dir"];
    const std::string mesh_dir = base_mesh_dir+sub_mesh_dir;

    int num_levels = Lx-Lx0+1;
    Eigen::VectorXi ndofs(num_levels);
    Eigen::VectorXd h_max(num_levels);
    Eigen::MatrixXd errSol(2,num_levels);
    errSol.setZero();

    auto testCase = make_waveO1_test_case(config);
    SparseGridsHandler *solver = nullptr;

    double ht_max, hx_max;
    Eigen::VectorXd errSol_(2);
    for (int k=0; k<num_levels; k++)
    {
        int sgLx = Lx0+k;
        solver = new SparseGridsHandler (config, testCase,
                                         mesh_dir,
                                         sgLx, sgL0x, sgL0t,
                                         load_init_mesh);

        std::tie(ndofs(k), ht_max, hx_max, errSol_)
                =  solver->compute_projection();
        h_max(k) = std::max(ht_max, hx_max);
        errSol.col(k) = errSol_;
        std::cout << "Level: " << sgLx
                  << ", ht_max: " << ht_max
                  << ", hx_max: " << hx_max
                  << ", ndofs: " << ndofs(k)
                  << std::endl;
        std::cout << "Error: " << errSol_.transpose()
                  << "\n" << std::endl;
        delete solver;
    }
    std::cout << "\n\nError:\n" << errSol << std::endl;

    write_json_file("waveSG", config, h_max, ndofs, errSol);
}

void run_localMeshRefinement
(const nlohmann::json& config, std::string base_mesh_dir,
 bool load_init_mesh=false)
{
    const int lx = config["level_x"];
    int deg = config["deg2_x"];

    int lx0 = 0;
    std::string sub_mesh_dir = config["mesh_dir"];
    const std::string mesh_dir = base_mesh_dir+sub_mesh_dir;
    const std::string mesh_file
            = mesh_dir+"/tri_mesh_l"
            +std::to_string(lx0)+".mesh";
    std::cout << "  Initial mesh file: "
              << mesh_file << std::endl;

    bool bool_nonConforming = true;
    if (config.contains("nonConforming")) {
        bool_nonConforming = config["nonConforming"];
    }

    auto xMesh
            = std::make_shared<Mesh>(mesh_file.c_str());

    // geometry
    std::shared_ptr<Polygon> lShaped
            = std::make_shared<LShaped>();
    Array<bool> refineFlags(1);
    Array<double> refineWeights(1);
    refineFlags[0] = true;
    refineWeights[0] = 1-2./3;
    lShaped->set_refine_flags(refineFlags);
    lShaped->set_refine_weights(refineWeights);

    // local refinement
    auto locMeshRef = std::make_unique
            <LocalMeshRefinement>(lShaped);

    double h = 1/std::pow(2, lx);
    locMeshRef->uniform(xMesh, h);
    locMeshRef->local(xMesh, deg, h, bool_nonConforming);

    std::string mesh_name;
    if (bool_nonConforming) {
        mesh_name = "../meshes/lShaped/rg_nc/deg"
                +std::to_string(deg)+"/mesh_l"
                +std::to_string(lx)+".mesh";
    }
    else {
        mesh_name = "../meshes/lShaped/rg_nc/deg"
                +std::to_string(deg)+"/mesh_l"
                +std::to_string(lx)+".mesh";
    }

    std::cout << mesh_name << std::endl;

    std::ofstream mesh_ofs(mesh_name.c_str());
    mesh_ofs.precision(12);
    xMesh->Print(mesh_ofs);
    mesh_ofs.close();
}

int main(int argc, char *argv[])
{   
    // Read config json
    auto config = get_global_config(argc, argv);
    const std::string host = config["host"];
    const std::string run = config["run"];
    std::string base_mesh_dir;

    // check if an initial mesh needs to be loaded
    bool load_init_mesh = false;
    if (config.contains("load_init_mesh")) {
        load_init_mesh = config["load_init_mesh"];
    }

    if (load_init_mesh) {
        base_mesh_dir.assign("../meshes/");
    } else {
        if (host == "local") {
            base_mesh_dir.assign(local_base_mesh_dir);
        }
        else if (host == "cluster") {
            base_mesh_dir.assign(cluster_base_mesh_dir);
        }
    }

    if (run == "simulationFG") {
        run_waveFG(config, base_mesh_dir, load_init_mesh);
    }
    else if (run == "simulationSG") {
        run_waveSG(config, base_mesh_dir, load_init_mesh);
    }
    else if (run == "projectionFG") {
        run_waveProjFG
                (config, base_mesh_dir, load_init_mesh);
    }
    else if (run == "projectionSG") {
        run_waveProjSG
                (config, base_mesh_dir, load_init_mesh);
    }
    else if (run == "convergenceFG") {
        run_convergence_waveFG
                (config, base_mesh_dir, load_init_mesh);
    }
    else if (run == "convergenceSG") {
        run_convergence_waveSG
                (config, base_mesh_dir, load_init_mesh);
    }
    else if (run == "convergenceProjFG") {
        run_convergence_waveProjFG
                (config, base_mesh_dir, load_init_mesh);
    }
    else if (run == "convergenceProjSG") {
        run_convergence_waveProjSG
                (config, base_mesh_dir, load_init_mesh);
    }
    else if (run == "localMeshRefinement") {
        run_localMeshRefinement
                (config, base_mesh_dir, load_init_mesh);
    }

    return 1;
}


// End of file
