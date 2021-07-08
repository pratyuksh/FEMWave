#include <gtest/gtest.h>

#include "mfem.hpp"
using namespace mfem;

#include <iostream>
#include <fstream>
#include <chrono>
#include <random>

#include <Eigen/Core>
#include <Eigen/Dense>

#include "../include/core/config.hpp"
#include "../include/mymfem/utilities.hpp"

using namespace std;


//! Test the point location algorithm implemented
//! in class PointLocator for meshes with local refinement
TEST(MfemUtil, pointLocator1)
{
    std::string input_dir
            = "../input/gammaShaped_bisecRefine/";

    int lx1 = 1;
    const std::string mesh_file1
            = input_dir+"mesh_lx"+to_string(lx1);

    int lx2 = 6;
    const std::string mesh_file2
            = input_dir+"mesh_lx"+to_string(lx2);

    //std::cout << mesh_file1 << std::endl;
    //std::cout << mesh_file2 << std::endl;

    Mesh mesh1(mesh_file1.c_str());
    Mesh mesh2(mesh_file2.c_str());

    const IntegrationRule *ir = nullptr;
    ir = &IntRules.Get(2, 3);

    std::random_device rd;     // initialise (seed) engine
    std::mt19937 rng(rd());    // random-number engine used
    std::uniform_int_distribution<int> uni(0,mesh1.GetNE()-1);

    auto elId1 = uni(rng);
    //std::cout << "Generate in element:"
    //          << elId1 << std::endl;
    //elId1 = 3;
    ElementTransformation *trans1
            = mesh1.GetElementTransformation(elId1);
    DenseMatrix true_x(mesh1.Dimension(), ir->GetNPoints());
    Vector xk;
    for (int k=0; k < true_x.NumCols(); k++) {
        true_x.GetColumnReference(k,xk);
        trans1->Transform(ir->IntPoint(k), xk);
    }

    Array <int> elIds1(true_x.NumCols());
    Array <IntegrationPoint> ips1(true_x.NumCols());

    Array <int> elIds2(true_x.NumCols());
    Array <IntegrationPoint> ips2(true_x.NumCols());

    auto start1 = std::chrono::high_resolution_clock::now();
    mesh2.FindPoints(true_x, elIds1, ips1);
    auto end1 = std::chrono::high_resolution_clock::now();

    auto start2 = std::chrono::high_resolution_clock::now();
    PointLocator point_locator(&mesh2);
    int init_elId = 0;
    for (int k=0; k < true_x.NumCols(); k++) {
        Vector xk;
        true_x.GetColumn(k, xk);
        std::tie (elIds2[k],ips2[k])
                = point_locator(xk, init_elId);
        init_elId = elIds2[k];
    }
    auto end2 = std::chrono::high_resolution_clock::now();

    ElementTransformation *trans2 = nullptr;
    Vector x1, x2, x;
    double TOL = 1E-8;
    for (int k=0; k < true_x.NumCols(); k++)
    {
        trans2 = mesh2.GetElementTransformation(elIds1[k]);
        trans2->Transform(ips1[k], x1);

        trans2 = mesh2.GetElementTransformation(elIds2[k]);
        trans2->Transform(ips2[k], x2);

        true_x.GetColumn(k,x);

        Vector xmx1(x.Size()), xmx2(x.Size());
        subtract(x, x1, xmx1);
        subtract(x, x2, xmx2);

        ASSERT_LE(xmx1.Norml1(), TOL);
        ASSERT_LE(xmx2.Norml1(), TOL);
    }

    auto duration1 = std::chrono::duration_cast
            <std::chrono::microseconds>(end1 - start1);
    auto duration2 = std::chrono::duration_cast
            <std::chrono::microseconds>(end2 - start2);

    std::cout << "Run times: "
              << "\t slow  " << duration1.count()
              << "\t fast  " << duration2.count()
              << std::endl;
}

//! Test the point location algorithm implemented
//! in class PointLocator for quasi-uniform meshes,
//! generated by _serial_ code.
//! The shared vertices table is generated,
//! but expected to be empty.
TEST(MfemUtil, pointLocator2)
{
    std::string input_dir
            = "../input/unitSquare_quasiUniform/";

    int lx1 = 1;
    const std::string mesh_file1
            = input_dir+"mesh_lx"+to_string(lx1);

    int lx2 = 6;
    const std::string mesh_file2
            = input_dir+"mesh_lx"+to_string(lx2);

    //std::cout << mesh_file1 << std::endl;
    //std::cout << mesh_file2 << std::endl;

    Mesh mesh1(mesh_file1.c_str());
    Mesh mesh2(mesh_file2.c_str());

    const IntegrationRule *ir = nullptr;
    ir = &IntRules.Get(2, 3);

    std::random_device rd;     // initialise (seed) engine
    std::mt19937 rng(rd());    // random-number engine used
    std::uniform_int_distribution<int> uni(0,mesh1.GetNE()-1);

    auto elId1 = uni(rng);
    ElementTransformation *trans1
            = mesh1.GetElementTransformation(elId1);
    DenseMatrix true_x(mesh1.Dimension(), ir->GetNPoints());
    Vector xk;
    for (int k=0; k < true_x.NumCols(); k++) {
        true_x.GetColumnReference(k,xk);
        trans1->Transform(ir->IntPoint(k), xk);
    }

    Array <int> elIds1(true_x.NumCols());
    Array <IntegrationPoint> ips1(true_x.NumCols());

    Array <int> elIds2(true_x.NumCols());
    Array <IntegrationPoint> ips2(true_x.NumCols());

    auto start1 = std::chrono::high_resolution_clock::now();
    mesh2.FindPoints(true_x, elIds1, ips1);
    auto end1 = std::chrono::high_resolution_clock::now();

    bool has_shared_vertices = true;
    auto start2 = std::chrono::high_resolution_clock::now();
    PointLocator point_locator(&mesh2, has_shared_vertices);
    int init_elId = 0;
    std::tie (elIds2, ips2) = point_locator(true_x,
                                            init_elId);
    auto end2 = std::chrono::high_resolution_clock::now();

    ElementTransformation *trans2 = nullptr;
    Vector x1, x2, x;
    double TOL = 1E-8;
    for (int k=0; k < true_x.NumCols(); k++)
    {
        trans2 = mesh2.GetElementTransformation(elIds1[k]);
        trans2->Transform(ips1[k], x1);

        trans2 = mesh2.GetElementTransformation(elIds2[k]);
        trans2->Transform(ips2[k], x2);

        true_x.GetColumn(k,x);

        Vector xmx1(x.Size()), xmx2(x.Size());
        subtract(x, x1, xmx1);
        subtract(x, x2, xmx2);

        ASSERT_LE(xmx1.Norml1(), TOL);
        ASSERT_LE(xmx2.Norml1(), TOL);
    }

    auto duration1 = std::chrono::duration_cast
            <std::chrono::microseconds>(end1 - start1);
    auto duration2 = std::chrono::duration_cast
            <std::chrono::microseconds>(end2 - start2);

    std::cout << "Run times: "
              << "\t slow  " << duration1.count()
              << "\t fast  " << duration2.count()
              << std::endl;
}
