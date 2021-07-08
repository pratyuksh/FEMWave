#ifndef WAVEO2_DISCRETISATION_HPP
#define WAVEO2_DISCRETISATION_HPP

#include "mfem.hpp"
using namespace mfem;

#include <iostream>

#include "../core/config.hpp"
#include "test_cases.hpp"
#include "coefficients.hpp"


class WaveO2XtFEM
{
public:
    //! Constructors
    WaveO2XtFEM (const nlohmann::json&,
                 std::shared_ptr<WaveO2TestCases>&);

    ~WaveO2XtFEM ();

    //! Sets-up the FE spaces and boundary conditions
    void set(std::shared_ptr<Mesh>&, std::shared_ptr<Mesh>&);

    //! Assembles the matrices
    void assemble_system();

    //! Builds the system matrix
    //void build_system_matrix();

    //void assemble_rhs(Vector *) const;
    //void assemble_source(Vector *) const;

private:
    //void apply_BCs(SparseMatrix&) const;
    //void apply_BCs(Vector&) const;

public:
    //! Returns test case
    inline std::shared_ptr<WaveO2TestCases>
    get_test_case() const {
        return m_testCase;
    }

    //! Returns temporal mesh
    inline Mesh* get_tMesh() const {
        return m_tFespace->GetMesh();
    }

    //! Returns spatial mesh
    inline Mesh* get_xMesh() const {
        return m_xFespace->GetMesh();
    }

    //! Returns FE collection in time
    inline FiniteElementCollection* get_tFec() const {
        return m_tFec;
    }

    //! Returns FE collection in space
    inline FiniteElementCollection* get_xFecs() const {
        return m_xFec;
    }

    //! Returns FEspace in time
    inline FiniteElementSpace* get_tFespace() const {
        return m_tFespace;
    }

    //! Returns FEspace in x
    inline FiniteElementSpace* get_xFespaces() const {
        return m_xFespace;
    }

    //! Returns essential boundary marker
    inline Array<int> get_xEss_bdr_marker() const {
        return m_xEss_bdr_marker;
    }

    //! Returns natural boundary marker
    inline Array<int> get_xNat_bdr_marker() const {
        return m_xNat_bdr_marker;
    }

    //! Returns system matrix
    inline SparseMatrix* get_wave_mat() const {
        return m_waveMat;
    }

    //! Prints info
    inline void print() const {
        std::cout << "Degree in x: "
                  << m_xDeg << std::endl;
        std::cout << "Degree in t: "
                  << m_tDeg << std::endl;
    }

private:
    const nlohmann::json& m_config;

    int m_tDeg;
    double m_endTime;

    int m_xndim, m_xDeg;

    std::shared_ptr<WaveO2TestCases> m_testCase;
    std::shared_ptr<WaveO2SqMediumCoeff> m_sqMed;

    std::unique_ptr<Mesh> m_tMesh0;

    FiniteElementCollection* m_tFec = nullptr;
    FiniteElementSpace* m_tFespace = nullptr;

    FiniteElementCollection* m_xFec = nullptr;
    FiniteElementSpace* m_xFespace = nullptr;

    SparseMatrix *m_tMass = nullptr;
    SparseMatrix *m_tStiff = nullptr;

    SparseMatrix *m_xMass = nullptr;
    SparseMatrix *m_xStiff = nullptr;

    SparseMatrix *m_waveMat = nullptr;

    mutable Array<int> m_xEss_bdr_marker;
    mutable Array<int> m_xNat_bdr_marker;
    Array<int> m_xEss_tdof_list;
    Array<int> m_ess_tdof_list;
};


#endif /// WAVEO2_DISCRETISATION_HPP
