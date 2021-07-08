#ifndef WAVEO1_DISCRETISATION_HPP
#define WAVEO1_DISCRETISATION_HPP

#include "mfem.hpp"
using namespace mfem;

#include <iostream>

#include "../core/config.hpp"
#include "test_cases.hpp"
#include "coefficients.hpp"


class WaveO1XtDG
{
public:
    WaveO1XtDG (const nlohmann::json&,
              std::shared_ptr<WaveO1TestCases>&);

    ~WaveO1XtDG ();
    
    void set(std::shared_ptr<Mesh>&,
             std::shared_ptr<Mesh>&);

    void assemble_system();
    void build_system_op();
    void build_system_matrix();

    void assemble_rhs_initial(BlockVector *) const;
    void assemble_rhs(const int, const BlockVector *, BlockVector *) const;
    void assemble_ics(BlockVector *) const;
    void assemble_source(const int, BlockVector *) const;
    void assemble_bcs(const int, BlockVector *) const;

    void set_vector(const int,
                    const Vector&,
                    Vector&) const;

    // for L2-projection
    void assemble_projector();
    void build_projector_matrix();
    void assemble_projection_rhs(int, BlockVector *) const;

private:
    void assemble_ics_xProjection(Vector&,
                                  Vector&) const;

    void assemble_source_xProjection(Vector&,
                                     Vector&,
                                     double) const;

    void assemble_bcs_xProjection(Vector&,
                                  Vector&,
                                  double t) const;

    void add_vector(const Array<int> &,
                    const Vector&,
                    Vector&) const;

    void assemble_xProjection_rhs(Vector&, Vector&,
                                  double) const;

public:
    inline std::shared_ptr<WaveO1TestCases> get_test_case() const {
        return m_testCase;
    }
    
    inline Mesh* get_tMesh() const {
        return m_tFespace->GetMesh();
    }
    
    inline Mesh* get_xMesh() const {
        return m_xFespaces[0]->GetMesh();
    }

    inline FiniteElementCollection* get_tFec() const {
        return m_tFec;
    }

    inline Array<FiniteElementCollection*> get_xFecs() const {
        return m_xFecs;
    }
    
    inline FiniteElementSpace* get_tFespace() const {
        return m_tFespace;
    }
    
    inline Array<FiniteElementSpace*> get_xFespaces() const {
        return m_xFespaces;
    }

    inline BlockOperator* get_wave_op() const {
        return m_waveOp;
    }

    inline Array<int> get_block_offsets() const {
        return m_block_offsets;
    }

    inline Array<int> get_slab_block_offsets() const {
        return m_slab_block_offsets;
    }

    inline Array<int> get_xEss_bdr_marker() const {
        return m_xEss_bdr_marker;
    }

    inline Array<int> get_xNat_bdr_marker() const {
        return m_xNat_bdr_marker;
    }

    inline SparseMatrix* get_wave_mat() const {
        return m_waveMat;
    }

    inline SparseMatrix* get_wave_projection_mat() const {
        return m_waveProjMat;
    }

    inline void print() const {
        std::cout << "Degree of pressure in x: "
                  << m_xDeg1 << std::endl;
        std::cout << "Degree of velocity in x: "
                  << m_xDeg2 << std::endl;
        std::cout << "Degree of pressure and velocity in t: "
                  << m_tDeg << std::endl;
    }

private:
    const nlohmann::json& m_config;
    
    int m_stabParamsType;

    int m_tDeg;
    double m_endTime;

    int m_xndim, m_xDeg1, m_xDeg2;

    std::shared_ptr<WaveO1TestCases> m_testCase;
    std::shared_ptr<WaveO1InvSqMediumCoeff> m_invSqMed;
    
    std::unique_ptr<Mesh> m_tMesh0;

    FiniteElementCollection* m_tFec = nullptr;
    FiniteElementSpace* m_tFespace = nullptr;
    FiniteElementSpace* m_tFespace0 = nullptr;
    
    Array<FiniteElementCollection*> m_xFecs;
    Array<FiniteElementSpace*> m_xFespaces;
    Array<int> m_block_offsets;
    Array<int> m_slab_block_offsets;
    
    SparseMatrix *m_slab_tMass = nullptr;
    SparseMatrix *m_slab_tGrad = nullptr;
    SparseMatrix *m_slab_tFluxp = nullptr;
    SparseMatrix *m_slab_tFluxm = nullptr;
    
    SparseMatrix *m_xMass0 = nullptr;
    SparseMatrix *m_xMass1 = nullptr;
    SparseMatrix *m_xGrad = nullptr;
    SparseMatrix *m_xDiv = nullptr;

    SparseMatrix *m_xFlux00 = nullptr;
    SparseMatrix *m_xFlux01 = nullptr;
    SparseMatrix *m_xFlux10 = nullptr;
    SparseMatrix *m_xFlux11 = nullptr;
    SparseMatrix *m_xDirichlet00 = nullptr;
    SparseMatrix *m_xDirichlet10 = nullptr;
    SparseMatrix *m_xNeumann01 = nullptr;
    SparseMatrix *m_xNeumann11 = nullptr;
    
    SparseMatrix *m_slab_rhsMat0 = nullptr;
    SparseMatrix *m_slab_rhsMat1 = nullptr;

    SparseMatrix *m_block00 = nullptr;
    SparseMatrix *m_block01 = nullptr;
    SparseMatrix *m_block10 = nullptr;
    SparseMatrix *m_block11 = nullptr;
    
    BlockOperator *m_waveOp = nullptr;
    SparseMatrix *m_waveMat = nullptr;
    SparseMatrix *m_waveProjMat = nullptr;
    
    mutable Array<int> m_xEss_bdr_marker;
    mutable Array<int> m_xNat_bdr_marker;
};


#endif /// WaveO1_DISCRETISATION_HPP
