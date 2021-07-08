#ifndef PARDISO_HPP
#define PARDISO_HPP

#include "config.hpp"

// PARDISO prototype.
#ifdef LIB_PARDISO
extern "C" void pardisoinit (void *, int    *, int *,
                             int  *, double *, int *);
extern "C" void pardiso     (void *, int    *, int *, int    *, int    *,
                             int  *, double *, int *, int    *, int    *,
                             int  *, int    *, int *, double *, double *,
                             int  *, double *);
extern "C" void pardiso_chkmatrix  (int *, int *, double *,
                                    int *, int *, int    *);
extern "C" void pardiso_chkvec     (int *, int *, double *, int *);
extern "C" void pardiso_residual (int *, int * ,
                                  double *, int *, int *,
                                  double *, double *, double *,
                                  double *, double *);
#else // MKL Pardiso
extern "C" void PARDISO     (void *, int    *, int *, int    *, int    *,
                             int  *, double *, int *, int    *, int    *,
                             int  *, int    *, int *, double *, double *,
                             int  *);
#endif


class PardisoSolver
{
public:
    PardisoSolver (int);

    PardisoSolver (const nlohmann::json& config)
    {
        m_verbose = config["pardiso_verbose"];

        m_mtype = 11; // real unsymemtric matrix
        m_solver = 0; // sparse direct solver
    }

    ~PardisoSolver ();

    void initialize (int, int *, int *, double *);
    void finalize();

    void factorize ();
    void solve (double *, double *);

private:
    int m_nnodes, m_nprocs;

    int m_mtype;
    int m_sizeA, m_nnz;
    int *m_rowPtrA, *m_colIdA; // CSR storage format
    double *m_dataA;

    int m_nrhs = 1;
    int m_solver;
    int m_maxfct, m_mnum, m_phase, m_error, m_verbose;

    int m_idum;
    double m_ddum;

    void *m_pt[64];
    int m_iparm[64];
    double m_dparm[64];

    bool m_finalized = true;
};

#endif /// PARDISO_HPP
