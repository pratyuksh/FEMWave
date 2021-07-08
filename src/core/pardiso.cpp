#include <../include/core/pardiso.hpp>
#include <iostream>

#define PARDISO_ERROR(err, errMsg)                 \
    if (err != 0) {                                \
        std::cout << errMsg << err << std::endl;   \
        exit(1);                                   \
    }

#define PARDISO_PRINT(val, msg)                 \
    {                                           \
        std::cout << msg << val << std::endl;   \
    }

/// Constructor
PardisoSolver :: PardisoSolver (int mtype) : m_mtype(mtype) 
{
    m_verbose = 0;
    m_solver = 0;
}

/// Destructor
PardisoSolver :: ~PardisoSolver ()
{
    if (!m_finalized) { finalize(); }
}

#ifdef LIB_PARDISO

/// Initialize Pardiso solver
void PardisoSolver :: initialize(int sizeA,
                                 int *rowPtrA,
                                 int *colIdA,
                                 double *dataA)
{
    m_sizeA = sizeA;
    m_rowPtrA = rowPtrA;
    m_colIdA = colIdA;
    m_dataA = dataA;
    m_nnz = m_rowPtrA[m_sizeA];

    /// Setup Pardiso control parameters
    /// initialize internal address pointers
    char *var;
    m_error = 0;

    pardisoinit(m_pt, &m_mtype, &m_solver,
                m_iparm, m_dparm, &m_error);

    if (m_error != 0)
    {
        if (m_error == -10 )
           std::cout << "No license file found."
                     << std::endl;
        if (m_error == -11 )
           std::cout << "License is expired."
                     << std::endl;
        if (m_error == -12 )
           std::cout << "Wrong username or hostname."
                     << std::endl;
    }

    var = getenv("OMP_NUM_THREADS");
    if(var != nullptr) {
        sscanf(var, "%d", &m_nprocs);
    } else {
        std::cout << "Set environment OMP_NUM_THREADS!"
                  << std::endl;
        exit(1);
    }
    m_iparm[2] = m_nprocs;

    m_maxfct = 1;  // Max number of numerical factorizations.
    m_mnum = 1;    // Type of factorization
    m_error  = 0;

    // Shift matrix index for Fortran 1-based index
    for (int i=0; i<=m_sizeA; i++) {
        m_rowPtrA[i] += 1;
    }
    for (int i=0; i<m_nnz; i++) {
        m_colIdA[i] += 1;
    }

    // Check matrix for consistency
    pardiso_chkmatrix  (&m_mtype, &m_sizeA,
                        m_dataA, m_rowPtrA, m_colIdA,
                        &m_error);
    PARDISO_ERROR(m_error,
                  "\nERROR in consistency of matrix: ")

    m_finalized = false;
}

/// Finalize
void PardisoSolver :: finalize()
{
    // Shift matrix index to C++ 0-based index
    for (int i=0; i<m_sizeA; i++) {
        m_rowPtrA[i] -= 1;
    }
    for (int i=0; i<m_nnz; i++) {
        m_colIdA[i] -= 1;
    }

    // Release internal memory
    m_phase = -1;
    pardiso(m_pt, &m_maxfct, &m_mnum, &m_mtype, &m_phase,
            &m_sizeA, m_dataA, m_rowPtrA, m_colIdA,
            &m_idum, &m_nrhs, m_iparm, &m_verbose,
            &m_ddum, &m_ddum, &m_error, m_dparm);

    m_finalized = true;
}

/// Factorize the initialized matrix
void PardisoSolver :: factorize()
{
    m_error = 0;

    // Re-ordering and Symbolic factorization
    m_phase = 11;
    pardiso(m_pt, &m_maxfct, &m_mnum, &m_mtype, &m_phase,
            &m_sizeA, m_dataA, m_rowPtrA, m_colIdA,
            &m_idum, &m_nrhs, m_iparm, &m_verbose,
            &m_ddum, &m_ddum, &m_error, m_dparm);
    PARDISO_ERROR(m_error,
                  "\nERROR during symbolic factorization: ")

    // Numerical factorization
    m_phase = 22;
    m_iparm[32] = 1; // compute determinant
    pardiso(m_pt, &m_maxfct, &m_mnum, &m_mtype, &m_phase,
            &m_sizeA, m_dataA, m_rowPtrA, m_colIdA,
            &m_idum, &m_nrhs, m_iparm, &m_verbose,
            &m_ddum, &m_ddum, &m_error, m_dparm);
    PARDISO_ERROR(m_error,
                  "\nERROR during numerical factorization: ")
}

/// Solve Ax = b
void PardisoSolver :: solve(double *b, double *x)
{
    m_error = 0;

    pardiso_chkvec(&m_sizeA, &m_nrhs, b, &m_error);
    PARDISO_ERROR(m_error, "\nERROR in rhs vector: ")

    // Back-substitution and iterative refinement
    m_phase = 33;
    m_iparm[7] = 1;
    pardiso(m_pt, &m_maxfct, &m_mnum, &m_mtype, &m_phase,
            &m_sizeA, m_dataA, m_rowPtrA, m_colIdA,
            &m_idum, &m_nrhs, m_iparm, &m_verbose,
            b, x, &m_error, m_dparm);
    PARDISO_ERROR(m_error, "\nERROR during solve: ")
#ifdef MYVERBOSE
    double normb, normr;
    double *y = static_cast<double *>
            (calloc(static_cast<size_t>(m_sizeA),
                    sizeof(double)));
    pardiso_residual (&m_mtype, &m_sizeA,
                      m_dataA, m_rowPtrA, m_colIdA,
                      b, x, y, &normb, &normr);
    PARDISO_PRINT(normr/normb, "Relative residual norm: ")
    free(y);
#endif
}

#else

/// Initialize Pardiso solver
void PardisoSolver :: initialize(int sizeA,
                                 int *rowPtrA,
                                 int *colIdA,
                                 double *dataA)
{
    m_sizeA = sizeA;
    m_rowPtrA = rowPtrA;
    m_colIdA = colIdA;
    m_dataA = dataA;
    m_nnz = m_rowPtrA[m_sizeA];
    
    /// Setup Pardiso control parameters
    for (int i = 0; i < 64; i++) { m_iparm[i] = 0; }
    m_iparm[34] = 1;  // Zero-based indexing, C-style
    m_iparm[0] = 1;   // No solver default
    m_iparm[1] = 2;   // Fill-in reordering from METIS
    m_iparm[3] = 0;   // No iterative-direct algorithm
    m_iparm[4] = 0;   // No user fill-in reducing permutation
    m_iparm[5] = 0;   // Write solution into x
    m_iparm[6] = 0;   // Not in use
    m_iparm[7] = 2;   // Max number of itr refinement steps
    m_iparm[8] = 0;   // Not in use
    m_iparm[9] = 8;   // Perturb the pivot elements with 1E-8
    m_iparm[10] = 1;  // Use nonsymmetric permutation and scaling MPS
    m_iparm[11] = 0;  // Conjugate transposed/transpose solve
    m_iparm[12] = 1;  // Max weighted matching algo is switched-on
    m_iparm[13] = 0;  // Output: Number of perturbed pivots
    m_iparm[14] = 0;  // Not in use
    m_iparm[15] = 0;  // Not in use
    m_iparm[16] = 0;  // Not in use
    m_iparm[17] = -1; // Output: Number of nonzeros in the factor LU
    m_iparm[18] = -1; // Output: Mflops for LU factorization
    m_iparm[19] = 0;  // Output: Numbers of CG Iterations
    m_maxfct = 1;     // Maximum number of numerical factorizations.
    m_mnum = 1;       // Which factorization to use.
    m_verbose = 0;    // Message print level
    
    // Initialize internal solver pointer
    for (int i = 0; i < 64; i++) {
        m_pt[i] = 0;
    }

    m_finalized = false;
}

/// Finalize
void PardisoSolver :: finalize()
{
    // Release internal memory
    m_phase = -1;
    PARDISO(m_pt, &m_maxfct, &m_mnum, &m_mtype, &m_phase,
            &m_sizeA, m_dataA, m_rowPtrA, m_colIdA,
            &m_idum, &m_nrhs, m_iparm, &m_verbose,
            &m_ddum, &m_ddum, &m_error);
    PARDISO_ERROR(m_error,
                  "\nERROR during memory release: ")

    m_finalized = true;
}

/// Factorize the initialized matrix
void PardisoSolver :: factorize()
{
    m_error = 0;
    
    // Re-ordering and Symbolic factorization
    m_phase = 11;
    PARDISO(m_pt, &m_maxfct, &m_mnum, &m_mtype, &m_phase,
            &m_sizeA, m_dataA, m_rowPtrA, m_colIdA,
            &m_idum, &m_nrhs, m_iparm, &m_verbose,
            &m_ddum, &m_ddum, &m_error);
    //PARDISO_PRINT(m_iparm[17],
    //              "Number of nonzeros in factors: ")
    //PARDISO_PRINT(m_iparm[18],
    //              "Number of factorization MFLOPS: ")
    PARDISO_ERROR(m_error,
                  "\nERROR during symbolic factorization: ")

    // Numerical factorization
    m_phase = 22;
    PARDISO(m_pt, &m_maxfct, &m_mnum, &m_mtype, &m_phase,
            &m_sizeA, m_dataA, m_rowPtrA, m_colIdA,
            &m_idum, &m_nrhs, m_iparm, &m_verbose,
            &m_ddum, &m_ddum, &m_error);
    //PARDISO_PRINT(m_iparm[29],
    //        "Number of zero or negative pivots: ")
    PARDISO_ERROR(m_error,
                  "\nERROR during numerical factorization: ")
}

/// Solve Ax = b
void PardisoSolver :: solve(double *b, double *x)
{
    m_error = 0;

    // Back-substitution and iterative refinement
    m_phase = 33;
    PARDISO(m_pt, &m_maxfct, &m_mnum, &m_mtype, &m_phase,
            &m_sizeA, m_dataA, m_rowPtrA, m_colIdA,
            &m_idum, &m_nrhs, m_iparm, &m_verbose,
            b, x, &m_error);
    PARDISO_ERROR(m_error, "\nERROR during solve: ")
}

#endif

#undef PARDISO_ERROR
#undef PARDISO_PRINT


// End of file
