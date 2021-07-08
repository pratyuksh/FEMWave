#ifndef WAVEO1_ASSEMBLY_HPP
#define WAVEO1_ASSEMBLY_HPP

#include "mfem.hpp"
using namespace mfem;

//#include "../mymfem/mymixedbf.hpp"
#include "../mymfem/mybfinteg.hpp"


namespace mymfem {

// Gradient Integrator in time; (du/dt, v)
class TimeGradientIntegrator
        : public BilinearFormIntegrator
{
public:
    TimeGradientIntegrator() {}

    void AssembleElementMatrix
    (const FiniteElement &, ElementTransformation &,
     DenseMatrix &) override;

    void AssembleElementMatrix2
    (const FiniteElement &,
     const FiniteElement &,
     ElementTransformation &, DenseMatrix &)
    override {}

private:
#ifndef MFEM_THREAD_SAFE
    Vector shape;
    DenseMatrix dshape;
#endif
};

// Flux plus integrator in time; u(0)*v(0)
class TimeFluxPlusIntegrator
        : public BilinearFormIntegrator
{
public:
    TimeFluxPlusIntegrator() {}

    void AssembleElementMatrix
    (const FiniteElement &, ElementTransformation &,
     DenseMatrix &) override;

    void AssembleElementMatrix2
    (const FiniteElement &,
     const FiniteElement &,
     ElementTransformation &, DenseMatrix &)
    override {}

private:
#ifndef MFEM_THREAD_SAFE
    Vector shape;
#endif
};

// Flux minus integrator in time; u(1)*v(0)
class TimeFluxMinusIntegrator
        : public BilinearFormIntegrator
{
public:
    TimeFluxMinusIntegrator() {}

    void AssembleElementMatrix
    (const FiniteElement &, ElementTransformation &,
     DenseMatrix &) override;

    void AssembleElementMatrix2
    (const FiniteElement &,
     const FiniteElement &,
     ElementTransformation &, DenseMatrix &)
    override {}

private:
#ifndef MFEM_THREAD_SAFE
    Vector shape1, shape2;
#endif
};

// Vector Gradient Integrator; (grad(u), v)
class VectorGradientIntegrator
        : public BilinearFormIntegrator
{
public:
    VectorGradientIntegrator() {}

    VectorGradientIntegrator(MatrixCoefficient *M_)
        : M(M_) {}

    VectorGradientIntegrator(MatrixCoefficient& M_)
        : M(&M_) {}

    VectorGradientIntegrator(Coefficient *q_)
        : q(q_) {}

    VectorGradientIntegrator(Coefficient& q_)
        : q(&q_) {}

    void AssembleElementMatrix
    (const FiniteElement &, ElementTransformation &,
     DenseMatrix &) override {};

    void AssembleElementMatrix2
    (const FiniteElement &,
     const FiniteElement &,
     ElementTransformation &, DenseMatrix &)
    override;

private:
    void AssembleBlock(const int dim,
                       const int test_ndofs,
                       const int trial_ndofs,
                       const Vector& shape,
                       const DenseMatrix& dshape,
                       DenseMatrix& outMat) const;

    MatrixCoefficient *M = nullptr;
    Coefficient *q = nullptr;
#ifndef MFEM_THREAD_SAFE
    Vector shape;
    DenseMatrix dshape;
#endif
};

// Flux Block at position 0,0 Integrator
class FluxBlock00Integrator
        : public BilinearFormIntegrator
{
public:
    explicit FluxBlock00Integrator(double hx, int stabParamsType=1)
        : m_hx(hx), m_stabParamsType(stabParamsType) {}

    void AssembleFaceMatrix (const FiniteElement &fe1,
                             const FiniteElement &fe2,
                             FaceElementTransformations &ftr,
                             DenseMatrix &elmat) override;

private:
    inline void AssembleBlock (const int row_ndofs,
                               const int col_ndofs,
                               const int row_offset,
                               const int col_offset,
                               const Vector &row_shape,
                               const Vector &col_shape,
                               const double sCoeff,
                               DenseMatrix &elmat) const;

#ifndef MFEM_THREAD_SAFE
    mutable Vector shape1, shape2;
#endif
    double m_hx;
    int m_stabParamsType;
};

// Flux Block at position 1,1 Integrator
class FluxBlock11Integrator
        : public BilinearFormIntegrator
{
public:
    explicit FluxBlock11Integrator(double hx, int stabParamsType=1)
        : m_hx(hx), m_stabParamsType(stabParamsType) {}

    void AssembleFaceMatrix (const FiniteElement &fe1,
                             const FiniteElement &fe2,
                             FaceElementTransformations &ftr,
                             DenseMatrix &elmat) override;

private:
    inline void AssembleBlock (const int dim,
                               const int row_ndofs,
                               const int col_ndofs,
                               const int row_offset,
                               const int col_offset,
                               const Vector &row_shape,
                               const Vector &col_shape,
                               const Vector& nor,
                               const double sCoeff,
                               DenseMatrix &elmat) const;

#ifndef MFEM_THREAD_SAFE
    mutable Vector shape1, shape2, nor;
#endif
    double m_hx;
    int m_stabParamsType;
};

// Flux Block at position 1,0 Integrator
class FluxBlock10Integrator
        : public mymfem::MyBilinearFormIntegrator
{
public:
    void MyAssembleFaceMatrix (const FiniteElement &trial_fe1,
                               const FiniteElement &trial_fe2,
                               const FiniteElement &test_fe1,
                               const FiniteElement &test_fe2,
                               FaceElementTransformations &ftr,
                               DenseMatrix &elmat) const override;

private:
    inline void AssembleBlock (const int dim,
                               const int row_ndofs,
                               const int col_ndofs,
                               const int row_offset,
                               const int col_offset,
                               const Vector &row_shape,
                               const Vector &col_shape,
                               const Vector& nor,
                               const double sCoeff,
                               DenseMatrix &elmat) const;

#ifndef MFEM_THREAD_SAFE
    mutable Vector nor;
    mutable Vector tr_shape1, tr_shape2;
    mutable Vector te_shape1, te_shape2;
#endif
};

// Flux Block at position 0,1 Integrator
class FluxBlock01Integrator
        : public mymfem::MyBilinearFormIntegrator
{
public:
    void MyAssembleFaceMatrix (const FiniteElement &trial_fe1,
                               const FiniteElement &trial_fe2,
                               const FiniteElement &test_fe1,
                               const FiniteElement &test_fe2,
                               FaceElementTransformations &ftr,
                               DenseMatrix &elmat) const override;

private:
    inline void AssembleBlock (const int dim,
                               const int row_ndofs,
                               const int col_ndofs,
                               const int row_offset,
                               const int col_offset,
                               const Vector &row_shape,
                               const Vector &col_shape,
                               const Vector& nor,
                               const double sCoeff,
                               DenseMatrix &elmat) const;

#ifndef MFEM_THREAD_SAFE
    mutable Vector nor;
    mutable Vector tr_shape1, tr_shape2;
    mutable Vector te_shape1, te_shape2;
#endif
};

// Dirichlet Block at position 1,0 Integrator
class DirichletBlock10Integrator
        : public mymfem::MyBilinearFormIntegrator
{
public:
    void MyAssembleFaceMatrix (const FiniteElement &trial_fe1,
                               const FiniteElement &trial_fe2,
                               const FiniteElement &test_fe1,
                               const FiniteElement &test_fe2,
                               FaceElementTransformations &ftr,
                               DenseMatrix &elmat) const override;

private:
    inline void AssembleBlock (const int dim,
                               const int row_ndofs,
                               const int col_ndofs,
                               const int row_offset,
                               const int col_offset,
                               const Vector &row_shape,
                               const Vector &col_shape,
                               const Vector& nor,
                               const double sCoeff,
                               DenseMatrix &elmat) const;

#ifndef MFEM_THREAD_SAFE
    mutable Vector nor;
    mutable Vector tr_shape, te_shape;
#endif
};

// Neumann Block at position 0,1 Integrator
class NeumannBlock01Integrator
        : public mymfem::MyBilinearFormIntegrator
{
public:
    void MyAssembleFaceMatrix (const FiniteElement &trial_fe1,
                               const FiniteElement &trial_fe2,
                               const FiniteElement &test_fe1,
                               const FiniteElement &test_fe2,
                               FaceElementTransformations &ftr,
                               DenseMatrix &elmat) const override;

private:
    inline void AssembleBlock (const int dim,
                               const int row_ndofs,
                               const int col_ndofs,
                               const int row_offset,
                               const int col_offset,
                               const Vector &row_shape,
                               const Vector &col_shape,
                               const Vector& nor,
                               const double sCoeff,
                               DenseMatrix &elmat) const;

#ifndef MFEM_THREAD_SAFE
    mutable Vector nor;
    mutable Vector tr_shape, te_shape;
#endif
};


// Dirichlet linear form block at position 0 Integrator
class DirichletLFBlock0Integrator
        : public LinearFormIntegrator
{
public:
    DirichletLFBlock0Integrator(Coefficient &Qf,
                                int stabParamsType = 1)
        : Q(Qf), m_stabParamsType (stabParamsType) {}

    void AssembleRHSElementVect  (const FiniteElement &el,
                                  FaceElementTransformations &ftr,
                                  Vector &elvect) override;

    void AssembleRHSElementVect  (const FiniteElement &el,
                                  ElementTransformation &tr,
                                  Vector &elvect) override {}

private:
    Coefficient &Q;
    Vector shape;
    int m_stabParamsType;
};

// Dirichlet linear form block at position 1 Integrator
class DirichletLFBlock1Integrator
        : public LinearFormIntegrator
{
public:
    DirichletLFBlock1Integrator(Coefficient &Qf, double s=1)
        : Q(Qf), factor(s) {}

    void AssembleRHSElementVect  (const FiniteElement &el,
                                  FaceElementTransformations &ftr,
                                  Vector &elvect) override;

    void AssembleRHSElementVect  (const FiniteElement &el,
                                  ElementTransformation &tr,
                                  Vector &elvect) override {}

private:
    Coefficient &Q;
    Vector nor, shape;
    double factor;
};


// Neumann linear form block at position 0 Integrator
class NeumannLFBlock0Integrator
        : public LinearFormIntegrator
{
public:
    NeumannLFBlock0Integrator(VectorCoefficient &Qf, double s=1)
        : Q(Qf), factor(s) {}

    void AssembleRHSElementVect  (const FiniteElement &el,
                                  FaceElementTransformations &ftr,
                                  Vector &elvect) override;

    void AssembleRHSElementVect  (const FiniteElement &el,
                                  ElementTransformation &tr,
                                  Vector &elvect) override {}

private:
    VectorCoefficient &Q;
    Vector nor, shape;
    double factor;
};

// Neumann linear form block at position 1 Integrator
class NeumannLFBlock1Integrator
        : public LinearFormIntegrator
{
public:
    NeumannLFBlock1Integrator(VectorCoefficient &Qf,
                              int stabParamsType = 1)
        : Q(Qf), m_stabParamsType (stabParamsType) {}

    void AssembleRHSElementVect  (const FiniteElement &el,
                                  FaceElementTransformations &ftr,
                                  Vector &elvect) override;

    void AssembleRHSElementVect  (const FiniteElement &el,
                                  ElementTransformation &tr,
                                  Vector &elvect) override {}

private:
    VectorCoefficient &Q;
    Vector nor, shape;
    int m_stabParamsType;
};

}

#endif /// WAVEO1_ASSEMBLY_HPP
