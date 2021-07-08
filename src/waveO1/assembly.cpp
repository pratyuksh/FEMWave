#include "../../include/waveO1/assembly.hpp"
#include <assert.h>


// Gradient Integrator in time, dim = 1
void mymfem::TimeGradientIntegrator
:: AssembleElementMatrix (const FiniteElement &fe,
                          ElementTransformation &Trans,
                          DenseMatrix &elmat)
{
    int dim = fe.GetDim();
    assert(dim == 1);

    int ndofs = fe.GetDof();
#ifdef MFEM_THREAD_SAFE
    Vector shape(ndofs);
    DenseMatrix dshape(ndofs, dim);
#else
    shape.SetSize(ndofs);
    dshape.SetSize(ndofs, dim);
#endif
    elmat.SetSize(ndofs, ndofs);

    const IntegrationRule *ir = IntRule;
    if (ir == nullptr)
    {
        int order = 2*fe.GetOrder()-1;
        ir = &IntRules.Get(fe.GetGeomType(), order);
    }

    elmat = 0.0;
    for (int i = 0; i < ir->GetNPoints(); i++)
    {
        const IntegrationPoint &ip = ir->IntPoint(i);
        Trans.SetIntPoint(&ip);

        fe.CalcShape(ip, shape);
        fe.CalcPhysDShape(Trans, dshape);
        double w = ip.weight*Trans.Weight();

        Vector t;
        Trans.Transform(ip,t);

        for (int l=0; l<ndofs; l++) {
            for (int k=0; k<ndofs; k++) {
                elmat(k,l) += w*shape(k)*dshape(l,0);
            }
        }
    }
}

// Flux plus integrator in time; u(0)*v(0)
void mymfem::TimeFluxPlusIntegrator
:: AssembleElementMatrix (const FiniteElement &fe,
                          ElementTransformation &,
                          DenseMatrix &elmat)
{
    int dim = fe.GetDim();
    assert(dim == 1);

    int ndofs = fe.GetDof();
#ifdef MFEM_THREAD_SAFE
    Vector shape(ndofs);
#else
    shape.SetSize(ndofs);
#endif
    elmat.SetSize(ndofs, ndofs);

    IntegrationPoint ip;
    ip.Set1w(0.0, 1.0);
    fe.CalcShape(ip, shape);

    for (int l=0; l<ndofs; l++) {
        for (int k=0; k<ndofs; k++) {
            elmat(k,l) = shape(k)*shape(l);
        }
    }
}

// Flux minus integrator in time; u(1)*v(0)
void mymfem::TimeFluxMinusIntegrator
:: AssembleElementMatrix (const FiniteElement &fe,
                          ElementTransformation &,
                          DenseMatrix &elmat)
{
    int dim = fe.GetDim();
    assert(dim == 1);

    int ndofs = fe.GetDof();
#ifdef MFEM_THREAD_SAFE
    Vector shape1(ndofs);
    Vector shape2(ndofs);
#else
    shape1.SetSize(ndofs);
    shape2.SetSize(ndofs);
#endif
    elmat.SetSize(ndofs, ndofs);
    IntegrationPoint ip1, ip2;
    ip1.Set1w(0.0, 1.0);
    ip2.Set1w(1.0, 1.0);
    fe.CalcShape(ip1, shape1);
    fe.CalcShape(ip2, shape2);

    for (int l=0; l<ndofs; l++) {
        for (int k=0; k<ndofs; k++) {
            elmat(k,l) = shape1(k)*shape2(l);
        }
    }
}

// Vector Gradient Integrator
void mymfem::VectorGradientIntegrator
:: AssembleElementMatrix2 (const FiniteElement &trial_fe,
                           const FiniteElement &test_fe,
                           ElementTransformation &Trans,
                           DenseMatrix &elmat)
{
    int dim = trial_fe.GetDim();
    int trial_ndofs = trial_fe.GetDof();
    int test_ndofs = test_fe.GetDof();
#ifdef MFEM_THREAD_SAFE
    DenseMatrix dshape(trial_ndofs, dim);
    Vector shape(test_ndofs);
#else
    dshape.SetSize(trial_ndofs, dim);
    shape.SetSize(test_ndofs);
#endif
    elmat.SetSize(test_ndofs*dim, trial_ndofs);

    const IntegrationRule *ir = IntRule;
    if (ir == nullptr)
    {
        int order
                = trial_fe.GetOrder()+test_fe.GetOrder()-1;
        ir = &IntRules.Get(trial_fe.GetGeomType(), order);
    }

    elmat = 0.0;
    for (int i = 0; i < ir->GetNPoints(); i++)
    {
        const IntegrationPoint &ip = ir->IntPoint(i);
        Trans.SetIntPoint(&ip);

        trial_fe.CalcPhysDShape(Trans, dshape);
        test_fe.CalcShape(ip, shape);
        double w = ip.weight*Trans.Weight();

        DenseMatrix tmpMat(dim*test_ndofs, trial_ndofs);
        if (M) { // matrix coefficient
            DenseMatrix matM(dim);
            DenseMatrix buf(trial_ndofs, dim);
            M->Eval(matM, Trans, ip);
            MultABt(dshape, matM, buf);
            AssembleBlock(dim, test_ndofs, trial_ndofs,
                          shape, buf, tmpMat);
        }
        else {
            AssembleBlock(dim, test_ndofs, trial_ndofs,
                          shape, dshape, tmpMat);
        }

        if (q) { // scalar coefficient
            double val = q->Eval(Trans, ip);
            w *= val;
        }
        elmat.Add(w,tmpMat);
    }
}

void mymfem::VectorGradientIntegrator
:: AssembleBlock(const int dim,
                 const int test_ndofs,
                 const int trial_ndofs,
                 const Vector& shape,
                 const DenseMatrix& dshape,
                 DenseMatrix& outMat) const
{
    for (int k=0; k<dim; k++)
        for (int i=0; i<test_ndofs; i++)
            for (int j=0; j<trial_ndofs; j++)
            {
                outMat(i + k*test_ndofs, j)
                        = shape(i)*dshape(j,k);
            }
}

// Flux Block at position 0,0 Integrator
// For stabParamsType 1 or 3, alpha = 1
// For stabParamsType 2 or 4, alpha = 1/hF
void mymfem::FluxBlock00Integrator
:: AssembleFaceMatrix (const FiniteElement &fe1,
                       const FiniteElement &fe2,
                       FaceElementTransformations &ftr,
                       DenseMatrix &elmat)
{
#ifdef MFEM_THREAD_SAFE
    Vector shape1, shape2;
#endif
    int ndofs1 = fe1.GetDof();
    int ndofs2 = (ftr.Elem2No >= 0) ? fe2.GetDof() : 0;
    int ndofs = ndofs1 + ndofs2;

    shape1.SetSize(ndofs1);
    if (ndofs2) { shape2.SetSize(ndofs2); }
    elmat.SetSize(ndofs, ndofs);

    const IntegrationRule *ir = IntRule;
    if (ir == nullptr)
    {
        int order = 2*fe1.GetOrder();
        ir = &IntRules.Get(ftr.FaceGeom, order);
    }

    elmat = 0.0;
    double hF = 0;
    for (int i = 0; i < ir->GetNPoints(); i++)
    {
        const IntegrationPoint &ip = ir->IntPoint(i);
        ftr.Face->SetIntPoint(&ip);
        if (i == 0) { hF = ftr.Face->Weight(); }

        IntegrationPoint eip1, eip2;
        ftr.Loc1.Transform(ip, eip1);
        fe1.CalcShape(eip1, shape1);
        if (ndofs2) {
            ftr.Loc2.Transform(ip, eip2);
            fe2.CalcShape(eip2, shape2);
        }

        double coeff = ip.weight;
        //if (m_stabParamsType == 1 || m_stabParamsType == 3) {
        //    coeff *= ftr.Face->Weight();
        //}

        // (1,1) block
        AssembleBlock(ndofs1, ndofs1,
                      0, 0,
                      shape1, shape1,
                      +coeff, elmat);
        if (ndofs2)
        {
            // (1,2) block
            AssembleBlock(ndofs1, ndofs2,
                          0, ndofs1,
                          shape1, shape2,
                          -coeff, elmat);

            // (2,1) block
            AssembleBlock(ndofs2, ndofs1,
                          ndofs1, 0,
                          shape2, shape1,
                          -coeff, elmat);

            // (2,2) block
            AssembleBlock(ndofs2, ndofs2,
                          ndofs1, ndofs1,
                          shape2, shape2,
                          +coeff, elmat);
        }
    }
    if (m_stabParamsType == 1 || m_stabParamsType == 3) {
        elmat *= hF;
    }
    if (m_stabParamsType == 5) {
        elmat *= m_hx;
    }
}

void mymfem::FluxBlock00Integrator
:: AssembleBlock (const int row_ndofs,
                  const int col_ndofs,
                  const int row_offset,
                  const int col_offset,
                  const Vector &row_shape,
                  const Vector &col_shape,
                  const double sCoeff,
                  DenseMatrix &elmat) const
{
    for (int j=0; j < col_ndofs; j++)
        for (int i=0; i < row_ndofs; i++)
        {
            elmat(i+row_offset, j+col_offset)
                    += sCoeff*row_shape(i)*col_shape(j);
        }
}

// Flux Block at position 1,1 Integrator
// For stabParamsType 1 or 2, beta = 1
// For stabParamsType 3 or 4, beta = hF
void mymfem::FluxBlock11Integrator
:: AssembleFaceMatrix (const FiniteElement &fe1,
                       const FiniteElement &fe2,
                       FaceElementTransformations &ftr,
                       DenseMatrix &elmat)
{
    int dim = fe1.GetDim();
#ifdef MFEM_THREAD_SAFE
    Vector shape1, shape2, nor;
#endif
    int ndofs1 = fe1.GetDof();
    int ndofs2 = (ftr.Elem2No >= 0) ? fe2.GetDof() : 0;
    int ndofs = ndofs1 + ndofs2;

    nor.SetSize(dim);
    shape1.SetSize(ndofs1);
    if (ndofs2) { shape2.SetSize(ndofs2); }
    elmat.SetSize(dim*ndofs, dim*ndofs);

    const IntegrationRule *ir = IntRule;
    if (ir == nullptr)
    {
        int order = 2*fe1.GetOrder();
        ir = &IntRules.Get(ftr.FaceGeom, order);
    }

    elmat = 0.0;
    //double hF = 0;
    for (int i = 0; i < ir->GetNPoints(); i++)
    {
        const IntegrationPoint &ip = ir->IntPoint(i);
        ftr.Face->SetIntPoint(&ip);
        CalcOrtho(ftr.Face->Jacobian(), nor);
        nor /= ftr.Face->Weight(); // unit normal
        //if (i == 0) { hF = ftr.Face->Weight(); }

        IntegrationPoint eip1, eip2;
        ftr.Loc1.Transform(ip, eip1);
        fe1.CalcShape(eip1, shape1);
        if (ndofs2) {
            ftr.Loc2.Transform(ip, eip2);
            fe2.CalcShape(eip2, shape2);
        }

        double coeff = ip.weight;//*ftr.Face->Weight();
        if (m_stabParamsType == 1 || m_stabParamsType == 2) {
            coeff *= ftr.Face->Weight();
        }
        else if (m_stabParamsType == 3 || m_stabParamsType == 4) {
            coeff *= (ftr.Face->Weight()*ftr.Face->Weight());
        }
        else if (m_stabParamsType == 5) {
            coeff *= (ftr.Face->Weight()*ftr.Face->Weight());
            coeff /= m_hx;
        }

        // (1,1) block
        AssembleBlock(dim, ndofs1, ndofs1,
                      0, 0,
                      shape1, shape1,
                      nor, +coeff, elmat);
        if (ndofs2)
        {
            // (1,2) block
            AssembleBlock(dim, ndofs1, ndofs2,
                          0, dim*ndofs1,
                          shape1, shape2,
                          nor, -coeff, elmat);

            // (2,1) block
            AssembleBlock(dim, ndofs2, ndofs1,
                          dim*ndofs1, 0,
                          shape2, shape1,
                          nor, -coeff, elmat);

            // (2,2) block
            AssembleBlock(dim, ndofs2, ndofs2,
                          dim*ndofs1, dim*ndofs1,
                          shape2, shape2,
                          nor, +coeff, elmat);
        }
    }
    //if (m_stabParamsType == 3 || m_stabParamsType == 4) {
    //    elmat *= hF;
    //}
}

void mymfem::FluxBlock11Integrator
:: AssembleBlock (const int dim,
                  const int row_ndofs,
                  const int col_ndofs,
                  const int row_offset,
                  const int col_offset,
                  const Vector &row_shape,
                  const Vector &col_shape,
                  const Vector &nor,
                  const double sCoeff,
                  DenseMatrix &elmat) const
{
    for (int l=0; l<dim; l++)
        for (int k=0; k<dim; k++)
        {
            double tmp = sCoeff*nor(k)*nor(l);
            for (int j=0; j < col_ndofs; j++)
                for (int i=0; i < row_ndofs; i++)
                {
                    elmat(i+k*row_ndofs+row_offset,
                          j+l*col_ndofs+col_offset)
                            += tmp*row_shape(i)*col_shape(j);
                }
        }
}

// Flux Block at position 1,0 Integrator
void mymfem::FluxBlock10Integrator
:: MyAssembleFaceMatrix (const FiniteElement &tr_fe1,
                         const FiniteElement &tr_fe2,
                         const FiniteElement &te_fe1,
                         const FiniteElement &te_fe2,
                         FaceElementTransformations &ftr,
                         DenseMatrix &elmat) const
{
    int dim = te_fe1.GetDim();
    int tr_ndofs1 = tr_fe1.GetDof();
    int tr_ndofs2 = tr_fe2.GetDof();
    int tr_ndofs = tr_ndofs1 + tr_ndofs2;

    int te_ndofs1 = te_fe1.GetDof();
    int te_ndofs2 = te_fe2.GetDof();
    int te_ndofs = te_ndofs1 + te_ndofs2;
#ifdef MFEM_THREAD_SAFE
    Vector nor(dim);
    Vector tr_shape1(tr_ndofs1), tr_shape2(tr_ndofs2);
    Vector te_shape1(te_ndofs1), te_shape2(te_ndofs2);
#else
    nor.SetSize(dim);
    tr_shape1.SetSize(tr_ndofs1);
    tr_shape2.SetSize(tr_ndofs2);
    te_shape1.SetSize(te_ndofs1);
    te_shape2.SetSize(te_ndofs2);
#endif

    elmat.SetSize(dim*te_ndofs, tr_ndofs);

    const IntegrationRule *ir = IntRule;
    if (ir == nullptr)
    {
        int order = tr_fe1.GetOrder()
                +te_fe1.GetOrder();
        ir = &IntRules.Get(ftr.FaceGeom, order);
    }

    elmat = 0.0;
    for (int i = 0; i < ir->GetNPoints(); i++)
    {
        const IntegrationPoint &ip = ir->IntPoint(i);
        ftr.Face->SetIntPoint(&ip);
        CalcOrtho(ftr.Face->Jacobian(), nor);

        IntegrationPoint eip1, eip2;
        ftr.Loc1.Transform(ip, eip1);
        ftr.Loc2.Transform(ip, eip2);
        tr_fe1.CalcShape(eip1, tr_shape1);
        tr_fe2.CalcShape(eip2, tr_shape2);
        te_fe1.CalcShape(eip1, te_shape1);
        te_fe2.CalcShape(eip2, te_shape2);

        double coeff = 0.5*ip.weight;
        // (1,1) block
        AssembleBlock(dim, te_ndofs1, tr_ndofs1,
                      0, 0,
                      te_shape1, tr_shape1,
                      nor, -coeff, elmat);

        // (1,2) block
        AssembleBlock(dim, te_ndofs1, tr_ndofs2,
                      0, tr_ndofs1,
                      te_shape1, tr_shape2,
                      nor, +coeff, elmat);

        // (2,1) block
        AssembleBlock(dim, te_ndofs2, tr_ndofs1,
                      dim*te_ndofs1, 0,
                      te_shape2, tr_shape1,
                      nor, -coeff, elmat);

        // (2,2) block
        AssembleBlock(dim, te_ndofs2, tr_ndofs2,
                      dim*te_ndofs1, tr_ndofs1,
                      te_shape2, tr_shape2,
                      nor, +coeff, elmat);
    }
}

void mymfem::FluxBlock10Integrator
:: AssembleBlock (const int dim,
                  const int row_ndofs,
                  const int col_ndofs,
                  const int row_offset,
                  const int col_offset,
                  const Vector &row_shape,
                  const Vector &col_shape,
                  const Vector &nor,
                  const double sCoeff,
                  DenseMatrix &elmat) const
{
    for (int k=0; k<dim; k++)
    {
        double tmp = sCoeff*nor(k);
        for (int j=0; j < col_ndofs; j++)
            for (int i=0; i < row_ndofs; i++)
            {
                elmat(i+k*row_ndofs+row_offset,
                      j+col_offset)
                        += tmp*row_shape(i)*col_shape(j);
            }
    }
}

// Flux Block at position 0,1 Integrator
void mymfem::FluxBlock01Integrator
:: MyAssembleFaceMatrix (const FiniteElement &tr_fe1,
                         const FiniteElement &tr_fe2,
                         const FiniteElement &te_fe1,
                         const FiniteElement &te_fe2,
                         FaceElementTransformations &ftr,
                         DenseMatrix &elmat) const
{
    int tr_ndofs1 = tr_fe1.GetDof();
    int tr_ndofs2 = tr_fe2.GetDof();
    int tr_ndofs = tr_ndofs1 + tr_ndofs2;

    int dim = tr_fe1.GetDim();
    int te_ndofs1 = te_fe1.GetDof();
    int te_ndofs2 = te_fe2.GetDof();
    int te_ndofs = te_ndofs1 + te_ndofs2;
#ifdef MFEM_THREAD_SAFE
    Vector nor(dim);
    Vector tr_shape1(tr_ndofs1), tr_shape2(tr_ndofs2);
    Vector te_shape1(te_ndofs1), te_shape2(te_ndofs2);
#else
    nor.SetSize(dim);
    tr_shape1.SetSize(tr_ndofs1);
    tr_shape2.SetSize(tr_ndofs2);
    te_shape1.SetSize(te_ndofs1);
    te_shape2.SetSize(te_ndofs2);
#endif
    elmat.SetSize(te_ndofs, dim*tr_ndofs);

    const IntegrationRule *ir = IntRule;
    if (ir == nullptr)
    {
        int order = tr_fe1.GetOrder()
                +te_fe1.GetOrder();
        ir = &IntRules.Get(ftr.FaceGeom, order);
    }

    elmat = 0.0;
    for (int i = 0; i < ir->GetNPoints(); i++)
    {
        const IntegrationPoint &ip = ir->IntPoint(i);
        ftr.Face->SetIntPoint(&ip);
        CalcOrtho(ftr.Face->Jacobian(), nor);

        IntegrationPoint eip1, eip2;
        ftr.Loc1.Transform(ip, eip1);
        ftr.Loc2.Transform(ip, eip2);
        tr_fe1.CalcShape(eip1, tr_shape1);
        tr_fe2.CalcShape(eip2, tr_shape2);
        te_fe1.CalcShape(eip1, te_shape1);
        te_fe2.CalcShape(eip2, te_shape2);

        double coeff = 0.5*ip.weight;
        // (1,1) block
        AssembleBlock(dim, te_ndofs1, tr_ndofs1,
                      0, 0,
                      te_shape1, tr_shape1,
                      nor, -coeff, elmat);

        // (1,2) block
        AssembleBlock(dim, te_ndofs1, tr_ndofs2,
                      0, dim*tr_ndofs1,
                      te_shape1, tr_shape2,
                      nor, +coeff, elmat);

        // (2,1) block
        AssembleBlock(dim, te_ndofs2, tr_ndofs1,
                      te_ndofs1, 0,
                      te_shape2, tr_shape1,
                      nor, -coeff, elmat);

        // (2,2) block
        AssembleBlock(dim, te_ndofs2, tr_ndofs2,
                      te_ndofs1, dim*tr_ndofs1,
                      te_shape2, tr_shape2,
                      nor, +coeff, elmat);
    }
}

void mymfem::FluxBlock01Integrator
:: AssembleBlock (const int dim,
                  const int row_ndofs,
                  const int col_ndofs,
                  const int row_offset,
                  const int col_offset,
                  const Vector &row_shape,
                  const Vector &col_shape,
                  const Vector &nor,
                  const double sCoeff,
                  DenseMatrix &elmat) const
{
    for (int k=0; k<dim; k++)
    {
        double tmp = sCoeff*nor(k);
        for (int j=0; j < col_ndofs; j++)
            for (int i=0; i < row_ndofs; i++)
            {
                elmat(i+row_offset,
                      j+k*col_ndofs+col_offset)
                        += tmp*row_shape(i)*col_shape(j);
            }
    }
}

// Dirichlet Block at position 1,0 Integrator
void mymfem::DirichletBlock10Integrator
:: MyAssembleFaceMatrix (const FiniteElement &tr_fe,
                         const FiniteElement &,
                         const FiniteElement &te_fe,
                         const FiniteElement &,
                         FaceElementTransformations &ftr,
                         DenseMatrix &elmat) const
{
    int dim = te_fe.GetDim();
#ifdef MFEM_THREAD_SAFE
    Vector nor;
    Vector tr_shape, te_shape;
#endif
    int tr_ndofs = tr_fe.GetDof();
    int te_ndofs = te_fe.GetDof();

    nor.SetSize(dim);
    tr_shape.SetSize(tr_ndofs);
    te_shape.SetSize(te_ndofs);

    elmat.SetSize(dim*te_ndofs, tr_ndofs);

    const IntegrationRule *ir = IntRule;
    if (ir == nullptr)
    {
        int order = tr_fe.GetOrder()
                +te_fe.GetOrder();
        ir = &IntRules.Get(ftr.FaceGeom, order);
    }

    elmat = 0.0;
    for (int i = 0; i < ir->GetNPoints(); i++)
    {
        const IntegrationPoint &ip = ir->IntPoint(i);
        ftr.Face->SetIntPoint(&ip);
        CalcOrtho(ftr.Face->Jacobian(), nor);

        IntegrationPoint eip;
        ftr.Loc1.Transform(ip, eip);
        tr_fe.CalcShape(eip, tr_shape);
        te_fe.CalcShape(eip, te_shape);

        double coeff = ip.weight;
        AssembleBlock(dim, te_ndofs, tr_ndofs,
                      0, 0,
                      te_shape, tr_shape,
                      nor, -coeff, elmat);
    }
}

void mymfem::DirichletBlock10Integrator
:: AssembleBlock (const int dim,
                  const int row_ndofs,
                  const int col_ndofs,
                  const int row_offset,
                  const int col_offset,
                  const Vector &row_shape,
                  const Vector &col_shape,
                  const Vector &nor,
                  const double sCoeff,
                  DenseMatrix &elmat) const
{
    for (int k=0; k<dim; k++)
    {
        double tmp = sCoeff*nor(k);
        for (int j=0; j < col_ndofs; j++)
            for (int i=0; i < row_ndofs; i++)
            {
                elmat(i+k*row_ndofs+row_offset,
                      j+col_offset)
                        += tmp*row_shape(i)*col_shape(j);
            }
    }
}

// Neumann Block at position 0,1 Integrator
void mymfem::NeumannBlock01Integrator
:: MyAssembleFaceMatrix (const FiniteElement &tr_fe,
                         const FiniteElement &,
                         const FiniteElement &te_fe,
                         const FiniteElement &,
                         FaceElementTransformations &ftr,
                         DenseMatrix &elmat) const
{
    int dim = tr_fe.GetDim();
#ifdef MFEM_THREAD_SAFE
    Vector nor;
    Vector tr_shape, te_shape;
#endif
    int tr_ndofs = tr_fe.GetDof();
    int te_ndofs = te_fe.GetDof();

    nor.SetSize(dim);
    tr_shape.SetSize(tr_ndofs);
    te_shape.SetSize(te_ndofs);

    elmat.SetSize(te_ndofs, dim*tr_ndofs);

    const IntegrationRule *ir = IntRule;
    if (ir == nullptr)
    {
        int order = tr_fe.GetOrder()
                +te_fe.GetOrder();
        ir = &IntRules.Get(ftr.FaceGeom, order);
    }

    elmat = 0.0;
    for (int i = 0; i < ir->GetNPoints(); i++)
    {
        const IntegrationPoint &ip = ir->IntPoint(i);
        ftr.Face->SetIntPoint(&ip);
        CalcOrtho(ftr.Face->Jacobian(), nor);

        IntegrationPoint eip;
        ftr.Loc1.Transform(ip, eip);
        tr_fe.CalcShape(eip, tr_shape);
        te_fe.CalcShape(eip, te_shape);

        double coeff = ip.weight;
        AssembleBlock(dim, te_ndofs, tr_ndofs,
                      0, 0,
                      te_shape, tr_shape,
                      nor, -coeff, elmat);
    }
}

void mymfem::NeumannBlock01Integrator
:: AssembleBlock (const int dim,
                  const int row_ndofs,
                  const int col_ndofs,
                  const int row_offset,
                  const int col_offset,
                  const Vector &row_shape,
                  const Vector &col_shape,
                  const Vector &nor,
                  const double sCoeff,
                  DenseMatrix &elmat) const
{
    for (int k=0; k<dim; k++)
    {
        double tmp = sCoeff*nor(k);
        for (int j=0; j < col_ndofs; j++)
            for (int i=0; i < row_ndofs; i++)
            {
                elmat(i+row_offset,
                      j+k*col_ndofs+col_offset)
                        += tmp*row_shape(i)*col_shape(j);
            }
    }
}

// Dirichlet linear form block at position 0 Integrator
void mymfem::DirichletLFBlock0Integrator
:: AssembleRHSElementVect(const FiniteElement &el,
                          FaceElementTransformations &ftr,
                          Vector &elvect)
{
    int ndofs = el.GetDof();

    shape.SetSize(ndofs);

    elvect.SetSize(ndofs);
    elvect = 0.0;

    const IntegrationRule *ir = IntRule;
    if (ir == nullptr)
    {
        int order = 2*el.GetOrder() + 2;
        ir = &IntRules.Get(ftr.FaceGeom, order);
    }

    for (int i = 0; i < ir->GetNPoints(); i++)
    {
        const IntegrationPoint &ip = ir->IntPoint(i);
        ftr.Face->SetIntPoint(&ip);

        IntegrationPoint eip;
        ftr.Loc1.Transform(ip, eip);
        el.CalcShape(eip, shape);

        double val = ip.weight*Q.Eval(*ftr.Face, ip);
        if (m_stabParamsType == 1 || m_stabParamsType == 3) {
            val *= ftr.Face->Weight();
        }
        elvect.Add(val, shape);
    }
}

// Dirichlet linear form block at position 1 Integrator
void mymfem::DirichletLFBlock1Integrator
:: AssembleRHSElementVect(const FiniteElement &el,
                          FaceElementTransformations &ftr,
                          Vector &elvect)
{
    int dim = el.GetDim();
    int ndofs = el.GetDof();

    nor.SetSize(dim);
    shape.SetSize(ndofs);

    elvect.SetSize(dim*ndofs);
    elvect = 0.0;

    const IntegrationRule *ir = IntRule;
    if (ir == nullptr)
    {
        int order = 2*el.GetOrder() + 2;
        ir = &IntRules.Get(ftr.FaceGeom, order);
    }

    for (int i = 0; i < ir->GetNPoints(); i++)
    {
        const IntegrationPoint &ip = ir->IntPoint(i);
        ftr.Face->SetIntPoint(&ip);
        CalcOrtho(ftr.Face->Jacobian(), nor);

        IntegrationPoint eip;
        ftr.Loc1.Transform(ip, eip);
        el.CalcShape(eip, shape);

        nor *= ip.weight*Q.Eval(*ftr.Face, ip);
        for (int k = 0; k < dim; k++) {
            for (int j = 0; j < ndofs; j++) {
                elvect(j+k*ndofs) += nor(k)*shape(j);
            }
        }
    }
    elvect *= factor;
}

// Neumann linear form block at position 0 Integrator
void mymfem::NeumannLFBlock0Integrator
:: AssembleRHSElementVect(const FiniteElement &el,
                          FaceElementTransformations &ftr,
                          Vector &elvect)
{
    int dim = el.GetDim();
    int ndofs = el.GetDof();

    nor.SetSize(dim);
    shape.SetSize(ndofs);

    elvect.SetSize(ndofs);
    elvect = 0.0;

    const IntegrationRule *ir = IntRule;
    if (ir == nullptr)
    {
        int order = 2*el.GetOrder() + 2;
        ir = &IntRules.Get(ftr.FaceGeom, order);
    }

    for (int i = 0; i < ir->GetNPoints(); i++)
    {
        const IntegrationPoint &ip = ir->IntPoint(i);
        ftr.Face->SetIntPoint(&ip);
        CalcOrtho(ftr.Face->Jacobian(), nor);

        IntegrationPoint eip;
        ftr.Loc1.Transform(ip, eip);
        el.CalcShape(eip, shape);

        Vector v;
        Q.Eval(v, *ftr.Face, ip);
        double val = ip.weight*(v*nor);
        elvect.Add(val, shape);
    }
    elvect *= factor;
}

// Neumann linear form block at position 1 Integrator
void mymfem::NeumannLFBlock1Integrator
:: AssembleRHSElementVect(const FiniteElement &el,
                          FaceElementTransformations &ftr,
                          Vector &elvect)
{
    int dim = el.GetDim();
    int ndofs = el.GetDof();

    nor.SetSize(dim);
    shape.SetSize(ndofs);

    elvect.SetSize(dim*ndofs);
    elvect = 0.0;

    const IntegrationRule *ir = IntRule;
    if (ir == nullptr)
    {
        int order = 2*el.GetOrder() + 2;
        ir = &IntRules.Get(ftr.FaceGeom, order);
    }

    for (int i = 0; i < ir->GetNPoints(); i++)
    {
        const IntegrationPoint &ip = ir->IntPoint(i);
        ftr.Face->SetIntPoint(&ip);
        CalcOrtho(ftr.Face->Jacobian(), nor);
        nor /= ftr.Face->Weight(); // unit normal

        IntegrationPoint eip;
        ftr.Loc1.Transform(ip, eip);
        el.CalcShape(eip, shape);

        Vector v;
        Q.Eval(v, *ftr.Face, ip);
        double coeff = ip.weight*ftr.Face->Weight()*(v*nor);
        if (m_stabParamsType == 3 || m_stabParamsType == 4) {
            coeff *= ftr.Face->Weight();
        }

        for (int k = 0; k < dim; k++) {
            for (int j = 0; j < ndofs; j++) {
                elvect(j + k*ndofs) += coeff*nor(k)*shape(j);
            }
        }
    }
}

// End of file
