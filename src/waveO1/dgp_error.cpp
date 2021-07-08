#include "../../include/waveO1/observer.hpp"
#include "../../include/waveO1/utilities.hpp"
#include "../../include/mymfem/utilities.hpp"

#include <fstream>
#include <iostream>
#include <Eigen/Core>


double WaveO1Observer
:: eval_dgpFtimeError (std::shared_ptr<GridFunction>& p,
                      std::shared_ptr<GridFunction>& v,
                      double t) const
{
    double errp=0;
    double errv=0;

    // assumes that all elements have the same geometry
    auto faceElGeomType = p->FESpace()->GetMesh()
                ->GetFaceGeometryType(0);

    // L2 error in pressure
    // set integration rules
    IntegrationRules rule1{};
    const IntegrationRule *ir1;
    int order1 = 2*p->FESpace()->GetFE(0)->GetOrder()+2;
    ir1 = &rule1.Get(faceElGeomType, order1);
    m_pE_coeff->SetTime(t);
    m_assembleDgpSmoothError->SetIntegrationRule(ir1);
    errp = m_assembleDgpSmoothError->Ftime(p.get(), m_pE_coeff.get());

    // L2 error in velocity
    // set integration rules
    IntegrationRules rule2{};
    const IntegrationRule *ir2;
    int order2 = 2*v->FESpace()->GetFE(0)->GetOrder()+2;
    ir2 = &rule2.Get(faceElGeomType, order2);
    m_vE_coeff->SetTime(t);
    m_assembleDgpSmoothError->SetIntegrationRule(ir2);
    errv = m_assembleDgpSmoothError->Ftime(v.get(), m_vE_coeff.get());

    return (errp + errv);
}

//! Computes DG error for a space-like face
//! at a specified time
double WaveO1Observer
:: eval_dgpFspaceError
(std::shared_ptr<GridFunction>& p,
 std::shared_ptr<GridFunction>& v,
 double t,
 std::shared_ptr<WaveO1InvSqMediumCoeff>&) const
{
    double errp=0;
    double errv=0;

    // assumes that all elements have the same geometry
    auto elGeomType = p->FESpace()
                ->GetFE(0)->GetGeomType();

    IntegrationRules rule1{}, rule2{};
    const IntegrationRule *ir1, *ir2;

    m_pE_coeff->SetTime(t);
    m_vE_coeff->SetTime(t);

    // pressure
    int order1 = 2*p->FESpace()->GetFE(0)->GetOrder()+2;
    ir1 = &rule1.Get(elGeomType, order1);
    m_assembleDgpSmoothError->SetIntegrationRule(ir1);
    errp = m_assembleDgpSmoothError->Fspace
            (p.get(), m_pE_coeff.get());
            //(p.get(), m_pE_coeff.get(), invSqMed.get());

    // velocity
    int order2 = 2*v->FESpace()->GetFE(0)->GetOrder()+2;
    ir2 = &rule2.Get(elGeomType, order2);
    m_assembleDgpSmoothError->SetIntegrationRule(ir2);
    m_vE->ProjectCoefficient(*m_vE_coeff);
    errv = m_assembleDgpSmoothError->Fspace(v.get(), m_vE_coeff.get());

    return (errp + errv);
}

//! Computes DG^{+} error
std::tuple <double, double> WaveO1Observer
:: eval_xtDgpError (BlockVector& W) const
{
    double errDgFtime=0;
    double errDgFspace=0;

    int Nt = m_tWspace->GetNE();
    int xdimW1 = m_xW1space->GetTrueVSize();
    int xdimW2 = m_xW2space->GetTrueVSize();

    Vector& p = W.GetBlock(0);
    Vector& v = W.GetBlock(1);

    Vector pSol(xdimW1);
    std::shared_ptr<GridFunction> pGSol
            = std::make_shared<GridFunction>
            (m_xW1space, pSol);

    Vector vSol(xdimW2);
    std::shared_ptr<GridFunction> vGSol
            = std::make_shared<GridFunction>
            (m_xW2space, vSol);

    // medium
    auto invSqMed = std::make_shared<WaveO1InvSqMediumCoeff>
            (m_testCase);

    ElementTransformation *tTrans = nullptr;
    const FiniteElement *tFe = nullptr;
    Vector tShape;
    Array<int> tVdofs;

    // DG^{+} error for time-like faces
    Eigen::VectorXd bufErrFtime(Nt);
    bufErrFtime.setZero();

    double errFtime;
    for (int n=0; n<Nt; n++)
    {
        m_tWspace->GetElementVDofs(n, tVdofs);
        tTrans = m_tWspace->GetElementTransformation(n);

        tFe = m_tWspace->GetFE(n);
        int tNdofs = tFe->GetDof();
        tShape.SetSize(tNdofs);

        int order = 2*tFe->GetOrder()+2;
        const IntegrationRule *ir
                = &IntRules.Get(tFe->GetGeomType(), order);

        for (int i = 0; i < ir->GetNPoints(); i++)
        {
            const IntegrationPoint &ip = ir->IntPoint(i);
            tTrans->SetIntPoint(&ip);
            tFe->CalcShape(ip, tShape);

            // build solution at time t
            Vector t;
            tTrans->Transform(ip, t);
            build_xSol_FG(p, tShape, tVdofs, pSol);
            build_xSol_FG(v, tShape, tVdofs, vSol);

            errFtime = eval_dgpFtimeError(pGSol, vGSol, t(0));

            double w = ip.weight*tTrans->Weight();
            bufErrFtime(n) += w*errFtime;
        }
    }
    errDgFtime = std::sqrt(bufErrFtime.sum());

    // DG^{+} error for space-like faces
    Eigen::VectorXd bufErrFspace(Nt);
    bufErrFspace.setZero();

    Vector t;
    IntegrationPoint ip0;
    ip0.Set1w(0, 1);

    for (int n=0; n<Nt; n++)
    {
        int np = n;

        // build solution at tn^{+}
        m_tWspace->GetElementVDofs(np, tVdofs);
        tTrans = m_tWspace->GetElementTransformation(np);
        tFe = m_tWspace->GetFE(np);
        int tNdofs = tFe->GetDof();
        tShape.SetSize(tNdofs);

        tTrans->SetIntPoint(&ip0);
        tFe->CalcShape(ip0, tShape);
        tTrans->Transform(ip0, t);
        build_xSol_FG(p, tShape, tVdofs, pSol);
        build_xSol_FG(v, tShape, tVdofs, vSol);

        // compute error
        bufErrFspace(n) = eval_dgpFspaceError
                (pGSol, vGSol, t(0), invSqMed);
    }
    errDgFspace = std::sqrt(2*bufErrFspace.sum());

    return {errDgFtime, errDgFspace};
}

// End of file
