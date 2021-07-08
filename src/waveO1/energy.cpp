#include "../../include/waveO1/observer.hpp"
#include "../../include/waveO1/utilities.hpp"
#include "../../include/mymfem/utilities.hpp"

#include <iostream>
#include <Eigen/Core>


//! Computes DG jumps on time-like face at a specified time
double WaveO1Observer
:: eval_dgFtimeJump (std::shared_ptr<GridFunction>& p,
                     std::shared_ptr<GridFunction>& v) const
{
    double jumpp=0;
    double jumpv=0;

    // assumes that all elements have the same geometry
    auto faceElGeomType = p->FESpace()->GetMesh()
                ->GetFaceGeometryType(0);

    // jump in pressure
    // set integration rules
    IntegrationRules rule1{};
    const IntegrationRule *ir1;
    int order1 = 2*p->FESpace()->GetFE(0)->GetOrder()+2;
    ir1 = &rule1.Get(faceElGeomType, order1);
    m_assembleDgJumps->SetIntegrationRule(ir1);
    jumpp = m_assembleDgJumps->FtimeScalar(p.get());

    // jump in velocity
    // set integration rules
    IntegrationRules rule2{};
    const IntegrationRule *ir2;
    int order2 = 2*v->FESpace()->GetFE(0)->GetOrder()+2;
    ir2 = &rule2.Get(faceElGeomType, order2);
    m_assembleDgJumps->SetIntegrationRule(ir2);
    jumpv = m_assembleDgJumps->FtimeVector(v.get());

    return (jumpp + jumpv);
}

//! Computes DG jump for a space-like face
//! at a specified time
double WaveO1Observer
:: eval_dgFspaceJump (int degp, int degv, Mesh* mesh,
                      double t, std::shared_ptr
                      <WaveO1InvSqMediumCoeff>& invSqMed)
const
{
    double tracep=0;
    double tracev=0;

    // assumes that all elements have the same geometry
    auto elGeomType = mesh->GetElementBaseGeometry(0);

    IntegrationRules rule1{}, rule2{};
    const IntegrationRule *ir1, *ir2;

    // pressure
    int order1 = 2*degp+4;
    ir1 = &rule1.Get(elGeomType, order1);
    m_assembleDgJumps->SetIntegrationRule(ir1);
    m_pE_coeff->SetTime(t);
    tracep = m_assembleDgJumps->FspaceScalar
            (mesh, m_pE_coeff.get(), invSqMed.get());

    // velocity
    int order2 = 2*degv+4;
    ir2 = &rule2.Get(elGeomType, order2);
    m_assembleDgJumps->SetIntegrationRule(ir2);
    tracev = m_assembleDgJumps->FspaceVector
            (mesh, m_vE_coeff.get());

    return (tracep + tracev);
}

double WaveO1Observer
:: eval_dgFspaceJump (std::shared_ptr<GridFunction>& p,
                      std::shared_ptr<GridFunction>& v,
                      std::shared_ptr
                      <WaveO1InvSqMediumCoeff>& invSqMed)
const
{
    double jumpp=0;
    double jumpv=0;

    // assumes that all elements have the same geometry
    auto elGeomType = p->FESpace()
                ->GetFE(0)->GetGeomType();

    IntegrationRules rule1{}, rule2{};
    const IntegrationRule *ir1, *ir2;

    // pressure
    int order1 = 2*p->FESpace()->GetFE(0)->GetOrder()+2;
    ir1 = &rule1.Get(elGeomType, order1);
    m_assembleDgJumps->SetIntegrationRule(ir1);
    jumpp = m_assembleDgJumps->FspaceScalar
            (p.get(), invSqMed.get());

    // velocity
    int order2 = 2*v->FESpace()->GetFE(0)->GetOrder()+2;
    ir2 = &rule2.Get(elGeomType, order2);
    m_assembleDgJumps->SetIntegrationRule(ir2);
    jumpv = m_assembleDgJumps->FspaceVector(v.get());

    return (jumpp + jumpv);
}

double WaveO1Observer
:: eval_dgFspaceJump (std::shared_ptr<GridFunction>& p1,
                      std::shared_ptr<GridFunction>& p2,
                      std::shared_ptr<GridFunction>& v1,
                      std::shared_ptr<GridFunction>& v2,
                      std::shared_ptr
                      <WaveO1InvSqMediumCoeff>& invSqMed)
const
{
    double jumpp=0;
    double jumpv=0;

    // assumes that all elements have the same geometry
    auto elGeomType = p1->FESpace()
                ->GetFE(0)->GetGeomType();

    IntegrationRules rule1{}, rule2{};
    const IntegrationRule *ir1, *ir2;

    // pressure
    int order1 = 2*p1->FESpace()->GetFE(0)->GetOrder()+2;
    ir1 = &rule1.Get(elGeomType, order1);
    m_assembleDgJumps->SetIntegrationRule(ir1);
    jumpp = m_assembleDgJumps->FspaceScalar
            (p1.get(), p2.get(), invSqMed.get());

    // velocity
    int order2 = 2*v1->FESpace()->GetFE(0)->GetOrder()+2;
    ir2 = &rule2.Get(elGeomType, order2);
    m_assembleDgJumps->SetIntegrationRule(ir2);
    jumpv = m_assembleDgJumps->FspaceVector
            (v1.get(), v2.get());

    return (jumpp + jumpv);
}

// Computes energy
Vector WaveO1Observer
:: eval_energy(BlockVector& W) const
{
    int Nt = m_tWspace->GetNE();
    int xdimW1 = m_xW1space->GetTrueVSize();
    int xdimW2 = m_xW2space->GetTrueVSize();

    Vector energy(Nt+1);

    Vector& p = W.GetBlock(0);
    Vector& v = W.GetBlock(1);

    Vector pSol1(xdimW1), pSol2(xdimW1);
    auto pGSol1 = std::make_shared<GridFunction>
            (m_xW1space, pSol1);
    auto pGSol2 = std::make_shared<GridFunction>
            (m_xW1space, pSol2);

    Vector vSol1(xdimW2), vSol2(xdimW2);
    auto vGSol1 = std::make_shared<GridFunction>
            (m_xW2space, vSol1);
    auto vGSol2 = std::make_shared<GridFunction>
            (m_xW2space, vSol2);

    // medium
    auto invSqMed = std::make_shared
            <WaveO1InvSqMediumCoeff>(m_testCase);

    // energy @t=0
    int degp = pGSol1->FESpace()->GetFE(0)->GetOrder();
    int degv = vGSol1->FESpace()->GetFE(0)->GetOrder();
    double initEnergy = 0.5*eval_dgFspaceJump
            (degp, degv, pGSol1->FESpace()->GetMesh(),
             0, invSqMed);

    ElementTransformation *tTrans = nullptr;
    const FiniteElement *tFe = nullptr;
    Vector tShape;
    Array<int> tVdofs;

    // for time-like faces
    Eigen::VectorXd bufEnergyFtime(Nt);
    bufEnergyFtime.setZero();

    double jumpFtime;
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
            build_xSol_FG(p, tShape, tVdofs, pSol1);
            build_xSol_FG(v, tShape, tVdofs, vSol1);

            jumpFtime = eval_dgFtimeJump(pGSol1, vGSol1);

            double w = ip.weight*tTrans->Weight();
            bufEnergyFtime(n) += w*jumpFtime;
        }
    }

    // for space-like faces
    Eigen::VectorXd bufEnergyFspace(Nt-1);
    bufEnergyFspace.setZero();

    IntegrationPoint ip0, ip1;
    ip0.Set1w(0, 1);
    ip1.Set1w(1, 1);

    // for 0 < t < T
    for (int n=1; n<Nt; n++)
    {
        int nm = n-1;
        int np = n;

        // build solution at tn^{-}
        m_tWspace->GetElementVDofs(nm, tVdofs);
        tTrans = m_tWspace->GetElementTransformation(nm);
        tFe = m_tWspace->GetFE(nm);
        int tNdofs = tFe->GetDof();
        tShape.SetSize(tNdofs);

        tTrans->SetIntPoint(&ip1);
        tFe->CalcShape(ip1, tShape);
        build_xSol_FG(p, tShape, tVdofs, pSol1);
        build_xSol_FG(v, tShape, tVdofs, vSol1);

        // build solution at tn^{+}
        m_tWspace->GetElementVDofs(np, tVdofs);
        tTrans = m_tWspace->GetElementTransformation(np);
        tFe = m_tWspace->GetFE(np);
        tNdofs = tFe->GetDof();
        tShape.SetSize(tNdofs);

        tTrans->SetIntPoint(&ip0);
        tFe->CalcShape(ip0, tShape);
        build_xSol_FG(p, tShape, tVdofs, pSol2);
        build_xSol_FG(v, tShape, tVdofs, vSol2);

        bufEnergyFspace(n-1) = 0.5*eval_dgFspaceJump
                (pGSol1, pGSol2, vGSol1, vGSol2, invSqMed);
    }

    energy(0) = initEnergy;
    energy(1) = energy(0)
            - bufEnergyFtime(0);
    for (int n=2; n<=Nt; n++) {
        energy(n) = energy(n-1)
                - bufEnergyFtime(n-1) - bufEnergyFspace(n-2);
    }

    return energy;
}

// End of file
