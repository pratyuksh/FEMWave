#include "../../include/waveO1/sparse_grids_handler.hpp"
#include "../../include/waveO1/utilities.hpp"

#include <iostream>


//! Computes the DG^{+} error for time-like faces,
//! given a collection of sub-solutions
double SparseGridsHandler
:: eval_dgpFtimeError(Array<GridFunction *> pColl,
                      Array<GridFunction *> vColl,
                      double t) const
{
    double errp=0;
    double errv=0;

    // set exact solution coefficients at time t
    m_pE_coeff->SetTime(t);
    m_vE_coeff->SetTime(t);

    // integration rule for faces
    auto faceElGeomType = m_xMesh_finest->GetFaceBaseGeometry(0);
    int order = 2*m_xW1space_finest->GetOrder(0)+2;
    const IntegrationRule irFace
            = IntRules.Get(faceElGeomType, order);
    int numPoints = irFace.GetNPoints();

    // tranformation operators
    FaceElementTransformations *ftr_finest = nullptr;
    FaceElementTransformations *ftr_coarse = nullptr;
    auto invTr1 = std::make_unique<InverseElementTransformation>();
    auto invTr2 = std::make_unique<InverseElementTransformation>();

    // integration rules for face neighbouring elements
    IntegrationRule irEl1(numPoints);
    IntegrationRule irEl2(numPoints);
    Array<IntegrationPoint> eips1(numPoints);
    Array<IntegrationPoint> eips2(numPoints);
    
    // auxiliary variables
    Vector pSol1(numPoints), pSol2(numPoints);
    Vector pSol1_tmp(numPoints), pSol2_tmp(numPoints);

    DenseMatrix vSol1(m_xndim, numPoints), vSol2(m_xndim, numPoints);
    DenseMatrix vSol1_tmp(m_xndim, numPoints);
    DenseMatrix vSol2_tmp(m_xndim, numPoints);
    Vector buf1, buf2;

    // loop over all interior faces on the finest space mesh,
    // search on coarser meshes, evaluate grid functions
    // and build solution.
    for (int i=0; i<m_xMesh_finest->GetNumFaces(); i++)
    {
        ftr_finest = m_xMesh_finest->GetInteriorFaceTransformations(i);
        if (ftr_finest != nullptr)
        {
            int i1 = ftr_finest->Elem1No;
            int i2 = ftr_finest->Elem2No;

            // integration points for face neighbours
            ftr_finest->Loc1.Transform(irFace, irEl1);
            ftr_finest->Loc2.Transform(irFace, irEl2);

            // use combination formula on the sub-solutions

            // evaluate sub-solution on the finest mesh
            pColl[0]->GetValues(i1, irEl1, pSol1);
            pColl[0]->GetValues(i2, irEl2, pSol2);
            pSol1 *= m_combCoeffs[0];
            pSol2 *= m_combCoeffs[0];

            vColl[0]->GetVectorValues(*ftr_finest->Elem1, irEl1, vSol1);
            vColl[0]->GetVectorValues(*ftr_finest->Elem2, irEl2, vSol2);
            vSol1 *= m_combCoeffs[0];
            vSol2 *= m_combCoeffs[0];

            // physical points to search on coarse meshes
            DenseMatrix points;
            Vector pointsCol;
            ftr_finest->Face->Transform(irFace, points);
            points.GetColumnReference(0, pointsCol);

            // search points and build solution from sub-sols
            Array<bool> boolEdgeFound(m_xMeshFaceLocators.Size());
            boolEdgeFound = false;
            for (int j=0; j<m_xMeshFaceLocators.Size(); j++)
            {
                int j1 = j+1;
                int j2 = j+m_numUniqueLevels;

                // search mesh faces for the physical point
                int faceElId = (*m_xMeshFaceLocators[j])(pointsCol);
                if (faceElId != -1)
                {
                    boolEdgeFound[j] = true;

                    ftr_coarse = m_xMeshes[j1]
                            ->GetFaceElementTransformations(faceElId);
                    invTr1->SetTransformation(*ftr_coarse->Elem1);
                    invTr2->SetTransformation(*ftr_coarse->Elem2);

                    int elId1 = ftr_coarse->Elem1No;
                    int elId2 = ftr_coarse->Elem2No;

                    // find reference points
                    for (int k=0; k<numPoints; k++)
                    {
                        points.GetColumnReference(j, pointsCol);
                        int info1 = invTr1->Transform
                                (pointsCol, eips1[k]);
                        assert(info1 == InverseElementTransformation::
                               TransformResult::Inside);

                        int info2 = invTr2->Transform
                                (pointsCol, eips2[k]);
                        assert(info2 == InverseElementTransformation::
                               TransformResult::Inside);
                    }

                    for (int k=0; k<numPoints; k++)
                    {
                        // eval pSol[j1] - pSol[j2]
                        double val1 = pColl[j1]->GetValue
                                (elId1, eips1[k]);
                        double val2 = pColl[j2]->GetValue
                                (elId1, eips1[k]);
                        pSol1_tmp(k) = m_combCoeffs[j1]*val1
                                + m_combCoeffs[j2]*val2;

                        val1 = pColl[j1]->GetValue
                                (elId2, eips2[k]);
                        val2 = pColl[j2]->GetValue
                                (elId2, eips2[k]);
                        pSol2_tmp(k) = m_combCoeffs[j1]*val1
                                + m_combCoeffs[j2]*val2;

                        // eval vSol[j1] - vSol[j2]
                        vSol1_tmp.GetColumnReference(k, buf1);
                        vColl[j1]->GetVectorValue
                                (elId1, eips1[k], buf1);
                        vColl[j2]->GetVectorValue
                                (elId1, eips1[k], buf2);
                        buf1 *= m_combCoeffs[j1];
                        buf1.Add(m_combCoeffs[j2], buf2);

                        vSol2_tmp.GetColumnReference(k, buf1);
                        vColl[j1]->GetVectorValue
                                (elId2, eips2[k], buf1);
                        vColl[j2]->GetVectorValue
                                (elId2, eips2[k], buf2);
                        buf1 *= m_combCoeffs[j1];
                        buf1.Add(m_combCoeffs[j2], buf2);
                    }
                    pSol1.Add(1, pSol1_tmp);
                    pSol2.Add(1, pSol2_tmp);
                    vSol1.Add(1, vSol1_tmp);
                    vSol2.Add(1, vSol2_tmp);
                }
            }

            // If some points were not found on the coarse mesh faces
            // search point on mesh faces of the finest mesh
            // in coarse mesh elements, and build solution
            Array<int> elIds(numPoints);
            Array<int> init_elIds(m_xPointLocators.Size());
            init_elIds = 0;
            for (int j=0; j<m_xPointLocators.Size(); j++)
            {
                if (boolEdgeFound[j]) { continue; }

                std::tie (elIds, eips1)
                        = (*m_xPointLocators[j])(points,init_elIds[j]);

                int j1 = j+1;
                int j2 = j+m_numUniqueLevels;

                for (int k=0; k<numPoints; k++)
                {
                    // eval pSol[j1] - pSol[j2]
                    double val1 = pColl[j1]->GetValue
                            (elIds[k], eips1[k]);
                    double val2 = pColl[j2]->GetValue
                            (elIds[k], eips1[k]);
                    pSol1_tmp(k) = m_combCoeffs[j1]*val1
                            + m_combCoeffs[j2]*val2;

                    // eval vSol[j1] - vSol[j2]
                    vSol1_tmp.GetColumnReference(k, buf1);
                    vColl[j1]->GetVectorValue
                            (elIds[k], eips1[k], buf1);
                    vColl[j2]->GetVectorValue
                            (elIds[k], eips1[k], buf2);
                    buf1 *= m_combCoeffs[j1];
                    buf1.Add(m_combCoeffs[j2], buf2);
                }
                pSol1.Add(1, pSol1_tmp);
                pSol2.Add(1, pSol1_tmp);
                vSol1.Add(1, vSol1_tmp);
                vSol2.Add(1, vSol1_tmp);
            }

            // evaluate errors
            double coeff1, coeff2;
            Vector nor(m_xndim);
            Vector errvLocal;
            for (int j=0; j<numPoints; j++)
            {
                const IntegrationPoint &ip = irFace.IntPoint(j);
                ftr_finest->Face->SetIntPoint(&ip);

                // pressure averages
                coeff1 = ip.weight;
                if (m_stabParamsType == 1
                        || m_stabParamsType == 2) {
                    coeff1 *= ftr_finest->Face->Weight();
                }
                double pE = m_pE_coeff->Eval(*ftr_finest->Face, ip);
                double val = pE - 0.5*(pSol1(j) + pSol2(j));
                errp += coeff1*val*val;

                // velocity averages
                CalcOrtho(ftr_finest->Face->Jacobian(), nor);
                nor /= ftr_finest->Face->Weight(); // unit normal

                coeff2 = ip.weight*ftr_finest->Face->Weight();
                if (m_stabParamsType == 2
                        || m_stabParamsType == 4) {
                    coeff2 *= ftr_finest->Face->Weight();
                }
                m_vE_coeff->Eval(errvLocal, *ftr_finest->Face, ip);
                Vector tmp1(vSol1.GetColumn(j), m_xndim);
                Vector tmp2(vSol2.GetColumn(j), m_xndim);
                errvLocal.Add(-0.5, tmp1);
                errvLocal.Add(-0.5, tmp2);
                errv += coeff2*(errvLocal*nor)*(errvLocal*nor);
            }
        }
    }

    // boundary markers
    Array<int> xEss_bdr_marker
            = m_discrs[0]->get_xEss_bdr_marker();
    Array<int> xNat_bdr_marker
            = m_discrs[0]->get_xNat_bdr_marker();

    // loop over all boundary faces on the finest space mesh,
    // search on coarser meshes, evaluate grid functions
    // and build solution.
    for (int i=0; i<m_xMesh_finest->GetNBE(); i++)
    {
        const int bdr_attr = m_xMesh_finest->GetBdrAttribute(i);
        ftr_finest = m_xMesh_finest -> GetBdrFaceTransformations(i);

        if (xNat_bdr_marker[bdr_attr-1] == 1) // Neumann bdry
        {
            int i1 = ftr_finest->Elem1No;

            // element integration points
            ftr_finest->Loc1.Transform(irFace, irEl1);

            // use combination formula on the sub-solutions

            // evaluate sub-solution on the finest mesh
            pColl[0]->GetValues(i1, irEl1, pSol1);
            pSol1 *= m_combCoeffs[0];

            // physical points to search on coarse meshes
            DenseMatrix points;
            Vector pointsCol;
            ftr_finest->Face->Transform(irFace, points);
            points.GetColumnReference(0, pointsCol);

            // search points and build solution from sub-sols
            for (int j=0; j<m_xMeshFaceLocators.Size(); j++)
            {
                int j1 = j+1;
                int j2 = j+m_numUniqueLevels;

                // search mesh faces for the physical point
                int faceElId = (*m_xMeshFaceLocators[j])(pointsCol);
                if (faceElId != -1)
                {
                    ftr_coarse = m_xMeshes[j1]
                            ->GetFaceElementTransformations(faceElId);
                    invTr1->SetTransformation(*ftr_coarse->Elem1);

                    int elId1 = ftr_coarse->Elem1No;

                    // find reference points
                    for (int k=0; k<numPoints; k++)
                    {
                        points.GetColumnReference(j, pointsCol);
                        int info1 = invTr1->Transform
                                (pointsCol, eips1[k]);
                        assert(info1 == InverseElementTransformation::
                               TransformResult::Inside);
                    }

                    for (int k=0; k<numPoints; k++)
                    {
                        // eval pSol[j1] - pSol[j2]
                        double val1 = pColl[j1]->GetValue
                                (elId1, eips1[k]);
                        double val2 = pColl[j2]->GetValue
                                (elId1, eips1[k]);
                        pSol1_tmp(k) = m_combCoeffs[j1]*val1
                                + m_combCoeffs[j2]*val2;
                    }
                    pSol1.Add(1, pSol1_tmp);
                }
            }

            // evaluate errors
            double coeff;
            Vector nor(m_xndim);
            for (int j=0; j<numPoints; j++)
            {
                const IntegrationPoint &ip = irFace.IntPoint(j);
                ftr_finest->Face->SetIntPoint(&ip);

                coeff = ip.weight;
                if (m_stabParamsType == 1
                        || m_stabParamsType == 2) {
                    coeff *= ftr_finest->Face->Weight();
                }

                // pressure averages
                double pE = m_pE_coeff->Eval(*ftr_finest->Face, ip);
                errp += coeff*(pE - pSol1(j))*(pE - pSol1(j));
            }
        }
        else // Dirichlet bdry
        {
            // element integration points
            ftr_finest->Loc1.Transform(irFace, irEl1);

            // use combination formula on the sub-solutions

            // evaluate sub-solution on the finest mesh
            vColl[0]->GetVectorValues(*ftr_finest->Elem1, irEl1, vSol1);
            vSol1 *= m_combCoeffs[0];

            // physical points to search on coarse meshes
            DenseMatrix points;
            Vector pointsCol;
            ftr_finest->Face->Transform(irFace, points);
            points.GetColumnReference(0, pointsCol);

            // search points and build solution from sub-sols
            for (int j=0; j<m_xMeshFaceLocators.Size(); j++)
            {
                int j1 = j+1;
                int j2 = j+m_numUniqueLevels;

                // search mesh faces for the physical point
                int faceElId = (*m_xMeshFaceLocators[j])(pointsCol);
                if (faceElId != -1)
                {
                    ftr_coarse = m_xMeshes[j1]
                            ->GetFaceElementTransformations(faceElId);
                    invTr1->SetTransformation(*ftr_coarse->Elem1);

                    int elId1 = ftr_coarse->Elem1No;

                    // find reference points
                    for (int k=0; k<numPoints; k++)
                    {
                        points.GetColumnReference(j, pointsCol);
                        int info1 = invTr1->Transform
                                (pointsCol, eips1[k]);
                        assert(info1 == InverseElementTransformation::
                               TransformResult::Inside);
                    }

                    for (int k=0; k<numPoints; k++)
                    {
                        // eval pSol[j1] - pSol[j2]
                        vSol1_tmp.GetColumnReference(k, buf1);
                        vColl[j1]->GetVectorValue
                                (elId1, eips1[k], buf1);
                        vColl[j2]->GetVectorValue
                                (elId1, eips1[k], buf2);
                        buf1 *= m_combCoeffs[j1];
                        buf1.Add(m_combCoeffs[j2], buf2);
                    }
                    vSol1.Add(1, vSol1_tmp);
                }
            }

            // evaluate errors
            double coeff;
            Vector nor(m_xndim);
            Vector errvLocal(m_xndim);
            for (int j=0; j<numPoints; j++)
            {
                const IntegrationPoint &ip = irFace.IntPoint(j);
                ftr_finest->Face->SetIntPoint(&ip);

                CalcOrtho(ftr_finest->Face->Jacobian(), nor);
                nor /= ftr_finest->Face->Weight(); // unit normal

                coeff = ip.weight*ftr_finest->Face->Weight();
                if (m_stabParamsType == 2
                        || m_stabParamsType == 4) {
                    coeff *= ftr_finest->Face->Weight();
                }

                // velocity averages
                m_vE_coeff->Eval(errvLocal, *ftr_finest->Face, ip);
                Vector tmp(vSol1.GetColumn(j), m_xndim);
                errvLocal -= tmp;
                errv += coeff*(errvLocal*nor)*(errvLocal*nor);
            }
        }
    }
    //std::cout << errp << "\t" << errv << std::endl;
    return (errp + errv);
}

//! Computes the DG error for space-like faces,
//! given a collection of sub-solutions
double SparseGridsHandler
:: eval_dgpFspaceError
(Array<GridFunction *> pColl, Array<GridFunction *> vColl,
 double t, std::shared_ptr<WaveO1InvSqMediumCoeff>& invSqMed) const
{
    double errp=0;
    double errv=0;

    // set exact solution coefficients at time t
    m_pE_coeff->SetTime(t);
    m_vE_coeff->SetTime(t);

    // transformation operators
    ElementTransformation *trans_finest = nullptr;
    ElementTransformation *trans_coarse = nullptr;

    const FiniteElement *fe = nullptr;
    Array<int> init_elIds(m_xPointLocators.Size());
    init_elIds = 0;

    // auxiliary variables
    Vector pSol, pSol_tmp;
    DenseMatrix vSol, vSol_tmp;
    Vector buf1, buf2;

    // loop over all elements on the finest space mesh,
    // search on coarser meshes, evaluate grid functions
    // and build solution.
    for (int i=0; i<m_xMesh_finest->GetNE(); i++)
    {
        fe = m_xW1space_finest->GetFE(i);
        int order = 2*fe->GetOrder()+2;
        const IntegrationRule *ir
                = &IntRules.Get(fe->GetGeomType(), order);
        int numPoints = ir->GetNPoints();
        trans_finest = m_xW1space_finest->
                GetElementTransformation(i);

        // use combination formula on the sub-solutions

        // evaluate sub-solution on the finest mesh
        pColl[0]->GetValues(i, *ir, pSol);
        vColl[0]->GetVectorValues(*trans_finest, *ir, vSol);
        pSol *= m_combCoeffs[0];
        vSol *= m_combCoeffs[0];

        pSol_tmp.SetSize(pSol.Size());
        vSol_tmp.SetSize(vSol.NumRows(), vSol.NumCols());

        // physical points to search on coarse meshes
        DenseMatrix points;
        trans_finest->Transform(*ir, points);

        // search points and build solution from sub-sols
        Array<int> elIds(numPoints);
        Array<IntegrationPoint> ips(numPoints);
        for (int j=0; j<m_xPointLocators.Size(); j++)
        {
            std::tie (elIds,ips)
                    = (*m_xPointLocators[j])(points, init_elIds[j]);

            int j1 = j+1;
            int j2 = j+m_numUniqueLevels;
            for (int k=0; k<ips.Size(); k++)
            {
                // eval pSol[j1] - pSol[j2]
                double val1 = pColl[j1]->GetValue
                        (elIds[k], ips[k]);
                double val2 = pColl[j2]->GetValue
                        (elIds[k], ips[k]);
                pSol_tmp(k) = m_combCoeffs[j1]*val1
                            + m_combCoeffs[j2]*val2;

                // eval vSol[j1] - vSol[j2]
                trans_coarse = m_xMeshes[j1]
                        ->GetElementTransformation(elIds[k]);
                trans_coarse->SetIntPoint(&ips[k]);

                vSol_tmp.GetColumnReference(k, buf1);
                vColl[j1]->GetVectorValue(elIds[k], ips[k],
                                          buf1);
                vColl[j2]->GetVectorValue(elIds[k], ips[k],
                                          buf2);
                buf1 *= m_combCoeffs[j1];
                buf1.Add(m_combCoeffs[j2], buf2);
            }
            pSol.Add(1, pSol_tmp);
            vSol.Add(1, vSol_tmp);
        }

        // evaluate errors
        double pE, w;
        Vector errvLocal(m_xndim);
        for (int j=0; j<numPoints; j++)
        {
            const IntegrationPoint &ip = ir->IntPoint(j);
            trans_finest->SetIntPoint(&ip);
            w = ip.weight*trans_finest->Weight();
            pE = m_pE_coeff->Eval(*trans_finest, ip);
            m_vE_coeff->Eval(errvLocal, *trans_finest, ip);

            // pressure trace at t_n^{-}
            errp += w*(pE - pSol(j))*(pE - pSol(j))
                    *invSqMed->Eval(*trans_finest, ip);

            // velocity trace at t_n^{-}
            Vector tmp(vSol.GetColumn(j), m_xndim);
            errvLocal -= tmp;
            errv += w*(errvLocal*errvLocal);
        }
    }

    return (errp + errv);
}

// Computes DG^{+} error
std::tuple <double, double>
SparseGridsHandler :: eval_xtDgpError() const
{
    double errDgpFtime=0;
    double errDgpFspace=0;

    // medium
    auto invSqMed = std::make_shared<WaveO1InvSqMediumCoeff>
            (m_testCase);

    // allocate memory for collection of solutions
    // pColl wraps pSubSols as GF
    int Nt = m_tWspace_finest->GetNE();
    Array<Vector *> pSubSols(m_numSubSols);
    Array<Vector *> vSubSols(m_numSubSols);
    Array<GridFunction *> pColl(m_numSubSols);
    Array<GridFunction *> vColl(m_numSubSols);
    for (int k=0; k<m_numSubSols; k++)
    {
        auto xW1space = (m_discrs[k]->get_xFespaces())[0];
        auto xW2space = (m_discrs[k]->get_xFespaces())[1];
        int xdimW1 = xW1space->GetTrueVSize();
        int xdimW2 = xW2space->GetTrueVSize();
        pSubSols[k] = new Vector(xdimW1);
        vSubSols[k] = new Vector(xdimW2);
        pColl[k] = new GridFunction(xW1space, *pSubSols[k]);
        vColl[k] = new GridFunction(xW2space, *vSubSols[k]);
    }

    Eigen::VectorXd bufErrDgpFtime(Nt);
    bufErrDgpFtime.setZero();

    IntegrationPoint tIpC;
    Array<int> tElIdsC(m_tPointLocators.Size());
    Array<int> init_tElIdsC(m_tPointLocators.Size());

    int idf = m_numUniqueLevels-1;
    init_tElIdsC = 0;

    Vector t, tShape;
    Array<int> tVdofsF, tVdofsC;
    ElementTransformation *tTransF = nullptr;
    const FiniteElement *tFeF = nullptr;
    const FiniteElement *tFeC = nullptr;

    // DG error for time-like faces
    for (int n=0; n<Nt; n++)
    {
        tFeF = m_tWspace_finest->GetFE(n);
        tTransF = m_tWspace_finest->
                GetElementTransformation(n);
        m_tWspace_finest->GetElementVDofs(n, tVdofsF);

        // integration points
        int order = 2*tFeF->GetOrder()+2;
        const IntegrationRule *ir
                = &IntRules.Get(tFeF->GetGeomType(), order);

        tShape.SetSize(tFeF->GetDof());
        for (int i = 0; i < ir->GetNPoints(); i++)
        {
            const IntegrationPoint &tIpf = ir->IntPoint(i);
            tTransF->SetIntPoint(&tIpf);
            tFeF->CalcShape(tIpf, tShape);
            tTransF->Transform(tIpf, t);
            //std::cout << n << "\t" << i << "\t"
            //          << t(0) << std::endl;

            // set uSubSols at time t
            build_xSol_FG(m_subSols[idf]->GetBlock(0),
                          tShape, tVdofsF, *pSubSols[idf]);
            build_xSol_FG(m_subSols[idf]->GetBlock(1),
                          tShape, tVdofsF, *vSubSols[idf]);
            for (int j=0; j<m_tPointLocators.Size(); j++)
            {
                std::tie (tElIdsC[j], tIpC)
                        = (*m_tPointLocators[j])
                        (t, init_tElIdsC[j]);
                tFeC = m_discrs[j]->
                        get_tFespace()->GetFE(tElIdsC[j]);
                tFeC->CalcShape(tIpC, tShape);
                m_discrs[j]->get_tFespace()
                        ->GetElementVDofs
                        (tElIdsC[j], tVdofsC);

                int j1 = j;
                int j2 = j+m_numUniqueLevels;
                build_xSol_FG(m_subSols[j1]->GetBlock(0),
                              tShape, tVdofsC,
                              *pSubSols[j1]);
                build_xSol_FG(m_subSols[j2]->GetBlock(0),
                              tShape, tVdofsC,
                              *pSubSols[j2]);
                build_xSol_FG(m_subSols[j1]->GetBlock(1),
                              tShape, tVdofsC,
                              *vSubSols[j1]);
                build_xSol_FG(m_subSols[j2]->GetBlock(1),
                              tShape, tVdofsC,
                              *vSubSols[j2]);
            }

            // eval error at time t
            double errFtime = eval_dgpFtimeError(pColl, vColl, t(0));
            double w = tIpf.weight*tTransF->Weight();
            bufErrDgpFtime(n) += w*errFtime;
        }
    }
    errDgpFtime = std::sqrt(bufErrDgpFtime.sum());

    // DG error for space-like faces
    Eigen::VectorXd bufErrDgpFspace(Nt);
    bufErrDgpFspace.setZero();

    IntegrationPoint ip0;
    ip0.Set1w(0, 1);

    for (int n=0; n<Nt; n++)
    {
        int np = n;

        // build solution at tn^{+}
        m_tWspace_finest->GetElementVDofs(np, tVdofsF);
        tTransF = m_tWspace_finest->GetElementTransformation(np);
        tFeF = m_tWspace_finest->GetFE(np);
        int tNdofsF = tFeF->GetDof();
        tShape.SetSize(tNdofsF);

        tTransF->SetIntPoint(&ip0);
        tFeF->CalcShape(ip0, tShape);
        tTransF->Transform(ip0, t);

        build_xSol_FG(m_subSols[idf]->GetBlock(0),
                      tShape, tVdofsF, *pSubSols[idf]);
        build_xSol_FG(m_subSols[idf]->GetBlock(1),
                      tShape, tVdofsF, *vSubSols[idf]);
        for (int j=0; j<m_tPointLocators.Size(); j++)
        {
            std::tie (tElIdsC[j], tIpC)
                    = (*m_tPointLocators[j])
                    (t, init_tElIdsC[j]);
            tFeC = m_discrs[j]->
                    get_tFespace()->GetFE(tElIdsC[j]);
            tFeC->CalcShape(tIpC, tShape);
            m_discrs[j]->get_tFespace()
                    ->GetElementVDofs(tElIdsC[j], tVdofsC);

            int j1 = j;
            int j2 = j+m_numUniqueLevels;
            build_xSol_FG(m_subSols[j1]->GetBlock(0),
                          tShape, tVdofsC,
                          *pSubSols[j1]);
            build_xSol_FG(m_subSols[j2]->GetBlock(0),
                          tShape, tVdofsC,
                          *pSubSols[j2]);
            build_xSol_FG(m_subSols[j1]->GetBlock(1),
                          tShape, tVdofsC,
                          *vSubSols[j1]);
            build_xSol_FG(m_subSols[j2]->GetBlock(1),
                          tShape, tVdofsC,
                          *vSubSols[j2]);
        }

        // eval error
        bufErrDgpFspace(n)
                = eval_dgpFspaceError(pColl, vColl, t(0), invSqMed);
    }
    errDgpFspace = std::sqrt(2*bufErrDgpFspace.sum());

    // free memory
    for (int k=0; k<m_numSubSols; k++) {
        delete pSubSols[k];
        delete vSubSols[k];
        delete pColl[k];
        delete vColl[k];
    }

    return {errDgpFtime, errDgpFspace};
}

// End of file
