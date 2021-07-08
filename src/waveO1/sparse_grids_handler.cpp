#include "../../include/waveO1/sparse_grids_handler.hpp"
#include "../../include/waveO1/utilities.hpp"

#include <iostream>
using namespace std;


// Constructors
SparseGridsHandler
:: SparseGridsHandler (const nlohmann::json& config)
    : m_config(config)
{
    m_endTime = config["end_time"];

    m_stabParamsType = 1;
    if (config.contains("stab_parameters_type")) {
        m_stabParamsType = config["stab_parameters_type"];
    }

    m_bool_error = false;
    if (config.contains("eval_error")) {
        m_bool_error = config["eval_error"];
    }

    m_errorType = "relL2";
    if (config.contains("error_type")) {
        m_errorType = config["error_type"];
    }
}

SparseGridsHandler
:: SparseGridsHandler
(const nlohmann::json& config,
 std::shared_ptr<WaveO1TestCases>& testCase)
    : SparseGridsHandler (config)
{
    m_testCase = testCase;
}

SparseGridsHandler
:: SparseGridsHandler (const nlohmann::json& config,
                       std::string mesh_dir,
                       const int Lx,
                       const int L0x,
                       const int L0t,
                       bool load_init_mesh)
    : SparseGridsHandler(config)
{
    set(mesh_dir, load_init_mesh);
    set(Lx, L0x, L0t);
}

SparseGridsHandler
:: SparseGridsHandler
(const nlohmann::json& config,
 std::shared_ptr<WaveO1TestCases>& testCase,
 std::string mesh_dir,
 const int Lx, const int L0x, const int L0t,
 bool load_init_mesh)
    : SparseGridsHandler(config, testCase)
{
    set(mesh_dir, load_init_mesh);
    set(Lx, L0x, L0t);
}

// Destructor
SparseGridsHandler
:: ~SparseGridsHandler()
{
    for (int k=0; k<m_xPointLocators.Size(); k++) {
        if (m_xPointLocators[k])
        { delete m_xPointLocators[k]; }
    }

    for (int k=0; k<m_tPointLocators.Size(); k++) {
        if (m_tPointLocators[k])
        { delete m_tPointLocators[k]; }
    }

    for (int k=0; k<m_xMeshFaceLocators.Size(); k++) {
        if (m_xMeshFaceLocators[k])
        { delete m_xMeshFaceLocators[k]; }
    }
}

void SparseGridsHandler :: set (std::string mesh_dir,
                                const bool load_init_mesh)
{
    m_load_init_mesh = load_init_mesh;
    m_mesh_dir = mesh_dir;
}

void SparseGridsHandler :: set (const int Lx)
{
    set(Lx, m_L0x, m_L0t);
}

void SparseGridsHandler :: set (const int Lx,
                                const int L0x,
                                const int L0t)
{
    m_Lx = Lx;
    m_L0x = L0x;
    m_L0t = L0t;

    // test case
    if (!m_testCase) {
        m_testCase = make_waveO1_test_case(m_config);
    }

    // number of detail- or sub-solutions
    m_numSubSols = 2*(Lx-L0x)+1;

    // set space-time levels
    // and linear combination coefficients
    // for sub-solutions
    m_levels.SetSize(m_numSubSols, 2);
    m_combCoeffs.SetSize(m_numSubSols);

    m_numUniqueLevels = Lx-L0x+1;
    for (int k=0; k<m_numUniqueLevels; k++)
    {
        m_levels(k,0) = Lx - k;
        m_levels(k,1) = L0t + k;
        m_combCoeffs[k] = +1;
    }

    for (int k=0; k<m_numUniqueLevels-1; k++)
    {
        int kk = k + m_numUniqueLevels;
        m_levels(kk,0) = Lx - 1 - k;
        m_levels(kk,1) = L0t + k;
        m_combCoeffs[kk] = -1;
    }
}

// sets minimum resolution levels
void SparseGridsHandler :: set_minLevels(const int L0x,
                                         const int L0t)
{
    m_L0x = L0x;
    m_L0t = L0t;
}

// reset memory
void SparseGridsHandler :: reset()
{
    for (int k=0; k<m_subSols.Size(); k++) {
        m_subSols[k].reset();
    }

    for (int k=0; k<m_tMeshes.Size(); k++) {
        m_tMeshes[k].reset();
    }

    for (int k=0; k<m_xMeshes.Size(); k++) {
        m_xMeshes[k].reset();
    }

    for (int k=0; k<m_discrs.Size(); k++) {
        m_discrs[k].reset();
    }

    for (int k=0; k<m_obsrvs.Size(); k++) {
        m_obsrvs[k].reset();
    }
}

std::tuple<int, double, double, Eigen::VectorXd>
SparseGridsHandler :: operator()()
{
    // compute sub-solutions
    run();

    // compute error
    auto errSol = compute_error();

    // max meshwidth in space and time
    double ht_max, hx_max=-1;
    std::tie(ht_max, hx_max) = get_meshChars();

    return {std::move(get_ndofs()),
                std::move(ht_max), std::move(hx_max),
                std::move(errSol) };
}

// Computes the sub-solutions
void SparseGridsHandler :: run()
{
    // resets memory
    reset();

    // sequence of sub-solutions
    m_subSols.SetSize(m_numSubSols);
    m_tMeshes.SetSize(m_numSubSols);
    m_xMeshes.SetSize(m_numSubSols);
    m_discrs.SetSize(m_numSubSols);
    m_obsrvs.SetSize(m_numSubSols);
#ifdef MYVERBOSE
    std::cout << "\nRun Sparse solver ..." << std::endl;
    m_levels.Print();
#endif
    WaveO1Solver *wave_solver = nullptr;
    for (int k=0; k<m_numSubSols; k++)
    {
        int lx = m_levels(k,0);
        int lt = m_levels(k,1);
        wave_solver = new WaveO1Solver(m_config,
                                     m_testCase,
                                     m_mesh_dir,
                                     lx, lt,
                                     m_load_init_mesh);
        m_tMeshes[k] = wave_solver->get_tMesh();
        m_xMeshes[k] = wave_solver->get_xMesh();
        m_discrs[k] = wave_solver->get_discretisation();
        m_obsrvs[k] = wave_solver->get_observer();
        (*wave_solver)(m_subSols[k]);
        delete wave_solver;
    }

    // set finest meshes in space and time
    m_tMesh_finest = m_tMeshes[m_numUniqueLevels-1];
    m_xMesh_finest = m_xMeshes[0];
}

// Set point locators for unique coarse xMeshes
void SparseGridsHandler :: set_xPointLocators()
{
    m_xPointLocators.SetSize(m_numUniqueLevels-1);
    for (int k=1; k<m_numUniqueLevels; k++)
    {
        //m_xMeshes[k]->PrintInfo();
        m_xPointLocators[k-1]
                = new PointLocator(m_xMeshes[k].get());
    }
}

// Set point locators for unique coarse tMeshes
void SparseGridsHandler :: set_tPointLocators()
{
    m_tPointLocators.SetSize(m_numUniqueLevels-1);
    for (int k=0; k<m_numUniqueLevels-1; k++)
    {
        //m_tMeshes[k]->PrintInfo();
        m_tPointLocators[k]
                = new PointLocator(m_tMeshes[k].get());
    }
}

// Sets mesh face locators for unique coarse xMeshes
void SparseGridsHandler :: set_xMeshFaceLocators()
{
    m_xMeshFaceLocators.SetSize(m_numUniqueLevels-1);
    for (int k=1; k<m_numUniqueLevels; k++)
    {
        m_xMeshFaceLocators[k-1]
                = new MeshFaceLocator(m_xMeshes[k].get());
    }
}


Eigen::VectorXd SparseGridsHandler :: compute_error()
{
    Eigen::VectorXd errSol(2);
    errSol.setZero();

    if (m_bool_error)
    {
        // initialize
        init_eval_error();
        set_xPointLocators();

        if (m_errorType == "L2") {
#ifdef MYVERBOSE
            std::cout << "  Computing error at end time..."
                      << std::endl;
#endif
            double errpL2, pEL2, errvL2, vEL2;
            std::tie (errpL2, pEL2, errvL2, vEL2)
                    = eval_xL2Error_at_endTime();
#ifdef MYVERBOSE
            std::cout << "\tpressure L2 error: "
                      << errpL2 << "\t" << pEL2 << std::endl;
            std::cout << "\tvelocity L2 error: "
                      << errvL2 << "\t" << vEL2 << std::endl;
#endif
            errSol(0) = errpL2;
            errSol(1) = errvL2;
        }
        else if (m_errorType == "relL2") {
#ifdef MYVERBOSE
            std::cout << "  Computing relative L2 error "
                         "at end time..."
                      << std::endl;
#endif
            double errpL2, pEL2, errvL2, vEL2;
            std::tie (errpL2, pEL2, errvL2, vEL2)
                    = eval_xL2Error_at_endTime();
#ifdef MYVERBOSE
            std::cout << "\tpressure L2 error: "
                      << errpL2 << "\t" << pEL2 << std::endl;
            std::cout << "\tvelocity L2 error: "
                      << errvL2 << "\t" << vEL2 << std::endl;
#endif
            errSol(0) = errpL2/pEL2;
            errSol(1) = errvL2/vEL2;
        }
        else if (m_errorType == "L2L2") {
#ifdef MYVERBOSE
            std::cout << "  Computing space-time L2 error..."
                      << std::endl;
#endif
            set_tPointLocators();

            double errpL2L2, pEL2L2, errvL2L2, vEL2L2;
            std::tie (errpL2L2, pEL2L2, errvL2L2, vEL2L2)
                    = eval_xtL2Error();
#ifdef MYVERBOSE
            std::cout << "\tpressure L2L2 error: "
                      << errpL2L2 << "\t"
                      << pEL2L2 << std::endl;
            std::cout << "\tvelocity L2L2 error: "
                 << errvL2L2 << "\t"
                 << vEL2L2 << std::endl;
#endif
            errSol(0) = errpL2L2;
            errSol(1) = errvL2L2;
        }
        else if (m_errorType == "relL2L2") {
#ifdef MYVERBOSE
            std::cout << "  Computing relative space-time "
                         "L2 error..."
                      << std::endl;
#endif
            set_tPointLocators();

            double errpL2L2, pEL2L2, errvL2L2, vEL2L2;
            std::tie (errpL2L2, pEL2L2, errvL2L2, vEL2L2)
                    = eval_xtL2Error();
#ifdef MYVERBOSE
            std::cout << "\tpressure L2L2 error: "
                      << errpL2L2 << "\t"
                      << pEL2L2 << std::endl;
            std::cout << "\tvelocity L2L2 error: "
                 << errvL2L2 << "\t"
                 << vEL2L2 << std::endl;
#endif
            errSol(0) = errpL2L2/pEL2L2;
            errSol(1) = errvL2L2/vEL2L2;
        }
        else if (m_errorType == "DG") {
#ifdef MYVERBOSE
            std::cout << "  Computing DG error ... "
                      << std::endl;
#endif
            set_tPointLocators();
            set_xMeshFaceLocators();

            double errDgFtime, errDgFspace;
            std::tie (errDgFtime, errDgFspace)
                    = eval_xtDgError();
#ifdef MYVERBOSE
            std::cout << "\tDG Fspace error: "
                      << errDgFspace << std::endl;
            std::cout << "\tDG Ftime error: "
                      << errDgFtime << std::endl;
#endif
            errSol(0) = errDgFspace;
            errSol(1) = errDgFtime;
        }
        else if (m_errorType == "DGPlus") {
#ifdef MYVERBOSE
            std::cout << "  Computing DG+ error ... "
                      << std::endl;
#endif
            set_tPointLocators();
            set_xMeshFaceLocators();

            double errDgpFtime, errDgpFspace;
            std::tie (errDgpFtime, errDgpFspace)
                    = eval_xtDgpError();
#ifdef MYVERBOSE
            std::cout << "\tDG+ Fspace error: "
                      << errDgpFspace << std::endl;
            std::cout << "\tDG+ Ftime error: "
                      << errDgpFtime << std::endl;
#endif
            errSol(0) = errDgpFspace;
            errSol(1) = errDgpFtime;
        }
        else {
            std::cout << "Unknown error type!" << std::endl;
            abort();
        }
    }

    return errSol;
}

// Initializes variables needed for error evaluation
void SparseGridsHandler :: init_eval_error()
{
    m_xndim = m_discrs[0]->get_xMesh()->Dimension();

    // set exact solution coefficients
    m_pE_coeff = std::make_unique<WaveO1ExactPressureCoeff>
            (m_testCase);
    m_vE_coeff = std::make_unique<WaveO1ExactVelocityCoeff>
            (m_testCase);

    // set FE spaces at the finest space and time levels
    m_tWspace_finest
            = m_discrs[m_numUniqueLevels-1]->get_tFespace();
    m_xW1space_finest = m_discrs[0]->get_xFespaces()[0];
    m_xW2space_finest = m_discrs[0]->get_xFespaces()[1];

    // exact solution on the finest space level
    m_pE = std::make_shared<GridFunction>(m_xW1space_finest);
    m_vE = std::make_shared<GridFunction>(m_xW2space_finest);
}

// Given a collection of sub-solutions
// computes the L2 error
// solution coefficients assume that time t
// has been preset in the handler functions
std::tuple <double, double, double, double>
SparseGridsHandler
:: eval_xL2Error(Array<GridFunction *> pColl,
                 Array<GridFunction *> vColl) const
{
    double errpL2=0, pEL2=0;
    double errvL2=0, vEL2=0;

    // loop over all elements on the finest space mesh,
    // search on coarser meshes, evaluate grid functions
    // and build solution.
    const FiniteElement *fe = nullptr;
    ElementTransformation *trans_finest = nullptr;
    ElementTransformation *trans_coarse = nullptr;
    Array<int> init_elIds(m_xPointLocators.Size());
    init_elIds = 0;
    for (int i=0; i<m_xMesh_finest->GetNE(); i++)
    {
        fe = m_xW1space_finest->GetFE(i);
        int order = 2*fe->GetOrder()+2;
        const IntegrationRule *ir
                = &IntRules.Get(fe->GetGeomType(), order);
        int numPoints = ir->GetNPoints();
        trans_finest = m_xW1space_finest->
                GetElementTransformation(i);

        /// use combination formula on the sub-solutions
        // evaluate sub-solution on the finest mesh
        Vector pSol;
        DenseMatrix vSol;
        pColl[0]->GetValues(i, *ir, pSol);
        vColl[0]->GetVectorValues(*trans_finest, *ir, vSol);
        pSol *= m_combCoeffs[0];
        vSol *= m_combCoeffs[0];

        // physical points to search on coarse meshes
        DenseMatrix points;
        trans_finest->Transform(*ir, points);

        // The meshes are hierarchical. So, for all points
        // in an element of the finest mesh,the point
        // locator yields the same parent element on a
        // coarse mesh, i.e. elements elIds[n] are same.
        Vector pSol_tmp(numPoints);
        DenseMatrix vSol_tmp(m_xndim, numPoints);
        Array<int> elIds(numPoints);
        Array<IntegrationPoint> ips(numPoints);
        Vector buf1, buf2;
        for (int j=0; j<m_xPointLocators.Size(); j++)
        {
            std::tie (elIds,ips)
                    = (*m_xPointLocators[j])(points,init_elIds[j]);

            int j1 = j+1;
            int j2 = j+m_numUniqueLevels;
            //std::cout << "\n";
            //std::cout << i << "\t"
            //          << j << "\t"
            //          << j1 << "\t"
            //          << j2 << "\t"
            //          << m_combCoeffs[j1] << "\t"
            //          << m_combCoeffs[j2] << std::endl;
            //elIds.Print();
            for (int k=0; k<ips.Size(); k++)
            {
                // eval pSol[j1] - pSol[j2]
                double val1 = pColl[j1]->GetValue(elIds[k],
                                                  ips[k]);
                double val2 = pColl[j2]->GetValue(elIds[k],
                                                  ips[k]);
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
            //std::cout << i << "\t" << j << "\t"
            //          << pSol_tmp.Normlinf() << std::endl;
            pSol.Add(1, pSol_tmp);
            vSol.Add(1, vSol_tmp);
        }

        /// evaluate errors
        double pE, w;
        Vector err_vE(m_xndim);
        for (int j=0; j<numPoints; j++)
        {
            const IntegrationPoint &ip = ir->IntPoint(j);
            trans_finest->SetIntPoint(&ip);
            w = ip.weight*trans_finest->Weight();
            pE = m_pE_coeff->Eval(*trans_finest, ip);
            m_vE_coeff->Eval(err_vE, *trans_finest, ip);

            // L2 error pressure
            errpL2 += w*(pE - pSol(j))*(pE - pSol(j));

            // L2 error velocity
            Vector tmp(vSol.GetColumn(j), m_xndim);
            err_vE -= tmp;
            errvL2 += w*(err_vE*err_vE);
        }
    }
    errpL2 = std::sqrt(errpL2);
    errvL2 = std::sqrt(errvL2);

    // exact solution norms
    ConstantCoefficient zero(0);
    ConstantCoefficient one(1.0);
    VectorFunctionCoefficient zeroVec(m_xndim, zeroFn);
    pEL2 = m_pE->ComputeL2Error(zero);
    vEL2 = m_vE->ComputeL2Error(zeroVec);

    return {errpL2, pEL2, errvL2, vEL2};
}

// Computes L2 and H1 error for temperature at end time
std::tuple <double, double, double, double>
SparseGridsHandler :: eval_xL2Error_at_endTime() const
{
    // set exact solution coefficients for end time
    m_pE_coeff->SetTime(m_endTime);
    m_vE_coeff->SetTime(m_endTime);
    m_pE->ProjectCoefficient(*m_pE_coeff);
    m_vE->ProjectCoefficient(*m_vE_coeff);

    // collect pressure, velocity sub-solutions at end time
    Array<GridFunction *> pColl(m_numSubSols);
    Array<GridFunction *> vColl(m_numSubSols);
    Array<int> tVdofs;
    const FiniteElement *tFe = nullptr;
    Vector tShape;
    IntegrationPoint tIp;
    tIp.Set1w(1.0, 1.0);
    for (int k=0; k<m_numSubSols; k++)
    {
        auto tWspace = m_discrs[k]->get_tFespace();
        auto xW1space = (m_discrs[k]->get_xFespaces())[0];
        auto xW2space = (m_discrs[k]->get_xFespaces())[1];
        pColl[k] = new GridFunction(xW1space);
        vColl[k] = new GridFunction(xW2space);

        tFe = tWspace->GetFE(m_tMeshes[k]->GetNE()-1);
        tWspace->GetElementVDofs(m_tMeshes[k]->GetNE()-1,
                                 tVdofs);
        tShape.SetSize(tVdofs.Size());
        tFe->CalcShape(tIp, tShape);
        build_xSol_FG(m_subSols[k]->GetBlock(0), tShape,
                      tVdofs, *pColl[k]);
        build_xSol_FG(m_subSols[k]->GetBlock(1), tShape,
                      tVdofs, *vColl[k]);
    }

    // only for visualization
    WaveO1Observer observer(m_config, -1);
    for (int k=0; k<m_numSubSols; k++) {
        observer(*pColl[k]);
    }
    //observer(*m_pE);

    auto [errpL2, pEL2, errvL2, vEL2]
            = eval_xL2Error(pColl,vColl);
    for (int k=0; k<m_numSubSols; k++) {
        delete pColl[k];
        delete vColl[k];
    }

    return {errpL2, pEL2, errvL2, vEL2};
}

// Computes L2L2 error for pressure and velocity
std::tuple <double, double, double, double>
SparseGridsHandler :: eval_xtL2Error() const
{
    double errpL2L2=0, pEL2L2=0;
    double errvL2L2=0, vEL2L2=0;

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

    Eigen::VectorXd bufErrpL2L2(Nt); bufErrpL2L2.setZero();
    Eigen::VectorXd bufErrvL2L2(Nt); bufErrvL2L2.setZero();
    Eigen::VectorXd bufpEL2L2(Nt); bufpEL2L2.setZero();
    Eigen::VectorXd bufvEL2L2(Nt); bufvEL2L2.setZero();
    double errpL2, pEL2, errvL2, vEL2;

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
    // loop over all elements of the finest time mesh
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
            //tElIdsC.Print();
            //startIds.Print();

            // set exact solution coefficients at time t
            m_pE_coeff->SetTime(t(0));
            m_vE_coeff->SetTime(t(0));
            m_pE->ProjectCoefficient(*m_pE_coeff);
            m_vE->ProjectCoefficient(*m_vE_coeff);

            // eval error at time t
            std::tie (errpL2, pEL2, errvL2, vEL2)
                    = eval_xL2Error(pColl, vColl);
            double w = tIpf.weight*tTransF->Weight();
            bufErrpL2L2(n) += w*errpL2*errpL2;
            bufErrvL2L2(n) += w*errvL2*errvL2;
            bufpEL2L2(n) += w*pEL2*pEL2;
            bufvEL2L2(n) += w*vEL2*vEL2;
        }
    }
    errpL2L2 = std::sqrt(bufErrpL2L2.sum());
    errvL2L2 = std::sqrt(bufErrvL2L2.sum());
    pEL2L2 = std::sqrt(bufpEL2L2.sum());
    vEL2L2 = std::sqrt(bufvEL2L2.sum());

    // free memory
    for (int k=0; k<m_numSubSols; k++) {
        delete pSubSols[k];
        delete vSubSols[k];
        delete pColl[k];
        delete vColl[k];
    }

    return {errpL2L2, pEL2L2, errvL2L2, vEL2L2};
}

//! Computes the sub-solution projections
std::tuple<int, double, double, Eigen::VectorXd>
SparseGridsHandler :: compute_projection()
{
    // resets memory
    reset();

    // sequence of sub-solutions
    m_subSols.SetSize(m_numSubSols);
    m_tMeshes.SetSize(m_numSubSols);
    m_xMeshes.SetSize(m_numSubSols);
    m_discrs.SetSize(m_numSubSols);
    m_obsrvs.SetSize(m_numSubSols);
#ifdef MYVERBOSE
    std::cout << "\nRun Sparse solver ..." << std::endl;
    m_levels.Print();
#endif
    WaveO1Solver *wave_solver = nullptr;
    for (int k=0; k<m_numSubSols; k++)
    {
        int lx = m_levels(k,0);
        int lt = m_levels(k,1);
        wave_solver = new WaveO1Solver
                (m_config, m_testCase, m_mesh_dir,
                 lx, lt, m_load_init_mesh);
        m_tMeshes[k] = wave_solver->get_tMesh();
        m_xMeshes[k] = wave_solver->get_xMesh();
        m_discrs[k] = wave_solver->get_discretisation();
        m_obsrvs[k] = wave_solver->get_observer();
        wave_solver->projection(m_subSols[k]);
        delete wave_solver;
    }

    // set finest meshes in space and time
    m_tMesh_finest = m_tMeshes[m_numUniqueLevels-1];
    m_xMesh_finest = m_xMeshes[0];

    // compute error
    auto errSol = compute_error();

    // max meshwidth in space and time
    double ht_max, hx_max=-1;
    std::tie(ht_max, hx_max) = get_meshChars();

    return {std::move(get_ndofs()),
                std::move(ht_max), std::move(hx_max),
                std::move(errSol) };
}

// End of file
