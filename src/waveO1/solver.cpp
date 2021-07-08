#include "../../include/waveO1/solver.hpp"
#include "../../include/waveO1/utilities.hpp"

#include <iostream>


// Constructor
WaveO1Solver :: WaveO1Solver (const nlohmann::json& config)
    : m_config(config)
{
    m_endTime = config["end_time"];

    m_bool_error = false;
    if (config.contains("eval_error")) {
        m_bool_error = config["eval_error"];
    }

    m_errorType = "relL2";
    if (config.contains("error_type")) {
        m_errorType = config["error_type"];
    }

    m_meshFormat = "mesh";
    if (config.contains("mesh_format")) {
        m_meshFormat = m_config["mesh_format"];
    }

    m_solOutputAt = "time_ids";
    if (config.contains("sol_output_at")) {
        m_solOutputAt = config["sol_output_at"];
    }

    m_bool_measureSignal = false;
    if (config.contains("measure_signal")) {
        m_bool_measureSignal = config["measure_signal"];

        if (m_bool_measureSignal) {
            std::vector<double> tmp
                    = config["measurement_point_in_x"];
            m_xMeasurementPoint.SetSize(int(tmp.size()));
            for (int i=0; i<m_xMeasurementPoint.Size(); i++)
            {
                m_xMeasurementPoint(i)
                        = tmp[static_cast<unsigned int>(i)];
            }
        }
    }

    m_bool_energy = false;
    if (config.contains("eval_energy")) {
        m_bool_energy = config["eval_energy"];
    }

    if (m_solOutputAt == "time_stamps") {
        std::vector<double> tmp
                = config["time_stamps"];
        m_timeStamps.SetSize(int(tmp.size()));
        for (int i=0; i<m_timeStamps.Size(); i++) {
            m_timeStamps(i)
                    = tmp[static_cast<unsigned int>(i)];
        }
    }
}

WaveO1Solver :: WaveO1Solver
(const nlohmann::json& config,
 std::shared_ptr<WaveO1TestCases>&
 testCase)
    : WaveO1Solver (config)
{
    m_testCase = testCase;
}

WaveO1Solver :: WaveO1Solver (const nlohmann::json& config,
                              std::string mesh_dir,
                              const int lx,
                              const int lt,
                              const bool load_init_mesh)
    : WaveO1Solver (config)
{
    set(mesh_dir, lx, lt, load_init_mesh);
}

WaveO1Solver :: WaveO1Solver (const nlohmann::json& config,
                              std::shared_ptr<WaveO1TestCases>&
                              testCase,
                              std::string mesh_dir,
                              const int lx,
                              const int lt,
                              const bool load_init_mesh)
    : WaveO1Solver (config)
{
    m_testCase = testCase;
    set(mesh_dir, lx, lt, load_init_mesh);
}

WaveO1Solver :: WaveO1Solver (const nlohmann::json& config,
                              std::string mesh_dir,
                              const int lx,
                              const bool load_init_mesh)
    : WaveO1Solver (config)
{
    set(mesh_dir, lx, -2, load_init_mesh);
}

WaveO1Solver :: WaveO1Solver (const nlohmann::json& config,
                              std::shared_ptr<WaveO1TestCases>&
                              testCase,
                              std::string mesh_dir,
                              const int lx,
                              const bool load_init_mesh)
    : WaveO1Solver (config)
{
    m_testCase = testCase;
    set(mesh_dir, lx, -2, load_init_mesh);
}

// Destructor
WaveO1Solver :: ~ WaveO1Solver ()
{
    if (!m_pardiso_finalized) {
        finalize();
    }
}

// Sets test case, meshes, discretization, observer
void WaveO1Solver :: set (std::string mesh_dir,
                          const int lx,
                          const int lt,
                          const bool load_init_mesh)
{
    set(mesh_dir, load_init_mesh);
#ifdef MYVERBOSE
    std::cout << "\nRun Wave solver ..." << std::endl;
#endif
    set(lx, lt);
}

void WaveO1Solver :: set (std::string mesh_dir,
                          const bool load_init_mesh)
{
    m_load_init_mesh = load_init_mesh;
    m_mesh_dir = mesh_dir;

    m_lx0 = 0;
    if (m_load_init_mesh) {
        if (m_config.contains("init_mesh_level")) {
            m_lx0 = m_config["init_mesh_level"];
        }
    }

    m_mesh_elem_type = "tri";
    if (m_config.contains("mesh_elem_type")) {
        m_mesh_elem_type = m_config["mesh_elem_type"];
    }
}

void WaveO1Solver :: set(int lx, int lt)
{
    // test case
    if (!m_testCase) {
        m_testCase = make_waveO1_test_case(m_config);
    }

    // meshes
    set_meshes(lx, lt);

    // space-time discretisation
    m_discr.reset();
    m_discr = std::make_shared<WaveO1XtDG>
            (m_config, m_testCase);
    m_discr->set(m_tMesh, m_xMesh);
    m_blockOffsets = m_discr->get_block_offsets();
    m_slab_blockOffsets = m_discr->get_slab_block_offsets();

    // observer
    m_observer.reset();
    m_observer = std::make_shared<WaveO1Observer>
            (m_config, lx);
    m_observer->set(m_testCase, m_discr);
}

// Sets meshes
void WaveO1Solver :: set_meshes(int lx, int lt)
{
    // mesh in space
    m_xMesh.reset();
    if (m_load_init_mesh) { // load initial mesh and refine
        int num_refinements = lx - m_lx0;
        const std::string mesh_file
                = m_mesh_dir+"/"+m_mesh_elem_type+"_mesh_l"
                +std::to_string(m_lx0)+"."+m_meshFormat;
//#ifdef MYVERBOSE
        std::cout << "  Initial mesh file: "
                  << mesh_file << std::endl;
        std::cout << "  Number of uniform refinements: "
             << num_refinements << std::endl;
//#endif
        m_xMesh = std::make_shared<Mesh>(mesh_file.c_str());
        for (int k=0; k<num_refinements; k++) {
            m_xMesh->UniformRefinement();
        }
    }
    else {
        const std::string mesh_file =
                m_mesh_dir+"/mesh_l"
                +std::to_string(lx)+"."+m_meshFormat;
#ifdef MYVERBOSE
        std::cout << "  Mesh file: "
                  << mesh_file << std::endl;
#endif
        m_xMesh = std::make_shared<Mesh>(mesh_file.c_str());
    }

    // mesh in time
    m_tMesh.reset();
    if (lt == -2) {
        double hx_min, hx_max, kappax_min, kappax_max;
        m_xMesh->GetCharacteristics(hx_min, hx_max,
                                    kappax_min, kappax_max);

        // set time-level acc. to spatial-level
        lt = int(std::ceil(-std::log2(hx_max/m_endTime)));
    }
    int Nt = static_cast<int>(std::pow(2,lt));
    m_tMesh = std::make_shared<Mesh>(Nt, m_endTime);
}

// Runs the solver; used for the full-tensor version
std::tuple<int, double, double, Eigen::VectorXd> WaveO1Solver
:: operator() ()
{
    // run solver
    run();

    // visualization at end time
    auto pressure = get_pressureSol_at_endTime(m_W.get());
    (*m_observer)(pressure);

    // measure signal
    if (m_bool_measureSignal) {
        measure_signal(m_W.get());
    }

    // dump mesh
    m_observer->dump_mesh(m_xMesh);

    if (m_solOutputAt == "time_ids") {
        // dump solution at time steps
        dump_sol_at_time_steps(m_W.get());
    }
    else if (m_solOutputAt == "time_stamps") {
        // dump solution at time stamps
        dump_sol_at_time_stamps(m_W.get());
    }

    // compute error
    auto errSol = compute_error(m_W.get());

    // compute energy
    compute_energy(m_W.get());

    // max meshwidth in space and time
    double ht_max, hx_max=-1;
    std::tie(ht_max, hx_max) = get_meshChars();

    return {std::move(get_ndofs()),
                std::move(ht_max), std::move(hx_max),
                std::move(errSol) };
}

// Runs the solver; used for the sparse-grids handler
void WaveO1Solver
:: operator()(std::unique_ptr<BlockVector>& W)
{
    // MFEM solution variables
    W = std::make_unique<BlockVector>(m_blockOffsets);

    // MFEM solution and rhs variables for a time slab
    m_slabW = std::make_unique<BlockVector>
            (m_slab_blockOffsets);
    m_slabB = std::make_unique<BlockVector>
            (m_slab_blockOffsets);

    // initialize
    init ();

    // solve
    solve (W.get());

    // visualization at end time
    auto pressure = get_pressureSol_at_endTime(W.get());
    (*m_observer)(pressure);
}

// Runs the solver
void WaveO1Solver
:: run()
{
    // MFEM solution variables
    m_W.reset();
    m_W = std::make_unique<BlockVector>(m_blockOffsets);

    // MFEM solution and rhs variables for a time slab
    m_slabW.reset();
    m_slabB.reset();
    m_slabW = std::make_unique<BlockVector>
            (m_slab_blockOffsets);
    m_slabB = std::make_unique<BlockVector>
            (m_slab_blockOffsets);

    // initialize
    init ();

    // solve
    solve (m_W.get());
}

// Initializes the solver
void WaveO1Solver :: init ()
{
    assemble_system();
}

// Assembles the system
void WaveO1Solver :: assemble_system()
{
    if (!m_pardiso_finalized) {
        finalize();
    }

    m_discr->assemble_system();

    m_discr->build_system_matrix();
    m_waveMat = m_discr->get_wave_mat();

    // set Pardiso solver
    int mtype = 11; // real indefinite
    m_pardisoSolver.reset();
    m_pardisoSolver = std::make_unique<PardisoSolver>(mtype);
    m_pardisoSolver->initialize(m_waveMat->Size(),
                                m_waveMat->GetI(),
                                m_waveMat->GetJ(),
                                m_waveMat->GetData());
    m_pardisoSolver->factorize();
    m_pardiso_finalized = false;
}

void WaveO1Solver :: solve (BlockVector *W)
{
    int Nt = m_tMesh->GetNE();
    (*W) = 0.;
    for (int n=0; n<Nt; n++) {
        solve_one_time_slab(n, W);
    }
}

void WaveO1Solver :: solve_one_time_slab (const int n,
                                        BlockVector *W)
{
    if (n == 0) { assemble_rhs_initial(); }
    else { assemble_rhs(n); }
    m_pardisoSolver->solve(m_slabB->GetData(),
                           m_slabW->GetData());
    copy_slabW_to_W(n, m_slabW.get(), W);
}

// Assembles the right-hand side in the first time slab
void WaveO1Solver :: assemble_rhs_initial()
{
    (*m_slabB) = 0.0;
    m_discr->assemble_rhs_initial(m_slabB.get());
}

void WaveO1Solver :: assemble_rhs(int n)
{
    (*m_slabB) = 0.0;
    m_discr->assemble_rhs(n, m_slabW.get(), m_slabB.get());
}

void WaveO1Solver :: copy_slabW_to_W(const int n,
                                   const BlockVector *slabW,
                                   BlockVector *W) const
{
    m_discr->set_vector(n, slabW->GetBlock(0),
                        W->GetBlock(0));
    m_discr->set_vector(n, slabW->GetBlock(1),
                        W->GetBlock(1));
}

// Releases memory
void WaveO1Solver :: finalize() const
{
    if (!m_pardiso_finalized) {
#ifdef MYVERBOSE
        std::cout << "\nRelease Pardiso Memory"
                  << std::endl;
#endif
        m_pardisoSolver->finalize();
        m_pardiso_finalized = true;
    }
}

// Computes the numerical solution error
Eigen::VectorXd WaveO1Solver
:: compute_error(BlockVector *W) const
{
    Eigen::VectorXd errSol(2);
    errSol.setZero();

    if (m_bool_error)
    {
        if (m_errorType == "L2") {
#ifdef MYVERBOSE
            std::cout << "  Computing L2 error "
                         "at end time..."
                      << std::endl;
#endif
            auto pressure = get_pressureSol_at_endTime(W);
            auto velocity = get_velocitySol_at_endTime(W);

            double errpL2, pEL2, errvL2, vEL2;
            std::tie (errpL2, pEL2, errvL2, vEL2)
                    = m_observer->eval_xL2Error
                    (pressure, velocity, m_endTime);
#ifdef MYVERBOSE
            std::cout << "\tpressure L2 error: "
                      << errpL2 << "\t"
                      << pEL2 << std::endl;
            std::cout << "\tvelocity L2 error: "
                      << errvL2 << "\t"
                      << vEL2 << std::endl;
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
            auto pressure = get_pressureSol_at_endTime(W);
            auto velocity = get_velocitySol_at_endTime(W);

            double errpL2, pEL2, errvL2, vEL2;
            std::tie (errpL2, pEL2, errvL2, vEL2)
                    = m_observer->eval_xL2Error
                    (pressure, velocity, m_endTime);
#ifdef MYVERBOSE
            std::cout << "\tpressure L2 error: "
                      << errpL2 << "\t"
                      << pEL2 << std::endl;
            std::cout << "\tvelocity L2 error: "
                      << errvL2 << "\t"
                      << vEL2 << std::endl;
#endif
            errSol(0) = errpL2/pEL2;
            errSol(1) = errvL2/vEL2;
        }
        else if (m_errorType == "L2L2") {
#ifdef MYVERBOSE
            std::cout << "  Computing space-time L2 error..."
                      << std::endl;
#endif
            double errpL2L2, pEL2L2, errvL2L2, vEL2L2;
            std::tie (errpL2L2, pEL2L2, errvL2L2, vEL2L2)
                    = m_observer->eval_xtL2Error(*W);
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
            double errpL2L2, pEL2L2, errvL2L2, vEL2L2;
            std::tie (errpL2L2, pEL2L2, errvL2L2, vEL2L2)
                    = m_observer->eval_xtL2Error(*W);
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
            std::cout << "  Computing space-time DG error..."
                      << std::endl;
#endif
            double errDgFtime, errDgFspace;
            std::tie (errDgFtime, errDgFspace)
                    = m_observer->eval_xtDgError(*W);
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
//#ifdef MYVERBOSE
            std::cout << "  Computing space-time DG+ error..."
                      << std::endl;
//#endif
            double errDgpFtime, errDgpFspace;
            std::tie (errDgpFtime, errDgpFspace)
                    = m_observer->eval_xtDgpError(*W);
//#ifdef MYVERBOSE
            std::cout << "\tDG+ Fspace error: "
                      << errDgpFspace << std::endl;
            std::cout << "\tDG+ Ftime error: "
                      << errDgpFtime << std::endl;
//#endif
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

//! Computes the energy for all t=t_n
void WaveO1Solver
:: compute_energy(BlockVector *W) const
{
    if (m_bool_energy)
    {
        // eval energy
        auto energy = m_observer->eval_energy(*W);

        // time-series
        int Nt = m_tMesh->GetNE();
        double dt = m_endTime/Nt;
        Vector timeSeries(Nt+1);
        for (int n=0; n<=Nt; n++) {
            timeSeries(n) = n*dt;
        }

        m_observer->dump_energy(timeSeries, energy);
    }
}

// Evaluates the pressure solution vector at end time T
std::shared_ptr <GridFunction> WaveO1Solver
:: get_pressureSol_at_endTime(BlockVector *W) const
{
    // auxiliary variables
    auto tFespace = m_discr->get_tFespace();
    auto xFespaces = m_discr->get_xFespaces();

    // dofs for the time element (T-h, T)
    Array<int> tVdofs;
    const FiniteElement *tFe
            = tFespace->GetFE(m_tMesh->GetNE()-1);
    tFespace->GetElementVDofs(m_tMesh->GetNE()-1, tVdofs);

    Vector tShape(tVdofs.Size());
    IntegrationPoint ip;
    ip.Set1w(1.0, 1.0);
    tFe->CalcShape(ip, tShape);

    // pressure solution at T
    Vector pSol;
    std::shared_ptr <GridFunction> pressure
            = std::make_shared<GridFunction>(xFespaces[0]);
    pressure->GetTrueDofs(pSol);
    build_xSol_FG(W->GetBlock(0), tShape, tVdofs, pSol);

    return pressure;
}

// Evaluates the velocity solution vector at end time T
std::shared_ptr <GridFunction> WaveO1Solver
:: get_velocitySol_at_endTime(BlockVector *W) const
{
    // auxiliary variables for visualization
    // and error computation at time T
    auto tFespace = m_discr->get_tFespace();
    auto xFespaces = m_discr->get_xFespaces();

    // dofs for the time element (T-h, T)
    Array<int> tVdofs;
    const FiniteElement *tFe
            = tFespace->GetFE(m_tMesh->GetNE()-1);
    tFespace->GetElementVDofs(m_tMesh->GetNE()-1, tVdofs);

    Vector tShape(tVdofs.Size());
    IntegrationPoint ip;
    ip.Set1w(1.0, 1.0);
    tFe->CalcShape(ip, tShape);

    // flux solution at T
    Vector vSol;
    std::shared_ptr <GridFunction> velocity
            = std::make_shared<GridFunction>(xFespaces[1]);
    velocity->GetTrueDofs(vSol);
    build_xSol_FG(W->GetBlock(1), tShape, tVdofs, vSol);

    return velocity;
}

// Sets the pressure solution vector at end time T
void WaveO1Solver
:: set_pressureSol_at_endTime()
{
    m_pressure = get_pressureSol_at_endTime(m_W.get());
}

// Sets the flux solution vector at end time T
void WaveO1Solver
:: set_velocitySol_at_endTime()
{
    m_velocity = get_velocitySol_at_endTime(m_W.get());
}

//! Dumps the solution at specified time steps
void WaveO1Solver
:: dump_sol_at_time_steps(BlockVector *W) const
{
    auto tFespace = m_discr->get_tFespace();
    auto xFespaces = m_discr->get_xFespaces();
    int Nt = m_tMesh->GetNE();

    // pressure solution variable
    Vector pSol;
    std::shared_ptr <GridFunction> pressure
            = std::make_shared<GridFunction>(xFespaces[0]);
    pressure->GetTrueDofs(pSol);

    // velocity solution variable
    Vector vSol;
    std::shared_ptr <GridFunction> velocity
            = std::make_shared<GridFunction>(xFespaces[1]);
    velocity->GetTrueDofs(vSol);

    Vector tShape;
    Array<int> tVdofs;
    IntegrationPoint ip;

    // write solution at t=0
    ip.Set1w(0.0, 1.0);
    const FiniteElement *tFe = tFespace->GetFE(0);
    tFespace->GetElementVDofs(0, tVdofs);
    tShape.SetSize(tVdofs.Size());
    tFe->CalcShape(ip, tShape);
    build_xSol_FG(W->GetBlock(0), tShape,
                  tVdofs, pSol);
    build_xSol_FG(W->GetBlock(1), tShape,
                  tVdofs, vSol);
    m_observer->dump_sol_at_time_steps
            (pressure, velocity, 0);

    // write solution at t=t_{n}
    ip.Set1w(1.0, 1.0);
    for (int n=0; n<Nt; n++)
    {
        // dofs for the time element n
        const FiniteElement *tFe = tFespace->GetFE(n);
        tFespace->GetElementVDofs(n, tVdofs);

        // build solution
        tShape.SetSize(tVdofs.Size());
        tFe->CalcShape(ip, tShape);
        build_xSol_FG(W->GetBlock(0), tShape,
                      tVdofs, pSol);
        build_xSol_FG(W->GetBlock(1), tShape,
                      tVdofs, vSol);

        // write
        m_observer->dump_sol_at_time_steps
                (pressure, velocity, n+1);
    }
}

//! Dumps the solution at specified time stamps
void WaveO1Solver
:: dump_sol_at_time_stamps(BlockVector *W) const
{
    auto tFespace = m_discr->get_tFespace();
    auto xFespaces = m_discr->get_xFespaces();

    // locate time stamps on the time mesh
    int M = m_timeStamps.Size();
    Array<int> tElIds(M);
    Array<IntegrationPoint> tIps(M);
    PointLocator point_locator(m_tMesh.get());
    for (int i=0; i<M; i++) {
        Vector t(1);
        t(0) = m_timeStamps(i);
        std::tie (tElIds[i], tIps[i])
                = point_locator(t, 0);
        std::cout << "Measurement point t = "
                  << t(0) << " "
                  << "is located in mesh element "
                  << tElIds[i] << std::endl;
    }

    // pressure solution variable
    Vector pSol;
    std::shared_ptr <GridFunction> pressure
            = std::make_shared<GridFunction>(xFespaces[0]);
    pressure->GetTrueDofs(pSol);

    // velocity solution variable
    Vector vSol;
    std::shared_ptr <GridFunction> velocity
            = std::make_shared<GridFunction>(xFespaces[1]);
    velocity->GetTrueDofs(vSol);

    Vector tShape;
    Array<int> tVdofs;
    for (int i=0; i<M; i++)
    {
        // dofs for the time element n
        const FiniteElement *tFe
                = tFespace->GetFE(tElIds[i]);
        tFespace->GetElementVDofs(tElIds[i], tVdofs);

        tShape.SetSize(tVdofs.Size());
        tFe->CalcShape(tIps[i], tShape);
        build_xSol_FG(W->GetBlock(0), tShape,
                      tVdofs, pSol);
        build_xSol_FG(W->GetBlock(1), tShape,
                      tVdofs, vSol);

        m_observer->dump_sol_at_time_stamps
                (pressure, velocity, m_timeStamps(i));
    }
}

//! Evaluates the velocity solution vector at end time T
void WaveO1Solver
:: measure_signal(BlockVector *W) const
{
    auto tFespace = m_discr->get_tFespace();
    auto xFespaces = m_discr->get_xFespaces();
    int Nt = m_tMesh->GetNE();
    int M = 10; // # measurement points per time element

    // locate measurement point in physical space
    int xElId;
    IntegrationPoint xIp;
    PointLocator pointLocator(m_xMesh.get());
    std::tie (xElId,xIp)
            = pointLocator(m_xMeasurementPoint, 0);
    std::cout << "Measurement point x = ("
              << m_xMeasurementPoint(0) << ", "
              << m_xMeasurementPoint(1) << ") "
              << "is located in mesh element "
              << xElId << std::endl;

    // set measurement points in time
    Array<IntegrationPoint> tIps(M);
    for (int i=0; i<M; i++) {
        tIps[i].Set1w((i+1)*(1./M), 1.0);
    }

    // signal and solution variables
    Vector timeSeries(M*Nt);
    Vector pressureSignal(M*Nt);
    Vector pSol;
    std::shared_ptr <GridFunction> pressure
            = std::make_shared<GridFunction>(xFespaces[0]);
    pressure->GetTrueDofs(pSol);

    // measure
    double dt = m_endTime/Nt;
    Array<int> tVdofs;
    Vector tShape;
    for (int n=0; n<Nt; n++)
    {
        // dofs for the time element n
        const FiniteElement *tFe = tFespace->GetFE(n);
        tFespace->GetElementVDofs(n, tVdofs);

        // evaluate pressure at measurement points in time
        tShape.SetSize(tVdofs.Size());
        for (int i=0; i<M; i++)
        {
            tFe->CalcShape(tIps[i], tShape);
            build_xSol_FG(W->GetBlock(0), tShape,
                          tVdofs, pSol);
            pressureSignal(i+n*M)
                    = pressure->GetValue(xElId, xIp);
            timeSeries(i+n*M) = n*dt + dt*tIps[i].x;
        }
    }
    m_observer->dump_pressure_signal
            (timeSeries, pressureSignal);
}

// Routines for computing L2-projection
std::tuple<int, double, double, Eigen::VectorXd> WaveO1Solver
:: projection ()
{
    // MFEM solution variable
    m_W = std::make_unique<BlockVector>(m_blockOffsets);

    // MFEM solution and rhs variables for a time slab
    m_slabW = std::make_unique<BlockVector>
            (m_slab_blockOffsets);
    m_slabB = std::make_unique<BlockVector>
            (m_slab_blockOffsets);

    // initialize_projection
    init_projection();

    // solve
    solve_projection (m_W.get());

    // visualization at end time
    auto pressure = get_pressureSol_at_endTime(m_W.get());
    (*m_observer)(pressure);

    auto velocity = get_velocitySol_at_endTime(m_W.get());
    (*m_observer)(velocity);

    // compute error
    auto errSol = compute_error(m_W.get());

    // max meshwidth in space and time
    double ht_max, hx_max;
    std::tie(ht_max, hx_max) = get_meshChars();

    return {std::move(get_ndofs()),
                std::move(ht_max), std::move(hx_max),
                std::move(errSol) };
}

void WaveO1Solver
:: projection (std::unique_ptr<BlockVector>& W)
{
    // MFEM solution variable
    W = std::make_unique<BlockVector>(m_blockOffsets);

    // MFEM solution and rhs variables for a time slab
    m_slabW = std::make_unique<BlockVector>
            (m_slab_blockOffsets);
    m_slabB = std::make_unique<BlockVector>
            (m_slab_blockOffsets);

    // initialize projection
    init_projection();

    // compute projection
    solve_projection(W.get());
}

void WaveO1Solver :: init_projection()
{
    assemble_projection_system();
}

void WaveO1Solver :: assemble_projection_system()
{
    if (!m_pardiso_finalized) {
        finalize();
    }

    m_discr->assemble_projector();
    m_discr->build_projector_matrix();
    m_waveProjMat = m_discr->get_wave_projection_mat();

    // set Pardiso solver
    int mtype = 2; // real symmetric positive definite
    m_pardisoSolver = std::make_unique<PardisoSolver>(mtype);
    m_pardisoSolver->initialize(m_waveProjMat->Size(),
                                m_waveProjMat->GetI(),
                                m_waveProjMat->GetJ(),
                                m_waveProjMat->GetData());
    m_pardisoSolver->factorize();
    m_pardiso_finalized = false;
}

void WaveO1Solver :: solve_projection (BlockVector *W)
{
    int Nt = m_tMesh->GetNE();

    (*W) = 0.0;
    for (int n=0; n<Nt; n++) {
        assemble_projection_rhs(n);
        m_pardisoSolver->solve(m_slabB->GetData(),
                               m_slabW->GetData());
        copy_slabW_to_W(n, m_slabW.get(), W);
    }
}

void WaveO1Solver :: assemble_projection_rhs(int n) const
{
    (*m_slabB) = 0.0;
    m_discr->assemble_projection_rhs(n, m_slabB.get());
}

// End of file
