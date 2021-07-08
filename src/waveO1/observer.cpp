#include "../../include/waveO1/observer.hpp"
#include "../../include/waveO1/utilities.hpp"
#include "../../include/mymfem/utilities.hpp"

#include <fstream>
#include <iostream>
#include <Eigen/Core>

#include <filesystem>
namespace fs = std::filesystem;


// Constructor
WaveO1Observer
:: WaveO1Observer (const nlohmann::json& config, int lx)
    : BaseObserver (config, lx)
{
    std::string problem_type = config["problem_type"];

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

    m_bool_energy = false;
    if (config.contains("eval_energy")) {
        m_bool_energy = config["eval_energy"];
    }
    
    std::string base_out_dir = "../output";
    if (config.contains("base_out_dir")) {
        base_out_dir = config["base_out_dir"];
    }

    std::string sub_out_dir = "wave_"+problem_type;
    if (config.contains("sub_out_dir")) {
        sub_out_dir = config["sub_out_dir"];
    }

    m_output_dir = base_out_dir+"/"+sub_out_dir+"/lx"
            +std::to_string(lx)+"/";
    m_solName_suffix = "_lx"+std::to_string(lx);

    if (m_bool_dumpOut) {
        m_dump_out_nsteps = 8;
        if (config.contains("dump_output_num_time_steps")) {
            m_dump_out_nsteps
                    = config["dump_output_num_time_steps"];
        }
    }

    m_bool_dumpSignal = false;
    if (config.contains("dump_signal")) {
        m_bool_dumpSignal = config["dump_signal"];
    }

    m_bool_dumpEnergy = false;
    if (config.contains("dump_energy")) {
        m_bool_dumpEnergy = config["dump_energy"];
    }

    if (m_bool_dumpOut
            || m_bool_dumpSignal
            || m_bool_dumpEnergy) {
        fs::create_directories(m_output_dir);
    }
}

void WaveO1Observer
:: set (std::shared_ptr<WaveO1TestCases>& testCase,
        std::shared_ptr<WaveO1XtDG>& discr)
{
    m_testCase = testCase;
    m_discr = discr;

    m_xndim = discr->get_xMesh()->Dimension();
    m_tWspace = m_discr->get_tFespace();
    m_xW1space = m_discr->get_xFespaces()[0];
    m_xW2space = m_discr->get_xFespaces()[1];

    if (m_bool_error || m_bool_energy)
    {
        m_pE_coeff
                = std::make_unique<WaveO1ExactPressureCoeff>
                (m_testCase);
        m_vE_coeff
                = std::make_unique<WaveO1ExactVelocityCoeff>
                (m_testCase);
    }

    if (m_bool_error)
    {
        m_pE = std::make_unique<GridFunction>(m_xW1space);
        m_vE = std::make_unique<GridFunction>(m_xW2space);

        if (m_errorType == "DG")
        {
            Array<int> xEss_bdr_marker
                    = m_discr->get_xEss_bdr_marker();
            Array<int> xNat_bdr_marker
                    = m_discr->get_xNat_bdr_marker();

            m_assembleDgError
                    = std::make_unique
                    <mymfem::AssembleDgError>
                    (xEss_bdr_marker,
                     xNat_bdr_marker,
                     m_stabParamsType);
        }
        else if (m_errorType == "DGPlus")
        {
            Array<int> xEss_bdr_marker
                    = m_discr->get_xEss_bdr_marker();
            Array<int> xNat_bdr_marker
                    = m_discr->get_xNat_bdr_marker();

            m_assembleDgpSmoothError
                    = std::make_unique
                    <mymfem::AssembleDgPlusSmoothError>
                    (xEss_bdr_marker,
                     xNat_bdr_marker,
                     m_stabParamsType);
        }
    }

    if (m_bool_energy)
    {
        Array<int> xEss_bdr_marker
                = m_discr->get_xEss_bdr_marker();
        Array<int> xNat_bdr_marker
                = m_discr->get_xNat_bdr_marker();

        m_assembleDgJumps
                = std::make_unique
                <mymfem::AssembleDgJumps>
                (xEss_bdr_marker,
                 xNat_bdr_marker,
                 m_stabParamsType);
    }
}

//! Dumps solution
void WaveO1Observer
:: dump_sol (std::shared_ptr<GridFunction>& p,
             std::shared_ptr<GridFunction>& v) const
{
    if (m_bool_dumpOut)
    {
        std::string sol_name1
                = m_output_dir+"pressure"+m_solName_suffix;
        std::string sol_name2
                = m_output_dir+"velocity"+m_solName_suffix;
        std::cout << sol_name1 << std::endl;
        std::cout << sol_name2 << std::endl;

        std::ofstream sol_ofs1(sol_name1.c_str());
        sol_ofs1.precision(m_precision);
        p->Save(sol_ofs1);

        std::ofstream sol_ofs2(sol_name2.c_str());
        sol_ofs2.precision(m_precision);
        v->Save(sol_ofs2);
    }
}

//! Dumps solution at given time ids
void WaveO1Observer
:: dump_sol_at_time_steps
(std::shared_ptr<GridFunction>& p,
 std::shared_ptr<GridFunction>& v,
 int time_step_id) const
{
    if (!m_bool_dumpOut) { return; }

    if (time_step_id % m_dump_out_nsteps == 0)
    {
        std::string solName_suffix = m_solName_suffix
                +"_tId"+std::to_string(time_step_id);

        std::string sol_name1
                = m_output_dir+"pressure"+solName_suffix;
        std::string sol_name2
                = m_output_dir+"velocity"+solName_suffix;
        std::cout << sol_name1 << std::endl;
        std::cout << sol_name2 << std::endl;

        std::ofstream sol_ofs1(sol_name1.c_str());
        sol_ofs1.precision(m_precision);
        p->Save(sol_ofs1);

        std::ofstream sol_ofs2(sol_name2.c_str());
        sol_ofs2.precision(m_precision);
        v->Save(sol_ofs2);
    }
}

//! Dumps solution at given time stamps
void WaveO1Observer
:: dump_sol_at_time_stamps
(std::shared_ptr<GridFunction>& p,
 std::shared_ptr<GridFunction>& v,
 double t) const
{
    if (!m_bool_dumpOut) { return; }

    {
        std::string solName_suffix = m_solName_suffix
                +"_t"+std::to_string(t);

        std::string sol_name1
                = m_output_dir+"pressure"+solName_suffix;
        std::string sol_name2
                = m_output_dir+"velocity"+solName_suffix;
        std::cout << sol_name1 << std::endl;
        std::cout << sol_name2 << std::endl;

        std::ofstream sol_ofs1(sol_name1.c_str());
        sol_ofs1.precision(m_precision);
        p->Save(sol_ofs1);

        std::ofstream sol_ofs2(sol_name2.c_str());
        sol_ofs2.precision(m_precision);
        v->Save(sol_ofs2);
    }
}

// Dumps pressure signal
void WaveO1Observer
:: dump_pressure_signal (Vector& time_series,
                         Vector& signal) const
{
    if (m_bool_dumpSignal)
    {
        std::string outfile
                = m_output_dir+"pressure_signal"
                +m_solName_suffix+".json";
        std::cout << outfile << std::endl;

        auto json = nlohmann::json{};
        for (int i=0; i<time_series.Size(); i++) {
            json["time"][static_cast<unsigned int>(i)]
                    = time_series(i);
            json["pressure"][static_cast<unsigned int>(i)]
                    = signal(i);
        }

        auto file = std::ofstream(outfile);
        assert(file.good());
        file << json.dump(2);
        file.close();
    }
}

//! Dumps energy
void WaveO1Observer
:: dump_energy (Vector& timeSeries,
                Vector& energy) const
{
    if (m_bool_dumpEnergy)
    {
        std::string outfile
                = m_output_dir+"energy_evolution"
                +m_solName_suffix+".json";
        std::cout << outfile << std::endl;

        auto json = nlohmann::json{};
        for (int i=0; i<energy.Size(); i++) {
            json["time"][static_cast<unsigned int>(i)]
                    = timeSeries(i);
            json["energy"][static_cast<unsigned int>(i)]
                    = energy(i);
        }

        auto file = std::ofstream(outfile);
        assert(file.good());
        file << json.dump(2);
        file.close();
    }
}

// Compute spatial L2 error for pressure and velocity
std::tuple <double, double, double, double> WaveO1Observer
:: eval_xL2Error (std::shared_ptr<GridFunction>& p,
                  std::shared_ptr<GridFunction>& v,
                  double t) const
{
    double errpL2=0, pEL2=0;
    double errvL2=0, vEL2=0;

    ConstantCoefficient zero(0);
    VectorFunctionCoefficient zeroVec(m_xndim, zeroFn);

    // assumes that all elements have the same geometry
    auto elGeomType = p->FESpace()
                ->GetFE(0)->GetGeomType();

    // L2 error in pressure
    // set integration rules
    IntegrationRules rule1{};
    const IntegrationRule
            *irs1[Geometry::Type::NUM_GEOMETRIES];
    int order1 = 2*p->FESpace()->GetFE(0)->GetOrder()+2;
    irs1[elGeomType] = &rule1.Get(elGeomType, order1);
    m_pE_coeff->SetTime(t);
    m_pE->ProjectCoefficient(*m_pE_coeff);
    pEL2 = m_pE->ComputeL2Error(zero, irs1);
    errpL2 = p->ComputeL2Error(*m_pE_coeff, irs1);

    // L2 error in velocity
    // set integration rules
    IntegrationRules rule2{};
    const IntegrationRule
            *irs2[Geometry::Type::NUM_GEOMETRIES];
    int order2 = 2*v->FESpace()->GetFE(0)->GetOrder()+2;
    irs2[elGeomType] = &rule2.Get(elGeomType, order2);
    m_vE_coeff->SetTime(t);
    m_vE->ProjectCoefficient(*m_vE_coeff);
    vEL2 = m_vE->ComputeL2Error(zeroVec, irs2);
    errvL2 = v->ComputeL2Error(*m_vE_coeff, irs2);

    return {errpL2, pEL2, errvL2, vEL2};
}

// Compute space-time L2 error for pressure and velocity
std::tuple <double, double, double, double> WaveO1Observer
:: eval_xtL2Error (BlockVector& W) const
{
    double errpL2L2=0, pEL2L2=0;
    double errvL2L2=0, vEL2L2=0;

    int Nt = m_tWspace->GetNE();
    int xdimW1 = m_xW1space->GetTrueVSize();
    int xdimW2 = m_xW2space->GetTrueVSize();
    
    Vector& p = W.GetBlock(0);
    Vector& v = W.GetBlock(1);

    Vector pSol(xdimW1);
    std::shared_ptr<GridFunction> pGSol
            = std::make_shared<GridFunction>(m_xW1space,
                                             pSol);
            
    Vector vSol(xdimW2);
    std::shared_ptr<GridFunction> vGSol
            = std::make_shared<GridFunction>(m_xW2space,
                                             vSol);

    ElementTransformation *tTrans = nullptr;
    const FiniteElement *tFe = nullptr;
    Vector tShape;

    Eigen::VectorXd bufErrpL2L2(Nt); bufErrpL2L2.setZero();
    Eigen::VectorXd bufErrvL2L2(Nt); bufErrvL2L2.setZero();
    Eigen::VectorXd bufpEL2L2(Nt); bufpEL2L2.setZero();
    Eigen::VectorXd bufvEL2L2(Nt); bufvEL2L2.setZero();

    Array<int> tVdofs;
    double errpL2, pEL2, errvL2, vEL2;
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

            std::tie (errpL2, pEL2, errvL2, vEL2)
                    = eval_xL2Error(pGSol, vGSol, t(0));

            double w = ip.weight*tTrans->Weight();
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

    bufErrpL2L2.setZero();
    bufErrvL2L2.setZero();
    bufpEL2L2.setZero();
    bufvEL2L2.setZero();

    return {errpL2L2, pEL2L2, errvL2L2, vEL2L2};
}

// End of file
