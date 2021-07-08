#include "../../include/waveO2/observer.hpp"
#include "../../include/waveO2/utilities.hpp"
#include "../../include/mymfem/utilities.hpp"

#include <fstream>
#include <iostream>
#include <Eigen/Core>

#include <filesystem>
namespace fs = std::filesystem;


//! Constructor
WaveO2Observer
:: WaveO2Observer (const nlohmann::json& config, int lx)
    : BaseObserver (config, lx)
{
    std::string problem_type = config["problem_type"];

    m_bool_error = false;
    if (config.contains("eval_error")) {
        m_bool_error = config["eval_error"];
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
        fs::create_directories(m_output_dir);
    }
}

//! Sets the exact solution coefficients
//! and auxiliary variables needed
void WaveO2Observer
:: set (std::shared_ptr<WaveO2TestCases>& testCase,
        std::shared_ptr<WaveO2XtFEM>& discr)
{
    m_testCase = testCase;
    m_discr = discr;

    m_xndim = discr->get_xMesh()->Dimension();
    m_tVspace = m_discr->get_tFespace();
    m_xVspace = m_discr->get_xFespaces();

    if (m_bool_error)
    {
        m_uECoeff = std::make_unique
                <WaveO2ExactSolutionCoeff>(m_testCase);
        m_dudtECoeff = std::make_unique
                <WaveO2ExactGradtSolutionCoeff>(m_testCase);
        m_dudxECoeff = std::make_unique
                <WaveO2ExactGradxSolutionCoeff>(m_testCase);

        m_uE = std::make_unique<GridFunction>(m_xVspace);
        m_dudtE = std::make_unique<GridFunction>(m_xVspace);
    }
}

//! Computes H1 error for solution
std::tuple <double, double, double, double> WaveO2Observer
:: eval_xH1Error (std::shared_ptr<GridFunction>& u,
                  double t) const
{
    double erruL2=0, uEL2=0;
    double erruH1=0, uEH1=0;

    ConstantCoefficient zero(0);
    ConstantCoefficient one(1.0);
    VectorFunctionCoefficient zeroVec(m_xndim, zeroFn);

    // L2 error in solution
    m_uECoeff->SetTime(t);
    m_uE->ProjectCoefficient(*m_uECoeff);
    uEL2 = m_uE->ComputeL2Error(zero);
    erruL2 = u->ComputeL2Error(*m_uECoeff);

    // H1 error in solution
    m_dudxECoeff->SetTime(t);
    uEH1 = m_uE->ComputeH1Error(&zero, &zeroVec, &one,
                                1.0, 1);
    erruH1 = u->ComputeH1Error(m_uECoeff.get(),
                               m_dudxECoeff.get(),
                               &one, 1.0, 1);

    return {std::move(erruL2), std::move(uEL2),
                std::move(erruH1), std::move(uEH1)};
}


// End of file
