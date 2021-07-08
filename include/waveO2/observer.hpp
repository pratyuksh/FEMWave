#ifndef WAVEO2_OBSERVER_HPP
#define WAVEO2_OBSERVER_HPP

#include "mfem.hpp"
using namespace mfem;

#include "../core/config.hpp"
#include "../core/base_observer.hpp"

#include "../mymfem/utilities.hpp"

#include "test_cases.hpp"
#include "coefficients.hpp"
#include "discretisation.hpp"
#include "utilities.hpp"


// Observer for WaveO2 system
class WaveO2Observer : public BaseObserver
{
public:
    //! Constructors
    WaveO2Observer () : BaseObserver () {}

    WaveO2Observer (const nlohmann::json&, int);

    //! Sets the exact solution coefficients
    //! and auxiliary variables needed
    void set (std::shared_ptr<WaveO2TestCases>& testCase,
              std::shared_ptr<WaveO2XtFEM>& discr);

    //void dump_sol (std::shared_ptr<GridFunction>&) const;

    //! Computes H1 error at a given time t
    std::tuple <double, double, double, double>
    eval_xH1Error (std::shared_ptr<GridFunction>&,
                   double t) const;

    //! Computes L2H1 error
    /*std::tuple <double, double, double, double>
    eval_xtL2H1Error (Vector&) const;*/

private:
    bool m_bool_error;
    std::string m_solName_prefix;
    std::string m_solName_suffix;

    int m_xndim;
    std::shared_ptr<WaveO2TestCases> m_testCase;
    std::shared_ptr<WaveO2XtFEM> m_discr;

    mutable std::unique_ptr
    <WaveO2MediumCoeff> m_medCoeff;
    mutable std::unique_ptr
    <WaveO2ExactSolutionCoeff> m_uECoeff;
    mutable std::unique_ptr
    <WaveO2ExactGradtSolutionCoeff> m_dudtECoeff;
    mutable std::unique_ptr
    <WaveO2ExactGradxSolutionCoeff> m_dudxECoeff;

    FiniteElementSpace* m_tVspace = nullptr;
    FiniteElementSpace* m_xVspace = nullptr;

    mutable std::unique_ptr<GridFunction> m_uE;
    mutable std::unique_ptr<GridFunction> m_dudtE;
};


#endif /// WAVEO2_OBSERVER_HPP
