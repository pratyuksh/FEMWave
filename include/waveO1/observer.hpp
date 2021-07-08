#ifndef WAVEO1_OBSERVER_HPP
#define WAVEO1_OBSERVER_HPP

#include "mfem.hpp"
using namespace mfem;

#include "../core/config.hpp"
#include "../core/base_observer.hpp"

#include "../mymfem/utilities.hpp"

#include "test_cases.hpp"
#include "coefficients.hpp"
#include "discretisation.hpp"
#include "utilities.hpp"


// Observer for WaveO1 system
class WaveO1Observer : public BaseObserver
{
public:
    WaveO1Observer () : BaseObserver () {}

    WaveO1Observer (const nlohmann::json&, int);

    void set (std::shared_ptr<WaveO1TestCases>& testCase,
              std::shared_ptr<WaveO1XtDG>& discr);
    
    //! Dumps solution
    void dump_sol (std::shared_ptr<GridFunction>&,
                   std::shared_ptr<GridFunction>&) const;

    //! Dumps solution at given time ids
    void dump_sol_at_time_steps
    (std::shared_ptr<GridFunction>&,
     std::shared_ptr<GridFunction>&,
     int time_step_id) const;

    //! Dumps solution at given time stamps
    void dump_sol_at_time_stamps
    (std::shared_ptr<GridFunction>&,
     std::shared_ptr<GridFunction>&,
     double t) const;

    //! Dumps pressure signal
    void dump_pressure_signal (Vector& time_series,
                               Vector& signal) const;

    //! Dumps energy
    void dump_energy (Vector& time_series,
                      Vector& energy) const;
    
    std::tuple <double, double, double, double>
    eval_xL2Error (std::shared_ptr<GridFunction>&,
                   std::shared_ptr<GridFunction>&,
                   double t) const;

    std::tuple <double, double, double, double>
    eval_xtL2Error (BlockVector&) const;

    //! Computes DG error
    double eval_dgFtimeError
    (std::shared_ptr<GridFunction>&,
     std::shared_ptr<GridFunction>&,
     double t) const;

    double eval_dgFspaceError
    (std::shared_ptr<GridFunction>&,
     std::shared_ptr<GridFunction>&,
     double t,
     std::shared_ptr<WaveO1InvSqMediumCoeff>&) const;

    double eval_dgFspaceError
    (std::shared_ptr<GridFunction>&,
     std::shared_ptr<GridFunction>&,
     std::shared_ptr<GridFunction>&,
     std::shared_ptr<GridFunction>&,
     double t,
     std::shared_ptr
     <WaveO1InvSqMediumCoeff>&) const;

    std::tuple <double, double>
    eval_xtDgError (BlockVector&) const;

    //! Computes energy functional
    double eval_dgFtimeJump
    (std::shared_ptr<GridFunction>&,
     std::shared_ptr<GridFunction>&) const;

    double eval_dgFspaceJump
    (int degp, int degv, Mesh* mesh, double t,
     std::shared_ptr<WaveO1InvSqMediumCoeff>&) const;

    double eval_dgFspaceJump
    (std::shared_ptr<GridFunction>&,
     std::shared_ptr<GridFunction>&,
     std::shared_ptr
     <WaveO1InvSqMediumCoeff>&) const;

    double eval_dgFspaceJump
    (std::shared_ptr<GridFunction>&,
     std::shared_ptr<GridFunction>&,
     std::shared_ptr<GridFunction>&,
     std::shared_ptr<GridFunction>&,
     std::shared_ptr
     <WaveO1InvSqMediumCoeff>&) const;

    Vector eval_energy(BlockVector&) const;

    //! Computes DG^{+} error
    double eval_dgpFtimeError
    (std::shared_ptr<GridFunction>&,
     std::shared_ptr<GridFunction>&,
     double t) const;

    double eval_dgpFspaceError
    (std::shared_ptr<GridFunction>&,
     std::shared_ptr<GridFunction>&,
     double t,
     std::shared_ptr<WaveO1InvSqMediumCoeff>&) const;

    std::tuple <double, double>
    eval_xtDgpError (BlockVector&) const;
    
private:
    bool m_bool_error;
    std::string m_errorType;

    bool m_bool_energy;
    bool m_bool_dumpEnergy;

    bool m_bool_dumpSignal;
    std::string m_solName_prefix;
    std::string m_solName_suffix;

    int m_dump_out_nsteps;
    int m_stabParamsType;

    int m_xndim;
    std::shared_ptr<WaveO1TestCases> m_testCase;
    std::shared_ptr<WaveO1XtDG> m_discr;

    std::unique_ptr<WaveO1ExactPressureCoeff> m_pE_coeff;
    std::unique_ptr<WaveO1ExactVelocityCoeff> m_vE_coeff;

    FiniteElementSpace* m_tWspace = nullptr;
    FiniteElementSpace* m_xW1space = nullptr;
    FiniteElementSpace* m_xW2space = nullptr;

    std::unique_ptr<GridFunction> m_pE;
    std::unique_ptr<GridFunction> m_vE;

    std::unique_ptr<mymfem::AssembleDgError>
    m_assembleDgError;

    std::unique_ptr<mymfem::AssembleDgPlusSmoothError>
    m_assembleDgpSmoothError;

    std::unique_ptr<mymfem::AssembleDgJumps>
    m_assembleDgJumps;
};


#endif /// WaveO1_OBSERVER_HPP
