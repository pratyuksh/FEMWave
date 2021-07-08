#ifndef WAVEO1_SOLVER_HPP
#define WAVEO1_SOLVER_HPP

#include "mfem.hpp"
using namespace mfem;

#include <Eigen/Core>

#include "../core/config.hpp"
#include "../core/pardiso.hpp"

#include "../mymfem/utilities.hpp"

#include "test_cases_factory.hpp"
#include "coefficients.hpp"
#include "discretisation.hpp"
#include "observer.hpp"


class WaveO1Solver
{
public:
    WaveO1Solver (const nlohmann::json& config);

    WaveO1Solver (const nlohmann::json& config,
                std::shared_ptr<WaveO1TestCases>& testCase);

    WaveO1Solver (const nlohmann::json& config,
                std::string mesh_dir,
                const int lx,
                const int lt,
                const bool load_init_mesh=false);

    WaveO1Solver (const nlohmann::json& config,
                std::shared_ptr<WaveO1TestCases>& testCase,
                std::string mesh_dir,
                const int lx,
                const int lt,
                const bool load_init_mesh=false);

    WaveO1Solver (const nlohmann::json& config,
                std::string mesh_dir,
                const int lx,
                const bool load_init_mesh=false);

    WaveO1Solver (const nlohmann::json& config,
                std::shared_ptr<WaveO1TestCases>& testCase,
                std::string mesh_dir,
                const int lx,
                const bool load_init_mesh=false);

    ~ WaveO1Solver ();

    // sets meshes, discretization, observer
    void set(std::string mesh_dir,
             const int lx,
             const int lt=-2,
             const bool load_init_mesh=false);

    void set(std::string mesh_dir,
             const bool load_init_mesh=false);
    void set(int lx, int lt=-2);
    void set_meshes(int lx, int lt=-2);

    // used for the full-tensor version
    std::tuple<int, double, double, Eigen::VectorXd>
    operator()();

    // used for the sparse-grids handler
    void operator()(std::unique_ptr<BlockVector>&);
    
    void run();
    void init ();
    void solve (BlockVector *W);

private:
    void assemble_system();

    void solve_one_time_slab (const int n, BlockVector *W);
    void assemble_rhs_initial();
    void assemble_rhs(int);
    void copy_slabW_to_W(const int n, const BlockVector *,
                         BlockVector *) const;

public:
    void finalize () const;

    // error computation
    Eigen::VectorXd compute_error(BlockVector *W) const;

    // energy computation
    void compute_energy(BlockVector *W) const;

    // solution at end-time
    std::shared_ptr <GridFunction>
    get_pressureSol_at_endTime (BlockVector *W) const;

    std::shared_ptr <GridFunction>
    get_velocitySol_at_endTime (BlockVector *W) const;

    void set_pressureSol_at_endTime ();

    void set_velocitySol_at_endTime ();

    void dump_sol_at_time_steps(BlockVector *W) const;

    void dump_sol_at_time_stamps(BlockVector *W) const;

    void measure_signal (BlockVector *W) const;

    // for L2-projection
    std::tuple<int, double, double, Eigen::VectorXd>
    projection();
    void projection(std::unique_ptr<BlockVector>&);

    void init_projection ();

private:
    void assemble_projection_system();
    void assemble_projection_rhs(int) const;

public:
    void solve_projection (BlockVector *W);

    inline std::shared_ptr<Mesh> get_tMesh() const {
        return m_tMesh;
    }

    inline std::shared_ptr<Mesh> get_xMesh() const {
        return m_xMesh;
    }

    inline std::shared_ptr<WaveO1TestCases>
    get_test_case() const {
        return m_testCase;
    }

    inline std::shared_ptr<WaveO1XtDG>
    get_discretisation() const {
        return m_discr;
    }

    inline std::shared_ptr<WaveO1Observer>
    get_observer() const {
        return m_observer;
    }

    inline std::pair<double, double>
    get_meshChars() const
    {
        double ht_min, ht_max, kappat_min, kappat_max;
        m_tMesh->GetCharacteristics(ht_min, ht_max,
                                    kappat_min, kappat_max);

        double hx_min, hx_max, kappax_min, kappax_max;
        m_xMesh->GetCharacteristics(hx_min, hx_max,
                                    kappax_min, kappax_max);

        return {std::move(ht_max), std::move(hx_max)};
    }

    inline int get_ndofs() const {
        return std::move(m_W->Size());
    }

private:
    const nlohmann::json& m_config;
    std::string m_mesh_dir;
    std::string m_mesh_elem_type;

    double m_endTime;
    std::shared_ptr<Mesh> m_tMesh;
    std::shared_ptr<Mesh> m_xMesh;
    bool m_load_init_mesh;
    int m_lx0;

    std::shared_ptr<WaveO1TestCases> m_testCase;
    std::shared_ptr<WaveO1XtDG> m_discr;
    std::shared_ptr<WaveO1Observer> m_observer;

    bool m_bool_error;
    std::string m_errorType;

    std::string m_meshFormat;

    std::string m_solOutputAt;
    Vector m_timeStamps;

    bool m_bool_measureSignal;
    Vector m_xMeasurementPoint;

    bool m_bool_energy;

    Array<int> m_blockOffsets;
    Array<int> m_slab_blockOffsets;

    BlockOperator *m_waveOp = nullptr;
    SparseMatrix *m_waveMat = nullptr;
    SparseMatrix *m_waveProjMat = nullptr;

    std::unique_ptr <BlockVector> m_W;
    mutable std::unique_ptr <BlockVector> m_slabW;
    mutable std::unique_ptr <BlockVector> m_slabB;

    std::shared_ptr <GridFunction> m_pressure;
    std::shared_ptr <GridFunction> m_velocity;

#ifdef PARDISO_HPP
    std::unique_ptr<PardisoSolver> m_pardisoSolver;
#endif
    mutable bool m_pardiso_finalized = true;
};


#endif /// WaveO1_SOLVER_HPP
