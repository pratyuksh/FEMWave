#ifndef WAVEO1_SPARSE_GRIDS_HANDLER_HPP
#define WAVEO1_SPARSE_GRIDS_HANDLER_HPP

#include "mfem.hpp"
using namespace mfem;

#include <Eigen/Core>

#include "../core/config.hpp"
#include "../mymfem/utilities.hpp"

#include "discretisation.hpp"
#include "observer.hpp"
#include "solver.hpp"

class SparseGridsHandler
{
public:
    SparseGridsHandler (const nlohmann::json& config);

    SparseGridsHandler
    (const nlohmann::json& config,
     std::shared_ptr<WaveO1TestCases>& testCase);

    SparseGridsHandler (const nlohmann::json& config,
                        std::string mesh_dir,
                        const int Lx,
                        const int L0x,
                        const int L0t,
                        bool load_init_mesh=false);

    SparseGridsHandler
    (const nlohmann::json& config,
     std::shared_ptr<WaveO1TestCases>& testCase,
     std::string mesh_dir,
     const int Lx, const int L0x, const int L0t,
     bool load_init_mesh=false);

    ~SparseGridsHandler ();

    void set(std::string mesh_dir,
             const bool load_init_mesh=false);
    void set(const int Lx);
    void set(const int Lx, const int L0x, const int L0t);

    // sets minimum resolution levels
    void set_minLevels(const int L0x, const int L0t);

    void reset();

    std::tuple<int, double, double, Eigen::VectorXd>
    operator() ();

    void run ();

    // point locators
    void set_xPointLocators();
    void set_tPointLocators();

    // mesh face
    void set_xMeshFaceLocators();

    // error computation
    Eigen::VectorXd compute_error();

    void init_eval_error();

    //! Computes L2 error
    std::tuple <double, double, double, double>
    eval_xL2Error(Array<GridFunction *> pColl,
                  Array<GridFunction *> vColl) const;

    std::tuple <double, double, double, double>
    eval_xL2Error_at_endTime() const;

    std::tuple <double, double, double, double>
    eval_xtL2Error() const;

    //! Computes L2-projection
    std::tuple<int, double, double, Eigen::VectorXd>
    compute_projection();

    //! Computes DG error
    double eval_dgFtimeError(Array<GridFunction *>,
                             Array<GridFunction *>,
                             double t) const;

    double eval_dgFspaceError
    (Array<GridFunction *>,
     Array<GridFunction *>,
     double t,
     std::shared_ptr<WaveO1InvSqMediumCoeff>&) const;

    double eval_dgFspaceError
    (Array<GridFunction *>,
     Array<GridFunction *>,
     Array<GridFunction *>,
     Array<GridFunction *>,
     double t,
     std::shared_ptr<WaveO1InvSqMediumCoeff>&) const;

    std::tuple <double, double>
    eval_xtDgError() const;

    //! Computes DG+ error
    double eval_dgpFtimeError(Array<GridFunction *>,
                             Array<GridFunction *>,
                             double t) const;

    double eval_dgpFspaceError
    (Array<GridFunction *>,
     Array<GridFunction *>,
     double t,
     std::shared_ptr<WaveO1InvSqMediumCoeff>&) const;

    std::tuple <double, double>
    eval_xtDgpError() const;

    inline std::pair<double, double>
    get_meshChars() const
    {
        double ht_min, ht_max, kappat_min, kappat_max;
        m_tMesh_finest->GetCharacteristics(ht_min,
                                           ht_max,
                                           kappat_min,
                                           kappat_max);

        double hx_min, hx_max, kappax_min, kappax_max;
        m_xMesh_finest->GetCharacteristics(hx_min,
                                           hx_max,
                                           kappax_min,
                                           kappax_max);

        return {std::move(ht_max), std::move(hx_max)};
    }

    inline int get_ndofs() const
    {
        int ndofs = 0;
        for (int k=0; k<m_numSubSols; k++) {
            ndofs += m_subSols[k]->Size();
        }
        return std::move(ndofs);
    }

private:
    const nlohmann::json& m_config;
    std::string m_mesh_dir;
    bool m_load_init_mesh;

    int m_stabParamsType;

    int m_xndim;
    double m_endTime;
    bool m_bool_error = false;
    std::string m_errorType;

    int m_numSubSols;
    int m_Lx, m_L0x, m_L0t;
    int m_numUniqueLevels;

    Array2D<int> m_levels;
    Array<int> m_combCoeffs;

    Array<std::unique_ptr<BlockVector>> m_subSols;
    Array<std::shared_ptr<Mesh>> m_tMeshes;
    Array<std::shared_ptr<Mesh>> m_xMeshes;
    Array<std::shared_ptr<WaveO1XtDG>> m_discrs;
    Array<std::shared_ptr<WaveO1Observer>> m_obsrvs;

    std::shared_ptr<WaveO1TestCases> m_testCase;
    std::unique_ptr<WaveO1ExactPressureCoeff> m_pE_coeff;
    std::unique_ptr<WaveO1ExactVelocityCoeff> m_vE_coeff;

    std::shared_ptr<Mesh> m_tMesh_finest;
    std::shared_ptr<Mesh> m_xMesh_finest;

    FiniteElementSpace* m_tWspace_finest = nullptr;
    FiniteElementSpace* m_xW1space_finest = nullptr;
    FiniteElementSpace* m_xW2space_finest = nullptr;

    std::shared_ptr<GridFunction> m_pE;
    std::shared_ptr<GridFunction> m_vE;

    Array<PointLocator *> m_xPointLocators;
    Array<PointLocator *> m_tPointLocators;

    Array<MeshFaceLocator *> m_xMeshFaceLocators;
};

#endif /// WaveO1_SPARSE_GRIDS_HANDLER_HPP
