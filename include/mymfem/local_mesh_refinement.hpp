#ifndef MYMFEM_LOCAL_MESH_REFINEMENT_HPP
#define MYMFEM_LOCAL_MESH_REFINEMENT_HPP

#include <memory>

#include "../../include/mymfem/geometry.hpp"

#include "mfem.hpp"
using namespace mfem;


class LocalMeshRefinement
{
public:
    //! Constructor
    LocalMeshRefinement (std::shared_ptr<Polygon>& polygon)
        : m_polygon(polygon)
    {
        m_distR0 = eval_R0();
    }

    //! Evaluates R0 distance
    //! for all corners of the polygon
    Array<double> eval_R0();

    //! uniform mesh refinement
    void uniform(std::shared_ptr<Mesh>&, double h0);

    //! local mesh refinement
    void local(std::shared_ptr<Mesh>&,
               int deg, double h0,
               bool nonConforming=true);

protected:
    std::shared_ptr<Polygon> m_polygon;
    double m_factor = 0.49; // must be < 0.5
    Array<double> m_distR0;
};


#endif /// MYMFEM_LOCAL_MESH_REFINEMENT_HPP
