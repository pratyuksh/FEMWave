#ifndef MYMFEM_ERROR_HPP
#define MYMFEM_ERROR_HPP

#include "mfem.hpp"
using namespace mfem;

#include <Eigen/Core>
#include <Eigen/Dense>

#include <memory>

#include "../mymfem/utilities.hpp"


class ComputeH1Error
{
public:
    std::tuple <double, double, double, double>
    operator()(const GridFunction&, const GridFunction&,
               bool has_shared_vertices = false);

    std::pair < Array <int>, Array <IntegrationPoint> >
    get_ref_points (ElementTransformation&, int init_cell_id2);

    std::tuple <double, double, double, double>
    test_slow(const GridFunction&, const GridFunction&);

    std::pair < Array <int>, Array <IntegrationPoint> >
    get_ref_points_slow (ElementTransformation&, Mesh&);

private:
    std::unique_ptr <const IntegrationRule> m_ir;
    std::unique_ptr <PointLocator> m_point_locator;
};


#endif /// MYMFEM_ERROR_HPP
