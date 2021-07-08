#ifndef MYMFEM_UTILITIES_HPP
#define MYMFEM_UTILITIES_HPP

#include "mfem.hpp"
using namespace mfem;

#include <memory>


//! Releases SparseMatrix memory and seta pointer to null
void clear (SparseMatrix* &mat);

//! Zero function
void zeroFn(const Vector&, Vector&);

//! Returns the distance of the mesh element i
//! to a given point
double get_element_distance_to_point
(std::shared_ptr<Mesh>& mesh, int elId, Vector& pt);

//! Returns the reference coordinates of a point on a line segment
//double get_pointRefCoordsInLineSegment();

//! Returns the upper-triangular, including the diagonal,
//! for the input matrix A
SparseMatrix& get_upper_triangle(const SparseMatrix& A);

//! Searches the mesh-elements where a given physical point lies
class PointLocator
{
public:
    PointLocator (Mesh *mesh, bool has_shared_vertices=false)
        : m_has_shared_vertices (has_shared_vertices),
          m_mesh (mesh)
    {
        m_invTr = std::make_unique
                <InverseElementTransformation>();

        if (!m_has_shared_vertices) {
            m_vToEl.reset(m_mesh->GetVertexToElementTable());
        }
        else {
            m_shared_vertices = std::make_unique<Table>();
            get_shared_vertices_table();

            m_vToEl = std::make_unique<Table>();
            get_vertices_to_elements_table();
        }
    }

    std::pair <Array<int>, Array<IntegrationPoint>>
    operator() (const DenseMatrix&, int&) const;

    std::pair <int, IntegrationPoint> operator()
    (const Vector&, const int) const;

    void get_shared_vertices_table () const;

    void get_vertices_to_elements_table () const;

    inline std::pair <int, IntegrationPoint>
    compute_ref_point (const Vector& x, const int elId) const
    {
        IntegrationPoint ip;
        m_invTr->SetTransformation
                (*m_mesh->GetElementTransformation(elId));
        auto info = m_invTr->Transform(x, ip);
        return {info, ip};
    }

private:
    double m_TOL = 1E-12;
    bool m_has_shared_vertices = false;

    Mesh *m_mesh = nullptr;
    std::unique_ptr<InverseElementTransformation> m_invTr;

    std::unique_ptr<Table> m_vToEl;
    std::unique_ptr<Table> m_shared_vertices;
};

//! Searches the mesh-faces where a given physical point lies
class MeshFaceLocator
{
public:
    MeshFaceLocator (Mesh *mesh) : m_mesh (mesh) {}

    int operator() (const Vector&) const;

private:
    double m_TOL = 1E-14;
    Mesh *m_mesh = nullptr;
};

#endif /// MYMFEM_UTILITIES_HPP
