#include "../../include/mymfem/utilities.hpp"

#include <assert.h>


//! Releases SparseMatrix memory and seta pointer to null
void clear (SparseMatrix* &mat) {
    delete mat;
    mat = nullptr;
}

//! Zero function
void zeroFn(const Vector&, Vector& f) {
   f = 0.0;
}

//! Returns the distance of the mesh element i
//! to a given point
double get_element_distance_to_point
(std::shared_ptr<Mesh>& mesh, int elId, Vector& pt)
{
    double dist = 1E+2;
    Array<int> v;

    mesh->GetElementVertices (elId, v);
    for (int i=0; i<v.Size(); i++)
    {
        double *data = mesh->GetVertex(v[i]);
        Vector coords(data, mesh->Dimension());
        dist = std::min(coords.DistanceTo(pt.GetData()),
                        dist);
    }

    return dist;
}

//! Returns the upper-triangular, including the diagonal,
//! for the input matrix A
SparseMatrix& get_upper_triangle(const SparseMatrix& A)
{
    int numRows = A.NumRows();
    int numCols = A.NumCols();
    const int *iA = A.GetI();
    const int *jA = A.GetJ();
    const double *dA = A.GetData();

    SparseMatrix* UA = new SparseMatrix(numRows, numCols);
    for (int i=0; i<numRows; i++) {
        for(int k=iA[i]; k<iA[i+1]; k++) {
            int j = jA[k];
            if (j >= i) { UA->_Add_(i, j, dA[k]); }
        }
    }
    UA->Finalize();

    return *UA;
}


//! Locates a given physical point x to a mesh element
//! using the initial guess init_elId
std::pair <Array<int>, Array<IntegrationPoint>>
PointLocator :: operator() (const DenseMatrix& X,
                            int& init_elId) const
{
    int numPoints = X.NumCols();
    Array<int> elIds(numPoints);
    Array<IntegrationPoint> ips(numPoints);

    Vector x(X.NumRows());
    for (int i=0; i<numPoints; i++)
    {
        X.GetColumn(i, x);
        std::tie (elIds[i], ips[i]) = (*this)(x, init_elId);
        init_elId = elIds[i];
    }

    return {elIds, ips};
}

//! Locates a given physical point x to a mesh element
//! using the initial guess init_elId
std::pair <int, IntegrationPoint>
PointLocator :: operator() (const Vector& x,
                            const int init_elId) const
{
    auto [info_, ip_] = compute_ref_point(x, init_elId);
    if (info_ == InverseElementTransformation::Inside) {
        return {init_elId, ip_};
    }

    bool found = false;
    int info;
    int elId = init_elId, old_elId = init_elId;
    IntegrationPoint ip;

    Array <int> vertices;
    Vector z (m_mesh->Dimension());
    double min_dist = std::numeric_limits<double>::max();

    while(!found)
    {
        /// find the element closest to point x amongst the
        /// neighbours of the vertices of the current element
        m_mesh->GetElementVertices(elId, vertices);
        for (int i=0; i < vertices.Size(); i++)
        {
            int v = vertices[i];
            int ne = m_vToEl->RowSize(v);
            const int* els = m_vToEl->GetRow(v);
            //std::cout << "\n\nFor vertex: " << i << std::endl;
            for (int j=0; j<ne; j++)
            {
                if (els[j] == elId) {continue;}

                m_mesh->GetElementTransformation(els[j])
                        ->Transform(Geometries.GetCenter
                         (m_mesh->GetElementBaseGeometry
                          (els[j])), z);
                double dist = z.DistanceTo(x.GetData());

                //std::cout << els[j] << "\t" << elId << "\t"
                //          << dist << "\t"
                //          << min_dist << std::endl;
                if (dist < min_dist)
                {
                    min_dist = dist;
                    elId = els[j];
                    //std::cout << "Minimum: "
                    //          << min_dist << "\t"
                    //          << elId << std::endl;
                }
            }
        }

        if (elId == old_elId) {
            /// if elId is not upadted, search
            /// all neighbours of its vertices
            for (int i=0; i < vertices.Size(); i++)
            {
                int v = vertices[i];
                int ne = m_vToEl->RowSize(v);
                const int* els = m_vToEl->GetRow(v);

                for (int j=0; j<ne; j++)
                {
                    if (els[j] == elId) {continue;}
                    //std::cout << "Search neighbours: "
                    //          << x(0) << "\t"
                    //          << x(1) << "\t"
                    //          << els[j] << std::endl;
                    std::tie (info, ip) = compute_ref_point(x, els[j]);
                    if (info == InverseElementTransformation::Inside) {
                        //std::cout << "Exception: "
                        //          << els[j] << std::endl;
                        return {els[j],ip};
                    }
                }
            }
            /// in case the neighbour search fails,
            /// loop over other elements
            //if (elId < m_mesh->GetNE()-1) { elId++;}
            //else { elId = 0; }
        }
        else {
            /// if elId is upadted, check if it contains point x
            std::tie (info, ip) = compute_ref_point(x, elId);
            if (info == InverseElementTransformation::Inside) {
                found = true;
            }
            old_elId = elId;
        }
    }

    return {elId, ip};
}

//! Generates a table of vertices,
//! which are shared by different processors.
//! Needed when a ParMesh is written to one file,
//! the same vertices at the partition interfaces
//! will have a different numbering for each processor
void PointLocator :: get_shared_vertices_table() const
{
    int nv = m_mesh->GetNV();
    int dim = m_mesh->Dimension();

    m_shared_vertices->MakeI (nv);
    for (int i=0; i<nv; i++)
    {
        Vector v1(m_mesh->GetVertex(i), dim);
        for (int j=i+1; j<nv; j++)
        {
            if (v1.DistanceTo(m_mesh->GetVertex(j)) < m_TOL) {
                m_shared_vertices->AddAColumnInRow(i);
                m_shared_vertices->AddAColumnInRow(j);
            }
        }
    }
    m_shared_vertices->MakeJ();

    for (int i=0; i<nv; i++)
    {
        Vector v1(m_mesh->GetVertex(i), dim);
        for (int j=i+1; j<nv; j++)
        {
            if (v1.DistanceTo(m_mesh->GetVertex(j)) < m_TOL) {
                m_shared_vertices->AddConnection(i, j);
                m_shared_vertices->AddConnection(j, i);
            }
        }
    }
    m_shared_vertices->ShiftUpI();
}

//! Generates the vertex to element table,
//! uses the shared vertices table
void PointLocator :: get_vertices_to_elements_table() const
{
    int ne = m_mesh->GetNE();
    int nv = m_mesh->GetNV();

    m_vToEl->MakeI(nv);

    Array<int> v, sv;
    for (int i=0; i<ne; i++)
    {
        m_mesh->GetElementVertices(i, v);
        for (int j=0; j<v.Size(); j++) {
            m_vToEl->AddAColumnInRow(v[j]);

            int nsv = m_shared_vertices->RowSize(v[j]);
            m_shared_vertices->GetRow(v[j], sv);
            for (int k=0; k<nsv; k++) {
                m_vToEl->AddAColumnInRow(sv[k]);
            }
        }
    }
    m_vToEl->MakeJ();

    for (int i=0; i<ne; i++)
    {
        m_mesh->GetElementVertices(i, v);
        for (int j=0; j<v.Size(); j++) {
            m_vToEl->AddConnection(v[j], i);

            int nsv = m_shared_vertices->RowSize(v[j]);
            m_shared_vertices->GetRow(v[j], sv);
            for (int k=0; k<nsv; k++) {
                m_vToEl->AddConnection(sv[k], i);
            }
        }
    }
    m_vToEl->ShiftUpI();
}

//! Locates a given physical point x to a mesh face
int MeshFaceLocator :: operator()
(const Vector& x) const
{
    int faceElId = -1;
    Array<int> faceVertices(2);
    for (int i=0; i<m_mesh->GetNumFaces(); i++)
    {
        m_mesh->GetFaceVertices(i, faceVertices);
        Vector v1(m_mesh->GetVertex(faceVertices[0]),2);
        Vector v2(m_mesh->GetVertex(faceVertices[1]),2);

        auto distV1ToX = v1.DistanceTo(x);
        auto distV2ToX = v2.DistanceTo(x);
        auto distV1ToV2 = v1.DistanceTo(v2);
        if (std::fabs(distV1ToX + distV2ToX - distV1ToV2) <= m_TOL)
        {
            faceElId = i;
            break;
        }
    }

    return faceElId;
}

// End of file
