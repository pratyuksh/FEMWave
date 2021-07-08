#include "../../include/core/error.hpp"


std::tuple <double, double, double, double>
ComputeH1Error :: operator() (const GridFunction& u1,
                              const GridFunction& u2,
                              bool has_shared_vertices)
{
    double errorL2 = 0;
    double errorH10 = 0;

    double u2L2 = 0;
    double u2H10 = 0;

    const FiniteElementSpace *fes1 = u1.FESpace();
    const FiniteElementSpace *fes2 = u2.FESpace();

    Mesh *mesh1 = fes1->GetMesh();
    Mesh *mesh2 = fes2->GetMesh();

    int geom_type = fes1->GetFE(0)->GetGeomType();
    int order = 2*fes1->GetOrder(0);
    m_ir = std::make_unique<IntegrationRule>
            (IntRules.Get(geom_type, order));

    m_point_locator = std::make_unique <PointLocator>
            (mesh2, has_shared_vertices);
    int init_cell_id2 = 0;

    int dim = mesh1->Dimension();
    int n = m_ir->GetNPoints();

    Vector gradu1(dim);
    Vector gradu2(dim);
    Vector gradu1mu2(dim);
    for (int i=0; i < mesh1->GetNE(); i++)
    {
        ElementTransformation *trans1
                = mesh1->GetElementTransformation(i);
        auto [elIds2, ips2] = get_ref_points(*trans1, init_cell_id2);
        init_cell_id2 = elIds2[0];

        for (int k=0; k < n; k++)
        {
            const IntegrationPoint &ip1 = m_ir->IntPoint(k);
            const IntegrationPoint &ip2 = ips2[k];

            //std::cout << "\n" << ip1.x << "\t" << ip1.y << std::endl;
            //std::cout << elIds2[k] << "\t"
            //          << ip2.x << "\t" << ip2.y << std::endl;

            double valu1 = u1.GetValue(i, ip1);
            trans1->SetIntPoint(&ip1);
            u1.GetGradient(*trans1, gradu1);

            double valu2 = u2.GetValue(elIds2[k], ip2);
            ElementTransformation *trans2
                    = mesh2->GetElementTransformation(elIds2[k]);
            trans2->SetIntPoint(&ip2);
            u2.GetGradient(*trans2, gradu2);

            double coeff = ip1.weight*trans1->Weight();
            errorL2 += coeff*(valu1-valu2)*(valu1-valu2);
            u2L2 += coeff*valu2*valu2;

            subtract (gradu1, gradu2, gradu1mu2);
            errorH10 += coeff*(gradu1mu2*gradu1mu2);
            u2H10 += coeff*(gradu2*gradu2);
        }
    }

    return {sqrt(errorL2), sqrt(u2L2), sqrt(errorH10), sqrt(u2H10)};
}

std::pair < Array <int>, Array <IntegrationPoint> >
ComputeH1Error :: get_ref_points(ElementTransformation& trans1,
                                 int init_cell_id2)
{
    int n = m_ir->GetNPoints();
    Array <int> elIds2(n);
    Array <IntegrationPoint> ips2(n);
    Vector z;

    for (int k=0; k < n; k++)
    {
        trans1.Transform(m_ir->IntPoint(k), z);
        std::tie (elIds2[k], ips2[k])
                = (*m_point_locator)(z, init_cell_id2);

        init_cell_id2 = elIds2[k];
    }

    return {elIds2, ips2};
}

std::tuple <double, double, double, double>
ComputeH1Error :: test_slow (const GridFunction& u1,
                             const GridFunction& u2)
{
    double errorL2 = 0;
    double errorH10 = 0;

    double u2L2 = 0;
    double u2H10 = 0;

    const FiniteElementSpace *fes1 = u1.FESpace();
    const FiniteElementSpace *fes2 = u2.FESpace();

    Mesh *mesh1 = fes1->GetMesh();
    Mesh *mesh2 = fes2->GetMesh();

    int geom_type = fes1->GetFE(0)->GetGeomType();
    int order = 2*fes1->GetOrder(0);
    m_ir = std::make_unique<IntegrationRule>
            (IntRules.Get(geom_type, order));

    int dim = mesh1->Dimension();
    int n = m_ir->GetNPoints();

    Vector gradu1(dim);
    Vector gradu2(dim);
    Vector gradu1mu2(dim);
    for (int i=0; i < mesh1->GetNE(); i++)
    {
        ElementTransformation *trans1
                = mesh1->GetElementTransformation(i);
        auto [elIds2, ips2] = get_ref_points_slow(*trans1, *mesh2);

        for (int k=0; k < n; k++)
        {
            const IntegrationPoint &ip1 = m_ir->IntPoint(k);
            const IntegrationPoint &ip2 = ips2[k];

            //std::cout << "\n" << ip1.x << "\t" << ip1.y << std::endl;
            //std::cout << elIds2[k] << "\t"
            //          << ip2.x << "\t" << ip2.y << std::endl;

            double valu1 = u1.GetValue(i, ip1);
            trans1->SetIntPoint(&ip1);
            u1.GetGradient(*trans1, gradu1);

            double valu2 = u2.GetValue(elIds2[k], ip2);
            ElementTransformation *trans2
                    = mesh2->GetElementTransformation(elIds2[k]);
            trans2->SetIntPoint(&ip2);
            u2.GetGradient(*trans2, gradu2);

            double coeff = ip1.weight*trans1->Weight();
            errorL2 += coeff*(valu1-valu2)*(valu1-valu2);
            u2L2 += coeff*valu2*valu2;

            subtract (gradu1, gradu2, gradu1mu2);
            errorH10 += coeff*(gradu1mu2*gradu1mu2);
            u2H10 += coeff*(gradu2*gradu2);
        }
    }

    return {sqrt(errorL2), sqrt(u2L2), sqrt(errorH10), sqrt(u2H10)};
}

std::pair < Array <int>, Array <IntegrationPoint> >
ComputeH1Error :: get_ref_points_slow(ElementTransformation& trans1,
                                      Mesh &mesh2)
{
    int dim = trans1.GetSpaceDim();
    int n = m_ir->GetNPoints();

    Vector pt;
    DenseMatrix pts1(dim, n);
    for (int k=0; k < n; k++)
    {
        pts1.GetColumnReference(k, pt);
        trans1.Transform(m_ir->IntPoint(k), pt);
    }

    Array <int> elIds2(n);
    Array <IntegrationPoint> ips2(n);
    mesh2.FindPoints(pts1, elIds2, ips2);

    for (int i=0; i<n; i++) {
        if (elIds2[i] == -1){
            std::cout << i << "\t"
                      << ips2[i].x << "\t"
                      << ips2[i].y << "\t"
                      << pts1(0,i) << "\t"
                      << pts1(1,i) << std::endl;
        }
    }

    return {elIds2, ips2};
}


// End of file
