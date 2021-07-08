#include "../../include/waveO1/utilities.hpp"


// build the solution in a space-time element
// tVdofs gives the indices for reading
// the space-time solution uXtSol
// for the space-time element required.
// assumes that uSol size is preset
void build_xSol_FG(const Vector& uXtSol,
                   const Vector& tShape,
                   Array<int> tVdofs,
                   Vector& uSol)
{
    int tNdofs = tShape.Size();
    int xdimV = uSol.Size();

    uSol = 0.0;
    for (int j=0; j<tNdofs; j++){
        int shift = tVdofs[j]*xdimV;
        for (int k=0; k<xdimV; k++) {
            uSol(k) += tShape(j)*uXtSol(k + shift);
        }
    }
}

// Assembles DG error for scalar FE space
// Time-like face
double mymfem::AssembleDgError
:: Ftime (GridFunction* u, Coefficient* uE)
{
    double errDg = 0;

    auto fes = u->FESpace();
    auto mesh = fes->GetMesh();

    FaceElementTransformations *ftr = nullptr;
    IntegrationPoint eip1, eip2;

    for (int i=0; i<mesh->GetNumFaces(); i++)
    {
        ftr = mesh -> GetInteriorFaceTransformations (i);
        if (ftr != nullptr)
        {
            int i1 = ftr->Elem1No;
            int i2 = ftr->Elem2No;
            for (int j=0; j<m_ir->GetNPoints(); j++)
            {
                const IntegrationPoint &ip
                        = m_ir->IntPoint(j);
                ftr->Face->SetIntPoint(&ip);

                ftr->Loc1.Transform(m_ir->IntPoint(j), eip1);
                ftr->Loc2.Transform(m_ir->IntPoint(j), eip2);

                double val = uE->Eval(*ftr->Face, ip)
                        - u->GetValue(i1, eip1);

                val -= uE->Eval(*ftr->Face, ip)
                        - u->GetValue(i2, eip2);

                double coeff = ip.weight;
                if (m_stabParamsType == 1
                        || m_stabParamsType == 3) {
                    coeff *= ftr->Face->Weight();
                }
                errDg += coeff*val*val;
            }
        }
    }

    // Dirichlet boundary faces
    int nbfaces = mesh->GetNBE();
    for (int i = 0; i < nbfaces; i++)
    {
        const int bdr_attr = mesh->GetBdrAttribute(i);
        if (m_xEss_bdr_marker[bdr_attr-1] == 0) { continue; }

        ftr = mesh -> GetBdrFaceTransformations(i);

        // face element 1
        int i1 = ftr->Elem1No;
        for (int j=0; j<m_ir->GetNPoints(); j++)
        {
            const IntegrationPoint &ip
                    = m_ir->IntPoint(j);
            ftr->Face->SetIntPoint(&ip);

            ftr->Loc1.Transform(m_ir->IntPoint(j), eip1);
            double val = uE->Eval(*ftr->Face, ip)
                    - u->GetValue(i1, eip1);

            double coeff = ip.weight;
            if (m_stabParamsType == 1
                    || m_stabParamsType == 3) {
                coeff *= ftr->Face->Weight();
            }
            errDg += coeff*val*val;
        }
    }

    return errDg;
}


// Assembles DG error for vector FE space
// Time-like face
double mymfem::AssembleDgError
:: Ftime (GridFunction* u, VectorCoefficient* uE)
{
    double errDg = 0;

    auto fes = u->FESpace();
    auto mesh = fes->GetMesh();
    int dim = mesh->Dimension();

    Vector uSol, uESol;
    DenseMatrix errVal;

    FaceElementTransformations *ftr = nullptr;
    IntegrationPoint eip1, eip2;
    Vector nor(dim);

    errVal.SetSize(dim, m_ir->GetNPoints());
    for (int i=0; i<mesh->GetNumFaces(); i++)
    {
        ftr = mesh -> GetInteriorFaceTransformations (i);
        if (ftr != nullptr)
        {
            int i1 = ftr->Elem1No;
            int i2 = ftr->Elem2No;

            Vector buf(dim);
            for (int j=0; j<m_ir->GetNPoints(); j++)
            {
                const IntegrationPoint &ip
                        = m_ir->IntPoint(j);
                ftr->Face->SetIntPoint(&ip);

                CalcOrtho(ftr->Face->Jacobian(), nor);
                nor /= ftr->Face->Weight(); // unit normal

                ftr->Loc1.Transform(m_ir->IntPoint(j), eip1);
                ftr->Loc2.Transform(m_ir->IntPoint(j), eip2);

                uE->Eval(uESol, *ftr->Face, ip);
                u->GetVectorValue(i1, eip1, uSol);
                for (int k=0; k<dim; k++)
                    buf(k) = uESol(k) - uSol(k);

                u->GetVectorValue(i2, eip2, uSol);
                for (int k=0; k<dim; k++)
                    buf(k) -= uESol(k) - uSol(k);

                double coeff = ip.weight*ftr->Face->Weight();
                if (m_stabParamsType == 3
                        || m_stabParamsType == 4) {
                    coeff *= ftr->Face->Weight();
                }
                errDg += coeff*(buf*nor)*(buf*nor);
            }
        }
    }

    // Neumann boundary faces
    int nbfaces = mesh->GetNBE();
    for (int i = 0; i < nbfaces; i++)
    {
        const int bdr_attr = mesh->GetBdrAttribute(i);
        if (m_xNat_bdr_marker[bdr_attr-1] == 0) { continue; }

        ftr = mesh -> GetBdrFaceTransformations(i);
        int i1 = ftr->Elem1No;

        Vector buf(dim);
        for (int j=0; j<m_ir->GetNPoints(); j++)
        {
            const IntegrationPoint &ip
                    = m_ir->IntPoint(j);
            ftr->Face->SetIntPoint(&ip);

            CalcOrtho(ftr->Face->Jacobian(), nor);
            nor /= ftr->Face->Weight(); // unit normal

            ftr->Loc1.Transform(m_ir->IntPoint(j), eip1);

            uE->Eval(uESol, *ftr->Face, ip);
            u->GetVectorValue(i1, eip1, uSol);
            for (int k=0; k<dim; k++)
                buf(k) = uESol(k) - uSol(k);

            double coeff = ip.weight*ftr->Face->Weight();
            if (m_stabParamsType == 3
                    || m_stabParamsType == 4) {
                coeff *= ftr->Face->Weight();
            }
            errDg += coeff*(buf*nor)*(buf*nor);
        }
    }

    return errDg;
}

// Assembles DG error for scalar FE space
// Time-like face
double mymfem::AssembleDgError
:: Fspace (GridFunction* u, Coefficient* uE,  Coefficient* medium)
{
    double errDg = 0;

    auto fes = u->FESpace();
    auto mesh = fes->GetMesh();

    double errVal;
    ElementTransformation *tr = nullptr;

    for (int i=0; i<mesh->GetNE(); i++)
    {
        tr = mesh->GetElementTransformation(i);

        for (int j=0; j<m_ir->GetNPoints(); j++)
        {
            const IntegrationPoint &ip = m_ir->IntPoint(j);
            tr->SetIntPoint(&ip);

            errVal = uE->Eval(*tr, ip) - u->GetValue(i, ip);
            errDg += ip.weight*tr->Weight()*medium->Eval(*tr, ip)
                    *(errVal*errVal);
        }
    }

    return errDg;
}

double mymfem::AssembleDgError
:: Fspace (GridFunction* u1, GridFunction* u2, Coefficient*,
           Coefficient* medium)
{
    double errDg = 0;

    auto fes = u1->FESpace();
    auto mesh = fes->GetMesh();

    double errVal;
    ElementTransformation *tr = nullptr;

    for (int i=0; i<mesh->GetNE(); i++)
    {
        tr = mesh->GetElementTransformation(i);

        for (int j=0; j<m_ir->GetNPoints(); j++)
        {
            const IntegrationPoint &ip = m_ir->IntPoint(j);
            tr->SetIntPoint(&ip);

            errVal = u1->GetValue(i, ip)
                    - u2->GetValue(i, ip);
            errDg += ip.weight*tr->Weight()*medium->Eval(*tr, ip)
                    *(errVal*errVal);
        }
    }

    return errDg;
}

double mymfem::AssembleDgError
:: Fspace (GridFunction* u, VectorCoefficient* uE)
{
    double errDg = 0;

    auto fes = u->FESpace();
    auto mesh = fes->GetMesh();
    int dim = mesh->Dimension();

    Vector uSol(dim), uESol(dim);
    Vector errVal(dim);
    ElementTransformation *tr = nullptr;

    for (int i=0; i<mesh->GetNE(); i++)
    {
        tr = mesh->GetElementTransformation(i);

        for (int j=0; j<m_ir->GetNPoints(); j++)
        {
            const IntegrationPoint &ip = m_ir->IntPoint(j);
            tr->SetIntPoint(&ip);

            uE->Eval(uESol, *tr, ip);
            u->GetVectorValue(i, ip, uSol);
            for (int k=0; k<dim; k++)
                errVal(k) = uESol(k) - uSol(k);
            errDg += ip.weight*tr->Weight()*(errVal*errVal);
        }
    }

    return errDg;
}

double mymfem::AssembleDgError
:: Fspace (GridFunction* u1, GridFunction* u2,
           VectorCoefficient*)
{
    double errDg = 0;

    auto fes = u1->FESpace();
    auto mesh = fes->GetMesh();
    int dim = mesh->Dimension();

    Vector uSol1(dim), uSol2(dim);
    Vector errVal(dim);
    ElementTransformation *tr = nullptr;

    for (int i=0; i<mesh->GetNE(); i++)
    {
        tr = mesh->GetElementTransformation(i);

        for (int j=0; j<m_ir->GetNPoints(); j++)
        {
            const IntegrationPoint &ip = m_ir->IntPoint(j);
            tr->SetIntPoint(&ip);

            u1->GetVectorValue(i, ip, uSol1);
            u2->GetVectorValue(i, ip, uSol2);
            for (int k=0; k<dim; k++)
                errVal(k) = uSol1(k) - uSol2(k);
            errDg += ip.weight*tr->Weight()*(errVal*errVal);
        }
    }

    return errDg;
}


// Assembles DG jumps for scalar FE space
// Time-like face
double mymfem::AssembleDgJumps
:: FtimeScalar (GridFunction* u)
{
    double tJump = 0;

    auto fes = u->FESpace();
    auto mesh = fes->GetMesh();

    FaceElementTransformations *ftr = nullptr;
    IntegrationPoint eip1, eip2;

    for (int i=0; i<mesh->GetNumFaces(); i++)
    {
        ftr = mesh -> GetInteriorFaceTransformations (i);
        if (ftr != nullptr)
        {
            int i1 = ftr->Elem1No;
            int i2 = ftr->Elem2No;
            for (int j=0; j<m_ir->GetNPoints(); j++)
            {
                const IntegrationPoint &ip
                        = m_ir->IntPoint(j);
                ftr->Face->SetIntPoint(&ip);

                ftr->Loc1.Transform(m_ir->IntPoint(j), eip1);
                ftr->Loc2.Transform(m_ir->IntPoint(j), eip2);

                double val = u->GetValue(i1, eip1)
                        - u->GetValue(i2, eip2);

                double coeff = ip.weight;
                if (m_stabParamsType == 1
                        || m_stabParamsType == 3) {
                    coeff *= ftr->Face->Weight();
                }
                tJump += coeff*val*val;
            }
        }
    }

    // Dirichlet boundary faces
    int nbfaces = mesh->GetNBE();
    for (int i = 0; i < nbfaces; i++)
    {
        const int bdr_attr = mesh->GetBdrAttribute(i);
        if (m_xEss_bdr_marker[bdr_attr-1] == 0) { continue; }

        ftr = mesh -> GetBdrFaceTransformations(i);

        int i1 = ftr->Elem1No;
        for (int j=0; j<m_ir->GetNPoints(); j++)
        {
            const IntegrationPoint &ip
                    = m_ir->IntPoint(j);
            ftr->Face->SetIntPoint(&ip);

            ftr->Loc1.Transform(m_ir->IntPoint(j), eip1);
            double val = u->GetValue(i1, eip1);

            double coeff = ip.weight;
            if (m_stabParamsType == 1
                    || m_stabParamsType == 3) {
                coeff *= ftr->Face->Weight();
            }
            tJump += coeff*val*val;
        }
    }

    return tJump;
}

// Assembles DG error for vector FE space
// Time-like face
double mymfem::AssembleDgJumps
:: FtimeVector (GridFunction* u)
{
    double tJump = 0;

    auto fes = u->FESpace();
    auto mesh = fes->GetMesh();
    int dim = mesh->Dimension();

    Vector uSol1, uSol2;

    FaceElementTransformations *ftr = nullptr;
    IntegrationPoint eip1, eip2;
    Vector nor(dim);

    for (int i=0; i<mesh->GetNumFaces(); i++)
    {
        ftr = mesh -> GetInteriorFaceTransformations (i);
        if (ftr != nullptr)
        {
            int i1 = ftr->Elem1No;
            int i2 = ftr->Elem2No;

            Vector buf(dim);
            for (int j=0; j<m_ir->GetNPoints(); j++)
            {
                const IntegrationPoint &ip
                        = m_ir->IntPoint(j);
                ftr->Face->SetIntPoint(&ip);

                CalcOrtho(ftr->Face->Jacobian(), nor);
                nor /= ftr->Face->Weight(); // unit normal

                ftr->Loc1.Transform(m_ir->IntPoint(j), eip1);
                ftr->Loc2.Transform(m_ir->IntPoint(j), eip2);

                u->GetVectorValue(i1, eip1, uSol1);
                u->GetVectorValue(i2, eip2, uSol2);
                for (int k=0; k<dim; k++)
                    buf(k) = uSol1(k) - uSol2(k);

                double coeff = ip.weight*ftr->Face->Weight();
                if (m_stabParamsType == 3
                        || m_stabParamsType == 4) {
                    coeff *= ftr->Face->Weight();
                }
                tJump += coeff*(buf*nor)*(buf*nor);
            }
        }
    }

    // Neumann boundary faces
    int nbfaces = mesh->GetNBE();
    for (int i = 0; i < nbfaces; i++)
    {
        const int bdr_attr = mesh->GetBdrAttribute(i);
        if (m_xNat_bdr_marker[bdr_attr-1] == 0) { continue; }

        ftr = mesh -> GetBdrFaceTransformations(i);
        int i1 = ftr->Elem1No;

        Vector buf(dim);
        for (int j=0; j<m_ir->GetNPoints(); j++)
        {
            const IntegrationPoint &ip
                    = m_ir->IntPoint(j);
            ftr->Face->SetIntPoint(&ip);

            CalcOrtho(ftr->Face->Jacobian(), nor);
            nor /= ftr->Face->Weight(); // unit normal

            ftr->Loc1.Transform(m_ir->IntPoint(j), eip1);

            u->GetVectorValue(i1, eip1, uSol1);

            double coeff = ip.weight*ftr->Face->Weight();
            if (m_stabParamsType == 3
                    || m_stabParamsType == 4) {
                coeff *= ftr->Face->Weight();
            }
            tJump += coeff*(uSol1*nor)*(uSol1*nor);
        }
    }

    return tJump;
}

// Assembles DG jumps for scalar FE space
// Space-like faces
double mymfem::AssembleDgJumps
:: FspaceScalar(Mesh *mesh,
                Coefficient* u,
                Coefficient* medium)
{
    double xJump = 0;

    double val;
    ElementTransformation *tr = nullptr;

    for (int i=0; i<mesh->GetNE(); i++)
    {
        tr = mesh->GetElementTransformation(i);
        for (int j=0; j<m_ir->GetNPoints(); j++)
        {
            const IntegrationPoint &ip = m_ir->IntPoint(j);
            tr->SetIntPoint(&ip);

            val = u->Eval(*tr, ip);
            xJump += ip.weight*tr->Weight()
                    *medium->Eval(*tr, ip)*(val*val);
        }
    }

    return xJump;
}

double mymfem::AssembleDgJumps
:: FspaceScalar(GridFunction* u, Coefficient* medium)
{
    double xJump = 0;

    auto fes = u->FESpace();
    auto mesh = fes->GetMesh();

    double val;
    ElementTransformation *tr = nullptr;

    for (int i=0; i<mesh->GetNE(); i++)
    {
        tr = mesh->GetElementTransformation(i);
        for (int j=0; j<m_ir->GetNPoints(); j++)
        {
            const IntegrationPoint &ip = m_ir->IntPoint(j);
            tr->SetIntPoint(&ip);

            val = u->GetValue(i, ip);
            xJump += ip.weight*tr->Weight()
                    *medium->Eval(*tr, ip)
                    *(val*val);
        }
    }

    return xJump;
}

double mymfem::AssembleDgJumps
:: FspaceScalar (GridFunction* u1, GridFunction* u2,
                 Coefficient* medium)
{
    double xJump = 0;

    auto fes = u1->FESpace();
    auto mesh = fes->GetMesh();

    double val;
    ElementTransformation *tr = nullptr;

    for (int i=0; i<mesh->GetNE(); i++)
    {
        tr = mesh->GetElementTransformation(i);
        for (int j=0; j<m_ir->GetNPoints(); j++)
        {
            const IntegrationPoint &ip = m_ir->IntPoint(j);
            tr->SetIntPoint(&ip);

            val = u1->GetValue(i, ip) - u2->GetValue(i, ip);
            xJump += ip.weight*tr->Weight()
                    *medium->Eval(*tr, ip)
                    *(val*val);
        }
    }

    return xJump;
}

// Assembles DG jumps for vector FE space
// Time-like face
double mymfem::AssembleDgJumps
:: FspaceVector(Mesh* mesh, VectorCoefficient* u)
{
    double xJump = 0;
    int dim = mesh->Dimension();

    Vector uSol(dim);
    ElementTransformation *tr = nullptr;

    for (int i=0; i<mesh->GetNE(); i++)
    {
        tr = mesh->GetElementTransformation(i);
        for (int j=0; j<m_ir->GetNPoints(); j++)
        {
            const IntegrationPoint &ip = m_ir->IntPoint(j);
            tr->SetIntPoint(&ip);

            u->Eval(uSol, *tr, ip);
            xJump += ip.weight*tr->Weight()*(uSol*uSol);
        }
    }

    return xJump;
}

double mymfem::AssembleDgJumps
:: FspaceVector(GridFunction* u)
{
    double xJump = 0;

    auto fes = u->FESpace();
    auto mesh = fes->GetMesh();
    int dim = mesh->Dimension();

    Vector uSol(dim);
    ElementTransformation *tr = nullptr;

    for (int i=0; i<mesh->GetNE(); i++)
    {
        tr = mesh->GetElementTransformation(i);
        for (int j=0; j<m_ir->GetNPoints(); j++)
        {
            const IntegrationPoint &ip = m_ir->IntPoint(j);
            tr->SetIntPoint(&ip);

            u->GetVectorValue(i, ip, uSol);
            xJump += ip.weight*tr->Weight()*(uSol*uSol);
        }
    }

    return xJump;
}

double mymfem::AssembleDgJumps
:: FspaceVector (GridFunction* u1, GridFunction* u2)
{
    double xJump = 0;

    auto fes = u1->FESpace();
    auto mesh = fes->GetMesh();
    int dim = mesh->Dimension();

    Vector uSol1(dim), uSol2(dim);
    Vector val(dim);
    ElementTransformation *tr = nullptr;

    for (int i=0; i<mesh->GetNE(); i++)
    {
        tr = mesh->GetElementTransformation(i);
        for (int j=0; j<m_ir->GetNPoints(); j++)
        {
            const IntegrationPoint &ip = m_ir->IntPoint(j);
            tr->SetIntPoint(&ip);

            u1->GetVectorValue(i, ip, uSol1);
            u2->GetVectorValue(i, ip, uSol2);
            for (int k=0; k<dim; k++)
                val(k) = uSol1(k) - uSol2(k);
            xJump += ip.weight*tr->Weight()*(val*val);
        }
    }

    return xJump;
}

// Assembles DG^{+} smooth error for scalar FE space
// Time-like face
double mymfem::AssembleDgPlusSmoothError
:: Ftime (GridFunction* u, Coefficient* uE)
{
    double errDgp = 0;

    auto fes = u->FESpace();
    auto mesh = fes->GetMesh();

    Vector errVal;
    ElementTransformation *tr = nullptr;
    FaceElementTransformations *ftr = nullptr;
    IntegrationPoint eip;

    errVal.SetSize(m_ir->GetNPoints());

    // interior faces
    for (int i=0; i<mesh->GetNumFaces(); i++)
    {
        ftr = mesh -> GetInteriorFaceTransformations (i);
        if (ftr != nullptr)
        {
            // face element 1
            int i1 = ftr->Elem1No;
            tr = ftr->Elem1;
            for (int j=0; j<m_ir->GetNPoints(); j++)
            {
                ftr->Loc1.Transform(m_ir->IntPoint(j), eip);
                errVal(j) = 0.5*(uE->Eval(*tr, eip)
                        - u->GetValue(i1, eip));
            }

            // face element 2
            int i2 = ftr->Elem2No;
            tr = ftr->Elem2;
            for (int j=0; j<m_ir->GetNPoints(); j++)
            {
                ftr->Loc2.Transform(m_ir->IntPoint(j), eip);
                errVal(j) += 0.5*(uE->Eval(*tr, eip)
                        - u->GetValue(i2, eip));
            }

            for (int j=0; j<m_ir->GetNPoints(); j++)
            {
                const IntegrationPoint &ip
                        = m_ir->IntPoint(j);
                ftr->Face->SetIntPoint(&ip);

                double coeff = ip.weight;
                if (m_stabParamsType == 1
                        || m_stabParamsType == 2) {
                    coeff *= ftr->Face->Weight();
                }
                errDgp += coeff*errVal(j)*errVal(j);
            }
        }
    }

    // Neumann boundary faces
    int nbfaces = mesh->GetNBE();
    for (int i = 0; i < nbfaces; i++)
    {
        const int bdr_attr = mesh->GetBdrAttribute(i);
        if (m_xNat_bdr_marker[bdr_attr-1] == 0) { continue; }

        ftr = mesh -> GetBdrFaceTransformations(i);

        // face element 1
        int i1 = ftr->Elem1No;
        tr = ftr->Elem1;
        for (int j=0; j<m_ir->GetNPoints(); j++)
        {
            ftr->Loc1.Transform(m_ir->IntPoint(j), eip);
            errVal(j) = uE->Eval(*tr, eip)
                    - u->GetValue(i1, eip);
        }

        for (int j=0; j<m_ir->GetNPoints(); j++)
        {
            const IntegrationPoint &ip
                    = m_ir->IntPoint(j);
            ftr->Face->SetIntPoint(&ip);

            double coeff = ip.weight;
            if (m_stabParamsType == 1
                    || m_stabParamsType == 2) {
                coeff *= ftr->Face->Weight();
            }
            errDgp += coeff*errVal(j)*errVal(j);
        }
    }

    return errDgp;
}


// Assembles DG^{+} smooth error for vector FE space
// Time-like face
double mymfem::AssembleDgPlusSmoothError
:: Ftime (GridFunction* u, VectorCoefficient* uE)
{
    double errDgp = 0;

    auto fes = u->FESpace();
    auto mesh = fes->GetMesh();
    int dim = mesh->Dimension();

    Vector uSol, uESol;
    DenseMatrix errVal;
    Vector nor(dim);

    ElementTransformation *tr = nullptr;
    FaceElementTransformations *ftr = nullptr;
    IntegrationPoint eip;

    errVal.SetSize(dim, m_ir->GetNPoints());

    // interior faces
    for (int i=0; i<mesh->GetNumFaces(); i++)
    {
        ftr = mesh -> GetInteriorFaceTransformations (i);
        if (ftr != nullptr)
        {
            // face element 1
            int i1 = ftr->Elem1No;
            tr = ftr->Elem1;
            for (int j=0; j<m_ir->GetNPoints(); j++)
            {
                ftr->Loc1.Transform(m_ir->IntPoint(j), eip);
                uE->Eval(uESol, *tr, eip);
                u->GetVectorValue(i1, eip, uSol);
                for (int k=0; k<dim; k++)
                    errVal(k,j) = 0.5*(uESol(k) - uSol(k));
            }

            // face element 2
            int i2 = ftr->Elem2No;
            tr = ftr->Elem2;
            for (int j=0; j<m_ir->GetNPoints(); j++)
            {
                ftr->Loc2.Transform(m_ir->IntPoint(j), eip);
                uE->Eval(uESol, *tr, eip);
                u->GetVectorValue(i2, eip, uSol);
                for (int k=0; k<dim; k++)
                    errVal(k,j) += 0.5*(uESol(k) - uSol(k));
            }

            Vector buf;
            for (int j=0; j<m_ir->GetNPoints(); j++)
            {
                const IntegrationPoint &ip
                        = m_ir->IntPoint(j);
                ftr->Face->SetIntPoint(&ip);

                CalcOrtho(ftr->Face->Jacobian(), nor);
                nor /= ftr->Face->Weight(); // unit normal

                double coeff = ip.weight*ftr->Face->Weight();
                if (m_stabParamsType == 2
                        || m_stabParamsType == 4) {
                    coeff *= ftr->Face->Weight();
                }
                errVal.GetColumnReference(j, buf);
                errDgp += coeff*(buf*nor)*(buf*nor);
            }
        }
    }

    // Dirichlet boundary faces
    int nbfaces = mesh->GetNBE();
    for (int i = 0; i < nbfaces; i++)
    {
        const int bdr_attr = mesh->GetBdrAttribute(i);
        if (m_xEss_bdr_marker[bdr_attr-1] == 0) { continue; }

        ftr = mesh -> GetBdrFaceTransformations(i);

        // face element 1
        int i1 = ftr->Elem1No;
        tr = ftr->Elem1;
        for (int j=0; j<m_ir->GetNPoints(); j++)
        {
            ftr->Loc1.Transform(m_ir->IntPoint(j), eip);
            uE->Eval(uESol, *tr, eip);
            u->GetVectorValue(i1, eip, uSol);
            for (int k=0; k<dim; k++)
                errVal(k,j) = uESol(k) - uSol(k);
        }

        Vector buf;
        for (int j=0; j<m_ir->GetNPoints(); j++)
        {
            const IntegrationPoint &ip
                    = m_ir->IntPoint(j);
            ftr->Face->SetIntPoint(&ip);

            CalcOrtho(ftr->Face->Jacobian(), nor);
            nor /= ftr->Face->Weight(); // unit normal

            double coeff = ip.weight*ftr->Face->Weight();
            if (m_stabParamsType == 2
                    || m_stabParamsType == 4) {
                coeff *= ftr->Face->Weight();
            }
            errVal.GetColumnReference(j, buf);
            errDgp += coeff*(buf*nor)*(buf*nor);
        }
    }

    return errDgp;
}

// Assembles DG^{+} smooth error for scalar FE space
// space-like face
double mymfem::AssembleDgPlusSmoothError
:: Fspace (GridFunction* u, Coefficient* uE)
{
    double errDgp = 0;

    auto fes = u->FESpace();
    auto mesh = fes->GetMesh();

    double errVal;
    ElementTransformation *tr = nullptr;

    for (int i=0; i<mesh->GetNE(); i++)
    {
        tr = mesh->GetElementTransformation(i);

        for (int j=0; j<m_ir->GetNPoints(); j++)
        {
            const IntegrationPoint &ip = m_ir->IntPoint(j);
            tr->SetIntPoint(&ip);

            errVal = uE->Eval(*tr, ip) - u->GetValue(i, ip);
            errDgp += ip.weight*tr->Weight()*(errVal*errVal);
        }
    }

    return errDgp;
}

double mymfem::AssembleDgPlusSmoothError
:: Fspace (GridFunction* u, VectorCoefficient* uE)
{
    double errDgp = 0;

    auto fes = u->FESpace();
    auto mesh = fes->GetMesh();
    int dim = mesh->Dimension();

    Vector uSol(dim), uESol(dim);
    Vector errVal(dim);
    ElementTransformation *tr = nullptr;

    for (int i=0; i<mesh->GetNE(); i++)
    {
        tr = mesh->GetElementTransformation(i);

        for (int j=0; j<m_ir->GetNPoints(); j++)
        {
            const IntegrationPoint &ip = m_ir->IntPoint(j);
            tr->SetIntPoint(&ip);

            uE->Eval(uESol, *tr, ip);
            u->GetVectorValue(i, ip, uSol);
            for (int k=0; k<dim; k++)
                errVal(k) = uESol(k) - uSol(k);
            errDgp += ip.weight*tr->Weight()*(errVal*errVal);
        }
    }

    return errDgp;
}

// End of file
