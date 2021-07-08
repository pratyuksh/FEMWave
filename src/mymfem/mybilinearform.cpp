#include "../../include/mymfem/mybilinearform.hpp"
#include "../../include/mymfem/utilities.hpp"


void mymfem::MyBilinearForm :: MyAssemble (int skip_zeros)
{
    // allocate memory
    if (mat == nullptr) {
        AllocMat();
    }

    // assemble domain integrators
    if (mydbfi.Size())
    {
        DenseMatrix elmat;
        ElementTransformation *eltrans = nullptr;

        for (int i = 0; i < fes -> GetNE(); i++)
        {
            fes->GetElementVDofs(i, vdofs);

            const FiniteElement &fe = *fes->GetFE(i);
            eltrans = fes->GetElementTransformation(i);
            for (int k = 0; k < mydbfi.Size(); k++)
            {
                mydbfi[k]->AssembleElementMatrix(fe, *eltrans, elmat);
                mat->AddSubMatrix(vdofs, vdofs, elmat, skip_zeros);
            }
        }
    }

    // assemble interior face integrators
    if (myifbfi.Size())
    {
        DenseMatrix elmat;
        FaceElementTransformations *ftr = nullptr;
        Array<int> vdofs2;
        Mesh *mesh = fes -> GetMesh();

        int nfaces = mesh->GetNumFaces();
        for (int i = 0; i < nfaces; i++)
        {
            ftr = mesh -> GetInteriorFaceTransformations (i);
            if (ftr != nullptr)
            {
                const FiniteElement &fe1 = *(fes->GetFE (ftr->Elem1No));
                const FiniteElement &fe2 = *(fes->GetFE (ftr->Elem2No));

                fes -> GetElementVDofs (ftr -> Elem1No, vdofs);
                fes -> GetElementVDofs (ftr -> Elem2No, vdofs2);
                vdofs.Append (vdofs2);
                for (int k = 0; k < myifbfi.Size(); k++)
                {
                    myifbfi[k] -> AssembleFaceMatrix
                            (fe1, fe2, *ftr, elmat);

                    mat -> AddSubMatrix (vdofs, vdofs, elmat, skip_zeros);
                }
            }
        }
    }

    // assemble boundary face integrators
    if (mybfbfi.Size())
    {
        DenseMatrix elmat;
        FaceElementTransformations *ftr = nullptr;
        Mesh *mesh = fes -> GetMesh();

        int nbfaces = mesh->GetNBE();
        for (int i = 0; i < nbfaces; i++)
        {
            ftr = mesh -> GetBdrFaceTransformations(i);
            const FiniteElement &fe = *(fes->GetFE (ftr->Elem1No));

            fes -> GetElementVDofs (ftr -> Elem1No, vdofs);
            for (int k = 0; k < mybfbfi.Size(); k++)
            {
                mybfbfi[k] -> AssembleFaceMatrix
                        (fe, fe, *ftr, elmat);
                mat -> AddSubMatrix (vdofs, vdofs, elmat, skip_zeros);
            }
        }
    }
}


// End of file
