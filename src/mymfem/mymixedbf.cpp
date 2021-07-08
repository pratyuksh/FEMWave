#include "../../include/mymfem/mymixedbf.hpp"


void mymfem::MyMixedBilinearForm :: MyAssemble (int skip_zeros)
{
    // allocate memory
    if (mat == nullptr) {
        mat = new SparseMatrix(height, width);
    }

    // assemble interior face integrators
    if (myifbfi.Size())
    {
        DenseMatrix elmat;
        FaceElementTransformations *ftr = nullptr;
        Array<int> trial_vdofs, test_vdofs;
        Array<int> trial_vdofs2, test_vdofs2;
        Mesh *mesh = trial_fes -> GetMesh();

        int nfaces = mesh->GetNumFaces();
        for (int i = 0; i < nfaces; i++)
        {
            ftr = mesh -> GetInteriorFaceTransformations (i);
            if (ftr != nullptr)
            {
                const FiniteElement &tr_fe1
                        = *(trial_fes->GetFE (ftr->Elem1No));
                const FiniteElement &tr_fe2
                        = *(trial_fes->GetFE (ftr->Elem2No));
                
                const FiniteElement &te_fe1
                        = *(test_fes->GetFE (ftr->Elem1No));
                const FiniteElement &te_fe2
                        = *(test_fes->GetFE (ftr->Elem2No));

                trial_fes -> GetElementVDofs (ftr -> Elem1No,
                                              trial_vdofs);
                trial_fes -> GetElementVDofs (ftr -> Elem2No,
                                              trial_vdofs2);
                trial_vdofs.Append (trial_vdofs2);
                
                test_fes -> GetElementVDofs (ftr -> Elem1No,
                                             test_vdofs);
                test_fes -> GetElementVDofs (ftr -> Elem2No,
                                             test_vdofs2);
                test_vdofs.Append (test_vdofs2);

                for (int k = 0; k < myifbfi.Size(); k++)
                {
                    myifbfi[k] -> MyAssembleFaceMatrix (tr_fe1, tr_fe2,
                                                        te_fe1, te_fe2,
                                                        *ftr, elmat);
                    //elmat.Print();
                    mat -> AddSubMatrix (test_vdofs, trial_vdofs,
                                         elmat, skip_zeros);
                }
            }
        }
    }

    // assemble boundary face integrators
    if (mybfbfi.Size())
    {
        DenseMatrix elmat;
        FaceElementTransformations *ftr = nullptr;
        Array<int> trial_vdofs, test_vdofs;
        Mesh *mesh = trial_fes -> GetMesh();

        // set boundary markers
        Array<int> bdr_attr_marker(mesh->bdr_attributes.Size() ?
                                   mesh->bdr_attributes.Max() : 0);
        bdr_attr_marker = 0;
        for (int k = 0; k < mybfbfi.Size(); k++)
        {
            if (mybfmarker[k] == nullptr)
            {
                bdr_attr_marker = 1;
                break;
            }
            Array<int> &bdr_marker = *mybfmarker[k];
            MFEM_ASSERT(bdr_marker.Size() == bdr_attr_marker.Size(),
                        "invalid boundary marker for boundary face "
                        "integrator #" << k << ", counting from zero");
            for (int i = 0; i < bdr_attr_marker.Size(); i++)
            {
                bdr_attr_marker[i] |= bdr_marker[i];
            }
        }

        int nbfaces = mesh->GetNBE();
        for (int i = 0; i < nbfaces; i++)
        {
            const int bdr_attr = mesh->GetBdrAttribute(i);
            if (bdr_attr_marker[bdr_attr-1] == 0) { continue; }

            ftr = mesh -> GetBdrFaceTransformations(i);
            const FiniteElement &tr_fe
                    = *(trial_fes->GetFE (ftr->Elem1No));
            const FiniteElement &te_fe
                    = *(test_fes->GetFE (ftr->Elem1No));

            trial_fes -> GetElementVDofs (ftr -> Elem1No,
                                          trial_vdofs);
            test_fes -> GetElementVDofs (ftr -> Elem1No,
                                         test_vdofs);

            for (int k = 0; k < mybfbfi.Size(); k++)
            {
                mybfbfi[k] -> MyAssembleFaceMatrix (tr_fe, tr_fe,
                                                    te_fe, te_fe,
                                                    *ftr, elmat);

                mat -> AddSubMatrix (test_vdofs, trial_vdofs,
                                     elmat, skip_zeros);
            }
        }
    }
}


// End of file
