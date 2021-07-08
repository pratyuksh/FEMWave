#ifndef MYMFEM_MIXEDBF_HPP
#define MYMFEM_MIXEDBF_HPP

#include "mfem.hpp"
using namespace mfem;

#include "mybfinteg.hpp"


namespace mymfem {

class MyMixedBilinearForm : public MixedBilinearForm
{
public:
    MyMixedBilinearForm(FiniteElementSpace *tr_fes,
                        FiniteElementSpace *te_fes)
        : MixedBilinearForm(tr_fes, te_fes) { }

    void MyAddInteriorFaceIntegrator (MyBilinearFormIntegrator * bfi) {
        myifbfi.Append (bfi);
    }

    void MyAddBoundaryFaceIntegrator (MyBilinearFormIntegrator * bfi) {
        mybfbfi.Append (bfi);
        mybfmarker.Append(nullptr);
    }

    void MyAddBoundaryFaceIntegrator (MyBilinearFormIntegrator * bfi,
                                      Array<int> &bdr_attr_marker) {
        mybfbfi.Append (bfi);
        mybfmarker.Append(&bdr_attr_marker);
    }

    void MyAssemble(int skip_zeros=1);

    void MyClear()
    {
        for (int k=0; k < myifbfi.Size(); k++) { delete myifbfi[k]; }
        for (int k=0; k < mybfbfi.Size(); k++) { delete mybfbfi[k]; }
    }

private:
    Array<MyBilinearFormIntegrator *> myifbfi;
    Array<MyBilinearFormIntegrator *> mybfbfi;
    Array<Array<int>*> mybfmarker;
};

}

#endif /// MYMFEM_MIXEDBF_HPP
