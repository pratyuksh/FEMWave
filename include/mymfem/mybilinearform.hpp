#ifndef MYMFEM_BILINEARFORM_HPP
#define MYMFEM_BILINEARFORM_HPP

#include "mfem.hpp"
using namespace mfem;

#include "../includes.hpp"


namespace mymfem {

class MyBilinearForm : public BilinearForm
{
public:
    MyBilinearForm(FiniteElementSpace *fes) : BilinearForm(fes) { }
    
    void MyAddDomainIntegrator (BilinearFormIntegrator * bfi) {
        mydbfi.Append (bfi);
    }

    void MyAddInteriorFaceIntegrator (BilinearFormIntegrator * bfi) {
        myifbfi.Append (bfi);
    }

    void MyAddBoundaryFaceIntegrator (BilinearFormIntegrator * bfi) {
        mybfbfi.Append (bfi);
        mybfmarker.Append(nullptr);
    }

    void MyAddBoundaryFaceIntegrator (BilinearFormIntegrator * bfi,
                                      Array<int> &bdr_attr_marker) {
        mybfbfi.Append (bfi);
        mybfmarker.Append(&bdr_attr_marker);
    }
    
    void MyAddFaceIntegrator (BilinearFormIntegrator * bfi) {
        MyAddInteriorFaceIntegrator(bfi);
        MyAddBoundaryFaceIntegrator(bfi);
    }
    
    void MyAssemble(int skip_zeros=1);
    
    void MyClear()
    {
        for (int k=0; k < mydbfi.Size(); k++) { delete mydbfi[k]; }
        for (int k=0; k < myifbfi.Size(); k++) { delete myifbfi[k]; }
        for (int k=0; k < mybfbfi.Size(); k++) { delete mybfbfi[k]; }
        for (int k=0; k < myfbfi.Size(); k++) { delete myfbfi[k]; }
    }

private:
    Array<BilinearFormIntegrator *> mydbfi;
    Array<BilinearFormIntegrator *> myifbfi;
    Array<BilinearFormIntegrator *> mybfbfi;
    Array<BilinearFormIntegrator *> myfbfi;
    Array<Array<int>*> mybfmarker;
};

}

#endif /// MYMFEM_BILINEARFORM_HPP
