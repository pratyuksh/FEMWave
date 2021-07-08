#ifndef MYMFEM_BFINTEG_HPP
#define MYMFEM_BFINTEG_HPP

#include "mfem.hpp"
using namespace mfem;


namespace mymfem {

/// Abstract base class MyBilinearFormIntegrator
class MyBilinearFormIntegrator
        : public BilinearFormIntegrator
{
public:
    virtual void MyAssembleFaceMatrix
    (const FiniteElement &trial_fe1,
     const FiniteElement &trial_fe2,
     const FiniteElement &test_fe1,
     const FiniteElement &test_fe2,
     FaceElementTransformations &ftr, DenseMatrix &elmat)
    const
    {
        MFEM_ABORT("MyAssembleFaceMatrix (mixed form) "
                   "is not implemented "
                   "for this Integrator class.");
    }
    
    virtual void MyAssembleFaceMatrix
    (const FiniteElement &trial_fe,
     const FiniteElement &test_fe,
     FaceElementTransformations &ftr, DenseMatrix &elmat)
    const
    {
        MFEM_ABORT("MyAssembleFaceMatrix (mixed form) "
                   "is not implemented "
                   "for this Integrator class.");
    }

    virtual ~MyBilinearFormIntegrator() {}
};

}

#endif /// MYMFEM_BFINTEG_HPP
