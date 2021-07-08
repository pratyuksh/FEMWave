#ifndef WAVEO2_UTILITIES_HPP
#define WAVEO2_UTILITIES_HPP

#include "mfem.hpp"
using namespace mfem;


void build_xSol_FG(const Vector& uXtSol,
                   const Vector& tShape,
                   Array<int> tVdofs,
                   Vector& uSol);

#endif /// WAVEO2_UTILITIES_HPP
