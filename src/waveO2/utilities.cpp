#include "../../include/waveO2/utilities.hpp"


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

// End of file
