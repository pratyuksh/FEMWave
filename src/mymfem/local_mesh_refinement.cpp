#include "../../include/mymfem/local_mesh_refinement.hpp"
#include "../../include/mymfem/utilities.hpp"


//! Evaluates R0 distance
//! for all corners of the polygon
Array<double> LocalMeshRefinement
:: eval_R0()
{
    Array<double> R0;

    auto cornersCoords = m_polygon->get_corners_coords();
    int nc = cornersCoords.NumCols();
    R0.SetSize(nc);

    Vector vi(cornersCoords.NumRows());
    Vector vj(cornersCoords.NumRows());
    for (int i=0; i<nc; i++) {
        double minDist = 1E+5;
        cornersCoords.GetColumn(i,vi);
        for (int j=0; j<nc; j++) {
            if (i != j) {
                cornersCoords.GetColumn(j,vj);
                double val = vi.DistanceTo(vj);
                minDist = std::min(val, minDist);
            }
        }
        R0[i] = m_factor*minDist;
    }

    return R0;
}

//! uniform mesh refinement
void LocalMeshRefinement
:: uniform (std::shared_ptr<Mesh>& mesh, double h0)
{
    double hMin, hMax, kappaMin, kappaMax;
    mesh->GetCharacteristics(hMin, hMax,
                             kappaMin, kappaMax);

    while (hMax > h0) {
        std::cout << hMax << "\t" << h0 << std::endl;
        mesh->UniformRefinement();
        mesh->GetCharacteristics(hMin, hMax,
                                 kappaMin, kappaMax);
    }
}

//! local mesh refinement
void LocalMeshRefinement
:: local (std::shared_ptr<Mesh>& mesh,
          int deg, double h0, bool nonConforming)
{
    auto cornerCoords = m_polygon->get_corners_coords();
    auto sCornersIds = m_polygon->get_singular_corners_ids();
    auto refineFlags = m_polygon->get_refine_flags();
    auto refineWeights = m_polygon->get_refine_weights();

    int nsc = sCornersIds.Size();

    Vector corner;
    for (int i=0; i<nsc; i++)
    {
        double zeta = refineWeights[i];
        if (refineFlags[i])
        {
            int ii = sCornersIds[i]; // corner id
            cornerCoords.GetColumn(ii, corner);

            double gamma = 1 - zeta;
            auto R0 = m_distR0[ii];
            auto K = std::ceil
                    (-std::log2(h0)*(deg+1)/gamma -1);
            int NLocRef = int(2*K+1);
            std::cout << h0 << "\t" << gamma << "\t"
                      << K << "\t" << NLocRef << "\t"
                      << R0 << std::endl;

            Array<int> markedEls;
            for (int j=0; j<NLocRef; j++)
            {
                double weight = 1 - gamma;
                double expo = -j*(deg+weight)/(2.*(deg+1));
                double hMin = h0*std::pow(2, expo);
                double RMax = R0*std::pow(2, -j/2.);

                int ne = mesh->GetNE();
                for (int k=0; k<ne; k++)
                {
                    double hEl = mesh->GetElementSize(k);
                    double dist
                            = get_element_distance_to_point
                            (mesh, k, corner);
                    if (dist < RMax && hEl > hMin) {
                        markedEls.Append(k);
                    }
                }
                mesh->GeneralRefinement(markedEls,
                                        nonConforming);
                markedEls.DeleteAll();
            }
        }
    }
}

// End of file
