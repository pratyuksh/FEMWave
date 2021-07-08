#ifndef MYMFEM_GEOMETRY_HPP
#define MYMFEM_GEOMETRY_HPP

#include "mfem.hpp"
using namespace mfem;


class Polygon
{
public:
    //! Checks if given vertex is a singular corner
    std::pair<bool, int> is_singular_corner(Vector& v)
    {
        int nsc = m_singularCornersIds.Size();
        if (nsc == 0) { return {false, -1}; };

        Vector coords(v.Size());
        for (int i=0; i<nsc; i++)
        {
            int ii = m_singularCornersIds[i];
            m_cornersCoords.GetColumn(ii, coords);
            if (coords.DistanceTo(v) <= m_tol) {
                return {true, i};
            }
        }

        return {false, -2};
    }

    //! Sets corners coords
    void set_corners_coords
    (DenseMatrix& cornersCoords) {
        m_cornersCoords = cornersCoords;
    }

    //! Returns corners coords
    DenseMatrix get_corners_coords() const {
        return m_cornersCoords;
    }

    //! Returns number of corners
    int get_num_corners() const {
        return m_cornersCoords.Size();
    }

    //! Sets singular corners indices
    void set_singular_corners_ids
    (Array<int>& singularCornerIds) {
        m_singularCornersIds = singularCornerIds;
    }

    //! Returns singular corners ids
    Array<int> get_singular_corners_ids() const {
        return m_singularCornersIds;
    }

    //! Sets refine flags
    void set_refine_flags
    (Array<bool>& refineFlags) {
        m_refineFlags = refineFlags;
    }

    //! Returns refine flags
    Array<bool> get_refine_flags() const {
        return m_refineFlags;
    }

    //! Sets refine weights
    void set_refine_weights
    (Array<double>& refineWeights) {
        m_refineWeights = refineWeights;
    }

    //! Returns refine weights
    Array<double> get_refine_weights() const {
        return m_refineWeights;
    }

protected:
    DenseMatrix m_cornersCoords;
    Array<int> m_singularCornersIds;
    Array<bool> m_refineFlags;
    Array<double> m_refineWeights;
    double m_tol = 1E-15;
};


//! L-shaped domain
class LShaped : public Polygon
{
public:
    //! Constructor
    LShaped ()
    {
        double a = 0.5;
        double left = -a;
        double right = a;
        double bottom = -a;
        double top = a;

        double xMid = 0.5*(left+right);
        double yMid = 0.5*(bottom+top);

        // add corners
        m_cornersCoords.SetSize(2, 6);

        m_cornersCoords(0,0) = left;
        m_cornersCoords(1,0) = bottom;

        m_cornersCoords(0,1) = right;
        m_cornersCoords(1,1) = bottom;

        m_cornersCoords(0,2) = right;
        m_cornersCoords(1,2) = yMid;

        m_cornersCoords(0,3) = xMid;
        m_cornersCoords(1,3) = yMid;

        m_cornersCoords(0,4) = xMid;
        m_cornersCoords(1,4) = top;

        m_cornersCoords(0,5) = left;
        m_cornersCoords(1,5) = top;

        // default singular corners
        m_singularCornersIds.SetSize(1);
        m_singularCornersIds[0] = 3;

        // default refine flags and weights
        m_refineFlags.SetSize(1);
        m_refineWeights.SetSize(1);
        m_refineFlags[0] = false;
        m_refineWeights = 0.5;
    }
};


#endif /// MYMFEM_GEOMETRY_HPP
