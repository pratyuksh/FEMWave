#include "../../include/waveO2/test_cases.hpp"
#include <fmt/format.h>


// Case: UnitSquare_Test1
// Smooth solution
// Homogeneous Dirichlet boundary conditions
// Homogeneous source
double WaveO2TestCase <UnitSquare_Test1>
:: exact_sol(const Vector& x, const double t) const
{
    double uSol
            = sin(M_PI*x(0))*sin(M_PI*x(1))
            *sin(std::sqrt(2)*M_PI*m_c*t);
    return uSol;
}

double WaveO2TestCase <UnitSquare_Test1>
:: exact_gradt_sol(const Vector& x, const double t) const
{
    double dudtSol
            = sin(M_PI*x(0))*sin(M_PI*x(1))
            *std::sqrt(2)*M_PI*m_c
            *cos(std::sqrt(2)*M_PI*m_c*t);
    return dudtSol;
}

Vector WaveO2TestCase <UnitSquare_Test1>
:: exact_gradx_sol (const Vector& x, const double t) const
{
    Vector dudxSol(m_dim);
    dudxSol(0) = cos(M_PI*x(0))*sin(M_PI*x(1));
    dudxSol(1) = sin(M_PI*x(0))*cos(M_PI*x(1));
    dudxSol *= M_PI*sin(std::sqrt(2)*M_PI*m_c*t);
    return dudxSol;
}

double WaveO2TestCase <UnitSquare_Test1>
:: source (const Vector&, const double) const
{
    double f = 0;
    return f;
}

void WaveO2TestCase <UnitSquare_Test1>
:: set_bdry_dirichlet(Array<int>& bdr_marker) const
{
    bdr_marker = 1;
}

// End of file
