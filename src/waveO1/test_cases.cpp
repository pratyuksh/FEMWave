#include "../../include/waveO1/test_cases.hpp"
#include <fmt/format.h>


// Case: UnitSquare_Test1
// Smooth solution
// Homogeneous Dirichlet boundary conditions
// Homogeneous source
double WaveO1TestCase <UnitSquare_Test1>
:: pressure_sol(const Vector& x, const double t) const
{
    double pressure
            = sin(M_PI*x(0))*sin(M_PI*x(1))
            *cos(std::sqrt(2)*M_PI*m_c*t);
    pressure *= std::sqrt(2)*M_PI*m_c;
    return pressure;
}

Vector WaveO1TestCase <UnitSquare_Test1>
:: velocity_sol (const Vector& x, const double t) const
{
    Vector v(m_dim);
    v(0) = cos(M_PI*x(0))*sin(M_PI*x(1));
    v(1) = sin(M_PI*x(0))*cos(M_PI*x(1));
    v *= -M_PI*sin(std::sqrt(2)*M_PI*m_c*t);
    return v;
}

double WaveO1TestCase <UnitSquare_Test1>
:: source (const Vector&, const double) const
{
    double f = 0;
    return f;
}

void WaveO1TestCase <UnitSquare_Test1>
:: set_bdry_dirichlet(Array<int>& bdr_marker) const
{
    bdr_marker = 1;
}


// Case: UnitSquare_Test2
// Smooth solution
// Homogeneous Dirichlet boundary conditions
// Non-homogeneous source
double WaveO1TestCase <UnitSquare_Test2>
:: pressure_sol(const Vector& x, const double t) const
{
    double pressure
            = sin(M_PI*x(0))*sin(M_PI*x(1))
            *cos(std::sqrt(3)*M_PI*m_c*t);
    pressure *= std::sqrt(3)*M_PI*m_c;
    return pressure;
}

Vector WaveO1TestCase <UnitSquare_Test2>
:: velocity_sol (const Vector& x, const double t) const
{
    Vector v(m_dim);
    v(0) = cos(M_PI*x(0))*sin(M_PI*x(1));
    v(1) = sin(M_PI*x(0))*cos(M_PI*x(1));
    v *= -M_PI*sin(std::sqrt(3)*M_PI*m_c*t);
    return v;
}

double WaveO1TestCase <UnitSquare_Test2>
:: source (const Vector& x, const double t) const
{
    double f = -M_PI*M_PI
            *sin(M_PI*x(0))*sin(M_PI*x(1))
            *sin(std::sqrt(3)*M_PI*m_c*t);
    return f;
}

void WaveO1TestCase <UnitSquare_Test2>
:: set_bdry_dirichlet(Array<int>& bdr_marker) const
{
    bdr_marker = 1;
}


// Case: UnitSquare_Test3
// Smooth solution
// Non-homogeneous Dirichlet boundary conditions
// Homogeneous source
double WaveO1TestCase <UnitSquare_Test3>
:: pressure_sol(const Vector& x, const double t) const
{
    double pressure
            = cos(M_PI*x(0))*cos(M_PI*x(1))
            *cos(std::sqrt(2)*M_PI*m_c*t);
    pressure *= std::sqrt(2)*M_PI*m_c;
    return pressure;
}

Vector WaveO1TestCase <UnitSquare_Test3>
:: velocity_sol (const Vector& x, const double t) const
{
    Vector v(m_dim);
    v(0) = -sin(M_PI*x(0))*cos(M_PI*x(1));
    v(1) = -cos(M_PI*x(0))*sin(M_PI*x(1));
    v *= -M_PI*sin(std::sqrt(2)*M_PI*m_c*t);
    return v;
}

double WaveO1TestCase <UnitSquare_Test3>
:: source (const Vector&, const double) const
{
    double f = 0;
    return f;
}

void WaveO1TestCase <UnitSquare_Test3>
:: set_bdry_dirichlet(Array<int>& bdr_marker) const
{
    bdr_marker = 1;
}


// Case: UnitSquare_Test4
// Smooth solution
// Homogeneous Dirichlet boundary conditions
// Non-homogeneous source
double WaveO1TestCase <UnitSquare_Test4>
:: pressure_sol(const Vector& x, const double t) const
{
    double pressure
            = sin(M_PI*x(0))*sin(M_PI*x(1));
    pressure *= M_PI*sin(2*M_PI*t);
    return pressure;
}

Vector WaveO1TestCase <UnitSquare_Test4>
:: velocity_sol (const Vector& x, const double t) const
{
    Vector v(m_dim);
    v(0) = cos(M_PI*x(0))*sin(M_PI*x(1));
    v(1) = sin(M_PI*x(0))*cos(M_PI*x(1));
    v *= -M_PI*sin(M_PI*t)*sin(M_PI*t);
    return v;
}

double WaveO1TestCase <UnitSquare_Test4>
:: source (const Vector& x, const double t) const
{
    double f = 2*M_PI*M_PI*(cos(2*M_PI*t)
                          + sin(M_PI*t)*sin(M_PI*t));
    f *= sin(M_PI*x(0))*sin(M_PI*x(1));
    return f;
}

void WaveO1TestCase <UnitSquare_Test4>
:: set_bdry_dirichlet(Array<int>& bdr_marker) const
{
    bdr_marker = 1;
}


// Case: GammaShaped_Test1
// Non-homogeneous Neumann boundary
// Non-homogeneous source
double WaveO1TestCase <GammaShaped_Test1>
:: radius(const Vector& x) const
{
    return x.Norml2();
}

double WaveO1TestCase <GammaShaped_Test1>
:: polar_angle(const Vector& x) const
{
    double theta = 2*M_PI*(x(1) < 0) + std::atan2(x(1), x(0));
    return theta;
}

double WaveO1TestCase <GammaShaped_Test1>
:: pressure_sol(const Vector& x, const double t) const
{
    double r = radius(x);
    double theta = polar_angle(x);

    double pressure
            = std::pow(r,m_gamma)*sin(m_gamma*theta)
            *cos(std::sqrt(2)*M_PI*t);
    pressure *= std::sqrt(2)*M_PI;
    return pressure;
}

Vector WaveO1TestCase <GammaShaped_Test1>
:: velocity_sol (const Vector& x, const double t) const
{
    double r = radius(x);
    double theta = polar_angle(x);

    double coeff = m_gamma-1;
    Vector v(m_dim);
    v(0) = sin(coeff*theta);
    v(1) = cos(coeff*theta);
    v *= -m_gamma*std::pow(r,coeff)*sin(std::sqrt(2)*M_PI*t);
    return v;
}

double WaveO1TestCase <GammaShaped_Test1>
:: source (const Vector& x, const double t) const
{
    double r = radius(x);
    double theta = polar_angle(x);

    double f = -2*M_PI*M_PI*std::pow(r,m_gamma)
            *sin(m_gamma*theta)*sin(std::sqrt(2)*M_PI*t);
    return f;
}

void WaveO1TestCase <GammaShaped_Test1>
:: set_bdry_dirichlet(Array<int>& bdr_marker) const
{
    bdr_marker = 0;
}


// Case: LShaped_Test1
// Homogeneous Dirichlet boundary
// Non-homogeneous source
double WaveO1TestCase <LShaped_Test1>
:: radius(const Vector& x) const
{
    return x.Norml2();
}

double WaveO1TestCase <LShaped_Test1>
:: polar_angle(const Vector& x) const
{
    double theta = std::atan2(-x[0],x[1])
                + (x[0] >= 0 && x[1] <= 0 ? 2*M_PI : 0);
    return theta;
}

/*double WaveO1TestCase <LShaped_Test1>
:: pressure_sol(const Vector& x, const double t) const
{
    double r = radius(x);
    double theta = polar_angle(x);

    double pressure
            = std::pow(r, m_gamma)*sin(m_gamma*theta)
            *(1./4. - x[0]*x[0])*(1./4. - x[1]*x[1])
            *cos(std::sqrt(2)*M_PI*t);
    pressure *= std::sqrt(2)*M_PI;
    return pressure;
}

Vector WaveO1TestCase <LShaped_Test1>
:: velocity_sol (const Vector& x, const double t) const
{
    double r = radius(x);
    double theta = polar_angle(x);

    double coeff = m_gamma-1;
    Vector v(m_dim);

    v(0) = -m_gamma*std::pow(r, coeff)*cos(coeff*theta)
                *(1./4. - x[0]*x[0]);
    v(0) += std::pow(r, m_gamma)*sin(m_gamma*theta)*(- 2*x[0]);
    v(0) *= (1./4. - x[1]*x[1]);

    v(1) = m_gamma*std::pow(r, coeff)*sin(coeff*theta)
            *(1./4. - x[1]*x[1]);
    v(1) += std::pow(r, m_gamma)*sin(m_gamma*theta)*(- 2*x[1]);
    v(1) *= (1./4. - x[0]*x[0]);

    v *= -sin(std::sqrt(2)*M_PI*t);
    return v;
}

double WaveO1TestCase <LShaped_Test1>
:: source (const Vector& x, const double t) const
{
    double r = radius(x);
    double theta = polar_angle(x);
    double coeff = m_gamma-1;

    double f = -2*M_PI*M_PI*std::pow(r, m_gamma)
            *sin(m_gamma*theta)
            *(1./4. - x[0]*x[0])*(1./4. - x[1]*x[1])
            *sin(std::sqrt(2)*M_PI*t);

    double f1 = 4*m_gamma*std::pow(r, coeff)
            *cos(coeff*theta)*x[0];
    f1 -= 2*pow(r,m_gamma)*sin(m_gamma*theta);
    f1 *= (1./4. - x[1]*x[1]);

    double f2 = -4*m_gamma*pow(r, coeff)
            *sin(coeff*theta)*x[1];
    f2 -= 2*pow(r,m_gamma)*sin(m_gamma*theta);
    f2 *= (1./4. - x[0]*x[0]);

    f -= (f1 + f2)*sin(std::sqrt(2)*M_PI*t);

    return f;
}*/

void WaveO1TestCase <LShaped_Test1>
:: set_bdry_dirichlet(Array<int>& bdr_marker) const
{
    bdr_marker = 1;
}

double WaveO1TestCase <LShaped_Test1>
:: pressure_sol(const Vector& x, const double t) const
{
    double r = radius(x);
    double theta = polar_angle(x);

    double pressure
            = std::pow(r, m_gamma)*sin(m_gamma*theta)
            *(1 - x[0]*x[0])*(1 - x[1]*x[1])
            *cos(std::sqrt(2)*M_PI*t);
    pressure *= std::sqrt(2)*M_PI;
    return pressure;
}

Vector WaveO1TestCase <LShaped_Test1>
:: velocity_sol (const Vector& x, const double t) const
{
    double r = radius(x);
    double theta = polar_angle(x);

    double coeff = m_gamma-1;
    Vector v(m_dim);

    v(0) = -m_gamma*std::pow(r, coeff)*cos(coeff*theta)
                *(1 - x[0]*x[0]);
    v(0) += std::pow(r, m_gamma)*sin(m_gamma*theta)*(- 2*x[0]);
    v(0) *= (1 - x[1]*x[1]);

    v(1) = m_gamma*std::pow(r, coeff)*sin(coeff*theta)
            *(1 - x[1]*x[1]);
    v(1) += std::pow(r, m_gamma)*sin(m_gamma*theta)*(- 2*x[1]);
    v(1) *= (1 - x[0]*x[0]);

    v *= -sin(std::sqrt(2)*M_PI*t);
    return v;
}

double WaveO1TestCase <LShaped_Test1>
:: source (const Vector& x, const double t) const
{
    double r = radius(x);
    double theta = polar_angle(x);
    double coeff = m_gamma-1;

    double f = -2*M_PI*M_PI*std::pow(r, m_gamma)
            *sin(m_gamma*theta)
            *(1 - x[0]*x[0])*(1 - x[1]*x[1])
            *sin(std::sqrt(2)*M_PI*t);

    double f1 = 4*m_gamma*std::pow(r, coeff)
            *cos(coeff*theta)*x[0];
    f1 -= 2*pow(r,m_gamma)*sin(m_gamma*theta);
    f1 *= (1 - x[1]*x[1]);

    double f2 = -4*m_gamma*pow(r, coeff)
            *sin(coeff*theta)*x[1];
    f2 -= 2*pow(r,m_gamma)*sin(m_gamma*theta);
    f2 *= (1 - x[0]*x[0]);

    f -= (f1 + f2)*sin(std::sqrt(2)*M_PI*t);

    return f;
}


// Case: LShaped_Test2
// Homogeneous Dirichlet boundary
// Non-homogeneous source
double WaveO1TestCase <LShaped_Test2>
:: radius(const Vector& x) const
{
    return x.Norml2();
}

double WaveO1TestCase <LShaped_Test2>
:: polar_angle(const Vector& x) const
{
    double theta = std::atan2(-x[0],x[1])
                + (x[0] >= 0 && x[1] <= 0 ? 2*M_PI : 0);
    return theta;
}

/*double WaveO1TestCase <LShaped_Test2>
:: pressure_sol(const Vector& x, const double t) const
{
    double r = radius(x);
    double theta = polar_angle(x);

    double pressure
            = std::pow(r, m_gamma)*sin(m_gamma*theta)
            *(1./4. - x[0]*x[0])*(1./4. - x[1]*x[1]);
    pressure *= M_PI*sin(2*M_PI*t);
    return pressure;
}

Vector WaveO1TestCase <LShaped_Test2>
:: velocity_sol (const Vector& x, const double t) const
{
    double r = radius(x);
    double theta = polar_angle(x);

    double coeff = m_gamma-1;
    Vector v(m_dim);

    v(0) = -m_gamma*std::pow(r, coeff)*cos(coeff*theta)
                *(1./4. - x[0]*x[0]);
    v(0) += std::pow(r, m_gamma)*sin(m_gamma*theta)*(- 2*x[0]);
    v(0) *= (1./4. - x[1]*x[1]);

    v(1) = m_gamma*std::pow(r, coeff)*sin(coeff*theta)
            *(1./4. - x[1]*x[1]);
    v(1) += std::pow(r, m_gamma)*sin(m_gamma*theta)*(- 2*x[1]);
    v(1) *= (1./4. - x[0]*x[0]);

    v *= -sin(M_PI*t)*sin(M_PI*t);
    return v;
}

double WaveO1TestCase <LShaped_Test2>
:: source (const Vector& x, const double t) const
{
    double r = radius(x);
    double theta = polar_angle(x);
    double coeff = m_gamma-1;

    double f = 2*M_PI*M_PI
            *std::pow(r, m_gamma)
            *sin(m_gamma*theta)
            *(1./4. - x[0]*x[0])*(1./4. - x[1]*x[1])
            *cos(2*M_PI*t);

    double f1 = 4*m_gamma*std::pow(r, coeff)
            *cos(coeff*theta)*x[0];
    f1 -= 2*pow(r,m_gamma)*sin(m_gamma*theta);
    f1 *= (1./4. - x[1]*x[1]);

    double f2 = -4*m_gamma*pow(r, coeff)
            *sin(coeff*theta)*x[1];
    f2 -= 2*pow(r,m_gamma)*sin(m_gamma*theta);
    f2 *= (1./4. - x[0]*x[0]);

    f -= (f1 + f2)*sin(M_PI*t)*sin(M_PI*t);

    return f;
}*/

void WaveO1TestCase <LShaped_Test2>
:: set_bdry_dirichlet(Array<int>& bdr_marker) const
{
    bdr_marker = 1;
}

double WaveO1TestCase <LShaped_Test2>
:: pressure_sol(const Vector& x, const double t) const
{
    double r = radius(x);
    double theta = polar_angle(x);

    double pressure
            = std::pow(r, m_gamma)*sin(m_gamma*theta)
            *(1 - x[0]*x[0])*(1 - x[1]*x[1]);
    pressure *= M_PI*sin(2*M_PI*t);
    return pressure;
}

Vector WaveO1TestCase <LShaped_Test2>
:: velocity_sol (const Vector& x, const double t) const
{
    double r = radius(x);
    double theta = polar_angle(x);

    double coeff = m_gamma-1;
    Vector v(m_dim);

    v(0) = -m_gamma*std::pow(r, coeff)*cos(coeff*theta)
                *(1 - x[0]*x[0]);
    v(0) += std::pow(r, m_gamma)*sin(m_gamma*theta)*(- 2*x[0]);
    v(0) *= (1 - x[1]*x[1]);

    v(1) = m_gamma*std::pow(r, coeff)*sin(coeff*theta)
            *(1 - x[1]*x[1]);
    v(1) += std::pow(r, m_gamma)*sin(m_gamma*theta)*(- 2*x[1]);
    v(1) *= (1 - x[0]*x[0]);

    v *= -sin(M_PI*t)*sin(M_PI*t);
    return v;
}

double WaveO1TestCase <LShaped_Test2>
:: source (const Vector& x, const double t) const
{
    double r = radius(x);
    double theta = polar_angle(x);
    double coeff = m_gamma-1;

    double f = 2*M_PI*M_PI
            *std::pow(r, m_gamma)
            *sin(m_gamma*theta)
            *(1 - x[0]*x[0])*(1 - x[1]*x[1])
            *cos(2*M_PI*t);

    double f1 = 4*m_gamma*std::pow(r, coeff)
            *cos(coeff*theta)*x[0];
    f1 -= 2*pow(r,m_gamma)*sin(m_gamma*theta);
    f1 *= (1 - x[1]*x[1]);

    double f2 = -4*m_gamma*pow(r, coeff)
            *sin(coeff*theta)*x[1];
    f2 -= 2*pow(r,m_gamma)*sin(m_gamma*theta);
    f2 *= (1 - x[0]*x[0]);

    f -= (f1 + f2)*sin(M_PI*t)*sin(M_PI*t);

    return f;
}

/*
// Case: LShaped_Test3
// Non-homogeneous Dirichlet boundary
// Non-homogeneous source
double WaveO1TestCase <LShaped_Test3>
:: radius(const Vector& x) const
{
    return x.Norml2();
}

double WaveO1TestCase <LShaped_Test3>
:: polar_angle(const Vector& x) const
{
    double theta = std::atan2(-x[0],x[1])
                + (x[0] >= 0 && x[1] <= 0 ? 2*M_PI : 0);
    return theta;
}

double WaveO1TestCase <LShaped_Test3>
:: pressure_sol(const Vector& x, const double t) const
{
    double r = radius(x);
    double theta = polar_angle(x);

    double pressure
            = std::pow(r, m_gamma)*sin(m_gamma*theta);
    pressure *= std::sqrt(2)*M_PI*cos(std::sqrt(2)*M_PI*t);
    return pressure;
}

Vector WaveO1TestCase <LShaped_Test3>
:: velocity_sol (const Vector& x, const double t) const
{
    double r = radius(x);
    double theta = polar_angle(x);

    double coeff = m_gamma-1;
    Vector v(m_dim);

    v(0) = -m_gamma*std::pow(r, coeff)*cos(coeff*theta);
    v(1) = m_gamma*std::pow(r, coeff)*sin(coeff*theta);

    v *= -sin(std::sqrt(2)*M_PI*t);
    return v;
}

double WaveO1TestCase <LShaped_Test3>
:: source (const Vector& x, const double t) const
{
    double r = radius(x);
    double theta = polar_angle(x);

    double f = std::pow(r, m_gamma)*sin(m_gamma*theta)*
            (-2*M_PI*M_PI)*sin(std::sqrt(2)*M_PI*t);

    return f;
}

void WaveO1TestCase <LShaped_Test3>
:: set_bdry_dirichlet(Array<int>& bdr_marker) const
{
    bdr_marker = 1;
}

// Case: LShaped_Test4
// Non-homogeneous Dirichlet boundary
// Non-homogeneous source
double WaveO1TestCase <LShaped_Test4>
:: radius(const Vector& x) const
{
    return x.Norml2();
}

double WaveO1TestCase <LShaped_Test4>
:: polar_angle(const Vector& x) const
{
    double theta = std::atan2(-x[0],x[1])
                + (x[0] >= 0 && x[1] <= 0 ? 2*M_PI : 0);
    return theta;
}

double WaveO1TestCase <LShaped_Test4>
:: pressure_sol(const Vector& x, const double t) const
{
    double r = radius(x);
    double theta = polar_angle(x);

    double pressure
            = std::pow(r, m_gamma)*sin(m_gamma*theta);
    pressure *= M_PI*sin(2*M_PI*t);
    return pressure;
}

Vector WaveO1TestCase <LShaped_Test4>
:: velocity_sol (const Vector& x, const double t) const
{
    double r = radius(x);
    double theta = polar_angle(x);

    double coeff = m_gamma-1;
    Vector v(m_dim);

    v(0) = -m_gamma*std::pow(r, coeff)*cos(coeff*theta);
    v(1) = m_gamma*std::pow(r, coeff)*sin(coeff*theta);

    v *= -sin(M_PI*t)*sin(M_PI*t);
    return v;
}

double WaveO1TestCase <LShaped_Test4>
:: source (const Vector& x, const double t) const
{
    double r = radius(x);
    double theta = polar_angle(x);

    double f = std::pow(r, m_gamma)*sin(m_gamma*theta)*
            2*M_PI*M_PI*cos(2*M_PI*t);

    return f;
}

void WaveO1TestCase <LShaped_Test4>
:: set_bdry_dirichlet(Array<int>& bdr_marker) const
{
    bdr_marker = 1;
}


// Case: LShaped_Test4
// Homogeneous Dirichlet boundary
// Non-homogeneous source
double WaveO1TestCase <LShaped_Test5>
:: pressure_sol(const Vector& x, const double t) const
{
    double pressure
            = sin(2*M_PI*x(0))*sin(2*M_PI*x(1));
    pressure *= M_PI*sin(2*M_PI*t);
    return pressure;
}

Vector WaveO1TestCase <LShaped_Test5>
:: velocity_sol (const Vector& x, const double t) const
{
    Vector v(m_dim);
    v(0) = cos(2*M_PI*x(0))*sin(2*M_PI*x(1));
    v(1) = sin(2*M_PI*x(0))*cos(2*M_PI*x(1));
    v *= -2*M_PI*sin(M_PI*t)*sin(M_PI*t);
    return v;
}

double WaveO1TestCase <LShaped_Test5>
:: source (const Vector& x, const double t) const
{
    double f = 2*M_PI*M_PI*(cos(2*M_PI*t)
                          + 4*sin(M_PI*t)*sin(M_PI*t));
    f *= sin(2*M_PI*x(0))*sin(2*M_PI*x(1));
    return f;
}

void WaveO1TestCase <LShaped_Test5>
:: set_bdry_dirichlet(Array<int>& bdr_marker) const
{
    bdr_marker = 1;
}
*/

// Case: SquareTwoShaped_Test1
double WaveO1TestCase <SquareTwoPiece_Test1>
:: medium(const Vector& x) const
{
    double c = 0;
    if (x(0) <= m_xs) {
        c = 1;
    }
    else {
        c = 3;
    }

    return c;
}

double WaveO1TestCase <SquareTwoPiece_Test1>
:: init_pressure(const Vector&) const
{
    return 0;
}

Vector WaveO1TestCase <SquareTwoPiece_Test1>
:: init_velocity(const Vector& x) const
{
    double r2 = x.DistanceSquaredTo(m_x0);
    double invLambda2 = 1./(m_lambda*m_lambda);

    Vector v(m_dim);
    v(0) = (x(0) - m_x0(0));
    v(1) = (x(1) - m_x0(1));
    v *= 2*exp(-r2*invLambda2)*invLambda2;

    return v;
}

double WaveO1TestCase <SquareTwoPiece_Test1>
:: source (const Vector&, const double) const
{
    return 0;
}

void WaveO1TestCase <SquareTwoPiece_Test1>
:: set_bdry_dirichlet(Array<int>& bdr_marker) const
{
    bdr_marker = 1;
}


// Case: SquareTwoShaped_Test2
double WaveO1TestCase <SquareTwoPiece_Test2>
:: medium(const Vector& x) const
{
    double c = 0;
    if (x(0) > m_xs && x(1) > 1) {
        c = 3;
    }
    else if (x(0) <= m_xs && x(1) > 1) {
        c = 1;
    }
    else if (x(0) <= m_xs && x(1) <= 1) {
        c = 3;
    }
    else {
        c = 1;
    }

    return c;
}

double WaveO1TestCase <SquareTwoPiece_Test2>
:: init_pressure(const Vector&) const
{
    return 0;
}

Vector WaveO1TestCase <SquareTwoPiece_Test2>
:: init_velocity(const Vector& x) const
{
    double r2 = x.DistanceSquaredTo(m_x0);
    double invLambda2 = 1./(m_lambda*m_lambda);

    Vector v(m_dim);
    v(0) = (x(0) - m_x0(0));
    v(1) = (x(1) - m_x0(1));
    v *= 2*exp(-r2*invLambda2)*invLambda2;

    return v;
}

double WaveO1TestCase <SquareTwoPiece_Test2>
:: source (const Vector&, const double) const
{
    return 0;
}

void WaveO1TestCase <SquareTwoPiece_Test2>
:: set_bdry_dirichlet(Array<int>& bdr_marker) const
{
    bdr_marker = 1;
}

// End of file
