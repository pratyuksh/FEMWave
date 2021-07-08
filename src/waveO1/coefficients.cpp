#include "../../include/waveO1/coefficients.hpp"


//-------------------------------//
//  Exact Solution Coefficients  //
//-------------------------------//

// WaveO1 exact pressure coefficient
double WaveO1ExactPressureCoeff
:: Eval (ElementTransformation& T,
         const IntegrationPoint& ip)
{
    double x[3];
    Vector transip(x, 3);
    T.Transform(ip, transip);
    return (m_testCase->pressure_sol(transip, GetTime()));
}

double WaveO1ExactPressureCoeff
:: Eval (ElementTransformation& T,
         const IntegrationPoint& ip,
         double t)
{
    double x[3];
    Vector transip(x, 3);
    T.Transform(ip, transip);
    return (m_testCase->pressure_sol(transip, t));
}

// WaveO1 exact velocity coefficient
void WaveO1ExactVelocityCoeff
:: Eval (Vector& v,
         ElementTransformation& T,
         const IntegrationPoint& ip)
{
    double x[3];
    Vector transip(x, 3);
    T.Transform(ip, transip);
    v = m_testCase->velocity_sol(transip, GetTime());
}

void WaveO1ExactVelocityCoeff
:: Eval (Vector& v,
         ElementTransformation& T,
         const IntegrationPoint& ip,
         double t)
{
    double x[3];
    Vector transip(x, 3);
    T.Transform(ip, transip);
    v = m_testCase->velocity_sol(transip, t);
}


//---------------------------------//
//  Initial Solution Coefficients  //
//---------------------------------//

// WaveO1 initial pressure coefficient
double WaveO1InitialPressureCoeff
:: Eval (ElementTransformation& T,
         const IntegrationPoint& ip)
{
    double x[3];
    Vector transip(x, 3);
    T.Transform(ip, transip);
    return (m_testCase->init_pressure(transip));
}

// WaveO1 initial velocity coefficient
void WaveO1InitialVelocityCoeff
:: Eval (Vector& v,
         ElementTransformation& T,
         const IntegrationPoint& ip)
{
    double x[3];
    Vector transip(x, 3);
    T.Transform(ip, transip);
    v = m_testCase->init_velocity(transip);
}


//------------------------//
//  Boundary Coefficient  //
//------------------------//

// WaveO1 boundary pressure coefficient
double WaveO1BdryPressureCoeff
:: Eval (ElementTransformation& T,
         const IntegrationPoint& ip)
{
    double x[3];
    Vector transip(x, 3);
    T.Transform(ip, transip);
    return (m_testCase->bdry_pressure(transip, GetTime()));
}

double WaveO1BdryPressureCoeff
:: Eval (ElementTransformation& T,
         const IntegrationPoint& ip,
         double t)
{
    double x[3];
    Vector transip(x, 3);
    T.Transform(ip, transip);
    return (m_testCase->bdry_pressure(transip, t));
}

// WaveO1 boundary velocity coefficient
void WaveO1BdryVelocityCoeff
:: Eval (Vector& v,
         ElementTransformation& T,
         const IntegrationPoint& ip)
{
    double x[3];
    Vector transip(x, 3);
    T.Transform(ip, transip);
    v = m_testCase->bdry_velocity(transip, GetTime());
}

void WaveO1BdryVelocityCoeff
:: Eval (Vector& v,
         ElementTransformation& T,
         const IntegrationPoint& ip,
         double t)
{
    double x[3];
    Vector transip(x, 3);
    T.Transform(ip, transip);
    v = m_testCase->bdry_velocity(transip, t);
}


//-----------------------//
//  Source Coefficients  //
//-----------------------//

// WaveO1 source coefficient
double WaveO1SourceCoeff :: Eval (ElementTransformation& T,
                                const IntegrationPoint& ip)
{
    double x[3];
    Vector transip(x, 3);
    T.Transform(ip, transip);
    return m_testCase->source(transip, GetTime());
}

double WaveO1SourceCoeff :: Eval (ElementTransformation& T,
                                const IntegrationPoint& ip,
                                double t)
{
    double x[3];
    Vector transip(x, 3);
    T.Transform(ip, transip);
    return m_testCase->source(transip, t);
}


//-----------------------//
//  Medium Coefficients  //
//-----------------------//

// WaveO1 medium coefficient
double WaveO1MediumCoeff :: Eval (ElementTransformation& T,
                                const IntegrationPoint& ip)
{
    double x[3];
    Vector transip(x, 3);
    T.Transform(ip, transip);
    return m_testCase->medium(transip);
}

// WaveO1 inverse square medium coefficient
double WaveO1InvSqMediumCoeff :: Eval (ElementTransformation& T,
                                const IntegrationPoint& ip)
{
    double x[3];
    Vector transip(x, 3);
    T.Transform(ip, transip);
    auto val = m_testCase->medium(transip);
    return 1./(val*val);
}

// End of file
