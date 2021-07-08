#include "../../include/waveO2/coefficients.hpp"


//-------------------------------//
//  Exact Solution Coefficients  //
//-------------------------------//

//! WaveO2 exact solution coefficient
double WaveO2ExactSolutionCoeff
:: Eval (ElementTransformation& T,
         const IntegrationPoint& ip)
{
    double x[3];
    Vector transip(x, 3);
    T.Transform(ip, transip);
    return (m_testCase->exact_sol(transip, GetTime()));
}

double WaveO2ExactSolutionCoeff
:: Eval (ElementTransformation& T,
         const IntegrationPoint& ip,
         double t)
{
    double x[3];
    Vector transip(x, 3);
    T.Transform(ip, transip);
    return (m_testCase->exact_sol(transip, t));
}

//! WaveO2 exact solution gradient in t coefficient
double WaveO2ExactGradtSolutionCoeff
:: Eval (ElementTransformation& T,
         const IntegrationPoint& ip)
{
    double x[3];
    Vector transip(x, 3);
    T.Transform(ip, transip);
    return (m_testCase->exact_gradt_sol(transip, GetTime()));
}

double WaveO2ExactGradtSolutionCoeff
:: Eval (ElementTransformation& T,
         const IntegrationPoint& ip,
         double t)
{
    double x[3];
    Vector transip(x, 3);
    T.Transform(ip, transip);
    return (m_testCase->exact_gradt_sol(transip, t));
}

//! WaveO2 exact solution gradient in x coefficient
void WaveO2ExactGradxSolutionCoeff
:: Eval (Vector& v,
         ElementTransformation& T,
         const IntegrationPoint& ip)
{
    double x[3];
    Vector transip(x, 3);
    T.Transform(ip, transip);
    v = m_testCase->exact_gradx_sol(transip, GetTime());
}

void WaveO2ExactGradxSolutionCoeff
:: Eval (Vector& v,
         ElementTransformation& T,
         const IntegrationPoint& ip,
         double t)
{
    double x[3];
    Vector transip(x, 3);
    T.Transform(ip, transip);
    v = m_testCase->exact_gradx_sol(transip, t);
}


//---------------------------------//
//  Initial Solution Coefficients  //
//---------------------------------//

//! WaveO2 initial solution coefficient
double WaveO2InitialSolutionCoeff
:: Eval (ElementTransformation& T,
         const IntegrationPoint& ip)
{
    double x[3];
    Vector transip(x, 3);
    T.Transform(ip, transip);
    return (m_testCase->init_sol(transip));
}

//! WaveO2 initial solution gradient in t coefficient
double WaveO2InitialGradtSolutionCoeff
:: Eval (ElementTransformation& T,
         const IntegrationPoint& ip)
{
    double x[3];
    Vector transip(x, 3);
    T.Transform(ip, transip);
    return (m_testCase->init_gradt_sol(transip));
}


//------------------------//
//  Boundary Coefficient  //
//------------------------//

//! WaveO2 boundary solution coefficient
double WaveO2BdrySolutionCoeff
:: Eval (ElementTransformation& T,
         const IntegrationPoint& ip)
{
    double x[3];
    Vector transip(x, 3);
    T.Transform(ip, transip);
    return (m_testCase->bdry_sol(transip, GetTime()));
}

double WaveO2BdrySolutionCoeff
:: Eval (ElementTransformation& T,
         const IntegrationPoint& ip,
         double t)
{
    double x[3];
    Vector transip(x, 3);
    T.Transform(ip, transip);
    return (m_testCase->bdry_sol(transip, t));
}


//-----------------------//
//  Source Coefficients  //
//-----------------------//

//! WaveO2 source coefficient
double WaveO2SourceCoeff
:: Eval (ElementTransformation& T,
         const IntegrationPoint& ip)
{
    double x[3];
    Vector transip(x, 3);
    T.Transform(ip, transip);
    return m_testCase->source(transip, GetTime());
}

double WaveO2SourceCoeff
:: Eval (ElementTransformation& T,
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

//! WaveO2 medium coefficient
double WaveO2MediumCoeff
:: Eval (ElementTransformation& T,
         const IntegrationPoint& ip)
{
    double x[3];
    Vector transip(x, 3);
    T.Transform(ip, transip);
    return m_testCase->medium(transip);
}

//! WaveO2 square medium coefficient
double WaveO2SqMediumCoeff
:: Eval (ElementTransformation& T,
         const IntegrationPoint& ip)
{
    double x[3];
    Vector transip(x, 3);
    T.Transform(ip, transip);
    auto val = m_testCase->medium(transip);
    return (val*val);
}

// End of file
