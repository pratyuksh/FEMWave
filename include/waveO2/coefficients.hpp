#ifndef WAVEO2_COEFFICIENTS_HPP
#define WAVEO2_COEFFICIENTS_HPP

#include "test_cases.hpp"


//-------------------------------//
//  Exact Solution Coefficients  //
//-------------------------------//

//! Class for WaveO2
//! exact solution coefficient
class WaveO2ExactSolutionCoeff : public Coefficient
{
public:
    WaveO2ExactSolutionCoeff
    (std::shared_ptr<WaveO2TestCases> testCase)
        : m_testCase (testCase) {}
    
    virtual double Eval(ElementTransformation &,
                        const IntegrationPoint &);
    
    double Eval(ElementTransformation &,
                const IntegrationPoint &,
                double);
    
protected:
    std::shared_ptr<WaveO2TestCases> m_testCase;
};

//! Class for WaveO2
//! exact solution gradient in t coefficient
class WaveO2ExactGradtSolutionCoeff : public Coefficient
{
public:
    WaveO2ExactGradtSolutionCoeff
    (std::shared_ptr<WaveO2TestCases> testCase)
        : m_testCase (testCase) {}

    virtual double Eval(ElementTransformation &,
                        const IntegrationPoint &);

    double Eval(ElementTransformation &,
                const IntegrationPoint &,
                double);

protected:
    std::shared_ptr<WaveO2TestCases> m_testCase;
};

//! Class for WaveO2
//! exact solution gradient in x coefficient
class WaveO2ExactGradxSolutionCoeff
        : public VectorCoefficient
{
public:
    WaveO2ExactGradxSolutionCoeff
    (std::shared_ptr<WaveO2TestCases> testCase)
        : VectorCoefficient (testCase->get_dim()),
          m_testCase (testCase) {}

    virtual void Eval(Vector&, ElementTransformation&,
                      const IntegrationPoint&);
    
    void Eval(Vector&, ElementTransformation&,
              const IntegrationPoint&, double);
    
private:
    std::shared_ptr<WaveO2TestCases> m_testCase;
};


//---------------------------------//
//  Initial Solution Coefficients  //
//---------------------------------//

//! Class for WaveO2
//! initial solution coefficient
class WaveO2InitialSolutionCoeff : public Coefficient
{
public:
    WaveO2InitialSolutionCoeff
    (std::shared_ptr<WaveO2TestCases> testCase)
        : m_testCase (testCase) {}

    virtual double Eval(ElementTransformation &,
                        const IntegrationPoint &);

protected:
    std::shared_ptr<WaveO2TestCases> m_testCase;
};

//! Class for WaveO2
//! initial solution gradient in t coefficient
class WaveO2InitialGradtSolutionCoeff : public Coefficient
{
public:
    WaveO2InitialGradtSolutionCoeff
    (std::shared_ptr<WaveO2TestCases> testCase)
        : m_testCase (testCase) {}

    virtual double Eval(ElementTransformation &,
                        const IntegrationPoint &);

protected:
    std::shared_ptr<WaveO2TestCases> m_testCase;
};


//------------------------//
//  Boundary Coefficient  //
//------------------------//

//! Class for WaveO2
//! boundary solution coefficient
class WaveO2BdrySolutionCoeff : public Coefficient
{
public:
    WaveO2BdrySolutionCoeff
    (std::shared_ptr<WaveO2TestCases> testCase)
        : m_testCase (testCase) {}

    virtual double Eval(ElementTransformation &,
                        const IntegrationPoint &);
    
    double Eval(ElementTransformation &,
                const IntegrationPoint &,
                double);

private:
    std::shared_ptr<WaveO2TestCases> m_testCase;
};


//-----------------------//
//  Source Coefficients  //
//-----------------------//

//! Class for WaveO2 source coefficient
class WaveO2SourceCoeff : public Coefficient
{
public:
    WaveO2SourceCoeff
    (std::shared_ptr<WaveO2TestCases> testCase)
        : m_testCase (testCase) {}

    virtual double Eval(ElementTransformation &,
                        const IntegrationPoint &);
    
    double Eval(ElementTransformation &,
                const IntegrationPoint &,
                double);

private:
    std::shared_ptr<WaveO2TestCases> m_testCase;
};


//-------------------------//
//  Material Coefficients  //
//-------------------------//

//! Class for WaveO2 medium coefficient
class WaveO2MediumCoeff : public Coefficient
{
public:
    WaveO2MediumCoeff
    (std::shared_ptr<WaveO2TestCases> testCase)
        : m_testCase (testCase) {}

    virtual double Eval(ElementTransformation &,
                        const IntegrationPoint &);

private:
    std::shared_ptr<WaveO2TestCases> m_testCase;
};

//! Class for WaveO2 square medium coefficient
class WaveO2SqMediumCoeff : public Coefficient
{
public:
    WaveO2SqMediumCoeff
    (std::shared_ptr<WaveO2TestCases> testCase)
        : m_testCase (testCase) {}

    virtual double Eval(ElementTransformation &,
                        const IntegrationPoint &);

private:
    std::shared_ptr<WaveO2TestCases> m_testCase;
};

#endif /// WAVEO2_COEFFICIENTS_HPP
