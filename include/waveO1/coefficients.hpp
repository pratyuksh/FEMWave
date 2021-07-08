#ifndef WAVEO1_COEFFICIENTS_HPP
#define WAVEO1_COEFFICIENTS_HPP

#include "test_cases.hpp"


//-------------------------------//
//  Exact Solution Coefficients  //
//-------------------------------//

// class for WaveO1
// exact pressure coefficient
class WaveO1ExactPressureCoeff : public Coefficient
{
public:
    WaveO1ExactPressureCoeff
    (std::shared_ptr<WaveO1TestCases> testCase)
        : m_testCase (testCase) {}
    
    virtual double Eval(ElementTransformation &,
                        const IntegrationPoint &);
    
    double Eval(ElementTransformation &,
                const IntegrationPoint &,
                double);
    
protected:
    std::shared_ptr<WaveO1TestCases> m_testCase;
};


// class for WaveO1
// exact velocity coefficient
class WaveO1ExactVelocityCoeff : public VectorCoefficient
{
public:
    WaveO1ExactVelocityCoeff (std::shared_ptr
                            <WaveO1TestCases> testCase)
        : VectorCoefficient (testCase->get_dim()),
          m_testCase (testCase) {}

    virtual void Eval(Vector&, ElementTransformation&,
                      const IntegrationPoint&);
    
    void Eval(Vector&, ElementTransformation&,
              const IntegrationPoint&, double);
    
private:
    std::shared_ptr<WaveO1TestCases> m_testCase;
};


//---------------------------------//
//  Initial Solution Coefficients  //
//---------------------------------//

// class for WaveO1
// initial pressure coefficient
class WaveO1InitialPressureCoeff : public Coefficient
{
public:
    WaveO1InitialPressureCoeff(std::shared_ptr
                             <WaveO1TestCases> testCase)
        : m_testCase (testCase) {}

    virtual double Eval(ElementTransformation &,
                        const IntegrationPoint &);

protected:
    std::shared_ptr<WaveO1TestCases> m_testCase;
};

// class for WaveO1
// initial velocity coefficient
class WaveO1InitialVelocityCoeff : public VectorCoefficient
{
public:
    WaveO1InitialVelocityCoeff (std::shared_ptr
                              <WaveO1TestCases> testCase)
        : VectorCoefficient (testCase->get_dim()),
          m_testCase (testCase) {}

    virtual void Eval(Vector&, ElementTransformation&,
                      const IntegrationPoint&);
    
private:
    std::shared_ptr<WaveO1TestCases> m_testCase;
};


//------------------------//
//  Boundary Coefficient  //
//------------------------//

// class for WaveO1
// boundary pressure coefficient
class WaveO1BdryPressureCoeff : public Coefficient
{
public:
    WaveO1BdryPressureCoeff (std::shared_ptr
                           <WaveO1TestCases> testCase)
        : m_testCase (testCase) {}

    virtual double Eval(ElementTransformation &,
                        const IntegrationPoint &);
    
    double Eval(ElementTransformation &,
                const IntegrationPoint &,
                double);

private:
    std::shared_ptr<WaveO1TestCases> m_testCase;
};

// class for WaveO1
// boundary velocity coefficient
class WaveO1BdryVelocityCoeff : public VectorCoefficient
{
public:
    WaveO1BdryVelocityCoeff (std::shared_ptr
                           <WaveO1TestCases> testCase)
        : VectorCoefficient (testCase->get_dim()),
          m_testCase (testCase) {}

    virtual void Eval(Vector&, ElementTransformation&,
                      const IntegrationPoint&);
    
    void Eval(Vector&, ElementTransformation&,
              const IntegrationPoint&, double);
    
private:
    std::shared_ptr<WaveO1TestCases> m_testCase;
};


//-----------------------//
//  Source Coefficients  //
//-----------------------//

// class for WaveO1 source coefficient
class WaveO1SourceCoeff : public Coefficient
{
public:
    WaveO1SourceCoeff (std::shared_ptr
                     <WaveO1TestCases> testCase)
        : m_testCase (testCase) {}

    virtual double Eval(ElementTransformation &,
                        const IntegrationPoint &);
    
    double Eval(ElementTransformation &,
                const IntegrationPoint &,
                double);

private:
    std::shared_ptr<WaveO1TestCases> m_testCase;
};


//-------------------------//
//  Material Coefficients  //
//-------------------------//

// class for WaveO1 medium coefficient
class WaveO1MediumCoeff : public Coefficient
{
public:
    WaveO1MediumCoeff (std::shared_ptr
                     <WaveO1TestCases> testCase)
        : m_testCase (testCase) {}

    virtual double Eval(ElementTransformation &,
                        const IntegrationPoint &);

private:
    std::shared_ptr<WaveO1TestCases> m_testCase;
};

// class for WaveO1 inverse square medium coefficient
class WaveO1InvSqMediumCoeff : public Coefficient
{
public:
    WaveO1InvSqMediumCoeff (std::shared_ptr
                          <WaveO1TestCases> testCase)
        : m_testCase (testCase) {}

    virtual double Eval(ElementTransformation &,
                        const IntegrationPoint &);

private:
    std::shared_ptr<WaveO1TestCases> m_testCase;
};

#endif /// WAVEO1_COEFFICIENTS_HPP
