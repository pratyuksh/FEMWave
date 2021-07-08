#ifndef WAVEO2_TEST_CASES_HPP
#define WAVEO2_TEST_CASES_HPP

#include "../core/config.hpp"
#include "mfem.hpp"

using namespace mfem;

//! Cases
enum {UnitSquare_Test1,
      UnitSquare_Test2,
      LShaped_Test1};


//! Abstract Base class for the WaveO2 test cases
class WaveO2TestCases
{
public:

    virtual ~WaveO2TestCases() = default;
    
    virtual double medium(const Vector&) const = 0;
    
    virtual double exact_sol(const Vector&,
                             const double) const = 0;
    
    virtual double exact_gradt_sol(const Vector&,
                                   const double) const = 0;
    
    virtual Vector exact_gradx_sol(const Vector&,
                                   const double) const = 0;
    
    virtual double init_sol(const Vector&) const = 0;

    virtual double init_gradt_sol(const Vector&) const = 0;
    
    virtual double bdry_sol(const Vector&,
                            const double) const = 0;
    
    virtual double source(const Vector&,
                          const double) const = 0;
    
    virtual void set_bdry_dirichlet(Array<int>&) const = 0;
    
    virtual int get_dim() const = 0;
};


//! class WaveO2TestCase, base template
template<int ProblemType>
class WaveO2TestCase;

//! Case: UnitSquare_Test1
//! Smooth solution
//! Homogeneous Dirichlet boundary
//! Homogeneous source
template<>
class WaveO2TestCase <UnitSquare_Test1>
        : public WaveO2TestCases
{
public:
    explicit WaveO2TestCase (const nlohmann::json& config)
        : m_config(config), m_dim(2) {
    }

    double exact_sol(const Vector&,
                     const double) const override;
    
    double exact_gradt_sol(const Vector&,
                           const double) const override;
    
    Vector exact_gradx_sol(const Vector&,
                           const double) const override;

    double source(const Vector&,
                  const double) const override;

    void set_bdry_dirichlet(Array<int>&) const override;

    double medium(const Vector& x) const override {
        return m_c;
    }

    double init_sol(const Vector& x) const override {
        return exact_sol(x,0);
    }

    double init_gradt_sol(const Vector& x) const override {
        return exact_gradt_sol(x,0);
    }

    double bdry_sol(const Vector& x,
                    const double t) const override {
        return exact_sol(x, t);
    }

    inline int get_dim() const override {
        return m_dim;
    }

private:
    const nlohmann::json& m_config;
    int m_dim;

    double m_c = 1;
};

#endif /// WAVEO2_TEST_CASES_HPP
