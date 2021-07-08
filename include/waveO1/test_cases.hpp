#ifndef WAVEO1_TEST_CASES_HPP
#define WAVEO1_TEST_CASES_HPP

#include "../core/config.hpp"
#include "mfem.hpp"

using namespace mfem;

// Cases
enum {UnitSquare_Test1,
      UnitSquare_Test2,
      UnitSquare_Test3,
      UnitSquare_Test4,
      GammaShaped_Test1,
      LShaped_Test1,
      LShaped_Test2,
      LShaped_Test3,
      LShaped_Test4,
      LShaped_Test5,
      SquareTwoPiece_Test1,
      SquareTwoPiece_Test2};


// Abstract Base class for the WaveO1 test cases
class WaveO1TestCases
{
public:

    virtual ~WaveO1TestCases() = default;
    
    virtual double medium(const Vector&) const = 0;
    
    virtual double pressure_sol(const Vector&,
                                const double) const = 0;

    virtual Vector velocity_sol(const Vector&,
                                const double) const = 0;
    
    virtual double init_pressure(const Vector&) const = 0;

    virtual Vector init_velocity(const Vector&) const = 0;
    
    virtual double bdry_pressure(const Vector&,
                                 const double) const = 0;
    
    virtual Vector bdry_velocity(const Vector&,
                                 const double) const = 0;

    virtual double source(const Vector&,
                          const double) const = 0;
    
    virtual void set_bdry_dirichlet(Array<int>&) const = 0;
    
    virtual int get_dim() const = 0;
};


// class WaveO1TestCase, base template
template<int ProblemType>
class WaveO1TestCase;


// Case: UnitSquare_Test1
// Smooth solution
// Homogeneous Dirichlet boundary
// Homogeneous source
template<>
class WaveO1TestCase <UnitSquare_Test1>
        : public WaveO1TestCases
{
public:
    explicit WaveO1TestCase (const nlohmann::json& config)
        : m_config(config), m_dim(2) {
    }

    double pressure_sol(const Vector&,
                        const double) const override;

    Vector velocity_sol(const Vector&,
                        const double) const override;

    double source(const Vector&,
                  const double) const override;

    void set_bdry_dirichlet(Array<int>&) const override;

    double medium(const Vector& x) const override {
        return m_c;
    }

    double init_pressure(const Vector& x) const override {
        return pressure_sol(x,0);
    }

    Vector init_velocity(const Vector& x) const override {
        return velocity_sol(x,0);
    }

    double bdry_pressure(const Vector& x,
                         const double t) const override {
        return pressure_sol(x, t);
    }

    Vector bdry_velocity(const Vector& x,
                         const double t) const override {
        return velocity_sol(x, t);
    }

    inline int get_dim() const override {
        return m_dim;
    }

private:
    const nlohmann::json& m_config;
    int m_dim;

    double m_c = 1;
};


// Case: UnitSquare_Test2
// Smooth solution
// Homogeneous Dirichlet boundary
// Non-homogeneous source
template<>
class WaveO1TestCase <UnitSquare_Test2>
        : public WaveO1TestCases
{
public:
    explicit WaveO1TestCase (const nlohmann::json& config)
        : m_config(config), m_dim(2) {
    }

    double pressure_sol(const Vector&,
                        const double) const override;

    Vector velocity_sol(const Vector&,
                        const double) const override;

    double source(const Vector&,
                  const double) const override;

    void set_bdry_dirichlet(Array<int>&) const override;

    double medium(const Vector& x) const override {
        return m_c;
    }

    double init_pressure(const Vector& x) const override {
        return pressure_sol(x,0);
    }

    Vector init_velocity(const Vector& x) const override {
        return velocity_sol(x,0);
    }

    double bdry_pressure(const Vector& x,
                         const double t) const override {
        return pressure_sol(x, t);
    }

    Vector bdry_velocity(const Vector& x,
                         const double t) const override {
        return velocity_sol(x, t);
    }

    inline int get_dim() const override {
        return m_dim;
    }

private:
    const nlohmann::json& m_config;
    int m_dim;

    double m_c = 1;
};


// Case: UnitSquare_Test3
// Smooth solution
// Non-homogeneous Dirichlet boundary
// Homogeneous source
template<>
class WaveO1TestCase <UnitSquare_Test3>
        : public WaveO1TestCases
{
public:
    explicit WaveO1TestCase (const nlohmann::json& config)
        : m_config(config), m_dim(2) {
    }

    double pressure_sol(const Vector&,
                        const double) const override;

    Vector velocity_sol(const Vector&,
                        const double) const override;

    double source(const Vector&,
                  const double) const override;

    void set_bdry_dirichlet(Array<int>&) const override;

    double medium(const Vector& x) const override {
        return m_c;
    }

    double init_pressure(const Vector& x) const override {
        return pressure_sol(x,0);
    }

    Vector init_velocity(const Vector& x) const override {
        return velocity_sol(x,0);
    }

    double bdry_pressure(const Vector& x,
                         const double t) const override {
        return pressure_sol(x, t);
    }

    Vector bdry_velocity(const Vector& x,
                         const double t) const override {
        return velocity_sol(x, t);
    }

    inline int get_dim() const override {
        return m_dim;
    }

private:
    const nlohmann::json& m_config;
    int m_dim;

    double m_c = 1;
};

// Case: UnitSquare_Test4
// Smooth solution
// Homogeneous Dirichlet boundary
// Non-homogeneous source
template<>
class WaveO1TestCase <UnitSquare_Test4>
        : public WaveO1TestCases
{
public:
    explicit WaveO1TestCase (const nlohmann::json& config)
        : m_config(config), m_dim(2) {
    }

    double pressure_sol(const Vector&,
                        const double) const override;

    Vector velocity_sol(const Vector&,
                        const double) const override;

    double source(const Vector&,
                  const double) const override;

    void set_bdry_dirichlet(Array<int>&) const override;

    double medium(const Vector& x) const override {
        return m_c;
    }

    double init_pressure(const Vector& x) const override {
        return pressure_sol(x,0);
    }

    Vector init_velocity(const Vector& x) const override {
        return velocity_sol(x,0);
    }

    double bdry_pressure(const Vector& x,
                         const double t) const override {
        return pressure_sol(x, t);
    }

    Vector bdry_velocity(const Vector& x,
                         const double t) const override {
        return velocity_sol(x, t);
    }

    inline int get_dim() const override {
        return m_dim;
    }

private:
    const nlohmann::json& m_config;
    int m_dim;

    double m_c = 1;
};

// Case: GammaShaped_Test1
// Singular solution
// Non-homogeneous Neumann boundary
// Non-homogeneous source
template<>
class WaveO1TestCase <GammaShaped_Test1>
        : public WaveO1TestCases
{
public:
    explicit WaveO1TestCase (const nlohmann::json& config)
        : m_config(config), m_dim(2) {
    }

    double pressure_sol(const Vector&,
                        const double) const override;

    Vector velocity_sol(const Vector&,
                        const double) const override;

    double source(const Vector&,
                  const double) const override;

    void set_bdry_dirichlet(Array<int>&) const override;

    double medium(const Vector& x) const override {
        return m_c;
    }

    double init_pressure(const Vector& x) const override {
        return pressure_sol(x,0);
    }

    Vector init_velocity(const Vector& x) const override {
        return velocity_sol(x,0);
    }

    double bdry_pressure(const Vector& x,
                         const double t) const override {
        return pressure_sol(x, t);
    }

    Vector bdry_velocity(const Vector& x,
                         const double t) const override {
        return velocity_sol(x, t);
    }

    inline int get_dim() const override {
        return m_dim;
    }

private:
    inline double radius(const Vector& x) const;
    inline double polar_angle(const Vector& x) const;

private:
    const nlohmann::json& m_config;
    int m_dim;

    double m_c = 1;
    double m_gamma = 2./3.;
};


// Case: LShaped_Test1
// Singular solution
// Homogeneous Dirichlet boundary
// Non-homogeneous source
template<>
class WaveO1TestCase <LShaped_Test1>
        : public WaveO1TestCases
{
public:
    explicit WaveO1TestCase (const nlohmann::json& config)
        : m_config(config), m_dim(2) {
    }

    double pressure_sol(const Vector&,
                        const double) const override;

    Vector velocity_sol(const Vector&,
                        const double) const override;

    double source(const Vector&,
                  const double) const override;

    void set_bdry_dirichlet(Array<int>&) const override;

    double medium(const Vector& x) const override {
        return m_c;
    }

    double init_pressure(const Vector& x) const override {
        return pressure_sol(x,0);
    }

    Vector init_velocity(const Vector& x) const override {
        return velocity_sol(x,0);
    }

    double bdry_pressure(const Vector& x,
                         const double t) const override {
        return pressure_sol(x, t);
    }

    Vector bdry_velocity(const Vector& x,
                         const double t) const override {
        return velocity_sol(x, t);
    }

    inline int get_dim() const override {
        return m_dim;
    }

private:
    inline double radius(const Vector& x) const;
    inline double polar_angle(const Vector& x) const;

private:
    const nlohmann::json& m_config;
    int m_dim;

    double m_c = 1;
    double m_gamma = 2./3.;
};

// Case: LShaped_Test2
// Singular solution
// Homogeneous Dirichlet boundary
// Non-homogeneous source
template<>
class WaveO1TestCase <LShaped_Test2>
        : public WaveO1TestCases
{
public:
    explicit WaveO1TestCase (const nlohmann::json& config)
        : m_config(config), m_dim(2) {
    }

    double pressure_sol(const Vector&,
                        const double) const override;

    Vector velocity_sol(const Vector&,
                        const double) const override;

    double source(const Vector&,
                  const double) const override;

    void set_bdry_dirichlet(Array<int>&) const override;

    double medium(const Vector& x) const override {
        return m_c;
    }

    double init_pressure(const Vector& x) const override {
        return pressure_sol(x,0);
    }

    Vector init_velocity(const Vector& x) const override {
        return velocity_sol(x,0);
    }

    double bdry_pressure(const Vector& x,
                         const double t) const override {
        return pressure_sol(x, t);
    }

    Vector bdry_velocity(const Vector& x,
                         const double t) const override {
        return velocity_sol(x, t);
    }

    inline int get_dim() const override {
        return m_dim;
    }

private:
    inline double radius(const Vector& x) const;
    inline double polar_angle(const Vector& x) const;

private:
    const nlohmann::json& m_config;
    int m_dim;

    double m_c = 1;
    double m_gamma = 2./3.;
};

/*
// Case: LShaped_Test3
// Singular solution
// Non-homogeneous Dirichlet boundary
// Non-homogeneous source
template<>
class WaveO1TestCase <LShaped_Test3>
        : public WaveO1TestCases
{
public:
    explicit WaveO1TestCase (const nlohmann::json& config)
        : m_config(config), m_dim(2) {
    }

    double pressure_sol(const Vector&,
                        const double) const override;

    Vector velocity_sol(const Vector&,
                        const double) const override;

    double source(const Vector&,
                  const double) const override;

    void set_bdry_dirichlet(Array<int>&) const override;

    double medium(const Vector& x) const override {
        return m_c;
    }

    double init_pressure(const Vector& x) const override {
        return pressure_sol(x,0);
    }

    Vector init_velocity(const Vector& x) const override {
        return velocity_sol(x,0);
    }

    double bdry_pressure(const Vector& x,
                         const double t) const override {
        return pressure_sol(x, t);
    }

    Vector bdry_velocity(const Vector& x,
                         const double t) const override {
        return velocity_sol(x, t);
    }

    inline int get_dim() const override {
        return m_dim;
    }

private:
    inline double radius(const Vector& x) const;
    inline double polar_angle(const Vector& x) const;

private:
    const nlohmann::json& m_config;
    int m_dim;

    double m_c = 1;
    double m_gamma = 2./3.;
};

// Case: LShaped_Test4
// Singular solution
// Non-homogeneous Dirichlet boundary
// Non-homogeneous source
template<>
class WaveO1TestCase <LShaped_Test4>
        : public WaveO1TestCases
{
public:
    explicit WaveO1TestCase (const nlohmann::json& config)
        : m_config(config), m_dim(2) {
    }

    double pressure_sol(const Vector&,
                        const double) const override;

    Vector velocity_sol(const Vector&,
                        const double) const override;

    double source(const Vector&,
                  const double) const override;

    void set_bdry_dirichlet(Array<int>&) const override;

    double medium(const Vector& x) const override {
        return m_c;
    }

    double init_pressure(const Vector& x) const override {
        return pressure_sol(x,0);
    }

    Vector init_velocity(const Vector& x) const override {
        return velocity_sol(x,0);
    }

    double bdry_pressure(const Vector& x,
                         const double t) const override {
        return pressure_sol(x, t);
    }

    Vector bdry_velocity(const Vector& x,
                         const double t) const override {
        return velocity_sol(x, t);
    }

    inline int get_dim() const override {
        return m_dim;
    }

private:
    inline double radius(const Vector& x) const;
    inline double polar_angle(const Vector& x) const;

private:
    const nlohmann::json& m_config;
    int m_dim;

    double m_c = 1;
    double m_gamma = 2./3.;
};


// Case: LShaped_Test5
// Smooth solution
// Homogeneous Dirichlet boundary
// Non-homogeneous source
template<>
class WaveO1TestCase <LShaped_Test5>
        : public WaveO1TestCases
{
public:
    explicit WaveO1TestCase (const nlohmann::json& config)
        : m_config(config), m_dim(2) {
    }

    double pressure_sol(const Vector&,
                        const double) const override;

    Vector velocity_sol(const Vector&,
                        const double) const override;

    double source(const Vector&,
                  const double) const override;

    void set_bdry_dirichlet(Array<int>&) const override;

    double medium(const Vector& x) const override {
        return m_c;
    }

    double init_pressure(const Vector& x) const override {
        return pressure_sol(x,0);
    }

    Vector init_velocity(const Vector& x) const override {
        return velocity_sol(x,0);
    }

    double bdry_pressure(const Vector& x,
                         const double t) const override {
        return pressure_sol(x, t);
    }

    Vector bdry_velocity(const Vector& x,
                         const double t) const override {
        return velocity_sol(x, t);
    }

    inline int get_dim() const override {
        return m_dim;
    }

private:
    const nlohmann::json& m_config;
    int m_dim;

    double m_c = 1;
    double m_gamma = 2./3.;
};
*/

// Case: SquareTwoPiece_Test1
template<>
class WaveO1TestCase <SquareTwoPiece_Test1>
        : public WaveO1TestCases
{
public:
    explicit WaveO1TestCase (const nlohmann::json& config)
        : m_config(config), m_dim(2)
    {
        m_x0.SetSize(2);
        m_x0 = 1.; // Gaussian wave center
        m_lambda = 0.01; // Gaussian wave width
    }

    double init_pressure(const Vector&) const override;

    Vector init_velocity(const Vector&) const override;

    double source(const Vector&,
                  const double) const override;

    void set_bdry_dirichlet(Array<int>&) const override;

    double medium(const Vector& x) const override;

    double pressure_sol(const Vector&, const double) const override {
        return 0;
    }

    Vector velocity_sol(const Vector&, const double) const override {
        return Vector(2);
    }

    double bdry_pressure(const Vector& x,
                         const double t) const override {
        return 0;
    }

    Vector bdry_velocity(const Vector& x,
                         const double t) const override {
        Vector v(2);
        v = 0.;
        return v;
    }

    inline int get_dim() const override {
        return m_dim;
    }

private:
    const nlohmann::json& m_config;
    int m_dim;

    double m_xs = 1.2;
    Vector m_x0;
    double m_lambda;
};


// Case: SquareTwoPiece_Test2
template<>
class WaveO1TestCase <SquareTwoPiece_Test2>
        : public WaveO1TestCases
{
public:
    explicit WaveO1TestCase (const nlohmann::json& config)
        : m_config(config), m_dim(2)
    {
        m_x0.SetSize(2);
        m_x0(0) = 1.; // Gaussian wave center
        m_x0(1) = 1.125;
        m_lambda = 0.01; // Gaussian wave width
    }

    double init_pressure(const Vector&) const override;

    Vector init_velocity(const Vector&) const override;

    double source(const Vector&,
                  const double) const override;

    void set_bdry_dirichlet(Array<int>&) const override;

    double medium(const Vector& x) const override;

    double pressure_sol(const Vector& x,
                        const double) const override {
        return init_pressure(x);
    }

    Vector velocity_sol(const Vector& x,
                        const double) const override {
        return init_velocity(x);
    }

    double bdry_pressure(const Vector& x,
                         const double t) const override {
        return 0;
    }

    Vector bdry_velocity(const Vector& x,
                         const double t) const override {
        Vector v(2);
        v = 0.;
        return v;
    }

    inline int get_dim() const override {
        return m_dim;
    }

private:
    const nlohmann::json& m_config;
    int m_dim;

    double m_xs = 1.2;
    Vector m_x0;
    double m_lambda;
};

#endif /// WAVEO1_TEST_CASES_HPP
