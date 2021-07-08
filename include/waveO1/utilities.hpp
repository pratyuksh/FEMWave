#ifndef WAVEO1_UTILITIES_HPP
#define WAVEO1_UTILITIES_HPP

#include "mfem.hpp"
using namespace mfem;


void build_xSol_FG(const Vector& uXtSol,
                   const Vector& tShape,
                   Array<int> tVdofs,
                   Vector& uSol);


namespace mymfem {

// DG error
class AssembleDgError
{
public:
    explicit AssembleDgError(Array<int>& xEss_bdr_marker,
                             Array<int>& xNat_bdr_marker,
                             int stabParamsType=1)
        : m_xEss_bdr_marker(xEss_bdr_marker),
          m_xNat_bdr_marker(xNat_bdr_marker),
          m_stabParamsType(stabParamsType) {}

    explicit AssembleDgError(IntegrationRule* ir,
                             Array<int>& xEss_bdr_marker,
                             Array<int>& xNat_bdr_marker,
                             int stabParamsType=1)
        : m_ir(ir),
          m_xEss_bdr_marker(xEss_bdr_marker),
          m_xNat_bdr_marker(xNat_bdr_marker),
          m_stabParamsType(stabParamsType) {}

    // time-like face
    double Ftime(GridFunction* u,
                 Coefficient* uE_coeff);

    double Ftime(GridFunction* u,
                 VectorCoefficient* uE_coeff);

    // space-like face
    double Fspace(GridFunction* u, Coefficient* uE,
                  Coefficient* medium=nullptr);
    double Fspace(GridFunction* u1, GridFunction* u2,
                  Coefficient* uE,
                  Coefficient* medium=nullptr);

    double Fspace(GridFunction* u, VectorCoefficient* uE);
    double Fspace(GridFunction* u1, GridFunction* u2,
                  VectorCoefficient* uE);

    void SetIntegrationRule (const IntegrationRule* ir)
    { m_ir = ir; }

private:
    const IntegrationRule *m_ir = nullptr;
    Array<int> m_xEss_bdr_marker;
    Array<int> m_xNat_bdr_marker;
    int m_stabParamsType;
};


// DG jumps
class AssembleDgJumps
{
public:
    explicit AssembleDgJumps(Array<int>& xEss_bdr_marker,
                             Array<int>& xNat_bdr_marker,
                             int stabParamsType=1)
        : m_xEss_bdr_marker(xEss_bdr_marker),
          m_xNat_bdr_marker(xNat_bdr_marker),
          m_stabParamsType(stabParamsType) {}

    explicit AssembleDgJumps(IntegrationRule* ir,
                             Array<int>& xEss_bdr_marker,
                             Array<int>& xNat_bdr_marker,
                             int stabParamsType=1)
        : m_ir(ir),
          m_xEss_bdr_marker(xEss_bdr_marker),
          m_xNat_bdr_marker(xNat_bdr_marker),
          m_stabParamsType(stabParamsType) {}

    // time-like face
    double FtimeScalar(GridFunction* u);
    double FtimeVector(GridFunction* u);

    // space-like face
    double FspaceScalar(Mesh*,
                        Coefficient* u,
                        Coefficient* medium=nullptr);
    double FspaceScalar(GridFunction* u,
                        Coefficient* medium=nullptr);
    double FspaceScalar(GridFunction* u1, GridFunction* u2,
                        Coefficient* medium=nullptr);

    double FspaceVector(Mesh*,
                        VectorCoefficient* u);
    double FspaceVector(GridFunction* u);
    double FspaceVector(GridFunction* u1, GridFunction* u2);

    void SetIntegrationRule (const IntegrationRule* ir)
    { m_ir = ir; }

private:
    const IntegrationRule *m_ir = nullptr;
    Array<int> m_xEss_bdr_marker;
    Array<int> m_xNat_bdr_marker;
    int m_stabParamsType;
};


// DG^{+} smooth error
class AssembleDgPlusSmoothError
{
public:
    explicit AssembleDgPlusSmoothError
    (Array<int>& xEss_bdr_marker,
     Array<int>& xNat_bdr_marker,
     int stabParamsType=1)
        : m_xEss_bdr_marker(xEss_bdr_marker),
          m_xNat_bdr_marker(xNat_bdr_marker),
          m_stabParamsType(stabParamsType) {}

    explicit AssembleDgPlusSmoothError
    (IntegrationRule* ir,
     Array<int>& xEss_bdr_marker,
     Array<int>& xNat_bdr_marker,
     int stabParamsType=1)
        : m_ir(ir),
          m_xEss_bdr_marker(xEss_bdr_marker),
          m_xNat_bdr_marker(xNat_bdr_marker),
          m_stabParamsType(stabParamsType) {}

    // time-like face
    double Ftime(GridFunction* u,
                 Coefficient* uE_coeff);

    double Ftime(GridFunction* u,
                 VectorCoefficient* uE_coeff);

    // space-like face
    double Fspace(GridFunction* u, Coefficient* uE);

    double Fspace(GridFunction* u, VectorCoefficient* uE);

    void SetIntegrationRule (const IntegrationRule* ir)
    { m_ir = ir; }

private:
    const IntegrationRule *m_ir = nullptr;
    Array<int> m_xEss_bdr_marker;
    Array<int> m_xNat_bdr_marker;
    int m_stabParamsType;
};

}

#endif /// WAVEO1_UTILITIES_HPP
