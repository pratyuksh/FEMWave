#include "../../include/waveO2/discretisation.hpp"
#include "../../include/waveO2/assembly.hpp"
#include "../include/mymfem/utilities.hpp"


//! Constructor
WaveO2XtFEM
:: WaveO2XtFEM (const nlohmann::json& config,
                std::shared_ptr<WaveO2TestCases>& testCase)
    : m_config (config),
      m_testCase(testCase)
{
    m_endTime = config["end_time"];
    m_tDeg = config["deg_t"];
    m_xDeg = config["deg_x"];

    m_sqMed = std::make_shared<WaveO2SqMediumCoeff>
            (m_testCase);
}

//! Destructor
WaveO2XtFEM :: ~WaveO2XtFEM()
{
    if (m_tFec) { delete m_tFec; }
    if (m_tFespace) { delete m_tFespace; }

    if (m_xFec) { delete m_xFec; }
    if (m_xFespace) { delete m_xFespace; }

    if (m_waveMat) { delete m_waveMat; }
}

//! Sets-up the FE spaces and boundary conditions
void WaveO2XtFEM :: set(std::shared_ptr<Mesh>& tMesh,
                        std::shared_ptr<Mesh>& xMesh)
{
    m_tFec = new H1_FECollection(m_tDeg, 1,
                                 BasisType::GaussLobatto);
    m_tFespace = new FiniteElementSpace(tMesh.get(), m_tFec);

    m_xndim = xMesh->Dimension();
    m_xFec  = new H1_FECollection(m_xDeg, m_xndim,
                                  BasisType::GaussLobatto);
    m_xFespace = new FiniteElementSpace(xMesh.get(), m_xFec);

    int tdimV = m_tFespace->GetTrueVSize();
    int xdimV = m_xFespace->GetTrueVSize();
#ifdef MYVERBOSE
    std::cout << "\n\tNumber of degrees of freedom "
                 "in tV_space: "
         << tdimV << std::endl;
    std::cout << "\tNumber of degrees of freedom "
                 "in xVspace: "
         << xdimV << std::endl;
#endif
    // Mark boundary dofs for xV_space
    if (xMesh->bdr_attributes.Size())
    {
        m_xEss_bdr_marker.SetSize
                (xMesh->bdr_attributes.Max());
        m_testCase->set_bdry_dirichlet(m_xEss_bdr_marker);

        m_xFespace->GetEssentialTrueDofs
                (m_xEss_bdr_marker, m_xEss_tdof_list);

        int xNEssDofs = m_xEss_tdof_list.Size();
        m_ess_tdof_list.SetSize(xNEssDofs*tdimV);
        for (int j=0; j<tdimV; j++)
            for (int i=0; i<xNEssDofs; i++)
                m_ess_tdof_list[i + j*xNEssDofs]
                        = j*xdimV + m_xEss_tdof_list[i];
    }
}

//! Assembles the matrices
void WaveO2XtFEM :: assemble_system()
{
    /// assemble matrices in space

    // xMass form, xVspace
    BilinearForm *xMass_form
            = new BilinearForm(m_xFespace);
    xMass_form->AddDomainIntegrator(new MassIntegrator);
    xMass_form->Assemble();
    xMass_form->Finalize();
    m_xMass = xMass_form->LoseMat();
    delete xMass_form;

    // xStiff form, xVspace
    BilinearForm *xStiff_form
            = new BilinearForm(m_xFespace);
    xStiff_form->AddDomainIntegrator
            (new DiffusionIntegrator(*m_sqMed));
    xStiff_form->Assemble();
    xStiff_form->Finalize();
    m_xStiff = xStiff_form->LoseMat();
    delete xStiff_form;
}

// End of file
