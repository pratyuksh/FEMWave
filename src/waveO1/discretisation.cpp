#include "../../include/waveO1/discretisation.hpp"
#include "../../include/waveO1/assembly.hpp"
#include "../include/mymfem/mymixedbf.hpp"
#include "../include/mymfem/utilities.hpp"

#include <fstream>


// Constructor
WaveO1XtDG
:: WaveO1XtDG (const nlohmann::json& config,
             std::shared_ptr<WaveO1TestCases>& testCase)
    : m_config (config),
      m_testCase(testCase)
{
    m_stabParamsType = 1;
    if (config.contains("stab_parameters_type")) {
        m_stabParamsType = config["stab_parameters_type"];
    }

    m_xDeg1 = config["deg1_x"];
    m_xDeg2 = config["deg2_x"];

    m_tDeg = config["deg_t"];
    m_endTime = config["end_time"];

    m_tMesh0 = std::make_unique<Mesh>(1, 1);

    m_invSqMed = std::make_shared<WaveO1InvSqMediumCoeff>
            (m_testCase);
}

// Destructor
WaveO1XtDG :: ~WaveO1XtDG()
{
    if (m_tFec) { delete m_tFec; }
    if (m_tFespace0) { delete m_tFespace0; }
    if (m_tFespace) { delete m_tFespace; }
    
    for (int i=0; i<m_xFecs.Size(); i++)
    {
        if (m_xFecs[i]) { delete m_xFecs[i]; }
        if (m_xFespaces[i]) { delete m_xFespaces[i]; }
    }

    if (m_waveOp) { delete m_waveOp; }
    if (m_waveMat) { delete m_waveMat; }
}

// Set-up the finite element spaces and boundary conditions
void WaveO1XtDG :: set(std::shared_ptr<Mesh>& tMesh,
                     std::shared_ptr<Mesh>& xMesh)
{
    m_tFec = new L2_FECollection(m_tDeg, 1,
                                 BasisType::GaussLegendre);
    m_tFespace0
            = new FiniteElementSpace(m_tMesh0.get(), m_tFec);
    m_tFespace
            = new FiniteElementSpace(tMesh.get(), m_tFec);
    
    m_xndim = xMesh->Dimension();
    FiniteElementCollection *xL2_coll1
            = new L2_FECollection(m_xDeg1, m_xndim,
                                  BasisType::GaussLegendre);
    auto *xW1_space
            = new FiniteElementSpace(xMesh.get(), xL2_coll1);
    
    FiniteElementCollection *xL2_coll2
            = new L2_FECollection(m_xDeg2, m_xndim,
                                  BasisType::GaussLegendre);
    auto *xW2_space
            = new FiniteElementSpace(xMesh.get(), xL2_coll2, m_xndim);

    int tdimW = m_tFespace->GetTrueVSize();
    int xdimW1 = xW1_space->GetTrueVSize();
    int xdimW2 = xW2_space->GetTrueVSize();
#ifdef MYVERBOSE
    std::cout << "\nNumber of degrees of freedom "
                 "in tW_space: "
              << tdimW << std::endl;
    std::cout << "Number of degrees of freedom "
                 "in xW1_space: "
              << xdimW1 << std::endl;
    std::cout << "Number of degrees of freedom"
                 " in xW2_space: "
              << xdimW2 << std::endl;
#endif
    m_xFecs.Append(xL2_coll1);
    m_xFecs.Append(xL2_coll2);
    
    m_xFespaces.Append(xW1_space);
    m_xFespaces.Append(xW2_space);

    /// Define the two BlockStructure of the problem.
    m_block_offsets.SetSize(3);
    m_block_offsets[0] = 0;
    m_block_offsets[1] = xdimW1*tdimW;
    m_block_offsets[2] = xdimW2*tdimW;
    m_block_offsets.PartialSum();
#ifdef MYVERBOSE
    std::cout << "Block offsets:" << "\t"
              << m_block_offsets[0] << "\t"
              << m_block_offsets[1] << "\t"
              << m_block_offsets[2] << std::endl;
#endif
    int tdimW_slab = m_tFespace0->GetTrueVSize();
    m_slab_block_offsets.SetSize(3);
    m_slab_block_offsets[0] = 0;
    m_slab_block_offsets[1] = xdimW1*tdimW_slab;
    m_slab_block_offsets[2] = xdimW2*tdimW_slab;
    m_slab_block_offsets.PartialSum();
    
    // Mark boundary attributes
    if (xMesh->bdr_attributes.Size())
    {
        m_xEss_bdr_marker.SetSize
                (xMesh->bdr_attributes.Max());
        m_testCase->set_bdry_dirichlet(m_xEss_bdr_marker);
        
        m_xNat_bdr_marker.SetSize
                (xMesh->bdr_attributes.Max());
        for (int i=0; i<m_xNat_bdr_marker.Size(); i++) {
            m_xNat_bdr_marker[i]
                    = (m_xEss_bdr_marker[i] > 0 ? 0 : 1);
        }
    }
}

/// Assemble matrices
void WaveO1XtDG :: assemble_system()
{   
    /// assemble matrices in time

    // tMass form
    double ht = m_endTime/m_tFespace->GetNE();
    BilinearForm *tMass0_form
            = new BilinearForm(m_tFespace0);
    tMass0_form->AddDomainIntegrator(new MassIntegrator);
    tMass0_form->Assemble();
    tMass0_form->Finalize();
    m_slab_tMass = tMass0_form->LoseMat();
    delete tMass0_form;
    (*m_slab_tMass) *= ht;

    // tGrad form
    BilinearForm *tGrad_form
            = new BilinearForm(m_tFespace0);
    tGrad_form->AddDomainIntegrator
            (new mymfem::TimeGradientIntegrator);
    tGrad_form->Assemble();
    tGrad_form->Finalize();
    m_slab_tGrad = tGrad_form->LoseMat();
    delete tGrad_form;

    // tFluxes
    BilinearForm *tFluxp_form
            = new BilinearForm(m_tFespace0);
    tFluxp_form->AddDomainIntegrator
            (new mymfem::TimeFluxPlusIntegrator);
    tFluxp_form->Assemble();
    tFluxp_form->Finalize();
    m_slab_tFluxp = tFluxp_form->LoseMat();
    delete tFluxp_form;

    BilinearForm *tFluxm_form
            = new BilinearForm(m_tFespace0);
    tFluxm_form->AddDomainIntegrator
            (new mymfem::TimeFluxMinusIntegrator);
    tFluxm_form->Assemble();
    tFluxm_form->Finalize();
    m_slab_tFluxm = tFluxm_form->LoseMat();
    delete tFluxm_form;

    /// assemble matrices in space

    // global mesh-size
    double dbuf, hxMax;
    get_xMesh()->GetCharacteristics(dbuf, hxMax, dbuf, dbuf);

    // xMass0 form, xW1_space
    BilinearForm *xMass0_form
            = new BilinearForm(m_xFespaces[0]);
    xMass0_form->AddDomainIntegrator(new MassIntegrator
                                     (*m_invSqMed));
    xMass0_form->Assemble();
    xMass0_form->Finalize();
    m_xMass0 = xMass0_form->LoseMat();
    delete xMass0_form;

    // xMass1 form, xW2_space
    BilinearForm *xMass1_form
            = new BilinearForm(m_xFespaces[1]);
    xMass1_form->AddDomainIntegrator
            (new VectorMassIntegrator);
    xMass1_form->Assemble();
    xMass1_form->Finalize();
    m_xMass1 = xMass1_form->LoseMat();
    delete xMass1_form;

    // xDiv form
    MixedBilinearForm *xDiv_form
            = new MixedBilinearForm(m_xFespaces[1],
                                    m_xFespaces[0]);
    xDiv_form->AddDomainIntegrator
            (new VectorDivergenceIntegrator);
    xDiv_form->Assemble();
    xDiv_form->Finalize();
    m_xDiv = xDiv_form->LoseMat();
    delete xDiv_form;

    // xGrad form
    MixedBilinearForm *xGrad_form
            = new MixedBilinearForm(m_xFespaces[0],
                                    m_xFespaces[1]);
    xGrad_form->AddDomainIntegrator
            (new mymfem::VectorGradientIntegrator);
    xGrad_form->Assemble();
    xGrad_form->Finalize();
    m_xGrad = xGrad_form->LoseMat();
    delete xGrad_form;

    // xFlux00 form
    BilinearForm *xFlux00_form
            = new BilinearForm(m_xFespaces[0]);
    xFlux00_form->AddInteriorFaceIntegrator
            (new mymfem::FluxBlock00Integrator(hxMax, m_stabParamsType));
    xFlux00_form->AddBdrFaceIntegrator
            (new mymfem::FluxBlock00Integrator(hxMax, m_stabParamsType),
             m_xEss_bdr_marker);
    xFlux00_form->Assemble();
    xFlux00_form->Finalize();
    m_xFlux00 = xFlux00_form->LoseMat();
    delete xFlux00_form;

    // xFlux11 form
    BilinearForm *xFlux11_form
            = new BilinearForm(m_xFespaces[1]);
    xFlux11_form->AddInteriorFaceIntegrator
            (new mymfem::FluxBlock11Integrator(hxMax, m_stabParamsType));
    xFlux11_form->AddBdrFaceIntegrator
            (new mymfem::FluxBlock11Integrator(hxMax, m_stabParamsType),
             m_xNat_bdr_marker);
    xFlux11_form->Assemble();
    xFlux11_form->Finalize();
    m_xFlux11 = xFlux11_form->LoseMat();
    delete xFlux11_form;

    // xFlux10 form
    mymfem::MyMixedBilinearForm *xFlux10_form
            = new mymfem::MyMixedBilinearForm
            (m_xFespaces[0], m_xFespaces[1]);
    xFlux10_form->MyAddInteriorFaceIntegrator
            (new mymfem::FluxBlock10Integrator);
    xFlux10_form->MyAddBoundaryFaceIntegrator
            (new mymfem::DirichletBlock10Integrator,
             m_xEss_bdr_marker);
    xFlux10_form->MyAssemble();
    xFlux10_form->Finalize();
    m_xFlux10 = xFlux10_form->LoseMat();
    xFlux10_form->MyClear();
    delete xFlux10_form;

    // xFlux01 form
    mymfem::MyMixedBilinearForm *xFlux01_form
            = new mymfem::MyMixedBilinearForm
            (m_xFespaces[1], m_xFespaces[0]);
    xFlux01_form->MyAddInteriorFaceIntegrator
            (new mymfem::FluxBlock01Integrator);
    xFlux01_form->MyAddBoundaryFaceIntegrator
            (new mymfem::NeumannBlock01Integrator,
             m_xNat_bdr_marker);
    xFlux01_form->MyAssemble();
    xFlux01_form->Finalize();
    m_xFlux01 = xFlux01_form->LoseMat();
    xFlux01_form->MyClear();
    delete xFlux01_form;

    /// compute blocks
    /*SparseMatrix *tmp1, *tmp2;

    // blocks at row 0, column 0 and row 1, column 1
    tmp1 = Add(*m_slab_tGrad, *m_slab_tFluxp);
    m_block00 = OuterProduct(*tmp1, *m_xMass1);
    m_block11 = OuterProduct(*tmp1, *m_xMass2);
    delete tmp1;
    // block at row 0, column 0
    tmp1 = Add(*m_xFlux00, *m_xDirichlet00);
    tmp2 = OuterProduct(*m_slab_tMass, *tmp1);
    m_block00->Add(1, *tmp2);
    delete tmp1;
    delete tmp2;
    // block at row 1, column 1
    tmp1 = Add(*m_xFlux11, *m_xNeumann11);
    tmp2 = OuterProduct(*m_slab_tMass, *tmp1);
    m_block11->Add(1, *tmp2);
    delete tmp1;
    delete tmp2;

    // block at row 0, column 1
    tmp1 = Add(*m_xDiv, *m_xFlux01);
    tmp1->Add(1, *m_xNeumann01);
    m_block01 = OuterProduct(*m_slab_tMass, *tmp1);
    delete tmp1;

    // block at row 1, column 0
    tmp1 = Add(*m_xGrad, *m_xFlux10);
    tmp1->Add(1, *m_xDirichlet10);
    m_block10 = OuterProduct(*m_slab_tMass, *tmp1);
    delete tmp1;*/

    /// compute blocks
    SparseMatrix *tmp1;
    SparseMatrix *bufBlock00, *bufBlock11;

    // blocks at row 0, column 0 and row 1, column 1
    tmp1 = Add(*m_slab_tGrad, *m_slab_tFluxp);
    bufBlock00 = OuterProduct(*tmp1, *m_xMass0);
    bufBlock11 = OuterProduct(*tmp1, *m_xMass1);
    delete tmp1;
    // block at row 0, column 0
    tmp1 = OuterProduct(*m_slab_tMass, *m_xFlux00);
    m_block00 = Add(*tmp1, *bufBlock00);
    delete bufBlock00;
    delete tmp1;
    // block at row 1, column 1
    tmp1 = OuterProduct(*m_slab_tMass, *m_xFlux11);
    m_block11 = Add(*tmp1, *bufBlock11);
    delete bufBlock11;
    delete tmp1;

    // block at row 0, column 1
    tmp1 = Add(*m_xDiv, *m_xFlux01);
    m_block01 = OuterProduct(*m_slab_tMass, *tmp1);
    delete tmp1;

    // block at row 1, column 0
    tmp1 = Add(*m_xGrad, *m_xFlux10);
    m_block10 = OuterProduct(*m_slab_tMass, *tmp1);
    delete tmp1;

    /// rhs matrices
    m_slab_rhsMat0 = OuterProduct(*m_slab_tFluxm, *m_xMass0);
    m_slab_rhsMat1 = OuterProduct(*m_slab_tFluxm, *m_xMass1);
}

// Build waveOp
void WaveO1XtDG :: build_system_op()
{
    if (m_waveOp) {
        delete m_waveOp;
        m_waveOp = nullptr;
    }
    m_waveOp = new BlockOperator(m_slab_block_offsets);
    m_waveOp->SetBlock(0,0, m_block00);
    m_waveOp->SetBlock(0,1, m_block01);
    m_waveOp->SetBlock(1,0, m_block10);
    m_waveOp->SetBlock(1,1, m_block11);
}

// Build system matrix from waveOp
void WaveO1XtDG :: build_system_matrix()
{
    if (m_waveMat) { clear(m_waveMat); }

    BlockMatrix *waveBlockMat
            = new BlockMatrix(m_slab_block_offsets);
    waveBlockMat->SetBlock(0,0, m_block00);
    waveBlockMat->SetBlock(0,1, m_block01);
    waveBlockMat->SetBlock(1,0, m_block10);
    waveBlockMat->SetBlock(1,1, m_block11);

    m_waveMat = waveBlockMat->CreateMonolithic();
    delete waveBlockMat;

    if (!m_waveMat->ColumnsAreSorted()) {
        m_waveMat->SortColumnIndices();
    }
}

// Assembles rhs at the first time-slab
void WaveO1XtDG
:: assemble_rhs_initial(BlockVector* slabB) const
{
    assemble_ics(slabB);
    assemble_source(0, slabB);
    assemble_bcs(0, slabB);
}

// Assembles rhs
void WaveO1XtDG
:: assemble_rhs(const int n, const BlockVector * slabW,
                BlockVector* slabB) const
{
    // contribution from previous time slab
    m_slab_rhsMat0->Mult(slabW->GetBlock(0),
                         slabB->GetBlock(0));
    m_slab_rhsMat1->Mult(slabW->GetBlock(1),
                         slabB->GetBlock(1));

    assemble_source(n, slabB);
    assemble_bcs(n, slabB);
}

// Assembles initial conditions
void WaveO1XtDG
:: assemble_ics(BlockVector* slabB) const
{
    int n=0;

    int xdimW1 = m_xFespaces[0]->GetTrueVSize();
    int xdimW2 = m_xFespaces[1]->GetTrueVSize();
    Vector bx1(xdimW1), bx2(xdimW2);

    const FiniteElement *tFe = nullptr;
    Vector elvec1, elvec2, tShape;
    {
        // compute elvec
        tFe = m_tFespace->GetFE(n);
        int tNdofs = tFe->GetDof();
        tShape.SetSize(tNdofs);
        elvec1.SetSize(tNdofs*xdimW1);
        elvec2.SetSize(tNdofs*xdimW2);

        IntegrationPoint ip;
        ip.Set1w(0.0, 1.0);
        tFe->CalcShape(ip, tShape);

        assemble_ics_xProjection(bx1, bx2);

        for (int j=0; j<tNdofs; j++) {
            double coeff = tShape(j);
            for (int k=0; k<xdimW1; k++) {
                elvec1(k + j*xdimW1) = coeff*bx1(k);
            }
            for (int k=0; k<xdimW2; k++) {
                elvec2(k + j*xdimW2) = coeff*bx2(k);
            }
        }
    }

    Array<int> slab_vdofs;
    m_tFespace0->GetElementVDofs (0, slab_vdofs);
    add_vector(slab_vdofs, elvec1, slabB->GetBlock(0));
    add_vector(slab_vdofs, elvec2, slabB->GetBlock(1));
}

void WaveO1XtDG
:: assemble_ics_xProjection(Vector& b1,
                            Vector& b2) const
{
    // initial pressure; FixMe for medium coefficient
    LinearForm lF1(m_xFespaces[0]);
    WaveO1InitialPressureCoeff pressure(m_testCase);
    pressure.SetTime(0);
    lF1.AddDomainIntegrator
            (new DomainLFIntegrator(pressure));
    lF1.Assemble();
    b1 = lF1.GetData();

    // initial velocity
    LinearForm lF2(m_xFespaces[1]);
    WaveO1InitialVelocityCoeff velocity(m_testCase);
    velocity.SetTime(0);
    lF2.AddDomainIntegrator
            (new VectorDomainLFIntegrator(velocity));
    lF2.Assemble();
    b2 = lF2.GetData();
}

// Assembles source
void WaveO1XtDG
:: assemble_source(const int n, BlockVector* slabB) const
{
    int xdimW1 = m_xFespaces[0]->GetTrueVSize();
    int xdimW2 = m_xFespaces[1]->GetTrueVSize();
    Vector bx1(xdimW1), bx2(xdimW2);

    ElementTransformation *trans = nullptr;
    const FiniteElement *fe = nullptr;
    Vector elvec1, elvec2, shape;
    {
        trans = m_tFespace->GetElementTransformation(n);

        // compute elvec
        fe = m_tFespace->GetFE(n);
        int tNdofs = fe->GetDof();
        shape.SetSize(tNdofs);
        elvec1.SetSize(tNdofs*xdimW1);
        elvec2.SetSize(tNdofs*xdimW2);

        int order = 2*fe->GetOrder()+2;
        const IntegrationRule *ir
                = &IntRules.Get(fe->GetGeomType(), order);

        elvec1 = 0.0;
        elvec2 = 0.0;
        for (int i = 0; i < ir->GetNPoints(); i++)
        {
            const IntegrationPoint &ip = ir->IntPoint(i);
            trans->SetIntPoint(&ip);
            fe->CalcShape(ip, shape);

            Vector t;
            trans->Transform(ip, t);
            assemble_source_xProjection(bx1, bx2, t(0));

            // elvec += w*shape*bx
            double w = ip.weight*trans->Weight();
            for (int j=0; j<tNdofs; j++) {
                double coeff = w*shape(j);
                for (int k=0; k<xdimW1; k++) {
                    elvec1(k + j*xdimW1) += coeff*bx1(k);
                }
                for (int k=0; k<xdimW2; k++) {
                    elvec2(k + j*xdimW2) += coeff*bx2(k);
                }
            }
        }
    }

    Array<int> slab_vdofs;
    m_tFespace0->GetElementVDofs (0, slab_vdofs);
    add_vector(slab_vdofs, elvec1, slabB->GetBlock(0));
    add_vector(slab_vdofs, elvec2, slabB->GetBlock(1));
}

void WaveO1XtDG
:: assemble_source_xProjection(Vector& b1,
                               Vector& b2,
                               double t) const
{
    LinearForm lForm(m_xFespaces[0]);
    WaveO1SourceCoeff source(m_testCase);
    source.SetTime(t);
    lForm.AddDomainIntegrator
            (new DomainLFIntegrator(source));
    lForm.Assemble();
    b1 = lForm.GetData();

    b2 = 0.;
}

// Assembles boundary conditions
void WaveO1XtDG
:: assemble_bcs(const int n, BlockVector* slabB) const
{
    int xdimW1 = m_xFespaces[0]->GetTrueVSize();
    int xdimW2 = m_xFespaces[1]->GetTrueVSize();
    Vector bx1(xdimW1), bx2(xdimW2);

    ElementTransformation *trans = nullptr;
    const FiniteElement *fe = nullptr;
    Vector elvec1, elvec2, shape;
    {
        trans = m_tFespace->GetElementTransformation(n);

        // compute elvec
        fe = m_tFespace->GetFE(n);
        int tNdofs = fe->GetDof();
        shape.SetSize(tNdofs);
        elvec1.SetSize(tNdofs*xdimW1);
        elvec2.SetSize(tNdofs*xdimW2);

        int order = 2*fe->GetOrder()+2;
        const IntegrationRule *ir
                = &IntRules.Get(fe->GetGeomType(), order);

        elvec1 = 0.0;
        elvec2 = 0.0;
        for (int i = 0; i < ir->GetNPoints(); i++)
        {
            const IntegrationPoint &ip = ir->IntPoint(i);
            trans->SetIntPoint(&ip);
            fe->CalcShape(ip, shape);

            Vector t;
            trans->Transform(ip, t);
            assemble_bcs_xProjection(bx1, bx2, t(0));

            // elvec += w*shape*bx
            double w = ip.weight*trans->Weight();
            for (int j=0; j<tNdofs; j++) {
                double coeff = w*shape(j);
                for (int k=0; k<xdimW1; k++) {
                    elvec1(k + j*xdimW1) += coeff*bx1(k);
                }
                for (int k=0; k<xdimW2; k++) {
                    elvec2(k + j*xdimW2) += coeff*bx2(k);
                }
            }
        }
    }

    Array<int> slab_vdofs;
    m_tFespace0->GetElementVDofs (0, slab_vdofs);
    add_vector(slab_vdofs, elvec1, slabB->GetBlock(0));
    add_vector(slab_vdofs, elvec2, slabB->GetBlock(1));
}

void WaveO1XtDG
:: assemble_bcs_xProjection (Vector& b1,
                             Vector& b2,
                             double t) const
{
    // Dirichlet boundary conditions
    WaveO1BdryPressureCoeff pressure(m_testCase);
    pressure.SetTime(t);

    LinearForm dForm1(m_xFespaces[0]);
    dForm1.AddBdrFaceIntegrator
            (new mymfem::DirichletLFBlock0Integrator
             (pressure, m_stabParamsType),
             m_xEss_bdr_marker);
    dForm1.Assemble();
    b1 = dForm1.GetData();

    LinearForm dForm2(m_xFespaces[1]);
    dForm2.AddBdrFaceIntegrator
            (new mymfem::DirichletLFBlock1Integrator(pressure, -1),
             m_xEss_bdr_marker);
    dForm2.Assemble();
    b2 = dForm2.GetData();

    // Neumann boundary conditions
    WaveO1BdryVelocityCoeff velocity(m_testCase);
    velocity.SetTime(t);

    Vector tmp(b1.Size());
    LinearForm nForm1(m_xFespaces[0]);
    nForm1.AddBdrFaceIntegrator
            (new mymfem::NeumannLFBlock0Integrator(velocity),
             m_xNat_bdr_marker);
    nForm1.Assemble();
    tmp = nForm1.GetData();
    b1.Add(-1, tmp);

    tmp.SetSize(b2.Size());
    LinearForm nForm2(m_xFespaces[1]);
    nForm2.AddBdrFaceIntegrator
            (new mymfem::NeumannLFBlock1Integrator(velocity),
             m_xNat_bdr_marker);
    nForm2.Assemble();
    tmp = nForm2.GetData();
    b2.Add(1, tmp);

    // Release memory
    tmp.Destroy();
}


// Adds elvec to b
void WaveO1XtDG :: add_vector(const Array<int> &vdofs,
                            const Vector& elvec,
                            Vector& b) const
{
    int ndofs = elvec.Size();
    int tNdofs = vdofs.Size();
    int xNdofs = ndofs/tNdofs;

    const bool use_dev
            = vdofs.UseDevice() || elvec.UseDevice();
    auto d_elvec = elvec.Read(use_dev);
    auto d_b = b.ReadWrite(use_dev);

    for(int i=0; i<tNdofs; i++)
    {
        const int j = vdofs[i];
        const int ii = i*xNdofs;
        const int jj = j*xNdofs;
        for (int k=0; k<xNdofs; k++)
            d_b[k + jj] += d_elvec[k + ii];
    }
}

void WaveO1XtDG :: set_vector(const int n,
                            const Vector& slabB,
                            Vector& B) const
{
    Array<int> slab_vdofs, vdofs;
    m_tFespace0->GetElementVDofs (0, slab_vdofs);
    m_tFespace->GetElementVDofs (n, vdofs);

    int ndofs = slabB.Size();
    int tNdofs = vdofs.Size();
    int xNdofs = ndofs/tNdofs;

    const bool use_dev
            = vdofs.UseDevice() || slabB.UseDevice();
    auto d_slabB = slabB.Read(use_dev);
    auto d_B = B.ReadWrite(use_dev);

    for(int i=0; i<tNdofs; i++)
    {
        const int j = vdofs[i];
        const int ii = slab_vdofs[i]*xNdofs;
        const int jj = j*xNdofs;
        for (int k=0; k<xNdofs; k++)
            d_B[k + jj] = d_slabB[k + ii];
    }
}

// Assembles matrices for L2-projection
void WaveO1XtDG :: assemble_projector()
{
    /// assemble matrices in time

    // tMass form
    double ht = m_endTime/m_tFespace->GetNE();
    BilinearForm *tMass0_form
            = new BilinearForm(m_tFespace0);
    tMass0_form->AddDomainIntegrator(new MassIntegrator);
    tMass0_form->Assemble();
    tMass0_form->Finalize();
    m_slab_tMass = tMass0_form->LoseMat();
    delete tMass0_form;
    (*m_slab_tMass)*=ht;

    /// assemble matrices in space

    // xMass0 form, xW1_space
    BilinearForm *xMass0_form
            = new BilinearForm(m_xFespaces[0]);
    xMass0_form->AddDomainIntegrator(new MassIntegrator);
    xMass0_form->Assemble();
    xMass0_form->Finalize();
    m_xMass0 = xMass0_form->LoseMat();
    delete xMass0_form;

    // xMass1 form, xW2_space
    BilinearForm *xMass1_form
            = new BilinearForm(m_xFespaces[1]);
    xMass1_form->AddDomainIntegrator
            (new VectorMassIntegrator);
    xMass1_form->Assemble();
    xMass1_form->Finalize();
    m_xMass1 = xMass1_form->LoseMat();
    delete xMass1_form;

    /// compute blocks
    m_block00 = OuterProduct(*m_slab_tMass, *m_xMass0);
    m_block11 = OuterProduct(*m_slab_tMass, *m_xMass1);
}

// Builds projection matrix
// stores only the upper triangle
void WaveO1XtDG :: build_projector_matrix()
{
    BlockMatrix *waveProjBlockMat
            = new BlockMatrix(m_slab_block_offsets);
    waveProjBlockMat->SetBlock(0,0, m_block00);
    waveProjBlockMat->SetBlock(1,1, m_block11);

    SparseMatrix *waveProjMatBuf
            = waveProjBlockMat->CreateMonolithic();
    m_waveProjMat = &get_upper_triangle(*waveProjMatBuf);
    delete waveProjBlockMat;
    delete waveProjMatBuf;

    if (!m_waveProjMat->ColumnsAreSorted()) {
        m_waveProjMat->SortColumnIndices();
    }
}

// Assembles projection rhs for slab n
void WaveO1XtDG
:: assemble_projection_rhs(int n, BlockVector* slabB) const
{
    int xdimW1 = m_xFespaces[0]->GetTrueVSize();
    int xdimW2 = m_xFespaces[1]->GetTrueVSize();
    Vector bx1(xdimW1), bx2(xdimW2);

    ElementTransformation *trans = nullptr;
    const FiniteElement *fe = nullptr;
    Vector elvec1, elvec2, shape;
    {
        trans = m_tFespace->GetElementTransformation(n);

        // compute elvec
        fe = m_tFespace->GetFE(n);
        int tNdofs = fe->GetDof();
        shape.SetSize(tNdofs);
        elvec1.SetSize(tNdofs*xdimW1);
        elvec2.SetSize(tNdofs*xdimW2);

        int order = 2*fe->GetOrder()+2;
        const IntegrationRule *ir
                = &IntRules.Get(fe->GetGeomType(), order);

        elvec1 = 0.0;
        elvec2 = 0.0;
        for (int i = 0; i < ir->GetNPoints(); i++)
        {
            const IntegrationPoint &ip = ir->IntPoint(i);
            trans->SetIntPoint(&ip);
            fe->CalcShape(ip, shape);

            Vector t;
            trans->Transform(ip, t);
            assemble_xProjection_rhs(bx1, bx2, t(0));

            // elvec += w*shape*bx
            double w = ip.weight*trans->Weight();
            for (int j=0; j<tNdofs; j++){
                double coeff = w*shape(j);
                for (int k=0; k<xdimW1; k++) {
                    elvec1(k + j*xdimW1) += coeff*bx1(k);
                }
                for (int k=0; k<xdimW2; k++) {
                    elvec2(k + j*xdimW2) += coeff*bx2(k);
                }
            }
        }
    }

    Array<int> slab_vdofs;
    m_tFespace0->GetElementVDofs (0, slab_vdofs);
    add_vector(slab_vdofs, elvec1, slabB->GetBlock(0));
    add_vector(slab_vdofs, elvec2, slabB->GetBlock(1));
}

void WaveO1XtDG
:: assemble_xProjection_rhs(Vector& b1, Vector& b2,
                            double t) const
{
    LinearForm pSol_form(m_xFespaces[0]);
    WaveO1ExactPressureCoeff pCoeff(m_testCase);
    pCoeff.SetTime(t);
    pSol_form.AddDomainIntegrator
            (new DomainLFIntegrator(pCoeff));
    pSol_form.Assemble();
    b1 = pSol_form.GetData();

    LinearForm vSol_form(m_xFespaces[1]);
    WaveO1ExactVelocityCoeff vCoeff(m_testCase);
    vCoeff.SetTime(t);
    vSol_form.AddDomainIntegrator
            (new VectorDomainLFIntegrator(vCoeff));
    vSol_form.Assemble();
    b2 = vSol_form.GetData();
}

// End of file
