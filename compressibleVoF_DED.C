/*---------------------------------------------------------------------------*\
  =========                 |
  \\      /  F ield         | OpenFOAM: The Open Source CFD Toolbox
   \\    /   O peration     | Website:  https://openfoam.org
    \\  /    A nd           | Copyright (C) 2022-2025 OpenFOAM Foundation
     \\/     M anipulation  |
-------------------------------------------------------------------------------
License
    This file is part of OpenFOAM.

\*---------------------------------------------------------------------------*/

#include "compressibleVoF_DED.H"
#include "localEulerDdtScheme.H"
#include "fvcDdt.H"
#include "fvcDiv.H"
#include "addToRunTimeSelectionTable.H"

// * * * * * * * * * * * * * * Static Data Members * * * * * * * * * * * * * //

namespace Foam
{
namespace solvers
{
    defineTypeNameAndDebug(compressibleVoF_DED, 0);
    addToRunTimeSelectionTable(solver, compressibleVoF_DED, fvMesh);
}
}


// * * * * * * * * * * * * * Protected Member Functions  * * * * * * * * * * //

bool Foam::solvers::compressibleVoF_DED::read()
{
    twoPhaseVoFSolver::read();

    const dictionary& alphaControls = mesh.solution().solverDict(alpha1.name());

    vDotResidualAlpha =
        alphaControls.lookupOrDefault("vDotResidualAlpha", 1e-4);

    readDEDProperties();

    return true;
}


// * * * * * * * * * * * * * * * * Constructors  * * * * * * * * * * * * * * //

Foam::solvers::compressibleVoF_DED::compressibleVoF_DED(fvMesh& mesh)
:
    twoPhaseVoFSolver
    (
        mesh,
        autoPtr<twoPhaseVoFMixture>(new compressibleTwoPhaseVoFMixture(mesh))
    ),

    mixture_
    (
        refCast<compressibleTwoPhaseVoFMixture>(twoPhaseVoFSolver::mixture)
    ),

    p(mixture_.p()),

    vDot
    (
        IOobject
        (
            "vDot",
            runTime.name(),
            mesh,
            IOobject::READ_IF_PRESENT,
            IOobject::AUTO_WRITE
        ),
        alpha1()*fvc::div(phi)()()
    ),

    pressureReference_
    (
        p,
        p_rgh,
        pimple.dict(),
        false
    ),

    alphaRhoPhi1
    (
        IOobject::groupName("alphaRhoPhi", alpha1.group()),
        fvc::interpolate(mixture_.thermo1().rho())*alphaPhi1
    ),

    alphaRhoPhi2
    (
        IOobject::groupName("alphaRhoPhi", alpha2.group()),
        fvc::interpolate(mixture_.thermo2().rho())*alphaPhi2
    ),

    K("K", 0.5*magSqr(U)),

    pRecoil_
    (
        IOobject
        (
            "pRecoil",
            runTime.name(),
            mesh,
            IOobject::READ_IF_PRESENT,
            IOobject::AUTO_WRITE
        ),
        mesh,
        dimensionedScalar(dimPressure, 0)
    ),

    fLiquid_
    (
        IOobject
        (
            "fLiquid",
            runTime.name(),
            mesh,
            IOobject::READ_IF_PRESENT,
            IOobject::AUTO_WRITE
        ),
        mesh,
        dimensionedScalar(dimless, 0)
    ),

    fLiquidPrev_
    (
        IOobject
        (
            "fLiquidPrev",
            runTime.name(),
            mesh,
            IOobject::READ_IF_PRESENT,
            IOobject::NO_WRITE
        ),
        fLiquid_
    ),

    marangoniForce_
    (
        IOobject
        (
            "marangoniForce",
            runTime.name(),
            mesh,
            IOobject::READ_IF_PRESENT,
            IOobject::AUTO_WRITE
        ),
        mesh,
        dimensionedVector
        (
            "zero",
            dimensionSet(1, -2, -2, 0, 0, 0, 0),
            vector::zero
        )
    ),

    marangoniForcef_
    (
        IOobject
        (
            "marangoniForcef",
            runTime.name(),
            mesh,
            IOobject::NO_READ,
            IOobject::NO_WRITE
        ),
        mesh,
        dimensionedScalar
        (
            "zero",
            dimensionSet(1, -2, -2, 0, 0, 0, 0),
            0
        )
    ),

    gasSpongeResistance_
    (
        IOobject
        (
            "gasSpongeResistance",
            runTime.name(),
            mesh,
            IOobject::NO_READ,
            IOobject::AUTO_WRITE
        ),
        mesh,
        dimensionedScalar
        (
            "zero",
            dimensionSet(1, -3, -1, 0, 0, 0, 0),
            0.0
        )
    ),

    mushyResistance_
    (
        IOobject
        (
            "mushyResistance",
            runTime.name(),
            mesh,
            IOobject::NO_READ,
            IOobject::AUTO_WRITE
        ),
        mesh,
        dimensionedScalar
        (
            "zero",
            dimensionSet(1, -3, -1, 0, 0, 0, 0),
            0.0
        )
    ),

    evaporationSink_
    (
        IOobject
        (
            "evaporationSink",
            runTime.name(),
            mesh,
            IOobject::READ_IF_PRESENT,
            IOobject::AUTO_WRITE
        ),
        mesh,
        dimensionedScalar
        (
            "zero",
            dimensionSet(1, -1, -3, 0, 0, 0, 0),
            0.0
        )
    ),

    mDotEvap_
    (
        IOobject
        (
            "mDotEvap",
            runTime.name(),
            mesh,
            IOobject::READ_IF_PRESENT,
            IOobject::AUTO_WRITE
        ),
        mesh,
        dimensionedScalar
        (
            "zero",
            dimensionSet(1, -3, -1, 0, 0, 0, 0),
            0.0
        )
    ),

    interfaceDelta_
    (
        IOobject
        (
            "interfaceDelta",
            runTime.name(),
            mesh,
            IOobject::NO_READ,
            IOobject::AUTO_WRITE
        ),
        mesh,
        dimensionedScalar("zero", dimless/dimLength, 0)
    ),

    Qlaser_
    (
        IOobject
        (
            "Qlaser",
            runTime.name(),
            mesh,
            IOobject::NO_READ,
            IOobject::AUTO_WRITE
        ),
        mesh,
        dimensionedScalar
        (
            "zero",
            dimensionSet(1, -1, -3, 0, 0, 0, 0),
            0
        )
    ),

    laserPower_(0),
    absorptivity_(0),
    laserSpeed_(0),
    xStart_(0),
    yStart_(0),
    laserRadius_(0),
    laserRampTime_(1e-3),
    forceRelaxation_(0.3),
    evaporationRelaxation_(0.3),
    laserTopSurfaceBias_(1.0),
    laserMaxDeltaT_(500),
    laserMinIntegralFraction_(0.05),
    PsatMax_(2e6),
    TmaxEval_(6000),
    boilSoftWidth_(50),
    alphaEvap_(1),
    limitEvapByLaser_(false),
    massConservingEvaporation_(true),
    nLatentCorrectors_(2),
    liquidFractionRelaxation_(0.7),
    laserAbsorptionLength_(0),
    laserBulkFraction_(0.25),
    solidusTemperature_(1700),
    liquidusTemperature_(1800),
    Tboil_(3000),
    Lv_(0),
    Mv_(0.056),
    Rs_(8.314462618/0.056),
    emissivity_(0.5),
    Tamb_(300),

    Lf_("Lf", dimEnergy/dimMass, 0),
    dSigma_dT_
    (
        "dSigma_dT",
        dimensionSet(1, 0, -2, -1, 0, 0, 0),
        0
    ),
    mushyConstant_
    (
        "mushyConstant",
        dimensionSet(1, -3, -1, 0, 0, 0, 0),
        1e8
    ),
    mushyEpsilon_(1e-6),

    spongeZStart_(0),
    spongeZEnd_(0),
    spongeTimeScale_(1e-5),

    momentumTransport
    (
        rho,
        U,
        phi,
        rhoPhi,
        alphaPhi1,
        alphaPhi2,
        alphaRhoPhi1,
        alphaRhoPhi2,
        mixture_
    ),

    thermophysicalTransport(momentumTransport),

    mixture(mixture_)
{
    read();
    updateInterfaceDelta();
    updateMushyAndSpongeResistance();

    if (correctPhi || mesh.topoChanging())
    {
        rAU = new volScalarField
        (
            IOobject
            (
                "rAU",
                runTime.name(),
                mesh,
                IOobject::READ_IF_PRESENT,
                IOobject::AUTO_WRITE
            ),
            mesh,
            dimensionedScalar(dimTime/dimDensity, 1)
        );
    }
}


Foam::solvers::compressibleVoF_DED::~compressibleVoF_DED()
{}


// * * * * * * * * * * * * * * Member Functions  * * * * * * * * * * * * * * //

void Foam::solvers::compressibleVoF_DED::prePredictor()
{
    readDEDProperties();

    twoPhaseVoFSolver::prePredictor();

    const volScalarField& rho1 = mixture_.thermo1().rho();
    const volScalarField& rho2 = mixture_.thermo2().rho();

    alphaRhoPhi1 = fvc::interpolate(rho1)*alphaPhi1;
    alphaRhoPhi2 = fvc::interpolate(rho2)*alphaPhi2;

    rhoPhi = alphaRhoPhi1 + alphaRhoPhi2;

    // Continuity residuals include Phase-C mass transfer metal→gas:
    //   metal: S1 = -mDotEvap
    //   gas:   S2 = +mDotEvap
    // contErr = ddt + div - S
    contErr1 =
    (
        fvc::ddt(alpha1, rho1)()() + fvc::div(alphaRhoPhi1)()()
      - (fvModels().source(alpha1, rho1)&rho1)()
      + mDotEvap_()
    );

    contErr2 =
    (
        fvc::ddt(alpha2, rho2)()() + fvc::div(alphaRhoPhi2)()()
      - (fvModels().source(alpha2, rho2)&rho2)()
      - mDotEvap_()
    );

    updateInterfaceDelta();
    updateMushyAndSpongeResistance();
}


void Foam::solvers::compressibleVoF_DED::momentumTransportPredictor()
{
    momentumTransport.predict();
}


void Foam::solvers::compressibleVoF_DED::thermophysicalTransportPredictor()
{
    thermophysicalTransport.predict();
}


void Foam::solvers::compressibleVoF_DED::momentumTransportCorrector()
{
    momentumTransport.correct();
}


void Foam::solvers::compressibleVoF_DED::thermophysicalTransportCorrector()
{
    thermophysicalTransport.correct();
}


void Foam::solvers::compressibleVoF_DED::postSolve()
{
    // Freeze time-level liquid fraction for next-step latent-heat ddt
    fLiquidPrev_ == fLiquid_;

    VoFSolver::postSolve();
}


// ************************************************************************* //
