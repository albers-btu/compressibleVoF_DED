/*---------------------------------------------------------------------------*\
  =========                 |
  \\      /  F ield         | OpenFOAM: The Open Source CFD Toolbox
   \\    /   O peration     | Website:  https://openfoam.org
    \\  /    A nd           | Copyright (C) 2022-2025 OpenFOAM Foundation
     \\/     M anipulation  |
-------------------------------------------------------------------------------
License
    This file is part of OpenFOAM.

Class
    Foam::solvers::compressibleVoF_DED

Description
    Compressible two-phase VOF solver extended for DED / laser melt-pool
    physics (Phase A–C):

    Base: OpenFOAM compressibleVoF + PIMPLE, dual compressible thermo.

    Phase A: stable laser, radiation linearisation, continuous Psat,
             under-relaxed evaporation / Marangoni / recoil.

    Phase B: shared CSF interface measure, face Marangoni + recoil,
             outer-corrector re-coupling, mushy/sponge fields.

    Phase C: Picard latent-heat treatment, mass-conserving evaporation
             (α + continuity), Beer–Lambert bulk laser absorption option.

    Phase D: continuum powder / material addition on free surface under the
             laser (mass + energy + α / continuity / pressure coupling).

SourceFiles
    compressibleVoF_DED.C
    DEDPhysics.C

\*---------------------------------------------------------------------------*/

#ifndef compressibleVoF_DED_H
#define compressibleVoF_DED_H

#include "twoPhaseVoFSolver.H"
#include "compressibleTwoPhaseVoFMixture.H"
#include "compressibleInterPhaseTransportModel.H"
#include "compressibleInterPhaseThermophysicalTransportModel.H"

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

namespace Foam
{
namespace solvers
{

/*---------------------------------------------------------------------------*\
                          Class compressibleVoF_DED Declaration
\*---------------------------------------------------------------------------*/

class compressibleVoF_DED
:
    public twoPhaseVoFSolver
{

protected:

    // Phase properties

        compressibleTwoPhaseVoFMixture& mixture_;


    // Thermophysical properties

        volScalarField& p;
        volScalarField::Internal vDot;
        scalar vDotResidualAlpha;


    // Pressure reference

        Foam::pressureReference pressureReference_;


    // Kinematic properties

        surfaceScalarField alphaRhoPhi1;
        surfaceScalarField alphaRhoPhi2;
        volScalarField K;


    // DED / melt-pool fields

        volScalarField pRecoil_;
        volScalarField fLiquid_;
        //- Liquid fraction at previous time level (for latent-heat ddt)
        volScalarField fLiquidPrev_;
        volVectorField marangoniForce_;
        surfaceScalarField marangoniForcef_;
        volScalarField gasSpongeResistance_;
        volScalarField mushyResistance_;
        volScalarField evaporationSink_;
        //- Volumetric mass evaporation rate metal→gas [kg/m3/s]
        volScalarField mDotEvap_;
        //- Volumetric powder mass addition rate (gas/metal → metal) [kg/m3/s]
        volScalarField mDotPowder_;
        //- Thermal power density from powder enthalpy [W/m3] (can be negative)
        volScalarField powderHeat_;
        volScalarField interfaceDelta_;
        //- Laser volumetric power [W/m3] (diagnostic)
        volScalarField Qlaser_;


    // Cached DED model coefficients

        scalar laserPower_;
        scalar absorptivity_;
        scalar laserSpeed_;
        scalar xStart_;
        scalar yStart_;
        scalar laserRadius_;
        scalar laserRampTime_;
        scalar forceRelaxation_;
        scalar evaporationRelaxation_;
        scalar laserTopSurfaceBias_;
        scalar laserMaxDeltaT_;
        scalar laserMinIntegralFraction_;
        scalar PsatMax_;
        scalar TmaxEval_;
        scalar boilSoftWidth_;
        scalar alphaEvap_;
        Switch limitEvapByLaser_;
        Switch massConservingEvaporation_;
        label nLatentCorrectors_;
        scalar liquidFractionRelaxation_;
        //- Beer–Lambert absorption length [m]; <=0 disables bulk absorption
        scalar laserAbsorptionLength_;
        //- Fraction of deposited power routed to bulk (rest surface); [0,1]
        scalar laserBulkFraction_;

        //- Laser deposition model: "geometric" | "rayTrace"
        word laserModel_;
        //- Rays per side for rayTrace (total rays ≈ nRaySide^2)
        label nRaySide_;
        //- Max surface reflections after first hit (0 = single hit only)
        label maxLaserReflections_;
        //- Ray march step as fraction of typical cell length
        scalar laserRayStepFraction_;
        //- Aperture / start plane z (absolute); < -GREAT/2 → domain max z
        scalar laserRayOriginZ_;
        //- Minimum |μ|=|d·n| for absorption (grazing cutoff)
        scalar laserMinCosIncidence_;

        //- Continuum powder feed (material addition)
        Switch powderEnabled_;
        //- Nominal powder mass feed rate [kg/s]
        scalar powderFeedRate_;
        //- Catch efficiency η ∈ [0,1] → deposited ṁ = η * feedRate
        scalar powderCatchEfficiency_;
        //- Enthalpy temperature of added metal [K]
        scalar powderTemperature_;
        //- Gaussian radius for powder distribution [m] (default: laserRadius)
        scalar powderRadius_;
        //- Minimum liquid fraction on free surface to catch powder
        scalar powderMinLiquidFraction_;
        //- If true, solid powder pays latent heat Lf on capture
        Switch powderSolid_;

        scalar solidusTemperature_;
        scalar liquidusTemperature_;
        scalar Tboil_;
        scalar Lv_;
        scalar Mv_;
        scalar Rs_;
        scalar emissivity_;
        scalar Tamb_;

        dimensionedScalar Lf_;
        dimensionedScalar dSigma_dT_;
        dimensionedScalar mushyConstant_;
        scalar mushyEpsilon_;

        scalar spongeZStart_;
        scalar spongeZEnd_;
        scalar spongeTimeScale_;


    // Transport

        compressibleInterPhaseTransportModel momentumTransport;
        compressibleInterPhaseThermophysicalTransportModel
            thermophysicalTransport;


    // Cached temporary fields

        tmp<volScalarField::Internal> contErr1;
        tmp<volScalarField::Internal> contErr2;


    // Protected Member Functions

        virtual bool read();
        void readDEDProperties();
        void updateInterfaceDelta();
        void updateMushyAndSpongeResistance();

        //- Equilibrium liquid fraction from temperature
        scalar liquidFraction(const scalar T) const;

        //- dfL/dT for apparent heat capacity [1/K]
        scalar dLiquidFraction_dT(const scalar T) const;

        //- Update fLiquid_ from T with under-relaxation
        void updateLiquidFraction(const volScalarField& T);

        //- Build laser power density into Qlaser_ (dispatches model)
        void updateLaserSource(const volScalarField& CvEff);

        //- Geometric free-surface + projected-area laser (legacy / fallback)
        void updateLaserSourceGeometric(const volScalarField& CvEff);

        //- Ray-traced laser with optional multi-reflection (CPU + OpenMP)
        void updateLaserSourceRayTrace(const volScalarField& CvEff);

        //- Cap Qlaser_ by thermal headroom and redistribute leftover power
        void redistributeLaserPower
        (
            const volScalarField& CvEff,
            const scalar depositedPower,
            const scalar dt,
            label& nLaserLimited
        );

        //- Evaporation energy sink + optional mass transfer rate
        void updateEvaporation
        (
            const volScalarField& T,
            const volScalarField& CvEff,
            const volScalarField& rho1
        );

        //- Continuum powder mass + enthalpy source on free surface under laser
        void updatePowderSource(const volScalarField& T);

        //- Marangoni volumetric + CSF face force; recoil pressure
        void updateSurfaceForces(const volScalarField& T);

        virtual const Foam::pressureReference& pressureReference() const
        {
            return pressureReference_;
        }

        virtual bool incompressible() const
        {
            return mixture.incompressible();
        }

        virtual bool divergent() const
        {
            return
               !incompressible()
             || fvModels().addsSupToField(mixture_.rho1().name())
             || fvModels().addsSupToField(mixture_.rho2().name())
             || massConservingEvaporation_
             || (powderEnabled_ && powderFeedRate_ > SMALL);
        }

        virtual tmp<volScalarField> psiByRho() const
        {
            return mixture.psiByRho();
        }

        virtual void alphaSuSp
        (
            tmp<volScalarField::Internal>& Su,
            tmp<volScalarField::Internal>& Sp
        );

        virtual tmp<fvVectorMatrix> divDevTau(volVectorField& U);

public:

        const compressibleTwoPhaseVoFMixture& mixture;

        TypeName("compressibleVoF_DED");

        compressibleVoF_DED(fvMesh& mesh);
        compressibleVoF_DED(const compressibleVoF_DED&) = delete;
        virtual ~compressibleVoF_DED();

        virtual void prePredictor();
        virtual void momentumTransportPredictor();
        virtual void thermophysicalTransportPredictor();
        virtual void momentumPredictor();
        virtual void thermophysicalPredictor();
        virtual void pressureCorrector();
        virtual void momentumTransportCorrector();
        virtual void thermophysicalTransportCorrector();
        virtual void postSolve();

        void operator=(const compressibleVoF_DED&) = delete;
};


// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

} // End namespace solvers
} // End namespace Foam

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

#endif

// ************************************************************************* //
