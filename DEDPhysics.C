/*---------------------------------------------------------------------------*\
  =========                 |
  \\      /  F ield         | OpenFOAM: The Open Source CFD Toolbox
   \\    /   O peration     | Website:  https://openfoam.org
    \\  /    A nd           | Copyright (C) 2022-2025 OpenFOAM Foundation
     \\/     M anipulation  |
-------------------------------------------------------------------------------
License
    This file is part of OpenFOAM.

    Phase-C DED physics helpers: interface, mushy/sponge, liquid fraction,
    laser (surface + Beer–Lambert bulk), evaporation, surface forces.

\*---------------------------------------------------------------------------*/

#include "compressibleVoF_DED.H"
#include "fvcGrad.H"
#include "fvcFlux.H"
#include "constants.H"

// * * * * * * * * * * * * * Protected Member Functions  * * * * * * * * * * //

void Foam::solvers::compressibleVoF_DED::readDEDProperties()
{
    IOdictionary laserProperties
    (
        IOobject
        (
            "laserProperties",
            runTime.constant(),
            mesh,
            IOobject::MUST_READ_IF_MODIFIED,
            IOobject::NO_WRITE
        )
    );

    laserPower_ = laserProperties.lookup<scalar>("laserPower");
    absorptivity_ = laserProperties.lookup<scalar>("absorptivity");
    laserSpeed_ = laserProperties.lookup<scalar>("laserSpeed");
    xStart_ = laserProperties.lookup<scalar>("xStart");
    yStart_ = laserProperties.lookup<scalar>("yStart");
    laserRadius_ = laserProperties.lookup<scalar>("laserRadius");

    laserRampTime_ =
        laserProperties.lookupOrDefault<scalar>("laserRampTime", 1e-3);
    forceRelaxation_ =
        laserProperties.lookupOrDefault<scalar>("forceRelaxation", 0.3);
    evaporationRelaxation_ =
        laserProperties.lookupOrDefault<scalar>
        (
            "evaporationRelaxation",
            forceRelaxation_
        );
    laserTopSurfaceBias_ =
        laserProperties.lookupOrDefault<scalar>("laserTopSurfaceBias", 1.0);
    laserMaxDeltaT_ =
        laserProperties.lookupOrDefault<scalar>("laserMaxDeltaT", 500.0);
    laserMinIntegralFraction_ =
        laserProperties.lookupOrDefault<scalar>
        (
            "laserMinIntegralFraction",
            0.05
        );
    PsatMax_ = laserProperties.lookupOrDefault<scalar>("PsatMax", 2.0e6);
    TmaxEval_ = laserProperties.lookupOrDefault<scalar>("TmaxEval", 6000.0);
    boilSoftWidth_ =
        laserProperties.lookupOrDefault<scalar>("boilSoftWidth", 50.0);
    alphaEvap_ = laserProperties.lookupOrDefault<scalar>("alphaEvap", 1.0);
    limitEvapByLaser_ =
        laserProperties.lookupOrDefault<Switch>("limitEvapByLaser", false);

    // Phase C
    massConservingEvaporation_ =
        laserProperties.lookupOrDefault<Switch>
        (
            "massConservingEvaporation",
            true
        );
    nLatentCorrectors_ =
        laserProperties.lookupOrDefault<label>("nLatentCorrectors", 2);
    nLatentCorrectors_ = max(nLatentCorrectors_, 1);
    liquidFractionRelaxation_ =
        laserProperties.lookupOrDefault<scalar>
        (
            "liquidFractionRelaxation",
            0.7
        );
    laserAbsorptionLength_ =
        laserProperties.lookupOrDefault<scalar>
        (
            "laserAbsorptionLength",
            0.0
        );
    laserBulkFraction_ =
        laserProperties.lookupOrDefault<scalar>("laserBulkFraction", 0.25);
    laserBulkFraction_ = min(max(laserBulkFraction_, 0.0), 1.0);

    solidusTemperature_ =
        laserProperties.lookup<scalar>("solidusTemperature");
    liquidusTemperature_ =
        laserProperties.lookup<scalar>("liquidusTemperature");
    Tboil_ = laserProperties.lookup<scalar>("Tboil");
    Lv_ = laserProperties.lookup<scalar>("Lv");
    Mv_ = laserProperties.lookup<scalar>("Mv");
    Rs_ = 8.314462618/Mv_;
    emissivity_ = laserProperties.lookup<scalar>("emissivity");
    Tamb_ = laserProperties.lookup<scalar>("Tamb");

    Lf_.value() = laserProperties.lookup<scalar>("Lf");
    dSigma_dT_.value() = laserProperties.lookup<scalar>("dSigma_dT");

    mushyConstant_.value() =
        laserProperties.lookupOrDefault<scalar>("mushyConstant", 1e8);
    mushyEpsilon_ =
        laserProperties.lookupOrDefault<scalar>("mushyEpsilon", 1e-6);

    IOdictionary gasSpongeProperties
    (
        IOobject
        (
            "gasSpongeProperties",
            runTime.constant(),
            mesh,
            IOobject::MUST_READ_IF_MODIFIED,
            IOobject::NO_WRITE
        )
    );

    spongeZStart_ = gasSpongeProperties.lookup<scalar>("zStart");
    spongeZEnd_ = gasSpongeProperties.lookup<scalar>("zEnd");
    spongeTimeScale_ =
        gasSpongeProperties.lookupOrDefault<scalar>("timeScale", 1e-5);
}


Foam::scalar Foam::solvers::compressibleVoF_DED::liquidFraction
(
    const scalar T
) const
{
    if (T <= solidusTemperature_)
    {
        return 0.0;
    }
    if (T >= liquidusTemperature_)
    {
        return 1.0;
    }

    return
        (T - solidusTemperature_)
       /max(liquidusTemperature_ - solidusTemperature_, SMALL);
}


Foam::scalar Foam::solvers::compressibleVoF_DED::dLiquidFraction_dT
(
    const scalar T
) const
{
    if (T <= solidusTemperature_ || T >= liquidusTemperature_)
    {
        return 0.0;
    }

    return 1.0/max(liquidusTemperature_ - solidusTemperature_, SMALL);
}


void Foam::solvers::compressibleVoF_DED::updateInterfaceDelta()
{
    const volVectorField gradAlpha(fvc::grad(alpha1));
    const volScalarField magGradAlpha(mag(gradAlpha));

    forAll(interfaceDelta_, celli)
    {
        const scalar a = min(max(alpha1[celli], 0.0), 1.0);
        interfaceDelta_[celli] = 2.0*a*(1.0 - a)*magGradAlpha[celli];
    }

    interfaceDelta_.correctBoundaryConditions();
}


void Foam::solvers::compressibleVoF_DED::updateMushyAndSpongeResistance()
{
    forAll(mushyResistance_, celli)
    {
        const scalar a = min(max(alpha1[celli], 0.0), 1.0);
        const scalar fL = min(max(fLiquid_[celli], 0.0), 1.0);

        mushyResistance_[celli] =
            mushyConstant_.value()
           *a
           *sqr(1.0 - fL)
           /(pow3(fL) + mushyEpsilon_);
    }
    mushyResistance_.correctBoundaryConditions();

    const scalar dz = max(spongeZEnd_ - spongeZStart_, SMALL);

    forAll(gasSpongeResistance_, celli)
    {
        gasSpongeResistance_[celli] = 0.0;

        const scalar z = mesh.C()[celli].z();
        const scalar a = min(max(alpha1[celli], 0.0), 1.0);
        const scalar gasFraction = 1.0 - a;

        if (z > spongeZStart_ && gasFraction > 1e-3)
        {
            const scalar ramp =
                min(max((z - spongeZStart_)/dz, 0.0), 1.0);

            gasSpongeResistance_[celli] =
                rho[celli]
               *sqr(gasFraction)
               *sqr(ramp)
               /max(spongeTimeScale_, SMALL);
        }
    }
    gasSpongeResistance_.correctBoundaryConditions();
}


void Foam::solvers::compressibleVoF_DED::updateLiquidFraction
(
    const volScalarField& T
)
{
    const scalar omega =
        min(max(liquidFractionRelaxation_, 0.0), 1.0);

    forAll(fLiquid_, celli)
    {
        const scalar fEq = liquidFraction(T[celli]);
        fLiquid_[celli] =
            omega*fEq + (1.0 - omega)*fLiquid_[celli];
    }

    fLiquid_.correctBoundaryConditions();
}


void Foam::solvers::compressibleVoF_DED::updateLaserSource
(
    const volScalarField& CvEff
)
{
    const volVectorField gradAlpha(fvc::grad(alpha1));
    const volScalarField magGradAlpha(mag(gradAlpha));

    scalar powerFactor = 1.0;
    if (laserRampTime_ > SMALL)
    {
        powerFactor = min(1.0, runTime.value()/laserRampTime_);
    }

    const scalar depositedPower = absorptivity_*laserPower_*powerFactor;
    const scalar xLaser = xStart_ + laserSpeed_*runTime.value();
    const scalar sigma = max(laserRadius_/3.0, SMALL);
    const scalar dt = max(runTime.deltaTValue(), SMALL);

    // -----------------------------------------------------------------
    // Laser deposition is LOCALISED to the free surface (and a shallow
    // layer below it).  Never heat the full metal column under the beam
    // — that produces a vertical "needle" melt instead of a pool.
    // -----------------------------------------------------------------

    // Free-surface height: top of metal (α > 0.5)
    scalar zFree = -GREAT;
    forAll(alpha1, celli)
    {
        if (alpha1[celli] > 0.5)
        {
            zFree = max(zFree, mesh.C()[celli].z());
        }
    }
    reduce(zFree, maxOp<scalar>());
    if (zFree < -0.5*GREAT)
    {
        zFree = 0.005;
    }

    // Max penetration below free surface for any volumetric deposition
    // (a few absorption lengths, or a few cells if no bulk model)
    const scalar cellLenTyp = cbrt(gAverage(mesh.V()));
    const scalar zPen =
        (laserAbsorptionLength_ > SMALL)
      ? max(3.0*laserAbsorptionLength_, 2.0*cellLenTyp)
      : 3.0*cellLenTyp;

    // --- Horizontal Gaussian footprint (xy only) ---
    volScalarField gaussian
    (
        IOobject
        (
            "gaussian",
            runTime.name(),
            mesh,
            IOobject::NO_READ,
            IOobject::NO_WRITE
        ),
        mesh,
        dimensionedScalar("zero", dimless, 0)
    );

    forAll(gaussian, celli)
    {
        const point& C = mesh.C()[celli];
        const scalar dx = C.x() - xLaser;
        const scalar dy = C.y() - yStart_;
        gaussian[celli] = exp(-(dx*dx + dy*dy)/(2.0*sigma*sigma));
    }

    // --- Surface weight: CSF free-surface measure only ---
    // δ = 2 α (1-α) |∇α|  peaks at the metal/gas interface
    volScalarField surfaceWeight
    (
        IOobject
        (
            "laserSurfaceWeight",
            runTime.name(),
            mesh,
            IOobject::NO_READ,
            IOobject::NO_WRITE
        ),
        mesh,
        dimensionedScalar("zero", dimless/dimLength, 0)
    );

    forAll(surfaceWeight, celli)
    {
        const scalar a = min(max(alpha1[celli], 0.0), 1.0);
        const scalar magGa = magGradAlpha[celli];
        const scalar z = mesh.C()[celli].z();

        // Only deposit near free surface (not deep substrate / not high gas)
        if (z < zFree - zPen || z > zFree + 2.0*cellLenTyp)
        {
            continue;
        }

        scalar delta = 2.0*a*(1.0 - a)*magGa;
        delta = max(delta, interfaceDelta_[celli]);

        // Sharp MULES interfaces: also use |∇α| in the free-surface band
        if (a > 0.02 && a < 0.98)
        {
            delta = max(delta, magGa);
        }

        if (delta <= VSMALL)
        {
            continue;
        }

        scalar topWeight = 1.0;
        if (laserTopSurfaceBias_ > SMALL && magGa > VSMALL)
        {
            const vector nMetalToGas = -gradAlpha[celli]/(magGa + VSMALL);
            // Soft upward bias — floor so weight cannot vanish
            const scalar up = max(nMetalToGas.z(), 0.0);
            topWeight =
                max
                (
                    0.25,
                    (1.0 - laserTopSurfaceBias_) + laserTopSurfaceBias_*up
                );
        }

        surfaceWeight[celli] = delta*topWeight;
    }

    // --- Shallow Beer–Lambert bulk: depth measured DOWN from zFree ---
    // Only for cells in metal below the free surface within zPen.
    // (Old model used depth ~ (α-0.5)/|∇α| which is ~uniform in bulk → needle)
    volScalarField bulkWeight
    (
        IOobject
        (
            "laserBulkWeight",
            runTime.name(),
            mesh,
            IOobject::NO_READ,
            IOobject::NO_WRITE
        ),
        mesh,
        dimensionedScalar("zero", dimless/dimLength, 0)
    );

    if (laserAbsorptionLength_ > SMALL && laserBulkFraction_ > SMALL)
    {
        const scalar mu = 1.0/laserAbsorptionLength_;

        forAll(bulkWeight, celli)
        {
            const scalar a = min(max(alpha1[celli], 0.0), 1.0);
            const scalar z = mesh.C()[celli].z();

            if (a < 0.05 || gaussian[celli] < 1e-8)
            {
                continue;
            }

            // Depth into the metal from free surface (downward)
            const scalar depth = zFree - z;
            if (depth < 0 || depth > zPen)
            {
                continue;
            }

            bulkWeight[celli] = a*mu*exp(-mu*depth);
        }
    }

    // Blend: mostly free-surface CSF, optional shallow bulk
    const scalar fb =
        (laserAbsorptionLength_ > SMALL ? max(laserBulkFraction_, 0.0) : 0.0);
    const scalar fs = max(1.0 - fb, 0.0);

    volScalarField weight
    (
        IOobject
        (
            "laserWeight",
            runTime.name(),
            mesh,
            IOobject::NO_READ,
            IOobject::NO_WRITE
        ),
        gaussian*(fs*surfaceWeight + fb*bulkWeight)
    );

    // If free-surface measure is too weak, fall back to a SHALLOW surface
    // layer under the beam — never the full metal column.
    {
        scalar surfInt = 0.0;
        forAll(weight, celli)
        {
            surfInt += weight[celli]*mesh.V()[celli];
        }
        reduce(surfInt, sumOp<scalar>());

        if (surfInt < SMALL)
        {
            weight = dimensionedScalar("zero", dimless/dimLength, 0);
            forAll(weight, celli)
            {
                const scalar a = min(max(alpha1[celli], 0.0), 1.0);
                const scalar z = mesh.C()[celli].z();
                if
                (
                    gaussian[celli] > 1e-6
                 && a > 0.05
                 && z >= zFree - zPen
                 && z <= zFree + 2.0*cellLenTyp
                )
                {
                    // Compact surface kernel ~ 1/cell near free surface
                    const scalar depth = max(zFree - z, 0.0);
                    weight[celli] =
                        gaussian[celli]*a
                       *exp(-depth/max(zPen, SMALL))
                       /max(cellLenTyp, SMALL);
                }
            }
        }
    }

    // Zero any residual weight deep in the substrate / high in gas
    forAll(weight, celli)
    {
        const scalar z = mesh.C()[celli].z();
        if (z < zFree - zPen || z > zFree + 2.0*cellLenTyp)
        {
            weight[celli] = 0.0;
        }
    }

    scalar sourceIntegral = 0.0;
    forAll(weight, celli)
    {
        sourceIntegral += weight[celli]*mesh.V()[celli];
    }
    reduce(sourceIntegral, sumOp<scalar>());

    // sourceIntegral = ∫ weight dV ~ area (δ ~ 1/length)
    // If weight exists, scale so full depositedPower is applied.
    // Only skip deposition when the kernel is essentially empty.
    scalar scaleVal = 0.0;
    if (sourceIntegral > SMALL)
    {
        scaleVal = depositedPower/sourceIntegral;
    }

    const dimensionedScalar scale
    (
        "scale",
        dimensionSet(1, 0, -3, 0, 0, 0, 0),
        scaleVal
    );

    Qlaser_ = scale*weight;

    // For logging
    const scalar safeIntegral = max(sourceIntegral, SMALL);

    // Per-cell power density cap
    const volScalarField& Cv2 = mixture.thermo2().Cv();
    label nLaserLimited = 0;

    forAll(Qlaser_, celli)
    {
        const scalar a = min(max(alpha1[celli], 0.0), 1.0);
        const scalar CvMix = a*CvEff[celli] + (1.0 - a)*Cv2[celli];
        const scalar maxQ =
            rho[celli]*max(CvMix, SMALL)*laserMaxDeltaT_/dt;

        if (Qlaser_[celli] > maxQ)
        {
            Qlaser_[celli] = maxQ;
            ++nLaserLimited;
        }
    }
    reduce(nLaserLimited, sumOp<label>());

    Qlaser_.correctBoundaryConditions();

    if (debug || pimple.finalIter())
    {
        scalar appliedPower = 0.0;
        forAll(Qlaser_, celli)
        {
            appliedPower += Qlaser_[celli]*mesh.V()[celli];
        }
        reduce(appliedPower, sumOp<scalar>());

        Info<< "DED laser:"
            << " powerFactor=" << powerFactor
            << " deposited=" << depositedPower
            << " applied=" << appliedPower
            << " safeIntegral=" << safeIntegral
            << " bulkLen=" << laserAbsorptionLength_
            << " nLimited=" << nLaserLimited
            << " xLaser=" << xLaser
            << endl;
    }
}


void Foam::solvers::compressibleVoF_DED::updateEvaporation
(
    const volScalarField& T,
    const volScalarField& CvEff,
    const volScalarField& rho1
)
{
    const scalar dt = max(runTime.deltaTValue(), SMALL);
    const volVectorField gradAlpha(fvc::grad(alpha1));
    const volScalarField magGradAlpha(mag(gradAlpha));

    // Free-surface height (open melt pool top)
    scalar zFree = -GREAT;
    forAll(alpha1, celli)
    {
        if (alpha1[celli] > 0.5)
        {
            zFree = max(zFree, mesh.C()[celli].z());
        }
    }
    reduce(zFree, maxOp<scalar>());
    if (zFree < -0.5*GREAT)
    {
        zFree = 0.005;
    }

    volScalarField evaporationSinkNew
    (
        IOobject
        (
            "evaporationSinkNew",
            runTime.name(),
            mesh,
            IOobject::NO_READ,
            IOobject::NO_WRITE
        ),
        mesh,
        dimensionedScalar
        (
            "zero",
            dimensionSet(1, -1, -3, 0, 0, 0, 0),
            0
        )
    );

    volScalarField mDotNew
    (
        IOobject
        (
            "mDotNew",
            runTime.name(),
            mesh,
            IOobject::NO_READ,
            IOobject::NO_WRITE
        ),
        mesh,
        dimensionedScalar
        (
            "zero",
            dimensionSet(1, -3, -1, 0, 0, 0, 0),
            0
        )
    );

    scalar evaporationPowerRaw = 0.0;
    scalar maxMdot = 0.0;
    scalar maxPsat = 0.0;
    label limitedCells = 0;
    label nEvapCells = 0;
    label nHotCells = 0;

    forAll(evaporationSinkNew, celli)
    {
        const scalar a = min(max(alpha1[celli], 0.0), 1.0);
        const scalar magGa = magGradAlpha[celli];
        const scalar cellLen = cbrt(mesh.V()[celli]);
        const scalar z = mesh.C()[celli].z();

        const scalar deltaCSF = 2.0*a*(1.0 - a)*magGa;
        const scalar delta = max(max(interfaceDelta_[celli], deltaCSF), magGa);

        const bool freeSurfaceBand =
            (delta > VSMALL)
         && (a > 0.02) && (a < 0.98);

        const bool openPoolTop =
            (a > 0.25)
         && (z > zFree - 2.5*cellLen)
         && (z < zFree + 1.5*cellLen);

        if (!freeSurfaceBand && !openPoolTop)
        {
            continue;
        }

        // Only evaporate from liquid metal free surface
        const scalar fL =
            max(liquidFraction(T[celli]), fLiquid_[celli]);
        if (fL < 1e-3)
        {
            continue;
        }

        const scalar Tc = T[celli];
        if (Tc < solidusTemperature_)
        {
            continue;
        }

        const scalar Tsafe = min(max(Tc, SMALL), TmaxEval_);

        // Continuous Clausius–Clapeyron saturation pressure
        scalar Psat =
            101325.0
           *exp(Lv_/Rs_*(1.0/Tboil_ - 1.0/Tsafe));
        Psat = min(max(Psat, 0.0), PsatMax_);
        maxPsat = max(maxPsat, Psat);

        // Hertz–Knudsen mass flux [kg/m2/s]
        const scalar mdotA =
            alphaEvap_*Psat
           /sqrt(2.0*constant::mathematical::pi*Rs_*Tsafe);

        maxMdot = max(maxMdot, mdotA);

        // Surface measure for open-pool top cells without mixed α
        const scalar deltaEff =
            freeSurfaceBand ? delta : (1.0/max(cellLen, SMALL));

        // Volumetric rates (scale by liquid fraction)
        scalar mDotVol = fL*deltaEff*mdotA;   // kg/m3/s
        scalar qVol = mDotVol*Lv_;            // W/m3

        // Local energy limiter (available enthalpy above solidus)
        const scalar Eavailable =
            a
           *rho1[celli]
           *CvEff[celli]
           *max(Tc - solidusTemperature_, 0.0)
           *mesh.V()[celli];

        const scalar maxQVol = (0.5*Eavailable/dt)/mesh.V()[celli];

        if (qVol > maxQVol)
        {
            const scalar scale = maxQVol/max(qVol, SMALL);
            qVol *= scale;
            mDotVol *= scale;
            ++limitedCells;
        }

        // Mass limiter: cannot remove more metal than present
        const scalar maxMDotVol = 0.5*a*rho1[celli]/dt;
        if (mDotVol > maxMDotVol)
        {
            const scalar scale = maxMDotVol/max(mDotVol, SMALL);
            mDotVol *= scale;
            qVol *= scale;
        }

        evaporationSinkNew[celli] = qVol;
        mDotNew[celli] = mDotVol;
        evaporationPowerRaw += qVol*mesh.V()[celli];
        ++nEvapCells;

        if (Tc > Tboil_)
        {
            ++nHotCells;
        }
    }

    reduce(evaporationPowerRaw, sumOp<scalar>());
    reduce(maxMdot, maxOp<scalar>());
    reduce(maxPsat, maxOp<scalar>());
    reduce(limitedCells, sumOp<label>());
    reduce(nHotCells, sumOp<label>());
    reduce(nEvapCells, sumOp<label>());

    if (limitEvapByLaser_)
    {
        scalar depositedPower = absorptivity_*laserPower_;
        if (laserRampTime_ > SMALL)
        {
            depositedPower *=
                min(1.0, runTime.value()/laserRampTime_);
        }

        const scalar availablePower = max(depositedPower, 0.0);
        if
        (
            evaporationPowerRaw > availablePower
         && evaporationPowerRaw > SMALL
        )
        {
            const scalar s = availablePower/(evaporationPowerRaw + SMALL);
            evaporationSinkNew *= s;
            mDotNew *= s;
            evaporationPowerRaw = availablePower;
        }
    }

    const scalar omegaE = min(max(evaporationRelaxation_, 0.0), 1.0);
    evaporationSink_ =
        omegaE*evaporationSinkNew
      + (1.0 - omegaE)*evaporationSink_;
    evaporationSink_.correctBoundaryConditions();

    // Mass-conserving phase change metal → gas (default on)
    if (massConservingEvaporation_)
    {
        mDotEvap_ =
            omegaE*mDotNew
          + (1.0 - omegaE)*mDotEvap_;
    }
    else
    {
        // Energy-only evaporation: still keep mDot for diagnostics = 0 mass
        mDotEvap_ = dimensionedScalar
        (
            "zero",
            dimensionSet(1, -3, -1, 0, 0, 0, 0),
            0
        );
    }
    mDotEvap_.correctBoundaryConditions();

    if (debug || pimple.finalIter())
    {
        scalar evaporationPower = 0.0;
        scalar massRate = 0.0;
        forAll(evaporationSink_, celli)
        {
            evaporationPower += evaporationSink_[celli]*mesh.V()[celli];
            massRate += mDotEvap_[celli]*mesh.V()[celli];
        }
        reduce(evaporationPower, sumOp<scalar>());
        reduce(massRate, sumOp<scalar>());

        Info<< "DED evapor:"
            << " power=" << evaporationPower
            << " massRate=" << massRate
            << " max(mdotA)=" << maxMdot
            << " max(Psat)=" << maxPsat
            << " nCells=" << nEvapCells
            << " limitedCells=" << limitedCells
            << " hotCells=" << nHotCells
            << " massTransfer=" << massConservingEvaporation_
            << endl;
    }
}


void Foam::solvers::compressibleVoF_DED::updateSurfaceForces
(
    const volScalarField& T
)
{
    const volVectorField gradAlpha(fvc::grad(alpha1));
    const volScalarField magGradAlpha(mag(gradAlpha));
    const volVectorField gradT(fvc::grad(T));

    volVectorField marangoniNew
    (
        IOobject
        (
            "marangoniNew",
            runTime.name(),
            mesh,
            IOobject::NO_READ,
            IOobject::NO_WRITE
        ),
        mesh,
        dimensionedVector
        (
            "zero",
            dimensionSet(1, -2, -2, 0, 0, 0, 0),
            vector::zero
        )
    );

    scalar maxMarangoniForce = 0.0;
    label nMarangoniCells = 0;

    // Free-surface height estimate (metal substrate top ~ max z with α>0.5)
    scalar zFree = -GREAT;
    forAll(alpha1, celli)
    {
        if (alpha1[celli] > 0.5)
        {
            zFree = max(zFree, mesh.C()[celli].z());
        }
    }
    reduce(zFree, maxOp<scalar>());
    if (zFree < -0.5*GREAT)
    {
        zFree = 0.005; // case default free surface
    }

    forAll(marangoniNew, celli)
    {
        const scalar a = min(max(alpha1[celli], 0.0), 1.0);
        const scalar magGa = magGradAlpha[celli];
        const scalar fL =
            max(liquidFraction(T[celli]), fLiquid_[celli]);

        if (fL <= 1e-3 || a < 0.05)
        {
            continue;
        }

        const scalar deltaCSF = 2.0*a*(1.0 - a)*magGa;
        const scalar delta = max(max(interfaceDelta_[celli], deltaCSF), magGa);

        const scalar cellLen = cbrt(mesh.V()[celli]);
        const scalar z = mesh.C()[celli].z();

        // Path A: classical free-surface band (mixed α + |∇α|)
        const bool freeSurfaceBand =
            (delta > VSMALL)
         && (a > 0.02) && (a < 0.98)
         && (magGa*cellLen > 1e-4);

        // Path B: open melt-pool top layer (liquid metal near free surface height)
        // Activates when melt is open but α is nearly pure metal at the top.
        const bool openPoolTop =
            (fL > 0.25)
         && (a > 0.25)
         && (z > zFree - 2.5*cellLen)
         && (z < zFree + 1.5*cellLen);

        if (!freeSurfaceBand && !openPoolTop)
        {
            continue;
        }

        vector nHat = vector(0, 0, 1); // default free-surface normal (up)
        if (magGa > VSMALL)
        {
            // Prefer metal→gas normal; flip if it points into the metal
            nHat = -gradAlpha[celli]/(magGa + VSMALL);
            if (nHat.z() < 0)
            {
                nHat = -nHat;
            }
        }

        const vector gradTangent =
            gradT[celli] - nHat*(nHat & gradT[celli]);

        if (mag(gradTangent) < VSMALL)
        {
            continue;
        }

        // Surface measure: CSF δ, or 1/cellLen for open-pool top cells
        const scalar deltaEff =
            freeSurfaceBand ? delta : (1.0/max(cellLen, SMALL));

        // F_M = f_L (dσ/dT) δ ∇_∥T   [N/m³]
        marangoniNew[celli] =
            fL*dSigma_dT_.value()*deltaEff*gradTangent;

        maxMarangoniForce =
            max(maxMarangoniForce, mag(marangoniNew[celli]));
        ++nMarangoniCells;
    }
    reduce(maxMarangoniForce, maxOp<scalar>());
    reduce(nMarangoniCells, sumOp<label>());

    // Mild under-relaxation for stability; keep responsive on melt pool
    const scalar omegaF = min(max(forceRelaxation_, 0.0), 1.0);
    marangoniForce_ =
        omegaF*marangoniNew
      + (1.0 - omegaF)*marangoniForce_;
    marangoniForce_.correctBoundaryConditions();

    // Face projection kept for diagnostics / optional face path;
    // primary coupling is volumetric (tangential) via HbyA.
    marangoniForcef_ = fvc::flux(marangoniForce_)/mesh.magSf();

    if (pimple.finalIter())
    {
        Info<< "DED Marangoni:"
            << " enabled cells=" << nMarangoniCells
            << " max|F_M|=" << maxMarangoniForce
            << " dSigma_dT=" << dSigma_dT_.value()
            << endl;
    }

    // ---------------------------------------------------------------
    // Recoil pressure (Anisimov): p_r = 0.54 Psat(T)
    // Continuous Psat; applied on free surface / open-pool top.
    // Coupled in pressure/momentum as p_r * snGrad(α) (CSF).
    // ---------------------------------------------------------------
    volScalarField pRecoilNew
    (
        IOobject
        (
            "pRecoilNew",
            runTime.name(),
            mesh,
            IOobject::NO_READ,
            IOobject::NO_WRITE
        ),
        mesh,
        dimensionedScalar("zero", dimPressure, 0)
    );

    label nRecoilCells = 0;
    scalar maxPsatRecoil = 0.0;

    forAll(pRecoilNew, celli)
    {
        const scalar a = min(max(alpha1[celli], 0.0), 1.0);
        const scalar magGa = magGradAlpha[celli];
        const scalar cellLen = cbrt(mesh.V()[celli]);
        const scalar z = mesh.C()[celli].z();
        const scalar fL =
            max(liquidFraction(T[celli]), fLiquid_[celli]);

        const scalar deltaCSF = 2.0*a*(1.0 - a)*magGa;
        const scalar delta = max(max(interfaceDelta_[celli], deltaCSF), magGa);

        const bool freeSurfaceBand =
            (delta > VSMALL)
         && (a > 0.02) && (a < 0.98);

        const bool openPoolTop =
            (fL > 0.25)
         && (a > 0.25)
         && (z > zFree - 2.5*cellLen)
         && (z < zFree + 1.5*cellLen);

        if ((!freeSurfaceBand && !openPoolTop) || fL < 1e-3)
        {
            continue;
        }

        const scalar Tc = T[celli];
        if (Tc < solidusTemperature_)
        {
            continue;
        }

        const scalar Tsafe = min(max(Tc, SMALL), TmaxEval_);

        scalar Psat =
            101325.0
           *exp(Lv_/Rs_*(1.0/Tboil_ - 1.0/Tsafe));
        Psat = min(max(Psat, 0.0), PsatMax_);
        maxPsatRecoil = max(maxPsatRecoil, Psat);

        // Mild soft factor near boiling (never fully kills continuous Psat)
        scalar soft = 1.0;
        if (boilSoftWidth_ > SMALL)
        {
            soft =
                0.5*(1.0 + tanh((Tc - Tboil_)/boilSoftWidth_));
            soft = min(max(soft, 0.05), 1.0);
        }

        // Anisimov recoil pressure
        pRecoilNew[celli] = 0.54*Psat*soft;
        ++nRecoilCells;
    }

    reduce(nRecoilCells, sumOp<label>());
    reduce(maxPsatRecoil, maxOp<scalar>());

    pRecoil_ =
        omegaF*pRecoilNew
      + (1.0 - omegaF)*pRecoil_;
    pRecoil_.correctBoundaryConditions();

    if (debug || pimple.finalIter())
    {
        Info<< "DED forces:"
            << " max|F_M|=" << maxMarangoniForce
            << " max(pRecoil)=" << gMax(pRecoil_)
            << " nRecoilCells=" << nRecoilCells
            << " max(Psat)=" << maxPsatRecoil
            << endl;
    }
}


// ************************************************************************* //
