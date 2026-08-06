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
#include "Pstream.H"
#include "meshSearch.H"

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

    laserModel_ =
        laserProperties.lookupOrDefault<word>("laserModel", "rayTrace");
    nRaySide_ =
        max(laserProperties.lookupOrDefault<label>("nRaySide", 32), 4);
    maxLaserReflections_ =
        max
        (
            laserProperties.lookupOrDefault<label>("maxLaserReflections", 3),
            0
        );
    laserRayStepFraction_ =
        laserProperties.lookupOrDefault<scalar>("laserRayStepFraction", 0.35);
    laserRayOriginZ_ =
        laserProperties.lookupOrDefault<scalar>
        (
            "laserRayOriginZ",
            -GREAT
        );
    laserMinCosIncidence_ =
        laserProperties.lookupOrDefault<scalar>
        (
            "laserMinCosIncidence",
            0.05
        );

    // Continuum powder / material addition
    powderEnabled_ =
        laserProperties.lookupOrDefault<Switch>("powderEnabled", false);
    powderFeedRate_ =
        laserProperties.lookupOrDefault<scalar>("powderFeedRate", 0.0);
    powderCatchEfficiency_ =
        min
        (
            max
            (
                laserProperties.lookupOrDefault<scalar>
                (
                    "powderCatchEfficiency",
                    0.7
                ),
                0.0
            ),
            1.0
        );
    powderTemperature_ =
        laserProperties.lookupOrDefault<scalar>("powderTemperature", 300.0);
    powderRadius_ =
        laserProperties.lookupOrDefault<scalar>("powderRadius", -1.0);
    if (powderRadius_ <= SMALL)
    {
        powderRadius_ = laserRadius_;
    }
    powderMinLiquidFraction_ =
        laserProperties.lookupOrDefault<scalar>
        (
            "powderMinLiquidFraction",
            0.1
        );
    powderSolid_ =
        laserProperties.lookupOrDefault<Switch>("powderSolid", true);

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
    if (laserModel_ == "rayTrace" || laserModel_ == "ray")
    {
        updateLaserSourceRayTrace(CvEff);
    }
    else
    {
        // "geometric" (default fallback name) or any other word
        updateLaserSourceGeometric(CvEff);
    }
}


void Foam::solvers::compressibleVoF_DED::redistributeLaserPower
(
    const volScalarField& CvEff,
    const scalar depositedPower,
    const scalar dt,
    label& nLaserLimited
)
{
    // Cap: Q ≤ ρ Cv ΔT_max / Δt
    // Leftover power after saturation is reassigned to unsaturated cells.
    const volScalarField& Cv2 = mixture.thermo2().Cv();

    scalarField maxQ(mesh.nCells(), 0.0);
    forAll(maxQ, celli)
    {
        const scalar a = min(max(alpha1[celli], 0.0), 1.0);
        const scalar CvMix = a*CvEff[celli] + (1.0 - a)*Cv2[celli];
        maxQ[celli] = rho[celli]*max(CvMix, SMALL)*laserMaxDeltaT_/dt;
    }

    // Treat current Qlaser_ as an unnormalized weight field (power density).
    // Convert to weights proportional to Q*V for redistribution of depositedPower.
    scalarField weight(mesh.nCells(), 0.0);
    forAll(weight, celli)
    {
        weight[celli] = max(Qlaser_[celli], 0.0);
    }

    Qlaser_ = dimensionedScalar
    (
        "zero",
        dimensionSet(1, -1, -3, 0, 0, 0, 0),
        0
    );

    boolList saturated(mesh.nCells(), false);
    const label nRedistrib = 8;

    for (label iter = 0; iter < nRedistrib; ++iter)
    {
        scalar powerOnSat = 0.0;
        scalar freeIntegral = 0.0;

        forAll(weight, celli)
        {
            if (weight[celli] <= VSMALL)
            {
                continue;
            }

            if (saturated[celli])
            {
                powerOnSat += maxQ[celli]*mesh.V()[celli];
            }
            else
            {
                freeIntegral += weight[celli]*mesh.V()[celli];
            }
        }
        reduce(powerOnSat, sumOp<scalar>());
        reduce(freeIntegral, sumOp<scalar>());

        const scalar remainingPower =
            max(depositedPower - powerOnSat, 0.0);

        if
        (
            remainingPower <= SMALL*max(depositedPower, 1.0)
         || freeIntegral <= SMALL
        )
        {
            if (freeIntegral <= SMALL)
            {
                forAll(weight, celli)
                {
                    if (!saturated[celli])
                    {
                        Qlaser_[celli] = 0;
                    }
                }
            }
            break;
        }

        const scalar scaleVal = remainingPower/freeIntegral;
        label nNewSat = 0;

        forAll(weight, celli)
        {
            if (weight[celli] <= VSMALL)
            {
                Qlaser_[celli] = 0;
                continue;
            }

            if (saturated[celli])
            {
                Qlaser_[celli] = maxQ[celli];
                continue;
            }

            const scalar qTry = scaleVal*weight[celli];
            if (qTry >= maxQ[celli]*(1.0 - 1e-12))
            {
                Qlaser_[celli] = maxQ[celli];
                saturated[celli] = true;
                ++nNewSat;
            }
            else
            {
                Qlaser_[celli] = qTry;
            }
        }
        reduce(nNewSat, sumOp<label>());

        if (nNewSat == 0)
        {
            break;
        }
    }

    nLaserLimited = 0;
    forAll(Qlaser_, celli)
    {
        if
        (
            weight[celli] > VSMALL
         && Qlaser_[celli] >= (1.0 - 1e-6)*maxQ[celli]
        )
        {
            ++nLaserLimited;
        }
    }
    reduce(nLaserLimited, sumOp<label>());

    Qlaser_.correctBoundaryConditions();
}


void Foam::solvers::compressibleVoF_DED::updateLaserSourceGeometric
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
    // Geometric free-surface laser (anti-donut / keyhole-ready):
    //   - Deposit on geometric free surface (interface + metal|gas faces)
    //     including crater floor — not only a global zFree CSF rim.
    //   - Projected-area weight: max(n_metal→gas · e_up, n_min)
    //     so crater bottom (n_z~1) gets heat; steep walls get less.
    //   - Optional shallow Beer–Lambert bulk from LOCAL free-surface height.
    //   - Per-cell ΔT cap with iterative redistribution of remaining power
    //     so applied ≈ deposited (no silent power loss).
    // -----------------------------------------------------------------

    // Global metal top (diagnostics + gas-side clamp only)
    scalar zMetalTop = -GREAT;
    forAll(alpha1, celli)
    {
        if (alpha1[celli] > 0.5)
        {
            zMetalTop = max(zMetalTop, mesh.C()[celli].z());
        }
    }
    reduce(zMetalTop, maxOp<scalar>());
    if (zMetalTop < -0.5*GREAT)
    {
        zMetalTop = 0.005;
    }

    const scalar cellLenTyp = cbrt(gAverage(mesh.V()));
    const scalar zPen =
        (laserAbsorptionLength_ > SMALL)
      ? max(3.0*laserAbsorptionLength_, 2.0*cellLenTyp)
      : 3.0*cellLenTyp;

    // Beam from +z downward; incidence uses upward unit e_up
    const vector eUp(0, 0, 1);
    // Minimum projected factor (walls keep a little heat; floor dominates)
    const scalar nMinProj = 0.05;
    // Optional soft blend with old top-bias (0 = pure projected area)
    const scalar topBias = min(max(laserTopSurfaceBias_, 0.0), 1.0);

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

    // --- Geometric free-surface marker ---
    // (A) mixed interface band, or (B) metal cell sharing a face with gas
    boolList isFreeSurface(mesh.nCells(), false);

    forAll(alpha1, celli)
    {
        const scalar a = min(max(alpha1[celli], 0.0), 1.0);
        const scalar magGa = magGradAlpha[celli];
        const scalar cellLen = cbrt(mesh.V()[celli]);

        if (a > 0.02 && a < 0.98 && magGa*cellLen > 1e-6)
        {
            isFreeSurface[celli] = true;
        }
    }

    {
        const labelUList& own = mesh.owner();
        const labelUList& nei = mesh.neighbour();

        forAll(own, facei)
        {
            const label c0 = own[facei];
            const label c1 = nei[facei];
            const scalar a0 = min(max(alpha1[c0], 0.0), 1.0);
            const scalar a1 = min(max(alpha1[c1], 0.0), 1.0);

            // Metal | gas face → both sides treated as geometric free surface
            // candidates; weight is applied only on metal-bearing cells below.
            if ((a0 > 0.5 && a1 < 0.5) || (a1 > 0.5 && a0 < 0.5))
            {
                isFreeSurface[c0] = true;
                isFreeSurface[c1] = true;
            }
        }

        // Coupled processor patches: use boundary α
        forAll(mesh.boundary(), patchi)
        {
            const fvPatch& patch = mesh.boundary()[patchi];
            if (!patch.coupled())
            {
                continue;
            }

            const fvPatchScalarField& alphaPf =
                alpha1.boundaryField()[patchi];
            const labelUList& fc = patch.faceCells();

            forAll(fc, facei)
            {
                const label c0 = fc[facei];
                const scalar a0 = min(max(alpha1[c0], 0.0), 1.0);
                const scalar aN = min(max(alphaPf[facei], 0.0), 1.0);

                if ((a0 > 0.5 && aN < 0.5) || (a0 < 0.5 && aN > 0.5))
                {
                    isFreeSurface[c0] = true;
                }
            }
        }
    }

    // Local free-surface height for bulk depth (only free-surface cells)
    // Under the beam: use nearest FS sample in xy (all ranks) for each bulk cell.
    DynamicList<label> fsCells(1024);
    DynamicList<vector> fsCentersLocal(1024);

    forAll(isFreeSurface, celli)
    {
        if (!isFreeSurface[celli] || gaussian[celli] < 1e-8)
        {
            continue;
        }

        const scalar a = min(max(alpha1[celli], 0.0), 1.0);
        // Keep metal-bearing free-surface cells (crater floor + walls + flat)
        if (a < 0.05)
        {
            continue;
        }

        fsCells.append(celli);
        fsCentersLocal.append(mesh.C()[celli]);
    }

    // Parallel: bulk depth needs FS samples from every rank under the beam
    List<vector> fsCenters;
    {
        List<List<vector>> gathered(Pstream::nProcs());
        gathered[Pstream::myProcNo()] = fsCentersLocal;
        Pstream::gatherList(gathered);
        Pstream::scatterList(gathered);

        label nAll = 0;
        forAll(gathered, proci)
        {
            nAll += gathered[proci].size();
        }
        fsCenters.setSize(nAll);
        label k = 0;
        forAll(gathered, proci)
        {
            forAll(gathered[proci], i)
            {
                fsCenters[k++] = gathered[proci][i];
            }
        }
    }

    // --- Surface weight: geometric FS × projected area ---
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

    forAll(fsCells, i)
    {
        const label celli = fsCells[i];
        const scalar a = min(max(alpha1[celli], 0.0), 1.0);
        const scalar magGa = magGradAlpha[celli];
        const scalar cellLen = cbrt(mesh.V()[celli]);
        const scalar z = mesh.C()[celli].z();

        // Drop pure gas and cells far above the metal domain
        if (a < 0.05 || z > zMetalTop + 2.0*cellLenTyp)
        {
            continue;
        }

        // Metal→gas normal; default upward on flat open metal
        vector nHat = eUp;
        if (magGa > VSMALL)
        {
            nHat = -gradAlpha[celli]/(magGa + VSMALL);
            // Ensure outward from metal (prefer +z when ambiguous)
            if ((nHat & eUp) < 0 && a > 0.5)
            {
                nHat = -nHat;
            }
        }

        // Projected area facing the vertical beam
        const scalar cosInc = max(nHat & eUp, nMinProj);

        // Soft optional bias (legacy knob); projected area is primary
        const scalar projWeight =
            (1.0 - topBias)*cosInc
          + topBias*max(cosInc, 0.25);

        // Interface area density [1/m]
        scalar delta = max(interfaceDelta_[celli], 2.0*a*(1.0 - a)*magGa);
        if (a > 0.02 && a < 0.98)
        {
            delta = max(delta, magGa);
        }
        // Pure-metal geometric FS (crater floor): compact 1/Δ kernel
        if (delta <= VSMALL || a >= 0.98)
        {
            delta = max(delta, 1.0/max(cellLen, SMALL));
        }

        surfaceWeight[celli] = delta*projWeight;
    }

    // --- Shallow Beer–Lambert bulk from LOCAL free-surface height ---
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

    if
    (
        laserAbsorptionLength_ > SMALL
     && laserBulkFraction_ > SMALL
     && fsCenters.size()
    )
    {
        const scalar mu = 1.0/laserAbsorptionLength_;
        // Search radius in xy for matching a free-surface sample
        const scalar rSearch = max(2.5*cellLenTyp, 0.5*laserRadius_);
        const scalar rSearchSqr = rSearch*rSearch;

        forAll(bulkWeight, celli)
        {
            const scalar a = min(max(alpha1[celli], 0.0), 1.0);
            if (a < 0.05 || gaussian[celli] < 1e-8)
            {
                continue;
            }

            // Skip pure free-surface cells (handled by surface kernel)
            if (isFreeSurface[celli] && a < 0.98)
            {
                continue;
            }

            const point& C = mesh.C()[celli];
            scalar bestD2 = GREAT;
            scalar zFs = C.z();

            forAll(fsCenters, i)
            {
                const vector d = fsCenters[i] - C;
                const scalar d2 = d.x()*d.x() + d.y()*d.y();
                if (d2 < bestD2)
                {
                    bestD2 = d2;
                    zFs = fsCenters[i].z();
                }
            }

            if (bestD2 > rSearchSqr)
            {
                continue;
            }

            // Depth into metal from the local free surface (downward only)
            const scalar depth = zFs - C.z();
            if (depth < 0 || depth > zPen)
            {
                continue;
            }

            bulkWeight[celli] = a*mu*exp(-mu*depth);
        }
    }

    const scalar fb =
        (laserAbsorptionLength_ > SMALL ? max(laserBulkFraction_, 0.0) : 0.0);
    const scalar fSurf = max(1.0 - fb, 0.0);

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
        gaussian*(fSurf*surfaceWeight + fb*bulkWeight)
    );

    // Fallback: compact geometric metal layer under the beam
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
                if (gaussian[celli] < 1e-6 || a < 0.05)
                {
                    continue;
                }

                if (isFreeSurface[celli] || a > 0.5)
                {
                    const scalar cellLen = cbrt(mesh.V()[celli]);
                    const scalar depth = max(zMetalTop - mesh.C()[celli].z(), 0.0);
                    if (depth <= zPen)
                    {
                        weight[celli] =
                            gaussian[celli]*a
                           *exp(-depth/max(zPen, SMALL))
                           /max(cellLen, SMALL);
                    }
                }
            }
        }
    }

    // Seed Qlaser_ with weight (power density shape); redistribute to
    // conserve depositedPower under per-cell thermal caps.
    Qlaser_ = dimensionedScalar
    (
        "zero",
        dimensionSet(1, -1, -3, 0, 0, 0, 0),
        0
    );
    forAll(Qlaser_, celli)
    {
        Qlaser_[celli] = weight[celli];
    }

    scalar sourceIntegral = 0.0;
    forAll(weight, celli)
    {
        sourceIntegral += weight[celli]*mesh.V()[celli];
    }
    reduce(sourceIntegral, sumOp<scalar>());
    const scalar safeIntegral = max(sourceIntegral, SMALL);

    label nLaserLimited = 0;
    redistributeLaserPower(CvEff, depositedPower, dt, nLaserLimited);

    if (debug || pimple.finalIter())
    {
        scalar appliedPower = 0.0;
        forAll(Qlaser_, celli)
        {
            appliedPower += Qlaser_[celli]*mesh.V()[celli];
        }
        reduce(appliedPower, sumOp<scalar>());

        Info<< "DED laser (geometric):"
            << " powerFactor=" << powerFactor
            << " deposited=" << depositedPower
            << " applied=" << appliedPower
            << " safeIntegral=" << safeIntegral
            << " bulkLen=" << laserAbsorptionLength_
            << " nLimited=" << nLaserLimited
            << " nFS=" << returnReduce(label(fsCells.size()), sumOp<label>())
            << " xLaser=" << xLaser
            << endl;
    }
}


void Foam::solvers::compressibleVoF_DED::updateLaserSourceRayTrace
(
    const volScalarField& CvEff
)
{
    // -----------------------------------------------------------------
    // Ray-traced laser — energy-conserving physics
    //
    // Beam: incident power P_inc = laserPower * ramp (NOT pre-scaled by A).
    // At each free-surface hit (opaque-metal partition):
    //   μ = max(−d·n, 0)   (n = metal→gas)
    //   A(θ) = A0 * μ_eff   with μ_eff = max(μ, μ_min)  → A(normal)=A0
    //   I_abs  = A * I      deposited in metal (Beer–Lambert skin)
    //   I_refl = (1−A) * I  continues for multi-reflection
    //
    // MPI: every rank marches geometrically, but ONLY the rank that owns
    // the global first hit (min path length) may deposit for that segment.
    // That prevents double-counting when partitions split free surface/bulk.
    //
    // Thermal ΔT caps (redistributeLaserPower) are numerical safety only.
    // -----------------------------------------------------------------

    scalar powerFactor = 1.0;
    if (laserRampTime_ > SMALL)
    {
        powerFactor = min(1.0, runTime.value()/laserRampTime_);
    }

    const scalar incidentPower = laserPower_*powerFactor;
    const scalar xLaser = xStart_ + laserSpeed_*runTime.value();
    const scalar sigma = max(laserRadius_/3.0, SMALL);
    const scalar dt = max(runTime.deltaTValue(), SMALL);
    const label myProci = Pstream::myProcNo();

    const volVectorField gradAlpha(fvc::grad(alpha1));
    const volScalarField magGradAlpha(mag(gradAlpha));

    // Global domain bounds
    const boundBox& bb = mesh.bounds();
    scalar zTop = bb.max().z();
    scalar zBot = bb.min().z();
    reduce(zTop, maxOp<scalar>());
    reduce(zBot, minOp<scalar>());

    scalar zOrigin = laserRayOriginZ_;
    if (zOrigin < -0.5*GREAT)
    {
        zOrigin = zTop - SMALL;
    }

    const scalar cellLenTyp = cbrt(gAverage(mesh.V()));
    const scalar ds =
        max(laserRayStepFraction_*cellLenTyp, 0.1*cellLenTyp);
    const scalar maxPath = 2.0*mag(zTop - zBot) + 4.0*laserRadius_;
    const label maxSteps = max(label(maxPath/ds) + 10, 50);

    const label nSide = nRaySide_;
    const label nRays = nSide*nSide;
    const scalar halfW = 3.0*sigma;

    List<scalar> rayW(nRays, 0.0);
    List<vector> rayOrig(nRays, vector::zero);
    const vector d0(0, 0, -1);

    scalar wSum = 0.0;
    label ri = 0;
    for (label j = 0; j < nSide; ++j)
    {
        const scalar fy = (nSide == 1) ? 0.5 : scalar(j)/(nSide - 1);
        const scalar y = yStart_ + (fy - 0.5)*2.0*halfW;

        for (label i = 0; i < nSide; ++i)
        {
            const scalar fx = (nSide == 1) ? 0.5 : scalar(i)/(nSide - 1);
            const scalar x = xLaser + (fx - 0.5)*2.0*halfW;
            const scalar dx = x - xLaser;
            const scalar dy = y - yStart_;
            const scalar G = exp(-(dx*dx + dy*dy)/(2.0*sigma*sigma));

            rayW[ri] = G;
            rayOrig[ri] = vector(x, y, zOrigin);
            wSum += G;
            ++ri;
        }
    }

    if (wSum <= SMALL || incidentPower <= SMALL)
    {
        Qlaser_ = dimensionedScalar
        (
            "zero",
            dimensionSet(1, -1, -3, 0, 0, 0, 0),
            0
        );
        return;
    }

    // Normalise weights → each ray carries a fraction of incident power
    forAll(rayW, i)
    {
        rayW[i] = incidentPower*(rayW[i]/wSum);
    }

    scalarField powerCell(mesh.nCells(), 0.0);
    const meshSearch& searchEng = meshSearch::New(mesh);

    label nHitsLocal = 0;
    label nReflectLocal = 0;
    label nSegments = 0;

    // One segment = one straight path (first arrival or after a reflection).
    // max segments per ray = 1 + maxLaserReflections_
    const label nSegMax = maxLaserReflections_ + 1;

    for (label rayI = 0; rayI < nRays; ++rayI)
    {
        scalar Iray = rayW[rayI];
        if (Iray <= SMALL)
        {
            continue;
        }

        point p0 = rayOrig[rayI];
        vector d = d0;

        for (label seg = 0; seg < nSegMax && Iray > SMALL; ++seg)
        {
            ++nSegments;

            // ---- Pass A: local first free-surface / gas→metal hit ----
            scalar localT = GREAT;
            label localCell = -1;
            point localHit = p0;
            vector localN(0, 0, 1);

            {
                point p = p0;
                scalar prevA = -1.0;
                scalar t = 0.0;

                for (label s = 0; s < maxSteps; ++s)
                {
                    p += d*ds;
                    t += ds;

                    if (p.z() < zBot - ds || p.z() > zTop + 2*ds)
                    {
                        break;
                    }

                    const label celli = searchEng.findCell(p);
                    if (celli < 0)
                    {
                        // Keep prevA (do not invent a new "first entry"
                        // when re-entering another partition's bulk).
                        continue;
                    }

                    const scalar a = min(max(alpha1[celli], 0.0), 1.0);
                    const scalar cellLen = cbrt(mesh.V()[celli]);
                    const scalar magGa = magGradAlpha[celli];

                    // True free-surface / gas→metal crossing only
                    const bool gasToMetal =
                        (prevA >= 0.0 && prevA < 0.4 && a >= 0.4);

                    const bool interfaceBand =
                        (a > 0.2 && a < 0.95)
                     && (magGa*cellLen > 1e-4)
                     && (prevA >= 0.0 && prevA < 0.5);

                    if (gasToMetal || interfaceBand)
                    {
                        localT = t;
                        localCell = celli;
                        localHit = p;

                        if (magGa > VSMALL)
                        {
                            localN = -gradAlpha[celli]/(magGa + VSMALL);
                            if ((localN & (-d)) < 0)
                            {
                                localN = -localN;
                            }
                        }
                        else
                        {
                            localN = -d;
                        }
                        break;
                    }

                    prevA = a;
                }
            }

            // ---- Global first-hit ownership (min path length) ----
            scalar globalT = localT;
            reduce(globalT, minOp<scalar>());

            // Tie-break: lowest rank among those with localT ≈ globalT
            label ownerProci = labelMax;
            if (localCell >= 0 && localT <= globalT*(1.0 + 1e-9) + SMALL)
            {
                ownerProci = myProci;
            }
            reduce(ownerProci, minOp<label>());

            const bool iOwnHit =
                (globalT < 0.5*GREAT)
             && (localCell >= 0)
             && (myProci == ownerProci)
             && (localT <= globalT*(1.0 + 1e-9) + SMALL);

            if (globalT >= 0.5*GREAT)
            {
                // Ray missed the free surface everywhere
                break;
            }

            // Reconstruct hit geometry on non-owner ranks for reflection
            // (owner broadcasts hit point, normal, and A via reductions)
            vector hitP = localHit;
            vector nHat = localN;
            if (!iOwnHit)
            {
                hitP = vector::zero;
                nHat = vector::zero;
            }
            reduce(hitP, sumOp<vector>());
            reduce(nHat, sumOp<vector>());
            if (mag(nHat) > SMALL)
            {
                nHat /= mag(nHat);
            }
            else
            {
                nHat = -d;
            }
            if (!iOwnHit)
            {
                // Non-owners zeroed hitP; after sum only owner contributed
                // hitP is correct globally.
            }

            // Angle-aware absorptivity: A(n) = A0 * μ_eff  (A = A0 at normal)
            const scalar mu = max((-d) & nHat, 0.0);
            const scalar muEff = max(mu, laserMinCosIncidence_);
            const scalar A =
                min(max(absorptivity_*muEff, 0.0), 0.99);

            const scalar Iabs = A*Iray;
            const scalar Irefl = (1.0 - A)*Iray;

            // ---- Deposit Iabs in metal (owner only) via Beer–Lambert ----
            // Integrates to ≤ Iabs; leftover stays in the hit cell so energy
            // is conserved on this rank (no silent loss, no extra creation).
            if (iOwnHit && Iabs > SMALL)
            {
                ++nHitsLocal;

                const scalar Labs =
                    (laserAbsorptionLength_ > SMALL)
                  ? laserAbsorptionLength_
                  : max(cellLenTyp, SMALL);
                const scalar muAbs = 1.0/Labs;

                // Collect local metal samples along the transmitted ray
                DynamicList<label> blCells(16);
                DynamicList<scalar> blW(16);
                scalar wSumBL = 0.0;

                // Always include the hit cell
                {
                    const scalar w0 = 1.0;
                    blCells.append(localCell);
                    blW.append(w0);
                    wSumBL += w0;
                }

                point pb = localHit;
                const label nBulk = max(label(4.0*Labs/ds), 3);
                for (label b = 0; b < nBulk; ++b)
                {
                    pb += d*ds;
                    const label cB = searchEng.findCell(pb);
                    if (cB < 0)
                    {
                        break;
                    }
                    const scalar aB = min(max(alpha1[cB], 0.0), 1.0);
                    if (aB < 0.15)
                    {
                        break;
                    }

                    // Weight ∝ α * exp(−μ s) * ds  (Beer–Lambert kernel)
                    const scalar s = (b + 1)*ds;
                    const scalar w =
                        aB*exp(-muAbs*s)*ds;
                    if (w > VSMALL)
                    {
                        blCells.append(cB);
                        blW.append(w);
                        wSumBL += w;
                    }
                }

                if (wSumBL > SMALL)
                {
                    forAll(blCells, k)
                    {
                        powerCell[blCells[k]] += Iabs*(blW[k]/wSumBL);
                    }
                }
                else
                {
                    powerCell[localCell] += Iabs;
                }
            }

            // ---- Multi-reflection: continue with Irefl ----
            if (Irefl <= SMALL || seg >= maxLaserReflections_)
            {
                break;
            }

            // Reflect on all ranks with the same (global) normal
            const scalar dn = d & nHat;
            d = d - 2.0*dn*nHat;
            d /= (mag(d) + VSMALL);
            p0 = hitP + d*(0.51*ds); // nudge off the surface
            Iray = Irefl;

            if (iOwnHit)
            {
                ++nReflectLocal;
            }
        }
    }

    // Power density shape from physical cell power [W]
    Qlaser_ = dimensionedScalar
    (
        "zero",
        dimensionSet(1, -1, -3, 0, 0, 0, 0),
        0
    );
    forAll(Qlaser_, celli)
    {
        if (powerCell[celli] > VSMALL)
        {
            Qlaser_[celli] = powerCell[celli]/mesh.V()[celli];
        }
    }

    scalar absorbedLocal = 0.0;
    forAll(powerCell, celli)
    {
        absorbedLocal += powerCell[celli];
    }
    reduce(absorbedLocal, sumOp<scalar>());

    // Physics: absorbed ≤ incident (equality only if A=1 everywhere).
    // Do NOT inflate to incidentPower — that would invent energy.
    // Thermal limiter may only reduce further (numerical safety).
    const scalar depositedBudget = max(absorbedLocal, 0.0);

    label nLaserLimited = 0;
    if (depositedBudget > SMALL)
    {
        redistributeLaserPower(CvEff, depositedBudget, dt, nLaserLimited);
    }
    else
    {
        Qlaser_ = dimensionedScalar
        (
            "zero",
            dimensionSet(1, -1, -3, 0, 0, 0, 0),
            0
        );
    }

    reduce(nHitsLocal, sumOp<label>());
    reduce(nReflectLocal, sumOp<label>());
    reduce(nSegments, sumOp<label>());

    if (debug || pimple.finalIter())
    {
        scalar appliedPower = 0.0;
        forAll(Qlaser_, celli)
        {
            appliedPower += Qlaser_[celli]*mesh.V()[celli];
        }
        reduce(appliedPower, sumOp<scalar>());

        Info<< "DED laser (rayTrace):"
            << " powerFactor=" << powerFactor
            << " incident=" << incidentPower
            << " rayAbsorbed=" << absorbedLocal
            << " applied=" << appliedPower
            << " abs/inc="
            << (incidentPower > SMALL ? absorbedLocal/incidentPower : 0)
            << " nRays=" << nRays
            << " nHits=" << nHitsLocal
            << " nReflect=" << nReflectLocal
            << " nLimited=" << nLaserLimited
            << " xLaser=" << xLaser
            << " z0=" << zOrigin
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


void Foam::solvers::compressibleVoF_DED::updatePowderSource
(
    const volScalarField& T
)
{
    // Continuum DED powder feed:
    //   ṁ_dep = η * powderFeedRate * ramp(t)
    //   distributed with Gaussian under the laser on geometric free surface
    //   that is hot/liquid enough to catch powder.
    // Mass: S1 = +mDotPowder  (added metal)
    // Energy: powderHeat ≈ mDot*[Cp*(Tp−T) − Lf_if_solid]

    mDotPowder_ = dimensionedScalar
    (
        "zero",
        dimensionSet(1, -3, -1, 0, 0, 0, 0),
        0
    );
    powderHeat_ = dimensionedScalar
    (
        "zero",
        dimensionSet(1, -1, -3, 0, 0, 0, 0),
        0
    );

    if (!powderEnabled_ || powderFeedRate_ <= SMALL)
    {
        mDotPowder_.correctBoundaryConditions();
        powderHeat_.correctBoundaryConditions();
        return;
    }

    scalar powerFactor = 1.0;
    if (laserRampTime_ > SMALL)
    {
        powerFactor = min(1.0, runTime.value()/laserRampTime_);
    }

    const scalar depositedMassRate =
        powderCatchEfficiency_*powderFeedRate_*powerFactor;

    if (depositedMassRate <= SMALL)
    {
        mDotPowder_.correctBoundaryConditions();
        powderHeat_.correctBoundaryConditions();
        return;
    }

    const scalar xLaser = xStart_ + laserSpeed_*runTime.value();
    const scalar rP = max(powderRadius_, SMALL);
    const scalar sigma = max(rP/3.0, SMALL);

    const volVectorField gradAlpha(fvc::grad(alpha1));
    const volScalarField magGradAlpha(mag(gradAlpha));

    // Geometric free-surface cells under powder Gaussian
    boolList isFreeSurface(mesh.nCells(), false);

    forAll(alpha1, celli)
    {
        const scalar a = min(max(alpha1[celli], 0.0), 1.0);
        const scalar magGa = magGradAlpha[celli];
        const scalar cellLen = cbrt(mesh.V()[celli]);

        if (a > 0.02 && a < 0.98 && magGa*cellLen > 1e-6)
        {
            isFreeSurface[celli] = true;
        }
    }

    {
        const labelUList& own = mesh.owner();
        const labelUList& nei = mesh.neighbour();

        forAll(own, facei)
        {
            const label c0 = own[facei];
            const label c1 = nei[facei];
            const scalar a0 = min(max(alpha1[c0], 0.0), 1.0);
            const scalar a1 = min(max(alpha1[c1], 0.0), 1.0);

            if ((a0 > 0.5 && a1 < 0.5) || (a1 > 0.5 && a0 < 0.5))
            {
                isFreeSurface[c0] = true;
                isFreeSurface[c1] = true;
            }
        }

        forAll(mesh.boundary(), patchi)
        {
            const fvPatch& patch = mesh.boundary()[patchi];
            if (!patch.coupled())
            {
                continue;
            }

            const fvPatchScalarField& alphaPf =
                alpha1.boundaryField()[patchi];
            const labelUList& fc = patch.faceCells();

            forAll(fc, facei)
            {
                const label c0 = fc[facei];
                const scalar a0 = min(max(alpha1[c0], 0.0), 1.0);
                const scalar aN = min(max(alphaPf[facei], 0.0), 1.0);

                if ((a0 > 0.5 && aN < 0.5) || (a0 < 0.5 && aN > 0.5))
                {
                    isFreeSurface[c0] = true;
                }
            }
        }
    }

    volScalarField weight
    (
        IOobject
        (
            "powderWeight",
            runTime.name(),
            mesh,
            IOobject::NO_READ,
            IOobject::NO_WRITE
        ),
        mesh,
        dimensionedScalar("zero", dimless/dimLength, 0)
    );

    label nCatchCells = 0;

    forAll(weight, celli)
    {
        if (!isFreeSurface[celli])
        {
            continue;
        }

        const scalar a = min(max(alpha1[celli], 0.0), 1.0);
        if (a < 0.05)
        {
            continue;
        }

        const point& C = mesh.C()[celli];
        const scalar dx = C.x() - xLaser;
        const scalar dy = C.y() - yStart_;
        const scalar G = exp(-(dx*dx + dy*dy)/(2.0*sigma*sigma));

        if (G < 1e-6)
        {
            continue;
        }

        // Catch only on liquid / mushy free surface under the beam
        const scalar fL =
            max(liquidFraction(T[celli]), fLiquid_[celli]);
        if (fL < powderMinLiquidFraction_ && T[celli] < solidusTemperature_)
        {
            continue;
        }

        const scalar magGa = magGradAlpha[celli];
        const scalar cellLen = cbrt(mesh.V()[celli]);

        scalar delta = max(interfaceDelta_[celli], 2.0*a*(1.0 - a)*magGa);
        if (a > 0.02 && a < 0.98)
        {
            delta = max(delta, magGa);
        }
        if (delta <= VSMALL || a >= 0.98)
        {
            delta = max(delta, 1.0/max(cellLen, SMALL));
        }

        // Prefer liquid fraction for catch distribution
        const scalar catchF = max(fL, powderMinLiquidFraction_);
        weight[celli] = G*delta*catchF;
        ++nCatchCells;
    }
    reduce(nCatchCells, sumOp<label>());

    scalar wIntegral = 0.0;
    forAll(weight, celli)
    {
        wIntegral += weight[celli]*mesh.V()[celli];
    }
    reduce(wIntegral, sumOp<scalar>());

    // Fallback: any metal free-surface under beam (ignore liquid gate)
    if (wIntegral < SMALL)
    {
        weight = dimensionedScalar("zero", dimless/dimLength, 0);
        nCatchCells = 0;

        forAll(weight, celli)
        {
            if (!isFreeSurface[celli])
            {
                continue;
            }

            const scalar a = min(max(alpha1[celli], 0.0), 1.0);
            if (a < 0.05)
            {
                continue;
            }

            const point& C = mesh.C()[celli];
            const scalar dx = C.x() - xLaser;
            const scalar dy = C.y() - yStart_;
            const scalar G = exp(-(dx*dx + dy*dy)/(2.0*sigma*sigma));
            if (G < 1e-6)
            {
                continue;
            }

            const scalar cellLen = cbrt(mesh.V()[celli]);
            weight[celli] = G/max(cellLen, SMALL);
            ++nCatchCells;
        }
        reduce(nCatchCells, sumOp<label>());

        wIntegral = 0.0;
        forAll(weight, celli)
        {
            wIntegral += weight[celli]*mesh.V()[celli];
        }
        reduce(wIntegral, sumOp<scalar>());
    }

    if (wIntegral <= SMALL)
    {
        mDotPowder_.correctBoundaryConditions();
        powderHeat_.correctBoundaryConditions();

        if (debug || pimple.finalIter())
        {
            Info<< "DED powder: depositedMassRate=" << depositedMassRate
                << " applied=0 (no catch surface) nCatch=0"
                << endl;
        }
        return;
    }

    const scalar scale = depositedMassRate/wIntegral;
    const volScalarField& Cv1 = mixture.thermo1().Cv();

    forAll(mDotPowder_, celli)
    {
        if (weight[celli] <= VSMALL)
        {
            continue;
        }

        const scalar md = scale*weight[celli]; // kg/m3/s
        mDotPowder_[celli] = md;

        // Enthalpy of added metal relative to local mixture temperature.
        // q > 0 heats cell if Tp > T; cold powder cools (q < 0).
        // Solid powder also requires melting energy Lf.
        const scalar Tc = T[celli];
        const scalar Cp = max(Cv1[celli], SMALL);
        scalar q = md*Cp*(powderTemperature_ - Tc);

        if (powderSolid_ && powderTemperature_ < solidusTemperature_)
        {
            // Cost to melt captured solid powder into the pool
            q -= md*Lf_.value();
        }

        powderHeat_[celli] = q;
    }

    mDotPowder_.correctBoundaryConditions();
    powderHeat_.correctBoundaryConditions();

    if (debug || pimple.finalIter())
    {
        scalar massRate = 0.0;
        scalar heatPower = 0.0;
        forAll(mDotPowder_, celli)
        {
            massRate += mDotPowder_[celli]*mesh.V()[celli];
            heatPower += powderHeat_[celli]*mesh.V()[celli];
        }
        reduce(massRate, sumOp<scalar>());
        reduce(heatPower, sumOp<scalar>());

        Info<< "DED powder:"
            << " feed=" << powderFeedRate_
            << " eta=" << powderCatchEfficiency_
            << " deposited=" << depositedMassRate
            << " applied=" << massRate
            << " heatPower=" << heatPower
            << " nCatch=" << nCatchCells
            << " Tp=" << powderTemperature_
            << " xLaser=" << xLaser
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
