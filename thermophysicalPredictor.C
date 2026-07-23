/*---------------------------------------------------------------------------*\
  =========                 |
  \\      /  F ield         | OpenFOAM: The Open Source CFD Toolbox
   \\    /   O peration     | Website:  https://openfoam.org
    \\  /    A nd           | Copyright (C) 2022-2023 OpenFOAM Foundation
     \\/     M anipulation  |
-------------------------------------------------------------------------------
License
    This file is part of OpenFOAM.

    Energy:
      1) stock compressible VoF T equation (transport + apparent Cp)
      2) operator-split DED sources (laser, radiation, evaporation)
         applied as explicit ΔT = q·Δt / (ρ Cv) so heating is local and
         power-conserving (avoids under-drive of correction-form matrix).

\*---------------------------------------------------------------------------*/

#include "compressibleVoF_DED.H"
#include "fvcMeshPhi.H"
#include "fvcDdt.H"
#include "fvcGrad.H"
#include "fvmDiv.H"
#include "fvmSup.H"
#include "fvmLaplacian.H"

// * * * * * * * * * * * * * * Member Functions  * * * * * * * * * * * * * * //

void Foam::solvers::compressibleVoF_DED::thermophysicalPredictor()
{
    const volScalarField& rho1(mixture.rho1());
    const volScalarField& rho2(mixture.rho2());

    const volScalarField& e1(mixture.thermo1().he());
    const volScalarField& e2(mixture.thermo2().he());

    const fvScalarMatrix e1Source(fvModels().source(alpha1, rho1, e1));
    const fvScalarMatrix e2Source(fvModels().source(alpha2, rho2, e2));

    volScalarField& T = mixture_.T();

    updateInterfaceDelta();

    const volScalarField& Cv1 = mixture.thermo1().Cv();
    const volScalarField& Cv2 = mixture.thermo2().Cv();

    const dimensionedScalar sigmaSB
    (
        "sigmaSB",
        dimensionSet(1, 0, -3, -4, 0, 0, 0),
        5.670374419e-8
    );

    const scalar dt = max(runTime.deltaTValue(), SMALL);

    // -------------------------------------------------------------------
    // Picard latent-heat / energy correctors
    // -------------------------------------------------------------------

    for (label Ecorr = 0; Ecorr < nLatentCorrectors_; ++Ecorr)
    {
        // Apparent heat capacity: Cv + Lf * dfL/dT
        volScalarField CvEff
        (
            IOobject
            (
                "CvEff",
                runTime.name(),
                mesh,
                IOobject::NO_READ,
                IOobject::NO_WRITE
            ),
            Cv1
        );

        forAll(CvEff, celli)
        {
            CvEff[celli] += Lf_.value()*dLiquidFraction_dT(T[celli]);
        }

        // Stock compressible VoF energy equation (transport + latent Cv only).
        // DED volumetric sources are operator-split once after Picard.
        fvScalarMatrix TEqn
        (
            correction
            (
                CvEff
               *(
                    fvm::ddt(alpha1, rho1, T) + fvm::div(alphaRhoPhi1, T)
                  - (
                        e1Source.hasDiag()
                      ? fvm::Sp(contErr1(), T) + fvm::Sp(e1Source.A(), T)
                      : fvm::Sp(contErr1(), T)
                    )
                )
              + mixture.thermo2().Cv()()
               *(
                    fvm::ddt(alpha2, rho2, T) + fvm::div(alphaRhoPhi2, T)
                  - (
                        e2Source.hasDiag()
                      ? fvm::Sp(contErr2(), T) + fvm::Sp(e2Source.A(), T)
                      : fvm::Sp(contErr2(), T)
                    )
                )
            )

          + fvc::ddt(alpha1, rho1, e1) + fvc::div(alphaRhoPhi1, e1)
          - contErr1()*e1
          + fvc::ddt(alpha2, rho2, e2) + fvc::div(alphaRhoPhi2, e2)
          - contErr2()*e2

          - fvm::laplacian(thermophysicalTransport.kappaEff(), T)

          + (
                mixture.totalInternalEnergy()
              ?
                fvc::div(fvc::absolute(phi, U), p)()()
              + (fvc::ddt(rho, K) + fvc::div(rhoPhi, K))()()
              - (U()&(fvModels().source(rho, U)&U)())
              - (contErr1() + contErr2())*K
              :
                p*fvc::div(fvc::absolute(phi, U))()()
            )
         ==
            (e1Source&e1)
          + (e2Source&e2)
        );

        TEqn.relax();
        fvConstraints().constrain(TEqn);
        TEqn.solve();
        fvConstraints().constrain(T);

        updateLiquidFraction(T);
    }

    // -------------------------------------------------------------------
    // Operator-split DED sources — only on final PIMPLE outer corrector
    // so a later TEqn.solve() cannot wipe the deposited heat.
    // Immediately correctThermo() so he = he(p,T) stores the new energy.
    // -------------------------------------------------------------------
    if (pimple.finalIter())
    {
        volScalarField CvEff
        (
            IOobject
            (
                "CvEffSplit",
                runTime.name(),
                mesh,
                IOobject::NO_READ,
                IOobject::NO_WRITE
            ),
            Cv1
        );

        forAll(CvEff, celli)
        {
            CvEff[celli] += Lf_.value()*dLiquidFraction_dT(T[celli]);
        }

        updateLaserSource(CvEff);
        updateEvaporation(T, CvEff, rho1);

        scalar maxDT = 0.0;
        label nHeated = 0;

        forAll(T, celli)
        {
            const scalar a = min(max(alpha1[celli], 0.0), 1.0);
            const scalar rhoCv =
                a*rho1[celli]*CvEff[celli]
              + (1.0 - a)*rho2[celli]*Cv2[celli];

            if (rhoCv < SMALL)
            {
                continue;
            }

            const scalar delta = interfaceDelta_[celli];
            const scalar Tc = max(T[celli], SMALL);
            const scalar qRad =
                (delta > VSMALL)
              ? emissivity_*sigmaSB.value()
               *(pow4(Tc) - pow4(Tamb_))*delta
              : 0.0;

            const scalar qNet =
                Qlaser_[celli]
              - qRad
              - evaporationSink_[celli];

            if (mag(qNet) < VSMALL)
            {
                continue;
            }

            const scalar dT = qNet*dt/rhoCv;
            const scalar dTlim =
                min(max(dT, -laserMaxDeltaT_), laserMaxDeltaT_);

            T[celli] += dTlim;
            maxDT = max(maxDT, mag(dTlim));
            ++nHeated;
        }

        reduce(maxDT, maxOp<scalar>());
        reduce(nHeated, sumOp<label>());

        forAll(T, celli)
        {
            T[celli] = max(T[celli], 200.0);
        }
        T.correctBoundaryConditions();

        updateLiquidFraction(T);

        // Freeze heated state into phase enthalpies before next time step
        mixture_.correctThermo();
        mixture_.correct();

        Info<< "DED heatSplit:"
            << " max|dT|=" << maxDT
            << " nHeated=" << nHeated
            << " max(T)=" << gMax(T)
            << " max(Qlaser)=" << gMax(Qlaser_)
            << endl;
    }
    else
    {
        // Non-final outer: still need thermo properties for pressure
        mixture_.correctThermo();
        mixture_.correct();
    }

    updateSurfaceForces(T);
    updateMushyAndSpongeResistance();

    scalar moltenVolume = 0.0;
    forAll(fLiquid_, celli)
    {
        if (fLiquid_[celli] > 0.5)
        {
            moltenVolume +=
                alpha1[celli]*fLiquid_[celli]*mesh.V()[celli];
        }
    }
    reduce(moltenVolume, sumOp<scalar>());

    if (pimple.finalIter())
    {
        Info<< "DED thermo:"
            << " max(T)=" << gMax(T)
            << " max(fL)=" << gMax(fLiquid_)
            << " moltenVolume=" << moltenVolume
            << " nLatentCorr=" << nLatentCorrectors_
            << endl;
    }
}


// ************************************************************************* //
