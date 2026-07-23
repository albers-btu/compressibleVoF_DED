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
#include "fvcDiv.H"

// * * * * * * * * * * * * * * Member Functions  * * * * * * * * * * * * * * //

void Foam::solvers::compressibleVoF_DED::alphaSuSp
(
    tmp<volScalarField::Internal>& tSu,
    tmp<volScalarField::Internal>& tSp
)
{
    const dimensionedScalar Szero(vDot.dimensions(), 0);

    tSp = volScalarField::Internal::New("Sp", mesh, Szero);
    tSu = volScalarField::Internal::New("Su", mesh, Szero);

    volScalarField::Internal& Sp = tSp.ref();
    volScalarField::Internal& Su = tSu.ref();

    if (fvModels().addsSupToField(mixture.rho1().name()))
    {
        const volScalarField::Internal alpha2ByRho1(alpha2()/mixture.rho1()());
        const fvScalarMatrix alphaRho1Sup
        (
            fvModels().sourceProxy(alpha1, mixture.rho1(), alpha1)
        );

        Su += alpha2ByRho1*alphaRho1Sup.Su();
        Sp += alpha2ByRho1*alphaRho1Sup.Sp();
    }

    if (fvModels().addsSupToField(mixture.rho2().name()))
    {
        const volScalarField::Internal alpha1ByRho2(alpha1()/mixture.rho2()());
        const fvScalarMatrix alphaRho2Sup
        (
            fvModels().sourceProxy(alpha2, mixture.rho2(), alpha2)
        );

        Su -= alpha1ByRho2*(alphaRho2Sup.Su() + alphaRho2Sup.Sp());
        Sp += alpha1ByRho2*alphaRho2Sup.Sp();
    }

    // Phase-C mass-conserving evaporation: metal (1) → gas (2)
    // Mass rates: S1 = -mDot, S2 = +mDot
    // Converted to α1 sources like compressible VoF rho-source mapping:
    //   Su += (α2/ρ1)*S1 - (α1/ρ2)*S2
    //       = -mDot*(α2/ρ1 + α1/ρ2)
    if (massConservingEvaporation_)
    {
        const volScalarField::Internal& rho1 = mixture.rho1()();
        const volScalarField::Internal& rho2 = mixture.rho2()();

        forAll(Su, celli)
        {
            const scalar md = mDotEvap_[celli];
            if (md <= SMALL)
            {
                continue;
            }

            const scalar r1 = max(rho1[celli], SMALL);
            const scalar r2 = max(rho2[celli], SMALL);
            const scalar a1 = min(max(alpha1[celli], 0.0), 1.0);
            const scalar a2 = min(max(alpha2[celli], 0.0), 1.0);

            // Explicit α1 destruction (stabilised)
            Su[celli] -= md*(a2/r1 + a1/r2);
        }
    }

    forAll(vDot, celli)
    {
        if (vDot[celli] > 0.0)
        {
            Sp[celli] -=
                vDot[celli]/max(1.0 - alpha1[celli], vDotResidualAlpha);
            Su[celli] +=
                vDot[celli]/max(1.0 - alpha1[celli], vDotResidualAlpha);
        }
        else if (vDot[celli] < 0.0)
        {
            Sp[celli] += vDot[celli]/max(alpha1[celli], vDotResidualAlpha);
        }
    }
}


// ************************************************************************* //
