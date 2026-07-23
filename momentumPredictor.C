/*---------------------------------------------------------------------------*\
  =========                 |
  \\      /  F ield         | OpenFOAM: The Open Source CFD Toolbox
   \\    /   O peration     | Website:  https://openfoam.org
    \\  /    A nd           | Copyright (C) 2023 OpenFOAM Foundation
     \\/     M anipulation  |
-------------------------------------------------------------------------------
License
    This file is part of OpenFOAM.

    OpenFOAM is free software: you can redistribute it and/or modify it
    under the terms of the GNU General Public License as published by
    the Free Software Foundation, either version 3 of the License, or
    (at your option) any later version.

    OpenFOAM is distributed in the hope that it will be useful, but WITHOUT
    ANY WARRANTY; without even the implied warranty of MERCHANTABILITY or
    FITNESS FOR A PARTICULAR PURPOSE.  See the GNU General Public License
    for more details.

    You should have received a copy of the GNU General Public License
    along with OpenFOAM.  If not, see <http://www.gnu.org/licenses/>.

\*---------------------------------------------------------------------------*/

#include "compressibleVoF_DED.H"
#include "fvmSup.H"
#include "fvmDiv.H"
#include "fvcSnGrad.H"
#include "fvcReconstruct.H"

// * * * * * * * * * * * * * * Member Functions  * * * * * * * * * * * * * * //

Foam::tmp<Foam::fvVectorMatrix> Foam::solvers::compressibleVoF_DED::divDevTau
(
    volVectorField& U
)
{
    return
        momentumTransport.divDevTau(U)
      - fvm::Sp(contErr1() + contErr2(), U);
}


void Foam::solvers::compressibleVoF_DED::momentumPredictor()
{
    volVectorField& U = U_;

    // Recompute mushy resistance with latest fLiquid (updated each outer
    // corrector after thermophysicalPredictor) and current α / sponge.
    updateMushyAndSpongeResistance();

    tUEqn =
    (
        fvm::ddt(rho, U)
      + fvm::div(rhoPhi, U)
      + MRF.DDt(rho, U)
      + divDevTau(U)
      + fvm::Sp(mushyResistance_, U)
      + fvm::Sp(gasSpongeResistance_, U)
     ==
        fvModels().source(rho, U)
    );
    fvVectorMatrix& UEqn = tUEqn.ref();

    UEqn.relax();

    fvConstraints().constrain(UEqn);

    if (pimple.momentumPredictor())
    {
        // Normal CSF forces on faces; Marangoni as volumetric body force
        solve
        (
            UEqn
         ==
            marangoniForce_
          + fvc::reconstruct
            (
                (
                    surfaceTensionForce()
                  + fvc::interpolate(pRecoil_)*fvc::snGrad(alpha1)
                  - buoyancy.ghf*fvc::snGrad(rho)
                  - fvc::snGrad(p_rgh)
                ) * mesh.magSf()
            )
        );

        fvConstraints().constrain(U);

        K = 0.5*magSqr(U);
    }

    if (pimple.finalIter())
    {
        Info<< "DED momentum:"
            << " max(mushyResistance)=" << gMax(mushyResistance_)
            << " max(gasSpongeResistance)=" << gMax(gasSpongeResistance_)
            << endl;
    }
}


// ************************************************************************* //
