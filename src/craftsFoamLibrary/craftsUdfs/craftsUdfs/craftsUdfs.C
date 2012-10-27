/*---------------------------------------------------------------------------*\
  =========                 |
  \\      /  F ield         | OpenFOAM: The Open Source CFD Toolbox
   \\    /   O peration     |
    \\  /    A nd           | Copyright held by original author
     \\/     M anipulation  |
-------------------------------------------------------------------------------
License
    This file is part of OpenFOAM.

    OpenFOAM is free software; you can redistribute it and/or modify it
    under the terms of the GNU General Public License as published by the
    Free Software Foundation; either version 2 of the License, or (at your
    option) any later version.

    OpenFOAM is distributed in the hope that it will be useful, but WITHOUT
    ANY WARRANTY; without even the implied warranty of MERCHANTABILITY or
    FITNESS FOR A PARTICULAR PURPOSE.  See the GNU General Public License
    for more details.

    You should have received a copy of the GNU General Public License
    along with OpenFOAM; if not, write to the Free Software Foundation,
    Inc., 51 Franklin St, Fifth Floor, Boston, MA 02110-1301 USA

\*---------------------------------------------------------------------------*/

#include "craftsUdfs.H"
#include "craftsModel.H"

// * * * * * * * * * * * * * * * * Constructors  * * * * * * * * * * * * * * //

template<int matrixSize>
Foam::craftsUdfs<matrixSize>::craftsUdfs
(
    craftsModel<matrixSize>& model,
    const word& hooksDictName
)
:
    model(model),
    runTime(model.runTime()),
    mesh(model.mesh()),
    admVars(model.admVars()),
    admCoeffs(model.admCoeffs()),
    admReacs(model.admReacs()),
    eqns(runTime.eqns()),
    
    admHooksDict_
    (
        IOobject
        (
            hooksDictName,
            runTime.constant(),
            runTime,
            IOobject::MUST_READ,
            IOobject::NO_WRITE
        )
    )
{}


template<int matrixSize>
Foam::craftsUdfs<matrixSize>::craftsUdfs
(
    const craftsUdfs& admc
)
:
    model(admc.model),
    runTime(admc.runTime),
    mesh(admc.mesh),
    admVars(admc.admVars),
    admCoeffs(admc.admCoeffs),
    admReacs(admc.admReacs),
    eqns(admc.eqns),
    admHooksDict_(admc.admHooksDict_)
{}

// * * * * * * * * * * * * * * Member Functions  * * * * * * * * * * * * * * //

template<int matrixSize>
void Foam::craftsUdfs<matrixSize>::applyVariableLimits()
{
        admVars.applyStandardLimits();
        admVars.applyImplicitLimits();
}


// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

#   include "newCraftsUdfs.C"

// ************************************************************************* //
