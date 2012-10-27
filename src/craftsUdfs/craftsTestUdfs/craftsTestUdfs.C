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

#include "craftsTestUdfs.H"
#include "addToRunTimeSelectionTable.H"
#include "craftsModel.H"

// * * * * * * * * * * * * * * * * Constructors  * * * * * * * * * * * * * * //

template<int matrixSize>
Foam::craftsTestUdfs<matrixSize>::craftsTestUdfs
(
    craftsModel<matrixSize>& model,
    const word& hooksDictName
)
:
    craftsUdfs<matrixSize>(model, hooksDictName)
{}


// * * * * * * * * * * * * * * Member Functions  * * * * * * * * * * * * * * //

template<int matrixSize>
void Foam::craftsTestUdfs<matrixSize>::initializeTimestep()
{}


template<int matrixSize>
void Foam::craftsTestUdfs<matrixSize>::initializeImplicitLoop()
{}


template<int matrixSize>
Foam::label Foam::craftsTestUdfs<matrixSize>::implicitLoop()
{
    admStandardVariable& messMeUp(admVars.lookupStandard("S_1"));
    const admStandardVariable& reference(admVars.lookupStandard("S_0"));
    messMeUp(0) = reference(0) * 2;
Info << reference(0) << "," << messMeUp(0) << endl;
    return 0;
}


template<int matrixSize>
void Foam::craftsTestUdfs<matrixSize>::finalizeImplicitLoop()
{}


template<int matrixSize>
void Foam::craftsTestUdfs<matrixSize>::finalizeTimestep()
{}


template<int matrixSize>
void Foam::craftsTestUdfs<matrixSize>::saveState(const label slot)
{}


template<int matrixSize>
void Foam::craftsTestUdfs<matrixSize>::clearState(const label slot)
{}


template<int matrixSize>
void Foam::craftsTestUdfs<matrixSize>::loadState(const label slot)
{}


template<int matrixSize>
Foam::label Foam::craftsTestUdfs<matrixSize>::nStates() const
{
    return -1;
}


template<int matrixSize>
bool Foam::craftsTestUdfs<matrixSize>::validState(const label slot) const
{
    return true;
}

// ************************************************************************* //
