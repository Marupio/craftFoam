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

#include "craftsEmptyUdfs.H"
#include "addToRunTimeSelectionTable.H"
#include "craftsModel.H"

// * * * * * * * * * * * * * * * * Constructors  * * * * * * * * * * * * * * //

template<int matrixSize>
Foam::craftsEmptyUdfs<matrixSize>::craftsEmptyUdfs
(
    craftsModel<matrixSize>& model,
    const word& hooksDictName
)
:
    craftsUdfs<matrixSize>(model, hooksDictName)
{}


// * * * * * * * * * * * * * * Member Functions  * * * * * * * * * * * * * * //

template<int matrixSize>
void Foam::craftsEmptyUdfs<matrixSize>::initializeTimestep()
{}


/*template<int matrixSize>
void Foam::craftsEmptyUdfs<matrixSize>::applyVariableLimits()
{}*/


template<int matrixSize>
void Foam::craftsEmptyUdfs<matrixSize>::initializeImplicitLoop()
{}


template<int matrixSize>
Foam::label Foam::craftsEmptyUdfs<matrixSize>::implicitLoop()
{
    return 0;
}


template<int matrixSize>
void Foam::craftsEmptyUdfs<matrixSize>::finalizeImplicitLoop()
{}


template<int matrixSize>
void Foam::craftsEmptyUdfs<matrixSize>::finalizeTimestep()
{}


template<int matrixSize>
void Foam::craftsEmptyUdfs<matrixSize>::saveState(const label slot)
{}


template<int matrixSize>
void Foam::craftsEmptyUdfs<matrixSize>::clearState(const label slot)
{}


template<int matrixSize>
void Foam::craftsEmptyUdfs<matrixSize>::loadState(const label slot)
{}


template<int matrixSize>
Foam::label Foam::craftsEmptyUdfs<matrixSize>::nStates() const
{
    return -1;
}


template<int matrixSize>
bool Foam::craftsEmptyUdfs<matrixSize>::validState(const label slot) const
{
    return true;
}

// ************************************************************************* //
