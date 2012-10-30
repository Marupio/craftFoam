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

#include "craftsGasModel.H"

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

namespace Foam
{

// * * * * * * * * * * * * * * * Member Functions  * * * * * * * * * * * * * //

template<int matrixSize>
inline const IOdictionary& craftsGasModel<matrixSize>::craftsGasModelDict()
const
{
    return craftsGasModelDict_;
}


template<int matrixSize>
inline const word& craftsGasModel<matrixSize>::surfacePatchName() const
{
    return surfacePatchName_;
}


template<int matrixSize>
inline const fvPatch& craftsGasModel<matrixSize>::surfacePatch() const
{
    return mesh_.boundary()[surfacePatchIndex_];
}


template<int matrixSize>
inline const UPtrList<admVariable>& craftsGasModel<matrixSize>::vars() const
{
    return vars_;
}


template<int matrixSize>
inline const UPtrList<admVariable>& craftsGasModel<matrixSize>::bcVars() const
{
    return bcVars_;
}


template<int matrixSize>
inline const IOdictionary craftsGasModel<matrixSize>::outputDict() const
{
    return outputDict_;
}


// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

} // End namespace Foam

// ************************************************************************* //
