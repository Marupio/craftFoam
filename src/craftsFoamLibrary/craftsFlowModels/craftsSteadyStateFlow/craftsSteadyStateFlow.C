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

#include "craftsSteadyStateFlow.H"
#include "addToRunTimeSelectionTable.H"

// * * * * * * * * * * * * * * Static Data Members * * * * * * * * * * * * * //

namespace Foam
{
    defineTypeNameAndDebug(craftsSteadyStateFlow, 0);

    addToRunTimeSelectionTable
    (
        craftsFlow,
        craftsSteadyStateFlow,
        nameConstructor
    );
}


// * * * * * * * * * * * * * * * * Constructors  * * * * * * * * * * * * * * //

Foam::craftsSteadyStateFlow::craftsSteadyStateFlow
(
    volVectorField& U,
    surfaceScalarField& phi,
    volScalarField& p,
    admReactionReader& model
)
:
    craftsFlow(U, phi, p, model)
{}


// * * * * * * * * * * * * * * Member Functions  * * * * * * * * * * * * * * //

Foam::label Foam::craftsSteadyStateFlow::step()
{
    // do nothing
    return 0;
}


void Foam::craftsSteadyStateFlow::storeFlowSolution()
{
    // do nothing
}


void Foam::craftsSteadyStateFlow::retagFlowSolution()
{
    // do nothing
}


void Foam::craftsSteadyStateFlow::calculateScales()
{
    if (outputFlowErrorScales_)
    {
        Info << "flowErrorScales: steadyState has no error scales" << endl;
    }
}


Foam::scalar Foam::craftsSteadyStateFlow::testConvergence() const
{
    if (outputFlowResiduals_)
    {
        Info << "flowResiduals: steadyState has no residuals" << endl;
    }
    return -VGREAT;
}


Foam::scalar Foam::craftsSteadyStateFlow::calculateNextBestDeltaT() const
{
    if (outputFlowTimestepEstimate_)
    {
        Info << "flowTimestepEstimate: steadyState has no time estimate"
            << endl;
    }
    return VGREAT;
}


void Foam::craftsSteadyStateFlow::saveState(const label slot)
{
    // do nothing
}
            

void Foam::craftsSteadyStateFlow::clearState(const label slot)
{
    // do nothing
}

            
void Foam::craftsSteadyStateFlow::loadState(const label slot)
{
    // do nothing
}
            

Foam::label Foam::craftsSteadyStateFlow::nStates() const
{
    return -1;
}


bool Foam::craftsSteadyStateFlow::validState(const label slot) const
{
    return true;
}

// ************************************************************************* //
