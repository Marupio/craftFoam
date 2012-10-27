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

Class
    Foam::craftsSteadyStateFlow

Description
    An inactive flow model: assumes the initial velocity is constant.
    
SourceFiles
    craftsSteadyStateFlow.C

Author
    David L. F. Gaden

\*---------------------------------------------------------------------------*/

#ifndef craftsSteadyStateFlow_H
#define craftsSteadyStateFlow_H

#include "admFlow.H"
#include "admReactionReader.H"

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

namespace Foam
{

// Forward declarations
//class admCoefficient;


/*---------------------------------------------------------------------------*\
                         Class craftsSteadyStateFlow Declaration
\*---------------------------------------------------------------------------*/

class craftsSteadyStateFlow
:
    public admFlow
{

public:

    //- Runtime type information
    TypeName("steadyState");


    // Constructors

        //- Construct from dictionary, given name
        craftsSteadyStateFlow
        (
            volVectorField& U,
            surfaceScalarField& phi,
            volScalarField& p,
            admReactionReader& model
        );

    // Destructor

        virtual ~craftsSteadyStateFlow()
        {}


    // Member functions

        // Flow model functions

            //- Perform solution for one timestep
            virtual label step();

            //- Store the current flow field into a temporary store area for
            //  later retagging.  Allows the algorithm to take the midpoint
            //  flow solution for the reactions when using sub-stepping
            virtual void storeFlowSolution();
            
            //- Moves the temporary stored flow solution into the globally
            //  defined flow variables
            virtual void retagFlowSolution();

            //- Calculate the error scaling factors for each variable
            virtual void calculateScales();

            //- Compare the current state against the coarse solution
            virtual scalar testConvergence() const;

            //- Estimate the next best delta T to use
            virtual scalar calculateNextBestDeltaT() const;


        // State control
        //  Since OpenFOAM does not easily support undoing a timestep, or
        //  stepping back several timesteps, the save state framework must be
        //  used.
        
            //- Save data to backup variables
            virtual void saveState(const label slot = 0);
            
            //- Erase backup variables
            virtual void clearState(const label slot = 0);
            
            //- Load backup variables into active variables
            virtual void loadState(const label slot = 0);
            
            //- Return number of spots that are occupied
            //  A value of -1 means it doesn't track nStates
            virtual label nStates() const;

            //- Returns true if this is a loadable save state slot
            virtual bool validState(const label slot) const;
};


// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

} // End namespace Foam

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

#endif

// ************************************************************************* //
