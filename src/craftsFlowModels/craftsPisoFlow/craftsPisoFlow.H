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
    Foam::craftsPisoFlow

Description
    An incompressible PISO flow solver with turbulence and optional source term
    for the anaerobic digestion model ADM-MDA.
    
SourceFiles
    craftsPisoFlow.C

Author
    David L. F. Gaden

\*---------------------------------------------------------------------------*/

#ifndef craftsPisoFlow_H
#define craftsPisoFlow_H

#include "craftsFlow.H"
#include "singlePhaseTransportModel.H"
#include "turbulenceModel.H"
#include "admReactionReader.H"

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

namespace Foam
{

// Forward declarations
//class admCoefficient;

/*---------------------------------------------------------------------------*\
                         Class craftsPisoFlow Declaration
\*---------------------------------------------------------------------------*/

class craftsPisoFlow
:
    public craftsFlow
{

private:

    // Member data that is normally at the global scope in a solver

        //- Transport model
        singlePhaseTransportModel transport_;

        //- Turbulence model
        autoPtr<incompressible::turbulenceModel> turbulence_;
        
        //- Algorithm data
        autoPtr<fvVectorMatrix> UEqnPtr_;

        //- PISO dictionary settings
        const dictionary& pisoDict_;
        int nCorr_;
        int nNonOrthCorr_;
        bool momentumPredictor_;
        label pRefCell_;
        scalar pRefValue_;
        
        //- Continuity error data members
        scalar sumLocalContErr;
        scalar globalContErr;
        scalar cumulativeContErr;

    // Source term data

        //- Source term names
        const word sourceTermXName_;
        const word sourceTermYName_;
        const word sourceTermZName_;
        
        //- A dimensioned coefficient for zero source term
        const admZeroCoefficient zeroSource_;
        
        //- Source term coefficients
        const admCoefficient& sourceTermX_;
        const admCoefficient& sourceTermY_;
        const admCoefficient& sourceTermZ_;


    // Save spots
    
        //- Turbulence model
        autoPtr<symmTensorField> coarseSaveSpotRPtr_;
        autoPtr<symmTensorField> fineSaveSpotRPtr_;


    // Adaptive timestepping variables
    
        //- True when these variables are being ignored
        bool atsIgnoreR_;

        //- Minimum error scaling factor
        symmTensor atsRMinScale_;

        //- Turbulence error scaling data
        symmTensor atsRScale_;
        
        //- Convergence criterion
        symmTensor atsRTolerance_;

    
    // Private member functions
    
        //- Read settings in the dictionary
        void readDict();
    
        //- Solve the U equation
        void solveUEquation();

        //- Solve the P equation
        void solvePEquation(int corr);
        
        //- Calculate continuityErrs
        void continuityErrs();
    
public:

    //- Runtime type information
    TypeName("piso");


    // Constructors

        //- Construct from dictionary, given name
        craftsPisoFlow
        (
            volVectorField& U,
            surfaceScalarField& phi,
            volScalarField& p,
            admReactionReader& model
        );

    // Destructor

        virtual ~craftsPisoFlow()
        {}


    // Member functions

        // Access
        
            //- Transport model
            inline const singlePhaseTransportModel& transport() const;
            inline singlePhaseTransportModel& transport();

            //- Turbulence model
            inline const incompressible::turbulenceModel& turbulence() const;
            inline incompressible::turbulenceModel& turbulence();

        // Flow model functions

            //- Perform solution for one timestep
            virtual label step();

            //- Store the current state as the coarse solution, for later
            //  comparison against the fine solution
            virtual void storeCoarseSolution();

            //- Store the current state as the fine solution.  Only valid
            //  when sub-stepping; otherwise this will already be in memory.
            virtual void storeFineSolution();

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

#include "craftsPisoFlowI.H"

#endif

// ************************************************************************* //
