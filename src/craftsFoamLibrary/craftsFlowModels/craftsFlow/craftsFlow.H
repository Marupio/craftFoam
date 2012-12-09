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
    Foam::craftsFlow

Description
    Abstract base class for flow models in the anaerobic digester model, ADM1.
    
SourceFiles
    craftsFlowI.H
    craftsFlow.C
    newAdmFlow.C

Author
    David L. F. Gaden

\*---------------------------------------------------------------------------*/

#ifndef craftsFlow_H
#define craftsFlow_H

#include "admTime.H"
#include "fvCFD.H"

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

namespace Foam
{

// Forward declarations

class admReactionReader;


/*---------------------------------------------------------------------------*\
                           Class craftsFlow Declaration
\*---------------------------------------------------------------------------*/

class craftsFlow
{

public:

    //- Sub-stepping type enumeration
    enum enumSubSteppingType
    {
        OFF,
        FIXED_N_STEPS,
        FIXED_TIMESTEP,
        ADAPTIVE_TIMESTEP
    };


protected:

    //- Model references
    const admReactionReader& model_;
    const admTime& runTime_;
    const fvMesh& mesh_;
    
    //- Fluid dynamic data references from the global scope
    volVectorField& UGlobal_;
    surfaceScalarField& phiGlobal_;
    volScalarField& pGlobal_;

    //- Flow model settings subDictionary
    const dictionary& settingsDict_;
    
    //- Coarse solution save spot
    scalarField coarseSaveSpotPhi_;
    scalarField coarseSaveSpotP_;

    
    // Adaptive timestepping (ATS) variables
    //  Adaptive timestepping can work with:
    //      - reactions only
    //      - reactions and flow model
    //  This section determines this behaviour.

        //- True means ATS is enabled with the flow model
        Switch atsEnableWithFlowModel_;
        
        //- List of variables to omit in ATS calculations
        wordList atsIgnore_;
        
        //- True when these variables are being ignored
        bool atsIgnorePhi_;
        bool atsIgnoreP_;
        
        //- ATS tolerance
        scalar atsPhiTolerance_;
        scalar atsPTolerance_;

        //- ATS minimum error scales (user-defined)
        scalar atsPhiMinScale_;
        scalar atsPMinScale_;

        //- ATS scales (actual values)
        scalar atsPhiScale_;
        scalar atsPScale_;


    // SubStepping
    //  When enabled, the flow model takes nSubSteps for each step of the
    //  reaction model.  Number of steps must always be even.

        //- Sub-stepping behaviour type:
        //      OFF
        //      FIXED_N_STEPS
        //      FIXED_TIMESTEP
        //      ADAPTIVE_TIMESTEP
        word subSteppingTypeWord_;
        enumSubSteppingType subSteppingType_;

        //- True when sub-stepping is enabled
        bool subStepping_;

        //- Maximum flow timestep
        //  Used with:
        //      FIXED_TIMESTEP
        //      ADAPTIVE_TIMESTEP
        scalar ssMaxDeltaT_;
        
        //- Minimum flow timestep
        //  Used with:
        //      FIXED_TIMESTEP
        //      ADAPTIVE_TIMESTEP
        scalar ssMinDeltaT_;
    
        //- Maximum number of subSteps
        //  Used with:
        //      FIXED_TIMESTEP
        //      ADAPTIVE_TIMESTEP
        label ssMaxNSubSteps_;
        
        //- Minimum number of subSteps
        //  Used with:
        //      FIXED_TIMESTEP
        //      ADAPTIVE_TIMESTEP
        label ssMinNSubSteps_;

        //- Initial number of subSteps
        //  Used with:
        //      FIXED_N_STEPS
        //      FIXED_TIMESTEP
        //      ADAPTIVE_TIMESTEP
        label ssInitialNSubSteps_;
        
        //- Current number of subSteps
        //  Used with all subStepping types
        mutable label ssCurrentNSubSteps_;

        //- Sub-stepping - active fluid dynamic data
        autoPtr<volVectorField> USubPtr_;
        autoPtr<surfaceScalarField> phiSubPtr_;
        autoPtr<volScalarField> pSubPtr_;

        //- Sub-stepping - temporary stored fluid dynamic data
        autoPtr<volVectorField> UTmpPtr_;
        autoPtr<surfaceScalarField> phiTmpPtr_;
        autoPtr<volScalarField> pTmpPtr_;
        
        //- Sub-stepping - fine save spots
        autoPtr<scalarField> fineSaveSpotPhiPtr_;
        autoPtr<scalarField> fineSaveSpotPPtr_;

        //- Reference to active flow field variables:
        //      if not subStepping - U_, phi_, p_; and
        //      if subStepping - USubPtr(), phiSubPtr(), pSubPtr()
        volVectorField& UActive_;
        surfaceScalarField& phiActive_;
        volScalarField& pActive_;


    // Output flags
    // These are yes/no settings that change what is output to the console
    // during a run.  They are contained in admSettingsDict.  Derived classes
    // can add to this subdictionary.

        //- Flow solver performance (default true)
        bool outputFlowSolverPerformance_;

        //- Flow residuals (default true)
        bool outputFlowResiduals_;

        //- Flow variable error scaling factors  (default false)
        bool outputFlowErrorScales_;

        //- Flow timestep estimate (default true)
        bool outputFlowTimestepEstimate_;

        //- Flow continuity errors (default true)
        bool outputFlowContinuityErrors_;


private:

    // Private member functions
    
        //- Read the flowModel subDictionary of the settings dictionary
        void readDict();
    
        //- Disallow default bitwise copy construct and assignment
        void operator=(const craftsFlow&);
        craftsFlow(const craftsFlow&);

public:

    //- Runtime type information
    TypeName("craftsFlow");


    // Declare run-time constructor selection tables

        declareRunTimeSelectionTable
        (
            autoPtr,
            craftsFlow,
            nameConstructor,
            (
                volVectorField& U,
                surfaceScalarField& phi,
                volScalarField& p,
                admReactionReader& model
            ),
            (U, phi, p, model)
        );


    // Constructors

        //- Construct from components
        craftsFlow
        (
            volVectorField& U,
            surfaceScalarField& phi,
            volScalarField& p,
            admReactionReader& model
        );

    // Selectors


        //- Return a reference to a new craftsFlow model
        //  *** Located in newAdmFlow.C ***
        static autoPtr<craftsFlow> New
        (
            volVectorField& U,
            surfaceScalarField& phi,
            volScalarField& p,
            admReactionReader& model,
            const word& flowModelName
        );


    // Destructor

        virtual ~craftsFlow()
        {}


    // Member functions


        // Access

            //- Main model, runTime and mesh
            inline const admReactionReader& model() const;
            inline const admTime& runTime() const;
            inline const fvMesh& mesh() const;
        
            //- Velocity (global)
            inline const volVectorField& U() const;
            inline volVectorField& U();
            
            //- Mass flux / rho (global)
            inline const surfaceScalarField& phi() const;
            inline surfaceScalarField& phi();

            //- Pressure (global)
            inline const volScalarField& p() const;
            inline volScalarField& p();


            // Adaptive timestepping

                //- True means ATS is enabled with the flow model
                inline bool atsEnableWithFlowModel() const;
                
                //- List of variables to omit in ATS calculations
                inline const wordList& atsIgnore() const;
                
                //- True when these variables are being ignored
                inline bool atsIgnorePhi() const;
                inline bool atsIgnoreP() const;
                
                //- ATS tolerance
                inline const scalar& atsPhiTolerance() const;
                inline const scalar& atsPTolerance() const;

                //- ATS minimum error scales (user-defined)
                inline const scalar& atsPhiMinScale() const;
                inline const scalar& atsPMinScale() const;

                //- ATS scales (actual values)
                inline const scalar& atsPhiScale() const;
                inline const scalar& atsPScale() const;


            // SubStepping

                //- Returns user-selected behaviour of subStepping:
                //      OFF
                //      FIXED_N_STEPS
                //      FIXED_TIMESTEP
                //      ADAPTIVE_TIMESTEP
                inline const enumSubSteppingType& subSteppingType() const;

                //- Returns true if subStepping is enabled
                inline const bool& subStepping() const;
                
                //- Current number of subSteps
                inline const label& nSubSteps() const;

                //- Set a new number of subSteps, applies max and min, even
                //  only, and returns the actual value set.
                //  const because ugly
                label setNSubSteps(label newNSubSteps) const;


            //- Output flags
            inline const bool& outputFlowSolverPerformance() const;
            inline const bool& outputFlowResiduals() const;
            inline const bool& outputFlowTimestepEstimate() const;


        // Flow model functions

            //- Perform solution for one timestep
            virtual label step();

            //- Store the current state as the coarse solution, for later
            //  comparison against the fine solution
            virtual void storeCoarseSolution();

            //- Store the current state as the fine solution.  Only valid
            //  when sub-stepping; otherwise this will already be in memory.
            virtual void storeFineSolution();
            
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

#include "craftsFlowI.H"

#endif

// ************************************************************************* //
