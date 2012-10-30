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
    Foam::craftsGasModel

Description
    A simple gas model for the anaerobic digester model, ADM1.  Calculates gas
    transfer through: diffusion across the surface, and bubble formation
    throughout the volume.  Intended to be called during the implicit routines.
    See the implementation of craftsGasModel::implicitLoop() for the overall
    pseudo-code for the behaviour of this object.
    
SourceFiles
    craftsGasModelI.H
    craftsGasModel.C

Author
    David L. F. Gaden

\*---------------------------------------------------------------------------*/

#ifndef craftsGasModel_H
#define craftsGasModel_H

#include "admCoefficientManager.H"
#include "admVariableManager.H"
#include "craftsGasFlowOde.H"

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

namespace Foam
{

// Forward declarations

template<int>
class craftsModel;

/*---------------------------------------------------------------------------*\
                        Class craftsGasModel Declaration
\*---------------------------------------------------------------------------*/

template<int matrixSize>
class craftsGasModel
{

    // Private member variables

        //- Quick references
        admCoefficientManager& admCoeffs_;
        admVariableManager& admVars_;
        const admTime& runTime_;
        const fvMesh& mesh_;
        craftsModel<matrixSize>& model_;

        //- Input dictionary (in constant)
        IOdictionary craftsGasModelDict_;

        // Dictionary settings (universal)

            //- If true, output stats to the console
            bool reportPerformance_;

            //- Name of odeSolver
            const word odeSolverName_;

            //- Maximum iterations for odeSolver
            label odeSolverMaxIter_;

            //- Accuracy of ODE solver
            const scalar eps_;
            
            //- Default mass convergence criteria
            const scalar defaultMassConvergence_;
            
            //- Default delta (change with previous iter) convergence criteria
            const scalar defaultDeltaConvergence_;

            //- Default near-zero mass error scale factor
            const scalar defaultNearZeroMassScale_;

            //- Default near-zero delta error scale factor
            const scalar defaultNearZeroDeltaScale_;

            //- Name of output dictionary
            const word outputDictName_;
            
            //- Name of gas transfer surface patch
            const word surfacePatchName_;
            
            //- Index of gas transfer surface patch
            const label surfacePatchIndex_;
            
            //- Gas volume coefficient - must be constant / uniform
            const admCoefficient& coeff_V_;

            //- R, gas constant coefficient - must be constant / uniform
            const admCoefficient& coeff_R_;
            
            //- Temperature variable
            const admVariable& var_T_;

            //- Overall gas pipe resistance coefficient - must be constant /
            //  uniform
            const admCoefficient& coeff_kp_;
            
            //- Logistic function sharpness - must be constant / uniform
            //  Used for smoothing zero-crossing cut-off for q_gas
            const admCoefficient& coeff_kGasCutOff_;
            
            //- Logistic function offset - must be constant / uniform
            //  Used for off-setting zero-crossing cut-off for q_gas
            const admCoefficient& coeff_pGasOverPressure_;
            
            //- Atmospheric pressure - must be uniform
            const admCoefficient& coeff_pAtm_;
            
            //- Bulk liquid vapour pressure
            const admCoefficient& coeff_pLiq_;
            
        // Dictionary settings (species)
        
            //- List of species names
            wordList names_;
            
            //- Associated variables
            UPtrList<admVariable> vars_;
            
            //- Variables that have the associated dissolvedGasTransfer
            //  boundary condition - a vars_ is derived, there must be another
            //  variable (and only one) that has a suitable boundary condition
            //  to accommodate the mass transfer - finding this variable is
            //  done automatically by searching Jacobian variables.
            UPtrList<admVariable> bcVars_;
            
            //- Holds the new internal fields of the liquid variables as they
            //  were at the start of a timestep.  Used for bubble formation -
            //  bubbles are relative to the initial value at a timestep.
            PtrList<scalarField> varLiquidStart_;
            
            // True when bubble formation is enabled (indexed by species)
            boolList formBubbles_;

            //- Associated conversion to moles (n from S = n Kh p) - must be
            //  constant / uniform
            UPtrList<const admCoefficient> coeffs_n_;
        
            //- Overall diffusion gas transfer coefficient
            UPtrList<const admCoefficient> coeffs_kla_;

            //- Associated Henry's constants
            UPtrList<const admCoefficient> coeffs_kh_;
            
            //- Associated bubble transfer coefficient
            UPtrList<const admCoefficient> coeffs_kb_;
            
            //- Upper and lower limits on gas concentration
            scalarField upperLimits_;
            scalarField lowerLimits_;

            //- Concentration relaxation factors
            scalarList relaxationFactors_;

            //- Mass convergence criteria
            scalarField massConvergence_;
            
            //- Delta convergence criteria
            scalarField deltaConvergence_;
            
            //- Default near-zero mass error scale factor
            scalarField nearZeroMassScale_;

            //- Default near-zero delta error scale factor
            scalarField nearZeroDeltaScale_;

        // For the output dictionary

            //- Output dictionary (in timeName)
            IOdictionary outputDict_;
        
            //- Temperature
            dimensionedScalar T_;

            //- Total pressure
            dimensionedScalar pTotal_;
            
            //- Partial pressures
            PtrList<dimensionedScalar> pPartial_;

            //- Bulk liquid vapour pressure
            dimensionedScalar pLiq_;
            
            //- Also we output pAtm, but it already exists as a coefficient
            
            //- Integral of R * T * sGas dt from t_old to t_new, by species, an
            //  interim result used by boundary condition
            PtrList<dimensionedScalar> integralRTSgDt_;
            
            //- Gas flow rate
            dimensionedScalar qGas_;

            //- Gas concentrations - only up to date after
            //  calculatePressuresAndFlow()
            PtrList<dimensionedScalar> sGas_;
            
        // ODE related variables
        
            //- Temperature from previous iteration
            dimensionedScalar Told_;
        
            //- Atmospheric pressure from previous iteration
            dimensionedScalar pAtmOld_;
            
            //- Bulk liquid vapour pressure from previous iteration
            dimensionedScalar pLiqOld_;
        
            //- Field of Henry's constant for previous timestep
            PtrList<scalarField> khOld_;
            
            //- Bubble source term
            scalarField bubbleSource_;
            
            //- Gas concentrations from previous timestep
            scalarField sGasCoeffsOld_;

            //- Gas concentrations from previous iteration
            scalarField sGasCoeffs_;

            //- Gas concentrations after ode solution
            scalarField sGasCoeffsNew_;
            
        // For testing convergence
        
            //- Mass transfer from the liquid
            mutable scalarField mLiq_;
            
        // State control
        //  Member variables that need to be remembered during a saveState()

            boolList savedState_;
            PtrList<dimensionedScalar> T_a;
            PtrList<dimensionedScalar> Told_a;
            PtrList<PtrList<dimensionedScalar> > pPartial_a;
            PtrList<dimensionedScalar> pAtmOld_a;
            PtrList<dimensionedScalar> pLiq_a;
            PtrList<dimensionedScalar> pLiqOld_a;
            PtrList<PtrList<dimensionedScalar> > integralRTSgDt_a;
            PtrList<PtrList<scalarField> > khOld_a;
            PtrList<scalarField> sGasCoeffs_a;
            PtrList<scalarField> sGasCoeffsOld_a;
            PtrList<dictionary> outputDict_a;
            PtrList<scalarField> mLiq_a;
            labelList mLiqTimeIndex_a;

    // Private member functions

        //- Allocate memory by filling pointer lists
        void allocateMemory();
        
        //- Read outputDict, if present
        void readDict();

        //- Check for uniform() and/or constant() on coefficients that require
        //  it
        void checkUniformAndConstant() const;

        //- Check dimensions of all variables and coefficients
        void checkDimensions() const;

        //- Perform bubble mass transfer calculations on given species
        void bubbleTransfer(const label speciesIndex);
        
        //- Calculate the temperature in the gas volume - currently a simple
        //  average of the gas transfer surface temperature
        const dimensionedScalar calculateNewT() const;

        //- Calculate the bulk fluid vapour pressure in the gas volume - since
        //  pLiq depends on T, and T varies throughout the fluid volume, pLiq
        //  also varies throughout the fluid volume, even though that doesn't
        //  make sense.  This function returns the average of the pLiq over the
        //  cell volumes next to the gas transfer surface.
        const dimensionedScalar calculateNewPLiq() const;

        // Convergence test functions
        //  These calculate:
        //      percentError_i - convergenceCriterion_i
        //  for all species, and return the largest value.  Percent error is
        //  mass conservation (liquid versus gas transfer) or delta (change
        //  between iterations).  If either values are close to zero, error is
        //  scaled by nearZeroScale
        
            //- Test convergence based on change in variables between
            //  iterations
            scalar testDeltaConvergence() const;
        
            //- Test convergence based on mass conservation
            scalar testMassConvergence() const;
        
        //- Load sGasCoeffs_ back into the sGas_, calculate pPartial_ and qGas_
        void calculatePressuresAndFlow
        (
            const craftsGasFlowOde<matrixSize>& gasOde
        );

        //- Apply limits to the gas concentrations (after the ode solution)
        void applyLimits();

        //- Apply the relaxation factors while updating the solution
        //  sGasCoeffs_ with sGasCoeffsNew_
        void relax();
        
        //- Update the values in outputDict_
        void updateOutputDictionary();

public:

    //- Runtime type information
    TypeName("craftsGasModel");


    // Constructors

        //- Construct from components
        craftsGasModel
        (
            admCoefficientManager& admCoeffs,
            admVariableManager& admVars,
            craftsModel<matrixSize>& model,
            const word& craftsGasModelDictName = "craftsGasModelDict"
        );


    // Destructor

        virtual ~craftsGasModel()
        {}


    // Member functions

        // Access

            //- The input dictionary
            inline const IOdictionary& craftsGasModelDict() const;

            //- The surface patch name
            inline const word& surfacePatchName() const;

            //- The surface fvPatch
            inline const fvPatch& surfacePatch() const;

            //- The dissolved gas concentration variables
            inline const UPtrList<admVariable>& vars() const;

            //- The variables with the gas transfer boundary condition - these
            //  only differ from vars if the variable is a derived variable
            inline const UPtrList<admVariable>& bcVars() const;

            //- Output dictionary (in timeName)
            inline const IOdictionary outputDict() const;

        // Searching
        
            //- True if speciesName exists
            bool found(const word& speciesName);
            
            //- Returns the index number of the speciesName, -1 if not found
            label lookup(const word& speciesName) const;

            //- Returns the surfacePatch of the speciesName
            const fvPatchScalarField& surfacePatchField
            (
                const word& speciesName
            ) const;
            
            //- Returns the surfacePatch of the speciesIndex
            const fvPatchScalarField& surfacePatchField
            (
                const label speciesIndex
            ) const;

        // Functions executed during function hooks
        
            //- Create estimate for initial sGas values at new timestep.
            void initializeTimestep();

            //- Calculate mass transferred from liquid volume, and store
            //  varLiqStart
            void initializeImplicitLoop();
        
            //- Perform diffusion & mass transfer, and gas flow calculations
            //  on all species.  Returns error level:
            //      0 = success
            //     -1 = total failure, retry at shorter timestep
            //     -2 = convergence failure, solve coupled matrix again
            //          mass transfer error (from liquid versus to gas) too
            //          great - recalculate liquid
            //     -3 = stabilization failure, repeat implicitLoop
            //          previous iteration too different from current results
            //     -4 = function hook internal failure
            //          this error is never produced
            int implicitLoop();

            //- Executed when implicit loop converges - currently does nothing
            void finalizeImplicitLoop();

            //- Executed when timestep converges - saves local temporally
            //  interpolated variables
            void finalizeTimestep();

        // State control
        
            //- Save data to backup variables
            void saveState(const label slot = 0);
            
            //- Erase backup variables
            void clearState(const label slot = 0);
            
            //- Load backup variables into active variables
            void loadState(const label slot = 0);
            
            //- Return number of spots that are occupied
            //  A value of -1 means it doesn't track nStates
            label nStates() const;

            //- Returns true if this is a loadable save state slot
            bool validState(const label slot) const;
};


// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

} // End namespace Foam

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

#ifdef NoRepository
#   include "craftsGasModel.C"
#endif

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

#include "craftsGasModelI.H"

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

#endif

// ************************************************************************* //
