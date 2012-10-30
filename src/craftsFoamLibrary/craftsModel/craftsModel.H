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
    craftsModel

Description
    This is a PDE object that OpenFOAM can solve using ODESolver.  It is an
    object-oriented implementation of the Anaerobic Digester Model No. 1 (ADM1)
    as presented in the International Water Association's (IWA) Scientific and
    Technical Report No. 13 (STR13).  The implementation is flexible and can
    accommodate many modifications.

SourceFiles
    craftsModelI.H
    craftsModel.C
    craftsModelCreateIndices.C
    craftsModelFillTerms.C
    craftsModelImplicitAutoSolve.C
    craftsModelPointedToFunctions.C

Author
    David L. F. Gaden

\*---------------------------------------------------------------------------*/

#ifndef craftsModel_H
#define craftsModel_H

#include "admReactionReader.H"
#include "craftsFlow.H"
#include "craftsUdfs.H"
#include "Tuple2.H"
#include "blockLduSolvers.H"

#include "blockMatrixTools.H"

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

namespace Foam
{

/*---------------------------------------------------------------------------*\
                            Class craftsModel Declaration
\*---------------------------------------------------------------------------*/

template<int matrixSize>
class craftsModel
:
    public admReactionReader
{
    typedef VectorN<scalar, matrixSize> vectorType;
    typedef TensorN<scalar, matrixSize> tensorType;


public:

    // Enumeration for the type of step being taken (timestep doubling)
    enum stepType
    {
        COARSE,
        FINEA,
        FINEB
    };


private:
    // Private member data

        //- Copied flow field data
        const volVectorField& U_;
        const surfaceScalarField& phi_;


        // Function hooks
        
            //- Name of function hooks object
            word functionHooksName_;
            
            //- Function hooks object
            autoPtr<craftsUdfs<matrixSize> > hooks_;


        //- Flow model
        autoPtr<craftsFlow> flow_;

            //- This flag is used with a subStepping flow solver.  Allows it to
            //  load from the temporary save states to continue the solution.
            bool transitionToNextTimestep_;

        // Save spots, used for temporary storage of standard or implicit
        //  variable internal fields

            //- Save spots for standard variables
            PtrList<scalarField> internalStandardSaveSpot_;
            PtrList<scalarField> coarseStepStandardSaveSpot_;
            
            //- Save spots for implicit variables
            PtrList<scalarField> coarseStepImplicitSaveSpot_;


        // Output flags
        // These are yes/no settings that change what is output to the console
        // during a run.  They are contained in admSettingsDict.  Other classes
        // can add to this subdictionary.

            //- Summarizes results from one full step (default true)
            bool outputOverallSummary_;

            //- Summarizes reaction variable residuals (default true)
            bool outputReactionResidualSummary_;

            //- All reaction variable residuals (default false)
            bool outputReactionResidualDetails_;

            //- All reaction variable error scaling factors  (default false)
            bool outputReactionErrorScales_;

            //- Reaction solver performance (default false)
            bool outputReactionSolverPerformance_;

            //- Reaction timestep estimate (default true)
            bool outputReactionTimestepEstimate_;
            
            //- Implicit autosolve performance (default false)
            bool outputAutoSolvePerformance_;
            
            //- Implicit loop performance (default false)
            bool outputImplicitLoopSummary_;
            
            //- All implicit loop residuals (default false)
            bool outputImplicitLoopDetails_;


        // Control
        
            //- Maximum outer loop iterations (standard variables)
            label outerLoopMaxIterations_;

            //- Maximum inner loop iterations (implicit variables)
            label innerLoopMaxIterations_;

            //- Standard variable convergence criteria
            scalarField standardConvergence_;

            //- Implicit variable convergence criteria
            scalarField implicitConvergence_;


        // Convergence scaling factors

            //- Scaling factor for normalizing standard error, indexed by local
            //  variable index
            scalarField standardScale_;
            
            //- Scaling factor for normalizing implicit error, indexed by local
            //  variable index
            scalarField implicitScale_;

        //- For accounting for changes to standard variables that  are directly
        //  affected by user-defined functions.  Uses this method:
        //
        //      (in the middle of the step() function)
        //      udfDelta has value from previous iteration or timestep
        //          or set to zero for first iteration
        //      solve the PDE (S changes)
        //          note: udfDelta is subtracted as a source term to the PDE
        //              matrix in order to approximate the effect of the UDF on
        //              S
        //      S = S + udfDelta        --> addUdfDelta()
        //      preUdf = S              --> savePreUdf()
        //      run implicit routines (S changes again)
        //      oldUdfDelta = udfDelta
        //      udfDelta = preUdf - S   --> calculateUdfDelta()
        //      compare udfDelta with oldUdfDelta for convergence
        //      if not converged, repeat
        PtrList<scalarField> preUdf_;
        PtrList<scalarField> udfDelta_;
        PtrList<scalarField> oldUdfDelta_;

        //- Implicit variables are not actually calculated in the main coupled
        //  reaction matrix; rather they are calculated in user-defined
        //  functions.  However, in order for the standard variables to "feel"
        //  the effect of their dependence on implicit variables, the implicit
        //  variables are linearly approximated in the main coupled reaction
        //  matrix.  implicitDdt is used for this linear approximation source
        //  term, and is also used for convergence testing.
        PtrList<scalarField> implicitDdt_;
        PtrList<scalarField> oldImplicitDdt_;

        //- owner and neighbour list for lduAddressing of the matrix
        labelList owner_;
        labelList neighbour_;

            
        // Non-zero index lists - these store the non-zero components to
        //  reduce computational effort for the derivatives and Jacobian
        
            //- Coupled reaction matrix indices
            //   The coupled reaction matrix (Aij) relates standard variables
            //   and implicit variables, according to:
            //      [A]{S} = {B}
            //   Where [A] is the coupled reaction matrix, shown below:
            //
            //        standard     implicit
            //    ----------------------------
            //    |\s            |           | 
            //   s|  \d    nss   |           |
            //   t|    \i        |           |
            //   a|      \a      |    nsi    |
            //   n|        \g    |           |
            //   d|   oss    \   |           |
            //    |            \ |           |
            //    |--------------------------|
            //    |              |\i         |
            //   i|              |  \d   nii |
            //   m|     osi      |    \i     |
            //   p|              |      \a   |
            //    |              | oii    \g |
            //    |              |          \|
            //    ----------------------------
            //
            //   Reaction source terms are linearized using Newton-Raphson
            //   linearization.  Each component has an implicit term, which
            //   goes into the matrix, and an explicit term, which goes into
            //   the source term.  They are given by:
            //      Aij = -Rij = -dRi/dSj  (implicit term)
            //      Bi = -sum(dRi/dSk * Sk) (explicit term)
            //   where
            //      Rij is described below.
            //   Reaction source term equations in the upper and lower triangle
            //   of the matrix are evaluated faster by applying an indexing to
            //   each equation term, and looping through the indices.  The
            //   equation is given by:
            //      Rij = d(R_i)/d(y_j) = nu_ik * d(rate_k)/d(y_j)
            //                            + d(nu_ik)/d(y_j) * rate_k
            //   where:
            //      i = row
            //      j = column
            //      k = reaction index (dummy index that is summed over)
            //   The first term (termI) is on the first line, the second term
            //   (termII) is on the second line.  Each index pair (i,j) have a
            //   corresponding upper triangle, lower triangle and diagonal
            //   value.  These are referred to as owners (lower) and neighbours
            //   (upper).  An empty k list indicates a term that is zero.  The
            //   indexing is given by:
            //      i and j are given by:
            //          for diags, i = j = m;
            //          for lower: i is owner_[m], j is neighbour_[m];
            //          for upper: i is neighbour_[m], j is owner_[m];
            //      k is a list given by:
            //          ossTermIk_[m], ossTermIIk_[m],
            //          nssTermIk_[m], nssTermIIk_[m],
            //          sdiagTermIk_[m], or sdiagTermIIk_[m]
            //      where m is the variable index.
            //
            //   Implicit variables are pulled out and solved seperately, but
            //   a first order approximation of them is left in the matrix.
            //   Therefore, idiag and the associated source term hold terms for
            //   this linear approximation, and osi, oii and nii are zero.

            PtrList<UPtrList<const admReaction> > ossTermIk_;
            PtrList<UPtrList<const admReaction> > ossTermIIk_;
            PtrList<UPtrList<const admReaction> > nssTermIk_;
            PtrList<UPtrList<const admReaction> > nssTermIIk_;
            PtrList<UPtrList<const admReaction> > nsiTermIk_;
            PtrList<UPtrList<const admReaction> > nsiTermIIk_;
            PtrList<UPtrList<const admReaction> > sDiagTermIk_;
            PtrList<UPtrList<const admReaction> > sDiagTermIIk_;
            
            //- The built-in implicitAutoSolve function solves a single
            //   implicit variable at a time, holding all others constant.
            PtrList<UPtrList<const admReaction> > dRdSAutoSolveTermIk_;
            PtrList<UPtrList<const admReaction> > dRdSAutoSolveTermIIk_;

            //- Mesh owner and neighbour indices - these are shortcuts used
            //   to speed up transferring terms from the transport matrix to
            //   the reaction matrix.  transportOwners_[cellIndex] and
            //   transportNeighbours_[cellIndex] are a list of internal owner
            //   and neighbour faces for cellIndex.  The associated cells are
            //   given by transportOwners/NeighboursCells_[cellIndex]
            labelListList transportOwners_;
            labelListList transportNeighbours_;
            labelListList transportOwnersCells_;
            labelListList transportNeighboursCells_;

        //- The latest residual calculation is stored in this list, indexed by
        //  variable number.  This is for standard variables.
        mutable scalarField standardResidual_;

        //- The latest residual calculation is stored in this list, indexed by
        //  variable number.  This is for implicit variables.
        mutable scalarField implicitResidual_;

        //- The latest residual calculation is stored in this list, indexed by
        //  variable number.  This is for implicitDdt.
        mutable scalarField implicitDdtResidual_;

        //- The latest residual calculation is stored in this list, indexed by
        //  variable number.  This is for udfDelta.
        mutable scalarList udfDeltaResidual_;
        

        // Adaptive timestepping

            //- True if using errorBased timestep adjustment
            Switch atsUseAts_;

            //- Adaptive timestep convergence factor - the convergence criteria
            //  are multipiled by this factor when checking for an acceptable
            //  time step size
            scalar atsConvergenceFactor_;

            //- Adaptive timestep overclock factor - the final timestep that
            //  the adaptive timestepping algorithm chooses is multiplied by
            //  this fudge factor.  Allows you to pressure the timestep
            //  upwards, while leaving the convergence criteria alone.
            scalar atsOverclockFactor_;

            //- Maximum allowed increase in the timestep per iteration
            scalar atsMaxIncreaseFactor_;

            //- Maximum allowed reduction in the timestep per iteration
            scalar atsMaxReductionFactor_;
            
            //- Minimum reduction in the timestep when redoing a timestep
            //  (e.g. prevents retrying with a timestep 0.99 x original)
            scalar atsMinReductionFactor_;


            // Performance feedback
            //  An optional component of adaptive timestepping.  This one uses
            //  performance stats to compare how long a single coarse step took
            //  versus two fine steps.  If the fine steps were faster,
            //  performance feedback can force the timestep to smaller values
            //  in the interest of speeding up the overall process.
            
                //- True if using performance stats for timestep adjustment
                Switch atsPfUsePf_;

                //- Performance metric 0=unknown, 1=iterations, 2=cpuTime
                label atsPfMeasure_;
            
                //- Performance curve bias - negative favours smaller
                //  timesteps, positive favours larger timesteps, zero is a
                //  linear blend
                scalar atsPfBias_;
            

            // Performance data
            //  Information is stored on the last time step() was run.  The
            //  values are used for performance based adaptive timestepping.
            //  For instance, if it was faster to do two mini-steps than one
            //  large step, reducing the timestep may be faster.  This check
            //  can be in addition to timestep doubling accuracy checks.
            
                //- Total number of iterations required in the last step
                //  Custom function hooks can modify this number if desired
                label lastStepIterations_;
                label coarseStepIterations_;
                label fineStepIterations_;

                //- Total cpu time spent in the last step
                scalar lastStepCpuTime_;
                scalar coarseStepCpuTime_;
                scalar fineStepCpuTime_;


    // Private member functions

        //- Disallow default bitwise copy construct and assignment
        craftsModel(const craftsModel&);
        void operator=(const craftsModel&);

        //- Read convergence criteria for each standard and implicit variables
        void readConvergenceCriteria();


        // Create indices and sublists

            //- Build / rebuild all non-zero index lists
            void initialize();

            //- Create lduMesh for lduAddressing of matrix
            //  *** Located in craftsModelCreateIndices.C ***
            void createLduMesh();

            //- Create transportOwner_ and transportNeighbour indices
            //  *** Located in craftsModelCreateIndices.C ***
            void createTransportFaceIndices();
            
            //- Create autoSolve indices
            //  *** Located in craftsModelCreateIndices.C ***
            void createAutoSolveIndices();


        // Fill terms
        
            //- Put transport-related terms into the coupled reaction matrix
            //  *** Located in craftsModelFillTerms.C ***
            void fillTransportTerms
            (
                BlockLduMatrix<vectorType>& blockM,
                Field<vectorType>& blockX,
                Field<vectorType>& blockB
            );

            //- Put diagonals and source terms into the coupled reaction matrix
            //  for an approximation of the implicit variables.
            //  *** Located in craftsModelFillTerms.C ***
            void fillImplicitTerms
            (
                BlockLduMatrix<vectorType>& blockM,
                Field<vectorType>& blockX,
                Field<vectorType>& blockB
            );

            //- Put Newton-Raphson diagonal reaction terms into the matrix
            //  *** Located in craftsModelFillTerms.C ***
            void fillSDiagReactionTerms
            (
                BlockLduMatrix<vectorType>& blockM,
                Field<vectorType>& blockB
            );
            
            //- Put off-diagonal reaction terms into the matrix
            //  *** Located in craftsModelFillTerms.C ***
            void fillOssNssNsiReactionTerms
            (
                BlockLduMatrix<vectorType>& blockM,
                Field<vectorType>& blockB
            );
            
            //- Put main reaction source term into matrix source term
            //  *** Located in craftsModelFillTerms.C ***
            void fillSourceReactionTerms
            (
                BlockLduMatrix<vectorType>& blockM,
                Field<vectorType>& blockB
            );


        // UdfDelta functions
        //  Support functions for linearly approximating the effect of user-
        //  defined functions that change standard variables.
  
            //- Add udfDelta fields to standard variables that have
            //  changedByUDF flag set
            void addUdfDelta();
            
            //- Store the value of a standard variable (with changedByUDF set)
            //  before running user-defined functions
            void savePreUdf();
            
            //- Calculate the change in a variable after passing through the
            //  user-defined functions
            void calculateUdfDelta();


        // ImplicitDdt functions
        //  implicitDdt is used both to linearly approximate the effect of
        //  a UDF function on the implicit variables for the coupled reaction
        //  matrix, and to measure convergence of the implicit routines.
        
            //- Calculate the implicitDdt_ scalarField
            void getImplicitDdt();

public:

    TypeName("craftsModel");


    // Constructors

        //- Construct, reading from dictionary
        craftsModel
        (
            admTime& runTime,
            fvMesh& mesh,
            surfaceScalarField& phi,
            volVectorField& U,
            volScalarField& p,
            const word& variableDictName = "admVariableDict",
            const word& coefficientDictName = "admCoefficientDict",
            const word& inhibitionDictName = "admInhibitionDict",
            const word& reactionDictName = "admReactionDict",
            const word& settingsDictName = "admSettingsDict",
            const word& hooksDictName = "admHooksDict"
        );


    // Destructor
    
        virtual ~craftsModel();


    // Member Functions
    
        // Access

            //- Velocity field
            inline const volVectorField& U() const;
            inline const surfaceScalarField& phi() const;

            //- Name of function hooks object
            inline const word& functionHooksName() const;
            
            //- Function hooks object
            inline autoPtr<craftsUdfs<matrixSize> >& hooks();
            inline const autoPtr<craftsUdfs<matrixSize> >& hooks() const;

            //- The flow model
            inline const craftsFlow& flow() const;
            inline craftsFlow& flow();


            // Output flags

                //- Output flag for results summary from one full step
                inline const bool& outputOverallSummary() const;

                //- Output flag for reaction variable residual summary
                inline const bool& outputReactionResidualSummary() const;

                //- Output flag for all reaction variable residuals
                inline const bool& outputReactionResidualDetails() const;

                //- Output flag for reaction solver performance
                inline const bool& outputReactionSolverPerformance() const;

                //- Output flag for reaction timestep estimate
                inline const bool& outputReactionTimestepEstimate() const;

                //- Output flag for implicit autosolve performance
                inline const bool& outputAutoSolvePerformance() const;

                //- Output flag for implicit loop performance
                inline const bool& outputImplicitLoopSummary() const;

                //- Output flag for all implicit loop residuals
                inline const bool& outputImplicitLoopDetails() const;

            //- Maximum iterations
            inline const label& outerLoopMaxIterations() const;
            inline const label& innerLoopMaxIterations() const;
            
            //- Convergence criteria
            inline const scalarField& standardConvergence() const;
            inline const scalarField& implicitConvergence() const;

            //- Convergence scaling factors
            inline const scalarField& standardScale() const;
            inline const scalarField& implicitScale() const;

            //- Latest residual calculations
            inline const scalarField& standardResidual() const;
            inline const scalarField& implicitResidual() const;
            inline const scalarField& implicitDdtResidual() const;
            inline const scalarField& udfDeltaResidual() const;

            //- Adaptive timestepping
            //  See in member variables above for descriptions
            inline const Switch& atsUseAts() const;
            inline const scalar& atsConvergenceFactor() const;
            inline const scalar& atsOverclockFactor() const;
            inline const scalar& atsMaxIncreaseFactor() const;
            inline const scalar& atsMaxReductionFactor() const;
            inline const scalar& atsMinReductionFactor() const;
            inline const Switch& atsPfUsePf() const;
            inline const label& atsPfMeasure() const;
            inline const scalar& atsPfBias() const;
            inline const label& lastStepIterations() const;
            inline label& lastStepIterations();
            inline const scalar& lastStepCpuTime() const;


        // Calculate the scaling factors for normalizing the error calculations

            //- Calculate the scaling factors for all variables
            void calculateScales();

            //- Standard variables scale factors
            void calculateStandardScale();

            //- Implicit variables scale factors
            void calculateImplicitScale();


        // Convergence testing
        
            //- Test the convergence for all model variables.  For each
            //  variable, calculates a volume-weighted root mean square error,
            //  scaled by their respective scaling factors, then subtracts the
            //  convergence criterion.  Returns the highest value encountered.
            //  A negative return indicates all variables have converged.
            scalar testConvergence() const;

            //- Get the residual for standard variables.  Compares the current
            //  variables' internal fields with those stored in saveSpot, and
            //  scales the residual using the standardScale_ list.  Stores the
            //  result for each variable in standardResidual_.  Returns a tuple
            //  with first giving the global index of the variable with the
            //  greatest rmsError, and the second containing a scalar with the
            //  maximum (scaledResidual - convergence * factor).  A negative
            //  indicates convergence criteria have been met.  factor is
            //  included for adaptive timestepping, which may multiply the
            //  convergence criteria by an adjustment factor.
            //
            //  reportResults prints residuals to console as:
            //      "title residuals test:"
            Tuple2<label, scalar> testStandardConvergence
            (
                const PtrList<scalarField>& saveSpot,
                const scalar& factor = 1.0,
                bool reportResults = false,
                word title = "Unknown"
            ) const;

            //- Get the residual for implicit variables.  Compares the current
            //  variables' internal fields with those stored in saveSpot, and
            //  scales the residual using the implicitScale_ list.  Stores the
            //  result for each variable in standardResidual_.  Returns a tuple
            //  with first giving the global index of the variable with the
            //  greatest rmsError, and the second containing a scalar with the
            //  indicates convergence criteria have been met.  factor is
            //  included for adaptive timestepping, which may multiply the
            //  convergence criteria by an adjustment factor.
            //
            //  reportResults prints residuals to console as:
            //      "title residuals test:"
            Tuple2<label, scalar> testImplicitConvergence
            (
                const PtrList<scalarField>& saveSpot,
                const scalar& factor = 1.0,
                bool reportResults = false,
                word title = "Unknown"
            ) const;
            
            //- Get the residual for implicit variables.  Compares the current
            //  implicitDdt fields with those stored in saveSpot, and
            //  scales the residual using the implicitScale_ list.  Stores the
            //  result for each variable in implicitResidual_.  Returns a tuple
            //  with first giving the global index of the variable with the
            //  greatest rmsError, and the second containing a scalar with the
            //  maximum (scaledResidual - convergence).  A negative indicates
            //  convergence criteria have been met.
            Tuple2<label, scalar> testImplicitDdtConvergence() const;

            //- Calculate the residual for the udfDelta.  Compares the current
            //  udfDelta with oldUdfDelta for each changedByUdf variable.
            //  Saves these to udfDeltaResidual_.  Returns a tuple with first
            //  giving the global index of the variable with the greatest
            //  rmsError, and the second containing a scalar with the maximum
            //  (scaledUdfDeltaResidual - convergence) for each variable.  A
            //  negative indicates convergence criteria have been met.
            Tuple2<label, scalar> testUdfDeltaConvergence() const;


        // Adaptive timestepping
        
            //- Returns the "best overall time" based on a comparison between
            //  coarse save spots and current values for standard and implicit
            //  variables.
            scalar calculateNextBestDeltaT() const;

            //- Use this to inform the model that the timestep was acceptable,
            //  and we are now moving to the next one.  When a substepping flow
            //  model is used, this tells it to load the temporary save states
            //  at the start points of the timestep.
            void transitionToNextTimestep();


        // Solution - Flow and reactions combined

            //- Performs step(), records coarse performance, stores coarse
            //  solution.  Increments runTime.  Returns error level:
            //      0 = success
            //     -1 = total failure, retry at shorter timestep
            label coarseStep();
            
            //- Performs two step(), records fine performance.  Increments
            //  runTime.  Returns error level:
            //      0 = success
            //     -1 = total failure, retry at shorter timestep
            label doubleFineStep();


        // Solution - Flow only

            //- This is a dispatch for flow_->step() that allows it to perform
            //  sub-stepping.  Increments runTime by a normal step size.
            //  Passes on return value from flow_->step().  Behaviour changes
            //  depending on the type of step, st, e.g.:
            //
            //    st     | runTime increment           | stores coarse solution
            //    -------+-----------------------------+-----------------------
            //    COARSE | runTime_.plusPlusNoOutput() | YES
            //    FINEA  | runTime_.plusPlusNoOutput() | NO
            //    FINEB  | runTime_++;                 | NO
            label stepFlowDispatch(stepType st);


        // Solution - Reactions only

            //- This is a dispatch for stepReactions().  Passes along return
            //  value from stepReactions().  This function changes behaviour
            //  based on the type of step, st:
            //    COARSE - record coarse performance
            //             reset fine performance
            //             save variables to coarse save spot
            //    FINEA, FINEB - record fine performance
            label stepReactionsDispatch(stepType st);

            //- Solve the PDAE for a single timestep.  Returns error level:
            //      0 = success
            //     -1 = total failure, retry at shorter timestep
            label stepReactions();
            
            //- Performs the implicit routines sequentially on each cell in the
            //  mesh, and repeats until they stabilize.  Returns error level:
            //      0 = success
            //     -1 = total failure, retry at shorter timestep
            //     -2 = convergence failure, solve coupled matrix again
            label runImplicitRoutines();


            // implicitAutoSolve
            //  The implicitAutoSolve is a built-in implicit solver for stiff
            //  variables.  It pulls these variables out of the PDE, and solves
            //  them separately, holding the other variables constant.  It is
            //  only valid for implicit variables marked with autoSolve.
            //  All implicitAutoSolve functions return the error level:
            //      0 = success
            //     -1 = total failure, retry at shorter timestep
            //     -2 = convergence failure, solve coupled matrix again

                //- autoSolve a given variable at a given cellIndex
                //  This function performs wasteful calculations.  Use
                //  implicitAutoSolveField whenever possible.
                //  *** Located in craftsModelImplicitAutoSolve.C ***
                label implicitAutoSolve
                (
                    admImplicitVariable& impVar,
                    label cellIndex
                );

                //- autoSolve a given variable over the entire mesh
                //  *** Located in craftsModelImplicitAutoSolve.C ***
                label implicitAutoSolveField
                (
                    admImplicitVariable& impVar
                );

                //- autoSolve all variables marked with autoSolve, one at a
                //  time, at a given cell index.  This function performs
                //  wasteful calculations.  Use implicitAutoSolveAllField
                //  whenever possible.
                //  *** Located in craftsModelImplicitAutoSolve.C ***
                label implicitAutoSolveAll
                (
                    label cellIndex
                );

                //- autoSolve all variables marked with autoSolve, one at a
                //  time, over the entire mesh.
                //  *** Located in craftsModelImplicitAutoSolve.C ***
                label implicitAutoSolveAllField();


        // State control
        //  Since OpenFOAM does not easily support undoing a timestep, or
        //  stepping back several timesteps, the save state framework must be
        //  used.

            //- Shortcut for saving flow and time data only
            void saveFlowState(const label slot);

            //- Shortcut for loading flow and time data only
            void loadFlowState(const label slot);
        
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


        // Debugging easy output
            
            //- Report scaling values
            void reportScales(Ostream& os);
};


// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

} // End namespace Foam

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

#include "craftsModelI.H"
#ifdef NoRepository
#   include "craftsModel.C"
#endif

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

#endif

// ************************************************************************* //
