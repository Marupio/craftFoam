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

//#include "IOcraftsModelReferencer.H"
#include "craftsModel.H"
#include "admFlow.H"
#include "IOstreams.H"
#include "token.H"

// * * * * * * * * * * * * * Private Member Functions  * * * * * * * * * * * //

template <int matrixSize>
void Foam::craftsModel<matrixSize>::readConvergenceCriteria()
{
    forAll(admVars_.standard(), varIndex)
    {
        const admStandardVariable& var(admVars_.standard(varIndex));
        standardConvergence_[varIndex] = var.convergence();
    }
    forAll(admVars_.implicit(), varIndex)
    {
        const admImplicitVariable& var(admVars_.implicit(varIndex));
        implicitConvergence_[varIndex] = var.convergence();
    }
}


template <int matrixSize>
void Foam::craftsModel<matrixSize>::initialize()
{
    // Create lduMesh used for lduAdressing of matrix 
    createLduMesh();
    
    // Create owner & neighbour face lists for mesh of entire geometry
    createTransportFaceIndices();

    // Create autoSolve nonZero indices
    createAutoSolveIndices();

    // Initialize implicitDdt, oldImplicitDdt, and dtEst
    forAll(admVars_.implicit(), varIndex)
    {
        implicitDdt_.set
        (
            varIndex,
            new scalarField(mesh_.nCells(), 0.0)
        );
        oldImplicitDdt_.set
        (
            varIndex,
            new scalarField(mesh_.nCells(), 0.0)
        );
    }
    forAll(udfDelta_, deltaIndex)
    {
        preUdf_.set
        (
            deltaIndex,
            new scalarField(mesh_.nCells(), 0.0)
        );
        udfDelta_.set
        (
            deltaIndex,
            new scalarField(mesh_.nCells(), 0.0)
        );
        oldUdfDelta_.set
        (
            deltaIndex,
            new scalarField(mesh_.nCells(), 0.0)
        );
    }
    
    // Output headings for output reported to the console as arrays
    if (outputReactionResidualDetails_)
    {
        Info << "reactionResidualDetails: Headings for numbers in the output "
            << "body are:" << endl;
        Info << "reactionResidualDetails: Standard res = "
            << admVars_.tocStandard() << endl;
        Info << "reactionResidualDetails: Implicit res = "
            << admVars_.tocImplicit() << endl;
    }
    
    if (outputReactionErrorScales_)
    {
        Info << "reactionErrorScales: Headings for numbers in the output body "
            << "are:" << endl;
        Info << "reactionErrorScales: Standard scales = "
            << admVars_.tocStandard() << endl;
        Info << "reactionErrorScales: Implicit scales = "
            << admVars_.tocImplicit() << endl;
    }
    
    if (flow_->subStepping())
    {
        // Create initial save states
        saveState(0); // reaction start of current timestep
        saveState(1); // reaction temporary
        saveFlowState(2); // flow start of current timestep
        saveFlowState(3); // flow temporary
    }
}


template <int matrixSize>
void Foam::craftsModel<matrixSize>::addUdfDelta()
{
    forAll(admVars_.changedByUdf(), i)
    {
        label varIndex(admVars_.changedByUdf(i).localIndex());
        admVars_.standard(varIndex)().internalField() += udfDelta_[i];
    }
}


template <int matrixSize>
void Foam::craftsModel<matrixSize>::savePreUdf()
{
    forAll(admVars_.changedByUdf(), i)
    {
        label varIndex(admVars_.changedByUdf(i).localIndex());
        preUdf_[i] = admVars_.standard(varIndex)().internalField();
    }
}


template <int matrixSize>
void Foam::craftsModel<matrixSize>::calculateUdfDelta()
{
    forAll(admVars_.changedByUdf(), i)
    {
        label varIndex(admVars_.changedByUdf(i).localIndex());
        udfDelta_[i] =
            preUdf_[i] - admVars_.standard(varIndex)().internalField();
    }
}


template <int matrixSize>
void Foam::craftsModel<matrixSize>::getImplicitDdt()
{
    scalar rDtime(1 / runTime_.deltaT().value());
    
    // Calculate ddt term
    forAll(admVars_.implicit(), varIndex)
    {
        implicitDdt_[varIndex] =
        (
            admVars_.implicit
            (
                varIndex
            )().internalField() -
            admVars_.implicit(varIndex)().oldTime().internalField()
        ) * rDtime;
    }
}

// * * * * * * * * * * * * * * * * Constructors  * * * * * * * * * * * * * * //

template<int matrixSize>
Foam::craftsModel<matrixSize>::craftsModel
(
    admTime& runTime,
    fvMesh& mesh,
    surfaceScalarField& phi,
    volVectorField& U,
    volScalarField& p,
    const word& variableDictName,
    const word& coefficientDictName,
    const word& inhibitionDictName,
    const word& reactionDictName,
    const word& settingsDictName,
    const word& hooksDictName
)
:
    admReactionReader
    (
        runTime,
        mesh,
        variableDictName,
        coefficientDictName,
        inhibitionDictName,
        reactionDictName,
        settingsDictName,
        hooksDictName
    ),
    
    U_(U),
    phi_(phi),

    functionHooksName_(admSettingsDict_.lookup("functionHooks")),
    hooks_
    (
        admFcHooks<matrixSize>::New
        (
            functionHooksName_,
            *this,
            hooksDictName
        )
    ),

    flow_
    (
        admFlow::New
        (
            U,
            phi,
            p,
            *this,
            admSettingsDict_.subDict("flowModel").lookup("type")
        )
    ),
    
    transitionToNextTimestep_(false),

    internalStandardSaveSpot_(0),
    coarseStepStandardSaveSpot_(0),
    coarseStepImplicitSaveSpot_(0),

    outputReactionResidualSummary_
    (
        outputFlagsDict_.found("reactionResidualSummary")
      ? bool(Switch(outputFlagsDict_.lookup("reactionResidualSummary")))
      : true
    ),
    outputReactionResidualDetails_
    (
        outputFlagsDict_.found("reactionResidualDetails")
      ? bool(Switch(outputFlagsDict_.lookup("reactionResidualDetails")))
      : false
    ),
    outputReactionErrorScales_
    (
        outputFlagsDict_.found("reactionErrorScales")
      ? bool(Switch(outputFlagsDict_.lookup("reactionErrorScales")))
      : false
    ),
    outputReactionSolverPerformance_
    (
        outputFlagsDict_.found("reactionSolverPerformance")
      ? bool(Switch(outputFlagsDict_.lookup("reactionSolverPerformance")))
      : false
    ),
    outputReactionTimestepEstimate_
    (
        outputFlagsDict_.found("reactionTimestepEstimate")
      ? bool(Switch(outputFlagsDict_.lookup("reactionTimestepEstimate")))
      : true
    ),
    outputAutoSolvePerformance_
    (
        outputFlagsDict_.found("implicitAutoSolvePerformance")
      ? bool(Switch(outputFlagsDict_.lookup("implicitAutoSolvePerformance")))
      : false
    ),
    outputImplicitLoopSummary_
    (
        outputFlagsDict_.found("implicitLoopSummary")
      ? bool(Switch(outputFlagsDict_.lookup("implicitLoopSummary")))
      : false
    ),
    outputImplicitLoopDetails_
    (
        outputFlagsDict_.found("implicitLoopDetails")
      ? bool(Switch(outputFlagsDict_.lookup("implicitLoopDetails")))
      : false
    ),

    outerLoopMaxIterations_
    (
        readLabel(admSettingsDict_.lookup("outerLoopMaxIterations"))
    ),

    innerLoopMaxIterations_
    (
        readLabel(admSettingsDict_.lookup("innerLoopMaxIterations"))
    ),

    standardConvergence_
    (
        admVars_.nStandard(),
        scalar(-1.0)
    ),

    implicitConvergence_
    (
        admVars_.nImplicit(),
        scalar(-1.0)
    ),

    standardScale_(admVars_.nStandard()),
    implicitScale_(admVars_.nImplicit()),

    preUdf_(admVars_.changedByUdf().size()),
    udfDelta_(admVars_.changedByUdf().size()),
    oldUdfDelta_(admVars_.changedByUdf().size()),

    implicitDdt_(admVars_.nImplicit()),
    oldImplicitDdt_(admVars_.nImplicit()),
    
    owner_(),
    neighbour_(),

    ossTermIk_(),
    ossTermIIk_(),
    nssTermIk_(),
    nssTermIIk_(),
    nsiTermIk_(),
    nsiTermIIk_(),
    sDiagTermIk_(admVars_.nStandard()),
    sDiagTermIIk_(admVars_.nStandard()),

    dRdSAutoSolveTermIk_(admVars_.implicitAutoSolve().size()),
    dRdSAutoSolveTermIIk_(admVars_.implicitAutoSolve().size()),

    transportOwners_(mesh_.nCells()),
    transportNeighbours_(mesh_.nCells()),
    transportOwnersCells_(mesh_.nCells()),
    transportNeighboursCells_(mesh_.nCells()),

    standardResidual_(admVars_.nStandard()),
    implicitResidual_(admVars_.nImplicit()),
    implicitDdtResidual_(admVars_.nImplicit()),
    udfDeltaResidual_(admVars_.changedByUdf().size()),

    atsUseAts_
    (
        admSettingsDict_.subDict("adaptiveTimeStepping")
            .lookup("useAdaptiveTimeStepping")
    ),
    atsConvergenceFactor_
    (
        atsUseAts_
      ? readScalar
        (
            admSettingsDict_.subDict("adaptiveTimeStepping")
                .lookup("convergenceFactor")
        )
      : -1.0
    ),
    atsOverclockFactor_
    (
        admSettingsDict_.subDict("adaptiveTimeStepping")
                .found("overclockFactor")
      ? readScalar(admSettingsDict_.subDict("adaptiveTimeStepping")
                .lookup("overclockFactor"))
      : 1.0
    ),
    atsMaxIncreaseFactor_
    (
        readScalar
        (
            admSettingsDict_.subDict("adaptiveTimeStepping")
                .lookup("maxIncreaseFactor")
        )
    ),
    atsMaxReductionFactor_
    (
      readScalar
        (
            admSettingsDict_.subDict("adaptiveTimeStepping")
                .lookup("maxReductionFactor")
        )
    ),
    atsMinReductionFactor_
    (
        readScalar
        (
            admSettingsDict_.subDict("adaptiveTimeStepping")
                .lookup("minReductionFactor")
        )
    ),
    atsPfUsePf_
    (
        admSettingsDict_.subDict("adaptiveTimeStepping")
            .subDict("performanceFeedback")
            .lookup("usePerformanceFeedback")
    ),
    atsPfMeasure_(0),
    atsPfBias_
    (
        atsPfUsePf_
      ? readScalar
        (
            admSettingsDict_.subDict("adaptiveTimeStepping")
                .subDict("performanceFeedback").lookup("bias")
        )
      : 0.0
    )
{
    if (atsPfUsePf_)
    {
        word atsPfMeasureWord
        (
            admSettingsDict_.subDict("adaptiveTimeStepping")
                .subDict("performanceFeedback").lookup("measure")
        );
        if (atsPfMeasureWord == "iterations")
        {
            atsPfMeasure_ = 1;
        }
        else if (atsPfMeasureWord == "cpuTime")
        {
            atsPfMeasure_ = 2;
        }
        else
        {
            FatalIOErrorIn
            (
                "craftsModel::craftsModel",
                admSettingsDict_.subDict("adaptiveTimeStepping")
                    .subDict("performanceFeedback")
            )
                << "Expected either 'iterations' or 'cpuTime', read '"
                << atsPfMeasureWord << "'"
                << exit(FatalIOError);
        }
    }
    readConvergenceCriteria();
    initialize();
}

// * * * * * * * * * * * * * * * * Destructor  * * * * * * * * * * * * * * * //

template <int matrixSize>
Foam::craftsModel<matrixSize>::~craftsModel()
{}

// * * * * * * * * * * * * * * Member Functions  * * * * * * * * * * * * * * //

template <int matrixSize>
void Foam::craftsModel<matrixSize>::calculateScales()
{
    calculateStandardScale();
    calculateImplicitScale();
    flow_->calculateScales();
    if (outputReactionErrorScales_)
    {
        Info << "reactionErrorScales: Standard scales = " << standardScale_
            << endl;
        Info << "reactionErrorScales: Implicit scales = " << implicitScale_
            << endl;
    }
}


template <int matrixSize>
void Foam::craftsModel<matrixSize>::calculateStandardScale()
{
    // Root mean square of each variable field
    forAll(admVars_.standard(), varIndex)
    {
        const admStandardVariable& var(admVars_.standard(varIndex));
        standardScale_[varIndex] = max
        (
            calculateRms
            (
                var().internalField(),
                mesh_.V().field()
            ),
            var.minErrorScale()
        );
    }
}


template <int matrixSize>
void Foam::craftsModel<matrixSize>::calculateImplicitScale()
{
    // Root mean square of each variable field
    forAll(admVars_.implicit(), varIndex)
    {
        const admImplicitVariable& var(admVars_.implicit(varIndex));
        implicitScale_[varIndex] = max
        (
            calculateRms
            (
                var().internalField(),
                mesh_.V().field()
            ),
            var.minErrorScale()
        );
    }
}


template <int matrixSize>
scalar Foam::craftsModel<matrixSize>::testConvergence() const
{
    Tuple2<label, scalar> standardResult
    (
        testStandardConvergence
        (
            coarseStepStandardSaveSpot_,
            atsConvergenceFactor_,
            outputReactionResidualDetails_,
            "reactionResidualDetails"
        )
    );

    Tuple2<label, scalar> implicitResult
    (
        testImplicitConvergence
        (
            coarseStepImplicitSaveSpot_,
            atsConvergenceFactor_,
            outputReactionResidualDetails_,
            "reactionResidualDetails"
        )
    );
    
    if (standardResult.second() < implicitResult.second())
    {
        standardResult = implicitResult;
    }

    scalar flowResult(flow_->testConvergence());
    
    // Output results to console
    if (outputReactionResidualSummary_)
    {
        const admVariable& worstVariable
        (
            admVars_.all(standardResult.first())
        );

        Info << "reactionResidualSummary: ";
        if (standardResult.second() > 0)
        {
            Info << "(fail), " << worstVariable.name() << " is "
            << standardResult.second();
            if (flowResult > standardResult.second())
            {
                Info << "; flow model is worse";
            }
            Info << "." << endl;
        }
        else
        {
            Info << "(pass)";
            if (flowResult > 0)
            {
                Info << ", but flow model fails";
            }
            Info << endl;
        }
    }

    return max(flowResult, standardResult.second());
}


template <int matrixSize>
Tuple2<label, scalar> Foam::craftsModel<matrixSize>::testStandardConvergence
(
    const PtrList<scalarField>& saveSpot,
    const scalar& factor,
    bool reportResults,
    word title
) const
{
    scalarList residuals(saveSpot.size());
    Tuple2<label, scalar> returnMe(-1, -VGREAT);
    forAll(standardResidual_, varIndex)
    {
        standardResidual_[varIndex] = calculateRmsError
        (
            admVars_.standard(varIndex)().internalField(),
            saveSpot[varIndex],
            mesh_.V().field(),
            standardScale_[varIndex]
        );
        residuals[varIndex] = standardResidual_[varIndex]
              - standardConvergence_[varIndex] * factor;
        if (returnMe.second() < residuals[varIndex])
        {
            returnMe.second() = residuals[varIndex];
            returnMe.first() = admVars_.standard(varIndex).globalIndex();
        }
    }
    if (reportResults)
    {
        Info << title << ": Standard res = " << residuals << endl;
    }
    return returnMe;
}


template <int matrixSize>
Tuple2<label, scalar> Foam::craftsModel<matrixSize>::testImplicitConvergence
(
    const PtrList<scalarField>& saveSpot,
    const scalar& factor,
    bool reportResults,
    word title
) const
{
    scalarList residuals(saveSpot.size());
    Tuple2<label, scalar> returnMe(-1, -VGREAT);
    forAll(implicitResidual_, varIndex)
    {
        implicitResidual_[varIndex] = calculateRmsError
        (
            admVars_.implicit(varIndex)().internalField(),
            saveSpot[varIndex],
            mesh_.V().field(),
            implicitScale_[varIndex]
        );
        residuals[varIndex] = implicitResidual_[varIndex]
            - implicitConvergence_[varIndex] * factor;
        if (returnMe.second() < residuals[varIndex])
        {
            returnMe.first() = admVars_.implicit(varIndex).globalIndex();
            returnMe.second() = residuals[varIndex];
        }
    }
    if (reportResults)
    {
        Info << title << ": Implicit res = " << residuals << endl;
    }
    return returnMe;
}


template <int matrixSize>
Tuple2<label, scalar> Foam::craftsModel<matrixSize>::
    testImplicitDdtConvergence() const
{
    scalarList residuals(implicitDdt_.size());
    Tuple2<label, scalar> returnMe(-1, -VGREAT);
    forAll(implicitDdt_, varIndex)
    {
        implicitDdtResidual_[varIndex] = calculateRmsError
        (
            implicitDdt_[varIndex],
            oldImplicitDdt_[varIndex],
            mesh_.V().field(),
            implicitScale_[varIndex]
        );
        residuals[varIndex] = implicitDdtResidual_[varIndex]
            - implicitConvergence_[varIndex];
        if (returnMe.second() < residuals[varIndex])
        {
            returnMe.first() = admVars_.implicit(varIndex).globalIndex();
            returnMe.second() = residuals[varIndex];
        }
    }
    if (outputImplicitLoopDetails_)
    {
        Info << "implicitLoopDetails: residuals = " << residuals << endl;
    }
    return returnMe;
}


template <int matrixSize>
Tuple2<label, scalar> Foam::craftsModel<matrixSize>::testUdfDeltaConvergence()
    const
{
    Tuple2<label, scalar> returnMe(-1, -VGREAT);
    forAll(udfDeltaResidual_, udfIndex)
    {
        label varIndex
        (
            admVars_.changedByUdf(udfIndex).localIndex()
        );
        udfDeltaResidual_[udfIndex] = calculateRmsError
        (
            udfDelta_[udfIndex],
            oldUdfDelta_[udfIndex],
            mesh_.V().field(),
            standardScale_[varIndex]
        );
        scalar testedValue
        (
            udfDeltaResidual_[udfIndex] - standardConvergence_[varIndex]
        );
        if (returnMe.second() < testedValue)
        {
            returnMe.first() = admVars_.changedByUdf(udfIndex).globalIndex();
            returnMe.second() = testedValue;
        }
    }
    return returnMe;
}


template <int matrixSize>
Foam::scalar Foam::craftsModel<matrixSize>::calculateNextBestDeltaT() const
{
    if (!atsUseAts_)
    {
        return runTime_.deltaT().value();
    }
    
    label limitingVariable(-1);
    scalar reactionNextDeltaT(VGREAT);
    forAll(admVars_.standard(), varIndex)
    {
        scalar dtEst
        (
            calculateNextDeltaT
            (
                coarseStepStandardSaveSpot_[varIndex],
                admVars_.standard(varIndex).evaluateField(),
                standardConvergence_[varIndex] * atsConvergenceFactor_,
                mesh_.V().field(),
                runTime_.deltaT().value(),
                standardScale_[varIndex]
            )
        );
        if (reactionNextDeltaT > dtEst)
        {
            reactionNextDeltaT = dtEst;
            limitingVariable = admVars_.standard(varIndex).globalIndex();
        }
    }
    forAll(admVars_.implicit(), varIndex)
    {
        scalar dtEst
        (
            calculateNextDeltaT
            (
                coarseStepImplicitSaveSpot_[varIndex],
                admVars_.implicit(varIndex).evaluateField(),
                implicitConvergence_[varIndex] * atsConvergenceFactor_,
                mesh_.V().field(),
                runTime_.deltaT().value(),
                implicitScale_[varIndex]
            )
        );
        if (reactionNextDeltaT > dtEst)
        {
            reactionNextDeltaT = dtEst;
            limitingVariable = admVars_.implicit(varIndex).globalIndex();
        }
    }

    // Ask the flow model for its suggestions
    scalar flowNextBestTime(flow_->calculateNextBestDeltaT());

    if (flow_->subStepping())
    {
        label nextBestSubSteps
        (
            flow_->nSubSteps() * reactionNextDeltaT / flowNextBestTime
//            reactionNextDeltaT / flowNextBestTime
        );
        label lastSubSteps(flow_->nSubSteps());
        label newSubSteps(flow_->setNSubSteps(nextBestSubSteps));
        flowNextBestTime *= newSubSteps / lastSubSteps;
        if (flow_->outputFlowTimestepEstimate())
        {
            Info << "flowTimestepEstimate: subSteps changed from "
                << lastSubSteps << " to " << newSubSteps << endl;
        }
    }
    else if (flow_->outputFlowTimestepEstimate())
    {
        Info << "flowTimestepEstimate: subStepping disabled" << endl;
    }

    scalar nextDeltaT(min(flowNextBestTime, reactionNextDeltaT));

    // Performance feedback
    scalar pFactor(1.0);
    if (atsPfUsePf_)
    {
        scalar coarse(0);
        scalar fine(0);
        if (atsPfMeasure_ == 1)
        {
            coarse = coarseStepIterations_;
            fine = fineStepIterations_;
        }
        else
        {
            coarse = coarseStepCpuTime_;
            fine = fineStepCpuTime_;
        }
        scalar gamma(1.0);
        if (coarse > 0)
        {
            gamma = fine / coarse;
        }
        if (gamma < 1.0)
        {
            if (mag(atsPfBias_) < SMALL)
            {
                pFactor = gamma;
            }
            else
            {
                pFactor = (1 - exp(-gamma * atsPfBias_))
                        / (1 - exp(-atsPfBias_));
            }
        }
    }

    // Time step growth / reduction limits
    // Ensure it doesn't grow too fast
    scalar maxDeltaT(runTime_.deltaT().value() * atsMaxIncreaseFactor_);

    // Ensure it doesn't drop too quickly
    scalar minDeltaT(runTime_.deltaT().value() * atsMaxReductionFactor_);

    // Apply the limits
    nextDeltaT = min(maxDeltaT, nextDeltaT);
    nextDeltaT = max(minDeltaT, nextDeltaT);

    // Apply performance factor
    nextDeltaT =
        pFactor * nextDeltaT
      + (1 - pFactor) * runTime_.deltaT().value() / 2;

    // Apply overclock factor
    nextDeltaT *= atsOverclockFactor_;

    // Ensure it won't exceed the end time ( divide by 2 because we take double
    // steps )
    maxDeltaT = min
    (
        maxDeltaT,
        (runTime_.endTime().value() - runTime_.value()) / 2
    );
    
    // Output results to console
    if (outputReactionTimestepEstimate_)
    {
        Info << "reactionTimestepEstimate: "
            << admVars_.all(limitingVariable).name() << " dtEst = "
            << reactionNextDeltaT;

        if (flowNextBestTime < reactionNextDeltaT)
        {
            Info << ", flow model dtEst = " << flowNextBestTime
                << " (bottleneck)";
        }
        else
        {
            Info << " (bottleneck), flow model dtEst = " << flowNextBestTime;
        }
        Info << ", pFactor = " << pFactor << endl;
    }

    return nextDeltaT;
}


template <int matrixSize>
void Foam::craftsModel<matrixSize>::transitionToNextTimestep()
{
    transitionToNextTimestep_ = true;
    if (flow_->subStepping())
    {
        // Copy state 3->2, 1->0 (currently on state 1)
        saveState(0);
        loadFlowState(3);
        saveFlowState(2);
        loadState(0);
    }
}


template <int matrixSize>
Foam::label Foam::craftsModel<matrixSize>::coarseStep()
{
    label returnMe(stepFlowDispatch(COARSE)); // increments runTime
    if (returnMe)
    {
        return returnMe;
    }

    return stepReactionsDispatch(COARSE);
}


template <int matrixSize>
Foam::label Foam::craftsModel<matrixSize>::doubleFineStep()
{
    // First step
    label returnMe(stepFlowDispatch(FINEA)); // increments runTime
    if (returnMe)
    {
        return returnMe;
    }

    returnMe = stepReactionsDispatch(FINEA);
    if (returnMe)
    {
        return returnMe;
    }

    // Second step
    returnMe = stepFlowDispatch(FINEB); // increments runTime
    if (returnMe)
    {
        return returnMe;
    }

    return (stepReactionsDispatch(FINEB));
}


template <int matrixSize>
Foam::label Foam::craftsModel<matrixSize>::stepFlowDispatch(stepType st)
{
    if (flow_->subStepping())
    {

        // Load subStepping timeframe, either slot 1 (restart), or slot 3
        // (continue)
        scalar oldDeltaT(runTime_.deltaT().value());
        scalar newDeltaT(oldDeltaT / flow_->nSubSteps());
        if (transitionToNextTimestep_ || st == FINEB)
        {
            loadFlowState(3);
        }
        else
        {
            loadFlowState(2);
        }
        runTime_.setDeltaT(newDeltaT);

        // Step to the middle
        for (label count(0); count < (flow_->nSubSteps() / 2); count++)
        {
            runTime_.plusPlusNoOutput();
            label returnMe(flow_->step()); //&&&
            if (returnMe)
            {
                runTime_.loadState(0);
                runTime_.setDeltaT(oldDeltaT);
                return returnMe;
            }
        }

        // Store the solution (static velocity field)
        flow_->storeFlowSolution();

        // Step to end
        for (label count(0); count < (flow_->nSubSteps() / 2); count++)
        {
            runTime_.plusPlusNoOutput();
            label returnMe(flow_->step()); //&&&
            if (returnMe)
            {
                runTime_.loadState(0);
                runTime_.setDeltaT(oldDeltaT);
                return returnMe;
            }
        }

        switch (st)
        {
            case COARSE:
                flow_->storeCoarseSolution();
            break;
            case FINEB:
                flow_->storeFineSolution();
                // carry through to next case
            case FINEA:
                saveFlowState(3);
        }        

        if (transitionToNextTimestep_)
        {
            loadState(1);
        }
        else
        {
            loadState(0);
        }
        runTime_.setDeltaT(oldDeltaT);

        if (st == FINEA)
        {
            transitionToNextTimestep_ = false;
        }

        // Increment runTime by a full (non-sub) step depending on st
        if (st != FINEB)
        {
            runTime_.plusPlusNoOutput();
        }
        else
        {
            runTime_++;
        }

        flow_->retagFlowSolution();
        return 0;
    }
    else
    {
        if (st != FINEB)
        {
            runTime_.plusPlusNoOutput();
        }
        else
        {
            runTime_++;
        }
        
        label returnMe(flow_->step());
        if (returnMe)
        {
            return returnMe;
        }
        
        flow_->storeFlowSolution();
        
        if (st == COARSE)
        {
            flow_->storeCoarseSolution();
        }
        return returnMe;
    }
}


template <int matrixSize>
Foam::label Foam::craftsModel<matrixSize>::stepReactionsDispatch(stepType st)
{
    label returnMe(stepReactions());
    
    if (st == COARSE)
    {
        // Store coarse performance, reset fine performance
        coarseStepIterations_ = lastStepIterations_;
        coarseStepCpuTime_ = lastStepCpuTime_;
        fineStepIterations_ = 0;
        fineStepCpuTime_ = 0.0;
        
        // Store interim results
        admVars_.saveStandardInternalFields(coarseStepStandardSaveSpot_);
        admVars_.saveImplicitInternalFields(coarseStepImplicitSaveSpot_);
    }
    else
    {
        // Store fine performance
        fineStepIterations_ += lastStepIterations_;
        fineStepCpuTime_ += lastStepCpuTime_;
        if (flow_->subStepping())
        {
            saveState(1);
        }
    }
    
    return returnMe;
}


template <int matrixSize>
Foam::label Foam::craftsModel<matrixSize>::stepReactions()
{
    // Algorithm pseudo code
    // Guess implicit ddt (equal to previous timestep or zero)
    // (Field loop begin)
    // while (standardVariables differ from previous)
        // Create coupled matrix, source array, variable array
        // Map variables into variable array
        // Fill in terms matrix & source
        // Solve the matrix
        // Map variables back from variable array
        // Apply variable limits
        // Update coeffs and derived
        // (Implicit loop begin)
        // For each cell in the mesh:
            // Call implicit function hook
            // Update coeffs and derived variables
        // Store old implicit ddt
        // Calculate new implicit ddt
        // Calculate implicit residual
/*fileName outVariables(runTime_.path()/"outVariables.csv");
OFstream osv(outVariables);
fileName matrixMFile(runTime_.path()/"matrixBlockM");
OFstream osm(matrixMFile);
fileName matrixXFile(runTime_.path()/"matrixBlockX");
OFstream osx(matrixXFile);
fileName matrixBFile(runTime_.path()/"matrixBlockB");
OFstream osb(matrixBFile);
reportVariablesHeader(osv);*/

    // Reset last step iterations and time
    lastStepIterations_ = 0;
    lastStepCpuTime_ = 0.0;
    cpuTime timer;

    // Call function hook 'initializeTimestep'
    hooks_->initializeTimestep();

    // Outer loop (repeat until standard variables stop changing)
    for
    (
        label outerIteration(0);
        outerIteration < outerLoopMaxIterations_;
        outerIteration++
    )
    {
        // Save all variable fields temporarily to saveSpot_
        admVars_.saveStandardInternalFields(internalStandardSaveSpot_);

/*if (debug)
{
    osv << "iteration start " << outerIteration << endl;
    for (label i(0); i < mesh_.nCells(); i++)
    {
        reportVariables(osv, i);
    }
}*/
        // Create coupled matrix, source array, variable array
        BlockLduMatrix<vectorType> blockM(mesh_);

        Field<tensorType>& d = blockM.diag().asSquare();
        d = tensorType::zero;

        Field<vectorType>& l = blockM.lower().asLinear();
        Field<vectorType>& u = blockM.upper().asLinear();
        l = vectorType::zero;
        u = vectorType::zero;

        Field<vectorType> blockX(mesh_.nCells(), vectorType::zero);
        Field<vectorType> blockB(mesh_.nCells(), vectorType::zero);
/*if (debug > 1)
{
    osm << "Initialized" << endl;
    osx << "Initialized" << endl;
    osb << "Initialized" << endl;
    osm << "blockM = [" << blockM << "]" << endl;
    osx << "blockX = [" << blockX << "]" << endl;
    osb << "blockB = [" << blockB << "]" << endl;
}*/

        // Fill in terms of the matrix
        fillTransportTerms(blockM, blockX, blockB);
/*if (debug > 1)
{
    osm << "fillTransportTerms" << endl;
    osx << "fillTransportTerms" << endl;
    osb << "fillTransportTerms" << endl;
    osm << "blockM = [" << blockM << "]" << endl;
    osx << "blockX = [" << blockX << "]" << endl;
    osb << "blockB = [" << blockB << "]" << endl;
}*/

        fillImplicitTerms(blockM, blockX, blockB);

/*if (debug > 1)
{
    osm << "fillImplicitTerms" << endl;
    osx << "fillImplicitTerms" << endl;
    osb << "fillImplicitTerms" << endl;
    osm << "blockM = [" << blockM << "]" << endl;
    osx << "blockX = [" << blockX << "]" << endl;
    osb << "blockB = [" << blockB << "]" << endl;
}*/

        fillSDiagReactionTerms(blockM, blockB);

/*if (debug > 1)
{
    osm << "fillSDiagReactionTerms" << endl;
    osx << "fillSDiagReactionTerms" << endl;
    osb << "fillSDiagReactionTerms" << endl;
    osm << "blockM = [" << blockM << "]" << endl;
    osx << "blockX = [" << blockX << "]" << endl;
    osb << "blockB = [" << blockB << "]" << endl;
}*/

        fillSourceReactionTerms(blockM, blockB);

/*if (debug > 1)
{
    osm << "fillSourceReactionTerms" << endl;
    osx << "fillSourceReactionTerms" << endl;
    osb << "fillSourceReactionTerms" << endl;
    osm << "blockM = [" << blockM << "]" << endl;
    osx << "blockX = [" << blockX << "]" << endl;
    osb << "blockB = [" << blockB << "]" << endl;
}*/

        fillOssNssNsiReactionTerms(blockM, blockB);

/*if (debug > 1)
{
    osm << "fillOssNssNsiReactionTerms" << endl;
    osx << "fillOssNssNsiReactionTerms" << endl;
    osb << "fillOssNssNsiReactionTerms" << endl;
    osm << "blockM = [" << blockM << "]" << endl;
    osx << "blockX = [" << blockX << "]" << endl;
    osb << "blockB = [" << blockB << "]" << endl;
}*/

        //- Block coupled solver call
        BlockSolverPerformance<vectorType> solverPerf =
            BlockLduSolver<vectorType>::New
            (
                word("blockVar"),
                blockM,
                mesh_.solver("blockVar")
            )->solve(blockX, blockB);

        if (outputReactionSolverPerformance_)
        {
            solverPerf.print();
        }
//osx << "\nAfter solution" << endl;
//osx << blockX << endl;
//osx << max(blockX) << endl;

        // Retrieve solution for standard variables
        forAll(admVars_.standard(), varIndex)
        {
            blockMatrixTools::blockRetrieve
            (
                varIndex,
                admVars_.standard(varIndex)().internalField(),
                blockX
            );
        }

        // Retrieve solution for implicit variables (&&& is this necessary?)
        forAll(admVars_.implicit(), varIndex)
        {
            blockMatrixTools::blockRetrieve
            (
                varIndex + admVars_.nStandard(),
                admVars_.implicit(varIndex)().internalField(),
                blockX
            );
        }

        // Ensure the variables are all still within bounds
        hooks_->applyVariableLimits();

        // Add the effect of implicit routines from the previous iteration to
        // the variables
        addUdfDelta();

/*if (debug)
{
    osv << "crm solved" << endl;
    for (label i(0); i < mesh_.nCells(); i++)
    {
        reportVariables(osv, i);
    }
}*/

        // Variables have changed - set milestone (forces lazy evaluation
        // objects to recalculate)
        setMilestone();

        // Save current state of standard variables marked with changedByUdf
        savePreUdf();

        // Update implicit variables using implicitAutoSolve and user-
        // defined functions, loop until they stabilize
        label impResult(runImplicitRoutines());
/*if (debug)
{
    osv << "implicitRoutines result " << label(impResult) << endl;
    for (label i(0); i < mesh_.nCells(); i++)
    {
        reportVariables(osv, i);
    }
}*/
        switch (impResult)
        {
            default:
                // Unknown value
                FatalErrorIn("craftsModel::stepReactions")
                    << "Received error level " << impResult << " from "
                    << "runImplicitRoutines. Valid error levels are:"
                    << token::NL << token::TAB
                    << " 0 = success"
                    << token::NL << token::TAB
                    << "-1 = total failure, reduce delta t"
                    << token::NL << token::TAB
                    << "-2 = convergence failure, repeat timestep"
                    << abort(FatalError);
            case 0:
                // success - keep going
                break;
            case -1:
                // total failure - return to solver level
                return impResult;
            case -2:
                // convergence failure - keep going, then repeat
                break;
        }

        // Save old udfDelta and calculate the new one
        oldUdfDelta_ = udfDelta_;
        calculateUdfDelta();

        // Check for convergence
        scalar convergenceTest
        (
            max
            (
                testStandardConvergence(internalStandardSaveSpot_).second(),
                testUdfDeltaConvergence().second()
            )
        );
        
        Info << "overallSummary: t = " << runTime_.timeName() << ", dt = "
            << runTime_.deltaT().value() << ", i = "
            << outerIteration << ", res = " << convergenceTest
            << ", imp = " << label(impResult) << endl;
/*if (debug > 1)
{
    FatalErrorIn("craftsModel") << "Debug stop point" << abort(FatalError);
}*/

        if ((convergenceTest <= 0) && (impResult == 0))
        {
            hooks_->finalizeTimestep();
            lastStepCpuTime_ = timer.cpuTimeIncrement();
            lastStepIterations_ += outerIteration;
            return 0;
        }
    }
    WarningIn("craftsModel::step")
        << "Exceeded max iterations of " << outerLoopMaxIterations_ << endl;
    return -1;
}


template <int matrixSize>
Foam::label Foam::craftsModel<matrixSize>::runImplicitRoutines()
{
    // Call function hook 'initializeImplicitLoop'
    hooks_->initializeImplicitLoop();

    // Implicit loop (inner)
    for
    (
        label implicitIter(0);
        implicitIter < innerLoopMaxIterations_;
        implicitIter++
    )
    {
        // Call function hook 'implicitLoop'
        label impResult(hooks_->implicitLoop());
        
        switch (impResult)
        {
            default:
                // Unknown value
                FatalErrorIn("craftsModel::runImplicitRoutines")
                    << "Received error level " << impResult << " from "
                    << "function hooks implicitLoop. Valid error levels are:"
                    << token::NL << token::TAB
                    << " 0 = success"
                    << token::NL << token::TAB
                    << "-1 = total failure, reduce delta t"
                    << token::NL << token::TAB
                    << "-2 = convergence failure, repeat timestep"
                    << token::NL << token::TAB
                    << "-3 = stabilization failure, repeat implicit loop"
                    << abort(FatalError);
            case  0: // success, keep going
            case -3: // stabilization failure - repeat implicit loop
                break;
            case -1: // total failure - return to solver level
            case -2: // convergence failure - repeat timestep
                if (outputImplicitLoopSummary_)
                {
                    Info << "implicitLoopSummary: exiting (fail) at n = "
                        << implicitIter << endl;
                }
                return impResult;
        }

        // Store the old implicitDdt, and calculate the new one
        oldImplicitDdt_ = implicitDdt_;
        getImplicitDdt();

        if ((testImplicitDdtConvergence().second() <= 0) && impResult == 0)
        {
            // Call function hook 'finalizeImplicitLoop'
            hooks_->finalizeImplicitLoop();

            if (outputImplicitLoopSummary_)
            {
                Info << "implicitLoopSummary: exiting (success) at n = "
                    << implicitIter << endl;
            }

            return 0;
        }
    } // end implicit loop

    // Max iterations exceeded: warn, and return total failure
    WarningIn("craftsModel::runImplicitRoutines")
        << "Exceeded max iterations of " << innerLoopMaxIterations_
        << endl;
    return -1;
}


template <int matrixSize>
void Foam::craftsModel<matrixSize>::saveFlowState(const label slot)
{
    flow_->saveState(slot);
    runTime_.saveState(slot);
}


template <int matrixSize>
void Foam::craftsModel<matrixSize>::loadFlowState(const label slot)
{
    flow_->loadState(slot);
    runTime_.loadState(slot);
}


template <int matrixSize>
void Foam::craftsModel<matrixSize>::saveState(const label slot)
{
    admReactionReader::saveState(slot);
    flow_->saveState(slot);
    hooks_->saveState(slot);
}


template <int matrixSize>
void Foam::craftsModel<matrixSize>::clearState(const label slot)
{
    admReactionReader::clearState(slot);
    hooks_->clearState(slot);
}


template <int matrixSize>
void Foam::craftsModel<matrixSize>::loadState(const label slot)
{
    admReactionReader::loadState(slot);
    flow_->loadState(slot);
    hooks_->loadState(slot);
}


template <int matrixSize>
Foam::label Foam::craftsModel<matrixSize>::nStates() const
{
    return admReactionReader::nStates();
}


template <int matrixSize>
bool Foam::craftsModel<matrixSize>::validState(const label slot) const
{
    return admReactionReader::validState(slot);
}


template <int matrixSize>
void Foam::craftsModel<matrixSize>::reportScales(Ostream& os)
{
    os << runTime_.value() << "," << "scale,";
    if (admVars_.nStandard())
    {
        os << standardScale_[0];
    }
    for (label varIndex(1); varIndex < admVars_.nStandard(); varIndex++)
    {
        os << "," << standardScale_[varIndex];
    }
    if (admVars_.nStandard() && admVars_.nImplicit())
    {
        os << ",";
    }
    if (admVars_.nImplicit())
    {
        os << implicitScale_[0];
    }
    for (label varIndex(1); varIndex < admVars_.nImplicit(); varIndex++)
    {
        os << "," << implicitScale_[varIndex];
    }
    os << endl;
}


// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

#include "craftsModelCreateIndices.C"
#include "craftsModelFillTerms.C"
#include "craftsModelImplicitAutoSolve.C"

// ************************************************************************* //
