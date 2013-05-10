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

#include "craftsFlow.H"
#include "admReactionReader.H"

// * * * * * * * * * * * * * * Static Data Members * * * * * * * * * * * * * //

namespace Foam
{
    defineTypeNameAndDebug(craftsFlow, 0);

    defineRunTimeSelectionTable(craftsFlow, nameConstructor);
}


// * * * * * * * * * * * * * Private Member Functions  * * * * * * * * * * * //


void Foam::craftsFlow::readDict()
{
    // Set subStepping type based on the word
    if (subSteppingTypeWord_ == "off")
    {
        subSteppingType_ = OFF;
    }
    else if (subSteppingTypeWord_ == "fixedNSteps")
    {
        subSteppingType_ = FIXED_N_STEPS;
        ssCurrentNSubSteps_ = readLabel
        (
            settingsDict_.subDict("subStepping").lookup("nSteps")
        );
        if (ssCurrentNSubSteps_ % 2)
        {
            WarningIn("craftsFlow::readDict")
                << "nSteps must be even.  Using "
                << ssCurrentNSubSteps_ + 1 << " in place of "
                << ssCurrentNSubSteps_ << endl;
            
            ssCurrentNSubSteps_++;
        }
    }
    else if (subSteppingTypeWord_ == "fixedTimestep")
    {
        subSteppingType_ = FIXED_TIMESTEP;
    }
    else if (subSteppingTypeWord_ == "adaptiveTimestep")
    {
        subSteppingType_ = ADAPTIVE_TIMESTEP;
        if (!atsEnableWithFlowModel_)
        {
            FatalIOErrorIn("craftsFlow::readDict", settingsDict_)
                << "enableWithFlowModel must be 'yes' when subStepping type "
                << "is 'adaptiveTimestep"
                << exit(FatalIOError);
        }
    }
    else
    {
        FatalIOErrorIn
        (
            "craftsFlow::readDict",
            settingsDict_.subDict("subStepping")
        )
            << subSteppingTypeWord_ << " unrecognized flow model subStepping "
            << "behaviour.  Expecting: 'off', 'fixedNSteps', 'fixedTimeStep', "
            << "'adaptiveTimeStep'."
            << exit(FatalIOError);
    }

    // Adaptive timestepping tolerance
    if (settingsDict_.subDict("adaptiveTimestepping").found("tolerance"))
    {
        const dictionary& toleranceDict
        (
            settingsDict_.subDict("adaptiveTimestepping").subDict("tolerance")
        );
        if (toleranceDict.found("phi")
        )
        {
            atsPhiTolerance_ = readScalar(toleranceDict.lookup("phi"));
        }
        if (toleranceDict.found("p"))
        {
            atsPTolerance_ = readScalar(toleranceDict.lookup("p"));
        }
    }

    // Adaptive timestepping minimum error scale
    if (settingsDict_.subDict("adaptiveTimestepping").found("minErrorScale"))
    {
        const dictionary& minErrorDict
        (
            settingsDict_
                .subDict("adaptiveTimestepping")
                .subDict("minErrorScale")
        );
        if (minErrorDict.found("phi"))
        {
            atsPhiMinScale_ = readScalar(minErrorDict.lookup("phi"));
        }
        if (minErrorDict.found("p"))
        {
            atsPMinScale_ = readScalar(minErrorDict.lookup("p"));
        }
    }

    // Search for ignored fields
    forAll(atsIgnore_, i)
    {
        if (atsIgnore_[i] == "phi")
        {
            atsIgnorePhi_ = true;
        }
        if (atsIgnore_[i] == "p")
        {
            atsIgnoreP_ = true;
        }
    }

    // SubStepping error checking
    if
    (
        (ssMaxDeltaT_ <= ssMinDeltaT_)
     || (ssTargetDeltaT_ > ssMaxDeltaT_)
     || (ssTargetDeltaT_ < ssMinDeltaT_)
    )
    {
        FatalIOErrorIn
        (
            "craftsFlow::readDict",
            settingsDict_.subDict("subStepping")
        )
            << "maxDeltaT (" << ssMaxDeltaT_ << ") must be >= "
            << "targetDeltaT (" << ssTargetDeltaT_ << ") must be >= "
            << "minDeltaT (" << ssMinDeltaT_ << ")"
            << exit(FatalIOError);
    }
    if (ssInitialNSubSteps_ % 2)
    {
        // value is odd
        WarningIn("craftsFlow::readDict")
            << "initialNSteps must be even.  Using " << ssInitialNSubSteps_ + 1
            << " in place of " << ssInitialNSubSteps_ << endl;
        
        ssInitialNSubSteps_++;
        ssCurrentNSubSteps_ = ssInitialNSubSteps_;
    }
    if ((ssMaxNSubSteps_ < 2) || (ssMaxNSubSteps_ % 2))
    {
        FatalIOErrorIn
        (
            "craftsFlow::readDict",
            settingsDict_.subDict("subStepping")
        )
            << "maxNSteps = " << ssMaxNSubSteps_ << ", must be an even "
            << "integer greater than 1."
            << exit(FatalIOError);
    }
    if ((ssMinNSubSteps_ < 2) || (ssMinNSubSteps_ % 2))
    {
        FatalIOErrorIn
        (
            "craftsFlow::readDict",
            settingsDict_.subDict("subStepping")
        )
            << "minNSteps = " << ssMinNSubSteps_ << ", must be an even "
            << "integer greater than 1."
            << exit(FatalIOError);
    }
    if (subStepping_)
    {
        IOWarningIn("craftsFlow::readDict", settingsDict_)
            << "subStepping is enabled. Some boundary conditions look for "
            << "'U' or 'phi' by default.  Direct flow-related boundary "
            << "conditions to 'USub' or 'phiSub', and reaction-related "
            << "boundary conditions to 'UGlobal' or 'phiGlobal' instead by "
            << "adding 'U' or 'phi' keywords in their dictionaries."
            << endl;
    }
    else
    {
        IOWarningIn("craftsFlow::readDict", settingsDict_)
            << "subStepping is disabled. Some boundary conditions look for "
            << "'U' or 'phi' by default.  Direct them all to 'UGlobal' or "
            << "'phiGlobal' instead by adding 'U' or 'phi' keywords in their "
            << "dictionaries."
            << endl;
    }

}


// * * * * * * * * * * * * * * * * Constructors  * * * * * * * * * * * * * * //

Foam::craftsFlow::craftsFlow
(
    volVectorField& U,
    surfaceScalarField& phi,
    volScalarField& p,
    admReactionReader& model
)
:
    model_(model),
    runTime_(model.runTime()),
    mesh_(model.mesh()),
    UGlobal_(U),
    phiGlobal_(phi),
    pGlobal_(p),
    gammaTurbulent_
    (
        IOobject
        (
            "gammaTurbulent",
            runTime_.timeName(),
            mesh_,
            IOobject::READ_IF_PRESENT,
            IOobject::AUTO_WRITE
        ),
        mesh_,
        dimensionedScalar("gammaTurbulent", dimensionSet(0, 2, -1, 0, 0), 0)
    ),
    settingsDict_
    (
        model.admSettingsDict()
            .subDict("flowModel")
    ),
    coarseSaveSpotPhi_(phiGlobal_.internalField().size(), 0.0),
    coarseSaveSpotP_(pGlobal_.internalField().size(), 0.0),

    atsEnableWithFlowModel_
    (
      settingsDict_.subDict("adaptiveTimestepping")
        .lookup("enableWithFlowModel")
    ),

    atsIgnore_
    (
        settingsDict_.subDict("adaptiveTimestepping").found("ignore")
      ? settingsDict_.subDict("adaptiveTimestepping").lookup("ignore")
      : wordList(0)
    ),

    atsIgnorePhi_(false),
    atsIgnoreP_(false),

    atsPhiTolerance_(-VGREAT),
    atsPTolerance_(-VGREAT),

    atsPhiMinScale_(SMALL),
    atsPMinScale_(SMALL),

    atsPhiScale_(VSMALL),
    atsPScale_(VSMALL),
    
    subSteppingTypeWord_
    (
        settingsDict_.subDict("subStepping").lookup("type")
    ),
    subSteppingType_(OFF),
    subStepping_(subSteppingTypeWord_ != "off"),
    ssMaxDeltaT_
    (
        settingsDict_.subDict("subStepping").found("maxDeltaT")
      ? readScalar(settingsDict_.subDict("subStepping").lookup("maxDeltaT"))
      : VGREAT
    ),
    ssMinDeltaT_
    (
        settingsDict_.subDict("subStepping").found("minDeltaT")
      ? readLabel(settingsDict_.subDict("subStepping").lookup("minDeltaT"))
      : scalar(0)
    ),
    ssTargetDeltaT_
    (
        settingsDict_.subDict("subStepping").found("targetDeltaT")
      ? readScalar(settingsDict_.subDict("subStepping").lookup("targetDeltaT"))
      : ssMaxDeltaT_ / 2.0 + ssMinDeltaT_ / 2.0
    ),
    ssMaxNSubSteps_
    (
        settingsDict_.subDict("subStepping").found("maxNSteps")
      ? readLabel(settingsDict_.subDict("subStepping").lookup("maxNSteps"))
      : 2000000000
    ),
    ssMinNSubSteps_
    (
        settingsDict_.subDict("subStepping").found("minNSteps")
      ? readLabel(settingsDict_.subDict("subStepping").lookup("minNSteps"))
      : label(2)
    ),
    ssInitialNSubSteps_
    (
        settingsDict_.subDict("subStepping").found("initialNSteps")
      ? readLabel(settingsDict_.subDict("subStepping").lookup("initialNSteps"))
      : label(2)
    ),
    ssCurrentNSubSteps_
    (
        subStepping_
      ? ssInitialNSubSteps_
      : label(1)
    ),

    USubPtr_
    (
        subStepping_
      ? new volVectorField("USub", UGlobal_)
      : 0
    ),
    phiSubPtr_
    (
        subStepping_
      ? new surfaceScalarField("phiSub", phiGlobal_)
      : 0
    ),
    pSubPtr_
    (
        subStepping_
      ? new volScalarField("pSub", pGlobal_)
      : 0
    ),

    UTmpPtr_
    (
        subStepping_
      ? new volVectorField("UTmp", UGlobal_)
      : 0
    ),
    phiTmpPtr_
    (
        subStepping_
      ? new surfaceScalarField("phiTmp", phiGlobal_)
      : 0
    ),
    pTmpPtr_
    (
        subStepping_
      ? new volScalarField("pTmp", pGlobal_)
      : 0
    ),

    fineSaveSpotPhiPtr_
    (
        subStepping_
      ? new scalarField(phiGlobal_.internalField().size(), 0.0)
      : 0
    ),
    fineSaveSpotPPtr_
    (
        subStepping_
      ? new scalarField(pGlobal_.internalField().size(), 0.0)
      : 0
    ),

    UActive_
    (
        subStepping_
      ? USubPtr_()
      : UGlobal_
    ),
    phiActive_
    (
        subStepping_
      ? phiSubPtr_()
      : phiGlobal_
    ),
    pActive_
    (
        subStepping_
      ? pSubPtr_()
      : pGlobal_
    ),

    outputFlowSolverPerformance_
    (
        model_.outputFlagsDict().found("flowSolverPerformance")
      ? bool(Switch(model_.outputFlagsDict().lookup("flowSolverPerformance")))
      : false
    ),
    outputFlowResiduals_
    (
        model_.outputFlagsDict().found("flowResiduals")
      ? bool(Switch(model_.outputFlagsDict().lookup("flowResiduals")))
      : true
    ),
    outputFlowErrorScales_
    (
        model_.outputFlagsDict().found("flowErrorScales")
      ? bool(Switch(model_.outputFlagsDict().lookup("flowErrorScales")))
      : false
    ),
    outputFlowTimestepEstimate_
    (
        model_.outputFlagsDict().found("flowTimestepEstimate")
      ? bool(Switch(model_.outputFlagsDict().lookup("flowTimestepEstimate")))
      : true
    ),
    outputFlowContinuityErrors_
    (
        model_.outputFlagsDict().found("flowContinuityErrors")
      ? bool(Switch(model_.outputFlagsDict().lookup("flowContinuityErrors")))
      : true
    )
{
    readDict();
}

// * * * * * * * * * * * * * * Member Functions  * * * * * * * * * * * * * * //

Foam::label Foam::craftsFlow::setNSubSteps(label newNSubSteps) const
{
    // Force to even
    newNSubSteps += newNSubSteps % 2;
    
    // Apply limits
    ssCurrentNSubSteps_ = max(ssMinNSubSteps_, newNSubSteps);
    ssCurrentNSubSteps_ = min(ssMaxNSubSteps_, ssCurrentNSubSteps_);

    return ssCurrentNSubSteps_;
}


Foam::label Foam::craftsFlow::step()
{
    notImplemented("craftsFlow::step");
    return -1;
}


void Foam::craftsFlow::storeCoarseSolution()
{
    coarseSaveSpotPhi_ = phiActive_.internalField();
    coarseSaveSpotP_ = pActive_.internalField();
}


void Foam::craftsFlow::storeFineSolution()
{
    if (!subStepping_)
    {
        FatalErrorIn("craftsFlow::storeFineSolution")
            << "This function is invalid when sub-stepping is disabled."
            << abort(FatalError);
    }
    fineSaveSpotPhiPtr_() = phiActive_.internalField();
    fineSaveSpotPPtr_() = pActive_.internalField();
}


void Foam::craftsFlow::storeFlowSolution()
{
    if (subStepping_)
    {
        UTmpPtr_() == USubPtr_();
        phiTmpPtr_() == phiSubPtr_();
	    pTmpPtr_() == pSubPtr_();
    }
}


void Foam::craftsFlow::retagFlowSolution()
{
    UGlobal_ == UTmpPtr_();
    phiGlobal_ == phiTmpPtr_();
    pGlobal_ == pTmpPtr_();
}


void Foam::craftsFlow::calculateScales()
{
    if (atsEnableWithFlowModel_)
    {
        if (subStepping_)
        {
            if (!atsIgnorePhi_)
            {
                atsPhiScale_ = max
                (
                    model_.calculateRms(fineSaveSpotPhiPtr_()), atsPhiMinScale_
                );
            }
            if (!atsIgnoreP_)
            {
                atsPScale_ = max
                (
                    model_.calculateRms(fineSaveSpotPPtr_()), atsPMinScale_
                );
            }
        }
        else
        {
            if (!atsIgnorePhi_)
            {
                atsPhiScale_ = max
                (
                    model_.calculateRms(phiActive_), atsPhiMinScale_
                );
            }
            if (!atsIgnoreP_)
            {
                atsPScale_ = max
                (
                    model_.calculateRms(pActive_), atsPMinScale_
                );
            }
        }
        if (outputFlowErrorScales_)
        {
            Info << "flowErrorScales: U = " << atsPhiScale_ << ", p = "
                << atsPScale_ << endl;
        }
    }
    else if (outputFlowErrorScales_)
    {
        Info << "flowErrorScales: n/a (flow adaptiveTimestepping disabled)"
            << endl;
    }
}


Foam::scalar Foam::craftsFlow::testConvergence() const
{
    scalar returnMe(-VGREAT);
    scalar phiTest(-VGREAT);
    scalar pTest(-VGREAT);

    if (atsEnableWithFlowModel_)
    {
        if (subStepping_)
        {
            if (!atsIgnorePhi_)
            {
                phiTest = model_.calculateRmsError
                (
                    fineSaveSpotPhiPtr_(),
                    coarseSaveSpotPhi_,
                    mesh_.magSf().field(),
                    atsPhiScale_
                ) - atsPhiTolerance_;
            }
            if (!atsIgnoreP_)
            {
                pTest = model_.calculateRmsError
                (
                    fineSaveSpotPPtr_(),
                    coarseSaveSpotP_,
                    mesh_.V().field(),
                    atsPScale_
                ) - atsPTolerance_;
            }
        }
        else
        {
            if (!atsIgnorePhi_)
            {
                phiTest = model_.calculateRmsError
                (
                    phiActive_.internalField(),
                    coarseSaveSpotPhi_,
                    mesh_.magSf().field(),
                    atsPhiScale_
                ) - atsPhiTolerance_;
            }
            if (!atsIgnoreP_)
            {
                pTest = model_.calculateRmsError
                (
                    pActive_.internalField(),
                    coarseSaveSpotP_,
                    mesh_.V().field(),
                    atsPScale_
                ) - atsPTolerance_;
            }
        }
        
        returnMe = max(phiTest, pTest);

        if (outputFlowResiduals_)
        {
            Info << "flowResiduals: U = " << phiTest << ", p = " << pTest;
            if (returnMe > 0)
            {
                Info << " (fail)" << endl;
            }
            else
            {
                Info << " (pass)" << endl;
            }
        }
    }
    else if (outputFlowResiduals_)
    {
        Info << "flowResiduals: n/a (flow adaptiveTimestepping disabled)"
            << endl;
    }
    return returnMe;
}


Foam::scalar Foam::craftsFlow::calculateNextBestDeltaT() const
{
    // Don't have to worry about subStepping here
    if (atsEnableWithFlowModel_)
    {
        scalar dtPhi(VGREAT);
        scalar dtP(VGREAT);

        if (!atsIgnorePhi_)
        {
            dtPhi = admReactionReader::calculateNextDeltaT
            (
                coarseSaveSpotPhi_,
                phiActive_.internalField(),
                atsPhiTolerance_,
                mesh_.magSf().field(),
                runTime_.deltaT().value(),
                atsPhiScale_
            );
        }
        
        if (!atsIgnoreP_)
        {
            dtP = admReactionReader::calculateNextDeltaT
            (
                coarseSaveSpotP_,
                pActive_.internalField(),
                atsPTolerance_,
                mesh_.V().field(),
                runTime_.deltaT().value(),
                atsPScale_
            );
        }
        
        if (outputFlowTimestepEstimate_)
        {
            if (dtPhi <= dtP)
            {
                Info << "flowTimestepEstimate: U dtEst = " << dtPhi
                    << " (bottleneck), p dtEst = " << dtP << endl;
            }
            else
            {
                Info << "flowTimestepEstimate: U dtEst = " << dtPhi
                    << ", p dtEst = " << dtP << " (bottleneck)" << endl;
            }
        }
        return min(dtPhi, dtP);
    }
    else if (outputFlowTimestepEstimate_)
    {
        Info << "flowTimestepEstimate: n/a (flow adaptive timestepping "
            << "disabled)" << endl;
    }
    return VGREAT;
}


void Foam::craftsFlow::adjustSubSteppingAdaptiveTimestep
(
    scalar reactionNextDeltaT,
    scalar& flowNextDeltaT
) const
{
    if (subSteppingType_ == ADAPTIVE_TIMESTEP)
    {
        label nextBestSubSteps
        (
            ssCurrentNSubSteps_ * reactionNextDeltaT / flowNextDeltaT
        );
        label lastSubSteps(ssCurrentNSubSteps_);
        label newSubSteps(setNSubSteps(nextBestSubSteps));
        flowNextDeltaT *= newSubSteps / lastSubSteps;
        if (outputFlowTimestepEstimate_)
        {
            Info << "flowTimestepEstimate: subSteps changed from "
                << lastSubSteps << " to " << newSubSteps << endl;
        }
    }
}


void Foam::craftsFlow::adjustSubSteppingFixedTimestep
(
    scalar nextDeltaT
) const
{
    OStringStream errorMessage;
    bool failed(false);

    if (subSteppingType_ == FIXED_TIMESTEP)
    {
        bool successful(false);
        scalar idealNSteps(nextDeltaT / ssTargetDeltaT_);
        label actualNSteps(idealNSteps + 0.5);        
        scalar actualDeltaT;

        while (!successful && !failed)
        {
            actualNSteps = setNSubSteps(actualNSteps);

            actualDeltaT = nextDeltaT / actualNSteps;
            if (actualDeltaT > ssMaxDeltaT_)
            {
                if (actualNSteps == ssMaxNSubSteps_)
                {
                    errorMessage << "maxNSteps and maxDeltaT limit exceeded";
                    failed = true;
                }
                else
                {
                    actualNSteps += 2;
                }
            }
            else if (actualDeltaT < ssMinDeltaT_)
            {
                if (actualNSteps == ssMinNSubSteps_)
                {
                    errorMessage << "minNSteps and minDeltaT limit exceeded";
                    failed = true;
                }
                else
                {
                    actualNSteps -= 2;
                }
            }
            else
            {
                successful = true;
            }
        }
    }
    
    scalar flowDeltaT(nextDeltaT / ssCurrentNSubSteps_);
    
    if (!failed)
    {
        if (flowDeltaT > ssMaxDeltaT_)
        {
            errorMessage << "maxDeltaT limit exceeded";
            failed = true;
        }
        else if (flowDeltaT < ssMinDeltaT_)
        {
            errorMessage << "minDeltaT limit exceeded";
            failed = true;
        }
    }

    if (failed)
    {
        FatalErrorIn("craftsFlow::adjustSubSteppingFixedTimestep")
            << "Flow model sub-stepping requirements cannot be met:"
            << token::NL << token::TAB
            << "Reason: " << errorMessage.str()
            << token::NL << token::TAB
            << "Reaction deltaT = " << nextDeltaT
            << token::NL << token::TAB
            << "Flow deltaT (min, target, max):"
            << token::NL << token::TAB << token::TAB
            << "(" << ssMinDeltaT_ << ", " << ssTargetDeltaT_ << ", "
            << ssMaxDeltaT_ << ")"
            << token::NL << token::TAB
            << "Sub-steps (min, max):"
            << token::NL << token::TAB << token::TAB
            << "(" << ssMinNSubSteps_ << ", " << ssMaxNSubSteps_ << ")"
            << token::NL << token::TAB
            << ssCurrentNSubSteps_ << " gives flow delta T of "
            << flowDeltaT
            << token::NL << token::NL << token::TAB
            << "Change sub-stepping limits."
            << abort(FatalError);
    }
}


void Foam::craftsFlow::saveState(const label slot)
{
    UGlobal_.saveState(slot);
    phiGlobal_.saveState(slot);
    pGlobal_.saveState(slot);
    if (subStepping_)
    {
        USubPtr_().saveState(slot);
        phiSubPtr_().saveState(slot);
        pSubPtr_().saveState(slot);
    }
}

void Foam::craftsFlow::clearState(const label slot)
{
    // do nothing
}


void Foam::craftsFlow::loadState(const label slot)
{
    UGlobal_.loadState(slot);
    phiGlobal_.loadState(slot);
    pGlobal_.loadState(slot);
    if (subStepping_)
    {
        USubPtr_().loadState(slot);
        phiSubPtr_().loadState(slot);
        pSubPtr_().loadState(slot);
    }
}


Foam::label Foam::craftsFlow::nStates() const
{
    return -1;
}


bool Foam::craftsFlow::validState(const label slot) const
{
    return true;
}

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

#   include "newCraftsFlow.C"

// ************************************************************************* //
