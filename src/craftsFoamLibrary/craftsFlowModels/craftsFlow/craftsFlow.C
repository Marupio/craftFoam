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
    settingsDict_
    (
        model.admSettingsDict()
            .subDict("flowModel")
    ),
    coarseSaveSpotPhi_(phiGlobal_.internalField().size(), 0.0),
    coarseSaveSpotP_(pGlobal_.internalField().size(), 0.0),
    
    phiMinScale_(SMALL),
    pMinScale_(SMALL),
    
    phiScale_(VSMALL),
    pScale_(VSMALL),
    
    phiConvergence_(-VGREAT),
    pConvergence_(-VGREAT),

    subStepping_
    (
        Switch(settingsDict_.subDict("subStepping").lookup("useSubStepping"))
    ),
    
    maxNSubSteps_
    (
        settingsDict_.subDict("subStepping").found("maxNSteps")
      ? readLabel(settingsDict_.subDict("subStepping").lookup("maxNSteps"))
      : FOAM_LABEL_MAX
    ),
    minNSubSteps_
    (
        settingsDict_.subDict("subStepping").found("minNSteps")
      ? readLabel(settingsDict_.subDict("subStepping").lookup("minNSteps"))
      : label(2)
    ),
    initialNSubSteps_
    (
        settingsDict_.subDict("subStepping").found("initialNSteps")
      ? readLabel(settingsDict_.subDict("subStepping").lookup("initialNSteps"))
      : label(2)
    ),
    nSubSteps_
    (
        subStepping_
      ? initialNSubSteps_
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
    if (settingsDict_.found("adaptiveTimeSteppingTolerance"))
    {
        const dictionary& atsDict
        (
            settingsDict_.subDict("adaptiveTimeSteppingTolerance")
        );
        if (atsDict.found("phi"))
        {
            phiConvergence_ = readScalar(atsDict.lookup("phi"));
            if (phiConvergence_ < 0)
            {
                FatalIOErrorIn("craftsFlow::craftsFlow", atsDict)
                    << "phi timestepping tolerance cannot be negative"
                    << exit(FatalIOError);
            }
        }
        if (atsDict.found("p"))
        {
            pConvergence_ = readScalar(atsDict.lookup("p"));
            if (pConvergence_ < 0)
            {
                FatalIOErrorIn("craftsFlow::craftsFlow", atsDict)
                    << "p convergence cannot be negative"
                    << exit(FatalIOError);
            }
        }
    }
    if (settingsDict_.found("minErrorScale"))
    {
        const dictionary& scaleDict(settingsDict_.subDict("minErrorScale"));
        if (scaleDict.found("phi"))
        {
            phiMinScale_ = readScalar(scaleDict.lookup("phi"));
        }
        if (scaleDict.found("p"))
        {
            pMinScale_ = readScalar(scaleDict.lookup("p"));
        }
    }
    
    // SubStepping error checking
    if (initialNSubSteps_ % 2)
    {
        // value is odd
        WarningIn("craftsFlow::craftsFlow")
            << "initialNSteps must be even.  Using " << initialNSubSteps_ + 1
            << " in place of " << initialNSubSteps_ << endl;
        
        initialNSubSteps_++;
        nSubSteps_ = initialNSubSteps_;
    }
    if ((maxNSubSteps_ < 2) || (maxNSubSteps_ % 2))
    {
        FatalIOErrorIn
        (
            "craftsFlow::craftsFlow",
            settingsDict_.subDict("subStepping")
        )
            << "maxNSteps must be an even integer greater than 1."
            << exit(FatalIOError);
    }
    if ((minNSubSteps_ < 2) || (minNSubSteps_ % 2))
    {
        FatalIOErrorIn
        (
            "craftsFlow::craftsFlow",
            settingsDict_.subDict("subStepping")
        )
            << "minNSteps must be an even integer greater than 1."
            << exit(FatalIOError);
    }
    if (subStepping_)
    {
        IOWarningIn("craftsFlow::craftsFlow", settingsDict_)
            << "subStepping is enabled. Some boundary conditions look for "
            << "'U' or 'phi' by default.  Direct flow-related boundary "
            << "conditions to 'USub' or 'phiSub', and reaction-related "
            << "boundary conditions to 'UGlobal' or 'phiGlobal' instead by "
            << "adding 'U' or 'phi' keywords in their dictionaries."
            << endl;
    }
    else
    {
        IOWarningIn("craftsFlow::craftsFlow", settingsDict_)
            << "subStepping is disabled. Some boundary conditions look for "
            << "'U' or 'phi' by default.  Direct them all to 'UGlobal' or "
            << "'phiGlobal' instead by adding 'U' or 'phi' keywords in their "
            << "dictionaries."
            << endl;
    }
}

// * * * * * * * * * * * * * * Member Functions  * * * * * * * * * * * * * * //

Foam::label Foam::craftsFlow::setNSubSteps(label newNSubSteps) const
{
    // Force to even
    if (newNSubSteps % 2)
    {
        newNSubSteps++;
    }
    
    // Apply limits
    nSubSteps_ = max(minNSubSteps_, newNSubSteps);
    nSubSteps_ = min(maxNSubSteps_, nSubSteps_);

    return nSubSteps_;
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
            << endl;
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
    if (subStepping_)
    {
        phiScale_ = max
        (
            model_.calculateRms(fineSaveSpotPhiPtr_()), phiMinScale_
        );
        pScale_ = max(model_.calculateRms(fineSaveSpotPPtr_()), pMinScale_);
    }
    else
    {
        phiScale_ = max(model_.calculateRms(phiActive_), phiMinScale_);
        pScale_ = max(model_.calculateRms(pActive_), pMinScale_);
    }
    if (outputFlowErrorScales_)
    {
        Info << "flowErrorScales: U = " << phiScale_ << ", p = " << pScale_
            << endl;
    }
}


Foam::scalar Foam::craftsFlow::testConvergence() const
{
    scalar returnMe(-VGREAT);
    scalar phiTest(-VGREAT);
    scalar pTest(-VGREAT);
    if (subStepping_)
    {
//Info << "pCoarse = [" << coarseSaveSpotP_ << "]" << endl;
//Info << "pFine = [" << fineSaveSpotPPtr_() << "]" << endl;
        phiTest = model_.calculateRmsError
        (
            fineSaveSpotPhiPtr_(),
            coarseSaveSpotPhi_,
            mesh_.magSf().field(),
            phiScale_
        ) - phiConvergence_;
        pTest = model_.calculateRmsError
        (
            fineSaveSpotPPtr_(),
            coarseSaveSpotP_,
            mesh_.V().field(),
            pScale_
        ) - pConvergence_;
    }
    else
    {
        phiTest = model_.calculateRmsError
        (
            phiActive_.internalField(),
            coarseSaveSpotPhi_,
            mesh_.magSf().field(),
            phiScale_
        ) - phiConvergence_;
        pTest = model_.calculateRmsError
        (
            pActive_.internalField(),
            coarseSaveSpotP_,
            mesh_.V().field(),
            pScale_
        ) - pConvergence_;
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
    
    return returnMe;
}


Foam::scalar Foam::craftsFlow::calculateNextBestDeltaT() const
{
    // Don't have to worry about subStepping here
    scalar dtPhi
    (
        admReactionReader::calculateNextDeltaT
        (
            coarseSaveSpotPhi_,
            phiActive_.internalField(),
            phiConvergence_,
            mesh_.magSf().field(),
            runTime_.deltaT().value(),
            phiScale_
        )
    );

    scalar dtP
    (
        admReactionReader::calculateNextDeltaT
        (
            coarseSaveSpotP_,
            pActive_.internalField(),
            pConvergence_,
            mesh_.V().field(),
            runTime_.deltaT().value(),
            pScale_
        )
    );
    
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
