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

#include "craftsPimpleFlow.H"
#include "addToRunTimeSelectionTable.H"

// * * * * * * * * * * * * * * Static Data Members * * * * * * * * * * * * * //

namespace Foam
{
    defineTypeNameAndDebug(craftsPimpleFlow, 0);

    addToRunTimeSelectionTable
    (
        craftsFlow,
        craftsPimpleFlow,
        nameConstructor
    );
}


// * * * * * * * * * * * * * Private Member Functions  * * * * * * * * * * * //

void Foam::craftsPimpleFlow::solveUEquation()
{
    // Shorthand
    volVectorField& U(UActive_);
    surfaceScalarField& phi(phiActive_);
    volScalarField& p(pActive_);

    // Solve the momentum equation

    volVectorField::DimensionedInternalField sourceTerm
    (
        IOobject
        (
            "sourceTerm",
            runTime_.timeName(),
            mesh_,
            IOobject::NO_READ,
            IOobject::NO_WRITE
        ),
        mesh_,
        dimLength / dimTime / dimTime
    );
    sourceTerm.replace(0, sourceTermX_.evaluateDimensionedField());
    sourceTerm.replace(1, sourceTermY_.evaluateDimensionedField());
    sourceTerm.replace(2, sourceTermZ_.evaluateDimensionedField());

    UEqnPtr_.clear();
    UEqnPtr_.set
    (
        new fvVectorMatrix
        (
            fvm::ddt(U)
          + fvm::div(phi, U)
          + turbulence_->divDevReff(U)
          + sourceTerm
        )
    );
    fvVectorMatrix& UEqn(UEqnPtr_());

    if (oCorr_ == nOuterCorr_ - 1)
    {
        UEqn.relax(1);
    }
    else
    {
        UEqn.relax();
    }

    rUAPtr_.clear();
    rUAPtr_.set
    (
        new volScalarField
        (
            "rUA",
            1.0/UEqn.A()
        )
    );
    volScalarField& rUA(rUAPtr_());

    if (momentumPredictor_)
    {
        (UEqn == -fvc::grad(p))().solve(mesh_.solver("U"));

        if (oCorr_ == nOuterCorr_ - 1)
        {
            solve(UEqn == -fvc::grad(p), mesh_.solver("UFinal"));
        }
        else
        {
            solve(UEqn == -fvc::grad(p), mesh_.solver("U"));
        }
    }
    else
    {
        U = rUA * (UEqn.H() - fvc::grad(p));
        U.correctBoundaryConditions();
    }
}


void Foam::craftsPimpleFlow::solvePEquation(int corr)
{
    // Shorthand
    volVectorField& U(UActive_);
    surfaceScalarField& phi(phiActive_);
    volScalarField& p(pActive_);
    fvVectorMatrix& UEqn(UEqnPtr_());
    volScalarField& rUA(rUAPtr_());
    
    U = rUA * UEqn.H();

    if (nCorr_ <= 1)
    {
        UEqnPtr_.clear();
    }

    phi = (fvc::interpolate(U) & mesh_.Sf())
        + fvc::ddtPhiCorr(rUA, U, phi);

    adjustPhi(phi, U, p);

    // Non-orthogonal pressure corrector loop
    for (int nonOrth=0; nonOrth<=nNonOrthCorr_; nonOrth++)
    {
        // Pressure corrector
        fvScalarMatrix pEqn
        (
            fvm::laplacian(rUA, p) == fvc::div(phi)
        );

        pEqn.setReference(pRefCell_, pRefValue_);

        if
        (
            oCorr_ == nOuterCorr_ - 1
         && corr == nCorr_ - 1
         && nonOrth == nNonOrthCorr_
        )
        {
            pEqn.solve(mesh_.solver("pFinal"));
        }
        else
        {
            pEqn.solve(mesh_.solver("p"));
        }

        if (nonOrth == nNonOrthCorr_)
        {
            phi -= pEqn.flux();
        }
    }

    continuityErrs();

    // Explicitly relax pressure for momentum corrector except for last
    // corrector
    if (oCorr_ != nOuterCorr_ - 1)
    {
        p.relax();
    }

    U -= rUA*fvc::grad(p);
    U.correctBoundaryConditions();
}


void Foam::craftsPimpleFlow::continuityErrs()
{
    // Shorthand
    surfaceScalarField& phi(phiActive_);

    volScalarField contErr = fvc::div(phi);

    sumLocalContErr = runTime_.deltaT().value()*
        mag(contErr)().weightedAverage(mesh_.V()).value();

    globalContErr = runTime_.deltaT().value()*
        contErr.weightedAverage(mesh_.V()).value();

    cumulativeContErr += globalContErr;

    if (outputFlowContinuityErrors_)
    {
        Info << "flowContinuityErrors: sum local = " << sumLocalContErr
            << ", global = " << globalContErr
            << ", cumulative = " << cumulativeContErr
            << endl;
    }
}


// * * * * * * * * * * * * * * * * Constructors  * * * * * * * * * * * * * * //

Foam::craftsPimpleFlow::craftsPimpleFlow
(
    volVectorField& U,
    surfaceScalarField& phi,
    volScalarField& p,
    admReactionReader& model
)
:
    craftsFlow(U, phi, p, model),
    
    transport_(UActive_, phiActive_),
    turbulence_
    (
        incompressible::turbulenceModel::New
        (
            UActive_, phiActive_, transport_
        ).ptr()
    ),

    UEqnPtr_(0),
    rUAPtr_(0),

    pimpleDict_(mesh_.solutionDict().subDict("PIMPLE")),
    nCorr_
    (
        readInt(pimpleDict_.lookup("nCorrectors"))
    ),
    nNonOrthCorr_
    (
        pimpleDict_.lookupOrDefault<int>("nNonOrthogonalCorrectors", 0)
    ),
    oCorr_(0),
    nOuterCorr_
    (
        pimpleDict_.lookupOrDefault<int>("nOuterCorrectors", 0)
    ),
    momentumPredictor_
    (
        pimpleDict_.lookupOrDefault<Switch>("momentumPredictor", true)
    ),
    pRefCell_(0),
    pRefValue_(0.0),

    sumLocalContErr(0.0),
    globalContErr(0.0),
    cumulativeContErr(0.0),
    
    sourceTermXName_
    (
        transport_.found("sourceTermX")
      ? transport_.lookup("sourceTermX")
      : word::null
    ),
    sourceTermYName_
    (
        transport_.found("sourceTermY")
      ? transport_.lookup("sourceTermY")
      : word::null
    ),
    sourceTermZName_
    (
        transport_.found("sourceTermZ")
      ? transport_.lookup("sourceTermZ")
      : word::null
    ),
    
    zeroSource_
    (
        model.admVars(),
        "flowModelZeroSource",
        dimLength / dimTime / dimTime
    ),
        
    sourceTermX_
    (
        sourceTermXName_ != word::null
      ? model.admCoeffs()(sourceTermXName_)
      : zeroSource_
    ),
    sourceTermY_
    (
        sourceTermYName_ != word::null
      ? model.admCoeffs()(sourceTermYName_)
      : zeroSource_
    ),
    sourceTermZ_
    (
        sourceTermZName_ != word::null
      ? model.admCoeffs()(sourceTermZName_)
      : zeroSource_
    ),

    coarseSaveSpotR_
    (
        mesh_.nCells(), pTraits<symmTensor>::zero
    ),
    fineSaveSpotRPtr_
    (
        subStepping_
      ? new symmTensorField(mesh_.nCells(), pTraits<symmTensor>::zero)
      : 0
    ),

    RMinScale_
    (
        SMALL, SMALL, SMALL,
               SMALL, SMALL,
                      SMALL
    ),
    RScale_
    (
        SMALL, SMALL, SMALL,
               SMALL, SMALL,
                      SMALL
    ),
    RConvergence_
    (
        symmTensor
        (
            settingsDict_.subDict("adaptiveTimeSteppingTolerance").lookup("R")
        )
    )
{
    // Override lduMatrix::debug
    if (outputFlowSolverPerformance_)
    {
        lduMatrix::debug = 1;
    }
    else
    {
        lduMatrix::debug = 0;
    }

    // Check convergence and error scaling
    if (phiConvergence_ < 0)
    {
        FatalErrorIn("craftsPimpleFlow::craftsPimpleFlow")
            << "Positive value expected for phi adaptive timestepping "
            << "tolerance. Set in admSettingsDict, flowModel subDictionary, "
            << "adaptiveTimesteppingTolerance subdictionary, keyword phi"
            << abort(FatalError);
    }
    if (settingsDict_.found("minErrorScale"))
    {
        if (settingsDict_.subDict("minErrorScale").found("R"))
        {
            RMinScale_ = symmTensor
            (
                settingsDict_.subDict("minErrorScale").lookup("R")
            );
        }
    }

    // Reference pressure
    setRefCell
    (
        pGlobal_,
        mesh_.solutionDict().subDict("PIMPLE"),
        pRefCell_,
        pRefValue_
    );

    // Dimension checking source term
    if
    (
        dimensionSet::debug
     && (
            (sourceTermX_.dimensions() != sourceTermY_.dimensions())
         || (sourceTermX_.dimensions() != sourceTermZ_.dimensions())
        )
    )
    {
        WarningIn("craftsPimpleFlow::craftsPimpleFlow")
            << "Source term dimensions do not match for all directions"
            << endl;
    }
    sourceTermX_.dimensions() = sourceTermY_.dimensions();
    sourceTermX_.dimensions() = sourceTermZ_.dimensions();
}


// * * * * * * * * * * * * * * Member Functions  * * * * * * * * * * * * * * //

Foam::label Foam::craftsPimpleFlow::step()
{
    volScalarField& p(pActive_);

    // --- Pressure-velocity PIMPLE corrector loop
    for (oCorr_ = 0; oCorr_ < nOuterCorr_; oCorr_++)
    {
        if (nOuterCorr_ != 1)
        {
            p.storePrevIter();
        }
        solveUEquation();

        // --- PIMPLE loop
        for (int corr=0; corr<nCorr_; corr++)
        {
            solvePEquation(corr);
        }

        turbulence_->correct();
    }
    
    return 0;
}


void Foam::craftsPimpleFlow::storeCoarseSolution()
{
    coarseSaveSpotR_ = turbulence_->R()().internalField();
    this->craftsFlow::storeCoarseSolution();
}


void Foam::craftsPimpleFlow::storeFineSolution()
{
    fineSaveSpotRPtr_() = turbulence_->R()().internalField();
    this->craftsFlow::storeFineSolution();
}


void Foam::craftsPimpleFlow::calculateScales()
{
    symmTensorField& R(coarseSaveSpotR_);
    RScale_.component(0) = max
    (
        model_.calculateRms(R.component(0)),
        RMinScale_.component(0)
    );
    RScale_.component(1) = max
    (
        model_.calculateRms(R.component(1)),
        RMinScale_.component(1)
    );
    RScale_.component(2) = max
    (
        model_.calculateRms(R.component(2)),
        RMinScale_.component(2)
    );
    RScale_.component(3) = max
    (
        model_.calculateRms(R.component(3)),
        RMinScale_.component(3)
    );
    RScale_.component(4) = max
    (
        model_.calculateRms(R.component(4)),
        RMinScale_.component(4)
    );
    RScale_.component(5) = max
    (
        model_.calculateRms(R.component(5)),
        RMinScale_.component(5)
    );
    this->craftsFlow::calculateScales();
    if (outputFlowErrorScales_)
    {
        Info << "flowErrorScales: R = " << RScale_ << endl;
    }
}


Foam::scalar Foam::craftsPimpleFlow::testConvergence() const
{
    // Turbulence residuals
    scalar turbError(-VGREAT);

    if (subStepping_)
    {
        scalarField Rerr(pTraits<symmTensor>::nComponents, -VGREAT);
        forAll(Rerr, comp)
        {
            Rerr[comp] = admReactionReader::calculateRmsError
            (
                fineSaveSpotRPtr_().component(comp),
                coarseSaveSpotR_.component(comp),
                mesh_.V().field(),
                RScale_.component(comp)
            ) - RConvergence_.component(comp);
            turbError = max(turbError, Rerr[comp]);
        }
    }
    else
    {
        const symmTensorField R(turbulence_->R()().internalField());
        scalarField Rerr(pTraits<symmTensor>::nComponents, -VGREAT);
        forAll(Rerr, comp)
        {
            Rerr[comp] = admReactionReader::calculateRmsError
            (
                R.component(comp),
                coarseSaveSpotR_.component(comp),
                mesh_.V().field(),
                RScale_.component(comp)
            ) - RConvergence_.component(comp);
            turbError = max(turbError, Rerr[comp]);
        }
    }

    // p and phi residuals
    scalar pPhiError(this->craftsFlow::testConvergence());
    scalar returnMe(max(turbError, pPhiError));

    if (outputFlowResiduals_)
    {
        Info << "flowResiduals: Turbulence, R: " << turbError;
        if (turbError > 0)
        {
            Info << " (fail)." << endl;
        }
        else
        {
            Info << " (pass)." << endl;
        }
    }

    return returnMe;
}


Foam::scalar Foam::craftsPimpleFlow::calculateNextBestDeltaT() const
{
    const symmTensorField R(turbulence_->R()().internalField());

    scalarField dtR(pTraits<symmTensor>::nComponents, VGREAT);
    scalar minDtR(VGREAT);
    forAll(dtR, comp)
    {
        dtR[comp] = admReactionReader::calculateNextDeltaT
        (
            coarseSaveSpotR_.component(comp),
            R.component(comp),
            RConvergence_.component(comp),
            mesh_.V().field(),
            runTime_.deltaT().value(),
            RScale_.component(comp)
        );
        minDtR = min(minDtR, dtR[comp]);
    }
    
    scalar dtPhiP(this->craftsFlow::calculateNextBestDeltaT());
    
    if (outputFlowTimestepEstimate_)
    {
        Info << "flowTimestepEstimate: Turbulence, R dtEst = "
            << minDtR;
        if (dtPhiP < minDtR)
        {
            Info << " (not a bottleneck)" << endl;
        }
        else
        {
            Info << " (bottleneck)" << endl;
        }
    }
    
    return min(dtPhiP,minDtR);
}


void Foam::craftsPimpleFlow::saveState(const label slot)
{
    turbulence_->saveState(slot);
    this->craftsFlow::saveState(slot);
}
            

void Foam::craftsPimpleFlow::clearState(const label slot)
{
    // do nothing
}

            
void Foam::craftsPimpleFlow::loadState(const label slot)
{
    turbulence_->loadState(slot);
    this->craftsFlow::loadState(slot);
}
            

Foam::label Foam::craftsPimpleFlow::nStates() const
{
    return -1;
}


bool Foam::craftsPimpleFlow::validState(const label slot) const
{
    return true;
}

// ************************************************************************* //
