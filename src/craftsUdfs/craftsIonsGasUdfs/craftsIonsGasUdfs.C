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

#include "craftsIonsGasUdfs.H"
#include "addToRunTimeSelectionTable.H"
#include "craftsModel.H"

// * * * * * * * * * * * * * Private Member Functions  * * * * * * * * * * * //

template<int matrixSize>
Foam::label Foam::craftsIonsGasUdfs<matrixSize>::getIons()
{
    // Gather all variables
    // ion concentrations must be implicit variables
    admImplicitVariable& ivarShp(admVars.lookupImplicit("S_h_p"));
    admImplicitVariable& ivarSvam(admVars.lookupImplicit("S_va_m"));
    admImplicitVariable& ivarSbum(admVars.lookupImplicit("S_bu_m"));
    admImplicitVariable& ivarSprom(admVars.lookupImplicit("S_pro_m"));
    admImplicitVariable& ivarSacm(admVars.lookupImplicit("S_ac_m"));
    admImplicitVariable& ivarShco3m(admVars.lookupImplicit("S_hco3_m"));
    admImplicitVariable& ivarSnh4p(admVars.lookupImplicit("S_nh4_p"));
    // other concentrations - it does not matter what type of variable
    const admVariable& varSva(admVars.lookup("S_va"));
    const admVariable& varSbu(admVars.lookup("S_bu"));
    const admVariable& varSpro(admVars.lookup("S_pro"));
    const admVariable& varSac(admVars.lookup("S_ac"));
    const admVariable& varSic(admVars.lookup("S_ic"));
    const admVariable& varSin(admVars.lookup("S_in"));
    const admVariable& varScat(admVars.lookup("S_cat"));
    const admVariable& varSan(admVars.lookup("S_an"));

    // Field references
    // direct access for ion concentrations: varName().internalField()
    scalarField& S_h_p(ivarShp().internalField());
    scalarField& S_va_m(ivarSvam().internalField());
    scalarField& S_bu_m(ivarSbum().internalField());
    scalarField& S_pro_m(ivarSprom().internalField());
    scalarField& S_ac_m(ivarSacm().internalField());
    scalarField& S_hco3_m(ivarShco3m().internalField());
    scalarField& S_nh4_p(ivarSnh4p().internalField());

    // const ref only for other variables
    const scalarField& S_va(varSva.evaluateField());
    const scalarField& S_bu(varSbu.evaluateField());
    const scalarField& S_pro(varSpro.evaluateField());
    const scalarField& S_ac(varSac.evaluateField());
    const scalarField& S_ic(varSic.evaluateField());
    const scalarField& S_in(varSin.evaluateField());
    const scalarField& S_cat(varScat.evaluateField());
    const scalarField& S_an(varSan.evaluateField());

    // Gather coefficients
    const admCoefficient& coeffKava(admCoeffs("k_a_va"));
    const admCoefficient& coeffKabu(admCoeffs("k_a_bu"));
    const admCoefficient& coeffKapro(admCoeffs("k_a_pro"));
    const admCoefficient& coeffKaac(admCoeffs("k_a_ac"));
    const admCoefficient& coeffKgmac(admCoeffs("kg_m_ac"));
    const admCoefficient& coeffKgmbu(admCoeffs("kg_m_bu"));
    const admCoefficient& coeffKgmpro(admCoeffs("kg_m_pro"));
    const admCoefficient& coeffKgmva(admCoeffs("kg_m_va"));
    const admCoefficient& coeffKaco2(admCoeffs("k_a_co2"));
    const admCoefficient& coeffKain(admCoeffs("k_a_in"));
    const admCoefficient& coeffKw(admCoeffs("k_w"));
    
    // Coefficient values
    // These ones we assume are uniform, therefore a single element will do
    const scalar k_a_va(coeffKava.evaluate(0));
    const scalar k_a_bu(coeffKabu.evaluate(0));
    const scalar k_a_pro(coeffKapro.evaluate(0));
    const scalar k_a_ac(coeffKaac.evaluate(0));
    const scalar kg_m_ac(coeffKgmac.evaluate(0));
    const scalar kg_m_bu(coeffKgmbu.evaluate(0));
    const scalar kg_m_pro(coeffKgmpro.evaluate(0));
    const scalar kg_m_va(coeffKgmva.evaluate(0));
    // These coefficients are non-uniform, therefore a full field is required
    const scalarField k_a_co2(coeffKaco2.evaluateField());
    const scalarField k_a_in(coeffKain.evaluateField());
    const scalarField k_w(coeffKw.evaluateField());

    label maxIter(ivarShp.autoSolveMaxIter());
    scalar convergence(ivarShp.autoSolveConvergence());
    scalar convergedTo(VGREAT);
    scalar scaleFactor(model.implicitScale()[ivarShp.localIndex()]);
    label iter(0);
    scalarField numerator(mesh.nCells(), 0.0);
    scalarField denominator(mesh.nCells(), 1.0);

    // Normally we'd start the loop here, but the first iteration is pulled out
    // to test for unacceptable loss of accuracy
    S_va_m = k_a_va * S_va / (k_a_va + S_h_p);
    S_bu_m = k_a_bu * S_bu / (k_a_bu + S_h_p);
    S_pro_m = k_a_pro * S_pro / (k_a_pro + S_h_p);
    S_ac_m = k_a_ac * S_ac / (k_a_ac + S_h_p);
    S_hco3_m = k_a_co2 * S_ic / (k_a_co2 + S_h_p);
    S_nh4_p = S_h_p * S_in / (k_a_in + S_h_p);

    denominator =
        1.0 + k_a_in * S_in / pow((k_a_in + S_h_p), 2.0)
        + k_a_co2 * S_ic / pow((k_a_co2 + S_h_p), 2.0)
        + k_a_ac * S_ac / kg_m_ac / pow((k_a_ac + S_h_p), 2.0)
        + k_a_pro * S_pro / kg_m_pro / pow((k_a_pro + S_h_p), 2.0)
        + k_a_bu * S_bu / kg_m_bu / pow((k_a_bu + S_h_p), 2.0)
        + k_a_va * S_va / kg_m_va / pow((k_a_va + S_h_p), 2.0)
        + k_w / pow(S_h_p, 2.0);
    
    scalarField epsilon(mag(denominator * convergence));
    scalarField test
    (
        S_cat + S_nh4_p + S_h_p + S_hco3_m + S_ac_m / kg_m_ac
      + S_pro_m / kg_m_pro + S_bu_m / kg_m_bu + S_va_m / kg_m_va
      + k_w / S_h_p - S_an
    );
    forAll(test, cellIndex)
    {
        if ((test[cellIndex] + epsilon[cellIndex]) == test[cellIndex])
        {
            WarningIn("craftsIonsGasUdfs::getIons")
            << "Machine epsilon is not small enough to achieve the requested "
            << "accuracy at cell index " << cellIndex << endl;

            // Add iterations to performance stats
            model.lastStepIterations() += 1;

            // Tell solver to repeat using a smaller timestep            
            return -1;
        }
    }

    numerator =
        S_cat + S_nh4_p + S_h_p - S_hco3_m - S_ac_m / kg_m_ac
      - S_pro_m / kg_m_pro - S_bu_m / kg_m_bu - S_va_m / kg_m_va
      - k_w / S_h_p - S_an;

    S_h_p = S_h_p - numerator / denominator;

    // Now the loop begins, starting at the second iteration
    iter = 1;
    for (label i(1); i < maxIter; i++)
    {
        iter++;
        S_va_m = k_a_va * S_va / (k_a_va + S_h_p);
        S_bu_m = k_a_bu * S_bu / (k_a_bu + S_h_p);
        S_pro_m = k_a_pro * S_pro / (k_a_pro + S_h_p);
        S_ac_m = k_a_ac * S_ac / (k_a_ac + S_h_p);
        S_hco3_m = k_a_co2 * S_ic / (k_a_co2 + S_h_p);
        S_nh4_p = S_h_p * S_in / (k_a_in + S_h_p);

        numerator =
            S_cat + S_nh4_p + S_h_p - S_hco3_m - S_ac_m / kg_m_ac
          - S_pro_m / kg_m_pro - S_bu_m / kg_m_bu - S_va_m / kg_m_va
          - k_w / S_h_p - S_an;

        denominator =
            1.0 + k_a_in * S_in / pow((k_a_in + S_h_p), 2.0)
            + k_a_co2 * S_ic / pow((k_a_co2 + S_h_p), 2.0)
            + k_a_ac * S_ac / kg_m_ac / pow((k_a_ac + S_h_p), 2.0)
            + k_a_pro * S_pro / kg_m_pro / pow((k_a_pro + S_h_p), 2.0)
            + k_a_bu * S_bu / kg_m_bu / pow((k_a_bu + S_h_p), 2.0)
            + k_a_va * S_va / kg_m_va / pow((k_a_va + S_h_p), 2.0)
            + k_w / pow(S_h_p, 2.0);

        // Incorporate denominator, reusing memory for numerator
        numerator /= denominator;
        S_h_p = S_h_p - numerator;
        convergedTo = (model.calculateRms(numerator) / scaleFactor)
            - convergence;
        if (convergedTo <= 0)
        {
            iter--;
            break;
        }
    }
    
    //Stabilise solution - apply limits to S_h_p
    admVars.applyLimits(ivarShp);

    // Add iterations to performance stats
    model.lastStepIterations() += iter;

    if (outputIonSolverPerformance_)
    {
        Info << "ionSolverPerformance: res = " << convergedTo
            << ", n = " << iter << endl;
    }

    if (iter == maxIter)
    {
        WarningIn("craftsIonsGasUdfs::getIons")
            << "Exceeded maximum iterations" << endl;

        // Tell solver to repeat using a smaller timestep            
        return -1;
    }

    return 0;
}


// * * * * * * * * * * * * * * * * Constructors  * * * * * * * * * * * * * * //

template<int matrixSize>
Foam::craftsIonsGasUdfs<matrixSize>::craftsIonsGasUdfs
(
    craftsModel<matrixSize>& model,
    const word& hooksDictName
)
:
    craftsUdfs<matrixSize>(model, hooksDictName),
    admGas_
    (
        admCoeffs,
        admVars,
        model
    ),
    outputIonSolverPerformance_
    (
        model.outputFlagsDict().found("ionSolverPerformance")
      ? bool(Switch(model.outputFlagsDict().lookup("ionSolverPerformance")))
      : false
    ),
    outputFunctionHooksSummary_
    (
        model.outputFlagsDict().found("functionHooksSummary")
      ? bool(Switch(model.outputFlagsDict().lookup("functionHooksSummary")))
      : false
    )
{}

// * * * * * * * * * * * * * * Member Functions  * * * * * * * * * * * * * * //

template<int matrixSize>
void Foam::craftsIonsGasUdfs<matrixSize>::initializeTimestep()
{
    admGas_.initializeTimestep();
}


template<int matrixSize>
void Foam::craftsIonsGasUdfs<matrixSize>::applyVariableLimits()
{
    admVars.massConservingApplyStandardLimits();
    admVars.massConservingApplyImplicitLimits();
}


template<int matrixSize>
void Foam::craftsIonsGasUdfs<matrixSize>::initializeImplicitLoop()
{
    admGas_.initializeImplicitLoop();
}


template<int matrixSize>
Foam::label Foam::craftsIonsGasUdfs<matrixSize>::implicitLoop()
{
    label returnMe(0);

    // Perform gas model calculations
    label gasError(admGas_.implicitLoop());
    switch (gasError)
    {
        default:
            // Unknown value
            FatalErrorIn("craftsIonsGasUdfs::implicitLoop")
                << "Received error level " << gasError << " from "
                << "gas model implicitLoop. Valid error levels are:"
                << token::NL << token::TAB
                << " 0 = success"
                << token::NL << token::TAB
                << "-1 = total failure, reduce delta t"
                << token::NL << token::TAB
                << "-2 = convergence failure, repeat timestep"
                << token::NL << token::TAB
                << "-3 = stabilization failure, repeat implicit loop"
                << token::NL << token::TAB
                << "-4 = function hook internal failure"
                << abort(FatalError);
        case 0: // success, keep going
            break;
        case -1: // total failure - return to solver level
        case -2: // convergence failure - repeat timestep
            if (outputFunctionHooksSummary_)
            {
                Info << "functionHooksSummary:"
                    << " gasModel=" << gasError
                    << " ionModel"
                    << " (fail)" << endl;
            }
            return gasError;
        case -3: // stabilization failure - repeat implicit loop
        case -4: // funtion hook internal failure - handled same as -3
            returnMe = -3;
            break;
    }
    
    // Set milestone, force lazy evaluating objects to recalculate
    model.setMilestone();
    
    // Newton-Raphson solver for ion concentrations
    label ionsError(getIons());
    switch (ionsError)
    {
        default:
            // Unknown value
            FatalErrorIn("craftsIonsGasUdfs::implicitLoop")
                << "Received error level " << ionsError << " from "
                << "getIons. Valid error levels are:"
                << token::NL << token::TAB
                << " 0 = success"
                << token::NL << token::TAB
                << "-1 = total failure, reduce delta t"
                << token::NL << token::TAB
                << "-2 = convergence failure, repeat timestep"
                << token::NL << token::TAB
                << "-3 = stabilization failure, repeat implicit loop"
                << token::NL << token::TAB
                << "-4 = function hook internal failure"
                << abort(FatalError);
        case 0: // success, keep going
            break;
        case -1: // total failure - return to solver level
        case -2: // convergence failure - repeat timestep
            if (outputFunctionHooksSummary_)
            {
                Info << "functionHooksSummary:"
                    << " gasModel=" << gasError
                    << " ionModel=" << ionsError
                    << " (fail)" << endl;
            }
            return ionsError;
        case -3: // stabilization failure - repeat implicit loop
        case -4: // funtion hook internal failure - handled same as -3
            if (outputFunctionHooksSummary_)
            {
                Info << "functionHooksSummary: gasModel [ionModel="
                    << ionsError << "]" << endl;
            }
            return ionsError;
            returnMe = -3;
            break;
    }

    if (outputFunctionHooksSummary_)
    {
        Info << "functionHooksSummary:"
            << " gasModel=" << gasError
            << " ionModel=" << ionsError;
        if (returnMe == 0)
        {
            Info << " (pass)" << endl;
        }
        else
        {
            Info << " (fail)" << endl;
        }
    }

    // Set milestone, force lazy evaluating objects to recalculate
    model.setMilestone();

    return returnMe;
}


template<int matrixSize>
void Foam::craftsIonsGasUdfs<matrixSize>::finalizeImplicitLoop()
{
    admGas_.finalizeImplicitLoop();
}


template<int matrixSize>
void Foam::craftsIonsGasUdfs<matrixSize>::finalizeTimestep()
{
    admGas_.finalizeTimestep();
}


template<int matrixSize>
void Foam::craftsIonsGasUdfs<matrixSize>::saveState(const label slot)
{
    admGas_.saveState(slot);
}


template<int matrixSize>
void Foam::craftsIonsGasUdfs<matrixSize>::clearState(const label slot)
{
    admGas_.clearState(slot);
}


template<int matrixSize>
void Foam::craftsIonsGasUdfs<matrixSize>::loadState(const label slot)
{
    admGas_.loadState(slot);
}


template<int matrixSize>
Foam::label Foam::craftsIonsGasUdfs<matrixSize>::nStates() const
{
    return admGas_.nStates();
}


template<int matrixSize>
bool Foam::craftsIonsGasUdfs<matrixSize>::validState(const label slot) const
{
    return admGas_.validState(slot);
}

// ************************************************************************* //
