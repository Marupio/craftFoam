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

#include "craftsSh2GasUdfs.H"
#include "addToRunTimeSelectionTable.H"
#include "craftsModel.H"

// * * * * * * * * * * * * * Private Member Functions  * * * * * * * * * * * //

template<int matrixSize>
Foam::label Foam::craftsSh2GasUdfs<matrixSize>::getShp()
{
    // Gather all variables
    // ion concentrations must be implicit variables
    admImplicitVariable& ivarShp(admVars.lookupImplicit("S_h_p"));

    // other concentrations - it does not matter what type of variable
    const admVariable& varSvam(admVars.lookup("S_va_m"));
    const admVariable& varSbum(admVars.lookup("S_bu_m"));
    const admVariable& varSprom(admVars.lookup("S_pro_m"));
    const admVariable& varSacm(admVars.lookup("S_ac_m"));
    const admVariable& varShco3m(admVars.lookup("S_hco3_m"));
    const admVariable& varSnh4p(admVars.lookup("S_nh4_p"));
    const admVariable& varScat(admVars.lookup("S_cat"));
    const admVariable& varSan(admVars.lookup("S_an"));
    
    // Field references
    // direct access for ion concentrations: varName().internalField()
    scalarField& S_h_p(ivarShp().internalField());
    
    // const ref only for other variables
    const scalarField& S_va_m(varSvam.evaluateField());
    const scalarField& S_bu_m(varSbum.evaluateField());
    const scalarField& S_pro_m(varSprom.evaluateField());
    const scalarField& S_ac_m(varSacm.evaluateField());
    const scalarField& S_hco3_m(varShco3m.evaluateField());
    const scalarField& S_nh4_p(varSnh4p.evaluateField());
    const scalarField& S_cat(varScat.evaluateField());
    const scalarField& S_an(varSan.evaluateField());

    // Gather coefficients
    const admCoefficient& coeffKgmac(admCoeffs("kg_m_ac"));
    const admCoefficient& coeffKgmbu(admCoeffs("kg_m_bu"));
    const admCoefficient& coeffKgmpro(admCoeffs("kg_m_pro"));
    const admCoefficient& coeffKgmva(admCoeffs("kg_m_va"));
    const admCoefficient& coeffKw(admCoeffs("k_w"));
    
    // Coefficient values
    // These ones we assume are uniform, therefore a single element will do
    const scalar kg_m_ac(coeffKgmac.evaluate(0));
    const scalar kg_m_bu(coeffKgmbu.evaluate(0));
    const scalar kg_m_pro(coeffKgmpro.evaluate(0));
    const scalar kg_m_va(coeffKgmva.evaluate(0));
    // These coefficients are non-uniform, therefore a full field is required
    const scalarField k_w(coeffKw.evaluateField());

    // Quadratic equation: Ax^2 + Bx + C = 0
    // A is unity
    scalarField B
    (
        S_cat + S_nh4_p - S_hco3_m - S_an
      - S_ac_m / kg_m_ac
      - S_pro_m / kg_m_pro
      - S_bu_m / kg_m_bu
      - S_va_m / kg_m_va
    );
    // C is Kw

    scalarField temp
    (
        sqr(B) - 4 * k_w
    );

    if (min(temp) < 0)
    {
        WarningIn("craftsSh2GasUdfs::getShp")
        << "Imaginary hydrogen ion concentration detected." << endl;

        // Tell solver to repeat using a smaller timestep            
        return -1;
    }    
    
    S_h_p = (-B + sqrt(sqr(B) - 4 * k_w)) / 2;

    //Stabilise solution - apply limits to S_h_p
    admVars.applyLimits(ivarShp);

    return 0;
}


// * * * * * * * * * * * * * * * * Constructors  * * * * * * * * * * * * * * //

template<int matrixSize>
Foam::craftsSh2GasUdfs<matrixSize>::craftsSh2GasUdfs
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
    )
{}

// * * * * * * * * * * * * * * Member Functions  * * * * * * * * * * * * * * //

template<int matrixSize>
void Foam::craftsSh2GasUdfs<matrixSize>::initializeTimestep()
{
    admGas_.initializeTimestep();
}


template<int matrixSize>
void Foam::craftsSh2GasUdfs<matrixSize>::applyVariableLimits()
{
    admVars.massConservingApplyStandardLimits();
    admVars.massConservingApplyImplicitLimits();
}


template<int matrixSize>
void Foam::craftsSh2GasUdfs<matrixSize>::initializeImplicitLoop()
{
    admGas_.initializeImplicitLoop();
}


template<int matrixSize>
Foam::label Foam::craftsSh2GasUdfs<matrixSize>::implicitLoop()
{
    label returnMe(0);

    // Quadratic equation solver for hydrogen ion concentration
    label ionsError(getShp());
    switch (ionsError)
    {
        default:
            // Unknown value
            FatalErrorIn("craftsSh2GasUdfs::implicitLoop")
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
        case 0:
            // success, keep going
            break;
        case -1:
            // total failure - return to solver level
            return ionsError;
        case -2:
            // convergence failure - repeat timestep
            return ionsError;
        case -3:
            // stabilization failure - repeat implicit loop
            returnMe = -3;
            break;
        case -4:
            // funtion hook internal failure - handled same as -3
            returnMe = -3;
            break;
    }

    // Set milestone, force lazy evaluating objects to recalculate
    model.setMilestone();

    // Perform gas model calculations    
    label gasError(admGas_.implicitLoop());
    switch (gasError)
    {
        default:
            // Unknown value
            FatalErrorIn("craftsSh2GasUdfs::implicitLoop")
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
        case 0:
            // success, keep going
            break;
        case -1:
            // total failure - return to solver level
            return gasError;
        case -2:
            // convergence failure - repeat timestep
            return gasError;
        case -3:
            // stabilization failure - repeat implicit loop
            returnMe = -3;
            break;
        case -4:
            // funtion hook internal failure - handled same as -3
            returnMe = -3;
            break;
    }
    
    // Set milestone, force lazy evaluating objects to recalculate
    model.setMilestone();
    
    // Run Newton-Raphson solver on implicit variables marked for autoSolve
    // This solves S_h2
    label autoSolveError(model.implicitAutoSolveAllField());
    switch (autoSolveError)
    {
        default:
            // Unknown value
            FatalErrorIn("craftsSh2GasUdfs::implicitLoop")
                << "Received error level " << autoSolveError << " from "
                << "implicitAutoSolveAllField. Valid error levels are:"
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
        case 0:
            // success, keep going
            break;
        case -1:
            // total failure - return to solver level
            return autoSolveError;
        case -2:
            // convergence failure - repeat timestep
            return autoSolveError;
        case -3:
            // stabilization failure - repeat implicit loop
            returnMe = -3;
            break;
        case -4:
            // funtion hook internal failure - handled same as -3
            returnMe = -3;
            break;
    }

    // Set milestone, force lazy evaluating objects to recalculate
    model.setMilestone();

    return returnMe;
}


template<int matrixSize>
void Foam::craftsSh2GasUdfs<matrixSize>::finalizeImplicitLoop()
{
    admGas_.finalizeImplicitLoop();
}


template<int matrixSize>
void Foam::craftsSh2GasUdfs<matrixSize>::finalizeTimestep()
{
    admGas_.finalizeTimestep();
}


template<int matrixSize>
void Foam::craftsSh2GasUdfs<matrixSize>::saveState(const label slot)
{
    admGas_.saveState(slot);
}


template<int matrixSize>
void Foam::craftsSh2GasUdfs<matrixSize>::clearState(const label slot)
{
    admGas_.clearState(slot);
}


template<int matrixSize>
void Foam::craftsSh2GasUdfs<matrixSize>::loadState(const label slot)
{
    admGas_.loadState(slot);
}


template<int matrixSize>
Foam::label Foam::craftsSh2GasUdfs<matrixSize>::nStates() const
{
    return admGas_.nStates();
}


template<int matrixSize>
bool Foam::craftsSh2GasUdfs<matrixSize>::validState(const label slot)
    const
{
    return admGas_.validState(slot);
}

// ************************************************************************* //
