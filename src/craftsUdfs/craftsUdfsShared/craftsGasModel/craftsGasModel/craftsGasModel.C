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

#include "craftsGasModel.H"
#include "ODESolver.H"

// * * * * * * * * * * * * * Private Member Functions  * * * * * * * * * * * //

template<int matrixSize>
void Foam::craftsGasModel<matrixSize>::allocateMemory()
{
    // Fill in lists
    forAll(names_, i)
    {
        const dictionary& speciesDict
        (
            craftsGasModelDict_.subDict("species").subDict(names_[i])
        );
        
        vars_.set
        (
            i,
            & admVars_.lookup(word(speciesDict.lookup("liquidVariable")))
        );
        
        if (vars_[i].varType() == admVariable::vtderived)
        {
            const admDerivedVariable& derVar
            (
                admVars_.derived(vars_[i].localIndex())
            );
            label jFound(0);
            forAll(derVar.jVars(), jIndex)
            {
                const fvPatchScalarField& fvp
                (
                    derVar.evaluateGeometricField().boundaryField()
                    [surfacePatchIndex_]
                );
                if (fvp.type() == "dissolvedGasTransfer")
                {
                    if (jFound == 0)
                    {
                        jFound = jIndex;
                    }
                    else
                    {
                        FatalErrorIn("craftsGasModel::allocateMemory")
                            << "Gas model for " << names_[i] << " involves "
                            << "derived variable " << vars_[i].name()
                            << ", but this variable has more than one non-"
                            << "zero Jacobian variable with a dissolvedGas"
                            << "Transfer boundary condition on patch "
                            << surfacePatchName_ << ".  Cannot handle that."
                            << abort(FatalError);
                    }
                }
            }
            bcVars_.set
            (
                i,
                & admVars_.all(derVar.jVars()[jFound].globalIndex())
            );
        }
        else
        {
            const fvPatchScalarField& fvp
            (
                vars_[i].evaluateGeometricField().boundaryField()
                [surfacePatchIndex_]
            );
            if (fvp.type() != "dissolvedGasTransfer")
            {
                FatalErrorIn("craftsGasModel::allocateMemory")
                    << "Gas model for " << names_[i] << " associated variable "
                    << vars_[i].name() << " requires a dissolvedGasTransfer"
                    << "boundary condition on patch " << surfacePatchName_
                    << "."
                    << abort(FatalError);
            }
            bcVars_.set
            (
                i,
                & vars_[i]
            );
        }
        
        varLiquidStart_.set
        (
            i,
            new scalarField(vars_[i].evaluateField())
        );

        coeffs_n_.set
        (
            i,
            & admCoeffs_.readCoefficient
            (
                speciesDict,
                word("conversionToMoles"),
                word
                (
                    surfacePatchName_ + "(" + names_[i] + "_conversionToMoles)"
                )
            )
        );

        coeffs_kla_.set
        (
            i,
            & admCoeffs_.readCoefficient
            (
                speciesDict,
                word("k_L_a"),
                word
                (
                    surfacePatchName_ + "(" + names_[i] + "_k_L_a)"
                )
            )
        );

        coeffs_kh_.set
        (
            i,
            & admCoeffs_.readCoefficient
            (
                speciesDict,
                word("k_H"),
                word
                (
                    surfacePatchName_ + "(" + names_[i] + "_k_H)"
                )
            )
        );

        coeffs_kb_.set
        (
            i,
            & admCoeffs_.readCoefficient
            (
                speciesDict,
                word("k_b"),
                word
                (
                    surfacePatchName_ + "(" + names_[i] + "_k_b)"
                )
            )
        );

        upperLimits_[i] =
            speciesDict.found("upperLimit")
          ? readScalar(speciesDict.lookup("upperLimit"))
          : VGREAT;

        lowerLimits_[i] =
            speciesDict.found("lowerLimit")
          ? readScalar(speciesDict.lookup("lowerLimit"))
          : scalar(0.0);

        relaxationFactors_[i] =
            speciesDict.found("relaxation")
          ? readScalar(speciesDict.lookup("relaxation"))
          : scalar(1.0);

        pPartial_.set
        (
            i,
            new dimensionedScalar
            (
                "p_" + names_[i],
                dimPressure,
                0.0
            )
        );
        integralRTSgDt_.set
        (
            i,
            new dimensionedScalar
            (
                "interimValue_intRTSgDt_" + names_[i],
                dimensionSet(1, -1, -1, -1, 0, 0, 0)
                    * coeffs_n_[i].dimensions(),
                0.0
            )
        );
        sGas_.set
        (
            i,
            new dimensionedScalar
            (
                "Sg_" + names_[i],
                dimMoles / dimVolume * coeffs_n_[i].dimensions(),
                0.0
            )
        );
        khOld_.set
        (
            i,
            new scalarField
            (
                coeffs_kh_[i].evaluateField()
            )
        );
    }
}


template<int matrixSize>
void Foam::craftsGasModel<matrixSize>::readDict()
{
    // read input dict
    forAll(names_, speciesIndex)
    {
        const dictionary& speciesDict
        (
            craftsGasModelDict_
                .subDict("species")
                .subDict(names_[speciesIndex])
        );
        
        // Mass convergence
        if (speciesDict.found("massConvergence"))
        {
            massConvergence_[speciesIndex] = readScalar
            (
                speciesDict.lookup("massConvergence")
            );
        }
        else
        {
            massConvergence_[speciesIndex] = defaultMassConvergence_;
        }
        if (massConvergence_[speciesIndex] < 0)
        {
            FatalIOErrorIn("craftsGasModel::readDict", craftsGasModelDict_)
                << "Negative massConvergence criteria. Define a non-negative "
                << "massConvergence criteria or default massConvergence "
                << "criterion."
                << exit(FatalIOError);
        }
        if (speciesDict.found("nearZeroMassScale"))
        {
            nearZeroMassScale_[speciesIndex] = readScalar
            (
                speciesDict.lookup("nearZeroMassScale")
            );
        }
        else
        {
            nearZeroMassScale_[speciesIndex] = defaultNearZeroMassScale_;
        }
        
        // Delta convergence
        if (speciesDict.found("deltaConvergence"))
        {
            deltaConvergence_[speciesIndex] = readScalar
            (
                speciesDict.lookup("deltaConvergence")
            );
        }
        else
        {
            deltaConvergence_[speciesIndex] = defaultDeltaConvergence_;
        }
        if (deltaConvergence_[speciesIndex] < 0)
        {
            FatalIOErrorIn("craftsGasModel::readDict", craftsGasModelDict_)
                << "Negative deltaConvergence criteria. Define a non-negative "
                << "deltaConvergence criteria or default deltaConvergence "
                << "criterion."
                << exit(FatalIOError);
        }
        if (speciesDict.found("nearZeroDeltaScale"))
        {
            nearZeroDeltaScale_[speciesIndex] = readScalar
            (
                speciesDict.lookup("nearZeroDeltaScale")
            );
        }
        else
        {
            nearZeroDeltaScale_[speciesIndex] = defaultNearZeroDeltaScale_;
        }
        
        formBubbles_[speciesIndex] = Switch(speciesDict.lookup("formBubbles"));
        /*if (formBubbles_[speciesIndex])
        {
            FatalIOErrorIn("craftsGasModel::readDict", speciesDict)
                << "Bubble model is not currently functional.  Please disable "
                << "this option."
                << exit(FatalIOError);
        }*/
    }

    // Read output dict
    { // scope { } used to clear dsTemp
        dimensionedScalar dsTemp(outputDict_.lookup("T"));
        if (dimensionSet::debug && (dsTemp.dimensions() != T_.dimensions()))
        {
            WarningIn("craftsGasModel::readDict")
                << "T variable in " << outputDictName_ << " dimensions "
                << dsTemp.dimensions() << ", must be equal to variable "
                << var_T_.name() << ", with dimensions "
                << var_T_.dimensions() << endl;
        }
        T_ = dsTemp;
    }
    { // scope { } used to clear dsTemp
        dimensionedScalar dsTemp(outputDict_.lookup("Q"));
        if (dimensionSet::debug && (dsTemp.dimensions() != qGas_.dimensions()))
        {
            WarningIn("craftsGasModel::readDict")
                << "Q variable in " << outputDictName_ << " dimensions "
                << dsTemp.dimensions() << ", must be equal to coefficient "
                << coeff_kp_.name() << " (with dimensions "
                << coeff_kp_.dimensions() << ") * coefficient "
                << coeff_pAtm_.name() << " (with dimensions "
                << coeff_pAtm_.dimensions() << ") * coefficient "
                << coeff_kGasCutOff_.name() << " (with dimensions "
                << coeff_kGasCutOff_.dimensions() << ")" << endl;
        }
        qGas_ = dsTemp;
    }
    { // scope { } used to clear dsTemp
        dimensionedScalar dsTemp(outputDict_.lookup("p_total"));
        if
        (
            dimensionSet::debug
         && (dsTemp.dimensions() != pTotal_.dimensions())
        )
        {
            WarningIn("craftsGasModel::readDict")
                << "p_total variable in " << outputDictName_
                << " dimensions " << dsTemp.dimensions() << ", must be equal "
                << "to coefficient " << coeff_pAtm_.name() << ", with "
                << "dimensions " << coeff_pAtm_.dimensions() << endl;
        }
        pTotal_ = dsTemp;
    }
    { // scope { } used to clear dsTemp
        dimensionedScalar dsTemp(outputDict_.lookup(coeff_pLiq_.name()));
        if
        (
            dimensionSet::debug
         && (dsTemp.dimensions() != pLiq_.dimensions())
        )
        {
            WarningIn("craftsGasModel::readDict")
                << coeff_pLiq_.name() << " variable in " << outputDictName_
                << " dimensions " << dsTemp.dimensions() << ", must be "
                << "equal to its associated coefficient, "
                << coeff_pLiq_.name() << ", with dimensions "
                << coeff_pAtm_.dimensions() << endl;
        }
        pLiq_ = dsTemp;
    }

    // Ignore pAtm as this is a coefficient, and we have its latest value

    // Read partial pressures and vapour pressures
    forAll(names_, speciesIndex)
    {
        word pname("p_" + names_[speciesIndex]);
        word sname("Sg_" + names_[speciesIndex]);
        word iname("interimValue_intRTSgDt_" + names_[speciesIndex]);
        dimensionedScalar dsp(outputDict_.lookup(pname));
        dimensionedScalar dss(outputDict_.lookup(sname));
        if
        (
            dimensionSet::debug
         && (dsp.dimensions() != dimPressure)
        )
        {
            WarningIn("craftsGasModel::readDict")
                << pname << " variable in " << outputDictName_
                << " dimensions " << dsp.dimensions() << ", should be equal "
                << "to pressure " << dimPressure << endl;
        }
        pPartial_[speciesIndex] = dsp;

        if (outputDict_.found(iname))
        {
            dimensionedScalar dsi(outputDict_.lookup(iname));
            if
            (
                dimensionSet::debug
             && (
                    dsi.dimensions()
                 != (
                        dimensionSet(1, -1, -1, -1, 0, 0, 0)
                      * coeffs_n_[speciesIndex].dimensions()
                    )
                )
            )
            {
                WarningIn("craftsGasModel::readDict")
                    << iname << " variable in " << outputDictName_
                    << " dimensions " << dsi.dimensions() << ", should be "
                    << "equal to mass / (time * temperature * length) * "
                    << "coefficient " << coeffs_n_[speciesIndex].name()
                    << endl;
            }
            integralRTSgDt_[speciesIndex] = dsi;
        }

        if
        (
            dimensionSet::debug
         && (
                dss.dimensions()
             != (dimMoles / dimVolume * coeffs_n_[speciesIndex].dimensions())
            )
        )
        {
            WarningIn("craftsGasModel::readDict")
                << sname << " variable in " << outputDictName_
                << " dimensions " << dss.dimensions() << ", should be equal "
                << "moles / volume * coefficient "
                << coeffs_n_[speciesIndex].name() << endl;
        }
        sGas_[speciesIndex] = dss;
        sGasCoeffs_[speciesIndex] = dss.value();
        sGasCoeffsOld_[speciesIndex] = dss.value();
    }

    if
    (
        dimensionSet::debug
     && (coeff_pGasOverPressure_.dimensions() != dimPressure)
    )
    {
        WarningIn("craftsGasModel::readDict")
            << "Coefficient " << coeff_pGasOverPressure_.name()
            << " has dimensions " << coeff_pGasOverPressure_.dimensions()
            << ", should be " << dimPressure << endl;
    }
    coeff_pGasOverPressure_.dimensions() = dimPressure;

}


template<int matrixSize>
void Foam::craftsGasModel<matrixSize>::checkUniformAndConstant() const
{
    if (!coeff_V_.uniform() || !coeff_V_.constant())
    {
        FatalErrorIn("craftsGasModel::craftsGasModel")
            << "In gas model for surface " << surfacePatchName_ << ", "
            << "gasVolume coefficient " << coeff_V_.name() << " must be "
            << "uniform and constant."
            << abort(FatalError);
    }
    if (!coeff_R_.uniform() || !coeff_R_.constant())
    {
        FatalErrorIn("craftsGasModel::craftsGasModel")
            << "In gas model for surface " << surfacePatchName_ << ", "
            << "R coefficient " << coeff_R_.name() << " must be "
            << "uniform and constant."
            << abort(FatalError);
    }
    if (!coeff_kp_.uniform() || !coeff_kp_.constant())
    {
        FatalErrorIn("craftsGasModel::craftsGasModel")
            << "In gas model for surface " << surfacePatchName_ << ", "
            << "k_p coefficient " << coeff_kp_.name()
            << " must be uniform and constant."
            << abort(FatalError);
    }
    if (!coeff_kGasCutOff_.uniform() || !coeff_kGasCutOff_.constant())
    {
        FatalErrorIn("craftsGasModel::craftsGasModel")
            << "In gas model for surface " << surfacePatchName_ << ", "
            << "kGasCutOff coefficient " << coeff_kGasCutOff_.name()
            << " must be uniform and constant."
            << abort(FatalError);
    }
    if
    (
        !coeff_pGasOverPressure_.uniform()
     || !coeff_pGasOverPressure_.constant()
    )
    {
        FatalErrorIn("craftsGasModel::craftsGasModel")
            << "In gas model for surface " << surfacePatchName_ << ", "
            << "pGasOverPressure coefficient "
            << coeff_pGasOverPressure_.name()
            << " must be uniform and constant."
            << abort(FatalError);
    }
    if (!coeff_pAtm_.uniform())
    {
        FatalErrorIn("craftsGasModel::craftsGasModel")
            << "In gas model for surface " << surfacePatchName_ << ", "
            << "p_atm coefficient " << coeff_pAtm_.name() << " must be "
            << "uniform."
            << abort(FatalError);
    }

    forAll(names_, i)
    {
        if (!coeffs_n_[i].uniform() || !coeffs_n_[i].constant())
        {
            FatalErrorIn("craftsGasModel::craftsGasModel")
                << "In gas model for surface " << surfacePatchName_ << ", "
                << "conversionToMoles coefficient " << coeffs_n_[i].name()
                << " for gas species " << names_[i] << " must be uniform "
                << "and constant."
                << abort(FatalError);
        }
        if (!coeffs_kla_[i].uniform() || !coeffs_kla_[i].constant())
        {
            FatalErrorIn("craftsGasModel::craftsGasModel")
                << "In gas model for surface " << surfacePatchName_ << ", "
                << "k_L_a coefficient " << coeffs_kla_[i].name()
                << " must be uniform and constant."
                << abort(FatalError);
        }
    }
}


template<int matrixSize>
void Foam::craftsGasModel<matrixSize>::checkDimensions() const
{
    forAll(names_, speciesIndex)
    {
        if (dimensionSet::debug)
        {
            if
            (
                vars_[speciesIndex].dimensions()
             != (
                    coeffs_n_[speciesIndex].dimensions()
                  * coeffs_kh_[speciesIndex].dimensions()
                  * pPartial_[speciesIndex].dimensions()
                )
            )
            {
                // Warn the user, as OpenFOAM's dimension error is too vague
                WarningIn("craftsGasModel::checkDimensions")
                    << "In gas model for surface " << surfacePatchName_ << ", "
                    << "for species " << names_[speciesIndex] << ", dimension "
                    << "error:\n\tconversionToMoles * k_H * [Pascals] must "
                    << "equal dimensions of liquidVariable." << endl;
                
                // Next line will throw the error
                vars_[speciesIndex].dimensions() =
                    coeffs_n_[speciesIndex].dimensions()
                  * coeffs_kh_[speciesIndex].dimensions()
                  * pPartial_[speciesIndex].dimensions();
            }
            else if
            (
                (coeffs_kla_[speciesIndex].dimensions() != (dimless / dimTime))
             || (coeffs_kb_[speciesIndex].dimensions() != (dimless / dimTime))
            )
            {
                // Warn the user, as OpenFOAM's dimension error is too vague
                WarningIn("craftsGasModel::checkDimensions")
                    << "In gas model for surface " << surfacePatchName_ << ", "
                    << "for species " << names_[speciesIndex] << ", dimension "
                    << "error:\n\tk_L_a and k_b dimensions must be [1/s]."
                    << endl;
                
                // Next lines will throw the error
                coeffs_kla_[speciesIndex].dimensions() = dimless / dimTime;
                coeffs_kb_[speciesIndex].dimensions() = dimless / dimTime;
            }
            else if
            (
                (
                    sGas_[speciesIndex].dimensions()
                  * coeff_R_.dimensions()
                  * var_T_.dimensions()
                  / coeffs_n_[speciesIndex].dimensions()
                ) != pPartial_[speciesIndex].dimensions()
            )
            {
                // Warn the user, as OpenFOAM's dimension error is too vague
                WarningIn("craftsGasModel::checkDimensions")
                    << "In gas model for surface " << surfacePatchName_ << ", "
                    << "for species " << names_[speciesIndex] << ", dimension "
                    << "error:\n\tR * T must be equal to [Pressure / moles]."
                    << endl;
                
                // Next line will throw the error
                pPartial_[speciesIndex].dimensions() =
                    sGas_[speciesIndex].dimensions()
                  * coeff_R_.dimensions()
                  * var_T_.dimensions()
                  / coeffs_n_[speciesIndex].dimensions();
            }
        }
    }
}


template<int matrixSize>
void Foam::craftsGasModel<matrixSize>::bubbleTransfer(const label speciesIndex)
{
    scalar dsGas(0.0);

    // R * sGas / n
    scalar rsgn
    (
        coeff_R_.evaluate(0)
      * sGasCoeffs_[speciesIndex]
      / coeffs_n_[speciesIndex].evaluate(0)
    );
    
    // Solubility limit
    /*scalarField sc
    (
        coeffs_kh_[speciesIndex].evaluateField()
      * rsgn
      * var_T_.evaluateField()
    );*/

    forAll(vars_[speciesIndex].evaluateField(), cellIndex)
    {
        // Local copy of liquid variable
        scalar sliq(varLiquidStart_[speciesIndex][cellIndex]);

        // Solubility limit
        scalar sc
        (
            coeffs_kh_[speciesIndex].evaluate(cellIndex)
          * rsgn
          * var_T_.evaluate(cellIndex)
        );
        
        // Only form bubbles when over solubility limit
        if (sliq >= sc)
        {
            scalar deltaSliq
            (
                (sc - sliq)
              * coeffs_kb_[speciesIndex].evaluate(cellIndex)
            );
            scalar newSliq
            (
                sliq + deltaSliq * runTime_.deltaT().value()
            );
            
            // Ensure we don't remove too much (due to large kb or timestep)
            if (newSliq < sc)
            {
                deltaSliq = (sc - sliq) / runTime_.deltaT().value();
                //sGas_[speciesIndex].value() += sliq - sc;
                //sliq = sc;
            }
            sliq += deltaSliq * runTime_.deltaT().value();
            vars_[speciesIndex].assign(sliq, cellIndex);
            dsGas -=
                deltaSliq * mesh_.V()[cellIndex]
              * runTime_.deltaT().value();
        }
    }
    bubbleSource_[speciesIndex]
        += dsGas / coeff_V_.evaluate(0);
}


template<int matrixSize>
const dimensionedScalar Foam::craftsGasModel<matrixSize>::calculateNewT() const
{
    const fvPatchScalarField& surfaceTempField
    (
        var_T_.evaluateGeometricField().boundaryField()[surfacePatchIndex_]
    );
    return dimensionedScalar
    (
        "T",
        var_T_.dimensions(),
        average(surfaceTempField)
    );
}


template<int matrixSize>
const dimensionedScalar Foam::craftsGasModel<matrixSize>::calculateNewPLiq()
    const
{
    // Get the average of coeff_pLiq_ on the boundary - but since coefficients
    // don't have boundaries, the best we can do for now is the average across
    // the cells next to the boundary.
    const fvPatch& sf(mesh_.boundary()[surfacePatchIndex_]);

    scalar avg(0.0);

    // This loop can be optimized as a reverse pointer loop, but I won't bother
    forAll(sf, faceI)
    {
        // cellIndex is the index of the neighbouring cell volume
        label cellIndex(sf.faceCells()[faceI]);
        avg += coeff_pLiq_.evaluate(cellIndex);
    }
    avg /= sf.size();
    
    return dimensionedScalar
    (
        coeff_pLiq_.name(),
        coeff_pLiq_.dimensions(),
        avg
    );
}


template<int matrixSize>
Foam::scalar Foam::craftsGasModel<matrixSize>::testMassConvergence() const
{
    scalar maxResid(-VGREAT);

    forAll(names_, speciesIndex)
    {
        // Calculate mass transfer to the gas
        const labelList& surfaceNeighbourIndices
        (
            vars_[speciesIndex].evaluateGeometricField()
                .boundaryField()[surfacePatchIndex_].patch().faceCells()
        );

        scalar surfaceVolume(0.0);
        forAll(surfaceNeighbourIndices, listI)
        {
            label cellIndex(surfaceNeighbourIndices[listI]);
            surfaceVolume += mesh_.V().field()[cellIndex];
        }

        scalar meanKh
        (
            gSum
            (
                mesh_.V() * coeffs_kh_[speciesIndex].evaluateField()
            ) / gSum(mesh_.V())
        );
        scalar meanKhOld
        (
            gSum
            (
                mesh_.V() * khOld_[speciesIndex]
            ) / gSum(mesh_.V())
        );

        scalar sLiqAve
        (
            vars_[speciesIndex].evaluateGeometricField().weightedAverage
            (
                mesh_.V()
            ).value()
        );
        scalar sLiqOldAve
        (
            vars_[speciesIndex].evaluateGeometricField().oldTime()
                .weightedAverage
                (
                    mesh_.V()
                ).value()
        );

        scalar mGas
        (
            coeffs_kla_[speciesIndex].evaluate(0)
          * surfaceVolume
          * (
                (sLiqAve + sLiqOldAve) / 2 * runTime_.deltaT().value()
              - (meanKh + meanKhOld) / 2
              * integralRTSgDt_[speciesIndex].value()
            )
        );

        // Finally, compare the two

        // Check if either is near zero.  If so, use nearZeroErrorScale
        if
        (
            (mag(mLiq_[speciesIndex]) < nearZeroMassScale_[speciesIndex])
         || (mag(mGas) < nearZeroMassScale_[speciesIndex])
        )
        {
            maxResid = max
            (
                maxResid,
                mag
                (
                    (mLiq_[speciesIndex] - mGas)
                    / nearZeroMassScale_[speciesIndex]
                ) - massConvergence_[speciesIndex]
            );
        }
        else
        {
            maxResid = max
            (
                maxResid,
                mag
                (
                    2 * (mLiq_[speciesIndex] - mGas)
                    / (mLiq_[speciesIndex] + mGas)
                ) - massConvergence_[speciesIndex]
            );
        }
    }
    return maxResid;
}


template<int matrixSize>
Foam::scalar Foam::craftsGasModel<matrixSize>::testDeltaConvergence() const
{
    scalar maxResid(-VGREAT);
    
    forAll(names_, speciesIndex)
    {
        if
        (
            (
                mag(sGasCoeffs_[speciesIndex])
              < nearZeroDeltaScale_[speciesIndex]
            )
         || (
                mag(sGasCoeffsNew_[speciesIndex])
              < nearZeroDeltaScale_[speciesIndex]
            )
        )
        {
            maxResid = max
            (
                maxResid,
                mag
                (
                    (sGasCoeffs_[speciesIndex] - sGasCoeffsNew_[speciesIndex])
                    / nearZeroDeltaScale_[speciesIndex]
                ) - deltaConvergence_[speciesIndex]
            );
        }
        else
        {
            maxResid = max
            (
                maxResid,
                mag
                (
                    (sGasCoeffs_[speciesIndex] - sGasCoeffsNew_[speciesIndex])
                    / sGasCoeffs_[speciesIndex]
                ) - deltaConvergence_[speciesIndex]
            );
        }
    }
    return maxResid;
}


template<int matrixSize>
void Foam::craftsGasModel<matrixSize>::calculatePressuresAndFlow
(
    const craftsGasFlowOde<matrixSize>& gasOde
)
{
    pTotal_.value() = 0;

    forAll(sGas_, speciesIndex)
    {
        sGas_[speciesIndex].value() = sGasCoeffs_[speciesIndex];
        pPartial_[speciesIndex].value() =
            sGas_[speciesIndex].value() * coeff_R_.evaluate(0)
          * T_.value() / coeffs_n_[speciesIndex].evaluate(0);
        pTotal_ += pPartial_[speciesIndex];
        
        // Output integral of R * T * sGas dt from t_old to t_new:
        integralRTSgDt_[speciesIndex].value() =
            gasOde.integralRTSgDt()[speciesIndex];
    }
    pTotal_ += pLiq_;

    scalar deltaP(pTotal_.value() - coeff_pAtm_.evaluate(0));
    scalar Hexp
    (
        -2 * coeff_kGasCutOff_.evaluate(0)
      * (deltaP - coeff_pGasOverPressure_.evaluate(0))
    );
    if (mag(Hexp) >= log(VGREAT / 2)) 
    {
        Hexp = pos(Hexp) * VGREAT / 2;
    }
    else
    {
        Hexp = exp(Hexp);
    }
    qGas_.value() = coeff_kp_.evaluate(0)
      * deltaP * pTotal_.value() / coeff_pAtm_.evaluate(0) / (1 + Hexp);
}


template<int matrixSize>
void Foam::craftsGasModel<matrixSize>::applyLimits()
{
    max
    (
        sGasCoeffsNew_,
        sGasCoeffsNew_,
        lowerLimits_
    );
    min
    (
        sGasCoeffsNew_,
        sGasCoeffsNew_,
        upperLimits_
    );
}


template<int matrixSize>
void Foam::craftsGasModel<matrixSize>::relax()
{
    sGasCoeffs_ =
        sGasCoeffs_ * (1 - relaxationFactors_)
      + sGasCoeffsNew_ * relaxationFactors_;
}


template<int matrixSize>
void Foam::craftsGasModel<matrixSize>::updateOutputDictionary()
{
    outputDict_.clear();
    outputDict_.set(T_.name(), T_);
    outputDict_.set(pTotal_.name(), pTotal_);
    forAll(pPartial_, speciesIndex)
    {
        outputDict_.set
        (
            pPartial_[speciesIndex].name(), pPartial_[speciesIndex]
        );
    }
    outputDict_.set(pLiq_.name(), pLiq_);
    outputDict_.set(coeff_pAtm_.name(), coeff_pAtm_.evaluateDimensioned(0));
    outputDict_.set(qGas_.name(), qGas_);
    forAll(sGas_, speciesIndex)
    {
        outputDict_.set
        (
            sGas_[speciesIndex].name(), sGas_[speciesIndex]
        );
    }
    forAll(integralRTSgDt_, speciesIndex)
    {
        outputDict_.set
        (
            integralRTSgDt_[speciesIndex].name(),
            integralRTSgDt_[speciesIndex]
        );
    }
}

// * * * * * * * * * * * * * * * * Constructors  * * * * * * * * * * * * * * //

template<int matrixSize>
Foam::craftsGasModel<matrixSize>::craftsGasModel
(
    admCoefficientManager& admCoeffs,
    admVariableManager& admVars,
    craftsModel<matrixSize>& model,
    const word& craftsGasModelDictName
)
:
    admCoeffs_(admCoeffs),
    admVars_(admVars),
    runTime_(admVars.runTime()),
    mesh_(admVars.mesh()),
    model_(model),

    craftsGasModelDict_
    (
        IOobject
        (
            craftsGasModelDictName,
            runTime_.constant(),
            mesh_,
            IOobject::MUST_READ,
            IOobject::NO_WRITE
        )
    ),

    reportPerformance_
    (
        model_.outputFlagsDict().found("gasSolverPerformance")
      ? bool(Switch(model_.outputFlagsDict().lookup("gasSolverPerformance")))
      : true
    ),

    odeSolverName_
    (
        craftsGasModelDict_.subDict("universal").lookup("odeSolver")
    ),
    
    odeSolverMaxIter_
    (
        craftsGasModelDict_
            .subDict("universal")
            .found("odeSolverMaxIterations")
      ? readLabel
        (
            craftsGasModelDict_.subDict("universal").lookup
            (
                "odeSolverMaxIterations"
            )
        )
      : label(10000)
    ),

    eps_
    (
        readScalar
        (
            craftsGasModelDict_.subDict("universal").lookup("epsilon")
        )
    ),

    defaultMassConvergence_
    (
        craftsGasModelDict_.subDict("universal").found
        (
            "defaultMassConvergence"
        )
      ? readScalar
        (
            craftsGasModelDict_.subDict("universal").lookup
            (
                "defaultMassConvergence"
            )
        )
      : -VGREAT
    ),

    defaultDeltaConvergence_
    (
        craftsGasModelDict_.subDict("universal").found
        (
            "defaultDeltaConvergence"
        )
      ? readScalar
        (
            craftsGasModelDict_.subDict("universal").lookup
            (
                "defaultDeltaConvergence"
            )
        )
      : -VGREAT
    ),

    defaultNearZeroMassScale_
    (
        craftsGasModelDict_.subDict("universal").found
        (
            "nearZeroMassScale"
        )
      ? readScalar
        (
            craftsGasModelDict_.subDict("universal").lookup
            (
                "nearZeroMassScale"
            )
        )
      : SMALL
    ),

    defaultNearZeroDeltaScale_
    (
        craftsGasModelDict_.subDict("universal").found
        (
            "nearZeroDeltaScale"
        )
      ? readScalar
        (
            craftsGasModelDict_.subDict("universal").lookup
            (
                "nearZeroDeltaScale"
            )
        )
      : SMALL
    ),

    outputDictName_
    (
        craftsGasModelDict_.subDict("universal").found("outputDictionary")
      ? craftsGasModelDict_.subDict("universal").lookup("outputDictionary")
      : word("gasModel")
    ),

    surfacePatchName_
    (
        craftsGasModelDict_.subDict("universal").lookup("surfacePatch")
    ),
    surfacePatchIndex_(mesh_.boundaryMesh().findPatchID(surfacePatchName_)),

    coeff_V_
    (
        admCoeffs_.readCoefficient
        (
            craftsGasModelDict_.subDict("universal"),
            word("gasVolume"),
            word(surfacePatchName_ + "(gasVolume)")
        )
    ),

    coeff_R_
    (
        admCoeffs_.readCoefficient
        (
            craftsGasModelDict_.subDict("universal"),
            word("R"),
            word(surfacePatchName_ + "(R)")
        )
    ),

    var_T_
    (
        admVars_.lookup
        (
            word
            (
                craftsGasModelDict_.subDict("universal").lookup("T")
            )
        )
    ),

    coeff_kp_
    (
        admCoeffs_.readCoefficient
        (
            craftsGasModelDict_.subDict("universal"),
            word("k_p"),
            word(surfacePatchName_ + "(k_p)")
        )
    ),

    coeff_kGasCutOff_
    (
        admCoeffs_.readCoefficient
        (
            craftsGasModelDict_.subDict("universal"),
            word("k_gasCutOff"),
            word(surfacePatchName_ + "(k_gasCutOff)")
        )
    ),

    coeff_pGasOverPressure_
    (
        admCoeffs_.readCoefficient
        (
            craftsGasModelDict_.subDict("universal"),
            word("p_gasOverPressure"),
            word(surfacePatchName_ + "(p_gasOverPressure)")
        )
    ),

    coeff_pAtm_
    (
        admCoeffs_.readCoefficient
        (
            craftsGasModelDict_.subDict("universal"),
            word("p_atm"),
            word(surfacePatchName_ + "(p_atm)")
        )
    ),

    coeff_pLiq_
    (
        admCoeffs_.readCoefficient
        (
            craftsGasModelDict_.subDict("universal"),
            word("p_liq"),
            word(surfacePatchName_ + "(p_liq)")
        )
    ),
    
    names_(craftsGasModelDict_.subDict("species").toc()),

    vars_(names_.size()),
    bcVars_(names_.size()),
    varLiquidStart_(names_.size()),
    formBubbles_(names_.size(), false),
    coeffs_n_(names_.size()),
    coeffs_kla_(names_.size()),
    coeffs_kh_(names_.size()),
    coeffs_kb_(names_.size()),
    upperLimits_(names_.size()),
    lowerLimits_(names_.size()),
    relaxationFactors_(names_.size()),
    massConvergence_(names_.size()),
    deltaConvergence_(names_.size()),
    nearZeroMassScale_(names_.size()),
    nearZeroDeltaScale_(names_.size()),
    
    outputDict_
    (
        IOobject
        (
            outputDictName_,
            runTime_.timeName(),
            mesh_,
            IOobject::MUST_READ,
            IOobject::AUTO_WRITE
        )
    ),

    T_
    (
        "T",
        var_T_.dimensions(),
        0.0
    ),

    pTotal_
    (
        "p_total",
        coeff_pAtm_.dimensions(),
        0.0
    ),

    pPartial_(names_.size()),

    pLiq_
    (
        coeff_pLiq_.name(),
        coeff_pLiq_.dimensions(),
        0.0
    ),
    
    integralRTSgDt_(names_.size()),

    qGas_
    (
        "Q",
        coeff_kp_.dimensions() * coeff_pAtm_.dimensions()
      * coeff_kGasCutOff_.dimensions(),
        0.0
    ),

    sGas_(names_.size()),

    Told_
    (
        "Told",
        var_T_.dimensions(),
        0.0
    ),

    pAtmOld_
    (
        coeff_pAtm_.name() + "Old",
        coeff_pAtm_.dimensions(),
        0
    ),

    pLiqOld_
    (
        coeff_pLiq_.name() + "Old",
        coeff_pLiq_.dimensions(),
        0.0
    ),
    
    khOld_(names_.size()),

    bubbleSource_(names_.size(), 0.0),
    sGasCoeffsOld_(names_.size(), 0.0),
    sGasCoeffs_(names_.size(), 0.0),
    sGasCoeffsNew_(names_.size(), 0.0),
    mLiq_(names_.size(), 0.0)
{
    // Fill the pointer lists
    allocateMemory();

    // Read values from the dictionary
    readDict();
    
    // Confirm coefficients meet required uniform() and constant() properties,
    // depending on coefficient
    checkUniformAndConstant();

    // Check dimensions of all coefficients and variables
    checkDimensions();
    
    Told_ = calculateNewT();
    pAtmOld_ = coeff_pAtm_.evaluateDimensioned(0);
    pLiqOld_ = calculateNewPLiq();
}


// * * * * * * * * * * * * * * Member Functions  * * * * * * * * * * * * * * //

template<int matrixSize>
bool Foam::craftsGasModel<matrixSize>::found(const word& speciesName)
{
    forAll(names_, i)
    {
        if (names_[i] == speciesName)
        {
            return true;
        }
    }
    return false;
}


template<int matrixSize>
Foam::label Foam::craftsGasModel<matrixSize>::lookup
(
    const word& speciesName
) const
{
    forAll(names_, i)
    {
        if (names_[i] == speciesName)
        {
            return i;
        }
    }
    return -1;
}


template<int matrixSize>
const Foam::fvPatchScalarField&
    Foam::craftsGasModel<matrixSize>::surfacePatchField
(
    const word& speciesName
) const
{
    label speciesIndex(lookup(speciesName));

    if (speciesIndex < 0)
    {
        FatalErrorIn("craftsGasModel::surfacePatchField")
            << speciesName << " not found in " << craftsGasModelDict_.name()
            << " dictionary."
            << abort(FatalError);
    }

    return surfacePatchField(speciesIndex);
}


template<int matrixSize>
const Foam::fvPatchScalarField&
    Foam::craftsGasModel<matrixSize>::surfacePatchField
(
    const label speciesIndex
) const
{
    return vars_[speciesIndex].evaluateGeometricField()
        .boundaryField()[surfacePatchIndex_];
}


template<int matrixSize>
void Foam::craftsGasModel<matrixSize>::initializeTimestep()
{
    // Temporarily store current value
    scalarField sGasTemp(sGasCoeffs_);
    
    // Extrapolate for next guess
    sGasCoeffs_ = sGasCoeffs_ +
        (sGasCoeffs_ - sGasCoeffsOld_)
      / runTime_.deltaT0().value()
      * (runTime_.deltaT0().value() + runTime_.deltaT().value());
    
    // Archive previous value
    sGasCoeffsOld_ = sGasTemp;
}


template<int matrixSize>
void Foam::craftsGasModel<matrixSize>::initializeImplicitLoop()
{
    forAll(names_, speciesIndex)
    {
        scalarField gradientField
        (
            bcVars_[speciesIndex].evaluateGeometricField().boundaryField()
                [surfacePatchIndex_].snGrad()
        );

        // Must account for ddy, if vars and bcVars are different
        labelList surfaceNeighbourIndices
        (
            vars_[speciesIndex].evaluateGeometricField()
                .boundaryField()[surfacePatchIndex_].patch().faceCells()
        );
        forAll(surfaceNeighbourIndices, i)
        {
            label cellIndex(surfaceNeighbourIndices[i]);
            gradientField[i] *=
                vars_[speciesIndex].ddy(bcVars_[speciesIndex], cellIndex);
        }

        const dimensionedScalar gamma
        (
            vars_[speciesIndex].lookup("diffusion")
        );
        const scalarField& magSf
        (
            vars_[speciesIndex].mesh().boundary()[surfacePatchIndex_].magSf()
        );
        
        mLiq_[speciesIndex] =
            sum(gradientField * magSf)
          * runTime_.deltaT().value() * -gamma.value();

        // Store varLiquidStart_
        varLiquidStart_[speciesIndex] = vars_[speciesIndex].evaluateField();
    }
}

template<int matrixSize>
int Foam::craftsGasModel<matrixSize>::implicitLoop()
{
    // Overall pseudo-code for the behaviour of this object
    // initializeTimestep()
    //  Set sGasOld equal to current value of sGas
    //  Also linearly interpolate for an initial guess of a new sGas
    //      sGasTemp = sGas
    //      sGas = linear interpolation:
    //          sGasOld is t-2
    //          sGas is t-1
    //          new sGas is t
    //      sGasOld = sGasTemp
    // initializeImplicitLoop()
    //  Calculate mass transfer from liquid volume
    //  *HERE* (label for pseudo-code "goto" below...)
    //  Reaction matrix is solved in craftsModel
    //      sLiq_i gradient is based on pPartial_i
    //  implicitLoop()
    //      calculate diffusion transfer
    //          add to bubbleSource
    //      calculate bubble transfer
    //          add to bubbleSource
    //          sLiq_i changes
    //              this is accommodated by the new 'hydrib variable' handling
    //      use an ode solver to get a new set of sGas's
    //          ode solver runs over the range
    //              tStart = previous timestep (runTime - runTime.deltaT)
    //              tEnd = current timestep (runTime)
    //          some ode parameters vary with time, are linearly interpolated:
    //              pAtm, Tgas, and pLiq
    //              interpolated between tStart and tEnd
    //          solution goes into sGasNew
    //      compare sGas with sGasNew
    //          convergence fail!
    //              relax sGas
    //              calculate partial pressures & gas flows
    //              update dictionary
    //              return -1 (repeat timestep --> return to *HERE* above)
    //          convergence sucess!
    //              sGas = sGasNew
    //              pAtmOld = pAtm
    //              TgasOld = tGas
    //              pLiqOld = pLiq
    //              calculate partial pressures & gas flows
    //              update dictionary

    // Calculate bubble transfer sources
    forAll(names_, speciesIndex)
    {
        if (formBubbles_[speciesIndex])
        {
            bubbleTransfer(speciesIndex);
        }
    }

    T_ = calculateNewT();
    pLiq_ = calculateNewPLiq();
    sGasCoeffsNew_ = sGasCoeffsOld_;

    scalarField n(names_.size());
    scalarField kla(names_.size());
    UPtrList<const volScalarField> slFields(names_.size());
    PtrList<scalarField> khEnd(names_.size());
    forAll(n, speciesIndex)
    {
        n[speciesIndex] = coeffs_n_[speciesIndex].evaluate(0);
        kla[speciesIndex] = coeffs_kla_[speciesIndex].evaluate(0);
        slFields.set
        (
            speciesIndex,
            & vars_[speciesIndex].evaluateGeometricField()
        );
        khEnd.set
        (
            speciesIndex,
            new scalarField
            (
                coeffs_kh_[speciesIndex].evaluateField()
            )
        );
    }

    craftsGasFlowOde<matrixSize> solveMe
    (
        model_,
        slFields,
        surfacePatchIndex_,
        runTime_.value() - runTime_.deltaT().value(),
        runTime_.value(),
        coeff_V_.evaluate(0),
        coeff_R_.evaluate(0),
        coeff_kp_.evaluate(0),
        coeff_kGasCutOff_.evaluate(0),
        coeff_pGasOverPressure_.evaluate(0),
        Told_.value(),
        T_.value(),
        pAtmOld_.value(),
        coeff_pAtm_.evaluate(0),
        pLiqOld_.value(),
        pLiq_.value(),
        kla,
        khOld_,
        khEnd,
        n,
        bubbleSource_,
        upperLimits_,
        lowerLimits_,
        sGasCoeffsNew_
    );

    autoPtr<ODESolver> gasOdeSolver
    (
        ODESolver::New
        (
            odeSolverName_,
            solveMe
        )
    );
    
    gasOdeSolver->maxStep() = odeSolverMaxIter_;

    label itersTaken(0);
    label odeError
    (
        gasOdeSolver->solve
        (
            runTime_.value() - runTime_.deltaT().value(),
            runTime_.value(),
            eps_,
            runTime_.deltaT().value(),
            itersTaken
        )
    );
    if (odeError)
    {
        // Raise error to calling function
        return -1;
    }

    applyLimits();

    scalar massResid(testMassConvergence());
    scalar deltaResid(testDeltaConvergence());
    if (reportPerformance_)
    {
        Info << "gasSolverPerformance (" << surfacePatchName_
            << "): n = " << itersTaken
            << ", massRes = " << massResid
            << ", deltaRes = " << deltaResid << endl;
    }

    // Add iterations to performance stats
    model_.lastStepIterations() += itersTaken;

    label returnMe(0);
    if (massResid > 0)
    {
        // Convergence failure - calculate CRM again
        returnMe = -2;
    }
    else if (deltaResid > 0)
    {
        // Stabilization failure - repeat implicitLoop
        returnMe = -3;
    }

    if (returnMe != 0)
    {
        // Failure - relax, update, and return
        relax();
        calculatePressuresAndFlow(solveMe);
        updateOutputDictionary();
        return returnMe;
    }

    // Success - update, and return
    sGasCoeffs_ = sGasCoeffsNew_;

    calculatePressuresAndFlow(solveMe);
    updateOutputDictionary();
    
    return 0;
}


template<int matrixSize>
void Foam::craftsGasModel<matrixSize>::finalizeImplicitLoop()
{
    // do nothing
}


template<int matrixSize>
void Foam::craftsGasModel<matrixSize>::finalizeTimestep()
{
    Told_ = T_;
    pAtmOld_ = coeff_pAtm_.evaluateDimensioned(0);
    pLiqOld_ = pLiq_;
    forAll(khOld_, speciesIndex)
    {
        khOld_[speciesIndex] = coeffs_kh_[speciesIndex].evaluateField();
    }
}


template<int matrixSize>
void Foam::craftsGasModel<matrixSize>::saveState(const label slot)
{
    if ((slot > nStates()) || (slot < 0))
    {
        FatalErrorIn("craftsGasModel::saveState")
            << "saveState slot " << slot << " out of range: 0, "
            << nStates() << ". Can only save over existing "
            << "slots or create one at end."
            << abort(FatalError);
    }
    if (slot == nStates())
    {
        // Create a new save state slot
        savedState_.setSize(slot + 1);
        savedState_[slot] = false;
        T_a.setSize(slot + 1);
        Told_a.setSize(slot + 1);
        pPartial_a.setSize(slot + 1);
        pAtmOld_a.setSize(slot + 1);
        pLiq_a.setSize(slot + 1);
        pLiqOld_a.setSize(slot + 1);
        integralRTSgDt_a.setSize(slot + 1);
        khOld_a.setSize(slot + 1);
        sGasCoeffs_a.setSize(slot + 1);
        sGasCoeffsOld_a.setSize(slot + 1);
        outputDict_a.setSize(slot + 1);
        mLiq_a.setSize(slot + 1);
    }
    if (!validState(slot))
    {
        // Allocate memory for slot
        T_a.set
        (
            slot,
            new dimensionedScalar(T_)
        );
        Told_a.set
        (
            slot,
            new dimensionedScalar(Told_)
        );
        pPartial_a.set
        (
            slot,
            new PtrList<dimensionedScalar>(names_.size())
        );
        forAll(pPartial_a[slot], i)
        {
            pPartial_a[slot].set
            (
                i,
                new dimensionedScalar
                (
                    "p_" + names_[i],
                    dimPressure,
                    0.0
                )
            );
        }
        pAtmOld_a.set
        (
            slot,
            new dimensionedScalar(pAtmOld_)
        );
        pLiq_a.set
        (
            slot,
            new dimensionedScalar(pLiq_)
        );
        pLiqOld_a.set
        (
            slot,
            new dimensionedScalar(pLiqOld_)
        );
        integralRTSgDt_a.set
        (
            slot,
            new PtrList<dimensionedScalar>(names_.size())
        );
        forAll(integralRTSgDt_a[slot], i)
        {
            integralRTSgDt_a[slot].set
            (
                i,
                new dimensionedScalar
                (
                    "interimValue_intRTSgDt_" + names_[i],
                    dimensionSet(1, -1, -1, -1, 0, 0, 0)
                        * coeffs_n_[i].dimensions(),
                    0.0
                )
            );
        }
        khOld_a.set
        (
            slot,
            new PtrList<scalarField>(names_.size())
        );
        forAll(khOld_a[slot], i)
        {
            khOld_a[slot].set
            (
                i,
                new scalarField
                (
                    coeffs_kh_[i].evaluateField()
                )
            );
        }
        sGasCoeffs_a.set
        (
            slot,
            new scalarField(sGasCoeffs_)
        );
        sGasCoeffsOld_a.set
        (
            slot,
            new scalarField(sGasCoeffsOld_)
        );
        outputDict_a.set
        (
            slot,
            new dictionary(outputDict_)
        );
        mLiq_a.set
        (
            slot,
            new scalarField(mLiq_)
        );
    }

    // Fill in state data
    savedState_[slot] = true;
    T_a[slot] = T_;
    Told_a[slot] = Told_;
    forAll(pPartial_a[slot], i)
    {
        pPartial_a[slot][i].value() = pPartial_[i].value();
        khOld_a[slot][i] = khOld_[i];
        integralRTSgDt_a[slot][i].value() = integralRTSgDt_[i].value();
    }
    pAtmOld_a[slot] = pAtmOld_;
    pLiq_a[slot] = pLiq_;
    pLiqOld_a[slot] = pLiqOld_;
    sGasCoeffs_a[slot] = sGasCoeffs_;
    sGasCoeffsOld_a[slot] = sGasCoeffsOld_;
    outputDict_a[slot].clear();
    outputDict_a[slot] <<= outputDict_;
    mLiq_a[slot] = mLiq_;
}


template<int matrixSize>
void Foam::craftsGasModel<matrixSize>::clearState(const label slot)
{
    if ((slot >= nStates()) || (slot < 0))
    {
        WarningIn("craftsGasModel::clearState")
            << "Attempting to clear non-existent slot " << slot << ". "
            << "Ignoring. Range = 0, " << nStates() << endl;
        return;
    }
    pPartial_a[slot].clear();
    khOld_a[slot].clear();
    integralRTSgDt_a[slot].clear();
    sGasCoeffs_a[slot].clear();
    sGasCoeffsOld_a[slot].clear();
    outputDict_a[slot].clear();
    savedState_[slot] = false;
    mLiq_a[slot].clear();
}


template<int matrixSize>
void Foam::craftsGasModel<matrixSize>::loadState(const label slot)
{
    if ((slot >= nStates()) || (slot < 0))
    {
        FatalErrorIn("craftsGasModel::loadState")
            << "Load slot " << slot << " out of range 0, " << nStates()
            << abort(FatalError);
    }
    if (!validState(slot))
    {
        FatalErrorIn("craftsGasModel::loadState")
            << "Attempting to load empty save state " << slot
            << abort(FatalError);
    }
    T_ = T_a[slot];
    Told_ = Told_a[slot];
    forAll(pPartial_, i)
    {
        pPartial_[i].value() = pPartial_a[slot][i].value();
        khOld_[i] = khOld_a[slot][i];
        integralRTSgDt_[i].value() = integralRTSgDt_a[slot][i].value();
    }
    pAtmOld_ = pAtmOld_a[slot];
    pLiq_ = pLiq_a[slot];
    pLiqOld_ = pLiqOld_a[slot];
    sGasCoeffs_ = sGasCoeffs_a[slot];
    sGasCoeffsOld_ = sGasCoeffsOld_a[slot];
    outputDict_.clear();
    outputDict_ <<= outputDict_a[slot];
    mLiq_ = mLiq_a[slot];
}


template<int matrixSize>
Foam::label Foam::craftsGasModel<matrixSize>::nStates() const
{
    return savedState_.size();
}


template<int matrixSize>
bool Foam::craftsGasModel<matrixSize>::validState(const label slot) const
{
    // Is slot out-of-range?
    if ((slot >= nStates()) || (slot < 0))
    {
        return false;
    }
    return savedState_[slot];
}

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

// ************************************************************************* //
