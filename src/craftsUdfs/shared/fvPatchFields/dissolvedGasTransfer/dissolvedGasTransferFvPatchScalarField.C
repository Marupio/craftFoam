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

#include "dissolvedGasTransferFvPatchScalarField.H"
#include "fvPatchFieldMapper.H"
#include "volFields.H"
#include "surfaceFields.H"
#include "addToRunTimeSelectionTable.H"
#include "Time.H"
#include "admReactionReader.H"

// * * * * * * * * * * * * * Private Member Functions  * * * * * * * * * * * //

bool Foam::dissolvedGasTransferFvPatchScalarField::findDictionaries()
{
    if 
    (
        !db().objectRegistry::foundObject<IOdictionary>(inputDictName_)
     || !db().objectRegistry::foundObject<IOdictionary>(outputDictName_)
    )
    {
        // Lookup failed
        if (failedLookups_ < 2)
        {
            failedLookups_++;
        }
        else
        {
            WarningIn("dissolvedGasTransferFvPatchScalarField::"
                "updateCoeffs()")
                << "Cannot find either " << inputDictName_ << " or "
                << outputDictName_ << " required by "
                << dimensionedInternalField().name() << " on patch "
                << patch().name() << ".  Assuming zero gradient."
                << endl;
        }
        return false;
    }

    inputDictPtr_ = & db().objectRegistry::lookupObject<IOdictionary>
    (
        inputDictName_
    );
    
    const dictionary& speciesDict
    (
        inputDictPtr_->subDict("species").subDict(species_)
    );

    outputDictPtr_ = & db().objectRegistry::lookupObject<IOdictionary>
    (
        outputDictName_
    );

    const admReactionReader& model
    (
        db().time().objectRegistry::lookupObject<admReactionReader>
        (
            "admReactionReader"
        )
    );

    varPtr_me_ =
        & model.admVars().lookup(dimensionedInternalField().name());

    coeffPtr_kLa_ =
        & model.admCoeffs()
        (
            word(speciesDict.lookup("k_L_a"))
        );
    
    coeffPtr_n_ =
        & model.admCoeffs()
        (
            word(speciesDict.lookup("conversionToMoles"))
        );
    
    coeffPtr_kH_ =
        & model.admCoeffs()
        (
            word(speciesDict.lookup("k_H"))
        );

    varPtr_sLiq_ =
        & model.admVars().lookup
        (
            word(speciesDict.lookup("liquidVariable"))
        );

    if (varPtr_me_->varType() == admVariable::vtderived)
    {
        FatalErrorIn("dissolvedGasTransferFvPatchScalarField"
            "::findDictionaries")
            << "volScalarField " << varPtr_me_->name() << " has a dissolved"
            << "GasTransfer boundary condition for species " << species_
            << ", but is a derived variable.  Assign this boundary condition "
            << "to one of the variables on which it depends."
            << abort(FatalError);
    }

    if (* varPtr_sLiq_ != * varPtr_me_)
    {
        // This is a hack for the case that sLiq is a derived variable
        if (varPtr_sLiq_->varType() != admVariable::vtderived)
        {
            FatalErrorIn("dissolvedGasTransferFvPatchScalarField"
                "::findDictionaries")
                << "volScalarField " << varPtr_me_->name() << " has a "
                << "dissolvedGasTransfer boundary condition for species "
                << species_ << ".  Species is associated with a different "
                << "liquid variable, " << varPtr_sLiq_->name() << ".  Only "
                << "valid if " << varPtr_sLiq_->name() << " is a derived "
                << "variable."
                << abort(FatalError);
        }
        else if (!varPtr_sLiq_->ddyNonZero(* varPtr_me_))
        {
            FatalErrorIn("dissolvedGasTransferFvPatchScalarField"
                "::findDictionaries")
                << "volScalarField " << varPtr_me_->name() << " has a "
                << "dissolvedGasTransfer boundary condition for species "
                << species_ << ", which is a derived variable. The derivative "
                << "of " << varPtr_sLiq_->name() << " w.r.t. "
                << varPtr_me_->name() << " must be non-zero."
                << abort(FatalError);
        }
        else
        {
            WarningIn("dissolvedGasTransferFvPatchScalarField"
                "::findDictionaries")
                << "volScalarField " << varPtr_me_->name() << " has a "
                << "dissolvedGasTransfer boundary condition for species "
                << species_ << ", which is a derived variable.  Derived "
                << "variables are okay, but they are implemented in a hard-"
                << "coded way, and are subject to these unchecked theoretical "
                << "limitations:"
                << token::NL << token::TAB
                << "* the derivative of " << varPtr_sLiq_->name()
                << " w.r.t. " << varPtr_me_->name() << " is non-zero and "
                << "spatially uniform;"
                << token::NL << token::TAB
                << "* all other variables that " << varPtr_sLiq_->name()
                << " depends on are zeroGradient at patch " << patch().name()
                << "." << token::NL << endl;
        }
    }
    
    dictionariesFound_ = true;
    return true;
}

// * * * * * * * * * * * * * * * * Constructors  * * * * * * * * * * * * * * //

Foam::dissolvedGasTransferFvPatchScalarField::
    dissolvedGasTransferFvPatchScalarField
(
    const fvPatch& p,
    const DimensionedField<scalar, volMesh>& iF
)
:
    fixedGradientFvPatchScalarField(p, iF),
    inputDictName_("craftsGasModelDict"),
    outputDictName_("gasModel"),
    failedLookups_(0),
    species_("default"),
    dictionariesFound_(false)
{}


Foam::dissolvedGasTransferFvPatchScalarField::
    dissolvedGasTransferFvPatchScalarField
(
    const dissolvedGasTransferFvPatchScalarField& ptf,
    const fvPatch& p,
    const DimensionedField<scalar, volMesh>& iF,
    const fvPatchFieldMapper& mapper
)
:
    fixedGradientFvPatchScalarField(ptf, p, iF, mapper),
    inputDictName_(ptf.inputDictName_),
    outputDictName_(ptf.outputDictName_),
    failedLookups_(ptf.failedLookups_),
    species_(ptf.species_),
    dictionariesFound_(ptf.dictionariesFound_)
{}


Foam::dissolvedGasTransferFvPatchScalarField::
    dissolvedGasTransferFvPatchScalarField
(
    const fvPatch& p,
    const DimensionedField<scalar, volMesh>& iF,
    const dictionary& dict
)
:
    fixedGradientFvPatchScalarField(p, iF),
    inputDictName_
    (
        dict.found("inputDictionary")
      ? dict.lookup("inputDictionary")
      : word("craftsGasModelDict")
    ),
    outputDictName_
    (
        dict.found("outputDictionary")
      ? dict.lookup("outputDictionary")
      : word("gasModel")
    ),
    failedLookups_(0),
    species_(dict.lookup("species")),
    dictionariesFound_(false)
{
    if (dict.found("gradient"))
    {
        gradient() = scalarField("gradient", dict, p.size());
        fixedGradientFvPatchScalarField::updateCoeffs();
        fixedGradientFvPatchScalarField::evaluate();
    }
    else
    {
        fvPatchField<scalar>::operator=(patchInternalField());
        gradient() = 0.0;
    }
}


Foam::dissolvedGasTransferFvPatchScalarField::
    dissolvedGasTransferFvPatchScalarField
(
    const dissolvedGasTransferFvPatchScalarField& wbppsf
)
:
    fixedGradientFvPatchScalarField(wbppsf),
    inputDictName_(wbppsf.inputDictName_),
    outputDictName_(wbppsf.outputDictName_),
    failedLookups_(wbppsf.failedLookups_),
    species_(wbppsf.species_),
    dictionariesFound_(wbppsf.dictionariesFound_)
{}


Foam::dissolvedGasTransferFvPatchScalarField::
    dissolvedGasTransferFvPatchScalarField
(
    const dissolvedGasTransferFvPatchScalarField& wbppsf,
    const DimensionedField<scalar, volMesh>& iF
)
:
    fixedGradientFvPatchScalarField(wbppsf, iF),
    inputDictName_(wbppsf.inputDictName_),
    outputDictName_(wbppsf.outputDictName_),
    failedLookups_(wbppsf.failedLookups_),
    species_(wbppsf.species_),
    dictionariesFound_(wbppsf.dictionariesFound_)
{}


// * * * * * * * * * * * * * * * Member Functions  * * * * * * * * * * * * * //

void Foam::dissolvedGasTransferFvPatchScalarField::
    updateCoeffs()
{
    if (updated())
    {
        return;
    }

    if (!dictionariesFound_ && !findDictionaries())
    {
        gradient() = 0.0;
        return;
    }

    dimensionedScalar gamma
    (
        varPtr_sLiq_->lookup("diffusion")
    );

    word iName("interimValue_intRTSgDt_" + species_);
    if (outputDictPtr_->found(iName))
    {
        // More accurate "integral" method available
        scalar integralRTSgDt
        (
            dimensionedScalar
            (
                outputDictPtr_->lookup(iName)
            ).value()
        );

        // Take average of dissolved gas
        scalar sLiqAve
        (
            varPtr_sLiq_->evaluateGeometricField().weightedAverage
            (
                dimensionedInternalField().mesh().V()
            ).value()
        );
        scalar sLiqOldAve
        (
            varPtr_sLiq_->evaluateGeometricField().oldTime().weightedAverage
            (
                dimensionedInternalField().mesh().V()
            ).value()
        );

        scalar meanKh
        (
            gSum
            (
                coeffPtr_kH_->evaluateField()
              * dimensionedInternalField().mesh().V()
              / gSum(dimensionedInternalField().mesh().V())
            )
        );

        // TODO This calculation needs non-orthogonality correction
        forAll(gradient(), faceIndex)
        {
            label cellIndex(patch().faceCells()[faceIndex]);
            gradient()[faceIndex] =
               -coeffPtr_kLa_->evaluate(cellIndex)
              * dimensionedInternalField().mesh().V().field()[cellIndex]
              / gamma.value() / patch().magSf()[faceIndex]
              * (
                    (sLiqAve + sLiqOldAve) / 2.0
                  - integralRTSgDt / db().time().deltaT().value()
                  * meanKh
                ) / varPtr_sLiq_->ddy(* varPtr_me_, cellIndex);
        }
    }
    else
    {
        // Resort to pressure method
        scalar p
        (
            dimensionedScalar
            (
                outputDictPtr_->lookup("p_" + species_)
            ).value()
        );

        scalar sLiqAve
        (
            varPtr_sLiq_->evaluateGeometricField().weightedAverage
            (
                dimensionedInternalField().mesh().V()
            ).value()
        );
        scalar sLiqOldAve
        (
            varPtr_sLiq_->evaluateGeometricField().oldTime().weightedAverage
            (
                dimensionedInternalField().mesh().V()
            ).value()
        );

        scalar meanKh
        (
            gSum
            (
                coeffPtr_kH_->evaluateField()
              * dimensionedInternalField().mesh().V()
              / gSum(dimensionedInternalField().mesh().V())
            )
        );

        // TODO This calculation needs non-orthogonality correction
        forAll(gradient(), faceIndex)
        {
            label cellIndex(patch().faceCells()[faceIndex]);
            gradient()[faceIndex] = -coeffPtr_kLa_->evaluate(cellIndex) *
            (
                (sLiqAve + sLiqOldAve) / 2.0
              - coeffPtr_n_->evaluate(cellIndex)
              * meanKh
              * p
            ) / gamma.value()
              * dimensionedInternalField().mesh().V().field()[cellIndex]
              / patch().magSf()[faceIndex]
              / varPtr_sLiq_->ddy(* varPtr_me_, cellIndex);
        }
    }
    fixedGradientFvPatchScalarField::updateCoeffs();
}


void Foam::dissolvedGasTransferFvPatchScalarField::
    write(Ostream& os) const
{
    fvPatchScalarField::write(os);
    if (inputDictName_ != "craftsGasModelDict")
    {
        os.writeKeyword("inputDictionary") << inputDictName_
            << token::END_STATEMENT << nl;
    }
    if (outputDictName_ != "gasModel")
    {
        os.writeKeyword("dictionary") << outputDictName_
            << token::END_STATEMENT << nl;
    }
    os.writeKeyword("species") << species_ << token::END_STATEMENT << nl;
    gradient().writeEntry("gradient", os);
}


// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

namespace Foam
{
    makePatchTypeField
    (
        fvPatchScalarField,
        dissolvedGasTransferFvPatchScalarField
    );
}

// ************************************************************************* //
