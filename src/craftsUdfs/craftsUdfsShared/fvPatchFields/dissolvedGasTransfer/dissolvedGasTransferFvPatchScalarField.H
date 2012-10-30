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
    Foam::dissolvedGasTransferFvPatchScalarField

Description
    Foam::dissolvedGasTransferFvPatchScalarField

SourceFiles
    dissolvedGasTransferFvPatchScalarField.C

\*---------------------------------------------------------------------------*/

#ifndef dissolvedGasTransferFvPatchScalarField_H
#define dissolvedGasTransferFvPatchScalarField_H

#include "fvPatchFields.H"
#include "fixedGradientFvPatchFields.H"
#include "admReactionReader.H"

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

namespace Foam
{

/*---------------------------------------------------------------------------*\
                           Class dissolvedGasTransferFvPatch Declaration
\*---------------------------------------------------------------------------*/

class dissolvedGasTransferFvPatchScalarField
:
    public fixedGradientFvPatchScalarField
{
    // Private data

        // Dictionary search
        //  This boundary condition depends on finding two dictionaries in the
        //  object registry.  These variables relate to that search.
    
            //- Name of the dictionary containing the coefficients and other
            //  gas model settings (in constant, default craftsGasModelDict)
            const word inputDictName_;
            
            //- Name of the dictionary containing the partial pressures (in
            //  runTime.name, default gasModel)
            const word outputDictName_;

            //- Number of times I failed to find both input and output
            //  dictionaries.  After 2 I start throwing warnings.
            label failedLookups_;

            //- Name of this boundary condition's gas species
            const word species_;

        // Pointers
        //  These variables are initialized only after the dictionaries have
        //  been found.
        
            //- True once dictionaries have been found
            bool dictionariesFound_;

            //- Pointer to input dictionary
            const dictionary * inputDictPtr_;

            //- Pointer to output dictionary
            const dictionary * outputDictPtr_;

            //- Pointer to the dissolved gas concentration variable
            const admVariable * varPtr_sLiq_;

            //- Pointer to the variable associated with this boundary condition
            //  (may differ from sLiq if this is a derived variabl)
            const admVariable * varPtr_me_;

            //- Pointer to overall diffusion transfer coefficient
            const admCoefficient * coeffPtr_kLa_;
            
            //- Pointer to the conversion to moles coefficient
            const admCoefficient * coeffPtr_n_;
            
            //- Pointer to Henry's equilibrium coefficient
            const admCoefficient * coeffPtr_kH_;

    // Private member variables
    
        //- Find and read the requisite dictionaries
        bool findDictionaries();

public:

    //- Runtime type information
    TypeName("dissolvedGasTransfer");


    // Constructors

        //- Construct from patch and internal field
        dissolvedGasTransferFvPatchScalarField
        (
            const fvPatch&,
            const DimensionedField<scalar, volMesh>&
        );

        //- Construct from patch, internal field and dictionary
        dissolvedGasTransferFvPatchScalarField
        (
            const fvPatch&,
            const DimensionedField<scalar, volMesh>&,
            const dictionary&
        );

        //- Construct by mapping given dissolvedGasTransferFvPatchScalarField onto
        //  a new patch
        dissolvedGasTransferFvPatchScalarField
        (
            const dissolvedGasTransferFvPatchScalarField&,
            const fvPatch&,
            const DimensionedField<scalar, volMesh>&,
            const fvPatchFieldMapper&
        );

        //- Construct as copy
        dissolvedGasTransferFvPatchScalarField
        (
            const dissolvedGasTransferFvPatchScalarField&
        );

        //- Construct and return a clone
        virtual tmp<fvPatchScalarField> clone() const
        {
            return tmp<fvPatchScalarField>
            (
                new dissolvedGasTransferFvPatchScalarField(*this)
            );
        }

        //- Construct as copy setting internal field reference
        dissolvedGasTransferFvPatchScalarField
        (
            const dissolvedGasTransferFvPatchScalarField&,
            const DimensionedField<scalar, volMesh>&
        );

        //- Construct and return a clone setting internal field reference
        virtual tmp<fvPatchScalarField> clone
        (
            const DimensionedField<scalar, volMesh>& iF
        ) const
        {
            return tmp<fvPatchScalarField>
            (
                new dissolvedGasTransferFvPatchScalarField(*this, iF)
            );
        }

        
    // Member functions

        //- Update the coefficients associated with the patch field
        virtual void updateCoeffs();

        //- Write
        virtual void write(Ostream&) const;

};


// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

} // End namespace Foam

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

#endif

// ************************************************************************* //
