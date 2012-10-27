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

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

// * * * * * * * * * * * * * Private Member Functions  * * * * * * * * * * * //

template <int matrixSize>
void Foam::craftsModel<matrixSize>::fillTransportTerms
(
    BlockLduMatrix<vectorType>& blockM,
    Field<vectorType>& blockX,
    Field<vectorType>& blockB
)
{
    // Add transport equations for standard variables.  Variables marked with
    // changedByUdf flag must be handled differently:
    // - must use Euler ddt
    // - must include udfDelta as source term
    // First, handle variables without the flag
    forAll(admVars_.unchangedByUdf(), i)
    {
        admStandardVariable& Svar(admVars_.unchangedByUdf(i));
        
        fvScalarMatrix SEqn
        (
            fvm::ddt(Svar())
          + fvm::div
            (
                phi_, Svar()
            )
          - fvm::laplacian
            (
                Svar.gamma(),
                Svar()
            )
        );

        blockMatrixTools::insertEquation
        (
            Svar.localIndex(),
            SEqn,
            blockM,
            blockX,
            blockB
        );
    }

    // Now for variables that are marked with changedByUdf flag
    
    // Create the hard coded Eulerian scheme
    dictionary ddtSchemeDict;
    ddtSchemeDict.set("EulerDdt", "Euler");
    tmp<fv::ddtScheme<scalar> > EulerDdt
    (
        fv::ddtScheme<scalar>::New
        (    
            mesh_,
            ddtSchemeDict.lookup("EulerDdt")
        )
    );

    forAll(admVars_.changedByUdf(), i)
    {
        admStandardVariable& Svar(admVars_.changedByUdf(i));

        DimensionedField<scalar, volMesh> dfUdfDelta
        (
            IOobject
            (
                "udfDelta",
                runTime_.timeName(),
                mesh_,
                IOobject::NO_READ,
                IOobject::NO_WRITE
            ),
            mesh_,
            Svar.dimensions() / dimTime,
            udfDelta_[i] / runTime_.deltaT().value()
        );

        fvScalarMatrix SEqn
        (
            EulerDdt->fvmDdt(Svar())
          + fvm::div
            (
                phi_, Svar()
            )
          - fvm::laplacian
            (
                Svar.gamma(),
                Svar()
            )
         == -dfUdfDelta
        );
        
        blockMatrixTools::insertEquation
        (
            Svar.localIndex(),
            SEqn,
            blockM,
            blockX,
            blockB
        );
    }
}


template <int matrixSize>
void Foam::craftsModel<matrixSize>::fillImplicitTerms
(
    BlockLduMatrix<vectorType>& blockM,
    Field<vectorType>& blockX,
    Field<vectorType>& blockB
)
{
    // The actual values for the implicit variables are calculated in implicit
    // functions, either implicitAutoSolve, or a user-defined function hook.
    // A first order linear approximation is created here for the main
    // coupledReactionMatrix.  This approximation is necessary so the standard
    // variables in the main matrix still get the effect of their dependence on
    // the implicit variables.

    // Create the hard code Eulerian scheme
    dictionary ddtSchemeDict;
    ddtSchemeDict.set("EulerDdt", "Euler");
    tmp<fv::ddtScheme<scalar> > EulerDdt
    (
        fv::ddtScheme<scalar>::New
        (    
            mesh_,
            ddtSchemeDict.lookup("EulerDdt")
        )
    );
    
    forAll(admVars_.implicit(), varIndex)
    {
        const admVariable& Svar(admVars_.implicit(varIndex));
        
        DimensionedField<scalar, volMesh> dfImplicitDdt
        (
            IOobject
            (
                "implicitDdt",
                runTime_.timeName(),
                mesh_,
                IOobject::NO_READ,
                IOobject::NO_WRITE
            ),
            mesh_,
            Svar.dimensions() / dimTime,
            implicitDdt_[varIndex]
        );
        
        fvScalarMatrix SiEqn
        (
            EulerDdt->fvmDdt(admVars_.implicit(varIndex)())
         == dfImplicitDdt
        );
        
        label matrixPosition(admVars_.nStandard() + varIndex);

        blockMatrixTools::insertEquation
        (
            matrixPosition,
            SiEqn,
            blockM,
            blockX,
            blockB
        );
    }
}


template <int matrixSize>
void Foam::craftsModel<matrixSize>::fillSDiagReactionTerms
(
    BlockLduMatrix<vectorType>& blockM,
    Field<vectorType>& blockB
)
{
    forAll(admVars_.standard(), varIndex)
    {
        const admVariable& baseVar(admVars_.standard(varIndex));
        
        // sDiagTermI:
        forAll(sDiagTermIk_[varIndex], termI)
        {
            const admReaction& reaction(sDiagTermIk_[varIndex][termI]);
/*scalar dvalue
(
    -mesh_.V()[0] * reaction.yield(baseVar).evaluate(0)
  * reaction.rate().ddy(baseVar, 0)
);
scalar svalue
(
    dvalue * baseVar.evaluate(0)
);
Info << varIndex << "," << termI << "," << dvalue << "," << svalue << endl;
svalue /= mesh_.V()[0];
dvalue /= mesh_.V()[0];
Info << varIndex << "," << termI << "," << dvalue << "," << svalue << endl;
Info << "Before: blockM = " << blockM.diag().asSquare()[0](varIndex, varIndex)
    << ", blockB = " << blockB[0](varIndex) << endl;*/
            blockMatrixTools::blockIncrement
            (
                varIndex,
                varIndex,
                scalarField
                (
                    -mesh_.V() * reaction.yield(baseVar).evaluateField()
                  * reaction.rate().ddyField(baseVar)
                ),
                blockM.diag().asSquare()
            );

            blockMatrixTools::blockIncrement
            (
                varIndex,
                scalarField
                (
                    -mesh_.V() * reaction.yield(baseVar).evaluateField()
                  * reaction.rate().ddyField(baseVar)
                  * baseVar.evaluateField()
                ),
                blockB
            );
/*Info << "After: blockM = " << blockM.diag().asSquare()[0](varIndex, varIndex)
    << ", blockB = " << blockB[0](varIndex) << endl;*/
        }
        
        // sDiagTermII:
        forAll(sDiagTermIIk_[varIndex], termII)
        {
            const admReaction& reaction(sDiagTermIIk_[varIndex][termII]);

            blockMatrixTools::blockIncrement
            (
                varIndex,
                varIndex,
                scalarField
                (
                    -mesh_.V() * reaction.yield(baseVar).ddyField(baseVar)
                  * reaction.rate().evaluateField()
                ),
                blockM.diag().asSquare()
            );

            blockMatrixTools::blockIncrement
            (
                varIndex,
                scalarField
                (
                    -mesh_.V() * reaction.yield(baseVar).ddyField(baseVar)
                  * reaction.rate().evaluateField()
                  * baseVar.evaluateField()
                ),
                blockB
            );
        } // end termII
    } // end forAll(vars)
}


template <int matrixSize>
void Foam::craftsModel<matrixSize>::fillOssNssNsiReactionTerms
(
    BlockLduMatrix<vectorType>& blockM,
    Field<vectorType>& blockB
)
{
    // oss terms
    forAll(ossTermIk_, m)
    {
        // Lower (owner) terms
        label wrtIndex(owner_[m]);
        label varIndex(neighbour_[m]);
        scalarField sum(mesh_.nCells(), 0.0);
        scalarField sourceSum(mesh_.nCells(), 0.0);
        forAll(ossTermIk_[m], termI)
        {
            const admVariable& wrtVarOss(admVars_.standard(wrtIndex));
            const admVariable& baseVarOss(admVars_.standard(varIndex));
            const admReaction& reaction(ossTermIk_[m][termI]);

            sum -=
                mesh_.V() * reaction.yield(baseVarOss).evaluateField()
              * reaction.rate().ddyField(wrtVarOss);
            
            sourceSum -=
                mesh_.V() * reaction.yield(baseVarOss).evaluateField()
              * reaction.rate().ddyField(wrtVarOss)
              * wrtVarOss.evaluateField();
        }

        forAll(ossTermIIk_[m], termII)
        {
            const admVariable& wrtVarOss(admVars_.standard(wrtIndex));
            const admVariable& baseVarOss(admVars_.standard(varIndex));
            const admReaction& reaction(ossTermIIk_[m][termII]);
            
            sum -=
                mesh_.V() * reaction.yield(baseVarOss).ddyField(wrtVarOss)
              * reaction.rate().evaluateField();
            
            sourceSum -=
                mesh_.V() * reaction.yield(baseVarOss).ddyField(wrtVarOss)
              * reaction.rate().evaluateField()
              * wrtVarOss.evaluateField();
        }
        
        // Add sum to term in lower triangle
        label col(owner_[m]);
        label row(neighbour_[m]);
        blockMatrixTools::blockIncrement
        (
            row,
            col,
            sum,
            blockM.diag().asSquare()
        );

        blockMatrixTools::blockIncrement
        (
            varIndex,
            sourceSum,
            blockB
        );

        // nss terms
        // varIndex and wrtIndex are now opposite - swap them
        Swap(varIndex, wrtIndex);
        
        sum = 0.0;
        sourceSum = 0.0;
        forAll(nssTermIk_[m], termI)
        {
            const admVariable& wrtVarNss(admVars_.standard(wrtIndex));
            const admVariable& baseVarNss(admVars_.standard(varIndex));
            const admReaction& reaction(nssTermIk_[m][termI]);

            sum -=
                mesh_.V() * reaction.yield(baseVarNss).evaluateField()
              * reaction.rate().ddyField(wrtVarNss);
            
            sourceSum -=
                mesh_.V() * reaction.yield(baseVarNss).evaluateField()
              * reaction.rate().ddyField(wrtVarNss)
              * wrtVarNss.evaluateField();
        }

        forAll(nssTermIIk_[m], termII)
        {
            const admVariable& wrtVarNss(admVars_.standard(wrtIndex));
            const admVariable& baseVarNss(admVars_.standard(varIndex));
            const admReaction& reaction(nssTermIIk_[m][termII]);
            sum -=
                mesh_.V() * reaction.yield(baseVarNss).ddyField(wrtVarNss)
              * reaction.rate().evaluateField();
            
            sourceSum -=
                mesh_.V() * reaction.yield(baseVarNss).ddyField(wrtVarNss)
              * reaction.rate().evaluateField()
              * wrtVarNss.evaluateField();
        }

        // nsi terms
        forAll(nsiTermIk_[m], termI)
        {
            const admVariable& baseVarNsi(admVars_.standard(varIndex));
            const admVariable& wrtVarNsi
            (
                admVars_.implicit(wrtIndex - admVars_.nStandard())
            );
            const admReaction& reaction(nsiTermIk_[m][termI]);
            
            sum -=
                mesh_.V() * reaction.yield(baseVarNsi).evaluateField()
              * reaction.rate().ddyField(wrtVarNsi);
            
            sourceSum -=
                mesh_.V() * reaction.yield(baseVarNsi).evaluateField()
                * reaction.rate().ddyField(wrtVarNsi)
              * wrtVarNsi.evaluateField();
        }

        forAll(nsiTermIIk_[m], termII)
        {
            const admVariable& baseVarNsi(admVars_.standard(varIndex));
            const admVariable& wrtVarNsi
            (
                admVars_.implicit(wrtIndex - admVars_.nStandard())
            );
            const admReaction& reaction(nsiTermIIk_[m][termII]);

            sum -=
                mesh_.V() * reaction.yield(baseVarNsi).ddyField(wrtVarNsi)
              * reaction.rate().evaluateField();
            
            sourceSum -=
                mesh_.V() * reaction.yield(baseVarNsi).ddyField(wrtVarNsi)
              * reaction.rate().evaluateField()
              * wrtVarNsi.evaluateField();
        }
        // Add sum to term in upper triangle
//        row = owner_[index];
//        col = neighbour_[index];
        // col and row are now in opposite positions as this is upper triangle
        blockMatrixTools::blockIncrement
        (
            col,
            row,
            sum,
            blockM.diag().asSquare()
        );

        blockMatrixTools::blockIncrement
        (
            varIndex,
            sourceSum,
            blockB
        );
    }
}


template <int matrixSize>
void Foam::craftsModel<matrixSize>::fillSourceReactionTerms
(
    BlockLduMatrix<vectorType>& blockM,
    Field<vectorType>& blockB
)
{
    // To enable dimension-checking for all reaction terms, it is only
    // necessary to dimension-check this function.  This is achieved using
    // dimensionedScalarFields for all quantities, including mesh volumes.
    dimensionedScalarField volume
    (
        IOobject
        (
            "volume",
            runTime_.timeName(),
            mesh_,
            IOobject::NO_READ,
            IOobject::NO_WRITE
        ),
        mesh_,
        dimVolume,
        mesh_.V()
    );

    forAll(admVars_.standard(), varIndex)
    {
        const admVariable& baseVar(admVars_.standard(varIndex));
        scalarField sum(mesh_.nCells(), 0.0);
        dimensionedScalarField dsfSum
        (
            IOobject
            (
                "sum_ddt(" + baseVar.name() + ")",
                runTime_.timeName(),
                mesh_,
                IOobject::NO_READ,
                IOobject::NO_WRITE
            ),
            mesh_,
            dimless,
            scalarField(mesh_.nCells(), 0.0)
        );

        // Instead of a forAll, we do the first reaction differently to
        // initialize dimensions - enables dimension checking
        if (nonZeroStandardYields_[varIndex].size())
        {
            const admReaction& reaction(nonZeroStandardYields_[varIndex][0]);

            dsfSum.dimensions().reset
            (
                dimVolume * reaction.yield(baseVar).evaluateDims()
              * reaction.rate().evaluateDims()
            );

            dsfSum.field() = mesh_.V()
                * reaction.yield(baseVar).evaluateField()
                * reaction.rate().evaluateField();
        }

        // Remaining reactions        
        for
        (
            label yieldIndex(1);
            yieldIndex < nonZeroStandardYields_[varIndex].size();
            yieldIndex++
        )
        {
            const admReaction& reaction
            (
                nonZeroStandardYields_[varIndex][yieldIndex]
            );

            dsfSum +=
                volume * reaction.yield(baseVar).evaluateDimensionedField()
              * reaction.rate().evaluateDimensionedField();
        }

        blockMatrixTools::blockIncrement
        (
            varIndex,
            dsfSum.field(),
            blockB
        );
    }
}

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

// ************************************************************************* //
