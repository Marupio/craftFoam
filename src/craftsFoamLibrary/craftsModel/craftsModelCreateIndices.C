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
void Foam::craftsModel<matrixSize>::createLduMesh()
{
    for (label row(0); row < admVars_.nStandard(); row++)
    {
        const admVariable& rowVar(admVars_.standard(row));
        
        for
        (
            label col(row + 1);
            col < admVars_.nStandard() + admVars_.nImplicit();
            col++
        )
        {
/*            const admVariable * colVarPtr;
            if (col < admVars_.nStandard())
            {
                colVarPtr = & admVars_.standard(col);
            }
            else
            {
                colVarPtr = & admVars_.implicit(col - admVars_.nStandard());
            }
            const admVariable& colVar(* colVarPtr);*/
            bool nonZero(false);
            
            UPtrList<const admReaction> ossTermIk;
            UPtrList<const admReaction> ossTermIIk;
            UPtrList<const admReaction> nssTermIk;
            UPtrList<const admReaction> nssTermIIk;
            UPtrList<const admReaction> nsiTermIk;
            UPtrList<const admReaction> nsiTermIIk;

            if (col < admVars_.nStandard()) // then we are in nss/oss quadrant
            {
                const admVariable& colVar(admVars_.standard(col));
                // Create oss indices
                forAll(nonZeroStandardYields_[col], yieldIndex)
                {
                    const admReaction& reaction
                    (
                        nonZeroStandardYields_[col][yieldIndex]
                    );
                    // Term I (yield already guaranteed non-zero)
                    if
                    (
                        reaction.rate().ddyNonZero(rowVar)
                    )
                    {
                        label newIndex(ossTermIk.size());
                        ossTermIk.setSize(newIndex + 1);
                        ossTermIk.set(newIndex, &reaction);
                        nonZero = true;
                    } // end term I
                    
                    // Term II
                    if
                    (
                        reaction.yield(colVar).ddyNonZero(rowVar)
                     && reaction.rate().evaluateNonZero()
                    )
                    {
                        label newIndex(ossTermIIk.size());
                        ossTermIIk.setSize(newIndex + 1);
                        ossTermIIk.set(newIndex, &reaction);
                        nonZero = true;
                    } // end term II
                } // end k loop

                // Create nss indices
                forAll(nonZeroStandardYields_[row], yieldIndex)
                {
                    const admReaction& reaction
                    (
                        nonZeroStandardYields_[row][yieldIndex]
                    );

                    // Term I (yield already guaranteed non-zero)
                    if
                    (
                        reaction.rate().ddyNonZero(colVar)
                    )
                    {
                        label newIndex(nssTermIk.size());
                        nssTermIk.setSize(newIndex + 1);
                        nssTermIk.set(newIndex, &reaction);
                        nonZero = true;
                    } // end term I
                    
                    // Term II
                    if
                    (
                        reaction.yield(rowVar).ddyNonZero(colVar)
                     && reaction.rate().evaluateNonZero()
                    )
                    {
                        label newIndex(nssTermIIk.size());
                        nssTermIIk.setSize(newIndex + 1);
                        nssTermIIk.set(newIndex, &reaction);
                        nonZero = true;
                    } // end term II
                } // end k loop

            } // end if col < nStandard
            else // then we are in the nsi/osi quadrant
            {
                const admVariable& colVar
                (
                    admVars_.implicit(col - admVars_.nStandard())
                );

                // Create nsi indices
                forAll(nonZeroStandardYields_[row], yieldIndex)
                {
                    const admReaction& reaction
                    (
                        nonZeroStandardYields_[row][yieldIndex]
                    );

                    // Term I (yield already guaranteed non-zero)
                    if
                    (
                        reaction.rate().ddyNonZero(colVar)
                    )
                    {
                        label newIndex(nsiTermIk.size());
                        nsiTermIk.setSize(newIndex + 1);
                        nsiTermIk.set(newIndex, &reaction);
                        nonZero = true;
                    } // end term I
                    
                    // Term II
                    if
                    (
                        reaction.yield(rowVar).ddyNonZero(colVar)
                     && reaction.rate().evaluateNonZero()
                    )
                    {
                        label newIndex(nsiTermIIk.size());
                        nsiTermIIk.setSize(newIndex + 1);
                        nsiTermIIk.set(newIndex, &reaction);
                        nonZero = true;
                    } // end term II
                } // end k loop
            } // end else col >= nStandard

            // Add to owner, neighbour, and nsi lists
            if (nonZero)
            {
                label newIndex(owner_.size());

                owner_.setSize(newIndex + 1);
                owner_[newIndex] = row;

                neighbour_.setSize(newIndex + 1);
                neighbour_[newIndex] = col;
                
                ossTermIk_.setSize(newIndex + 1);
                ossTermIk_.set
                (
                    newIndex,
                    new UPtrList<const admReaction>(ossTermIk)
                );
                
                ossTermIIk_.setSize(newIndex + 1);
                ossTermIIk_.set
                (
                    newIndex,
                    new UPtrList<const admReaction>(ossTermIIk)
                );
                
                nssTermIk_.setSize(newIndex + 1);
                nssTermIk_.set
                (
                    newIndex,
                    new UPtrList<const admReaction>(nssTermIk)
                );
                
                nssTermIIk_.setSize(newIndex + 1);
                nssTermIIk_.set
                (
                    newIndex,
                    new UPtrList<const admReaction>(nssTermIIk)
                );

                nsiTermIk_.setSize(newIndex + 1);
                nsiTermIk_.set
                (
                    newIndex,
                    new UPtrList<const admReaction>(nsiTermIk)
                );
                
                nsiTermIIk_.setSize(newIndex + 1);
                nsiTermIIk_.set
                (
                    newIndex,
                    new UPtrList<const admReaction>(nsiTermIIk)
                );
            } // end if nonZero
        } // end col for loop
    } // end row for loop
    
    // Create sDiag indices
    forAll(admVars_.standard(), ij)
    {
        sDiagTermIk_.set
        (
            ij,
            new UPtrList<const admReaction>(0)
        );
        sDiagTermIIk_.set
        (
            ij,
            new UPtrList<const admReaction>(0)
        );

        const admVariable& ijVar(admVars_.standard(ij));

        forAll(nonZeroStandardYields_[ij], yieldIndex)
        {
            const admReaction& reaction
            (
                nonZeroStandardYields_[ij][yieldIndex]
            );
            
            // Term I (yield already guaranteed non-zero)
            if
            (
                reaction.rate().ddyNonZero(ijVar)
            )
            {
                label newIndex(sDiagTermIk_[ij].size());
                sDiagTermIk_[ij].setSize(newIndex + 1);
                sDiagTermIk_[ij].set(newIndex, &reaction);
            }
        
            // Term II
            if
            (
                reaction.yield(ijVar).ddyNonZero(ijVar)
             && reaction.rate().evaluateNonZero()
            )
            {
                label newIndex(sDiagTermIIk_[ij].size());
                sDiagTermIIk_[ij].setSize(newIndex + 1);
                sDiagTermIIk_[ij].set(newIndex, &reaction);
            } // end termII
        } // end reactions loop 
    } // ij loop
}


template <int matrixSize>
void Foam::craftsModel<matrixSize>::createTransportFaceIndices()
{
    for (label cellIndex(0); cellIndex < mesh_.nCells(); cellIndex++)
    {
        // Owner faces
        label i(0);
        while
        (
            (i < mesh_.owner().size())
         && (mesh_.owner()[i] != cellIndex)
        )
        {
            i++;
        }
        while
        (
            (i < mesh_.owner().size())
         && (mesh_.owner()[i] == cellIndex)
        )
        {
            label newIndex(transportOwners_[cellIndex].size());
            transportOwners_[cellIndex].setSize(newIndex + 1);
            transportOwners_[cellIndex][newIndex] = i;
            
            transportOwnersCells_[cellIndex].setSize(newIndex + 1);
            transportOwnersCells_[cellIndex][newIndex] =
                mesh_.neighbour()[i];
            i++;
        }
        
        // Neighbour faces
        forAll(mesh_.neighbour(), i)
        {
            if (mesh_.neighbour()[i] == cellIndex)
            {
                label newIndex(transportNeighbours_[cellIndex].size());
                transportNeighbours_[cellIndex].setSize(newIndex + 1);
                transportNeighbours_[cellIndex][newIndex] = i;
                transportNeighboursCells_[cellIndex].setSize(newIndex + 1);
                transportNeighboursCells_[cellIndex][newIndex] =
                    mesh_.owner()[i];
            }
        }
    }
}


template <int matrixSize>
void Foam::craftsModel<matrixSize>::createAutoSolveIndices()
{
    // Create autoSolve dRdS indices
    forAll(admVars_.implicitAutoSolve(), impIndex)
    {
        dRdSAutoSolveTermIk_.set
        (
            impIndex,
            new UPtrList<const admReaction>(0)
        );
        dRdSAutoSolveTermIIk_.set
        (
            impIndex,
            new UPtrList<const admReaction>(0)
        );

        const admVariable& impVar(admVars_.implicitAutoSolve(impIndex));

        forAll(nonZeroImplicitYields_[impIndex], yieldIndex)
        {
            const admReaction& reaction
            (
                nonZeroImplicitYields_[impIndex][yieldIndex]
            );

            // Term I (yield already guaranteed non-zero)
            if
            (
                reaction.rate().ddyNonZero(impVar)
            )
            {
                label newIndex(dRdSAutoSolveTermIk_[impIndex].size());
                dRdSAutoSolveTermIk_[impIndex].setSize(newIndex + 1);
                dRdSAutoSolveTermIk_[impIndex].set(newIndex, &reaction);
            }

            // Term II
            if
            (
                reaction.yield(impVar).ddyNonZero(impVar)
             && reaction.rate().evaluateNonZero()
            )
            {
                label newIndex(dRdSAutoSolveTermIIk_[impIndex].size());
                dRdSAutoSolveTermIIk_[impIndex].setSize(newIndex + 1);
                dRdSAutoSolveTermIIk_[impIndex].set(newIndex, &reaction);
            } // end termII
        } // end reactions loop 
    } // impIndex loop
}

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

// ************************************************************************* //
