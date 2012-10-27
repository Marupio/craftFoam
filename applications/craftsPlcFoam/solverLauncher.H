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

#define craftsCaseAndSolveCommand(matrixSize)                                 \
        case matrixSize:                                                      \
            solve<matrixSize>(runTime, mesh, control);                        \
            break;

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

#   include "createAdmTime.H"
#   include "createMesh.H"

    // Find how many standard and implicit variables are defined
    int nVars(findTheMatrixSize(runTime));

    // Create crafts framework
    switch (nVars)
    {

#       define VectorNMacro craftsCaseAndSolveCommand
#       include "VectorNRunMacro.H"

        default:
#           define VectorNToc(n) \
                tocList.setSize(tocList.size() + 1); \
                tocList[tocList.size() - 1] = word("VectorN_" #n);
                
            wordList tocList(0);
#           define VectorNMacro VectorNToc
#           include "VectorNRunMacro.H"

            FatalErrorIn("craftsPlcFoam::solverLauncher")
                << "Matrix size (standard + implicit) of " << nVars
                << " requires a VectorN_" << nVars << " instantiation in the "
                << "VectorN library.  Available sizes are: " << tocList
                << "Ensure your matrix size meets one of these, or define the "
                << "required VectorN_## in the src/OpenFOAM/include/"
                << "customUserDefines.H and recompile all components."
                << abort(FatalError);
#           undef VectorNToc
            break;
    }
