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

Application
    craftsAverages

Description
    Amalgamates all 'averages' dictionary files from all timestep directories
    into a single gnuplot-ready file in the case root.

Author
    David L. F. Gaden

\*---------------------------------------------------------------------------*/

#include "timeSelector.H"
#include "fvCFD.H"

using namespace Foam;

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

int main(int argc, char *argv[])
{
    timeSelector::addOptions(false, true);
    argList::validOptions.insert("ds", "");
    argList::validOptions.insert("dict", "dictionary name");

#   include "setRootCase.H"
#   include "createTime.H"
    instantList timeDirs = timeSelector::select0(runTime, args);
#   include "createMesh.H"

    word dictName("averages");
    bool dsFormat(false);
    if (args.optionFound("dict"))
    {
        dictName = args.options()["dict"];
    }
    if (args.optionFound("ds"))
    {
        dsFormat = true;
    }

    runTime.setTime(timeDirs[0], 0);
    Info << "Time = " << runTime.timeName() << endl;
    
//    fileName averagesDictFile(runTime.timeName()/"averages");
//    dictionary averagesDict(averagesDictFile);
    wordList variables;
    scalarListList averages;
    {
        IOdictionary averagesDict
        (
            IOobject
            (
                dictName,
                runTime.timeName(),
                mesh,
                IOobject::MUST_READ,
                IOobject::NO_WRITE
            )
        );
    
        wordList allVariables(averagesDict.toc());
        forAll(allVariables, i)
        {
            if (allVariables[i](0,12) != word("interimValue"))
            {
                label newIndex(variables.size());
                variables.setSize(newIndex + 1);
                variables[newIndex] = allVariables[i];
            }
        }
        sort(variables);
        averages.setSize(variables.size());
    
        if (dsFormat)
        {
            forAll(variables, varIndex)
            {
                averages[varIndex].setSize(timeDirs.size());
                averages[varIndex][0] = dimensionedScalar
                (
                    averagesDict.lookup(variables[varIndex])
                ).value();
            }
        }
        else
        {
            forAll(variables, varIndex)
            {
                averages[varIndex].setSize(timeDirs.size());
                averages[varIndex][0] =
                    readScalar(averagesDict.lookup(variables[varIndex]));
            }
        }
    }

    for (label timeI = 1; timeI < timeDirs.size(); timeI++)
    {
        runTime.setTime(timeDirs[timeI], timeI);
        Info<< "Time = " << runTime.timeName() << endl;

        IOdictionary averagesDict
        (
            IOobject
            (
                dictName,
                runTime.timeName(),
                mesh,
                IOobject::MUST_READ,
                IOobject::NO_WRITE
            )
        );
        
        wordList averagesToc(averagesDict.toc());
        wordList averagesVars;
        forAll(averagesToc, i)
        {
            if (averagesToc[i](0,12) != word("interimValue"))
            {
                label newIndex(averagesVars.size());
                averagesVars.setSize(newIndex + 1);
                averagesVars[newIndex] = averagesToc[i];
            }
        }
        sort(averagesVars);

        if (averagesVars != variables)
        {
            WarningIn("craftsAverages::main")
                << dictName << " file content differs between timestep "
                << timeDirs[timeI].name() << ", and "
                << timeDirs[timeI - 1].name() << "." << endl;

            forAll(variables, varIndex)
            {
                averages[varIndex][timeI] = 0.0;
            }
            continue;
        }
        
        if (dsFormat)
        {
            forAll(variables, varIndex)
            {
                averages[varIndex][timeI] = dimensionedScalar
                (
                    averagesDict.lookup(variables[varIndex])
                ).value();
            }
        }
        else
        {
            forAll(variables, varIndex)
            {
                averages[varIndex][timeI] =
                    readScalar(averagesDict.lookup(variables[varIndex]));
            }
        }

        Info<< endl;
    }

    fileName outputFile(runTime.path()/"compiled_" + dictName);
    OFstream os(outputFile);
    
    os << "Time" << " " << variables[0];
    for (label varIndex = 1; varIndex < variables.size(); varIndex++)
    {
        os << " " << variables[varIndex];
    }

    os << endl;

    forAll(timeDirs, timeI)
    {
        os << timeDirs[timeI].name() << " " << averages[0][timeI];
        for (label varIndex = 1; varIndex < variables.size(); varIndex++)
        {
            os << " " << averages[varIndex][timeI];
        }
        os << endl;
    }
    Info<< "End\n" << endl;

    return 0;
}

// ************************************************************************* //
