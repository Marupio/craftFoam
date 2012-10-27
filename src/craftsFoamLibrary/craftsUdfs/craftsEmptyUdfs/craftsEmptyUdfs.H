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
    Foam::craftsEmptyUdfs

Description
    A class containing fullCoupled function hooks:
    
    @verbatim
    
    * initializeTimestep (function hook)
        - Do nothing
    Field loop
    {
        Solve coupled reaction matrix
        * applyVariableLimits (function hook)
            - Use default
        * initializeImplicitLoop (function hook)
            - Do nothing
        Implicit loop
        {
            * implicitLoop (function hook)
            - Do nothing
        }
        * finalizeImplicitLoop (function hook)
            - Do nothing
    }
    * finalizeTimestep (function hook)
        - Do nothing
    
    @endverbatim

    If your function hooks object has its own data that must be recalled should
    the solver return to a previous timestep, you must define how to store and
    recover that data in the functions:
    
    @verbatim
    
    * saveState(slotNumber) - copy all necessary data into backup variables
    * clearState(slotNumber) - delete data from the backup
    * loadState(slotNumber) - copy the data back to the active variables
    * nStates() - return the number of slots available (-1 if it doesn't track)
    * validState(slotNumber) - true if slotNumber can be loaded
    
    @endverbatim

    There can be any number of saved states, indexed by slot number.  You must
    start at slot zero, and increase by one at a time to create a new slot.
    You do not need to save data that can be derived from standard/implicit/
    derived variables.
        
SourceFiles
    craftsEmptyUdfs.C

Author
    David L. F. Gaden

\*---------------------------------------------------------------------------*/

#ifndef craftsEmptyUdfs_H
#define craftsEmptyUdfs_H

#include "craftsUdfs.H"

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

namespace Foam
{

/*---------------------------------------------------------------------------*\
                      Class craftsEmptyUdfs Declaration
\*---------------------------------------------------------------------------*/

template<int matrixSize>
class craftsEmptyUdfs
:
    public craftsUdfs<matrixSize>
{

private:

protected:
    using craftsUdfs<matrixSize>::runTime;
    using craftsUdfs<matrixSize>::mesh;
    using craftsUdfs<matrixSize>::model;
    using craftsUdfs<matrixSize>::admReacs;
    using craftsUdfs<matrixSize>::admVars;
    using craftsUdfs<matrixSize>::admCoeffs;
    using craftsUdfs<matrixSize>::eqns;

public:

    //- Runtime type information
    TypeName("empty");

    // Constructors

        //- Construct from dictionary, given name
        craftsEmptyUdfs
        (
            craftsModel<matrixSize>& model,
            const word& hooksDictName
        );

    // Destructor

        virtual ~craftsEmptyUdfs()
        {}


    // Function hooks

        //- Executed at a new timestep, before the main solver loop begins
        virtual void initializeTimestep();
        
        //- Executed after solving the coupled reaction matrix, used to apply
        //  limits to the variables.  Leave this commented out to implement the
        //  default: a simple max / min, but this may interfere with mass
        //  conservation if variables stray too far out of limits.
        //virtual void applyVariableLimits();

        //- Executed before the implicit loop begins
        virtual void initializeImplicitLoop();
        
        //- Executed in the innermost loop, the implicit loop.  Return value
        //  indicates error level:
        //      0 = success
        //     -1 = total failure, retry at shorter timestep
        //     -2 = convergence failure, solve coupled matrix again
        //     -3 = stabilization failure, repeat implicitLoop
        virtual label implicitLoop();
        
        //- Executed after the implicit loop reaches solution
        virtual void finalizeImplicitLoop();

        //- Executed after a timestep reaches solution
        virtual void finalizeTimestep();

    // State control
    
        //- Save data to backup variables
        virtual void saveState(const label slot = 0);
        
        //- Erase backup variables
        virtual void clearState(const label slot = 0);
        
        //- Load backup variables into active variables
        virtual void loadState(const label slot = 0);
        
        //- Return number of spots that are occupied
        //  A value of -1 means it doesn't track nStates
        virtual label nStates() const;

        //- Returns true if this is a loadable save state slot
        virtual bool validState(const label slot) const;
};


// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

} // End namespace Foam

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

#ifdef NoRepository
#   include "craftsEmptyUdfs.C"
#endif

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

#endif

// ************************************************************************* //
