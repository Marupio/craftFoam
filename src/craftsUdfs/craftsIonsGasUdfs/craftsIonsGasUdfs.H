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
    Foam::craftsIonsGasUdfs

Description
    A class containing fullCoupled function hooks for the DAE implementation
    of ADM1.  The function hooks are:
    
    @verbatim
    
    * initializeTimestep (function hook)
        - Do nothing
    Field loop
    {
        Solve coupled reaction matrix
        * applyVariableLimits (function hook)
            - Use proportional mass-conserving correction
        * initializeImplicitLoop (function hook)
            - Do nothing
        Implicit loop
        {
            * implicitLoop (function hook)
                - admGas_.implicitLoop - solves the gas model
                - getIons - Newton Raphson ion model solver
                - autoSolveAll - solve S_h2
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
    craftsIonsGasUdfs.C

Author
    David L. F. Gaden

\*---------------------------------------------------------------------------*/

#ifndef craftsIonsGasUdfs_H
#define craftsIonsGasUdfs_H

#include "craftsGasModel.H"
#include "craftsUdfs.H"

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

namespace Foam
{

/*---------------------------------------------------------------------------*\
                      Class craftsIonsGasUdfs Declaration
\*---------------------------------------------------------------------------*/

template<int matrixSize>
class craftsIonsGasUdfs
:
    public craftsUdfs<matrixSize>
{
    // Private member data

        //- A simple diffusion / bubble gas mass transfer model
        craftsGasModel<matrixSize> admGas_;

        // Output flags
        // These are yes/no settings that change what is output to the console
        // during a run.  They are contained in admSettingsDict.  Other classes
        // can add to this subdictionary.

            //- Ion solver performance (default false)
            bool outputIonSolverPerformance_;
            
            //- Function hooks summary (default false)
            bool outputFunctionHooksSummary_;


protected:
    using craftsUdfs<matrixSize>::runTime;
    using craftsUdfs<matrixSize>::mesh;
    using craftsUdfs<matrixSize>::model;
    using craftsUdfs<matrixSize>::admReacs;
    using craftsUdfs<matrixSize>::admVars;
    using craftsUdfs<matrixSize>::admCoeffs;
    using craftsUdfs<matrixSize>::eqns;

private:
    // Private member functions

        //- Use Newton-Raphson linearization to find ion implicit values
        // Returns error level, same as for implicitLoop()
        label getIons();

public:

    //- Runtime type information
    TypeName("sh2IonsGas");

    // Constructors

        //- Construct from dictionary, given name
        craftsIonsGasUdfs
        (
            craftsModel<matrixSize>& model,
            const word& hooksDictName
        );

    // Destructor

        virtual ~craftsIonsGasUdfs()
        {}

    // Function hooks

        //- Executed at a new timestep, before the main solver loop begins
        virtual void initializeTimestep();
        
        //- Executed after solving the coupled reaction matrix, used to apply
        //  limits to the variables.  Uses a mass conserving limiter.  For gas
        //  transfer variables, the cells at the transfer surface do not get
        //  the mass conserving correction applied, but they do get limited.
        virtual void applyVariableLimits();

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
#   include "craftsIonsGasUdfs.C"
#endif

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

#endif

// ************************************************************************* //
