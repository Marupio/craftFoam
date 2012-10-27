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
    Foam::craftsUdfs

Description
    Abstract base class for custom function hooks in the anaerobic digester
    model (ADM1).  These function hooks are necessary when using implicit
    variables.  The solver layout is:
    
    @verbatim
    
    * initializeTimestep (function hook)
    Field loop
    {
        Solve coupled reaction matrix
        * applyVariableLimits (function hook)
        * initializeImplicitLoop (function hook)
        Implicit loop
        {
            * implicitLoop (function hook)
        }
        * finalizeImplicitLoop (function hook)
    }
    * finalizeTimestep (function hook)

    Loop descriptions:
        Field loop
        - This loop solves the full field and compares it with the previous
          iterations.
        Implicit loop
        - This loop performs all calculations for adm implicit variables.
    
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
    craftsUdfs.C
    newcraftsUdfs.C

Author
    David L. F. Gaden

\*---------------------------------------------------------------------------*/

#ifndef craftsUdfs_H
#define craftsUdfs_H

#include "autoPtr.H"
#include "admTime.H"
#include "scalarMatrices.H"
#include "fvMesh.H"

#include "VectorNFieldTypes.H"
#include "volVectorNFields.H"
#include "blockVectorNMatrices.H"

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

// This macro is currently obsolete
#define admTypeName(NameString)                                               \
static const char* typeName_()                                                \
    {                                                                         \
        std::ostringstream buf;                                               \
        buf << NameString "(" << matrixSize << ")";                           \
        return buf.str().c_str();                                             \
    }                                                                         \
    static const ::Foam::word typeName;                                       \
    static int debug;                                                         \
    virtual const word& type() const                                          \
    {                                                                         \
        return typeName;                                                      \
    }

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

namespace Foam
{

// Forward declaration of friend functions and operators

template<int matrixSize> class craftsModel;
class admVariableManager;
class admCoefficientManager;
class admReactionManager;

/*---------------------------------------------------------------------------*\
                        Class craftsUdfs Declaration
\*---------------------------------------------------------------------------*/

template<int matrixSize>
class craftsUdfs
{
protected:
    //- Reference to craftsModel - access to runTime, admVars, admCoeffs, etc..
    craftsModel<matrixSize>& model;

    //- Quick reference
    admTime& runTime;
    const fvMesh& mesh;
    admVariableManager& admVars;
    admCoefficientManager& admCoeffs;
    admReactionManager& admReacs;
    IOEquationReader& eqns;

    //- Dictionary for user-defined settings
    IOdictionary admHooksDict_;

public:

    //- Runtime type information

//#include "base.H"
    TypeName("craftsUdfs")

    // Declare run-time constructor selection tables

        declareRunTimeSelectionTable
        (
            autoPtr,
            craftsUdfs,
            fullCoupledRef,
            (
                craftsModel<matrixSize>& model,
                const word& hooksDictName
            ),
            (model, hooksDictName)
        );


    // Constructors

        //- Construct given craftsModel reference
        craftsUdfs
        (
            craftsModel<matrixSize>& model,
            const word& hooksDictName
        );

        //- Construct copy
        craftsUdfs(const craftsUdfs&);

    // Selectors


        //- Return a reference to a new craftsUdfs
        //  *** Located in newcraftsUdfs.C ***
        static autoPtr<craftsUdfs> New
        (
            word craftsUdfsName,
            craftsModel<matrixSize>& model,
            const word& hooksDictName
        );


    // Destructor

        virtual ~craftsUdfs()
        {}


    // Function hooks

        //- Executed at a new timestep, before the main solver loop begins
        virtual void initializeTimestep()
        {}

        //- Executed after solving the coupled reaction matrix, used to apply
        //  limits to the variables (default is simple max / min, but this may
        //  interfere with mass conservation if variables stray too far out of
        //  limits)
        virtual void applyVariableLimits();
        
        //- Executed before the implicit loop begins
        virtual void initializeImplicitLoop()
        {}
        
        //- Executed in the innermost loop, the implicit loop.  Return value
        //  indicates error level:
        //      0 = success
        //     -1 = total failure, retry at shorter timestep
        //     -2 = convergence failure, solve coupled matrix again
        //     -3 = stabilization failure, repeat implicitLoop
        virtual label implicitLoop()
        {
            return 0;
        }
        
        //- Executed after the implicit loop reaches solution
        virtual void finalizeImplicitLoop()
        {}

        //- Executed after a timestep reaches solution
        virtual void finalizeTimestep()
        {}
    
    // State control
    
        //- Save data to backup variables
        virtual void saveState(const label slot = 0)
        {}
        
        //- Erase backup variables
        virtual void clearState(const label slot = 0)
        {}
        
        //- Load backup variables into active variables
        virtual void loadState(const label slot = 0)
        {}
        
        //- Return number of spots that are occupied
        //  A value of -1 means it doesn't track nStates
        virtual label nStates() const
        {
            return -1;
        }

        //- Returns true if this is a loadable save state slot
        virtual bool validState(const label slot) const
        {
            return true;
        }
};


// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

} // End namespace Foam

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

#define MakeFullCoupledHooksBaseDefinitions(matrixSize)                       \
    defineTemplateTypeNameAndDebugWithName                                    \
    (                                                                         \
        craftsUdfs<matrixSize>,                                      \
        "craftsUdfs",                                                \
        0                                                                     \
    );                                                                        \
    defineTemplatedRunTimeSelectionTable                                      \
    (                                                                         \
        craftsUdfs, fullCoupledRef, matrixSize                       \
    );

#define MakeFullCoupledHooksDerivedDefinitions(matrixSize, lookupName, Name)  \
    defineTemplateTypeNameAndDebugWithName                                    \
    (                                                                         \
        Name<matrixSize>,                                                     \
        #Name,                                                                \
        0                                                                     \
    );                                                                        \
    addNamedTemplatedToRunTimeSelectionTable                                  \
    (                                                                         \
        craftsUdfs,                                                  \
        Name,                                                                 \
        matrixSize,                                                           \
        fullCoupledRef,                                                       \
        lookupName                                                            \
    );

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

#ifdef NoRepository
#   include "craftsUdfs.C"
#endif

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

#endif

// ************************************************************************* //
