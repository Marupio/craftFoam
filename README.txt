Source code from the main solver from my PhD work:

A general coupled reacting flow solver for incompressible flows with user-
defined implicit routines.  Solves partial differential algebraic equations
(PDAE) using a point-implicit coupled PDE solver, and user-defined implicit
routines, which are automatically coupled back into the PDE solver through
its source term.

This application depends on other repositories:
* equationReader
* plcEmulator
* multiSolver
* reactionReader
