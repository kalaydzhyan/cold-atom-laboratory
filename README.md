# Simulations of Bose-Einstein Condensates in microgravity

For description  of the problem and preliminary results read **writeup.pdf**

## latekick.m + earlykick.m

Matlab routine for simulation of BEC released from a magnetic trap. Delta-kick is applied at 100 ms (late) and 10 ms (early). In the current version of the code one can choose arbitrary frequencies of the trap, as well as cubic anharmonicities.
Modified GPELab toolbox is used to solve Gross-Pitaevskii equation numerically.

## adiabatic.m

Simulation of the adiabatic cooling in harmonic approximation. Initial trap frequency 1kHz, going down to about 1Hz.

## Mathematica files

**CALdkcZanh.nb**: Simulation of the delta-kick cooling for the thermal cloud trapped by the CAL chip. Z-wire trap used.

**CALdkcdimpleharm.nb**: Simulation of the delta-kick cooling for the thermal cloud and BEC in harmonic approximation of the trap potential. Z-wire + dimple wire.

**adiabatic.nb**: Adiabatic cooling applied to the condensate in harmonic approximation of the cooling potential. Preliminary calculations for the thermal cloud.
