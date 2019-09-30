Liquid Rocket Engine Performance Code
Stanford Student Space Initiative
Palo Alto, California

Open Source as of 9/2019

----------------------------------------------
ABOUT
This code is meant to simulate the performance of small liquid rocket engines.
Currently, the code is set to handle the properties of nitrous oxide as an
oxidizer. However, the properties of any fuel can be taken as input. 

Features:
 - Two-phase nitrous oxide injector flow simulation using the Homogeneous 
   Equilibrium Model
 - Ability to import real combustion data varying with oxidizer to fuel 
   ratio and combustion chamber pressure.
 - Capability of simulating a self-pressurized oxidizer tank, an oxidizer
   tank with pressurant gas loaded in the ullage, or an oxidizer tank with an
   external pressurant gas supply tank with a regulator. 
 - Over fifteen plotting outputs allowing extensive insight into rocket 
   motor dynamics.

----------------------------------------------
MAIN SCRIPTS
RunPerformanceCode.m
This script defines the characteristics of the liquid rocket engine to be
simulated and runs the simulation and output plotting. Options include
plotting outputs with test data in a specific format. This can run in a mode
that simulations only outflow of oxidizer (cold flow) or the full operation
of the rocket engine with combustion (hot fire). 

----------------------------------------------
TEST CASES
Test cases are meant to easily verify that parts of the code are working for
development purposes. However, they may serve useful purposes.

PlotN2OProperties.m
This script tests the generation of nitrous oxide oxidizer properties. It
also provides insight into the important properties of the chemical that
drive the rocket motor dynamics, such as the vapor pressure and its
dependence on temperature.

TestCodeAccuracy.m
This script tests the numerical convergence of the simulation and can be used
to determine the numerical error as a function of time step. 

TestNozzleCalc.m
This script tests the compressible quasi-1D flow calculator that is used
for nozzle calculations and injector gas flow calculations. 

TestOxMassFlux.m
This script tests the calculation of two-phase nitrous oxide flow
calculation. The functions tested here are used for the prediction of
injector mass flow rate. 

----------------------------------------------
GENERAL USAGE NOTES
 - Scripts should always be run from the directory in which they are located.
 - Combustion data is derived from the program RPA Lite v1.2. This program
   generates a text file which can be parsed with the provided function
   CombustionDataProcess.m
 - Nitrous oxide properties are derived from data from the National Institute 
   of Standards and Technology's online chemical webbook.
 - The "Test Data" directory is meant to hold test data to plot against
   performance code data. This data can be referenced by 
   RunHybridPerformance.m
 - The RunHybridPerformance.m file can be used to store all information 
   about a rocket motor design.

----------------------------------------------
DEVELOPMENT NOTES
 - For ease of calculation, all units within the simulation are metric. This
   includes the following units: m, s, kg, K, mol, J, N. However, outputs or
   inputs may be defined in other convenient units. 
 - Output variable recording is set up to be independent of the state
   calculation for run-time efficiency. 

----------------------------------------------
THEORETICAL BASIS
Tank Dynamics:
 - The liquid and gas within the oxidizer tank is assumed to be in thermal
   equilibrium. The van der Waal's equation is used to model the state of the
   gas and while not completely as accurate as empirical saturation state data,
   allows the calculation of non-saturated states. 
Pressurant Gas:
 - Pressurant gas can be modeled as loaded into the ullage volume or supplied 
   from an external tank. In either case, all gas within the tank is assumed to
   be in thermal equilibrium with the nitrous oxide. Isentropic expansion is
   simulated in the external pressurant tank. 
Injector Dynamics:
 - Two types of injector flow are modeled. While liquid remains in the tank, 
   two-phase flow is assumed, modeled using the Homogeneous Equilibrium Model.
   With only gas in the oxidizer tank, isentropic compressible quasi-1D flow 
   modelling is used. However, no frictional losses are modeled apart from a
   discharge coefficient factor for the injector orifices. 
Combustion:
 - The combustion is modeled by predicting the exhaust properties within the
   combustion chamber. An efficiency factor, or c-star efficiency, is applied
   to the exhaust gas temperature. 
Nozzle:
 - Quasi-1D isentropic flow is used to model the nozzle flow, including the
   cases of subsonic flow, supersonic flow, and normal shocks in the nozzle. 
   A nozzle exhaust thermal efficiency is applied to the exhaust energy. A
   correction factor accounting for the divergence factor of the nozzle is
   also applied. 
Integration:
 - Euler's method is used for integration, using a constant time step.
