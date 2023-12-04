# NSR
This file provides a description of the columns in CutoffWall.csv and Dam.csv as part of "Shoreline barriers may amplify coastal groundwater hazards with sea-level rise" by Su et al. submitted to Nature Scientific Report.

Model parameters
=====================
Model ID: A unique model identifier run for the various scenarios (e.g., 0-755)

Hydraulic Conductivity (K): Homogeneous and isotropic hydraulic conductivity with units of meters per day (m/day) used in the model ranging from 0.01 ~ 10 m/day (0.01, 0.1, 1, 10 m/day)

Surface Slope (S): topographic slope in non-dimensional units (e.g., meters per meter) ranging from 0.005 to 0.01 (0.005, 0.007, 0.01)

Barrier Depth: number indicating the count of layers opening, with total layers of 20. Each layer depth was set to approximately 2.5 m, vary slightly depends on surface topography. Also represents the opening size ranging from 0.1 to 0.9 (see Equation 1)

D*: relative opening calculated with barrier depth (D_opening)and aquifer depth (D_aquifer) using Equation 1.

Barrier Location: location of the barrier from the current sea level, moving inland. Ranging from 0 to 300 m inland, represent the relative location of 0 to 0.3.

L*: relative location calculated with the barrier construction location at current sea level and the entire distance of the disigned aquifer using Equation 2. 

Model results
=======================
R_intrusionUZF: UZF simulations, unitless ratio describing amount of protection provided by the barrier relative to a simulation with no barrier, (w/o barrier saltwedge length - w barrier saltwedge length) / w/o barrier saltwedge length (see Equation 3)

R_intrusionCHD: CHD simulations, unitless ratio describing amount of protection providedby the barrier relative to a simulation with no barrier, (w/o barrier saltwedge length - w barrier saltwedge length) / w/o barrier saltwedge length (see Equation 3)

R_flowUZF: UZF simulations, unitless ratio describing the proportion of groundwater flow that discharges inland of the barrier, total flow at barrier to inland aquifer/ total UZF recharge (see Equation 4)

R_flowCHD: CHD simulations, unitless ratio describing the proportion of groundwater flow that discharges inland of the barrier, total flow at barrier to inland aquifer/ total UZF recharge (see Equation 4)

R_emergence: unitless ratio describing the amount of land inland of the barrier with groundwater emergence, calculated from water table and model topography. Counted emergence cell numebr / inland cell number (see Equation 5)
