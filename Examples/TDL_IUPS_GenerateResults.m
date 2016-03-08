% TDL_IUPS_GenerateResults
%
% This example script runs the simulations used to generate the results
% presented at IUPS 2009 by Doel et al:
%
%   Assessing the relative influence of vascular branching structure on pulmonary blood flow distribution via computational models
%   T Doel, KS Burrowes, JP Whiteley, V Grau, DJ Gavaghan
%   IUPS, Journal of Physiological Sciences 59, 224 (2009)
%
% After generating the results by running this script, you can generate the
% graphs by running the following script:
%   TDL_IUPS_GenerateGraphs.m
%
%     Author: Tom Doel www.tomdoel.com
%     Part of TreeSolve. http://github.com/tomdoel/TreeSolve
%     Distributed under the GNU GPL v3 licence. Please see LICENSE file.
%    



% Clear everything
clear all; clc; close all;
TDL_AddPaths

parameters.dx = 1; % mm Space step
parameters.rho = 1.05*10^-6; % kg.mm^-3 density
parameters.nu = 3.2; % mm^2.s^-1 Kinematic viscosity
parameters.alpha = 1.1; % dimensionless, Flow profile parameter
parameters.G_0 = 5; % kPa
parameters.beta = 3.2; % dimensionless, Wall elasticity constant
parameters.dt = 0.0001; % s Time step

parameters.C = 1.85; % For poiseuille flow scheme
parameters.L_a = 10^-9; % Parameters for Navier-Stokes blood flow solution scheme
parameters.L_b = 10^-9; % Parameters for Navier-Stokes blood flow solution scheme
parameters.L_c = 10^-9; % Parameters for Navier-Stokes blood flow solution scheme

% For Navier Stokes scheme, equations to solve at bifurcation
parameters.use_newton_solver = false; % Only applies to flow conservation solver

% Pressure boundary conditions
parameters.inlet_applied_pressure = 2; % kPa
parameters.outlet_applied_pressure = 1.3; % kPa

% Impedance boundary conditions
parameters.inlet_applied_flow = 12000; % mm^3.s^-1
parameters.outlet_applied_impedance = 10^-6; % N.s.mm^-5

% Gravity
parameters.gravity_vector = [0 0 -1];

% Paramaters relating to plotting and solution
parameters.p_plot_range = [-1 5];
parameters.R_plot_range = [0 1.5];
parameters.V_plot_range = [-200 200];
parameters.max_time = 1;
parameters.draw_skip = 10; % How often to draw graphs during solve



% Generate lung model 
vessel_tree_initial = TDL_GenerateLung(4, parameters, pi/5, pi/8);

% Now we run the simulation a number of times, with different flow
% equations, conservation laws, boundary conditions and the presence of
% gravity

disp('*** Poiseuille Flow, gravity, pressure bifurcations, pressure boundaries ***');
parameters.solution_scheme = @TDL_SolvePoiseuilleFlow;
parameters.boundary_equation = @TDL_BoundaryEquationPressure;
parameters.gravity = 9800; % mm.s^-1
parameters.bifurcation_function = @TDL_SolvePressureBifurcation;
parameters.momentum_conservation = false; % True = use momentum instead of pressure conservation
vessel_tree = TDL_SolveVesselTree(vessel_tree_initial, parameters);
TDL_SaveVesselResults('iups_pois_grav_prbif_prbcs_', vessel_tree, parameters);

disp('*** Poiseuille Flow, gravity, pressure bifurcations, impedance boundaries ***'); 
parameters.solution_scheme = @TDL_SolvePoiseuilleFlow;
parameters.boundary_equation = @TDL_BoundaryEquationImpedance;
parameters.gravity = 9800; % mm.s^-1
parameters.bifurcation_function = @TDL_SolvePressureBifurcation;
parameters.momentum_conservation = false; % True = use momentum instead of pressure conservation
vessel_tree = TDL_SolveVesselTree(vessel_tree_initial, parameters);
TDL_SaveVesselResults('iups_pois_grav_prbif_impbcs_', vessel_tree, parameters);

disp('*** Poiseuille Flow, gravity, momentum bifurcations, impedance boundaries ***');
parameters.solution_scheme = @TDL_SolvePoiseuilleFlow;
parameters.boundary_equation = @TDL_BoundaryEquationImpedance;
parameters.gravity = 9800; % mm.s^-1
parameters.bifurcation_function = @TDL_SolveSimpleMomentumBifurcation;
parameters.momentum_conservation = true; % True = use momentum instead of pressure conservation
vessel_tree = TDL_SolveVesselTree(vessel_tree_initial, parameters);
TDL_SaveVesselResults('iups_pois_grav_mombif_impbcs_', vessel_tree, parameters);

disp('***  Navier Stokes, gravity, pressure bifurcations, impedance boundaries ***');
parameters.solution_scheme = @TDL_SolveNavierStokes;
parameters.boundary_equation = @TDL_BoundaryEquationImpedance;
parameters.gravity = 9800; % mm.s^-1
parameters.bifurcation_function = @TDL_SolvePressureBifurcation;
parameters.momentum_conservation = false; % True = use momentum instead of pressure conservation
vessel_tree = TDL_SolveVesselTree(vessel_tree_initial, parameters);
TDL_SaveVesselResults('iups_navs_grav_prbif_impbcs_', vessel_tree, parameters);

disp('***  Navier Stokes no gravity, pressure bifircations, impedance boundaries ***');
parameters.solution_scheme = @TDL_SolveNavierStokes;
parameters.boundary_equation = @TDL_BoundaryEquationImpedance;
parameters.gravity = 0; % mm.s^-1
parameters.bifurcation_function = @TDL_SolvePressureBifurcation;
parameters.momentum_conservation = false; % True = use momentum instead of pressure conservation
vessel_tree = TDL_SolveVesselTree(vessel_tree_initial, parameters);
TDL_SaveVesselResults('iups_navs_ng_prbif_impbcs_', vessel_tree, parameters);

TDL_IUPS_GenerateGraphs;
