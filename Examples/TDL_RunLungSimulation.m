% TDL_RunLungSimulation Script for running a simulation of blood flow
%
% You can choose how the simulation is run by varying the parameters below.
% 
%
%
%     Author: Tom Doel www.tomdoel.com
%     Part of TreeSolve. http://github.com/tomdoel/TreeSolve
%     Distributed under the GNU GPL v3 licence. Please see LICENSE file.
%   


% Clear everything
clear all; clc; close all;
TDL_AddPaths

%% 
clc; clear all;

parameters.dx = 2; % mm Space step
parameters.rho = 1.05*10^-6; % kg.mm^-3 density
parameters.nu = 3.2; % mm^2.s^-1 Kinematic viscosity
parameters.alpha = 1.1; % dimensionless, Flow profile parameter
parameters.G_0 = 5; % kPa
parameters.beta = 3.2; % dimensionless, Wall elasticity constant
parameters.dt = 0.0001; % s Time step

parameters.C = 1.85; % For poiseuille flow scheme
parameters.L_a = 10^-9; % Parameters for Nic Smith's blood flow solution scheme
parameters.L_b = 10^-9; % Parameters for Nic Smith's blood flow solution scheme
parameters.L_c = 10^-9; % Parameters for Nic Smith's blood flow solution scheme

% Governing Equations
parameters.solution_scheme = @TDL_SolvePoiseuilleFlow;
%parameters.solution_scheme = @TDL_SolveNavierStokes;

% For Navier Stokes scheme, equations to solve at bifurcation
%parameters.bifurcation_function = @TDL_SolvePressureBifurcation;
%parameters.bifurcation_function = @TDL_SolveWhiteleyBifurcation;
parameters.bifurcation_function = @TDL_SolveSimpleMomentumBifurcation;
parameters.use_newton_solver = false; % Only applies to JW scheme

% For Poiseuille Flow scheme, equations to solve at bifurcation
parameters.momentum_conservation = false; % True = use momentum instead of pressure conservation

% Boundary Conditions
%parameters.boundary_equation = @TDL_BoundaryEquationImpedance;
parameters.boundary_equation = @TDL_BoundaryEquationPressure;
%parameters.boundary_equation = @TDL_BoundaryEquationVariablePressure;

% Pressure boundary conditions
parameters.inlet_applied_pressure = 2; %2 used in Burrowes paper,  10.6 used in steady state %kPa
parameters.outlet_applied_pressure = 1.25; %1.25 used in Burrowes paper, 5.6 used in steady state; %kPa

% Impedance boundary conditions
parameters.inlet_applied_flow = 3500; % mm^3.s^-1
parameters.outlet_applied_impedance = 10^-3; % 1.44*10^-7; % N.s.mm^-5

% Gravity
%parameters.gravity_vector = [-1 0 0];
parameters.gravity_vector = [0 0 -1];
parameters.gravity = 9800; % mm.s^-1
%parameters.gravity = 0; % mm.s^-1

% Paramaters relating to plotting and solution
parameters.p_plot_range = [-1 5];
parameters.R_plot_range = [0 1.5];
parameters.V_plot_range = [-200 200];
parameters.max_time = 1;
parameters.draw_skip = 10; % How often to draw graphs during solve

% Generate lung model 
vessel_tree = TDL_GenerateLung(4, parameters, pi/6, pi/5);

% Solve
vessel_tree = TDL_SolveVesselTree(vessel_tree, parameters);

% Choose a name for the output files
results_folder_name = 'my_simulation';

% Save results
TDL_OutputToCmgui(fullfile('Results', results_folder_name), vessel_tree);
save(fullfile('Results', results_folder_name), 'vessel_tree');


