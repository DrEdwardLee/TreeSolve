function vessel_tree = TDL_SolvePoiseuilleFlow(vessel_tree, parameters)
    % TDL_SolvePoiseuilleFlow Runs a Poiseuille flow simulation of flow 
    %
    % This is a novel implementation adapted in part from Poiseuille flow
    % equations by Pedley et al 1970.
    %
    %
    %     Author: Tom Doel www.tomdoel.com
    %     Part of TreeSolve. http://github.com/tomdoel/TreeSolve
    %     Distributed under the GNU GPL v3 licence. Please see LICENSE file.
    %
    
    
    x_0 = InitialConditionsFromVesselTree(vessel_tree, parameters);
    if (imag(x_0) ~= zeros(size(x_0)))
        error('Imaginary input');
    end
    options = optimset('Display','notify', 'FunValCheck', 'on');
   
    [x_solved, fval, exitflag, output] = fsolve(@PoiseuilleSystem, x_0, options, vessel_tree, parameters);
    vessel_tree = StoreResultsInVesselTree(vessel_tree, x_solved, parameters);
end

function x_0 = InitialConditionsFromVesselTree(vessel_tree, parameters)
    number_of_vessels = length(vessel_tree);
    input_pressures = zeros(number_of_vessels, 1);
    output_pressures = zeros(number_of_vessels, 1);
    velocities = zeros(number_of_vessels, 1);

    for vessel_index = 1 : number_of_vessels
        vessel = vessel_tree(vessel_index);
        input_pressures(vessel_index) = vessel.p(1);
        output_pressures(vessel_index) = vessel.p(end);
        velocities(vessel_index) = mean(vessel.V);
    end
    x_0 = PutVariablesInVector(input_pressures, output_pressures, velocities);
end

function vessel_tree = StoreResultsInVesselTree(vessel_tree, x, parameters)
    p_vessels = GeneratePVessels(x, vessel_tree, parameters);
    
    for vessel_index = 1 : length(vessel_tree)
        vessel = vessel_tree(vessel_index);
        p_vessel = p_vessels(vessel_index);
        one_vector = ones(size(vessel.p));
        
        pressure_step = (p_vessel.output_pressure - p_vessel.input_pressure)/(length(vessel.p) - 1);
        pressure_vector = (p_vessel.input_pressure : pressure_step : p_vessel.output_pressure)';
        vessel_tree(vessel_index).p = pressure_vector;
        vessel_tree(vessel_index).R = TDL_CalculateRadiusFromPressure(pressure_vector, vessel.R_unstretched, parameters);
        vessel_tree(vessel_index).V = one_vector*p_vessel.velocity;
    end
end

function output_vector = PoiseuilleSystem(x, vessel_tree, parameters)
    
    output_vector = zeros(size(x));
    number_of_vessels = length(vessel_tree);
    p_vessels = GeneratePVessels(x, vessel_tree, parameters);
        
    output_index = 1;

    for vessel_index = 1 : number_of_vessels
        parent = p_vessels(vessel_index);

        % Poiseuille Flow Equations (one for each vessel)
        output_vector(output_index) = parent.input_pressure - parent.output_pressure - parent.flow*parent.resistance;
        
        output_index = output_index + 1;
        
        if (vessel_index == 1)
            % Boundary condition for inlet
            output_vector(output_index) = parameters.boundary_equation(100, parent.input_pressure, parent.flow, parameters, 1);
            output_index = output_index + 1;
        end
        
        child_vessels = vessel_tree(vessel_index).connected_to;
        if (isempty(child_vessels))
            % Boundary conditions for outlet
            
            output_vector(output_index) = parameters.boundary_equation(100, parent.output_pressure, parent.flow, parameters, -1);
            output_index = output_index + 1;
            
        else
            % Bifurcation equations: 3 for each bifurcation
            child_index_1 = child_vessels(1);
            child_index_2 = child_vessels(2);

            child_1 = p_vessels(child_index_1);
            child_2 = p_vessels(child_index_2);
            
            % conservation of flow
            output_vector(output_index) = parent.flow - child_1.flow - child_2.flow;
            output_index = output_index + 1;

           if (~parameters.momentum_conservation)
               % conservation of pressure
               output_vector(output_index) = child_1.input_pressure - parent.output_pressure;
               output_index = output_index + 1;
 
               output_vector(output_index) = child_2.input_pressure - parent.output_pressure;
               output_index = output_index + 1;
               
           else
               % conservation of momentum
                output_vector(output_index) = - parent.output_pressure*parent.radius^2 + child_1.input_pressure*child_1.radius^2*cos(child_1.angle) ...
                    + child_2.input_pressure*child_2.radius^2*cos(child_2.angle);
                output_index = output_index + 1;

                output_vector(output_index) = child_1.input_pressure*child_1.radius^2*sin(child_1.angle) + child_2.input_pressure*child_2.radius^2*sin(child_2.angle);
                output_index = output_index + 1;
           end

            
            
        end
    end
end

function p_vessels = GeneratePVessels(x, vessel_tree, parameters)
    [input_pressures output_pressures velocities] = GetVariablesFromVector(x);

    p_vessels = struct;
    
    for vessel_index = 1 : length(vessel_tree)
        vessel = vessel_tree(vessel_index);
        input_pressure = input_pressures(vessel_index);
        output_pressure = output_pressures(vessel_index);
        velocity = velocities(vessel_index);
            
        if (vessel_index == 1)
            p_vessels = NewPVessel(input_pressure, output_pressure, velocity, vessel, parameters);
        else
            p_vessels(vessel_index) = NewPVessel(input_pressure, output_pressure, velocity, vessel, parameters);            
        end
    end
end

function p_vessel = NewPVessel(input_pressure, output_pressure, velocity, vessel, parameters)

    R_unstretched = mean(vessel.R_unstretched); % Note we use mean
    
    mean_pressure = (input_pressure + output_pressure)/2;
    radius = TDL_CalculateRadiusFromPressure(mean_pressure, R_unstretched, parameters);
    flow = pi*radius^2 * velocity;
    
    resistance = PoiseuilleResistance(radius, velocity, vessel.length, parameters);
    
    p_vessel = struct('input_pressure', input_pressure, 'output_pressure', output_pressure, 'velocity', velocity, 'radius', radius, 'flow', flow, ...
        'length', vessel.length, 'resistance', resistance, 'angle', vessel.angle);
end

function Rr = PoiseuilleResistance(radius, velocity, vessel_length, parameters)
    mu = parameters.rho*parameters.nu; % Dynamic viscosity = kinematic viscosity x density
    d = 2*radius; % diameter
    Rp = 128*mu*vessel_length/(pi*d^4);
    Re = CalculateReynoldsNumber(radius, velocity, parameters);
    Rr = (Rp*parameters.C /(4*sqrt(2)))*(Re*d/vessel_length)^(1/2);
end

function Re = CalculateReynoldsNumber(radius, velocity, parameters)
    diameter = 2*radius;
    mean_velocity=abs(velocity);
    Re = mean_velocity.*diameter/parameters.nu;
end

function [input_pressures output_pressures velocities] = GetVariablesFromVector(x)
    number_of_vessels = length(x)/3;
    input_pressures = x(1 : number_of_vessels);
    output_pressures = x(number_of_vessels + 1 : 2*number_of_vessels);
    velocities = x(2*number_of_vessels + 1 : end);
end

function x = PutVariablesInVector(input_pressures, output_pressures, velocities)
    number_of_vessels = length(input_pressures);
    x = zeros(number_of_vessels*3, 1);
    x(1 : number_of_vessels) = input_pressures;
    x(number_of_vessels + 1 : 2*number_of_vessels) = output_pressures;
    x(2*number_of_vessels + 1 : end) = velocities;
end
