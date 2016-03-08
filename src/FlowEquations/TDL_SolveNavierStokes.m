function vessel_tree = TDL_SolveNavierStokes(vessel_tree, parameters)
    % TDL_SolveNavierStokes Runs a Navier-Stokes simulation of flow using a
    % system of equations reduced to 1D
    %
    % This is a novel implementation adapted in part from a scheme by 
    % Burrowes at Tawhai, 2006, adapted from Smith et al., 2002.
    %
    %
    %     Author: Tom Doel www.tomdoel.com
    %     Part of TreeSolve. http://github.com/tomdoel/TreeSolve
    %     Distributed under the GNU GPL v3 licence. Please see LICENSE file.
    %    

    % Set up tree plotting (caching handles saves drawing time)
    handle_matrix = TDL_SetUpVesselPlot(length(vessel_tree), parameters);

    max_time = parameters.max_time;
    time = 0;
    skip = 0;
    steady_state_count = 0;
    max_steady_state = 100;
    
    % We don't look for a steady state at each step, to speed up the solution
    matches_skipped = 0;    
    skip_before_matching = 1000;

    steady_state_to_match = vessel_tree;

    while time < max_time
        vessel_tree = NextTimeStep(vessel_tree, parameters, time);
        matches_skipped = matches_skipped + 1;
        if (matches_skipped > skip_before_matching)
            closeness = VesselTreeValuesMatch(steady_state_to_match, vessel_tree);
            if (closeness < 0.001)
                steady_state_count = steady_state_count + 1;
                if (steady_state_count > max_steady_state)
                    disp('Steady state reached');

                     % Plot the final results and fit to axes
                     for vessel_number = 1:length(vessel_tree)
                         TDL_PlotVessel(vessel_tree, vessel_number, parameters, 'b', handle_matrix, true);
                      end

                    return;
                end
            else
                steady_state_count = 0;
                steady_state_to_match = vessel_tree;
                matches_skipped = 0;
            end
        end

        if (mod(skip,parameters.draw_skip) == 0)
            TDL_PlotVesselTree(vessel_tree, parameters, 'b', handle_matrix, false);
            pause(0.001);
        end
        skip = skip + 1;
        time = time + parameters.dt;

        % Show current time
        %time

    end
    
    error('Maximum time exceeded without steady state');


function closeness = VesselTreeValuesMatch(tree_1, tree_2)
    for (vessel_index = 1 : length(tree_1))
        norm_p = norm(tree_1(vessel_index).p - tree_2(vessel_index).p);
        norm_R = norm(tree_1(vessel_index).R - tree_2(vessel_index).R);
        norm_V = norm(tree_1(vessel_index).V - tree_2(vessel_index).V);
    end
    closeness = max(max(norm_p, norm_R), norm_V);

    
function vessel_tree_t1 = NextTimeStep(vessel_tree_t0, parameters, time)

    vessel_tree = vessel_tree_t0;
    num_vessels = length(vessel_tree_t0);
    
    
    % Calculate intermediate and internal points for all vessels; these are required before
    % we solve the equations at the bifurcations
    for vessel_index = 1 : num_vessels
        vessel = vessel_tree(vessel_index);
        R_unstretched = vessel.R_unstretched;
        angle = vessel.angle_to_vertical;

        [iV_half, iR_half, ip_half, ipint_t1, iVint_t1, iRint_t1] = ...
            TDL_LWStep(vessel.V, vessel.R, vessel.p, R_unstretched, angle, parameters);
       
        % These are values at points t+dt/2, x+dx/2
        vessel_tree(vessel_index).V_half = iV_half;
        vessel_tree(vessel_index).R_half = iR_half;
        vessel_tree(vessel_index).p_half = ip_half;
        
        % The values at t+dt for internal points in the vessels
        vessel_tree(vessel_index).V_int_t1 = iVint_t1;
        vessel_tree(vessel_index).R_int_t1 = iRint_t1;
        vessel_tree(vessel_index).p_int_t1 = ipint_t1;       
    end

    
    % Now determine the end points of each vessel; these are either determined
    % by boundary conditions or by solving the conservation equations
    % across bifurcations
    for vessel_index = 1 : num_vessels
        vessel = vessel_tree(vessel_index);
         
        % Assume the first vessel is the inlet so impose boundary conditions
        if (vessel_index == 1)
            vessel_tree(vessel_index).firstpoint = VesselBoundaryStart(vessel, time, parameters);
        end
        
        % Now examine whether there are any vessels branching from the end
        % of this vessel
        if (isempty(vessel.connected_to))
            % No vessels branch from this one, so impose an outlet boundary
            % condition
            vessel_tree(vessel_index).lastpoint = VesselBoundaryEnd(vessel, time, parameters);

        else
            % Vessels branch from this one, so solve for the 3 vessel endpoints at the bifurcation
            child_index_1 = vessel.connected_to(1);
            child_index_2 = vessel.connected_to(2);

            % Calculate points at bifurcations
            [parent_t1, child_1_t1, child_2_t1] = parameters.bifurcation_function(vessel_tree(vessel_index), ...
                vessel_tree(child_index_1), vessel_tree(child_index_2), parameters);
    
            % Update boundary points at bifurcation
            vessel_tree(vessel_index).lastpoint = parent_t1;
            vessel_tree(child_index_1).firstpoint = child_1_t1;
            vessel_tree(child_index_2).firstpoint = child_2_t1;
        end
    end


    % Update vessel tree with the new values calculated above
    vessel_tree_t1 = vessel_tree_t0;
    for vessel_index = 1:num_vessels
        % Update the internal points
        vessel_tree_t1(vessel_index).p(2:end-1) = vessel_tree(vessel_index).p_int_t1;
        vessel_tree_t1(vessel_index).R(2:end-1) = vessel_tree(vessel_index).R_int_t1;
        vessel_tree_t1(vessel_index).V(2:end-1) = vessel_tree(vessel_index).V_int_t1;
        
        % Update the boundary points
        vessel_tree_t1(vessel_index).p(1) = vessel_tree(vessel_index).firstpoint.p;
        vessel_tree_t1(vessel_index).R(1) = vessel_tree(vessel_index).firstpoint.R;
        vessel_tree_t1(vessel_index).V(1) = vessel_tree(vessel_index).firstpoint.V;
    
        vessel_tree_t1(vessel_index).p(end) = vessel_tree(vessel_index).lastpoint.p;
        vessel_tree_t1(vessel_index).R(end) = vessel_tree(vessel_index).lastpoint.R;
        vessel_tree_t1(vessel_index).V(end) = vessel_tree(vessel_index).lastpoint.V;
    end
    
    
    
% Boundary conditions
function point_t1 = VesselBoundaryStart(vessel_t0, time, parameters)
    options = optimset('Display', 'off');
    [p1_t1, fval, exitflag, output] = fsolve(@VesselBoundaryStartForPressureEquation, vessel_t0.p(1), options, vessel_t0, time, parameters);
    if (exitflag < 1)
        error('TDL_SolveNavierStokes:VesselBoundaryStart: Solver Failed');
    end
    point_t1 = TDL_VesselBoundaryStartFromPressure(p1_t1, vessel_t0, parameters);

% Boundary Conditions
function point_t1 = VesselBoundaryEnd(vessel_t0, time, parameters)    
    options = optimset('Display', 'off');
    [pN_t1, fval, exitflag, output] = fsolve(@VesselBoundaryEndForPressureEquation, vessel_t0.p(end), options, vessel_t0, time, parameters);
    if (exitflag < 1)
        error('TDL_SolveNavierStokes:VesselBoundaryEnd: Solver Failed');
    end

    point_t1 = TDL_VesselBoundaryEndFromPressure(pN_t1, vessel_t0, parameters);

    
function remainder = VesselBoundaryStartForPressureEquation(pressure, vessel_t0, time, parameters)
    point_t1 = TDL_VesselBoundaryStartFromPressure(pressure, vessel_t0, parameters);
    flow = pi*point_t1.V*(point_t1.R)^2;
    remainder = parameters.boundary_equation(time, pressure, flow, parameters, 1);
      
function remainder = VesselBoundaryEndForPressureEquation(pressure, vessel_t0, time, parameters)
    point_t1 = TDL_VesselBoundaryEndFromPressure(pressure, vessel_t0, parameters);
    flow = pi*point_t1.V*(point_t1.R)^2;
    remainder = parameters.boundary_equation(time, pressure, flow, parameters, -1);
    