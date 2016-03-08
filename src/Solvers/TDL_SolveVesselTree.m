function vessel_tree = TDL_SolveVesselTree(vessel_tree, parameters)
    % TDL_SolveVesselTree Performs a fluid flow simulation
    %
    %     Author: Tom Doel www.tomdoel.com
    %     Part of TreeSolve. http://github.com/tomdoel/TreeSolve
    %     Distributed under the GNU GPL v3 licence. Please see LICENSE file.
    %    
    
    % Set up angles etc.
    vessel_tree = TDL_InitialiseVesselTree(vessel_tree, parameters);
    TDL_DrawVesselTree(vessel_tree);

    vessel_tree = parameters.solution_scheme(vessel_tree, parameters);
    
    % Set up tree plotting (caching hanles saves drawing time)
    handle_matrix = TDL_SetUpVesselPlot(length(vessel_tree), parameters);

    % Plot the final results and fit to axes
    TDL_PlotVesselTree(vessel_tree, parameters, 'b', handle_matrix, true);
end