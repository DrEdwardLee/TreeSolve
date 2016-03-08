function vessel_tree = TDL_ConvertVesselTreeToAQ(vessel_tree_input)
    % TDL_ConvertVesselTreeToAQ Helper function for TreeSolve
    %
    %
    %     Author: Tom Doel www.tomdoel.com
    %     Part of TreeSolve. http://github.com/tomdoel/TreeSolve
    %     Distributed under the GNU GPL v3 licence. Please see LICENSE file.
    % 

    vessel_tree = rmfield(vessel_tree_input, {'p', 'R', 'V'});
    
    for vessel_index = 1 : length(vessel_tree)
       vessel_tree(vessel_index).A = pi*((vessel_tree_input(vessel_index).R).^2); 
       vessel_tree(vessel_index).Q = vessel_tree(vessel_index).A.*vessel_tree_input(vessel_index).V;         
    end
