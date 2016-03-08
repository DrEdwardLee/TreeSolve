function vessel_tree = TDL_ConvertVesselTreeFromAQ(vessel_tree_input, parameters)
    % TDL_ConvertVesselTreeFromAQ Helper function for TreeSolve
    %
    %
    %     Author: Tom Doel www.tomdoel.com
    %     Part of TreeSolve. http://github.com/tomdoel/TreeSolve
    %     Distributed under the GNU GPL v3 licence. Please see LICENSE file.
    % 
    
    vessel_tree = rmfield(vessel_tree_input, {'A', 'Q'});
    
    for vessel_index = 1 : length(vessel_tree)
       
       vessel_tree(vessel_index).R = sqrt(vessel_tree_input(vessel_index).A/pi); 
       vessel_tree(vessel_index).V = vessel_tree_input(vessel_index).Q./vessel_tree_input(vessel_index).A;
       R_unstretched = vessel_tree(vessel_index).R_unstretched;
       vessel_tree(vessel_index).p = CalculateFhat(vessel_tree_input(vessel_index).A, R_unstretched, parameters);
    end
