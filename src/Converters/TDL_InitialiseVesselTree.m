function vessel_tree = TDL_InitialiseVesselTree(vessel_tree, parameters)
    % TDL_InitialiseVesselTree Helper function for TreeSolve
    %
    %
    %     Author: Tom Doel www.tomdoel.com
    %     Part of TreeSolve. http://github.com/tomdoel/TreeSolve
    %     Distributed under the GNU GPL v3 licence. Please see LICENSE file.
    % 
    
    % Intial conditions
    for vessel_index = 1 : length(vessel_tree)
       vessel_tree(vessel_index).dx = parameters.dx; % spatial step
       vessel_tree(vessel_index).N = 1 + ceil(vessel_tree(vessel_index).length/vessel_tree(vessel_index).dx); % Number of grid points
       vessel_tree(vessel_index).R_unstretched = ones(vessel_tree(vessel_index).N, 1)*vessel_tree(vessel_index).R_unstretched_average;
       vessel_tree(vessel_index).p = zeros(vessel_tree(vessel_index).N, 1);
       vessel_tree(vessel_index).R = TDL_CalculateRadiusFromPressure(vessel_tree(vessel_index).p, ...
           vessel_tree(vessel_index).R_unstretched, parameters);
       vessel_tree(vessel_index).V = zeros(vessel_tree(vessel_index).N, 1);   
    end

    
    % Calculate angles - this is a bit of a hack
    vessel_tree(1).angle = 0;
    for element_index = 1 : length(vessel_tree)
        this_vessel = vessel_tree(element_index);
        this_element_vector = this_vessel.end_point - this_vessel.start_point;
        
        angle_to_vertical = TDL_CalculateAngleBetweenVectors(this_element_vector, parameters.gravity_vector);
        vessel_tree(element_index).angle_to_vertical = angle_to_vertical;
        
        child_vessels = this_vessel.connected_to;
        if (~isempty(child_vessels))
            child_1 = child_vessels(1); 
            child_2 = child_vessels(2); 
            child_vessel_1 = vessel_tree(child_1);
            child_vessel_2 = vessel_tree(child_2);
            child1_vector = child_vessel_1.end_point - child_vessel_1.start_point;
            child2_vector = child_vessel_2.end_point - child_vessel_2.start_point;
            angle_1 = TDL_CalculateAngleBetweenVectors(this_element_vector, child1_vector);
            angle_2 = TDL_CalculateAngleBetweenVectors(this_element_vector, child2_vector);
            angle_2 = - angle_2;
            vessel_tree(child_1).angle = angle_1;
            vessel_tree(child_2).angle = angle_2;
        end
    end
end


