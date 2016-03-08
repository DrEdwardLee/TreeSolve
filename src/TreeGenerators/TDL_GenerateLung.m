function vessel_tree = TDL_GenerateLung(generations, parameters, angle1, angle2)
    % TDL_SolveVesselTree Generates an artifical vessel tree for use with TreeSolve
    %
    % The angle parameter allows you to specify an asymmetric tree
    %
    %     Author: Tom Doel www.tomdoel.com
    %     Part of TreeSolve. http://github.com/tomdoel/TreeSolve
    %     Distributed under the GNU GPL v3 licence. Please see LICENSE file.
    %    
    
    initial_radius = 3;
    initial_length = 20;
    max_tree_level = generations;

    length_modifier = 0.7;
    rad_modifier = 0.7;

    angle_change1 = angle1;
    angle_change2 = -angle2;

    first_coordinate = [0, 0, 0];

    start_point = first_coordinate;
    current_angle = pi/2;
    end_point = start_point + initial_length*[sin(current_angle), 0, cos(current_angle)];

    vessel_tree = NewVessel(initial_length, initial_radius, start_point, end_point, parameters);
    vessel_tree.angle2d = current_angle;

    if (generations > 0)
        vessel_tree = GenerateBranchesFor([1], 1, max_tree_level, 2, vessel_tree, parameters, ...
            length_modifier, rad_modifier, angle_change1, angle_change2);
    end

    % The angle2d field is only used for generating the tree
    vessel_tree = rmfield(vessel_tree, 'angle2d');
end



function vessel_tree = GenerateBranchesFor(...
    list_of_vessel_indices, tree_level, max_tree_level, next_vessel_index, ...
    vessel_tree, parameters, length_modifier, rad_modifier, angle_change1, angle_change2)
    next_set_of_indices = [];
    for vessel_index = list_of_vessel_indices
        vessel = vessel_tree(vessel_index);
        
        new_vessel_1_index = next_vessel_index;
        new_vessel_2_index = next_vessel_index + 1;
        next_vessel_index = next_vessel_index + 2;
        
        current_angle = vessel.angle2d;
        current_unstretched_radius = vessel.R_unstretched_average;

        next_radius = current_unstretched_radius*rad_modifier;
        next_length = vessel.length*length_modifier;
        next_start_point = vessel.end_point;
        next_angle_1 = current_angle + angle_change1;
        next_angle_2 = current_angle + angle_change2;
        next_end_point_1 = next_start_point + next_length*[sin(next_angle_1), 0, cos(next_angle_1)];
        next_end_point_2 = next_start_point + next_length*[sin(next_angle_2), 0, cos(next_angle_2)];
            
        child_1 = NewVessel(next_length, next_radius, next_start_point, next_end_point_1, parameters);
        child_1.angle2d = next_angle_1;
        vessel_tree(new_vessel_1_index) = child_1;
        
        child_2 = NewVessel(next_length, next_radius, next_start_point, next_end_point_2, parameters);
        child_2.angle2d = next_angle_2;
        vessel_tree(new_vessel_2_index) = child_2;

        vessel.connected_to = [new_vessel_1_index, new_vessel_2_index];
        vessel_tree(vessel_index) = vessel;
        
        next_set_of_indices = [next_set_of_indices, new_vessel_1_index, new_vessel_2_index];      %#ok<AGROW>
    end
    
    if (tree_level < max_tree_level)
        tree_level = tree_level + 1;
        vessel_tree = GenerateBranchesFor(next_set_of_indices, tree_level, max_tree_level, ...
            next_vessel_index, vessel_tree, parameters, length_modifier, rad_modifier, angle_change1, angle_change2);
    end
end
    
function vessel = NewVessel(length, R_unstretched, start_point, end_point, parameters)
    vessel = struct('length', length, 'start_point', start_point, 'end_point', end_point, 'R_unstretched_average', R_unstretched, 'connected_to', []);
end