function vessel_tree = TDL_VesselTreeFromMesh(element_nodes, node_coordinates, rad_coordinates, parameters)
    % TDL_VesselTreeFromMesh Generates an TreeSolve artifical tree using
    % the supplied parameters
    %
    %     Author: Tom Doel www.tomdoel.com
    %     Part of TreeSolve. http://github.com/tomdoel/TreeSolve
    %     Distributed under the GNU GPL v3 licence. Please see LICENSE file.
    %    

    number_of_elements = size(element_nodes, 1);
    number_of_nodes = size(node_coordinates, 1);
    
    vessel_tree = [];
    
    % Determine connectivity between elements
    nodes_to_elements = -1 * ones(number_of_nodes, 2);
    for element_index = 1 : number_of_elements
        node_1 = element_nodes(element_index, 1);
        if (nodes_to_elements(node_1, 1) == -1)
            nodes_to_elements(node_1, 1) = element_index;
        elseif (nodes_to_elements(node_1, 2) == -1)
            nodes_to_elements(node_1, 2) = element_index;
        else
           error('Unexpected number of elements for this node'); 
        end
    end
    
    for element_index = 1 : number_of_elements
        node_1 = element_nodes(element_index, 1);
        node_2 = element_nodes(element_index, 2);
        node_1_coords = node_coordinates(node_1, :);
        node_2_coords = node_coordinates(node_2, :);
        
        vessel_length = norm(node_1_coords - node_2_coords);
        vessel_R_unstretched = rad_coordinates{node_1}(1);
            
        new_vessel = NewVessel(vessel_length, vessel_R_unstretched, node_1_coords, node_2_coords, parameters);

        vessels_connected_to = nodes_to_elements(node_2, :);
        if (vessels_connected_to == [-1 -1])
            vessels_connected_to = [];
        end
        new_vessel.connected_to = vessels_connected_to;

        if (element_index == 1)
            vessel_tree = new_vessel;
        else
            vessel_tree(element_index) = new_vessel;
        end
    end
end

function vessel = NewVessel(length, R_unstretched, start_point, end_point, parameters)
    vessel = struct('length', length, 'start_point', start_point, 'end_point', end_point, 'R_unstretched_average', R_unstretched, 'connected_to', []);
end

