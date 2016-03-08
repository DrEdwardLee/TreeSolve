function vessel_tree = TDL_CreateTreeFromParameterVectors(lengths, radii, angles, parameters)
    % TDL_CreateTreeFromParameterVectors Generates a vessel tree using the supplied parameters
    %
    %
    %     Licence
    %     -------
    %     Part of the TD Pulmonary Toolkit. https://github.com/tomdoel/pulmonarytoolkit
    %     Author: Tom Doel, 2009 www.tomdoel.com
    %     Distributed under the GNU GPL v3 licence. Please see website for details.
    %    
    
    num_vessels = length(lengths);
    vessel_tree = [];

    for vessel_index = 1 : num_vessels
        vessel_tree(vessel_index).length = lengths(vessel_index);
        vessel_tree(vessel_index).R_unstretched = radii(vessel_index);
        vessel_tree(vessel_index).R_unstretched_average = radii(vessel_index);
        vessel_tree(vessel_index).angle = angles(vessel_index);
        vessel_tree(vessel_index).N = 1 + ceil(vessel_tree(vessel_index).length/parameters.dx); % Number of grid points
        vessel_tree(vessel_index).connected_to = [];
    end

    vessel_tree(1).start_point = [0 0 0];
    vessel_tree(1).end_point = vessel_tree(1).start_point + vessel_tree(1).length*[sin(0), 0, cos(0)];
    
    switch num_vessels
        case 1
            vessel_tree;
        case 3
            vessel_tree(1).connected_to = [2 3];
            
            vessel_tree(2).start_point = vessel_tree(1).end_point;
            vessel_tree(2).end_point = vessel_tree(2).start_point + vessel_tree(2).length*[sin(vessel_tree(2).angle), 0, cos(vessel_tree(2).angle)];
            vessel_tree(3).start_point = vessel_tree(1).end_point;
            vessel_tree(3).end_point = vessel_tree(3).start_point + vessel_tree(3).length*[sin(vessel_tree(3).angle), 0, cos(vessel_tree(3).angle)];
        case 7
            vessel_tree(1).connected_to = [2 3];
            vessel_tree(2).connected_to = [4 5];
            vessel_tree(3).connected_to = [6 7];
        otherwise
            error('unknown number of vessels');
    end
    
    vessel_tree = TDL_InitialiseVesselTree(vessel_tree, parameters);
