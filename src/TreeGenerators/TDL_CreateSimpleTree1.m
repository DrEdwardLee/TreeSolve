function vessel_tree = TDL_CreateSimpleTree1(parameters)
    % Generates a simple single-vessel tree
    %
    %
    %     Licence
    %     -------
    %     Part of the TD Pulmonary Toolkit. https://github.com/tomdoel/pulmonarytoolkit
    %     Author: Tom Doel, 2009 www.tomdoel.com
    %     Distributed under the GNU GPL v3 licence. Please see website for details.
    %    
    

    vessel_1.length = 50; %mm
    vessel_1.N = 1 + ceil(vessel_1.length/parameters.dx); % Number of grid points
    vessel_1.R_unstretched = 1; % Initial unstretched radius
    vessel_1.p = parameters.outlet_applied_pressure*ones(vessel_1.N, 1);
    vessel_1.R = 1 * ones(vessel_1.N, 1);
    vessel_1.V = zeros(vessel_1.N, 1);
    vessel_1.connected_to = [];
    
    vessel_tree(1) = vessel_1;
