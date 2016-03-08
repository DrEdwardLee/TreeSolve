function vessel_tree = TDL_CreateSimpleTree(parameters)
    % Generates a simple vessel tree representing a simple bifurcation
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
    vessel_1.p = 1.2*ones(vessel_1.N, 1);
    vessel_1.R = 1 * ones(vessel_1.N, 1);
    vessel_1.V = ones(vessel_1.N, 1);
    vessel_1.connected_to = [2 3];
    vessel_1.alpha_2 = pi/4;
    vessel_1.alpha_3 = pi/4;
    
    vessel_2.length = 50; %mm
    vessel_2.N = 1 + ceil(vessel_2.length/parameters.dx); % Number of grid points
    vessel_2.R_unstretched = 0.7; % Initial unstretched radius
    vessel_2.p = 1.2*ones(vessel_2.N, 1);
    vessel_2.R = 0.7 * ones(vessel_2.N, 1);
    vessel_2.V = ones(vessel_2.N, 1);
    vessel_2.connected_to = [];
    vessel_2.alpha_2 = pi/4;
    vessel_2.alpha_3 = pi/4;

    vessel_3.length = 50; %mm
    vessel_3.N = 1 + ceil(vessel_3.length/parameters.dx); % Number of grid points
    vessel_3.R_unstretched = 0.7; % Initial unstretched radius
    vessel_3.p = 1.2*ones(vessel_3.N, 1);
    vessel_3.R = 0.7 * ones(vessel_3.N, 1);
    vessel_3.V = ones(vessel_3.N, 1);
    vessel_3.connected_to = [];
    vessel_3.alpha_2 = pi/4;
    vessel_3.alpha_3 = pi/4;

    vessel_tree(1) = vessel_1;
    vessel_tree(2) = vessel_2;
    vessel_tree(3) = vessel_3;