function vessel_tree = TDL_CreateSimpleTree3(parameters)
    % Generates a simple 3-generation vessel tree
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
    vessel_1.p = 0*ones(vessel_1.N, 1);
    vessel_1.R = 1 * ones(vessel_1.N, 1);
    vessel_1.V = ones(vessel_1.N, 1);
    vessel_1.connected_to = [2 3];
    
    vessel_2.length = 50; %mm
    vessel_2.N = 1 + ceil(vessel_2.length/parameters.dx); % Number of grid points
    vessel_2.R_unstretched = 0.7; % Initial unstretched radius
    vessel_2.p = 1*ones(vessel_2.N, 1);
    vessel_2.R = 1 * ones(vessel_2.N, 1);
    vessel_2.V = ones(vessel_2.N, 1);
    vessel_2.connected_to = [4 5];

    vessel_3.length = 50; %mm
    vessel_3.N = 1 + ceil(vessel_3.length/parameters.dx); % Number of grid points
    vessel_3.R_unstretched = 0.7; % Initial unstretched radius
    vessel_3.p = 1*ones(vessel_3.N, 1);
    vessel_3.R = 1 * ones(vessel_3.N, 1);
    vessel_3.V = ones(vessel_3.N, 1);
    vessel_3.connected_to = [6 7];

    vessel_4.length = 50; %mm
    vessel_4.N = 1 + ceil(vessel_4.length/parameters.dx); % Number of grid points
    vessel_4.R_unstretched = 0.5; % Initial unstretched radius
    vessel_4.p = 1*ones(vessel_4.N, 1);
    vessel_4.R = 1 * ones(vessel_4.N, 1);
    vessel_4.V = ones(vessel_4.N, 1);
    vessel_4.connected_to = [];
    
    vessel_5.length = 50; %mm
    vessel_5.N = 1 + ceil(vessel_5.length/parameters.dx); % Number of grid points
    vessel_5.R_unstretched = 0.5; % Initial unstretched radius
    vessel_5.p = 1*ones(vessel_5.N, 1);
    vessel_5.R = 1 * ones(vessel_5.N, 1);
    vessel_5.V = ones(vessel_5.N, 1);
    vessel_5.connected_to = [];
        
    vessel_6.length = 50; %mm
    vessel_6.N = 1 + ceil(vessel_6.length/parameters.dx); % Number of grid points
    vessel_6.R_unstretched = 0.5; % Initial unstretched radius
    vessel_6.p = 1*ones(vessel_6.N, 1);
    vessel_6.R = 1 * ones(vessel_6.N, 1);
    vessel_6.V = ones(vessel_6.N, 1);
    vessel_6.connected_to = [];
    
    vessel_7.length = 50; %mm
    vessel_7.N = 1 + ceil(vessel_7.length/parameters.dx); % Number of grid points
    vessel_7.R_unstretched = 0.5; % Initial unstretched radius
    vessel_7.p = 1*ones(vessel_7.N, 1);
    vessel_7.R = 1 * ones(vessel_7.N, 1);
    vessel_7.V = ones(vessel_7.N, 1);
    vessel_7.connected_to = [];

    vessel_tree(1) = vessel_1;
    vessel_tree(2) = vessel_2;
    vessel_tree(3) = vessel_3;
    vessel_tree(4) = vessel_4;
    vessel_tree(5) = vessel_5;
    vessel_tree(6) = vessel_6;
    vessel_tree(7) = vessel_7;
    
    