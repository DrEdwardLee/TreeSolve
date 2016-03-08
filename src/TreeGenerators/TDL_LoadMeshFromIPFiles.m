function vessel_tree = TDL_LoadMeshFromIPFiles(folder_name, parameters)
    % TDL_LoadMeshFromIPFiles Generates an TreeSolve artifical tree using
    % the ipelem and ipnode files in the given directory
    %
    %     Author: Tom Doel www.tomdoel.com
    %     Part of TreeSolve. http://github.com/tomdoel/TreeSolve
    %     Distributed under the GNU GPL v3 licence. Please see LICENSE file.
    %    

    path = ['../Meshes/' folder_name '/'];
    [element_nodes node_coordinates rad_coordinates] = LoadMeshFromIPFiles(path, folder_name);
    
    vessel_tree = TDL_VesselTreeFromMesh(element_nodes, node_coordinates, rad_coordinates, parameters);
end


function [element_nodes node_coordinates rad_coordinates] = LoadMeshFromIPFiles(path, folder_name)
    element_nodes = ReadElementNodes([path folder_name '.ipelem']);
    node_coordinates = ReadNodesFile([path folder_name '.ipnode']);
    rad_coordinates = ReadRadiusFile([path 'radius_' folder_name '.ipfiel']);
end    


function element_nodes = ReadElementNodes(filename)
    elem_lines = TDL_ReadTextFile(filename);
    line_no = 0;
    
    % Header (2 lines + blank line)
    line_no = line_no + 1; line_text = elem_lines{line_no};
    if (~strcmp(line_text, 'CMISS Version 2.1  ipelem File Version 2'))
        ReportError(filename, line_no, line_text, 'Unexpected header string');
    end
    line_no = line_no + 1; line_text = elem_lines{line_no};
    if (~strcmp(line_text, 'Heading: parent'))
        ReportError(filename, line_no, line_text, 'Unexpected header string');
    end
    line_no = line_no + 1; line_text = elem_lines{line_no};
    if (~strcmp(line_text, ''))
        ReportError(filename, line_no, line_text, 'Unexpected header string');
    end
    
    % Read number of elements
    line_no = line_no + 1; line_text = elem_lines{line_no};
    number_of_elements = sscanf(line_text, 'The number of elements is [1]: %d');
    if (ischar(number_of_elements))
        ReportError(filename, line_no, line_text, 'Unable to read number of elements');
    end
    
    element_nodes = zeros(number_of_elements, 2);
    
    % Now read in elements
    for element_index = 1 : number_of_elements
        % First, look for empty line
        line_no = line_no + 1; line_text = elem_lines{line_no};
        if (~strcmp(line_text, ''))
            ReportError(filename, line_no, line_text, 'Was expecting empty line');
        end
    
        % Read in element number (the first value in the vector will be the
        % default value which we can ignore)
        line_no = line_no + 1; line_text = elem_lines{line_no};
        this_element_number_vec = sscanf(line_text, 'Element number [%d]: %d');
        this_element_number = this_element_number_vec(2);
        if (ischar(this_element_number))
          ReportError(filename, line_no, line_text, 'Unable to read element number');
        end

        % Read in number of geometric coordinates- we expect this to be 3
        line_no = line_no + 1; line_text = elem_lines{line_no};
        num_geometric_coords = sscanf(line_text, 'The number of geometric Xj-coordinates is [3]: %d');
        if (num_geometric_coords ~= 3)
            ReportError(filename, line_no, line_text, 'Unexpected number of geometric cordinates');            
        end

        % Read in basic functional types- we expect these to be 1
        for fun_type_index = 1 : 3
            line_no = line_no + 1; line_text = elem_lines{line_no};
            fun_type_vec = sscanf(line_text, 'The basis function type for geometric variable %d is [1]:  %d');
            if (fun_type_vec(2) ~= 1)
                ReportError(filename, line_no, line_text, 'Unexpected functional type');
            end
        end

        % Read in global numbers for this element
        line_no = line_no + 1; line_text = elem_lines{line_no};
        this_basis_global_numbers = sscanf(line_text, 'Enter the 2 global numbers for basis 1: %d %d');
        if (ischar(this_basis_global_numbers))
          ReportError(filename, line_no, line_text, 'Unable to read global numbers for basis');
        end
        
        element_nodes(this_element_number, :) = this_basis_global_numbers';
    end
end

function node_coordinates = ReadNodesFile(filename)
    elem_lines = TDL_ReadTextFile(filename);
    line_no = 0;

    % Header (2 lines + blank line)
    line_no = line_no + 1; line_text = elem_lines{line_no};
    if (~strcmp(line_text, 'CMISS Version 2.1  ipnode File Version 2'))
        ReportError(filename, line_no, line_text, 'Unexpected header string');
    end
    line_no = line_no + 1; line_text = elem_lines{line_no};
    if (~strcmp(line_text, 'Heading: parent'))
        ReportError(filename, line_no, line_text, 'Unexpected header string');
    end
    line_no = line_no + 1; line_text = elem_lines{line_no};
    if (~strcmp(line_text, ''))
        ReportError(filename, line_no, line_text, 'Unexpected header string');
    end
    
    % Read number of nodes
    line_no = line_no + 1; line_text = elem_lines{line_no};
    number_of_nodes_vec = sscanf(line_text, 'The number of nodes is [%d]: %d');
    if (ischar(number_of_nodes_vec))
        ReportError(filename, line_no, line_text, 'Unable to read number of nodes');
    end
    number_of_nodes = number_of_nodes_vec(2);
    
    % Read number of coordinates
    line_no = line_no + 1; line_text = elem_lines{line_no};
    number_of_coords_vec = sscanf(line_text, 'Number of coordinates [%d]: %d');
    if (ischar(number_of_coords_vec))
        ReportError(filename, line_no, line_text, 'Unable to read number of coordinates');
    end
    number_of_coords = number_of_coords_vec(2);

    % Check prompting for different versions
    for coord_index = 1 : number_of_coords
        line_no = line_no + 1; line_text = elem_lines{line_no};
        prompt = sscanf(line_text, ['Do you want prompting for different versions of nj=' num2str(coord_index) ' [N]? %s']);
        if (~strcmp(prompt, 'N'))
            ReportError(filename, line_no, line_text, 'Prompting for versions not as expected');
        end
    end
    
    % Check number derivatives
    for coord_index = 1 : number_of_coords
        line_no = line_no + 1; line_text = elem_lines{line_no};
        num_derivs = sscanf(line_text, ['The number of derivatives for coordinate ' num2str(coord_index) ' is [0]: %d']);
        if (num_derivs ~= 0)
            ReportError(filename, line_no, line_text, 'Number of derivates not as expected');
        end
    end
    
    node_coordinates = zeros(number_of_nodes, number_of_coords);
    
    % Now read in nodes
    for node_index = 1 : number_of_nodes

        % First, look for empty line
        line_no = line_no + 1; line_text = elem_lines{line_no};
        if (~strcmp(line_text, ''))
            ReportError(filename, line_no, line_text, 'Was expecting empty line');
        end
    
        % Read in node number (the first value in the vector will be the
        % default value which we can ignore)
        line_no = line_no + 1; line_text = elem_lines{line_no};
        this_node_number_vec = sscanf(line_text, 'Node number [%d]: %d');
        this_node_number = this_node_number_vec(2);
        if (ischar(this_node_number))
          ReportError(filename, line_no, line_text, 'Unable to read node number');
        end
        
        % Read in coordinates
        for coord_index = 1 : number_of_coords
            line_no = line_no + 1; line_text = elem_lines{line_no};
            coords_vec = sscanf(line_text, ['The Xj(' num2str(coord_index) ') coordinate is [%f]: %f']);
            coord = coords_vec(2);
            if (ischar(coord))
                ReportError(filename, line_no, line_text, 'Unable to read in coordinates');
            end
            node_coordinates(this_node_number, coord_index) = coord';
        end
    end
end

function radius_coordinates = ReadRadiusFile(filename)
    elem_lines = TDL_ReadTextFile(filename);
    line_no = 0;

    % Header (2 lines + blank line)
    line_no = line_no + 1; line_text = elem_lines{line_no};
    if (~strcmp(line_text, 'CMISS Version 2.1  ipfiel File Version 3'))
        ReportError(filename, line_no, line_text, 'Unexpected header string');
    end
    line_no = line_no + 1; line_text = elem_lines{line_no};
    if (~strcmp(line_text, 'Heading: parent'))
        ReportError(filename, line_no, line_text, 'Unexpected header string');
    end
    line_no = line_no + 1; line_text = elem_lines{line_no};
    if (~strcmp(line_text, ''))
        ReportError(filename, line_no, line_text, 'Unexpected header string');
    end
    
    % Read number of field variables
    line_no = line_no + 1; line_text = elem_lines{line_no};
    number_of_field_vars = sscanf(line_text, 'The number of field variables is [3]: %d');
    if (number_of_field_vars ~= 1)
        ReportError(filename, line_no, line_text, 'Number of field variables is not 1');
    end

    % Read number of nodes
    line_no = line_no + 1; line_text = elem_lines{line_no};
    number_of_nodes_vec = sscanf(line_text, 'The number of nodes is [%d]: %d');
    if (ischar(number_of_nodes_vec))
        ReportError(filename, line_no, line_text, 'Unable to read number of nodes');
    end
    number_of_nodes = number_of_nodes_vec(2);
    
    % Check for prompting of versions of field variable
    line_no = line_no + 1; line_text = elem_lines{line_no};
    prompt = sscanf(line_text, 'Do you want prompting for different versions of field variable 1 [N]? %s');
    if (~strcmp(prompt, 'Y'))
        ReportError(filename, line_no, line_text, 'Prompting for versions not as expected');
    end
    
    % Read number of derivatives for field variable
    line_no = line_no + 1; line_text = elem_lines{line_no};
    num_derivs = sscanf(line_text, 'The number of derivatives for field variable 1 is [0]: %d');
    if (num_derivs ~= 0)
        ReportError(filename, line_no, line_text, 'Number of derivatives not as expected');
    end

    radius_coordinates = {};
    
    % Read in radius values
    for node_index = 1 : number_of_nodes
        % First, look for empty line
        line_no = line_no + 1; line_text = elem_lines{line_no};
        if (~strcmp(line_text, ''))
            ReportError(filename, line_no, line_text, 'Was expecting empty line');
        end
    
        % Read in node number (the first value in the vector will be the
        % default value which we can ignore)
        line_no = line_no + 1; line_text = elem_lines{line_no};
        this_node_number_vec = sscanf(line_text, 'Node number [%d]: %d');
        this_node_number = this_node_number_vec(2);
        if (ischar(this_node_number))
          ReportError(filename, line_no, line_text, 'Unable to read node number');
        end
        
        % Read in number of versions for this node
        line_no = line_no + 1; line_text = elem_lines{line_no};
        num_versions = sscanf(line_text, 'The number of versions for field variable 1 is [1]: %d');
        if (ischar(num_versions))
            ReportError(filename, line_no, line_text, 'Unable to read number of versions');
        end
        
        radius_list = zeros(num_versions, 1);
        
        % Read field variables for each version
        for version_index = 1 : num_versions
            if (num_versions > 1)
                % Check version number
                line_no = line_no + 1; line_text = elem_lines{line_no};
                ver_number = sscanf(line_text, 'For version number %d:');
                if (ver_number ~= version_index)
                    ReportError(filename, line_no, line_text, 'Version number not as expected');
                end
            end
            
            % Field variable value
            line_no = line_no + 1; line_text = elem_lines{line_no};
            value_vec = sscanf(line_text, 'The field variable 1 value is [%f]: %f');
            if (ischar(value_vec))
                ReportError(filename, line_no, line_text, 'Unable to read radius value');
            end
            
            radius_value = value_vec(2);
            radius_list(version_index) = radius_value;
        end
        radius_coordinates{node_index} = radius_list;
    end    
end



function ReportError(filename, line_no, line_text, error_message)
    disp(['***** TDL_LoadMeshFromIPFiles : Error: Unable to parse ' line_text]);
    error(['Error reading file ' filename ' on line ' int2str(line_no) ' : ' error_message]);
end