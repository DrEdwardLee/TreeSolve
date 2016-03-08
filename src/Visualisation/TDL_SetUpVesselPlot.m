function handle_matrix = TDL_SetUpVesselPlot(num_vessels, parameters)
    % TDL_SetUpVesselPlot Sets up a subplot for vessel tree visualisation
    %
    %     Author: Tom Doel www.tomdoel.com
    %     Part of TreeSolve. http://github.com/tomdoel/TreeSolve
    %     Distributed under the GNU GPL v3 licence. Please see LICENSE file.
    %    
    
    num_generations = log2(num_vessels + 1);
    handle_matrix = zeros(3, num_vessels);
    
    figure;
    maximize;

    % Group same generations together in plots
    number_of_columns = num_generations;

    for vessel_index = 1 : num_vessels    
        current_generation = 1 + floor(log2(vessel_index));

        subplot_index = current_generation;

        handle_matrix(1, vessel_index) = TDL_subplot(3, number_of_columns, subplot_index);
        title(['Radius, generation ' int2str(current_generation)]);

        handle_matrix(2, vessel_index) = TDL_subplot(3, number_of_columns, subplot_index + number_of_columns);
        title(['Velocity, generation ' int2str(current_generation)]);

        handle_matrix(3, vessel_index) = TDL_subplot(3, number_of_columns, subplot_index + 2*number_of_columns);
        title(['Pressure, generation ' int2str(current_generation)]);
    end
end

