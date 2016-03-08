function TDL_PlotCompareFlows(results_1, results_2, legend_1, legend_2, filename)
    % TDL_PlotCompareFlows Draw a graph comparing flow results
    %
    %     Author: Tom Doel www.tomdoel.com
    %     Part of TreeSolve. http://github.com/tomdoel/TreeSolve
    %     Distributed under the GNU GPL v3 licence. Please see LICENSE file.
    %    
    

    figure;
    hold on;

    PlotFlows(results_1.vessel_tree, results_1.parameters, 'r');
    PlotFlows(results_2.vessel_tree, results_2.parameters, 'b--');


    set(0,'DefaultAxesColorOrder',[1 0 0;0 1 0;0 0 1],...
          'DefaultAxesLineStyleOrder','-|--|:')

    page_width = 16;
    ratio = 6/3;
    font_size = 9;
    %font_size = 15;
    legend_size = 6;
    axis tight;
    
    % We want y to start at zero
    axis_scale = axis;
    axis_scale(3) = 0;
    axis(axis_scale);

    h_legend = legend(legend_1, legend_2, 'Location', 'NorthEast');
    set(h_legend,'FontSize', legend_size);

    %set(h_legend, 'Box', 'off')
    xlabel('Distance from inlet (mm)');
    ylabel('Flow (mm^3/s)');
    
    
    text(10, 13000, 'Generation 1', 'HorizontalAlignment', 'center');
    text(27, 13000, '2', 'HorizontalAlignment', 'center');
    text(39, 13000, '3', 'HorizontalAlignment', 'center');
    text(46, 13000, '4', 'HorizontalAlignment', 'center');
    text(53, 13000, '5', 'HorizontalAlignment', 'center');

    saveaspngandeps(-1, filename, page_width, ratio, font_size);
    
   

end



function PlotFlows(vessel_tree, parameters, color)
    num_vessels = length(vessel_tree);
    starting_coordinate = 0;
    starting_prev_coords{1} = [];
    starting_prev_values{1} = [];
    
    for vessel_index = 1 : num_vessels
        vessel = vessel_tree(vessel_index);
        x_values = 0 : parameters.dx : parameters.dx*(vessel.N - 1);
        x_values = x_values + starting_coordinate(vessel_index);
%        x_values = [starting_prev_coords{vessel_index} x_values];
        flows = pi*(vessel.V).*((vessel.R).^2);
%        flows = [starting_prev_values{vessel_index} ; flows];

        hobj = plot(x_values, flows, color);
        
        if (vessel_index > 1)
            hAnnotation = get(hobj,'Annotation');
            hLegendEntry = get(hAnnotation','LegendInformation');
            set(hLegendEntry,'IconDisplayStyle','off')
        end

        children = vessel.connected_to;
        if (~isempty(children))
           child_1 = children(1); 
           child_2 = children(2); 

           starting_coordinate(child_1) = starting_coordinate(vessel_index) + vessel.length;
           starting_coordinate(child_2) = starting_coordinate(vessel_index) + vessel.length;
            
           starting_prev_coords{child_1} = x_values;
           starting_prev_values{child_1} = flows;
           starting_prev_coords{child_2} = x_values;
           starting_prev_values{child_2} = flows;

        end        
    end
end
