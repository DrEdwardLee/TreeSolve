function TDL_PlotVesselTree(vessel_tree, parameters, colour, handle_matrix, fit_to_axis)
    % TDL_PlotVesselTree Draw graphs showing parameter values in each vessel of a TreeSolve tree
    %
    %     Author: Tom Doel www.tomdoel.com
    %     Part of TreeSolve. http://github.com/tomdoel/TreeSolve
    %     Distributed under the GNU GPL v3 licence. Please see LICENSE file.
    %    
    
    num_vessels = length(vessel_tree);
    
    prev_current_generation = 0;
    
    for vessel_index = 1 : num_vessels
        
        current_generation = 1 + floor(log2(vessel_index));
        
        vessel = vessel_tree(vessel_index);
        x_values = [0 : parameters.dx : parameters.dx*(vessel.N - 1)];
        
        subplot(handle_matrix(1, vessel_index));
        HoldOnOff(prev_current_generation, current_generation);

        plot(x_values, vessel.R);
        if (isfield(parameters, 'R_plot_range'))
            axis([0 vessel.length parameters.R_plot_range]);
        end

        if (fit_to_axis)
            axis tight;
        end

        subplot(handle_matrix(2, vessel_index));
        HoldOnOff(prev_current_generation, current_generation);

        plot(x_values, vessel.V);
        if (isfield(parameters, 'V_plot_range'))
            axis([0 vessel.length parameters.V_plot_range]);
        end

        if (fit_to_axis)
            axis tight;
        end

        subplot(handle_matrix(3, vessel_index));
        HoldOnOff(prev_current_generation, current_generation);

        plot(x_values, vessel.p);
        if (isfield(parameters, 'p_plot_range'))
            axis([0 vessel.length parameters.p_plot_range]);
        end

        if (fit_to_axis)
            axis tight;
        end

        prev_current_generation = current_generation;
    end
end

function HoldOnOff(prev_current_generation, current_generation)
    if (prev_current_generation ~= current_generation)
        hold off;
    else
        hold on;
        hold all;
    end
end
            
