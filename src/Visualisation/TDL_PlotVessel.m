function TDL_PlotVessel(vessel_tree, vessel_number, parameters, colour, handle_matrix, fit_to_axis)
    % TDL_PlotVessel Updates individual vessel graphs in a TreeSolve tree visualisation
    %
    %     Author: Tom Doel www.tomdoel.com
    %     Part of TreeSolve. http://github.com/tomdoel/TreeSolve
    %     Distributed under the GNU GPL v3 licence. Please see LICENSE file.
    %    
    
    vessel = vessel_tree(vessel_number);
    x_values = [0 : parameters.dx : parameters.dx*(vessel.N - 1)];
    num_vessels = length(vessel_tree);

    subplot(handle_matrix(1, vessel_number));
    plot(x_values, vessel.R, colour);
    title(['Radius, vessel ' int2str(vessel_number)]);
    if (isfield(parameters, 'R_plot_range'))
        axis([0 vessel.length parameters.R_plot_range]);
    end

    if (fit_to_axis)
        axis tight;
    end

    subplot(handle_matrix(2, vessel_number));
    plot(x_values, vessel.V, colour);
    title(['Velocity, vessel ' int2str(vessel_number)]);
    if (isfield(parameters, 'V_plot_range'))
        axis([0 vessel.length parameters.V_plot_range]);
    end

    if (fit_to_axis)
        axis tight;
    end
   
    subplot(handle_matrix(3, vessel_number));
    plot(x_values, vessel.p, colour);
    title(['Pressure, vessel ' int2str(vessel_number)]);
    if (isfield(parameters, 'p_plot_range'))
        axis([0 vessel.length parameters.p_plot_range]);
    end

    if (fit_to_axis)
        axis tight;
    end
