function TDL_DrawVesselTree(vessel_tree)
    % TDL_DrawVesselTree Visualise a TreeSolve tree
    %
    %     Author: Tom Doel www.tomdoel.com
    %     Part of TreeSolve. http://github.com/tomdoel/TreeSolve
    %     Distributed under the GNU GPL v3 licence. Please see LICENSE file.
    %    
    
    figure;
    for i = 1 : length(vessel_tree)
        this_vessel = vessel_tree(i);
        start_point = this_vessel.start_point;
        end_point = this_vessel.end_point;
        radius = this_vessel.R_unstretched_average;

        line_handle = line(...
            [start_point(1) end_point(1)], [start_point(2) end_point(2)], [start_point(3) end_point(3)], ...
            'Color', 'r', 'LineWidth', 3*radius);
    end

    camorbit(30, 70)
    xlabel('x'); ylabel('y'); zlabel('z');
    axis equal;
end

