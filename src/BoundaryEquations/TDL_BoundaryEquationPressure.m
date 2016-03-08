function remainder = TDL_BoundaryEquationPressure(time, pressure, flow, parameters, which_direction)
    % TDL_BoundaryEquationPressure Pressure boundary condition function for use with TDL_SolveVesselTree
    %
    %     Author: Tom Doel www.tomdoel.com
    %     Part of TreeSolve. http://github.com/tomdoel/TreeSolve
    %     Distributed under the GNU GPL v3 licence. Please see LICENSE file.
    %    

    % Steady applied pressures at inlet and outlet. . Solve for remainder=0
    if (which_direction == 1) 
        remainder = pressure - parameters.inlet_applied_pressure;
    else
        remainder = pressure - parameters.outlet_applied_pressure;
    end
end