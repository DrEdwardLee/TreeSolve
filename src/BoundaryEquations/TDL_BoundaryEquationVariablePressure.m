function remainder = TDL_BoundaryEquationVariablePressure(time, pressure, flow, parameters, which_direction)
    % TDL_BoundaryEquationVariablePressure Variable pressure boundary condition function for use with TDL_SolveVesselTree
    %
    %     Author: Tom Doel www.tomdoel.com
    %     Part of TreeSolve. http://github.com/tomdoel/TreeSolve
    %     Distributed under the GNU GPL v3 licence. Please see LICENSE file.
    %    

    % Applied pressures at inlet and outlet increasing from zero to fixed values. Solve for remainder=0
    % Some solution schemes are more stable if these variable pressure boundary
    % conditions are used
    if (which_direction == 1) 
        remainder = pressure - (min(time, 0.1)/0.1).*(parameters.inlet_applied_pressure);
    else
        remainder = pressure - (min(time, 0.1)/0.1).*(parameters.outlet_applied_pressure);
    end
end