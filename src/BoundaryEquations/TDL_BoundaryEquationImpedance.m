function remainder = TDL_BoundaryEquationImpedance(time, pressure, flow, parameters, which_direction)
    % TDL_BoundaryEquationImpedance Impedance boundary condition function for use with TDL_SolveVesselTree
    %
    %     Author: Tom Doel www.tomdoel.com
    %     Part of TreeSolve. http://github.com/tomdoel/TreeSolve
    %     Distributed under the GNU GPL v3 licence. Please see LICENSE file.
    %    

    
    % Steady applied impedance at outlet and flow at inlet. Solve for remainder=0
    if (which_direction == 1) 
        remainder = flow - parameters.inlet_applied_flow;
    else
        remainder = pressure - flow*parameters.outlet_applied_impedance;
    end
end