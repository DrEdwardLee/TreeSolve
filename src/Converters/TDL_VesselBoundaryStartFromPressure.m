function point_t1 = TDL_VesselBoundaryStartFromPressure(p1_t1, vessel_t0, parameters)
    % TDL_VesselBoundaryStartFromPressure Helper function for TreeSolve
    %
    %
    %     Author: Tom Doel www.tomdoel.com
    %     Part of TreeSolve. http://github.com/tomdoel/TreeSolve
    %     Distributed under the GNU GPL v3 licence. Please see LICENSE file.
    % 
    
    point1_t0.p = vessel_t0.p(1); point1_t0.V = vessel_t0.V(1); point1_t0.R = vessel_t0.R(1);
    point2_t0.p = vessel_t0.p(2); point2_t0.V = vessel_t0.V(2); point2_t0.R = vessel_t0.R(2);
    point2_t1.p = vessel_t0.p_int_t1(1); point2_t1.V = vessel_t0.V_int_t1(1); point2_t1.R = vessel_t0.R_int_t1(1);
   
    point_t1.p = p1_t1;
    point_t1.R = TDL_CalculateRadiusFromPressure(p1_t1, vessel_t0.R_unstretched(1), parameters);
    point_t1.V = TDL_CalculateBoundaryVelocity(...
        point1_t0, point2_t0, point2_t1, point_t1.p, point_t1.R, vessel_t0.R_half(1), vessel_t0.V_half(1), vessel_t0.angle_to_vertical, parameters, 1);
