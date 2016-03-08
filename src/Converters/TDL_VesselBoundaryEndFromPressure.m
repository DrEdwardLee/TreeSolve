function point_t1 = TDL_VesselBoundaryEndFromPressure(pN_t1, vessel_t0, parameters)
    % TDL_VesselBoundaryEndFromPressure Helper function for TreeSolve
    %
    %
    %     Author: Tom Doel www.tomdoel.com
    %     Part of TreeSolve. http://github.com/tomdoel/TreeSolve
    %     Distributed under the GNU GPL v3 licence. Please see LICENSE file.
    % 
    

    point1_t0.p = vessel_t0.p(end); point1_t0.V = vessel_t0.V(end); point1_t0.R = vessel_t0.R(end);
    point2_t0.p = vessel_t0.p(end-1); point2_t0.V = vessel_t0.V(end-1); point2_t0.R = vessel_t0.R(end-1);
    point2_t1.p = vessel_t0.p_int_t1(end); point2_t1.V = vessel_t0.V_int_t1(end); point2_t1.R = vessel_t0.R_int_t1(end);

    point_t1.p = pN_t1;
    point_t1.R = TDL_CalculateRadiusFromPressure(pN_t1, vessel_t0.R_unstretched(end), parameters);
    point_t1.V = TDL_CalculateBoundaryVelocity( ...
        point1_t0, point2_t0, point2_t1, pN_t1, point_t1.R, vessel_t0.R_half(end), vessel_t0.V_half(end), vessel_t0.angle_to_vertical, parameters, -1);



