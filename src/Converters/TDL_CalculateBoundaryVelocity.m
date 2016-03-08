function V1_t1 = TDL_CalculateBoundaryVelocity(point1_t0, point2_t0, point2_t1, p1_t1, R1_t1, R_dash, V_dash, angle, parameters, which_direction)
    % TDL_CalculateBoundaryVelocity Helper function for TreeSolve
    %
    %
    %     Author: Tom Doel www.tomdoel.com
    %     Part of TreeSolve. http://github.com/tomdoel/TreeSolve
    %     Distributed under the GNU GPL v3 licence. Please see LICENSE file.
    %    
    
    dx = (which_direction)*parameters.dx;
    rho = parameters.rho;
    nu = parameters.nu;
    alpha = parameters.alpha;
    dt = parameters.dt;
    g = parameters.gravity;
    
    gravity_term = 2*dx*g*cos(angle);
        
    % Time k
    p1_k = point1_t0.p;
    R1_k = point1_t0.R;
    p2_k = point2_t0.p;
    V2_k = point2_t0.V;
    R2_k = point2_t0.R;
    
    % Time k+1
    p1_kp1 = p1_t1;
    R1_kp1 = R1_t1;
    p2_kp1 = point2_t1.p;
    R2_kp1 = point2_t1.R;
    
    V1_t1 = V2_k + (dt/(2*dx))*(...
        2*alpha*(V_dash.^2)*(R2_kp1 + R2_k - R1_kp1 - R1_k)/R_dash ...
        - (p2_kp1 + p2_k - p1_kp1 - p1_k)/rho - gravity_term ...
        ) - 2*dt*nu*(alpha/(alpha-1))*V_dash/(R_dash.^2) ...
        + (1/R_dash)*(dx/dt + 2*alpha*V_dash)*(R1_kp1 + R2_kp1 - R1_k - R2_k) ...
        - 2*V_dash*(R1_kp1 - R2_k)/R_dash;
    