function [V_half, R_half, p_half, p_int_t1, V_int_t1, R_int_t1] = TDL_LWStep(V_t0, R_t0, p_t0, R_unstretched, angle, parameters) 
    % TDL_LWStep Lax_Wendroff step for solving Navier-Stokes equations
    %
    % This is a novel implementation adapted in part from a scheme by 
    % Burrowes at Tawhai, 2006, adapted from Smith et al., 2002.
    % Equation numbers in the comments refer to Smith et al. 2002.
    %
    %
    %     Author: Tom Doel www.tomdoel.com
    %     Part of TreeSolve. http://github.com/tomdoel/TreeSolve
    %     Distributed under the GNU GPL v3 licence. Please see LICENSE file.
    %    
 

    % Uses the Lax-Wendroff scheme to calculate values of p, V, R at intermediate points (at time
    % t+dt/2 and space x+dx/2) and the internal points at time t+dt (i.e.
    % excluding the boundary points at each end).
    [V_half, R_half, p_half] = CalculateIntermediatePoints(V_t0, R_t0, p_t0, R_unstretched, angle, parameters);
    [p_int_t1, V_int_t1, R_int_t1] = CalculateInternalPoints(V_t0, R_t0, V_half, R_half, p_half, R_unstretched, angle, parameters);

    
%LW step 1
function [V_half, R_half, p_half] = CalculateIntermediatePoints(V, R, p, R_0, angle, parameters)
    rho = parameters.rho;
    alpha = parameters.alpha;
    %beta = parameters.beta;
    nu = parameters.nu;
    dt = parameters.dt;
    dx = parameters.dx;
    g = parameters.gravity;
    
    V_plus  = V(2 : end, :) + V(1 : end-1, :);
    V_minus = V(2 : end, :) - V(1 : end-1, :);
    R_plus  = R(2 : end, :) + R(1 : end-1, :);
    R_minus = R(2 : end, :) - R(1 : end-1, :);
    p_minus = p(2 : end, :) - p(1 : end-1, :);
    R0_plus = (1/2)*(R_0(2 : end) + R_0(1 : end-1));


    gravity_term = dx*g*cos(angle);
    
    %Eqn 2.31
    V_half = V_plus/2 - (dt/(2*dx))*(...
        (2*alpha-1)*V_plus.*V_minus/2 + (alpha-1)*(V_plus.^2).*R_minus./R_plus + p_minus/rho + gravity_term ...
        ) - 2*dt*nu*alpha*(V_plus./(R_plus.^2))/(alpha-1);

    %Eqn 2.32
    R_half = R_plus/2 - (dt/(2*dx))*(R_plus.*V_minus/4 + V_plus.*R_minus/2);

    %Eqn 2.33
    p_half = TDL_CalculatePressureFromRadius(R_half, R0_plus, parameters);

    
% LW step 2
function [pint_t1, Vint_t1, Rint_t1] = CalculateInternalPoints(V_t0, R_t0, V_half, R_half, p_half, R_initial, angle, parameters)
    rho = parameters.rho;
    alpha = parameters.alpha;
    nu = parameters.nu;
    dt = parameters.dt;
    dx = parameters.dx;
    g = parameters.gravity;

    V_hplus  = V_half(2:end) + V_half(1:end-1);
    V_hminus = V_half(2:end) - V_half(1:end-1);
    R_hplus  = R_half(2:end) + R_half(1:end-1);
    R_hminus = R_half(2:end) - R_half(1:end-1);
    p_hminus = p_half(2:end) - p_half(1:end-1);

    gravity_term = dx*g*cos(angle);
        
    % eqn 2.34
    Vint_t1 = V_t0(2:end-1) - (dt/dx)*(...
            (2*alpha-1)*V_hplus.*V_hminus/2 + (alpha-1)*(V_hplus.^2).*R_hminus./R_hplus + p_hminus/rho + gravity_term ...
        ) - (4*dt*nu*alpha/(alpha-1))*V_hplus./(R_hplus.^2);

    % eqn 2.35
    Rint_t1 = R_t0(2:end-1) - (dt/dx)*(R_hplus.*V_hminus/4 + V_hplus.*R_hminus/2);

    % eqn 2.36
    pint_t1 = TDL_CalculatePressureFromRadius(Rint_t1, R_initial(2:end-1), parameters);
    
