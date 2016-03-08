function [new_parent_point, new_child_point_1, new_child_point_2] = TDL_SolvePressureBifurcation(...
    % TDL_SolvePressureBifurcation Funtion for use with TDL_SolveVesselTree to solve solve bifurcation equations using pressure conservation equations
    %
    %     Author: Tom Doel www.tomdoel.com
    %     Part of TreeSolve. http://github.com/tomdoel/TreeSolve
    %     Distributed under the GNU GPL v3 licence. Please see LICENSE file.
    %    

    parent_vessel, child_vessel_1, child_vessel_2, parameters)
    
    % Starting guess: use previous values
    p_vector_guess = [parent_vessel.p(end); child_vessel_1.p(1); child_vessel_2.p(1)];

    
  p_vector_k1 = p_vector_guess;  
    
    iteration = 0;
    tol = 0.001;
    
    
    step_norm = 1;
    while (true)
 
        new_parent_point = TDL_VesselBoundaryEndFromPressure(p_vector_k1(1), parent_vessel, parameters);
        new_child_point_1 = TDL_VesselBoundaryStartFromPressure(p_vector_k1(2), child_vessel_1, parameters);
        new_child_point_2 = TDL_VesselBoundaryStartFromPressure(p_vector_k1(3), child_vessel_2, parameters);
        
        if (iteration > 1000)
            error(['Warning: solution did not converge; iteration: ' num2str(iteration) ', norm:', num2str(step_norm)]);
            break
        end

        if (step_norm < tol)
            break
        end

        iteration = iteration + 1;
        
        [p_vector_k1, step_norm] = SolveSmithBifurcation( ...
            new_parent_point, new_child_point_1, new_child_point_2, ...
            parent_vessel, child_vessel_1, child_vessel_2, parameters);
  
    end




function [p_vector_k1, step_norm] = SolveSmithBifurcation(...
    new_parent_point, new_child_point_1, new_child_point_2, ...
    parent_vessel, child_vessel_1, child_vessel_2, parameters)

    parent_point_t0.p = parent_vessel.p(end); parent_point_t0.R = parent_vessel.R(end); parent_point_t0.V = parent_vessel.V(end);
    child_point_1_t0.p = child_vessel_1.p(1); child_point_1_t0.R = child_vessel_1.R(1); child_point_1_t0.V = child_vessel_1.V(1);
    child_point_2_t0.p = child_vessel_2.p(1); child_point_2_t0.R = child_vessel_2.R(1); child_point_2_t0.V = child_vessel_2.V(1);

    R_vector_unstretched = [parent_vessel.R_unstretched; child_vessel_1.R_unstretched; child_vessel_2.R_unstretched];
    V_prime_vector_k = [parent_vessel.V_half(end); child_vessel_1.V_half(1); child_vessel_2.V_half(1)];
    R_prime_vector_k = [parent_vessel.R_half(end); child_vessel_1.R_half(1); child_vessel_2.R_half(1)];

    p_vector_k1 = [new_parent_point.p; new_child_point_1.p; new_child_point_2.p];
    
    A = NewtonMatrix(new_parent_point, new_child_point_1, new_child_point_2, ...
        V_prime_vector_k, R_prime_vector_k, R_vector_unstretched, parameters);
    b = NewtonRHS(parent_point_t0, child_point_1_t0, child_point_2_t0, ...
        new_parent_point, new_child_point_1, new_child_point_2, parameters);
    
    
    
    s = A\b;
    
    %s = cgs(A, b, 10e-8, 100, diag(diag(A)));
    
    
    
    
    step_norm = norm(s);
    
    p_vector_k1 = p_vector_k1 + s;
    
    %disp(['Condition number of A:' num2str(cond(A))]);

    
    
    
    


function A = NewtonMatrix(a1_t1, b1_t1, c1_t1, V_prime_vector_k, R_prime_vector_k, R_vector_unstretched, parameters)
    A = zeros(3);

    L_a = parameters.L_a;
    L_b = parameters.L_b;
    L_c = parameters.L_c;
    dt = parameters.dt;

    V_prime_a = V_prime_vector_k(1);
    V_prime_b = V_prime_vector_k(2);
    V_prime_c = V_prime_vector_k(3);
    
    R_prime_a = R_prime_vector_k(1);
    R_prime_b = R_prime_vector_k(2);
    R_prime_c = R_prime_vector_k(3);
        
    dFdp_a1_k1 = BoundarydFdp(a1_t1, V_prime_a, R_prime_a, R_vector_unstretched(1), parameters, -1);
    dFdp_b1_k1 = BoundarydFdp(b1_t1, V_prime_b, R_prime_b, R_vector_unstretched(2), parameters, 1);
    dFdp_c1_k1 = BoundarydFdp(c1_t1, V_prime_c, R_prime_c, R_vector_unstretched(3), parameters, 1);
    
    A(1, 1) = 1 - 2*L_a*dFdp_a1_k1/dt;
    A(1, 2) = -1 - 2*L_b*dFdp_b1_k1/dt;
    A(1, 3) = 0;
    A(2, 1) = 1 - 2*L_a*dFdp_a1_k1/dt;
    A(2, 2) = 0;
    A(2, 3) = -1 - 2*L_c*dFdp_c1_k1/dt;
    A(3, 1) = dFdp_a1_k1;
    A(3, 2) = -dFdp_b1_k1;
    A(3, 3) = -dFdp_c1_k1;
    


function rhs = NewtonRHS(a1_t0, b1_t0, c1_t0, a1_t1, b1_t1, c1_t1, parameters)
    dt = parameters.dt;
    L_a = parameters.L_a;
    L_b = parameters.L_b;
    L_c = parameters.L_c;

    % Pressure values at time k+1
    p_a1_k1 = a1_t1.p;
    p_b1_k1 = b1_t1.p;
    p_c1_k1 = c1_t1.p;
    
    % Pressure values at time k
    p_a1_k = a1_t0.p;
    p_b1_k = b1_t0.p;
    p_c1_k = c1_t0.p;

    F_a1_k1 = Flow(a1_t1);
    F_b1_k1 = Flow(b1_t1);
    F_c1_k1 = Flow(c1_t1);

    F_a1_k = Flow(a1_t0);
    F_b1_k = Flow(b1_t0);
    F_c1_k = Flow(c1_t0);

    rhs = zeros(3, 1);
    rhs(1) = (2*L_a/dt)*(F_a1_k1-F_a1_k) + (2*L_b/dt)*(F_b1_k1-F_b1_k) ...
        + p_b1_k1 + p_b1_k - p_a1_k1 - p_a1_k;
    rhs(2) = (2*L_a/dt)*(F_a1_k1-F_a1_k) + (2*L_c/dt)*(F_c1_k1-F_c1_k) ...
        + p_c1_k1 + p_c1_k - p_a1_k1 - p_a1_k;
    rhs(3) = F_c1_k1 + F_b1_k1 - F_a1_k1;
    
    
function F = Flow(point)
   F = pi*(point.R.^2).*point.V;


function dFdp = BoundarydFdp(point_t1, V_prime, R_prime, R_unstretched, parameters, vessel_direction)
% Smith et al. Eqn 2.64
    dRdp = CalculatedRdpFromPressure(point_t1.p, R_unstretched, parameters);
    dVdp = CalculatedVdp(dRdp, V_prime, R_prime, parameters, vessel_direction);
    dFdp = pi*(point_t1.R.^2).*dVdp + 2*pi*point_t1.R.*dRdp.*point_t1.V;



function dVdp = CalculatedVdp(dRdp, V_prime, R_prime, parameters, vessel_direction)
% Based on Smith et al eqn. 2.65
% vessel_direction = 1 for considering points at start of vessel
% and -1 for points at the end of vessel
    
    alpha = parameters.alpha;
    dx = vessel_direction*parameters.dx;
    dt = parameters.dt;
    
    dVdp = (dt/(2*dx))*(1/parameters.rho - 2*alpha*(V_prime.^2).*dRdp./R_prime) ...
        + (1./R_prime).*(dx/dt + 2*alpha*V_prime).*dRdp - 2*V_prime.*dRdp./R_prime;


function dRdp = CalculatedRdpFromPressure(p, R_unstretched, parameters)
%CALCULATEDRDPFROMPRESSURE Based on eqn. 2.66 in Smith et al 2002

    dRdp = (R_unstretched./(parameters.beta.*parameters.G_0)).*((p./parameters.G_0 + 1).^(1/parameters.beta - 1));

    
    

