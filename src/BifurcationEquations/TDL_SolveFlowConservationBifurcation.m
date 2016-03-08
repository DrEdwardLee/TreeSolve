function [new_parent_point, new_child_point_1, new_child_point_2] = TDL_SolveFlowConservationBifurcation(parent_vessel, child_vessel_1, child_vessel_2, parameters)
    % TDL_SolveFlowConservationBifurcation Funtion for use with TDL_SolveVesselTreeWithFlowConservation to solve solve bifurcation equations using flow conservation equations
    %
    % Implementation of a scheme proposed by Jonathan Whiteley, Kelly
    % Burrowes, Kathryn Gillow, David Gavaghan (2008), University of Oxford.
    %
    %     Author: Tom Doel www.tomdoel.com
    %     Part of TreeSolve. http://github.com/tomdoel/TreeSolve
    %     Distributed under the GNU GPL v3 licence. Please see LICENSE file.
    %    
    
    % Starting guess: use previous values
    p_vector_guess = [parent_vessel.p(end); child_vessel_1.p(1); child_vessel_2.p(1)];

    
    if (~parameters.use_newton_solver)    
        options = optimset('Display', 'off');


        [p_vector_solution, fval, exitflag, output] = fsolve(@NewtonRHS2, p_vector_guess, options, ...
            parent_vessel, child_vessel_1, child_vessel_2, parameters);

        if (exitflag < 1)
            disp('*** Solver failed');
            fval
            exitflag
            output
            error('Solver failed');
        end

        new_parent_point = TDL_VesselBoundaryEndFromPressure(p_vector_solution(1), parent_vessel, parameters);
        new_child_point_1 = TDL_VesselBoundaryStartFromPressure(p_vector_solution(2), child_vessel_1, parameters);
        new_child_point_2 = TDL_VesselBoundaryStartFromPressure(p_vector_solution(3), child_vessel_2, parameters);
    
    else

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

            [p_vector_k1, step_norm] = SolveWhiteleyBifurcation( ...
                new_parent_point, new_child_point_1, new_child_point_2, ...
                parent_vessel, child_vessel_1, child_vessel_2, parameters);

        end
    end






function [p_vector_k1, step_norm] = SolveWhiteleyBifurcation(...
    new_parent_point, new_child_point_1, new_child_point_2, ...
    parent_vessel, child_vessel_1, child_vessel_2, parameters)

    parent_point_t0.p = parent_vessel.p(end); parent_point_t0.R = parent_vessel.R(end); parent_point_t0.V = parent_vessel.V(end);
    child_point_1_t0.p = child_vessel_1.p(1); child_point_1_t0.R = child_vessel_1.R(1); child_point_1_t0.V = child_vessel_1.V(1);
    child_point_2_t0.p = child_vessel_2.p(1); child_point_2_t0.R = child_vessel_2.R(1); child_point_2_t0.V = child_vessel_2.V(1);

    R_vector_unstretched = [parent_vessel.R_unstretched(end); child_vessel_1.R_unstretched(1); child_vessel_2.R_unstretched(1)];
    V_prime_vector_k = [parent_vessel.V_half(end); child_vessel_1.V_half(1); child_vessel_2.V_half(1)];
    R_prime_vector_k = [parent_vessel.R_half(end); child_vessel_1.R_half(1); child_vessel_2.R_half(1)];

    p_vector_k1 = [new_parent_point.p; new_child_point_1.p; new_child_point_2.p];
    
    A = NewtonMatrix(new_parent_point, new_child_point_1, new_child_point_2, ...
        V_prime_vector_k, R_prime_vector_k, R_vector_unstretched, ...
        child_vessel_1.angle, child_vessel_2.angle, parameters);
        
    b = NewtonRHS(new_parent_point, new_child_point_1, new_child_point_2, ...
        child_vessel_1.angle, child_vessel_2.angle, R_vector_unstretched, parameters);
    
    
    
    s = A\b;
    
    
    p_vector_k1 = p_vector_k1 + s;

    step_norm = norm(s/p_vector_k1);


function A = NewtonMatrix(a1_t1, b1_t1, c1_t1, V_prime_vector_k, R_prime_vector_k, R_vector_unstretched, ...
    alpha_2, alpha_3, parameters)
    A = zeros(3);

    V_prime_a = V_prime_vector_k(1);
    V_prime_b = V_prime_vector_k(2);
    V_prime_c = V_prime_vector_k(3);
    
    R_prime_a = R_prime_vector_k(1);
    R_prime_b = R_prime_vector_k(2);
    R_prime_c = R_prime_vector_k(3);
        
    dFdp_a1_k1 = BoundarydFdp(a1_t1, V_prime_a, R_prime_a, R_vector_unstretched(1), parameters, -1);
    dFdp_b1_k1 = BoundarydFdp(b1_t1, V_prime_b, R_prime_b, R_vector_unstretched(2), parameters, 1);
    dFdp_c1_k1 = BoundarydFdp(c1_t1, V_prime_c, R_prime_c, R_vector_unstretched(3), parameters, 1);
    
    dRdp_a1_t1 = CalculatedRdpFromPressure(a1_t1.p, R_vector_unstretched(1), parameters);
    dRdp_b1_t1 = CalculatedRdpFromPressure(b1_t1.p, R_vector_unstretched(2), parameters);
    dRdp_c1_t1 = CalculatedRdpFromPressure(c1_t1.p, R_vector_unstretched(3), parameters);
    
    dVdp_a1_t1 = CalculatedVdp(dRdp_a1_t1, V_prime_a, R_prime_a, parameters, -1);
    dVdp_b1_t1 = CalculatedVdp(dRdp_b1_t1, V_prime_b, R_prime_b, parameters, 1);
    dVdp_c1_t1 = CalculatedVdp(dRdp_c1_t1, V_prime_c, R_prime_c, parameters, 1);
    
    dh1dp_a1_k1 = Calculatedh1dp(a1_t1, dFdp_a1_k1, dRdp_a1_t1, dVdp_a1_t1, R_vector_unstretched(1), parameters);
    dh1dp_b1_k1 = Calculatedh1dp(b1_t1, dFdp_b1_k1, dRdp_b1_t1, dVdp_b1_t1, R_vector_unstretched(2), parameters);
    dh1dp_c1_k1 = Calculatedh1dp(c1_t1, dFdp_c1_k1, dRdp_c1_t1, dVdp_c1_t1, R_vector_unstretched(3), parameters);
    
    A(1, 1) = -dh1dp_a1_k1;
    A(1, 2) = dh1dp_b1_k1*cos(alpha_2);
    A(1, 3) = dh1dp_c1_k1*cos(alpha_3);
    A(2, 1) = 0;
    A(2, 2) = dh1dp_b1_k1*sin(alpha_2);
    A(2, 3) = dh1dp_c1_k1*sin(alpha_3); % Note positive
    A(3, 1) = -dFdp_a1_k1;
    A(3, 2) = dFdp_b1_k1;
    A(3, 3) = dFdp_c1_k1;
    
function h = h1(point, R_0, parameters)
    A = pi*point.R.^2;
    Q = A.*point.V;
    f_hat = point.p;
    g = Calculateg(A, R_0, parameters);
    
    h = parameters.rho*parameters.alpha*(Q.^2)./A + A.*f_hat - g;
    
    f_hat_test = Calculatef_hat(A, R_0, parameters);
    if (abs(f_hat - f_hat_test) > 0.00001)
        disp(['no match: f_hat ' num2str(f_hat) ' ' num2str(f_hat_test)]);
    end
    
    
function p = Calculatef_hat(A, R_0, parameters)
    p = parameters.G_0*((A./(pi*R_0.^2)).^(parameters.beta/2) - 1);
    
function dh1dp = Calculatedh1dp(point, dQdp, dRdp, dVdp, R_0, parameters)
    rho = parameters.rho;
    beta = parameters.beta;
    alpha = parameters.alpha;
    G_0 = parameters.G_0;
    Q = Flow(point);
    R = point.R;
    V = point.V;
    A = pi*R.^2;
    
    dh1dp_test = 2*rho*alpha*Q*(dQdp - (Q/R)*dRdp)/(pi*R.^2) + pi*R.^2;
    
    dh1dA = -rho*alpha*Q.^2/A.^2 + G_0*(beta/2)*(A/(pi*R_0^2)).^(beta/2);
    dh1dQ = 2*rho*alpha*Q./A;
    dAdp = 2*pi*R*dRdp;
    
    dQdp = V*dAdp + A*dVdp;
    
    dh1dp = 2*pi*R*dh1dA*dRdp + dh1dQ*dQdp;
    
    if (abs(dh1dp - dh1dp_test) > 0.001)
        disp(['no match: dh1_dp ' num2str(dh1dp) ' ' num2str(dh1dp_test)]);
    end

function g = Calculateg(A, R_0, parameters)
    beta = parameters.beta;
    G_0 = parameters.G_0;
       
    g = (G_0*A^(beta/2+1))/(((pi*R_0^2)^(beta/2))*(beta/2+1)) - G_0*A + G_0*pi*R_0^2*beta/(beta+2);
    
    
    
    g_test = G_0*( (1/(1 + beta/2))*(A/(pi*R_0^2))^(beta/2)*A - A + beta*pi*R_0^2/(2 + beta) );

    if (abs(g - g_test) > 0.0001)
        error('failure in g');
    end
    
    
function rhs = NewtonRHS2(p_vector_t1, parent_vessel_t0, child1_vessel_t0, child2_vessel_t0, parameters)
    
    parent_point_t1 = TDL_VesselBoundaryEndFromPressure(p_vector_t1(1), parent_vessel_t0, parameters);
    child1_point_t1 = TDL_VesselBoundaryStartFromPressure(p_vector_t1(2), child1_vessel_t0, parameters);
    child2_point_t1 = TDL_VesselBoundaryStartFromPressure(p_vector_t1(3), child2_vessel_t0, parameters);

    F_a1_t1 = Flow(parent_point_t1);
    F_b1_t1 = Flow(child1_point_t1);
    F_c1_t1 = Flow(child2_point_t1);
    
    h1_a1 = h1(parent_point_t1, parent_vessel_t0.R_unstretched(end), parameters);
    h1_b1 = h1(child1_point_t1, child1_vessel_t0.R_unstretched(1), parameters);
    h1_c1 = h1(child2_point_t1, child2_vessel_t0.R_unstretched(1), parameters);

    alpha_2 = child1_vessel_t0.angle;
    alpha_3 = child2_vessel_t0.angle;
    
    rhs = zeros(3, 1);
    rhs(1) = h1_a1 - h1_b1*cos(alpha_2) - h1_c1*cos(alpha_3);
    rhs(2) = h1_b1*sin(alpha_2) + h1_c1*sin(alpha_3);
    rhs(3) = - F_a1_t1 + F_b1_t1 + F_c1_t1;
    
    
    
    

function rhs = NewtonRHS(a1_t1, b1_t1, c1_t1, alpha_2, alpha_3, R_vector_unstretched, parameters)
    F_a1_k1 = Flow(a1_t1);
    F_b1_k1 = Flow(b1_t1);
    F_c1_k1 = Flow(c1_t1);
    
    h1_a1 = h1(a1_t1, R_vector_unstretched(1), parameters);
    h1_b1 = h1(b1_t1, R_vector_unstretched(2), parameters);
    h1_c1 = h1(c1_t1, R_vector_unstretched(3), parameters);

    rhs = zeros(3, 1);
    rhs(1) = h1_a1 - h1_b1*cos(alpha_2) - h1_c1*cos(alpha_3);
    rhs(2) = - h1_b1*sin(alpha_2) - h1_c1*sin(alpha_3);
    rhs(3) = F_a1_k1 - F_b1_k1 - F_c1_k1;    
    
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

    
    
    
  