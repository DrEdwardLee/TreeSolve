function F = TDL_FlowConservationFunction(big_vector, vessel_tree, parameters, time)
    % TDL_FlowConservationFunction Solver for use with TDL_SolveVesselTreeWithFlowConservation
    %
    % Implementation of a scheme proposed by Jonathan Whiteley, Kelly
    % Burrowes, Kathryn Gillow, David Gavaghan (2008), University of Oxford.
    %
    %     Author: Tom Doel www.tomdoel.com
    %     Part of TreeSolve. http://github.com/tomdoel/TreeSolve
    %     Distributed under the GNU GPL v3 licence. Please see LICENSE file.
    %    

    F = zeros(size(big_vector));
    
    equation_vector_index = 1;
    num_vessels = length(vessel_tree);
    
    for vessel_index = 1 : num_vessels
        vessel = vessel_tree(vessel_index);
        A_t0 = vessel_tree(vessel_index).A;
        Q_t0 = vessel_tree(vessel_index).Q;
        
        R_unstretched_a = vessel_tree(vessel_index).R_unstretched;
                
        [A0_t1, Q0_t1, AN_t1, QN_t1] = GetVesselVariables(big_vector, vessel_index);
        
        A_t1_internal = vessel_tree(vessel_index).A_int_t1;
        Q_t1_internal = vessel_tree(vessel_index).Q_int_t1;
        A_t1 = [A0_t1; A_t1_internal; AN_t1];
        Q_t1 = [Q0_t1; Q_t1_internal; QN_t1];

        f1 = Calculatef1(A_t0, Q_t0, A_t1, Q_t1, parameters);
        F(equation_vector_index) = f1;
        equation_vector_index = equation_vector_index + 1;
      
        f2 = Calculatef2(vessel, A0_t1, Q0_t1, AN_t1, QN_t1, parameters);
        F(equation_vector_index) = f2;
        equation_vector_index = equation_vector_index + 1;
 
        
        if (~isempty(vessel_tree(vessel_index).connected_to))
            vessel_index_child_1 = vessel_tree(vessel_index).connected_to(1);
            vessel_index_child_2 = vessel_tree(vessel_index).connected_to(2);
            
            R_unstretched_b = vessel_tree(vessel_index_child_1).R_unstretched;
            R_unstretched_c = vessel_tree(vessel_index_child_2).R_unstretched;
            
            alpha_b = vessel_tree(vessel_index_child_1).angle;
            alpha_c = vessel_tree(vessel_index_child_2).angle;
            
            [A0_b_t1, Q0_b_t1, AN_b_t1, QN_b_t1] = GetVesselVariables(big_vector, vessel_index_child_1);
            [A0_c_t1, Q0_c_t1, AN_c_t1, QN_c_t1] = GetVesselVariables(big_vector, vessel_index_child_2);
        
            f3 = - Q_t1(end) + Q0_b_t1 + Q0_c_t1;  
            F(equation_vector_index) = f3;
            equation_vector_index = equation_vector_index + 1;
 
            h1_aN = Calculateh1(AN_t1, QN_t1, R_unstretched_a, parameters);
            h1_b1 = Calculateh1(A0_b_t1, Q0_b_t1, R_unstretched_b, parameters);
            h1_c1 = Calculateh1(A0_c_t1, Q0_c_t1, R_unstretched_c, parameters);
            
            f4 = -h1_aN + h1_b1*cos(alpha_b) + h1_c1*cos(alpha_c);
            F(equation_vector_index) = f4;                   
            equation_vector_index = equation_vector_index + 1;
 
            % Note: angles are measured in opposite directions in Whiteley
            % et al., whereas we define them in same direciton, so sign
            % difference
            f5 = h1_b1*sin(alpha_b) + h1_c1*sin(alpha_c);
            F(equation_vector_index) = f5;
            equation_vector_index = equation_vector_index + 1;
 
        else
            % Impose a boundary condition
            outlet_pressure = parameters.outlet_applied_pressure;
            G_0 = parameters.G_0;
            beta = parameters.beta;
            outlet_area = pi*R_unstretched_a^2*(outlet_pressure/G_0 + 1)^(2/beta);
            F(equation_vector_index) = outlet_area - AN_t1;
            equation_vector_index = equation_vector_index + 1;            
        end
    end
    
    % Finally, impose inlet bounday condition
    [A0_t1, Q0_t1, AN_t1, QN_t1] = GetVesselVariables(big_vector, 1);

    inlet_pressure = parameters.boundary_function(time, parameters, 1); 
    G_0 = parameters.G_0;
    beta = parameters.beta;
    R_unstretched = vessel_tree(1).R_unstretched;
    inlet_area = pi*R_unstretched^2*(inlet_pressure/G_0 + 1)^(2/beta);
    F(equation_vector_index) = inlet_area - A0_t1;

    equation_vector_index = equation_vector_index + 1;            

    if (equation_vector_index ~= 4*length(vessel_tree) + 1)
        error('Matrix is of incorrect size');
    end    
    
function f1 = Calculatef1(A_t0, Q_t0, A_t1, Q_t1, parameters)
    h = parameters.dx;
    dt = parameters.dt;
    A_int_diff_sum = sum(A_t1(2:end-1) - A_t0(2:end-1));

    f1 = (h/dt)*(0.5*(A_t1(1) - A_t0(1) + A_t1(end) - A_t0(end)) + A_int_diff_sum) ...
        - Q_t1(1) + Q_t1(end);
 
function f2 = Calculatef2(vessel, A0_t1, Q0_t1, AN_t1, QN_t1, parameters)
    h = parameters.dx;
    dt = parameters.dt;
    alpha = parameters.alpha;
    nu = parameters.nu;
    beta = parameters.beta;
    rho = parameters.rho;
    G_0 = parameters.G_0;
    
    Q0_t0 = vessel.Q(1);
    A0_t0 = vessel.A(1);
    QN_t0 = vessel.Q(end);
    AN_t0 = vessel.A(end);
    
    QA_int_sum = sum(vessel.Q_int_t1./vessel.A_int_t1);
    Q_int_diff_sum = sum(vessel.Q_int_t1 - vessel.Q(2:end-1));
    
    p0_t1 = PressureFromArea(A0_t1, vessel.R_unstretched, parameters);
    pN_t1 = PressureFromArea(AN_t1, vessel.R_unstretched, parameters);
    
    f2 = (h/dt)*(0.5*(Q0_t1 - Q0_t0 + QN_t1 - QN_t0) + Q_int_diff_sum) ...
        - alpha*Q0_t1^2/A0_t1 - (beta/(rho*(beta+2)))*(p0_t1 + G_0)*A0_t1 ...
        + alpha*QN_t1^2/AN_t1 + (beta/(rho*(beta+2)))*(pN_t1 + G_0)*AN_t1 ...
        + (pi*nu*alpha*h/(alpha-1))*(Q0_t1/A0_t1 + QN_t1/AN_t1 + 2*QA_int_sum);    
    
    
    
function h1 = Calculateh1(A, Q, R_unstretched, parameters)
    rho = parameters.rho;
    alpha = parameters.alpha;
    
    f_hat = PressureFromArea(A, R_unstretched, parameters);
    g = CalculateG(A, R_unstretched, parameters);
    
    h1 = rho*alpha*Q.^2./A + A.*f_hat - g;
    
function p = PressureFromArea(A, R_unstretched, parameters)
    p = parameters.G_0*((A./(pi*R_unstretched.^2)).^(parameters.beta/2) - 1);
    
function g = CalculateG(A, R_unstretched, parameters)
    beta = parameters.beta;
    G_0 = parameters.G_0;
    
    g = G_0*((A.^(beta/2+1))./(((pi*R_unstretched^2).^(beta/2)).*(beta/2+1)) - A + beta*pi*R_unstretched^2/(beta+2));

    
function [A0_index, Q0_index, AN_index, QN_index] = GetVariableIndexes(vessel_index)
    A0_index = (vessel_index - 1)*4 + 1;
    Q0_index = (vessel_index - 1)*4 + 2;
    AN_index = (vessel_index - 1)*4 + 3;
    QN_index = (vessel_index - 1)*4 + 4;

function [A0, Q0, AN, QN] = GetVesselVariables(big_vector, vessel_index)
    [A0_index, Q0_index, AN_index, QN_index] = GetVariableIndexes(vessel_index);
    A0 = big_vector(A0_index);
    Q0 = big_vector(Q0_index);
    AN = big_vector(AN_index);
    QN = big_vector(QN_index);  
    
