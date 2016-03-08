function vessel_tree = TDL_SolveVesselTreeWithFlowConservation(vessel_tree, parameters)
    % TDL_SolveVesselTreeWithFlowConservation Performs a fluid flow
    % simulation using a momentum conservation scheme that takes account of
    % branch angles
    %
    % Implementation of a scheme proposed by Jonathan Whiteley, Kelly
    % Burrowes, Kathryn Gillow, David Gavaghan (2008), University of Oxford.
    %
    % Please note: this implementation is experimental. If you have
    % difficulties, consider using TDL_SolveVesselTree.
    %
    % Use the parameter parameters.use_newton_solver to choose the solver
    %
    %     Author: Tom Doel www.tomdoel.com
    %     Part of TreeSolve. http://github.com/tomdoel/TreeSolve
    %     Distributed under the GNU GPL v3 licence. Please see LICENSE file.
    %    

    vessel_tree_1 = TDL_ConvertVesselTreeToAQ(vessel_tree);
    SolveVesselTree(vessel_tree_1, parameters);
    vessel_tree = TD_ConvertVesselTreeFromAQ(vessel_tree_1, parameters);

    
function vessel_tree = SolveVesselTree(vessel_tree, parameters)

    max_time = parameters.max_time;
    time = 0;
    skip = 0;
    steady_state_count = 0;
    
    steady_state_to_match = vessel_tree;

    while time < max_time
        vessel_tree = NextTimeStep(vessel_tree, parameters, time);
        closeness = VesselTreeValuesMatch(steady_state_to_match, vessel_tree);
        if (closeness < 0.001)
            steady_state_count = steady_state_count + 1;
            if (steady_state_count > 1000)
                disp('Steady state reached');
                return;
            end
        else
            steady_state_count = 0;
            steady_state_to_match = vessel_tree;
        end
        time = time + parameters.dt;

        if (mod(skip,10) == 0)
            tree_to_plot = TD_ConvertVesselTreeFromAQ(vessel_tree, parameters);
            for vessel_number = 1:length(vessel_tree)
                TDL_PlotVessel(tree_to_plot, vessel_number, parameters, 'b');
            end
            pause(0.001);
        end
        skip = skip + 1;
    end


function closeness = VesselTreeValuesMatch(tree_1, tree_2)
    for (vessel_index = 1 : length(tree_1))
        norm_p = norm(tree_1(vessel_index).A - tree_2(vessel_index).A);
        norm_R = norm(tree_1(vessel_index).Q - tree_2(vessel_index).Q);
    end
    closeness = min(norm_p, norm_R);

function vessel_tree_t1 = NextTimeStep(vessel_tree_t0, parameters, time)

    vessel_tree_intermediate = vessel_tree_t0;
    num_vessels = length(vessel_tree_t0);
    
    
    % Calculate intermediate and internal points for all vessels; these are required before
    % we solve the equations at the bifurcations
    for vessel_index = 1 : num_vessels
        vessel = vessel_tree_t0(vessel_index);
        R_unstretched_a = vessel.R_unstretched;

        A_t0 = vessel.A;
        Q_t0 = vessel.Q;
           
        [A_half Q_half] = LWStep1(A_t0, Q_t0, R_unstretched_a, parameters);
        [A_t1_internal, Q_t1_internal] = LWStep2(A_t0, Q_t0, A_half, Q_half, R_unstretched_a, parameters);
       
        % These are values at points t+dt/2, x+dx/2
        vessel_tree_intermediate(vessel_index).A_half = A_half;
        vessel_tree_intermediate(vessel_index).Q_half = Q_half;
        
        % The values at t+dt for internal points in the vessels
        vessel_tree_intermediate(vessel_index).Q_int_t1 = Q_t1_internal;
        vessel_tree_intermediate(vessel_index).A_int_t1 = A_t1_internal;      
    end

    if (parameters.use_newton_solver)
        vessel_tree_intermediate = SolveBifurcationsOld(vessel_tree_t0, vessel_tree_intermediate, parameters, time);
    else
        vessel_tree_intermediate = SolveBifurcations(vessel_tree_t0, vessel_tree_intermediate, parameters, time);
    end
    vessel_tree_t1 = UpdateVesselTree(vessel_tree_t0, vessel_tree_intermediate);


function vessel_tree_t1 = UpdateVesselTree(vessel_tree_t0, vessel_tree_intermediate)
    % Update vessel tree with the new values calculated above
    vessel_tree_t1 = vessel_tree_t0;
    for vessel_index = 1:length(vessel_tree_t0)
        % Update the internal points
        vessel_tree_t1(vessel_index).Q(2:end-1) = vessel_tree_intermediate(vessel_index).Q_int_t1;
        vessel_tree_t1(vessel_index).A(2:end-1) = vessel_tree_intermediate(vessel_index).A_int_t1;
        
        % Update the boundary points
        vessel_tree_t1(vessel_index).Q(1) = vessel_tree_intermediate(vessel_index).firstpoint.Q;
        vessel_tree_t1(vessel_index).A(1) = vessel_tree_intermediate(vessel_index).firstpoint.A;
    
        vessel_tree_t1(vessel_index).Q(end) = vessel_tree_intermediate(vessel_index).lastpoint.Q;
        vessel_tree_t1(vessel_index).A(end) = vessel_tree_intermediate(vessel_index).lastpoint.A;
    end

    
function vessel_tree_intermediate = SolveBifurcations(vessel_tree_t0, vessel_tree_intermediate, parameters, time)
    % Starting guess: use previous values
    big_vector = zeros(4*length(vessel_tree_t0), 1);
    for vessel_index = 1 : length(vessel_tree_t0)
        vessel_t0 = vessel_tree_t0(vessel_index);
        [A0_index, Q0_index, AN_index, QN_index] = GetVariableIndexes(vessel_index);
        big_vector(A0_index) = vessel_t0.A(1);
        big_vector(Q0_index) = vessel_t0.Q(1);
        big_vector(AN_index) = vessel_t0.A(end);
        big_vector(QN_index) = vessel_t0.Q(end);
    end
    
    options = optimset('Display', 'off');
    optiond = optimset;
    
    [big_vector_solved, fval,exitflag,output] = fsolve(@TDL_FlowConservationFunction, big_vector, options, vessel_tree_intermediate, parameters, time);
    
    if (exitflag < 1)
        disp('*** Solver failed');
        fval
        exitflag
        output
        error('Solver failed');
    else
       disp('Solver succeeded'); 
    end
        
    big_vector = big_vector_solved;
    
     % Write out results into vessel tree
    for vessel_index = 1 : length(vessel_tree_intermediate)
        [A0, Q0, AN, QN] = GetVesselVariables(big_vector, vessel_index);
        firstpoint.Q = Q0;
        firstpoint.A = A0;
        lastpoint.Q = QN;
        lastpoint.A = AN;
        vessel_tree_intermediate(vessel_index).firstpoint = firstpoint;
        vessel_tree_intermediate(vessel_index).lastpoint = lastpoint;
    end


function vessel_tree = SolveBifurcationsOld(vessel_tree_t0, vessel_tree, parameters, time)

    % Starting guess: use previous values
    big_vector = zeros(4*length(vessel_tree), 1);
    for vessel_index = 1 : length(vessel_tree)
        vessel = vessel_tree(vessel_index);
        [A0_index, Q0_index, AN_index, QN_index] = GetVariableIndexes(vessel_index);
        big_vector(A0_index) = vessel.A(1);
        big_vector(Q0_index) = vessel.Q(1);
        big_vector(AN_index) = vessel.A(end);
        big_vector(QN_index) = vessel.Q(end);
    end
    
    iteration = 0;
    tol = 0.1;

    while (true)
        [M b] = CalculateNewtonMatrixVector(vessel_tree, big_vector, parameters, time);
        s = M\b;
        
        if (imag(s) ~= 0)
            disp('Warning: s produced imaginary result');
           % s = real(s);
        end
        
        big_vector = big_vector + 0.01*s;
        step_norm = norm(s./big_vector);
        
        if (iteration > 100)
            error(['Warning: solution did not converge; iteration: ' num2str(iteration) ', norm:', num2str(step_norm)]);
            break
        end

        if (step_norm < tol)
            break
        end

        iteration = iteration + 1;  
    end
    
    % Write out results into vessel tree
    for vessel_index = 1 : length(vessel_tree)
        [A0, Q0, AN, QN] = GetVesselVariables(big_vector, vessel_index);
        firstpoint.Q = Q0;
        firstpoint.A = A0;
        lastpoint.Q = QN;
        lastpoint.A = AN;
        vessel_tree(vessel_index).firstpoint = firstpoint;
        vessel_tree(vessel_index).lastpoint = lastpoint;
    end
    
        

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
    

function [M b] = CalculateNewtonMatrixVector(vessel_tree, big_vector, parameters, time)
    M = zeros(length(big_vector));
    b = zeros(size(big_vector));
    
    h = parameters.dx;
    dt = parameters.dt;
    
    equation_vector_index = 1;
    num_vessels = length(vessel_tree);
    
    for vessel_index = 1 : num_vessels        
        A_t0 = vessel_tree(vessel_index).A;
        Q_t0 = vessel_tree(vessel_index).Q;
        
        A0_t0 = A_t0(1);
        Q0_t0 = Q_t0(1);
        AN_t0 = A_t0(end);
        QN_t0 = Q_t0(end);
        
        R_unstretched_a = vessel_tree(vessel_index).R_unstretched;
                
        [A0_index, Q0_index, AN_index, QN_index] = GetVariableIndexes(vessel_index);
        [A0_t1, Q0_t1, AN_t1, QN_t1] = GetVesselVariables(big_vector, vessel_index);
        
        A_t1_internal = vessel_tree(vessel_index).A_int_t1;
        Q_t1_internal = vessel_tree(vessel_index).Q_int_t1;
        A_t1 = [A0_t1; A_t1_internal; AN_t1];
        Q_t1 = [Q0_t1; Q_t1_internal; QN_t1];

        f1 = Calculatef1(A_t0, Q_t0, A_t1, Q_t1, parameters);
        
        df1dA0 = h/(2*dt);
        df1dQ0 = -1;
        df1dAN = h/(2*dt);
        df1dQN = 1;
        
        temp_mul_factor = 1;%1/1000;
        
        b(equation_vector_index) = - temp_mul_factor*f1;
        M(equation_vector_index, A0_index) =  temp_mul_factor*df1dA0;
        M(equation_vector_index, Q0_index) =  temp_mul_factor*df1dQ0;
        M(equation_vector_index, AN_index) =  temp_mul_factor*df1dAN;
        M(equation_vector_index, QN_index) =  temp_mul_factor*df1dQN;
        
        equation_vector_index = equation_vector_index + 1;
      
        f2 = Calculatef2(A_t0, Q_t0, A_t1, Q_t1, R_unstretched_a, parameters);
        
        temp_mul_factor = 1;%/10000;
        
        b(equation_vector_index) = - temp_mul_factor*f2;
        M(equation_vector_index, A0_index) = temp_mul_factor*Calculatedf2dA0(A_t0, Q_t0, A_t1, Q_t1, R_unstretched_a, parameters);
        M(equation_vector_index, Q0_index) = temp_mul_factor*Calculatedf2dQ0(A_t0, Q_t0, A_t1, Q_t1, R_unstretched_a, parameters);
        M(equation_vector_index, AN_index) = temp_mul_factor*Calculatedf2dAN(A_t0, Q_t0, A_t1, Q_t1, R_unstretched_a, parameters);
        M(equation_vector_index, QN_index) = temp_mul_factor*Calculatedf2dQN(A_t0, Q_t0, A_t1, Q_t1, R_unstretched_a, parameters);
        
        
        equation_vector_index = equation_vector_index + 1;
 
        
        if (~isempty(vessel_tree(vessel_index).connected_to))
            vessel_index_child_1 = vessel_tree(vessel_index).connected_to(1);
            vessel_index_child_2 = vessel_tree(vessel_index).connected_to(2);
            
            R_unstretched_b = vessel_tree(vessel_index_child_1).R_unstretched;
            R_unstretched_c = vessel_tree(vessel_index_child_2).R_unstretched;
            
            [bA0_index, bQ0_index, bAN_index, bQN_index] = GetVariableIndexes(vessel_index_child_1);
            [A0_b_t1, Q0_b_t1, AN_b_t1, QN_b_t1] = GetVesselVariables(big_vector, vessel_index_child_1);

            [cA0_index, cQ0_index, cAN_index, cQN_index] = GetVariableIndexes(vessel_index_child_2);
            [A0_c_t1, Q0_c_t1, AN_c_t1, QN_c_t1] = GetVesselVariables(big_vector, vessel_index_child_2);
        
            
            f3 = Q_t1(end) - Q0_b_t1 - Q0_c_t1;
            
            b(equation_vector_index) = -f3;
            M(equation_vector_index, QN_index) = 1;
            M(equation_vector_index, bQ0_index) = -1;
            M(equation_vector_index, cQ0_index) = -1;
        
            equation_vector_index = equation_vector_index + 1;
 
            h1_aN = Calculateh1(AN_t1, QN_t1, R_unstretched_a, parameters);
            h1_b1 = Calculateh1(A0_c_t1, Q0_c_t1, R_unstretched_b, parameters);
            h1_c1 = Calculateh1(A0_c_t1, Q0_c_t1, R_unstretched_c, parameters);
            
            alpha_2 = vessel_tree(vessel_index_child_1).angle;
            alpha_3 = vessel_tree(vessel_index_child_2).angle;
            
            
            f4 = - h1_aN + h1_b1*cos(alpha_2) + h1_c1*cos(alpha_3);
            
            dh1_dA_aN = CalculateDh1dA(AN_t1, QN_t1, R_unstretched_a, parameters);
            dh1_dQ_aN = CalculateDh1dQ(QN_t1, QN_t1, parameters);

            dh1_dA_b1 = CalculateDh1dA(A0_b_t1, Q0_b_t1, R_unstretched_b, parameters);
            dh1_dQ_b1 = CalculateDh1dQ(A0_b_t1, Q0_b_t1, parameters);

            dh1_dA_c1 = CalculateDh1dA(A0_c_t1, Q0_c_t1, R_unstretched_c, parameters);
            dh1_dQ_c1 = CalculateDh1dQ(A0_c_t1, Q0_c_t1, parameters);

            b(equation_vector_index) = - f4;
            M(equation_vector_index, AN_index) = dh1_dA_aN;
            M(equation_vector_index, QN_index) = dh1_dQ_aN;
            
            M(equation_vector_index, bA0_index) = -dh1_dA_b1*cos(alpha_2);
            M(equation_vector_index, bQ0_index) = -dh1_dQ_b1*cos(alpha_2);

            M(equation_vector_index, cA0_index) = -dh1_dA_c1*cos(alpha_3);
            M(equation_vector_index, cQ0_index) = -dh1_dQ_c1*cos(alpha_3);
            
            
                   
            equation_vector_index = equation_vector_index + 1;
 
            f5 = h1_b1*sin(alpha_2) + h1_c1*sin(alpha_3); % Note sign change
                          
            
            b(equation_vector_index) = - f5;
                        
            M(equation_vector_index, bA0_index) = dh1_dA_b1*sin(alpha_2);
            M(equation_vector_index, bQ0_index) = dh1_dQ_b1*sin(alpha_2);

            M(equation_vector_index, cA0_index) = -dh1_dA_c1*sin(alpha_3);
            M(equation_vector_index, cQ0_index) = -dh1_dQ_c1*sin(alpha_3);
            
            
            equation_vector_index = equation_vector_index + 1;
 
        else
            % Impose a boundary condition
            outlet_pressure = parameters.outlet_applied_pressure;
            G_0 = parameters.G_0;
            beta = parameters.beta;
            outlet_area = pi*R_unstretched_a^2*(outlet_pressure/G_0 + 1)^(2/beta);
            b(equation_vector_index) = outlet_area - AN_t1;
            
            M(equation_vector_index, AN_index) = 1;
           
            equation_vector_index = equation_vector_index + 1;            
        end
    end
    
    % Finally, impose inlet bounday condition
    [A0_index, Q0_index, AN_index, QN_index] = GetVariableIndexes(1);
    [A0_t1, Q0_t1, AN_t1, QN_t1] = GetVesselVariables(big_vector, 1);

    inlet_pressure = parameters.boundary_function(time, parameters, 1); %parameters.inlet_applied_pressure;
    G_0 = parameters.G_0;
    beta = parameters.beta;
    R_unstretched = vessel_tree(1).R_unstretched;
    inlet_area = pi*R_unstretched^2*(inlet_pressure/G_0 + 1)^(2/beta);
    b(equation_vector_index) = inlet_area - A0_t1;

    M(equation_vector_index, A0_index) = 1;

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
 
function f2 = Calculatef2(A_t0, Q_t0, A_t1, Q_t1, R_unstretched, parameters)
    h = parameters.dx;
    dt = parameters.dt;
    alpha = parameters.alpha;
    nu = parameters.nu;
    beta = parameters.beta;
    rho = parameters.rho;
    G_0 = parameters.G_0;
    
    Q_int_diff_sum = sum(Q_t1(2:end-1) - Q_t0(2:end-1));
    QA_int_sum = sum(Q_t1(2:end-1)./A_t1(2:end-1));

    
    p_t1_x0 = CalculateFhat(A_t1(1), R_unstretched, parameters);
    p_t1_xN = CalculateFhat(A_t1(end), R_unstretched, parameters);
    
    f2 = (h/dt)*(0.5*(Q_t1(1) - Q_t0(1) + Q_t1(end) - Q_t0(end)) + Q_int_diff_sum) ...
        - alpha*Q_t1(1)^2/A_t1(1) - (beta/(rho*(beta+2)))*(p_t1_x0 + G_0)*A_t1(1) ...
        + alpha*Q_t1(end)^2/A_t1(end) + (beta/(rho*(beta+2)))*(p_t1_xN + G_0)*A_t1(end) ...
        + (pi*nu*alpha*h/(alpha-1))*(Q_t1(1)/A_t1(1) + Q_t1(end)/A_t1(end) + 2*QA_int_sum);    
    
    
function df2dA0 = Calculatedf2dA0(A_t0, Q_t0, A_t1, Q_t1, R_unstretched, parameters)
    h = parameters.dx;
    alpha = parameters.alpha;
    nu = parameters.nu;
    beta = parameters.beta;
    rho = parameters.rho;
    G_0 = parameters.G_0;
    
    p_t1_x0 = CalculateFhat(A_t1(1), R_unstretched, parameters);

    df2dA0 = alpha*((Q_t1(1)/A_t1(1))^2) - (beta/(2*rho))*(p_t1_x0 + G_0) ...
        - (pi*nu*alpha*h/(alpha-1))*Q_t1(1)/(A_t1(1)^2);  
    

function df2dQ0 = Calculatedf2dQ0(A_t0, Q_t0, A_t1, Q_t1, R_unstretched, parameters)
    h = parameters.dx;
    dt = parameters.dt;
    alpha = parameters.alpha;
    nu = parameters.nu;
    
    df2dQ0 = h/(2*dt) - 2*alpha*Q_t1(1)/A_t1(1) + pi*nu*alpha*h/((alpha-1)*A_t1(1));
    
    
    
function df2dAN = Calculatedf2dAN(A_t0, Q_t0, A_t1, Q_t1, R_unstretched, parameters)
    h = parameters.dx;
    alpha = parameters.alpha;
    nu = parameters.nu;
    beta = parameters.beta;
    rho = parameters.rho;
    G_0 = parameters.G_0;
    
    p_t1_xN = CalculateFhat(A_t1(end), R_unstretched, parameters);
    
    df2dAN = - alpha*((Q_t1(end)/A_t1(end))^2) + (beta/(2*rho))*(p_t1_xN + G_0) ...
        - (pi*nu*alpha*h/(alpha-1))*Q_t1(end)/(A_t1(end)^2);  
    
    
function df2dQN = Calculatedf2dQN(A_t0, Q_t0, A_t1, Q_t1, R_unstretched, parameters)
    h = parameters.dx;
    dt = parameters.dt;
    alpha = parameters.alpha;
    nu = parameters.nu;
    
    df2dQN = h/(2*dt) + 2*alpha*Q_t1(end)/A_t1(end) + pi*nu*alpha*h/((alpha-1)*A_t1(end));
        
function [A_thalf Q_thalf] = LWStep1(A_t0, Q_t0, R_unstretched, parameters)
    dt = parameters.dt;
    dx = parameters.dx;
    nu = parameters.nu;
    alpha = parameters.alpha;
    rho = parameters.rho;
    
    A_i0 = A_t0(1:end-1);
    Q_i0 = Q_t0(1:end-1);
    A_i1 = A_t0(2:end);
    Q_i1 = Q_t0(2:end);
    
    h1_i0 = Calculateh1(A_i0, Q_i0, R_unstretched, parameters);
    h1_i1 = Calculateh1(A_i1, Q_i1, R_unstretched, parameters);
    
    A_thalf = (1/2)*(A_t0(1:end-1) + A_t0(2:end)) - (dt/(2*dx))*(Q_t0(2:end) - Q_t0(1:end-1));

    Q_thalf = (1/2)*(Q_i0 + Q_i1) - (dt/(2*rho*dx))*(h1_i1 - h1_i0) ...
        - pi*nu*alpha*dt/(alpha-1)*(Q_i0 + Q_i1)./(A_i0 + A_i1);

    
function [A_t1_internal, Q_t1_internal] = LWStep2(A_t0, Q_t0, A_half, Q_half, R_unstretched, parameters)
    dt = parameters.dt;
    dx = parameters.dx;
    nu = parameters.nu;
    alpha = parameters.alpha;
    rho = parameters.rho;
    
    h1_minus = Calculateh1(A_half(1:end-1), Q_half(1:end-1), R_unstretched, parameters);
    h1_plus  = Calculateh1(A_half(2:end), Q_half(2:end), R_unstretched, parameters);
    
    A_t1_internal = A_t0(2:end-1) - (dt/dx)*(Q_half(2:end) - Q_half(1:end-1));

    Q_t1_internal = Q_t0(2:end-1) - (dt/(rho*dx))*(h1_plus - h1_minus) ...
        - 2*pi*nu*alpha*dt/(alpha - 1)*Q_t0(2:end-1)./A_t0(2:end-1);
    
function h1 = Calculateh1(A, Q, R_unstretched, parameters)
    rho = parameters.rho;
    alpha = parameters.alpha;
    
    f_hat = CalculateFhat(A, R_unstretched, parameters);
    g = CalculateG(A, R_unstretched, parameters);
    
    h1 = rho*alpha*Q.^2./A + A.*f_hat - g;

    
function dh1dQ = CalculateDh1dQ(A, Q, parameters)
    dh1dQ = 2*parameters.rho*parameters.alpha*Q./A;    

function dh1dA = CalculateDh1dA(A, Q, R_unstretched, parameters)
    rho = parameters.rho;
    alpha = parameters.alpha;
    beta = parameters.beta;
    G_0 = parameters.G_0;
    
    dh1dA = -rho*alpha*(Q./A).^2 + (beta/2)*G_0*(A/(pi*R_unstretched^2)).^(beta/2);
    
    
function f_hat = CalculateFhat(A, R_unstretched, parameters)
    f_hat = parameters.G_0*((A./(pi*R_unstretched.^2)).^(parameters.beta/2) - 1);
    
function dpdA = CalculatedpdA(A, R_unstretched, parameters)
    beta = parameters.beta;
    G_0 = parameters.G_0;
    
    dpdA = (G_0*beta/2)*((A.^(beta/2-1))./((pi*R_unstretched.^2).^(beta/2)));


function g = CalculateG(A, R_unstretched, parameters)
    beta = parameters.beta;
    G_0 = parameters.G_0;
    
    g = G_0*((A.^(beta/2+1))./(((pi*R_unstretched^2).^(beta/2)).*(beta/2+1)) - A + beta*pi*R_unstretched^2/(beta+2));
