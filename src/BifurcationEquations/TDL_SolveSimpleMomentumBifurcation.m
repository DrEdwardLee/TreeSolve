function [new_parent_point, new_child_point_1, new_child_point_2] = TDL_SolveSimpleMomentumBifurcation(parent_vessel, child_vessel_1, child_vessel_2, parameters)
    % TDL_SolveSimpleMomentumBifurcation Funtion for use with TDL_SolveVesselTree to solve solve bifurcation equations using momentum conservation equations
    %
    %     Author: Tom Doel www.tomdoel.com
    %     Part of TreeSolve. http://github.com/tomdoel/TreeSolve
    %     Distributed under the GNU GPL v3 licence. Please see LICENSE file.
    %    
    
      
    % Starting guess: use previous values
    p_vector_guess = [parent_vessel.p(end); child_vessel_1.p(1); child_vessel_2.p(1)];

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
end    

function rhs = NewtonRHS2(p_vector_t1, parent_vessel_t0, child1_vessel_t0, child2_vessel_t0, parameters)
    
    parent_point_t1 = TDL_VesselBoundaryEndFromPressure(p_vector_t1(1), parent_vessel_t0, parameters);
    child1_point_t1 = TDL_VesselBoundaryStartFromPressure(p_vector_t1(2), child1_vessel_t0, parameters);
    child2_point_t1 = TDL_VesselBoundaryStartFromPressure(p_vector_t1(3), child2_vessel_t0, parameters);

    F_a1_t1 = Flow(parent_point_t1);
    F_b1_t1 = Flow(child1_point_t1);
    F_c1_t1 = Flow(child2_point_t1);
    
    h1_a1 = parent_point_t1.p * parent_point_t1.R^2;
    h1_b1 = child1_point_t1.p * child1_point_t1.R^2;
    h1_c1 = child2_point_t1.p * child2_point_t1.R^2;
    
    alpha_2 = child1_vessel_t0.angle;
    alpha_3 = child2_vessel_t0.angle;
    
    rhs = zeros(3, 1);
    rhs(1) = h1_a1 - h1_b1*cos(alpha_2) - h1_c1*cos(alpha_3);
    rhs(2) = h1_b1*sin(alpha_2) + h1_c1*sin(alpha_3);
    rhs(3) = - F_a1_t1 + F_b1_t1 + F_c1_t1;
end
 
function F = Flow(point)
   F = pi*(point.R.^2).*point.V;
end

    
    
    
  
