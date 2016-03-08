function flows = TDL_GetBoundaryFlows(vessel_tree)
    % TDL_GetBoundaryFlows Helper function for TreeSolve
    %
    %
    %     Author: Tom Doel www.tomdoel.com
    %     Part of TreeSolve. http://github.com/tomdoel/TreeSolve
    %     Distributed under the GNU GPL v3 licence. Please see LICENSE file.
    % 
    
    flows = zeros(1,3);
    Va = vessel_tree(1).V(1);
    Ra = vessel_tree(1).R(1);
    flows(1) = Va*pi*Ra^2;

    Vb = vessel_tree(2).V(end);
    Rb = vessel_tree(2).R(end);
    flows(2) = Vb*pi*Rb^2;

    Vc = vessel_tree(3).V(end);
    Rc = vessel_tree(3).R(end);
    flows(3) = Vc*pi*Rc^2;