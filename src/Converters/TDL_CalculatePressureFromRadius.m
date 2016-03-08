function p = TDL_CalculatePressureFromRadius(R, R_unstretched, parameters)
    % TDL_CalculatePressureFromRadius Helper function for TreeSolve
    %
    %
    %     Author: Tom Doel www.tomdoel.com
    %     Part of TreeSolve. http://github.com/tomdoel/TreeSolve
    %     Distributed under the GNU GPL v3 licence. Please see LICENSE file.
    % 
    
    beta = parameters.beta;
    G_0 = parameters.G_0;
    p = G_0.*(((R./R_unstretched).^beta) - 1);
end