function R = TDL_CalculateRadiusFromPressure(p, R_unstretched, parameters)
    % TDL_CalculateRadiusFromPressure Helper function for TreeSolve
    %
    %
    %     Author: Tom Doel www.tomdoel.com
    %     Part of TreeSolve. http://github.com/tomdoel/TreeSolve
    %     Distributed under the GNU GPL v3 licence. Please see LICENSE file.
    % 

    % Derived from eqn. 2.36 in Smith et. al. 2002 - see eqn 2.46
    R = R_unstretched .*(((p/parameters.G_0) + 1).^(1/parameters.beta));
