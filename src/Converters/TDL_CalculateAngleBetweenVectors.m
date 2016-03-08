function angle = TDL_CalculateAngleBetweenVectors(vector_1, vector_2)
    % TDL_CalculateAngleBetweenVectors Helper function for TreeSolve
    %
    %
    %     Author: Tom Doel www.tomdoel.com
    %     Part of TreeSolve. http://github.com/tomdoel/TreeSolve
    %     Distributed under the GNU GPL v3 licence. Please see LICENSE file.
    %    

    % Returns angle in range 0-pi
    angle = atan2(norm(cross(vector_1, vector_2)),dot(vector_1, vector_2));
end
