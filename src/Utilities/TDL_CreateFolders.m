function TDL_CreateFolders(path_to_create)
    % TDL_CreateFolders Create a directory if it does not exist
    %
    %     Author: Tom Doel www.tomdoel.com
    %     Part of TreeSolve. http://github.com/tomdoel/TreeSolve
    %     Distributed under the GNU GPL v3 licence. Please see LICENSE file.
    %    

    if ~(7 == exist(path_to_create, 'dir'))
        mkdir(path_to_create);
    end
end