function solution_dir = TDL_SaveVesselResults(file_prefix, vessel_tree, parameters)
    % TDL_SaveVesselResults Saves results from a vessel simulation
    %
    % Every set of results is saved in its own directory, to ensure all the
    % files are kept together. A unique folder is created on each run,
    % with the name Results/<dir_name>, where <dir_name> is the file_prefix plus 
    % a number to ensure the directory is unique.
    %
    % A .mat file is saved containing the vessel results, the parameters and
    % a version number to ensure old results aren't mixed with new results.
    %
    % Additionally, files suitable for viewing in CMGUi are generated from the
    % results.
    %
    %     Author: Tom Doel www.tomdoel.com
    %     Part of TreeSolve. http://github.com/tomdoel/TreeSolve
    %     Distributed under the GNU GPL v3 licence. Please see LICENSE file.
    %    



% This schema number should be incremented when the file format changes,
% or when using a newer version of the simulation code.
schema_number = 1; %#ok<NASGU>


if ~exist('Results', 'dir')
    mkdir Results;
end

% Get a unique name for the solution directory
suffix_number = 0;
solution_dir = [fullfile('Results', file_prefix) int2str(suffix_number)];
while (exist(solution_dir, 'dir'))
    suffix_number = suffix_number + 1;
    solution_dir = [fullfile('Results', file_prefix) int2str(suffix_number)];
end

% Create new solution directory
mkdir(solution_dir);

% Save the vessel results and parameter values into a .MAT file
save([fullfile(solution_dir, file_prefix) '.mat'], 'vessel_tree', 'parameters', 'schema_number');

% Export results to CMGUI format
TDL_OutputToCmgui(fullfile(solution_dir, file_prefix), vessel_tree);
