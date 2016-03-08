function maximize(fig)
    % maximize Displays a figure full-screen
    %
    %     Author: Tom Doel www.tomdoel.com
    %     Part of TreeSolve. http://github.com/tomdoel/TreeSolve
    %     Distributed under the GNU GPL v3 licence. Please see LICENSE file.
    %    
    
if nargin < 1
    fig = gcf;
end

old_units = get(fig, 'units');
set(fig, 'units', 'normalized', 'outerposition', [0, 0, 1, 1]);
set(fig, 'units', old_units);
