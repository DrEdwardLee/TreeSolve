function handle = TDL_subplot(rows, cols, number)
    % TDL_subplot customised subplot
    %
    %     Author: Tom Doel www.tomdoel.com
    %     Part of TreeSolve. http://github.com/tomdoel/TreeSolve
    %     Distributed under the GNU GPL v3 licence. Please see LICENSE file.
    %    

    xp = mod(number - 1, cols);
    yp = fix((number - 1)/cols);

    % Horizontally, graphs are placed together with spacing around the
    % borders, but vertically we also have space between graphs
    x_spacing = 0.05;
    width = (1.0 - 2*x_spacing)/cols;

    y_spacing = 0.05;
    height = (1.0 - y_spacing)/rows;
    plot_height = height - y_spacing;
    
    xc = xp*width + x_spacing;
    yc = (rows - 1 - yp)*height + y_spacing;
    
    handle = subplot('position', [xc, yc, width, plot_height]);
end

