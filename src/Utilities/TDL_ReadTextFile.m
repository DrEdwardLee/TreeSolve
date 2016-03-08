function file_contents = TDL_ReadTextFile(file_name)
    % TDL_ReadTextFile
    %
    % Returns an array of strings from the text file, one for each line,
    % with the leading and trailing white space and line breaks removed 
    % from each line
    %
    %     Author: Tom Doel www.tomdoel.com
    %     Part of TreeSolve. http://github.com/tomdoel/TreeSolve
    %     Distributed under the GNU GPL v3 licence. Please see LICENSE file.
    %    


    file_contents = {};    
    fid = fopen(file_name, 'rt');
    
    next_line = fgetl(fid);
    while (ischar(next_line))% ~= -1)
        file_contents{end + 1} = strtrim(next_line);
        next_line = fgetl(fid);
    end
        
    fclose(fid);
end

