imageName = 'C001';
base_dir = join([getenv("HOME"), "Desktop/Cole_Modelling/Triangle_Modelling/Mathias_Made"], '/');
img_dir = join([base_dir, "img"], '/');
mesh_dir = join([base_dir, "mesh"], '/');


% Define input and output file names
inputFileName = fullfile(mesh_dir, [imageName '_BC.poly']);
outputFileName = fullfile(mesh_dir, [imageName '_BC_2.poly']);
% Open input .poly file for reading
inputFile = fopen(inputFileName, 'r');
if inputFile == -1
    error('Failed to open input .poly file');
end

% Open output file for writing
outputFile = fopen(outputFileName, 'w');
if outputFile == -1
    fclose(inputFile);
    error('Failed to create output file');
end

% Read lines from input .poly file
line = fgetl(inputFile);
while ischar(line)
    % Check if the line represents the start of the segment list section
    if numel(strsplit(line)) == 2
        % Write the line as is to the output file
        fprintf(outputFile, '%s\n', line);
        
        % Move to the next line
        line = fgetl(inputFile);
        
        % Process lines until the end of the point list section
while ischar(line)
    % Check if the line represents the start of the segment list section
    if numel(strsplit(line)) == 2
        % Write the line as is to the output file
        fprintf(outputFile, '%s\n', line);
        
        % Move to the next line
        line = fgetl(inputFile);
        
        % Process lines until the end of the segment list section
        while ischar(line) && numel(strsplit(line)) == 4
            % Split the line into an array of strings
            lineSplit = strsplit(line, ' ');
        
            % Replace the value in the fourth column with '0'
            line = strrep(line, [' ' strtrim(lineSplit{4})], ' 0');
        
            % Write modified line to output file
            fprintf(outputFile, '%s\n', line);
        
            % Move to the next line
            line = fgetl(inputFile);
        end
    end
    
    % Write the line to the output file (including lines outside the segment list section)
    if ischar(line)
        fprintf(outputFile, '%s\n', line);
    end
    
    % Read next line
    line = fgetl(inputFile);
end
    end
    
    % Write the line to the output file (including lines outside the segment list section)
    if ischar(line)
        fprintf(outputFile, '%s\n', line);
    end
    
    % Read next line
    line = fgetl(inputFile);
end

% Close files
fclose(inputFile);
fclose(outputFile);

disp('Boundary conditions modified and written to output file.');