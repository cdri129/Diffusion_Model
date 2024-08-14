imageName = 'C001';
base_dir = join([getenv("HOME"), "Desktop/Cole_Modelling/Triangle_Modelling/Mathias_Made"], '/');
img_dir = join([base_dir, "img"], '/');
mesh_dir = join([base_dir, "mesh"], '/');


% Define input and output file names
inputFileName = fullfile(mesh_dir, [imageName '.poly']);
outputFileName = fullfile(mesh_dir, [imageName '_BC.poly']);
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
    % Split the line into columns
    columns = strsplit(line);
    
    % Check if the line represents a segment and the boundary marker is non-zero
    if numel(columns) == 5 && str2double(columns{5}) ~= 0
        % Change the boundary marker to 0
        columns{5} = '0';
        
        % Join the columns back into a line
        modifiedLine = strjoin(columns, ' ');
        
        % Write modified line to output file
        fprintf(outputFile, '%s\n', modifiedLine);
    else
        % Write unmodified line to output file
        fprintf(outputFile, '%s\n', line);
    end
    
    % Read next line
    line = fgetl(inputFile);
end

% Close files
fclose(inputFile);
fclose(outputFile);

disp('Boundary conditions modified and written to output file.');