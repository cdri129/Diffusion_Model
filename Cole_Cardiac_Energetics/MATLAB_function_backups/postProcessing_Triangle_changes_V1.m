% Define the path to the original file
original_file_path = '/home/cdri129/Desktop/Cole_Modelling/Triangle_Modelling/Mathias_Made/mesh/binary_2D2_2.poly';

% Read the contents of the original file
fid = fopen(original_file_path, 'r');
if fid == -1
    error('Failed to open the file: %s', original_file_path);
end

% Read the content
content = textscan(fid, '%d %d %d %d %d', 'delimiter', '\t');
fclose(fid);

% Modify the content to replace the 5th column with 0s
content{5} = zeros(size(content{5}));

% Define the path for the new file
new_file_path = '/home/cdri129/Desktop/Cole_Modelling/Triangle_Modelling/Mathias_Made/mesh/binary_2D2_2_bound0.poly';

% Write the modified content to the new file
fid = fopen(new_file_path, 'w');
if fid == -1
    error('Failed to open the file for writing: %s', new_file_path);
end

% Write the content to the new file
for i = 1:numel(content{1})
    fprintf(fid, '%d\t%d\t%d\t%d\t%d\n', content{1}(i), content{2}(i), content{3}(i), content{4}(i), content{5}(i));
end

fclose(fid);

disp(['File successfully created: ' new_file_path]);
