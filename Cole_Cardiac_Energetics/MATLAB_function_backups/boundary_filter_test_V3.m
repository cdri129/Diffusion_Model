imageName = 'C001';
base_dir = fullfile(getenv("HOME"), "Desktop/Cole_Modelling/Triangle_Modelling/Mathias_Made");

img_dir = fullfile(base_dir, "img");
mesh_dir = fullfile(base_dir, "mesh");

disp("Processing image " + imageName); 

desired_border_type = 2;

original_node_file = strcat(imageName, '.1.node');
original_file_path = fullfile(mesh_dir, imageName, original_node_file);

% Check if the original node file exists
if exist(original_file_path, 'file') ~= 2
    error('Original node file does not exist.');
end

fid = fopen(original_file_path, 'r');

if fid == -1
    error('Failed to open the original node file for reading.');
end

% Read the contents of the original node file
nodes = textscan(fid, '%d %f %f %d %d', 'HeaderLines', 1);
fclose(fid);

% Extract node attributes
node_ids = nodes{1};
x_coords = nodes{2};
y_coords = nodes{3};
mitochondria_status = nodes{4};
border_type = nodes{5};

% Filter nodes based on criteria
filtered_nodes = node_ids(mitochondria_status == 1 & border_type == desired_border_type);
filtered_x_coords = x_coords(mitochondria_status == 1 & border_type == desired_border_type);
filtered_y_coords = y_coords(mitochondria_status == 1 & border_type == desired_border_type);

% Write filtered nodes to a new node file
new_node_file = strcat('filtered_border_', num2str(desired_border_type),'_', imageName, '.node');

new_node_file_path = fullfile(mesh_dir, new_node_file);  % Full path of the new file

fid = fopen(new_node_file_path, 'w');
if fid == -1
    error('Failed to open the new node file for writing.');
end
% Write the number of filtered nodes
fprintf(fid, '%d\n', numel(filtered_nodes));

% Write each filtered node to the new node file
for i = 1:numel(filtered_nodes)
    fprintf(fid, '%d\t%f\t%f\t%d\t%d\n', filtered_nodes(i), filtered_x_coords(i), filtered_y_coords(i), 1, desired_border_type);
end

fclose(fid);

disp('Filtered node file has been successfully created at:');
disp(new_node_file_path);  % Display the full path of the new file