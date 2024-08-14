% Load data from the .poly file
file_path = '/home/cdri129/Desktop/Cole_Modelling/Triangle_Modelling/Mathias_Made/mesh/C001.poly';
data = importdata(file_path);

% Extract points from the data
points = data.data(:, 2:3);  % Assuming the x and y coordinates are in columns 2 and 3

% Plot the points
scatter(points(:, 1), points(:, 2), '.');  % Scatter plot of points
xlabel('X');
ylabel('Y');
title('Plot of points from .poly file');
axis equal;  % Set equal aspect ratio
