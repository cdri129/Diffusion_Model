clear all
% Copied from other script (not needed if running without clear all)

imageName = 'CA002';
%imageName = 'binary_2D2_2';
%imageName = 'small_blocky_mito_myo';
%imageName = 'offset_glyco_centered_mito';
%imageName = 'small_blocky_other';
%imageName = 'myo_far_glyco';
%imageName = 'far_glyco_centered_mito';
%imageName = 'myo_centered_mito';
%imageName = 'myo_far_mito';
%imageName = 'test_far_myo_far_glyco_mid';
%imageName = 'test_touch_far_myo_far_glyco';

base_dir = join([getenv("HOME"), "Desktop/Cole_Modelling/Triangle_Modelling/Mathias_Made"], '/');
img_dir = join([base_dir, "img"], '/');
mesh_dir = join([base_dir, "mesh"], '/');
%%%% MROE734 EDITS add semicolon to imageName and a disp statement
disp("Processing image " + imageName); 

file_name = strcat(imageName,'.1.csv');
file_path = join([mesh_dir, file_name], '/');
%%%% MROE734 EDITS using readtable instead of csvread
fnode = readtable(file_path);
Y = fnode.Var2;
X = fnode.Var3;
%%%% MROE734 EDITS picked up labels as well
labels = fnode.Var4;

%%%% MROE734 EDITS normalised coordinates
% Normalise coordinates for plotting between 0 and 1
% Y = Y./max(Y);
% X = X./max(X);

%%%% ColeD EDITS
% added 'figure;' in an attempt to fix the plotting issue
% Create a new figure window
figure;
% ind = find(Y >= 1000);
% plot(Y(ind), -X(ind), '.g');
%%%% MROE734 EDITS plot based on colours
% Plot individual based on labels
% Swap X and Y coordinates and negate Y to rotate 90 degrees clockwise
% Noticed that the image seems to always be rotated 90 degrees
% counterclockwise, so this should correct that
plot(Y(labels == 10), -X(labels == 10), '.g'); hold on;
plot(Y(labels == 1), -X(labels == 1), '.r');
plot(Y(labels == 2), -X(labels == 2), '.b');
plot(Y(labels == 3), -X(labels == 3), '.k');

%%%ColeD EDITS added hold off in an attempt to fix the plotting issue
hold off;
