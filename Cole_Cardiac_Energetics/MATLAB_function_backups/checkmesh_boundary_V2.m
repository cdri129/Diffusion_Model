
imageName = 'myo_centered_mito';
base_dir = join([getenv("HOME"), "Desktop/Cole_Modelling/Triangle_Modelling/Mathias_Made"], '/');
img_dir = join([base_dir, "img"], '/');
mesh_dir = join([base_dir, "mesh"], '/');
file_dir = join([mesh_dir, imageName], '/');

disp("Processing image " + imageName);

file_name = strcat(imageName,'.1.csv');
%file_path = join([file_dir, file_name], '/');

%This following lines is for Riti's images or other files just put in the
%mesh folder and not in an organized subfolder
file_path = join([mesh_dir, file_name], '/');

fnode = readtable(file_path);
Y = fnode.Var2;
X = fnode.Var3;
labels = fnode.Var4;
conditions_1 = fnode.Var4 == 10 & fnode.Var5 == 1; % Magenta
conditions_2 = fnode.Var4 == 1 & fnode.Var5 == 2;  % Blue
conditions_3 = fnode.Var4 == 10 & fnode.Var5 == 3; % Yellow
conditions_4 = (fnode.Var4 == 1) & fnode.Var5 == 0; % Red

figure;


% Plot nodes based on conditions
hold on;
plot(Y(conditions_1), -X(conditions_1), '.m'); % Magenta
plot(Y(conditions_2), -X(conditions_2), '.b'); % Blue
plot(Y(conditions_3), -X(conditions_3), '.y'); % Yellow
plot(Y(conditions_4), -X(conditions_4), '.r'); %Red
% % Plot other nodes
% plot(Y(labels == 10 & ~conditions_1), -X(labels == 10 & ~conditions_1), '.g');
% plot(Y(labels == 1 & ~conditions_2), -X(labels == 1 & ~conditions_2), '.r');
% plot(Y(labels == 2 & ~conditions_2), -X(labels == 2 & ~conditions_2), '.b');
% plot(Y(labels == 3 & ~conditions_3), -X(labels == 3 & ~conditions_3), '.k');

hold off;
