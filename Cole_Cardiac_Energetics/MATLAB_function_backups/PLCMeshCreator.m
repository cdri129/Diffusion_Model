% % %Author: Shouryadipta Ghosh 
% %A program to extract the edge of the cell boundary to create a PLC file.
% %%% IMOD -> Draw Model -> Model View -> Goto Objects -> Toggle between objects -> Take snapshots ->
% %%% -> I                 mport snapshots to Fiji -> Make Binary -> Fill Holes -> images for MATLAB
clear all
close all
clc
%% initialise these variables
% spatres = 0.011375; %the spatial resolution of the image stack being read in.

%spatres = 12.775; %
spatres =  1.0; %

%%%% MROE734 EDITS harmonised indentation
%%%% MROE734 EDITS created directory variables
% This should be the base of everything, HOME is /home/your_upi 
% change Documents/cole to whatever you need to be
imageName = 'myo_centered_mito';
base_dir = join([getenv("HOME"), "Desktop/Cole_Modelling/Triangle_Modelling/Mathias_Made"], '/');
img_dir = join([base_dir, "img"], '/');
mesh_dir = join([base_dir, "mesh"], '/');

%%%% MROE734 EDITS commented next three lines
% total_edges = 0;
% total_sl_edges = 0;
% total_sl_nodes = 0;
total_mito_edges = 0;
total_ryr_nodes = 0;
total_sl_mito_nodes = 0;
total_cell_nodes = 0;
%% 1,8 2,1

%%%% MROE734 EDITS moved image naming up here
%%%% MROE734 EDITS changed the full image path method
%%%% MROE734 EDITS add semicolon to imageName and a disp statement
disp("Processing image " + imageName); 
imageFileName = strcat(imageName, '.tif');
fullImagePath = join([img_dir, imageFileName], '/');

%%%% MROE734 EDITS remove loops that did nothing
Resolution = 13.7;
upscale_factor = 1;

% Cell_name = strcat(parent,'STZ',int2str(j),'\classified9',sprintf('%.1f',STZ(j).Names(i)),'.tif');

%%%% MROE734 EDITS renamed A to img and Cell to resized_img for clarity
img = imread(fullImagePath);
img = imresize(img,upscale_factor);
%%%% MROE734 EDITS removed resizing but left conversion to double and added
%%%% border
img_in_double = double(padarray(img, [5, 5], 'both'));
%img_in_double = double(padarray(img, [100, 100], 'both'));
% resized_img = zeros(size(img)+1200);
% resized_img(1201:end, 1201:end) = img;

%%%% MROE734 EDITS wrote the identifyRegions function to clean up this
%%%% section
% Identify regions based on image color
[mito_size, mito_bin] = identifyRegions(img_in_double, 85);
[myo_size, myo_bin] = identifyRegions(img_in_double, 155);
[glyco_size, glyco_bin] = identifyRegions(img_in_double, 190);
[tt_size, tt_bin] = identifyRegions(img_in_double, 162);

cell_size = mito_size + myo_size + glyco_size + tt_size;
cell_bin = mito_bin + myo_bin + glyco_bin + tt_bin;

% Visualize regions in color
%%
RGB = zeros(size(img_in_double));
RGB = RGB + cat(3, (tt_bin+myo_bin+glyco_bin)*1, (tt_bin+glyco_bin)*1, (tt_bin+mito_bin)*1);
%RGB = imresize(RGB,Resolution/11);
imshow(RGB);

%%%% MROE734 EDITS commented out because unused
% % Calculate volume fraction of regions
% P(1, 1) = (myo_size)/cell_size;
% P(2, 1) = (mito_size)/cell_size;
% P(3, 1) = (glyco_size)/cell_size;
% P(4, 1) = (sr)/size_cell;
% P(5, 1) = (tt_size)/cell_size;


%% create Cell Boundary PLC

% Rescale the image, else it takes edges to run
cell_img = cell_bin;
% cell_img = imresize(cell_bin,Resolution/10);
cell_img = imbinarize(cell_img);

%%%ColeD EDITS noticed there seems to be a white border to the mesh, but
%%%that is located inside of the image instead of at the border. I guess
%%%the issue comes from this section, but removing this section causes the
%%% following code to not work
%add a border of blank area
cell_img (1:size(cell_img,1),1:2) = 0;
cell_img (1:size(cell_img,1),size(cell_img,2)-1:size(cell_img,2)) = 0;
cell_img (1:2,1:size(cell_img,2)) = 0;
cell_img (size(cell_img,1)-1:size(cell_img,1),1:size(cell_img,2)) = 0;

% Without this block, P Kovesis's algorithm picks up small areas as complete cell
cellSL_binary = edge(cell_img,'log','thinning');
cellSL_binary = imfill(cellSL_binary,'holes');
se =  strel('disk',4);
cellSL_binary = imerode(cellSL_binary,se);
cellSL_binary = bwareaopen(cellSL_binary,5000);
cellSL_binary = edge(cellSL_binary,'log','thinning');


cell_boundary_indices = find(cellSL_binary);
[celly,cellx] = ind2sub(size(cell_img),cell_boundary_indices);
centre = [mean(celly),mean(cellx)];

%using P. Kovesi's code to get an ordered list of edges via an edge-linking
%algorithm

[cell_edgelist,cell_edgeim] = edgelink(cellSL_binary);
cell_edgepoints = cell2mat(cell_edgelist(1,1));

% Calculate the node points and edges of the cell membrane
cell_edgepoints = (cell_edgepoints-ones(size(cell_edgepoints)))*spatres*1; %why ??
%cell_edgepoints(:,3) = zeros(1,size(cell_edgepoints,1));
cell_edgepoints = cell_edgepoints(1:end-1,:);

% Plot edge of 'cell'
% figure(2)
% plot(cell_edgepoints(:,1),cell_edgepoints(:,2))

%create list of edges with node indices
total_edges = size(cell_edgepoints,1);
%%%% MROE734 EDITS removed useless for loop and replaced with vectors
tmp_edges = (1:total_edges)';
edges = [tmp_edges, tmp_edges, tmp_edges+1];
edges(end,3) = 1;
total_sl_edges = total_edges;
all_nodes = cell_edgepoints;
all_labels = 10*ones(size(cell_edgepoints)); %%%% MROE734 EDITS fix typo in label


total_sl_nodes = size(all_nodes,1);
holes = [];
IMS = [];


%% create T-tutbule PLC
%%%% MROE734 EDITS removed progressbar and replaced with print
disp("Creating T-tubule PLC")

%%%% MROE734 EDITS wrote a function to create PLC for different regions
[total_edges, ...
    TT_nodes, ...
    TT_labels, ...
    TT_edges, ...
    TT_holes] = createTTPLC(tt_bin, total_edges, spatres);

all_nodes = [all_nodes; TT_nodes];
all_labels = [all_labels; TT_labels];
edges = [edges; TT_edges];
holes = [holes; TT_holes];

total_TT_edges = total_edges-total_sl_edges;



%% create Mitochondria PLC
disp("Creating mitochondria PLC")
[total_edges, ...
    MT_nodes, ...
    MT_labels, ...
    MT_edges, ...
    MT_holes] = createMTPLC(mito_bin, total_edges, spatres);

all_nodes = [all_nodes; MT_nodes];
all_labels = [all_labels; MT_labels];
edges = [edges; MT_edges];
holes = [holes; MT_holes];

total_MT_edges = total_edges-total_sl_edges-total_TT_edges;


%% create Myofibrils PLC
disp("Creating myofibrils PLC")

[total_edges, ...
    MY_nodes, ...
    MY_labels, ...
    MY_edges, ...
    MY_holes] = createMYPLC(myo_bin, total_edges, spatres);

all_nodes = [all_nodes; MY_nodes];
all_labels = [all_labels; MY_labels];
edges = [edges; MY_edges];
holes = [holes; MY_holes];

total_MY_edges = total_edges-total_sl_edges-total_TT_edges-total_MT_edges;


%% create Glycogen PLC
disp("Creating glycogen PLC")

[total_edges, ...
    GL_nodes, ...
    GL_labels, ...
    GL_edges, ...
    GL_holes] = createGLPLC(glyco_bin, total_edges, spatres);

all_nodes = [all_nodes; GL_nodes];
all_labels = [all_labels; GL_labels];
edges = [edges; GL_edges];
holes = [holes; GL_holes];

total_GL_edges = total_edges-total_sl_edges-total_TT_edges-total_MT_edges-total_MY_edges;

%% write out PLC nodes
disp("Writing out PLC nodes")
%TRIANGLE = 'triangle_folder';
%outfolder = strcat(TRIANGLE,'\output\MESH\');
%create out folder
%mkdir(outfolder);
%print out .poly file for testing with triangle
%%%% MROE734 EDITS changed save location to mesh directory
mesh_name = strcat(imageName,'.poly');
mesh_path = join([mesh_dir, mesh_name], '/');
[fid,msg] = fopen(mesh_path,'w','native');
numNodes = size(all_nodes,1);


%first write out nodes
fprintf (fid,strcat(num2str(total_edges),' 2 1 1\n'));
for k = 1:total_sl_nodes
    fprintf (fid,strcat(num2str(k),'\t',num2str(all_nodes(k,2)),'\t',num2str(all_nodes(k,1)),'\t',num2str(10),'\t',num2str(1),'\n'));
end
%
for k = (total_sl_nodes+1):(total_sl_nodes+total_TT_edges)
    fprintf (fid,strcat(num2str(k),'\t',num2str(all_nodes(k,2)),'\t',num2str(all_nodes(k,1)),'\t',num2str(3),'\t',num2str(1),'\n'));
end

for k = (total_sl_nodes+total_TT_edges+1):(total_sl_nodes+total_TT_edges+total_MT_edges)
    fprintf (fid,strcat(num2str(k),'\t',num2str(all_nodes(k,2)),'\t',num2str(all_nodes(k,1)),'\t',num2str(1),'\t',num2str(2),'\n'));
end

for k = (total_sl_nodes+total_TT_edges+total_MT_edges+1):(total_sl_nodes+total_TT_edges+total_MT_edges+total_MY_edges)
    fprintf (fid,strcat(num2str(k),'\t',num2str(all_nodes(k,2)),'\t',num2str(all_nodes(k,1)),'\t',num2str(10),'\t',num2str(3),'\n'));
end

for k = (total_sl_nodes+total_TT_edges+total_MT_edges+total_MY_edges+1):(total_sl_nodes+total_TT_edges+total_MT_edges+total_MY_edges+total_GL_edges)
    fprintf (fid,strcat(num2str(k),'\t',num2str(all_nodes(k,2)),'\t',num2str(all_nodes(k,1)),'\t',num2str(2),'\t',num2str(4),'\n'));
end

%% write out PLC edges
disp("Writing out PLC edges")
%next write out edge segments
%total edges = edges for sarcolemma + edges for each mito


fprintf (fid,strcat(num2str(total_edges),' 0\n'));

for k = 1:(total_sl_edges)
    fprintf (fid,strcat(num2str(edges(k,1)),'\t',num2str(edges(k,2)),'\t',num2str(edges(k,3)),'\t',num2str(1),'\n'));
end

for k = (total_sl_edges+1):(total_sl_edges+total_TT_edges)
    fprintf (fid,strcat(num2str(edges(k,1)),'\t',num2str(edges(k,2)),'\t',num2str(edges(k,3)),'\t',num2str(2),'\n'));
end

for k = (total_sl_edges+total_TT_edges+1):(total_sl_edges+total_TT_edges+total_MT_edges)
    fprintf (fid,strcat(num2str(edges(k,1)),'\t',num2str(edges(k,2)),'\t',num2str(edges(k,3)),'\t',num2str(3),'\n'));
end

for k = (total_sl_edges+total_TT_edges+total_MT_edges+1):(total_sl_edges+total_TT_edges+total_MT_edges+total_MY_edges)
    fprintf (fid,strcat(num2str(edges(k,1)),'\t',num2str(edges(k,2)),'\t',num2str(edges(k,3)),'\t',num2str(4),'\n'));
end

for k = (total_sl_edges+total_TT_edges+total_MT_edges+total_MY_edges+1):(total_sl_edges+total_TT_edges+total_MT_edges+total_MY_edges+total_GL_edges)
    fprintf (fid,strcat(num2str(edges(k,1)),'\t',num2str(edges(k,2)),'\t',num2str(edges(k,3)),'\t',num2str(4),'\n'));
end

%% specify holes

% numholes = 0;
% fprintf (fid,strcat(num2str(numholes),'\n'));
% fprintf (fid,strcat(num2str(numholes),'\n'));

%now specify holes
numholes = size(holes,1);


fprintf (fid,strcat(num2str(numholes),'\n'));
for k = 1:numholes
    fprintf (fid,strcat(num2str(k),'\t',num2str(holes(k,2)),'\t',num2str(holes(k,1)),'\t',num2str(2),'\n'));
end
fclose(fid);



%%%% MROE734 EDITS removed MATLAB usage of triangle

% %% Use Triangle
% 'Triangulation Starts'
% progressbar([],[], [], [], [], 1/3--0.0001, []);
% 
% %terminal = 'gnome-terminal';
% TRIANGLE = '/home/cdri129/Desktop/Cole_Modelling/Triangle_Modelling/Riti_Made/Matlab/Triangle_stuff';
% 
% outputFolder = fullfile(TRIANGLE, 'Triangle_outputs');
% 
% % Create the subfolder if it doesn't exist
% if ~isfolder(outputFolder)
%     mkdir(outputFolder);
% else
%     disp([outputFolder, ' folder already exists']);
% end
% 
% % Change the working directory to the subfolder
% cd(outputFolder);
% 
% % Specify the input and output file paths
% inputFile = strcat(imageName,'.poly');
% %TRI_OUTFILE = strcat('D001.poly');
% %TRI_EXE = strcat('triangle.exe');
% 
% % Specify file paths within the subfolder
% TRI_OUTFILE = fullfile(outputFolder, inputFile);
% TRI_EXE = fullfile(outputFolder, 'triangle.exe');
% 
% %OUTFILE = strcat('2D2_2.poly');
% %EXEFILE = strcat('triangle');
% progressbar([],[], [], [], [], 2/3--0.0001, []);
% 
% outputFile = fullfile(outputFolder, imageName);
% 
% % Specify the full path to the triangle executable
% triangleExecutable = fullfile(TRIANGLE, 'triangle');  % Assuming triangle.exe is in the same folder
% % Run Triangle directly in MATLAB
% %cmd = ['"', TRI_EXE, '" "', TRI_OUTFILE, '" -j -a10 -q20 -C'];
% cmd = ['"', triangleExecutable, '" "', inputFile, '" -j -a10 -q20 -C'];
% status = system(cmd);
% 
% if status ~= 0
%     error('Triangle command failed with status code %d', status);
% end
% %system([terminal,' -l -c "',EXEFILE,' ',OUTFILE,' -j -a10 -q20"']);
% 
% %system([terminal,' ls "',EXEFILE,' ',OUTFILE,' -j -a10 -q20"']);
% 
% 'Triangulation Finish'
% progressbar([],[], [], [], [], 3/3--0.0001, []);
% progressbar(9/11-0.0001,[], [], [], [], [], []);

%%%% MROE734 EDITS removed postprocessing and created a new script
