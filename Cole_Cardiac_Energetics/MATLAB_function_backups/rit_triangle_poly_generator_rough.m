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

parent = '.\';

total_edges = 0;
total_sl_edges = 0;
total_sl_nodes = 0;
total_mito_edges = 0;
total_ryr_nodes = 0;
total_sl_mito_nodes = 0;
total_cell_nodes = 0;
%% 1,8 2,1
for j = 1
for i = 1
   % Cell_name = strcat(parent,'STZ',int2str(j),'\classified9',sprintf('%.1f',STZ(j).Names(i)),'.tif');
    Resolution = 13.7;
    
    A = imread('binary_2D2_2.tif');
    Cell = zeros(size(A)+1200);    
    Cell(1201:end, 1201:end) = A;
%    Cell = imread('binary_1C2_3.tif');
    
    % Identify regions based on image color
    mito(i) = sum(Cell(:)==85); 
    myo(i) = sum(Cell(:)==155);  
    %sr(i) = sum(Cell(:)==155);  
    glyco(i) = sum(Cell(:)==163);  
    white(i) = sum(Cell(:)==190);  
    
    size_cell(i)= (myo(i)+mito(i)+white(i)+glyco(i));
    Cell_size(i)= size_cell(i)*Resolution*Resolution/1100/1100;
    
    
    % Create duplicates
    mito_bin = Cell;
    myo_bin = Cell;
    glyco_bin = Cell;
    %sr_bin = Cell;
    tt_bin = Cell;

    %Assign duplicates to invidual regions
    mito_bin(mito_bin ~= 85) = 0;
    myo_bin(myo_bin ~= 155) = 0;
    glyco_bin(glyco_bin ~= 163) = 0;
    %sr_bin(sr_bin ~= 155) = 0;
    tt_bin(tt_bin ~= 190) = 0;

    mito_bin = imbinarize(mito_bin);
    myo_bin = imbinarize(myo_bin);
    glyco_bin = imbinarize(glyco_bin);
    %sr_bin = imbinarize(sr_bin);
    tt_bin = imbinarize(tt_bin);
    cell_bin = mito_bin + myo_bin + glyco_bin + tt_bin;

% Visualize regions in color
   %% 
    RGB = zeros(size(Cell));
    RGB = RGB + cat(3, (tt_bin+myo_bin+glyco_bin)*1, (tt_bin+glyco_bin)*1, (tt_bin+mito_bin)*1);
    RGB = imresize(RGB,Resolution/11);
    imshow(RGB);
    
    mito(i) = sum(mito_bin(:)==1);  
    
    
% Calculate volume fraction of regions  
    P(1,i) = (myo(i))/size_cell(i);
    P(2,i) = (mito(i))/size_cell(i);
    P(3,i) = (glyco(i))/size_cell(i);
    %P(4,i) = (sr(i))/size_cell(i);
    P(5,i) = (white(i))/size_cell(i);  
end
end


%% create Cell Boundary PLC

% Rescale the image, else it takes edges to run
cell = imresize(cell_bin,Resolution/10);
cell = imbinarize(cell);

%add a border of blank area 
cell (1:size(cell,1),1:2) = 0;
cell (1:size(cell,1),size(cell,2)-1:size(cell,2)) = 0;
cell (1:2,1:size(cell,2)) = 0;
cell (size(cell,1)-1:size(cell,1),1:size(cell,2)) = 0;

% Without this block, P Kovesis's algorithm picks up small areas as complete cell
cellSL_binary = edge(cell,'log','thinning');
cellSL_binary = imfill(cellSL_binary,'holes');
 se =  strel('disk',4);
 cellSL_binary = imerode(cellSL_binary,se);
 cellSL_binary = bwareaopen(cellSL_binary,5000);
cellSL_binary = edge(cellSL_binary,'log','thinning');
   
   
cell_boundary_indices = find(cellSL_binary);
[celly,cellx] = ind2sub(size(cell),cell_boundary_indices);
centre = [mean(celly),mean(cellx)];

%using P. Kovesi's code to get an ordered list of edges via an edge-linking
%algorithm

[cell_edgelist,cell_edgeim] = edgelink(cellSL_binary);
cell_edgepoints = cell2mat(cell_edgelist(1,1));

% Calculate the node points and edges of the cell membrane
cell_edgepoints = (cell_edgepoints-ones(size(cell_edgepoints)))*spatres*1; %why ??
%cell_edgepoints(:,3) = zeros(1,size(cell_edgepoints,1));
cell_edgepoints = cell_edgepoints(1:end-1,:);
%create list of edges with node indices
total_edges = size(cell_edgepoints,1);
for k = 1:(total_edges-1)
    edges(k,1) = k;
    edges(k,2) = k;
    edges(k,3) = k+1;
end
edges(k+1,1) = k+1;
edges(k+1,2) = k+1;
edges(k+1,3) = 1;
total_sl_edges = total_edges;
all_nodes = [cell_edgepoints];
all_lables = 10*ones(size(cell_edgepoints));

total_sl_nodes = size(all_nodes,1);
holes = [];
IMS = [];

progressbar('Creating Mesh','T-tubules','Mitochondria','Myofibrils','Glycogen','./Triangle','Assigning Attributes') 
progressbar([1/11-0.0001],[], [], [], [], [], []);       

%% create T-tutbule PLC
TT = imresize(tt_bin,Resolution/1);

 se =  strel('disk',4);
 % If you see ERROR: Blank / No mesh : hole has been placed outside a TT
 % Try increasing the se value to 10
 TT1 = imerode(TT,se);
 TT1 = bwareaopen(TT1,5000); %Para 2

 TT = imdilate(TT1,se);

%add a border
TT (1:size(TT,1),1:2) = 0; 
TT (1:size(TT,1),size(TT,2)-1:size(TT,2)) = 0;
TT (1:2,1:size(TT,2)) = 0;
TT (size(TT,1)-1:size(TT,1),1:size(TT,2)) = 0;

% TT = edge(TT,'log','thinning');
% TT = imfill(TT,'holes');

[TT_labelled,TT_nums] = bwlabel(TT,8);
% vislabels(TT_labelled);

    for label = 1:TT_nums 
        [r,c] = find(TT_labelled==label);   %find row and columns beloonging to label
        tt = bwselect(TT_labelled,c,r,8);
        tt = imfill(tt,'holes');
        tt = imresize(tt,1/10);
        TT_boundary = edge(tt,'log','thinning');
        allTT_edgepoints = [];
        [TT_edgelist,TT_edgeim] = edgelink(TT_boundary);
        TT_edgepoints = cell2mat(TT_edgelist(1,1));
        TT_edgepoints = TT_edgepoints(1:end-1,:);
        TT_edgepoints_spatcoor = (TT_edgepoints-ones(size(TT_edgepoints)))*spatres;
        %mito_edgepoints(:,3) = label*ones(1,size(mito_edgepoints,1));
        allTT_edgepoints = [allTT_edgepoints;TT_edgepoints_spatcoor];
        
        %add mitochondrial edges to edges list
             for k = (total_edges+1):(total_edges+size(TT_edgepoints,1)-1)
            edges(k,1) = k;
            edges(k,2) = k;
            edges(k,3) = k+1;
            end
        edges(k+1,1) = k+1;
        edges(k+1,2) = k+1;
        edges(k+1,3) = total_edges+1;
        total_edges = total_edges+size(TT_edgepoints,1);
            all_nodes = [all_nodes;allTT_edgepoints];
            all_lables = [all_lables;label*ones(size(allTT_edgepoints))];
              
            se =  strel('disk',1);  %% Para 3
            BW = imerode(tt,se);
            [r,c] = find(BW);    
%             pixel = [[],[]];
            pixel = [r(round(size(r,1)/2)),c(round(size(c,1)/2))];
            points = (pixel-ones(size(pixel)))*spatres;
            holes = [holes;points];
    
%      label
     progressbar([],label/TT_nums-0.0001, [], [], [], [], []);       
    end
    total_TT_edges = total_edges-total_sl_edges;
progressbar([2/11-0.0001],[], [], [], [], [], []);   

 %% create Mitochondria PLC
MT = imresize(mito_bin,Resolution/1);
MT = bwareaopen(MT,5000); %Para 2

 se =  strel('disk',50);% Para 1
 MT = imerode(MT,se);

%add a border
MT (1:size(MT),1:2) = 0;
MT (1:size(MT,1),size(MT,2)-1:size(MT,2)) = 0;
MT (1:2,1:size(MT,2)) = 0;
MT (size(MT,1)-1:size(MT,1),1:size(MT,2)) = 0;

% TT = edge(TT,'log','thinning');
% TT = imfill(TT,'holes');

[MT_labelled,MT_nums] = bwlabel(MT,8);
% vislabels(MT_labelled);

    for label = 1:MT_nums 
        [r,c] = find(MT_labelled==label);   %find row and columns beloonging to label
        MT = bwselect(MT_labelled,c,r,8);
        MT = imdilate(MT,se);
        MT = imfill(MT,'holes');
        MT = imresize(MT,1/10);
        MT_boundary = edge(MT,'log','thinning');
        allMT_edgepoints = [];
        [MT_edgelist,MT_edgeim] = edgelink(MT_boundary);
        MT_edgepoints = cell2mat(MT_edgelist(1,1));
        MT_edgepoints = MT_edgepoints(1:end-1,:);
        MT_edgepoints_spatcoor = (MT_edgepoints-ones(size(MT_edgepoints)))*spatres;
        %mito_edgepoints(:,3) = label*ones(1,size(mito_edgepoints,1));
        allMT_edgepoints = [allMT_edgepoints;MT_edgepoints_spatcoor];
        
        %add mitochondrial edges to edges list
             for k = (total_edges+1):(total_edges+size(MT_edgepoints,1)-1)
            edges(k,1) = k;
            edges(k,2) = k;
            edges(k,3) = k+1;
            end
        edges(k+1,1) = k+1;
        edges(k+1,2) = k+1;
        edges(k+1,3) = total_edges+1;
        total_edges = total_edges+size(MT_edgepoints,1);
            all_nodes = [all_nodes;allMT_edgepoints];
            all_lables = [all_lables;label*ones(size(allMT_edgepoints))];
            
%                  label
                 progressbar([],[], label/MT_nums-0.0001, [], [], [], []);       
    end
    total_MT_edges = total_edges-total_sl_edges-total_TT_edges;       
    progressbar([3/11-0.0001],[], [], [], [], [], []);  
    
     %% create Myofibrils PLC
MY = imresize(myo_bin,Resolution/1);
MY = bwareaopen(MY,5000); %Para 2

 se =  strel('disk',40);% Para 1
 MY = imerode(MY,se);


%add a border
MY (1:size(MY),1:2) = 0;
MY (1:size(MY,1),size(MY,2)-1:size(MY,2)) = 0;
MY (1:2,1:size(MY,2)) = 0;
MY (size(MY,1)-1:size(MY,1),1:size(MY,2)) = 0;

% TT = edge(TT,'log','thinning');
% TT = imfill(TT,'holes');

[MY_labelled,MY_nums] = bwlabel(MY,8);
% vislabels(MY_labelled);

    for label = 1:MY_nums 
        [r,c] = find(MY_labelled==label);   %find row and columns belonging to label
        MY = bwselect(MY_labelled,c,r,8);
        MY = imdilate(MY,se);        
        MY = imfill(MY,'holes');
        MY = imresize(MY,1/10);
        MY_boundary = edge(MY,'log','thinning');
        allMY_edgepoints = [];
        [MY_edgelist,MY_edgeim] = edgelink(MY_boundary);
        MY_edgepoints = cell2mat(MY_edgelist(1,1));
        MY_edgepoints = MY_edgepoints(1:end-1,:);
        MY_edgepoints_spatcoor = (MY_edgepoints-ones(size(MY_edgepoints)))*spatres;
        %mito_edgepoints(:,3) = label*ones(1,size(mito_edgepoints,1));
        allMY_edgepoints = [allMY_edgepoints;MY_edgepoints_spatcoor];
        
        %add mitochondrial edges to edges list
             for k = (total_edges+1):(total_edges+size(MY_edgepoints,1)-1)
            edges(k,1) = k;
            edges(k,2) = k;
            edges(k,3) = k+1;
            end
        edges(k+1,1) = k+1;
        edges(k+1,2) = k+1;
        edges(k+1,3) = total_edges+1;
        total_edges = total_edges+size(MY_edgepoints,1);
            all_nodes = [all_nodes;allMY_edgepoints];
            all_lables = [all_lables;label*ones(size(allMY_edgepoints))];
            
%                              label
                             progressbar([],[],[],label/MY_nums-0.0001, [], [], []);       
    end
    total_MY_edges = total_edges-total_sl_edges-total_TT_edges-total_MT_edges;                   
    progressbar([4/11-0.0001],[], [], [], [], [], []);       

%% create Glycogen PLC
GL = imresize(glyco_bin,Resolution/1);
GL = bwareaopen(GL,5000); %Para 2

 se =  strel('disk',10);% Para 1
 GL = imerode(GL,se);


%add a border
GL (1:size(GL),1:2) = 0;
GL (1:size(GL,1),size(GL,2)-1:size(GL,2)) = 0;
GL (1:2,1:size(GL,2)) = 0;
GL (size(GL,1)-1:size(GL,1),1:size(GL,2)) = 0;

% TT = edge(TT,'log','thinning');
% TT = imfill(TT,'holes');

[GL_labelled,GL_nums] = bwlabel(GL,8);
% vislabels(GL_labelled);

    for label = 1:GL_nums 
        [r,c] = find(GL_labelled==label);   %find row and columns beloonging to label
        GL = bwselect(GL_labelled,c,r,8);
        GL = imdilate(GL,se);        
        GL = imfill(GL,'holes');
        GL = imresize(GL,1/10);
        GL_boundary = edge(GL,'log','thinning');
        allGL_edgepoints = [];
        [GL_edgelist,GL_edgeim] = edgelink(GL_boundary);
        GL_edgepoints = cell2mat(GL_edgelist(1,1));
        GL_edgepoints = GL_edgepoints(1:end-1,:);
        GL_edgepoints_spatcoor = (GL_edgepoints-ones(size(GL_edgepoints)))*spatres;
        %mito_edgepoints(:,3) = label*ones(1,size(mito_edgepoints,1));
        allGL_edgepoints = [allGL_edgepoints;GL_edgepoints_spatcoor];
        
        %add mitochondrial edges to edges list
             for k = (total_edges+1):(total_edges+size(GL_edgepoints,1)-1)
            edges(k,1) = k;
            edges(k,2) = k;
            edges(k,3) = k+1;
            end
        edges(k+1,1) = k+1;
        edges(k+1,2) = k+1;
        edges(k+1,3) = total_edges+1;
        total_edges = total_edges+size(GL_edgepoints,1);
            all_nodes = [all_nodes;allGL_edgepoints];
            all_lables = [all_lables;label*ones(size(allGL_edgepoints))];
            
%                              label
                             progressbar([],[],[], [],label/GL_nums-0.0001, [], []);       
    end
    total_GL_edges = total_edges-total_sl_edges-total_TT_edges-total_MT_edges-total_MY_edges;                   
%   
    progressbar([5/11-0.0001],[], [], [], [], [], []);       


%% write out PLC nodes
%TRIANGLE = 'triangle_folder';
%outfolder = strcat(TRIANGLE,'\output\MESH\');
%create out folder
%mkdir(outfolder);
%print out .poly file for testing with triangle
outfile = strcat('2D2_2','.poly');
[fid,msg] = fopen(outfile,'w','native');
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
    progressbar([6/11-0.0001],[], [], [], [], [], []);       

%% write out PLC edges
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
    progressbar([7/11-0.0001],[], [], [], [], [], []);       

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
    progressbar([8/11-0.0001],[], [], [], [], [], []);       

%% Use Triangle
'Triangulation Starts'
    progressbar([],[], [], [], [], 1/3--0.0001, []);       

terminal = 'gnome-terminal';
%TRIANGLE = 'triangle_folder/';
OUTFILE = strcat('2D2_2.poly');
EXEFILE = strcat('triangle.exe');
    progressbar([],[], [], [], [], 2/3--0.0001, []);       

%system([terminal,' -l -c "',EXEFILE,' ',OUTFILE,' -j -a10 -q20"']);

system([terminal,' ls "',EXEFILE,' ',OUTFILE,' -j -a10 -q20"']);

'Triangulation Finish'
    progressbar([],[], [], [], [], 3/3--0.0001, []);       
    progressbar([9/11-0.0001],[], [], [], [], [], []);       

%% Post processing - we reassign each traingle node to correct attrbutes based on regions
NodePath = strcat('2D2_2','.1.node');
%NodePath = ('2COMP.1.node')
NodeID = fopen(NodePath,'r');
NODES = str2num(fgetl(NodeID));
NUMBER_OF_NODES = NODES(1);
ATRBT = zeros(1, NUMBER_OF_NODES);
tot_myo_nodes = 0;

cell_nodes = zeros(size(cell_bin));
mito_nodes = zeros(size(cell_bin));
myo_nodes = zeros(size(cell_bin));
glyco_nodes = zeros(size(cell_bin));
tt_nodes = zeros(size(cell_bin));
%sr_nodes = zeros(size(cell_bin));

 for k = 1:NUMBER_OF_NODES
   NODES = str2num(fgetl(NodeID));
   XX(k) = NODES(2);
   YY(k) = NODES(3);
   NODEX(k) = round((NODES(2)+1)*10/Resolution);
   NODEY(k) = round((NODES(3)+1)*10/Resolution);
   ATRBT(k) = NODES (4);
   BC(k) = NODES (5);
   cell_nodes(NODEY(k),NODEX(k))=1;
% % 
   if(mito_bin(NODEY(k),NODEX(k))==1 && BC(k)==0)
         ATRBT(k) =  1.0;
         mito_nodes(NODEY(k)-1:NODEY(k)+1,NODEX(k)-1:NODEX(k)+1)=1;
   elseif (glyco_bin(NODEY(k),NODEX(k))==1 && BC(k)==0) 
         ATRBT(k) =  2.0;
         glyco_nodes(NODEY(k)-1:NODEY(k)+1,NODEX(k)-1:NODEX(k)+1)=1; 
  elseif (tt_bin(NODEY(k),NODEX(k))==1 && BC(k)==0) 
         ATRBT(k) =  3.0;
         tt_nodes(NODEY(k)-1:NODEY(k)+1,NODEX(k)-1:NODEX(k)+1)=1;   
  elseif (myo_bin(NODEY(k),NODEX(k))==1 && BC(k)==0) 
         ATRBT(k) =  10.0;
         myo_nodes(NODEY(k)-1:NODEY(k)+1,NODEX(k)-1:NODEX(k)+1)=1;    
%   elseif (sr_bin(NODEY(k),NODEX(k))==1 && BC(k)==0) 
%          ATRBT(k) =  10.0;
%          myo_nodes(NODEY(k)-1:NODEY(k)+1,NODEX(k)-1:NODEX(k)+1)=1;            
   end
   NODEX(j) = NODES(2);
   NODEY(j) = NODES(3);
       progressbar([],[], [], [], [], [], j/NUMBER_OF_NODES*0.5--0.0001);       
 end
   fclose (NodeID);   
       progressbar([10/11-0.0001],[], [], [], [], [], []);       

NodeID = fopen(NodePath,'w');
fprintf(NodeID,'%d\t 2\t 1\t 1\n',NUMBER_OF_NODES);
 for k = 1:NUMBER_OF_NODES
fprintf(NodeID,'%d\t%d\t%d\t%d\t%d\t \n ',k,XX(k),YY(k),ATRBT(k),BC(k));
    progressbar([],[], [], [], [], [], 0.5+j/NUMBER_OF_NODES*0.5--0.0001);       
 end
fclose('all');
progressbar([11/11-0.0001],[], [], [], [], [], []);       
progressbar(1) ; 

    figure;
    RGB = zeros(size(cell_nodes));
    RGB = RGB + cat(3,  (mito_nodes+tt_nodes)*1, (myo_nodes+tt_nodes)*1, (tt_nodes+glyco_nodes)*1);
    %RGB = RGB + cat(3, (myo_nodes+tt_nodes)*1, (mito_nodes+tt_nodes)*1, (tt_nodes+glyco_nodes)*1);
    imshow(RGB);
    size(RGB)
