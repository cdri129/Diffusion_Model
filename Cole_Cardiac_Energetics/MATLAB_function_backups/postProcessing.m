%% Post processing - we reassign each traingle node to correct attrbutes based on regions
%%%% Lines to be ran in terminal:
% cd /home/cdri129/Desktop/Cole_Modelling/Triangle_Modelling/Mathias_Made/mesh
% triangle imageName.poly -j -a10 -q20 -O


imageName = 'myo_centered_mito';
% Copied from other script (not needed if running without clear all)
base_dir = join([getenv("HOME"), "Desktop/Cole_Modelling/Triangle_Modelling/Mathias_Made"], '/');
img_dir = join([base_dir, "img"], '/');
mesh_dir = join([base_dir, "mesh"], '/');
%%%% MROE734 EDITS add semicolon to imageName and a disp statement
disp("Processing image " + imageName); 

node_name = strcat(imageName,'.1.node');
node_path = join([mesh_dir, node_name], '/');
csv_name = strcat(imageName, '.1.csv');
csv_path = join([mesh_dir, csv_name], '/');
NodeID = fopen(node_path,'r');
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

%%%% MROE734 EDITS
% Pre-allocate space for node information
XX = zeros(1, NUMBER_OF_NODES);
YY = XX;
NODEX = XX;
NODEY = XX;
ATRBT = XX;
BC = XX;

for k = 1:NUMBER_OF_NODES
    NODES = str2num(fgetl(NodeID));
    XX(k) = NODES(2);
    YY(k) = NODES(3);
    NODEX(k) = round(NODES(2)+1);%round((NODES(2)+1)*10/Resolution);
    NODEY(k) = round(NODES(3)+1);%round((NODES(3)+1)*10/Resolution);
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
    NODEX(1) = NODES(2);
    NODEY(1) = NODES(3);
end
fclose (NodeID);

NodeID = fopen(node_path,'w');
fprintf(NodeID,'%d\t 2\t 1\t 1\n ',NUMBER_OF_NODES); %%%% Added a space at the end
for k = 1:NUMBER_OF_NODES
    fprintf(NodeID,'%d\t%d\t%d\t%d\t%d\t \n ',k,XX(k),YY(k),ATRBT(k),BC(k));
end
fclose(NodeID);

%%%% MROE734 EDITS Automatically writes the csv file now.
NodeID = fopen(csv_path,'w');
for k = 1:NUMBER_OF_NODES
    fprintf(NodeID,'%d\t%d\t%d\t%d\t%d\t \n',k,XX(k),YY(k),ATRBT(k),BC(k));
end
fclose(NodeID);

figure;
%%%% MROE734 EDITS renamed RGB_buthoole to RGB_mesh
RGB_mesh = zeros(size(cell_nodes));
RGB_mesh = RGB_mesh + cat(3,  (mito_nodes+tt_nodes)*1, (myo_nodes+tt_nodes)*1, (tt_nodes+glyco_nodes)*1);
%RGB = RGB + cat(3, (myo_nodes+tt_nodes)*1, (mito_nodes+tt_nodes)*1, (tt_nodes+glyco_nodes)*1);
imshow(RGB_mesh);
size(RGB_mesh)
