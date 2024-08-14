function [total_edges, GL_nodes, GL_labels, GL_edges, holes] = createGLPLC(glyco_bin, total_edges, spatres)
%CREATEGLPLC Creates PLC for the GL region
%
%   Args:
%   - glyco_bin, MxN logical -- binarised image of the region to analyse.
%   - total_edges, double -- number of edges.
%   - spatres, double -- spatial resolution.
%
%   Returns:
%   - total_edges, double -- number of edges and number of edges for the region.
%   - GL_nodes, Kx2 double -- coordinates of nodes for the region.
%   - GL_labels, Kx2 double -- labels for the nodes of the region.
%   - GL_edges, Kx2 double -- coordinates of edges for the region.
%   - holes, Kx2 double -- number of holes for the region.
holes = [];
GL_labels = [];
GL_nodes = [];
GL_edges = [];
%%%***ColeD changes following:
% I don't know why the value of 5000 was chosen as the area to remove
% pixels, but I found it larger than I liked for my images. I decided to
% lessen this to smaller amounts via experimentation
% I made glycogen's 'removal area' particularly small
%glyco_bin = bwareaopen(glyco_bin,5000); 
glyco_bin = bwareaopen(glyco_bin,2);%Para 2
se =  strel('disk',10);% Para 1
glyco_bin = imerode(glyco_bin,se);

%add a border
glyco_bin (1:size(glyco_bin),1:2) = 0;
glyco_bin (1:size(glyco_bin,1),size(glyco_bin,2)-1:size(glyco_bin,2)) = 0;
glyco_bin (1:2,1:size(glyco_bin,2)) = 0;
glyco_bin (size(glyco_bin,1)-1:size(glyco_bin,1),1:size(glyco_bin,2)) = 0;

[GL_labelled,GL_nums] = bwlabel(glyco_bin,8);

for label = 1:GL_nums
    [r,c] = find(GL_labelled==label);   %find row and columns beloonging to label
    glyco_bin = bwselect(GL_labelled,c,r,8);
    imshow(glyco_bin)
    glyco_bin = imdilate(glyco_bin,se);
    glyco_bin = imfill(glyco_bin,'holes');
    GL_boundary = edge(glyco_bin,'log','thinning');
    [GL_edgelist,~] = edgelink(GL_boundary);
    GL_edgepoints = cell2mat(GL_edgelist(1,1));
    GL_edgepoints = GL_edgepoints(1:end-1,:);
    GL_edgepoints_spatcoor = (GL_edgepoints-ones(size(GL_edgepoints)))*spatres;

    %add mitochondrial edges to edges list    
    %%%% MROE734 EDITS removed useless for loop and replaced with vectors
    tmp_edges = ((total_edges+1):(total_edges+size(GL_edgepoints,1)))';
    edges = [tmp_edges, tmp_edges, tmp_edges+1];
    edges(end,3) = total_edges+1;
    GL_edges = [GL_edges; edges];
    
    total_edges = total_edges+size(GL_edgepoints,1);
    %%%% MROE734 EDITS renamed all_nodes and all_labels for clarity
    GL_nodes = [GL_nodes;GL_edgepoints_spatcoor];
    GL_labels = [GL_labels;label*ones(size(GL_edgepoints_spatcoor))];
end
end