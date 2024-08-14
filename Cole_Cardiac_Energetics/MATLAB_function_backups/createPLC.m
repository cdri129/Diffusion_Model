function [total_edges, edges, nodes, labels, holes] = createPLC(img_bin, edges, total_edges, spatres)
%CREATEPLC Creates PLC for a given region
%
%   Args:
%   - img_bin, MxN logical -- binarised image of the region to analyse.
%   - total_edges, double -- number of edges.
%   - spatres, double -- spatial resolution.
%
%   Returns:
%   - total_edges, double -- number of edges and number of edges for the region.
%   - nodes, Kx2 double -- coordinates of nodes for the region.
%   - labels, Kx2 double -- labels for the nodes of the region
%   - holes, Kx2 double -- number of holes for the region
nodes = [];
labels = [];
holes = [];
se =  strel('disk',4);

% If you see ERROR: Blank / No mesh : hole has been placed outside a TT
% Try increasing the se value to 10
eroded_region = imerode(img_bin,se);
opened_region = bwareaopen(img_bin,5000); %Para 2
region = imdilate(opened_region,se);

%add a border
region(1:size(region),1:2) = 0;
region(1:size(region,1),size(region,2)-1:size(region,2)) = 0;
region(1:2,1:size(region,2)) = 0;
region(size(region,1)-1:size(region,1),1:size(region,2)) = 0;

[labelled, nums] = bwlabel(region,8);

for label = 1:nums
    [r,c] = find(labelled==label);   %find row and columns belonging to label
    ROI = bwselect(labelled,c,r,8);
    ROI = imdilate(ROI, se);
    ROI = imfill(ROI,'holes');
    boundary = edge(ROI,'log','thinning');
    [edgelist,~] = edgelink(boundary);
    edgepoints = cell2mat(edgelist(1,1));
    edgepoints = edgepoints(1:end-1,:);
    edgepoints_spatcoor = (edgepoints-ones(size(edgepoints)))*spatres;

    %add mitochondrial edges to edges list
    for k = (total_edges+1):(total_edges+size(edgepoints,1)-1)
        edges(k,1) = k;
        edges(k,2) = k;
        edges(k,3) = k+1;
    end

    edges(k+1,1) = k+1;
    edges(k+1,2) = k+1;
    edges(k+1,3) = total_edges+1;
    total_edges = total_edges+size(edgepoints,1);
    nodes = [nodes;edgepoints_spatcoor];
    labels = [labels;label*ones(size(edgepoints_spatcoor))];
    se =  strel('disk',1);  %% Para 3
    BW = imerode(ROI,se);
    [r,c] = find(BW);

    % Check if there are any points found
    if ~isempty(r) && ~isempty(c)
        pixel = [r(round(size(r, 1) / 2)), c(round(size(c, 1) / 2))];
        points = (pixel - ones(size(pixel))) * spatres;
        holes = [holes; points];
    else
        % Handle the case where no points are found (e.g., set a default value)
        % You might want to add specific handling based on your requirements
        disp('No points found.');
    end

    progressbar([],label/nums-0.0001, [], [], [], [], []);
end

progressbar(2/11-0.0001,[], [], [], [], [], []);

end