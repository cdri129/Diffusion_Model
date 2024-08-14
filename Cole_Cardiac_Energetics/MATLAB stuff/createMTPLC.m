function [total_edges, MT_nodes, MT_labels, MT_edges, holes] = createMTPLC(mito_bin, total_edges, spatres)
%CREATEMTPLC Creates PLC for the MT region
%
%   Args:
%   - mito_bin, MxN logical -- binarised image of the region to analyse.
%   - total_edges, double -- number of edges.
%   - spatres, double -- spatial resolution.
%
%   Returns:
%   - total_edges, double -- number of edges and number of edges for the region.
%   - MT_nodes, Kx2 double -- coordinates of nodes for the region.
%   - MT_labels, Kx2 double -- labels for the nodes of the region.
%   - MT_edges, Kx2 double -- coordinates of edges for the region.
%   - holes, Kx2 double -- number of holes for the region.
holes = [];
MT_labels = [];
MT_nodes = [];
MT_edges = [];


%%%***ColeD changes following:
% I don't know why the value of 5000 was chosen as the area to remove
% pixels, but I found it larger than I liked for my images. 
%mito_bin = bwareaopen(mito_bin,5000); 

mito_bin = bwareaopen(mito_bin,100); %Para 2
% figure
% imshow(mito_bin)


% %***ColeD changes start:
% % Here begins a section for debugging the effects of 'se' and how it effects the image
% 
% Noticed previous code seemed to make borders too large or small. Your
% results may vary, so be sure to test and check
%
%%% Debug code follows

% % Before opening
% figure;
% close;
% imshow(mito_bin, 'InitialMagnification', 'fit'); % Display the original binary image (before)

% Perform morphological opening test 1
% % % se = strel('line', 10, 45);
% % % mito_bin_after = imopen(mito_bin, se); % Apply morphological opening

% Perform morphological opening test 2
% sed = strel('disk',4);
% sel = strel('line', 8, 45);
% mito_bin_after = bwareaopen(mito_bin,5);
% mito_bin_after = imdilate(mito_bin_after,sed);
% mito_bin_after = imopen(mito_bin_after,sed);
% mito_bin_after = imerode(mito_bin_after,sel);
% mito_bin_after = bwareaopen(mito_bin_after,5);
% 
% 
% % Overlay the before and after images
% hold on;
% h = imshow(cat(3, mito_bin, zeros(size(mito_bin)), zeros(size(mito_bin))), 'InitialMagnification', 'fit'); % Display before image in red channel
% set(h, 'AlphaData', 0.5); % Set transparency for before image
% h2 = imshow(cat(3, zeros(size(mito_bin_after)), mito_bin_after, zeros(size(mito_bin_after))), 'InitialMagnification', 'fit'); % Display after image in green channel
% set(h2, 'AlphaData', 0.5); % Set transparency for after image
% hold off;
% 
% %***ColeD changes end

% se = strel('line', 10, 45);
% mito_bin = imopen(mito_bin, se);
% figure
% imshow(mito_bin)

% %***ColeD changes start
% % I think this previous disk version for erosion was removing too many
% % narrow sections of structures. Line is now being used
% %
%%%Removed code that was previously used:
% 50 seems to be too aggressive and leaves the resulting mito-bin image
% incredibly small
% %se =  strel('disk',50);% Para 1
% se =  strel('disk',10);
% mito_bin = imopen(mito_bin, se);
% figure
% imshow(mito_bin)

%mito_bin = imerode(mito_bin,se);
%imshow(mito_bin)
% %***ColeD changes end


sed = strel('disk',6);
sel = strel('line', 8, 45);
sed2 = strel('disk',5);
mito_bin = bwareaopen(mito_bin,5);
mito_bin = imdilate(mito_bin,sed);
mito_bin = imerode(mito_bin,sed2);
mito_bin = imopen(mito_bin,sed);
mito_bin = imerode(mito_bin,sel);
mito_bin = bwareaopen(mito_bin,20);

% figure;
% imshow(mito_bin);

%add a border
mito_bin (1:size(mito_bin),1:2) = 0;
mito_bin (1:size(mito_bin,1),size(mito_bin,2)-1:size(mito_bin,2)) = 0;
mito_bin (1:2,1:size(mito_bin,2)) = 0;
mito_bin (size(mito_bin,1)-1:size(mito_bin,1),1:size(mito_bin,2)) = 0;

[MT_labelled,MT_nums] = bwlabel(mito_bin,8);

for label = 1:MT_nums
    [r,c] = find(MT_labelled==label);   %find row and columns beloonging to label
    mito_bin = bwselect(MT_labelled,c,r,8);
    % figure;
    % imshow(mito_bin)

    % %***ColeD changes start:
    % % attempting to debug the exact effects of the imdilate
    % % trying to see how line vs disk works and which is preferable
    % %
    %%%Debug code follows:
    %
    % mito_pre_dilate = bwselect(MT_labelled,c,r,8);
    % %se = strel('line', 5, 45);
    % se = strel('disk', 1);
    % mito_dilate = imdilate(mito_pre_dilate,se);
    % 
    % % Overlay the before and after images
    % hold on;
    % h = imshow(cat(3, mito_pre_dilate, zeros(size(mito_pre_dilate)), zeros(size(mito_pre_dilate))), 'InitialMagnification', 'fit'); % Display before image in red channel
    % set(h, 'AlphaData', 0.5); % Set transparency for before image
    % h2 = imshow(cat(3, zeros(size(mito_dilate)), mito_dilate, zeros(size(mito_dilate))), 'InitialMagnification', 'fit'); % Display after image in green channel
    % set(h2, 'AlphaData', 0.5); % Set transparency for after image
    % hold off;
    % 
    % %***ColeD changes end

    
    se = strel('disk', 1);
    mito_bin = imdilate(mito_bin,se);
    %imshow(mito_bin)

    mito_bin = imfill(mito_bin,'holes');
    %imshow(mito_bin)

    MT_boundary = edge(mito_bin,'log','thinning');
    %imshow(MT_boundary)

    %***ColeD changes start:
    % noticed that after thinning, small disconnected
    % edges are leftover, added bwareaopen to remove any leftover. However, the
    % amount of area removed may need to be changed depending on one's own
    % images, so a new change was made
    %
    % Now attempting to remove this bwareaopen segment in favor of code
    % that will instead find all edges, see if they are disconnected, and
    % then remove all but the largest edge. This will hopefully work better
    % to remove these disconnected edges (which break the code if present)
    % without deleting singular edges that are smaller than the arbitrarily
    % chosen area
    % Removed bwareaopen code: MT_boundary = bwareaopen(MT_boundary, 100);
    %imshow(MT_boundary)

    %***New code starts here
    % Perform connected component analysis to identify separate edges
    CC = bwconncomp(MT_boundary);
    numEdges = CC.NumObjects;
    
    % If multiple edges are found, keep only the largest one
    if numEdges > 1
        % Calculate areas of all edges
        edgeAreas = cellfun(@numel, CC.PixelIdxList);
    
        % Find the index of the largest edge
        [~, largestEdgeIdx] = max(edgeAreas);
    
        % Create a binary image containing only the largest edge
        MT_boundary = false(size(MT_boundary));
        MT_boundary(CC.PixelIdxList{largestEdgeIdx}) = true;
    end
    % figure
    % imshow(MT_boundary)
    %***ColeD changes ends here:


    %***ColeD changes start:
    % checking all ~50+ labels can be tedious, especially when I check
    % 1-20, make a small change for 21+. Hopefully this will catch any
    % scenarios where I make a change to label 21's erosion if that erosion
    % proves to be too much for the previously checked labels.
    % Check if the image is all black (i.e., all pixels are zero)
    if sum(MT_boundary(:)) == 0
        msgbox(['Label ', num2str(label), ' is all black (all pixels are zero).'], 'All Black Label');
    end
    %***ColeD changes ends here:
    
    %***ColeD changes start:
    % When debugging I found that keeping figures open during each loop got
    % tedious. Here is some code to close the last couple of figures
    % Get the list of all figure handles
    fig_handles = get(groot, 'Children');
    
    % Determine the number of figures
    num_figures = numel(fig_handles);
    
    % Close the last X figures to prevent clutter
    num_last_figures_to_keep = 0;
    if num_figures > num_last_figures_to_keep
        close(fig_handles(1:num_figures - num_last_figures_to_keep));
    end
    %***ColeD changes end
    

    [MT_edgelist, ~] = edgelink(MT_boundary);
    MT_edgepoints = cell2mat(MT_edgelist(1,1));
    MT_edgepoints = MT_edgepoints(1:end-1,:);
    MT_edgepoints_spatcoor = (MT_edgepoints-ones(size(MT_edgepoints)))*spatres;

    %add mitochondrial edges to edges list
    %%%% MROE734 EDITS removed useless for loop and replaced with vectors
    tmp_edges = ((total_edges+1):(total_edges+size(MT_edgepoints,1)))';
    edges = [tmp_edges, tmp_edges, tmp_edges+1];
    edges(end,3) = total_edges+1;
    MT_edges = [MT_edges; edges];

    total_edges = total_edges+size(MT_edgepoints,1);
    %%%% MROE734 EDITS renamed all_nodes and all_labels for clarity
    MT_nodes = [MT_nodes;MT_edgepoints_spatcoor];
    MT_labels = [MT_labels;label*ones(size(MT_edgepoints_spatcoor))];
end
end