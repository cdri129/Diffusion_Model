function [total_edges, TT_nodes, TT_labels, TT_edges, holes] = createTTPLC(tt_bin, total_edges, spatres)
%CREATETTPLC Creates PLC for the TT region
%
%   Args:
%   - tt_bin, MxN logical -- binarised image of the region to analyse.
%   - total_edges, double -- number of edges.
%   - spatres, double -- spatial resolution.
%
%   Returns:
%   - total_edges, double -- number of edges and number of edges for the region.
%   - TT_nodes, Kx2 double -- coordinates of nodes for the region.
%   - TT_labels, Kx2 double -- labels for the nodes of the region.
%   - TT_edges, Kx2 double -- coordinates of edges for the region.
%   - holes, Kx2 double -- number of holes for the holes = [];
holes = [];
TT_labels = [];
TT_nodes = [];
TT_edges = [];


%***ColeD changes: This was used in original code
%se =  strel('disk',4);
% TT1 = imerode(tt_bin,se);
%***ColeD changes: end

% If you see ERROR: Blank / No mesh : hole has been placed outside a TT
% Try increasing the se value to 10

% %These following lines are purely for debugging and visualizing. Comment
% %out when not in use
%imshow(tt_bin);

% %These following lines are purely for debugging and visualizing. Comment
% %out when not in use
% % Display the binary image
% figure

%%%***ColeD changes following:
% I don't know why the value of 5000 was chosen as the area to remove
% pixels, but I found it larger than I liked for my images. I decided to
% lessen this to smaller amounts via experimentation
% TT1 = bwareaopen(TT1,5000); 



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
% imshow(tt_bin, 'InitialMagnification', 'fit'); % Display the original binary image (before)
% 
% % Perform morphological opening
% sed = strel('disk',4);
% sel = strel('line', 8, 45);
% tt_bin_after = bwareaopen(tt_bin,5)
% tt_bin_after = imdilate(tt_bin_after,sed);
% tt_bin_after = imopen(tt_bin_after,sed);
% tt_bin_after = imerode(tt_bin_after,sel);
% tt_bin_after = bwareaopen(tt_bin_after,5);
% % % tt_bin_after = imopen(tt_bin_after, se); % Apply morphological opening
% 
% % Overlay the before and after images
% hold on;
% h = imshow(cat(3, tt_bin, zeros(size(tt_bin)), zeros(size(tt_bin))), 'InitialMagnification', 'fit'); % Display before image in red channel
% set(h, 'AlphaData', 0.5); % Set transparency for before image
% h2 = imshow(cat(3, zeros(size(tt_bin_after)), tt_bin_after, zeros(size(tt_bin_after))), 'InitialMagnification', 'fit'); % Display after image in green channel
% set(h2, 'AlphaData', 0.5); % Set transparency for after image
% hold off;

% %***ColeD changes end

sed = strel('disk',3);
sel = strel('line', 8, 45);
tt_bin = bwareaopen(tt_bin,5);
tt_bin = imdilate(tt_bin,sed);
tt_bin = imopen(tt_bin,sed);
tt_bin = imerode(tt_bin,sel);
tt_bin = bwareaopen(tt_bin,5);


%add a border
tt_bin (1:size(tt_bin,1),1:2) = 0;
tt_bin (1:size(tt_bin,1),size(tt_bin,2)-1:size(tt_bin,2)) = 0;
tt_bin (1:2,1:size(tt_bin,2)) = 0;
tt_bin (size(tt_bin,1)-1:size(tt_bin,1),1:size(tt_bin,2)) = 0;

[TT_labelled,TT_nums] = bwlabel(tt_bin,8);
% Display the image
% figure;
% imshow(TT_labelled);

%***///ColeD changes: Don't forget to set label to 1 when done debugging
for label = 1:TT_nums
    [r,c] = find(TT_labelled==label);   %find row and columns beloonging to label
    tt = bwselect(TT_labelled,c,r,8);
    % figure;
    % imshow(tt)
    
    % %***ColeD changes start:
    % % attempting to debug the exact effects of the imdilate
    % % trying to see how line vs disk works and which is preferable
    % %
    %%%Debug code follows:
    % tt_pre_dilate = bwselect(TT_labelled,c,r,8);
    % %se = strel('line', 5, 45);
    % se = strel('disk', 1);
    % tt_dilate = imdilate(tt_pre_dilate,se);
    % 
    % % Overlay the before and after images
    % hold on;
    % h = imshow(cat(3, tt_pre_dilate, zeros(size(tt_pre_dilate)), zeros(size(tt_pre_dilate))), 'InitialMagnification', 'fit'); % Display before image in red channel
    % set(h, 'AlphaData', 0.5); % Set transparency for before image
    % h2 = imshow(cat(3, zeros(size(tt_dilate)), tt_dilate, zeros(size(tt_dilate))), 'InitialMagnification', 'fit'); % Display after image in green channel
    % set(h2, 'AlphaData', 0.5); % Set transparency for after image
    % hold off;

    % %***ColeD changes end

    se = strel('disk', 1);
    tt = imdilate(tt, se);
    % figure
    % imshow(tt);

    tt = imfill(tt,'holes');
    % figure;
    % imshow(tt)
    
    TT_boundary = edge(tt,'log','thinning');
    % figure;
    % imshow(TT_boundary)
    %***ColeD changes start:
    % I have found that when the previous 'edge' line finishes running,
    % some boundaries will have multiple smaller disconnected
    % boundaries/edges present
    % bwareaopen was being used to clean and discard these
    % floating edges. However, the area cleaned in the command seemed
    % arbitrary between images, leading to some disconnected edges
    % persisting, which causes the code to break completely and enter an
    % infinite loop at the 'edgelink' line. However, increasing the area to
    % be erased still misses some edges and erases entire, singular edges
    % sometimes.

    % To address this, I thought of instead adding code to find all edges,
    % calculate their area, and make a new image that keeps only the
    % largest area. This should leave the images with the smaller, singular
    % edges intact while still fixing the issue of having multiple
    % disconnected edges.

    %Removed bwareaopen line: TT_boundary = bwareaopen(TT_boundary, 25);
    % figure;
    % imshow(TT_boundary);

    %***New code starts here

    % Perform connected component analysis to identify separate edges
    CC = bwconncomp(TT_boundary);
    numEdges = CC.NumObjects;
    
    % If multiple edges are found, keep only the largest one
    if numEdges > 1
        % Calculate areas of all edges
        edgeAreas = cellfun(@numel, CC.PixelIdxList);
    
        % Find the index of the largest edge
        [~, largestEdgeIdx] = max(edgeAreas);
    
        % Create a binary image containing only the largest edge
        TT_boundary = false(size(TT_boundary));
        TT_boundary(CC.PixelIdxList{largestEdgeIdx}) = true;
    end
    % figure
    % imshow(TT_boundary)
    %***ColeD changes ends here:


    %***ColeD changes start:
    % checking all ~50+ labels can be tedious, especially when I check
    % 1-20, make a small change for 21+. Hopefully this will catch any
    % scenarios where I make a change to label 21's erosion if that erosion
    % proves to be too much for the previously checked labels.
    % Check if the image is all black (i.e., all pixels are zero)
    if sum(TT_boundary(:)) == 0
        msgbox(['Label ', num2str(label), ' is all black (all pixels are zero).'], 'All Black Label');
    end
    %***ColeD changes ends here:
    
    [TT_edgelist,~] = edgelink(TT_boundary);
    TT_edgepoints = cell2mat(TT_edgelist(1,1));
    TT_edgepoints = TT_edgepoints(1:end-1,:);
    TT_edgepoints_spatcoor = (TT_edgepoints-ones(size(TT_edgepoints)))*spatres;

    %add mitochondrial edges to edges list
    %%%% MROE734 EDITS removed useless for loop and replaced with vectors
    tmp_edges = ((total_edges+1):(total_edges+size(TT_edgepoints,1)))';
    edges = [tmp_edges, tmp_edges, tmp_edges+1];
    edges(end,3) = total_edges+1;
    TT_edges = [TT_edges; edges];    
    total_edges = total_edges+size(TT_edgepoints,1);
    %%%% MROE734 EDITS renamed all_nodes and all_labels for clarity
    TT_nodes = [TT_nodes;TT_edgepoints_spatcoor];
    TT_labels = [TT_labels;label*ones(size(TT_edgepoints_spatcoor))];
    % figure;
    % imshow(tt);
    se =  strel('disk',1);  %% Para 3
    BW = imerode(tt,se);
    % figure;
    % imshow(BW);
    % hold on;
    [r,c] = find(BW);
    

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
end