function [total_edges, MY_nodes, MY_labels, MY_edges, holes] = createMYPLC(myo_bin, total_edges, spatres)
%CREATEMYPLC Creates PLC for the MY region
%
%   Args:
%   - myo_bin, MxN logical -- binarised image of the region to analyse.
%   - total_edges, double -- number of edges.
%   - spatres, double -- spatial resolution.
%
%   Returns:
%   - total_edges, double -- number of edges and number of edges for the region.
%   - MY_nodes, Kx2 double -- coordinates of nodes for the region.
%   - MY_labels, Kx2 double -- labels for the nodes of the region.
%   - MY_edges, Kx2 double -- coordinates of edges for the region.
%   - holes, Kx2 double -- number of holes for the region.
holes = [];
MY_labels = [];
MY_nodes = [];
MY_edges = [];

%%%***ColeD changes following:
% I don't know why the value of 5000 was chosen as the area to remove
% pixels, but I found it larger than I liked for my images. 
%myo_bin = bwareaopen(myo_bin,5000); 

myo_bin = bwareaopen(myo_bin,20);%Para 2
%imshow(myo_bin)


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
% imshow(myo_bin, 'InitialMagnification', 'fit'); % Display the original binary image (before)
% 
% % Perform morphological opening
% se = strel('line', 5, 45);
% %se = strel('disk', 5);
% myo_bin_after = imopen(myo_bin, se); % Apply morphological opening
% 
% % Overlay the before and after images
% hold on;
% h = imshow(cat(3, myo_bin, zeros(size(myo_bin)), zeros(size(myo_bin))), 'InitialMagnification', 'fit'); % Display before image in red channel
% set(h, 'AlphaData', 0.5); % Set transparency for before image
% h2 = imshow(cat(3, zeros(size(myo_bin_after)), myo_bin_after, zeros(size(myo_bin_after))), 'InitialMagnification', 'fit'); % Display after image in green channel
% set(h2, 'AlphaData', 0.5); % Set transparency for after image
% hold off;
% 
%***ColeD changes end

se = strel('line', 8, 45);
myo_bin = imopen(myo_bin, se);
%imshow(myo_bin)

%add a border
myo_bin (1:size(myo_bin),1:2) = 0;
myo_bin (1:size(myo_bin,1),size(myo_bin,2)-1:size(myo_bin,2)) = 0;
myo_bin (1:2,1:size(myo_bin,2)) = 0;
myo_bin (size(myo_bin,1)-1:size(myo_bin,1),1:size(myo_bin,2)) = 0;

[MY_labelled,MY_nums] = bwlabel(myo_bin,8);

for label = 1:MY_nums
    [r,c] = find(MY_labelled==label);   %find row and columns belonging to label
    myo_bin = bwselect(MY_labelled,c,r,8);
    % figure
    % imshow(myo_bin)
    
    % %***ColeD changes start:
    % % attempting to debug the exact effects of the imdilate
    % % trying to see how line vs disk works and which is preferable
    % %
    %%%Debug code follows:
    % %
    % myo_pre_dilate = bwselect(MY_labelled,c,r,8);
    % %se = strel('line', 5, 45);
    % se = strel('disk', 1);
    % myo_dilate = imdilate(myo_pre_dilate,se);
    % 
    % % Overlay the before and after images
    % hold on;
    % h3 = imshow(cat(3, zeros(size(myo_pre_dilate)), myo_pre_dilate, zeros(size(myo_pre_dilate))), 'InitialMagnification', 'fit'); % Display before image in green channel
    % set(h3, 'AlphaData', 0.5); % Set transparency for before image
    % h4 = imshow(cat(3, myo_dilate, zeros(size(myo_dilate)), zeros(size(myo_dilate))), 'InitialMagnification', 'fit'); % Display after image in red channel
    % set(h4, 'AlphaData', 0.5); % Set transparency for after image
    % hold off;

    % 
    % %***ColeD changes end

    se = strel('disk', 1);
    myo_bin = imdilate(myo_bin,se);
    % figure
    % imshow(myo_bin)

    %***ColeD changes start:
    % Getting rid of filling holes as the myofibrils are supposed to curve
    % around many structures
    % myo_bin = imfill(myo_bin,'holes');
    % figure
    % imshow(myo_bin)
    %***ColeD changes end

    % %***ColeD changes start:
    % % Running into errors with a tiny edge, attempting to delete edges that
    % % are too small
    % % Calculate connected components in the binary image
    % CC = bwconncomp(myo_bin);
    % 
    % % Iterate over each connected component
    % for label = 1:CC.NumObjects
    %     % Calculate the area of the current connected component
    %     area = numel(CC.PixelIdxList{label});
    % 
    %     % If the area is less than 10 pixels, set the entire image to black
    %     if area < 15
    %         myo_bin(:) = 0;  % Set all elements to 0 (black)
    %         break;  % Exit the loop since we've found a small area
    %     end
    % end
    % %***ColeD changes end

    MY_boundary = edge(myo_bin,'log','thinning');
    % figure
    % imshow(MY_boundary)

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
    % Removed bwareaopne code: MY_boundary = bwareaopen(MY_boundary, 100);
    %imshow(MY_boundary)

    %***New code starts here
    % Perform connected component analysis to identify separate edges
    CC = bwconncomp(MY_boundary);
    numEdges = CC.NumObjects;
    
    % If multiple edges are found, keep only the largest one
    if numEdges > 1
        % Calculate areas of all edges
        edgeAreas = cellfun(@numel, CC.PixelIdxList);
    
        % Find the index of the largest edge
        [~, largestEdgeIdx] = max(edgeAreas);
    
        % Create a binary image containing only the largest edge
        MY_boundary = false(size(MY_boundary));
        MY_boundary(CC.PixelIdxList{largestEdgeIdx}) = true;
    end
    % figure
    % imshow(MY_boundary)

    %***ColeD changes ends here:

    %***ColeD changes start here:
    % Attempting to overlay a specific label over the entire image to help
    % with debugging.
    %
    %***Debug code starts here
    % % Overlay the edge of label X over the original image
    % for X = 1:MY_nums
    %     close
    %     myo_pre_dilate = MY_labelled == 1;  % Image with label 1
    % 
    %     % Overlay the image containing label 1 with the current label (X)
    %     MY_edge_label_X = MY_labelled == X;
    % 
    %     % Plot the image with label 1 in green and the current label in red
    %     hold on;
    %     h5 = imshow(cat(3, zeros(size(myo_pre_dilate)), myo_pre_dilate, zeros(size(myo_pre_dilate))), 'InitialMagnification', 'fit'); % Display before edge in green channel
    %     set(h5, 'AlphaData', 0.5); % Set transparency for before edge
    % 
    %     h6 = imshow(cat(3, MY_edge_label_X, zeros(size(MY_edge_label_X)), zeros(size(MY_edge_label_X))), 'InitialMagnification', 'fit'); % Display edge of label X in red channel
    %     set(h6, 'AlphaData', 0.5); % Set transparency for edge of label X
    %     hold off;
    % end
    %***ColeD changes end


    %***ColeD changes start:
    % checking all ~50+ labels can be tedious, especially when I check
    % 1-20, make a small change for 21+. Hopefully this will catch any
    % scenarios where I make a change to label 21's erosion if that erosion
    % proves to be too much for the previously checked labels.
    % Check if the image is all black (i.e., all pixels are zero)
    if sum(MY_boundary(:)) == 0
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


    [MY_edgelist,~] = edgelink(MY_boundary);
    MY_edgepoints = cell2mat(MY_edgelist(1,1));
    MY_edgepoints = MY_edgepoints(1:end-1,:);
    MY_edgepoints_spatcoor = (MY_edgepoints-ones(size(MY_edgepoints)))*spatres;

    %add myofibril edges to edges list
    %%%% MROE734 EDITS removed useless for loop and replaced with vectors
    tmp_edges = ((total_edges+1):(total_edges+size(MY_edgepoints,1)))';
    edges = [tmp_edges, tmp_edges, tmp_edges+1];
    edges(end,3) = total_edges+1;
    MY_edges = [MY_edges; edges];

    total_edges = total_edges+size(MY_edgepoints,1);
    %%%% MROE734 EDITS renamed all_nodes and all_labels for clarity
    MY_nodes = [MY_nodes;MY_edgepoints_spatcoor];
    MY_labels = [MY_labels;label*ones(size(MY_edgepoints_spatcoor))];
end