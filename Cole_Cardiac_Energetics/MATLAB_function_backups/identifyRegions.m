function [region_size, region_bin] = identifyRegions(img, val)
%IDENTIFYREGIONS Identify regions in an image based on the intensity value
%
%   Args:
%   - img, MxN double -- image to analyse.
%   - val, double -- value to use for region segmentation.
%
%   Returns:
%   - region_size, double -- number of pixel in the region.
%   - region_bin, MxN double -- binarised version of the image where 1 is
%   the desired region and 0 is the rest of the image.
region_size = sum(img(:) == val);
region_bin = img; % Duplicate the image
% Binarise using the value as a threshold
region_bin(region_bin ~= val) = 0;
region_bin(region_bin == val) = 1;
end