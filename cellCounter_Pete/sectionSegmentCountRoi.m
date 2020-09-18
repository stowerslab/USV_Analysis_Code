function [overlayAll, Ibw, CC] = sectionSegmentCountRoi(Isegment)
%% sectionSegmentCountRoi: function to read image and segment in order to quantify fluorescence
% Jason Keller, 13th Feb 2017
%
% First "setSectionRoi" script does all the manual steps first, so that time-consuming processing can be done all at once
% 
% TO DO: 
%       reinstate using "regionprops" to get rid of regions < 50 pixel area
%       try a couple different methods (one oversegment, one undersegment) and average the results
%
% inputs:
%       [I] image alread cropped to ROI & masked, and converted to normalized double; for better machine-readability, consider doing the following:
%           (1) Use soma/nuclear-expressing line or virus - then you could even use the control/GFP channel to compute cellular ROIs, no counterstain required
%           (2) Use NeuN or Nissl counterstain to look at only neurons
%           (3) Always normalize image exposure, for example, by saturating a certain % of pixels in each channel
%           (4) Image stacks are beter to get entire cell outlines, but can make cell overlap ambiguous (which is possible to clean up with watershed, but messy)
%           (5) use a scope with an auto-focus feature and advanced image stitching capabilities to avoid bleaching lines, uneven section focus, etc. for larger images
%
% outputs:
% (A)   threshold
% (B)   BW segmented images with only regional maxima for fluorescence
% (C)   connected components on the BW image


%% Image segmentation
% now clean up / preprocess and then segment cell bodies, only for nuclear stain:
% FILTERING
%first use a median filter in [m n] neighborhood to remove salt-and-pepper noise:
cellSize = [10 10]; % use [m n] of ~cell size
Ifilt = medfilt2(Isegment, cellSize); 
H = fspecial('gaussian', cellSize, 1);
Ifilt = imfilter(Ifilt, H, 'replicate');  %also smooth out a bit to help with oversegmentation w/ regionalmax
% figure; 
% subplot(2,2,1); imagesc(Isegment); title('original')
% subplot(2,2,2); imagesc(Ifilt); title('median/gaussian filtered')

% THRESHOLD & MORPHOLOGICAL PREPROCESS
% now convert to binary using the interactive "threshTool" (needs 8-bit input) from MATLAB Central to manually choose a global threshold
imageToThresh = Ifilt;
% imageToThresh = Isegment;
% subimage = uint8(imageToThresh .* 255);  %this needs to be scaled
% thresh = threshToolMod(subimage) / 255 * (maxPixelValue/maxImageBitValue);  % thresh must be normalized between 0 & 1 for 16-bit image
% thresh = 31/255; %this was manually chosen using tool
thresh = graythresh(imageToThresh); % use default Otsu method
Ithresh = im2bw(imageToThresh, thresh);
% figure; imshow(Ithresh);

% Find cells &  MORPHOLOGICAL CLEAN-UP using regional max
Ithresh = bwareaopen(Ithresh,10,8); % removes from a binary image all connected components (objects) that have fewer than X pixels; 10 pixels seems to be nice for 10x 512px
Ifilt(~Ithresh) = 0;
% figure; imagesc(Ifilt); colormap(gray);
Ibw = imregionalmax(Ifilt); %assume 1 regional max per cell %clean up bwdist output before watershed to prevent oversegmentation
se = strel('disk',1); 
Ibw = imdilate(Ibw,se); %expand single regional max points so we can see them

figure; %QC - plot found points on top of original image
overlayR = Isegment; %first make grayscale part
overlayG = Isegment;
overlayB = Isegment;
overlayDot = [1 1 0];  % use yellow
overlayR(Ibw) = repmat(overlayDot(1), length(find(Ibw)), 1); %now color overlay
overlayG(Ibw) = repmat(overlayDot(2), length(find(Ibw)), 1);
overlayB(Ibw) = repmat(overlayDot(3), length(find(Ibw)), 1);
overlayAll = cat(3, overlayR, overlayG, overlayB);
imshow(overlayAll);

Ibw = bwmorph(Ibw,'shrink', Inf); %shrink to individual points
CC = bwconncomp(Ibw);

end %end of function

