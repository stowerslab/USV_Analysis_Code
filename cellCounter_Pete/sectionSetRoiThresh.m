function [roiPatchArray, xMin, xMax, yMin, yMax, thresh, Ibw, cellCoordinates, roiCoord, reducedim, boundary] = sectionSetRoiThresh(maxImageBitValue, maxPixelValue, IchooseRoi, Isegment, roiCoord, ROIlength, ROIwidth, thresh)
%% sectionSetRoiThresh: function to read image and pick ROIs and set thresholds in order to quantify fluorescence
% Jason Keller, 5th Feb 2016
%
% This script does all the manual steps first (picking ROIs and thresholds), so that time-consuming processing can be done all at once
% 
% TO DO: 
%       reinstate using "regionprops" to get rid of regions < 50 pixel area
%       try a couple different methods (one oversegment, one undersegment) and average the results
%
% inputs:
%       [Ir Ig Ib] image including your region of interest; for better machine-readability, consider doing the following:
%           (1) Use soma/nuclear-expressing line or virus - then you could even use the control/GFP channel to compute cellular ROIs, no counterstain required
%           (2) Use NeuN or Nissl counterstain to look at only neurons
%           (3) Always normalize image exposure, for example, by saturating a certain % of pixels in each channel
%           (4) Image stacks are beter to get entire cell outlines, but can make cell overlap ambiguous (which is possible to clean up with watershed, but messy)
% %           (5) use a scope with an auto-focus feature and advanced image stitching capabilities to avoid bleaching lines, uneven section focus, etc. for larger images
%
% outputs:
% (A)   manual ROI masks for area of interest
% (B)   threshold
% (C)   BW segmented images with only regional maxima for fluorescence

%% ROI selection
% manually perform selection of ROIs using selected channel:
numRoi = 1;

roiPhrases = {'Select ROI #1 (PMC):'};
figure; 
hIm = imagesc(IchooseRoi); %create image handle and show image
set(gcf, 'Name',' Double-click to set ROI');
hAxes = gca;  % get handle to axes for 'waitfor' funtion below
colormap(gray); %set to grayscale
daspect([1 1 1])

% set(gcf, 'Position', get(0,'Screensize')); % Maximize figure - commented since this can stretch figure if not at same aspect ratio as monitor
% zoom on
% waitfor(hAxes,'CameraPosition');  % wait until zoom changes - note that you only have one chance to change ;-)

%% pete: RECTANGULAR ROI. change dimension of rectangle here.

if length(roiCoord) < 1
    h = imrect(gca, [0 0 ROIwidth ROIlength]); % [1600 1600] dimensions of the rectangle [xposition, yposition, width, length]
    setResizable(h,0); % not resizable
    rect = wait(h);
    % now move to appropriate position
    % command line blocked until rectangle is double clicked

%     I2 = imcrop(IchooseRoi, rect); % cropped image saved in I2, show with image(I2)
    % imagesc(I2);
    % hAxes = gca;
    % hIm = I2;
    % zoom off
    xposROI = round(rect(1)); % x position of ROI
    yposROI = round(rect(2));
else
    xposROI = roiCoord(1);
    yposROI = roiCoord(2);
    h = imrect(gca, [roiCoord(1) roiCoord(2) ROIwidth ROIlength]); % [1600 1600] dimensions of the rectangle [xposition, yposition, width, length]
    % comment these 4 lines if you don't want to confirm each separate ROI
%     setResizable(h,0); % not resizable
%     rect = wait(h);
%     xposROI = round(rect(1)); % x position of ROI
%     yposROI = round(rect(2));

end
roiCoord = [xposROI, yposROI];
%%

roiXLim = get(hAxes, 'XLim'); %save axis limits for later plotting
roiYLim = get(hAxes, 'YLim'); %save axis limits for later plotting

% PETE: do not change to uint8 if image is 8bit
% xMin = round(uint16(roiXLim(1))); 
% xMax = round(uint16(roiXLim(2)));
% yMin = round(uint16(roiYLim(1)));
% yMax = round(uint16(roiYLim(2)));
xMin = round(double(roiXLim(1))); 
xMax = round(double(roiXLim(2)));
yMin = round(double(roiYLim(1)));
yMax = round(double(roiYLim(2)));

cmap =  [0         0    1.0000;...
         0    0.5000         0;...
    1.0000         0         0;...
         0    0.7500    0.7500]; %colormap to ID different ROIs, first four lines from colormap(lines)

roiPatchArray = struct(   ...    %structure array containing data for each ROI selected in the image        
                                'hPatch',{},...          %Handles to ROI patch objects
                                'patchPosition',{},...   %Positions of patches, specified in pixel coordinates
                                'patchMask',{},...       %binary masks for each ROI patch (pixels inside patch = 1)                                                                                
                                'patchColor',{});        %color value of patch as index into 'lines' colomap

for i=1:numRoi
    set(gcf, 'Name',roiPhrases{i});
    % User selects ROI with mouse:
%     hRoi = imfreehand(gca);
    % PETE: rectangle from cropping is roi,
    hRoi = h;
    hMask = createMask(hRoi,hIm);
    pos = hRoi.getPosition;
    color = cmap(i,:);
    hPatch = patch('Parent',gca,'Visible','on',...  %draw colored patch to visualize ROI
                   'XData',pos(:,1),'YData',pos(:,2),...
                   'EdgeColor',color,'LineWidth',2,'FaceColor','none');%,...
%                  'FaceAlpha',0.15,'FaceColor',color);  %NOTE: something about specifying FaceAlpha really slows down computer, graphics rendering                

    %save ROIs:   
    roiPatchArray(end+1).hPatch = hPatch;
    roiPatchArray(end).patchPosition = pos;
    roiPatchArray(end).patchMask = hMask;  %access ex: imshow(roiPatchArray(1,1).patchMask)
    roiPatchArray(end).patchColor = i; %store color for later use
%     delete(hRoi); %delete imfreehand object
end
close(gcf); %close figure

%% Image segmentation
% now clean up / preprocess and then segment cell bodies, only for nuclear stain:
% FILTERING
%first use a median filter in [m n] neighborhood to remove salt-and-pepper noise:
Ifilt = Isegment;

num_boundaries_to_return = 3;

% Isegment(~roiPatchArray(1,1).patchMask) = 0;
% Isegment(~Isegment) = 0;

boundary = detectBoundary(Isegment, num_boundaries_to_return);

% Ifilt = medfilt2(Isegment, [3 3]); % use [m n] of ~cell size
% H = fspecial('gaussian', [3 3], 1); % [4 4]
% Ifilt = imfilter(Ifilt, H, 'replicate');  %also smooth out a bit to help with oversegmentation w/ regionalmax

% figure;
% imagesc(Ifilt)
% title("Test");

% figure; 
% middle = rot90(round(mean(roiPatchArray(1,1).patchPosition))); % as a spot check, plot a couple hundred pixels around the patch
% subplot(2,2,1); imagesc(Isegment(middle(1)-100:middle(1)+100,middle(2)-100:middle(2)+100)); title('original')
% subplot(2,2,2); imagesc(Ifilt(middle(1)-100:middle(1)+100,middle(2)-100:middle(2)+100)); title('median/gaussian filtered')

% THRESHOLD & MORPHOLOGICAL PREPROCESS
% now convert to binary using the interactive "threshTool" (needs 8-bit input) from MATLAB Central to manually choose a global threshold
imageToThresh = Ifilt;

% imageToThresh = Isegment;
% subimage = uint8(imageToThresh(middle(1)-100:middle(1)+100,middle(2)-100:middle(2)+100) .* (255/maxPixelValue));  %this needs to be scaled by max 16-bit value since the threshold will be eventually applied as such
% thresh = threshToolMod(subimage) / 255 * (maxPixelValue/maxImageBitValue);  % thresh must be normalized between 0 & 1 for 16-bit image
% thresh = 31/255; %this was manually chosen using tool

%% pete: change threshold here, but when you decide, use the SAME method for the rest of the dataset


% thresh = graythresh(imageToThresh) * 3; % use default Otsu method 




% %pete: use thresh gui to select thresh manually
% thresh = thresh_tool(imageToThresh) / 255; % divide by 255 to get percentage thresh

% thresh = 0.0; % pete: manual set percentage threshold to count, 0 to 1




%%
Ithresh = im2bw(imageToThresh, thresh);

% WATERSHED & MORPHOLOGICAL CLEAN-UP
% use distance transform and regional max
Ithresh(~roiPatchArray(1,1).patchMask) = 0; %only count within ROI!
Ithresh = bwareaopen(Ithresh,5,8); % pete: use 1 for low res, 5 usually. (using 3 for ec7-15) just check by eye and then use the same within dataset. removes from a binary image all connected components (objects) that have fewer than X pixels; 10 pixels seems to be nice for 10x 512px
% Idist = bwdist(~Ithresh); %compute distance transform (effectively makes small grayscale 'basins' around each cell)
Ifilt(~Ithresh) = 0;
% Ifilt(~Ithresh) = 0;
figure; 
reducedim = imagesc(Ifilt);
title('reduced');



%clean up bwdist output before watershed to prevent oversegmentation:
% Ibw = imregionalmax(Idist);  %assume 1 regional max per cell
Ibw = imregionalmax(Ifilt);

se = strel('disk',1); % pete: this is for dilating to prevent oversegmentation. I use 1. original = 3, this is the radius of the point
Ibw = imdilate(Ibw,se);
figure;
imagesc(Ibw);
title('dilated');
Ibw = bwmorph(Ibw,'shrink', Inf); %shrink to individual points

% figure; 
% imagesc(Ibw);
CC = bwconncomp(Ibw);
% CC.PixelIdxList;
cellCoordinates = {}; % 

for i = 1:length(CC.PixelIdxList)
    
    currCoord = CC.PixelIdxList{i};
    if length(currCoord) > 1 % sometimes will be list of connected components of an object. only count 1 of object by getting 1st idx
%         for j = 1:length(currCoord)
       [row,col] = ind2sub([yMax-1,xMax-1], currCoord(1));

       cellCoordinates{end+1} = [col - xposROI,row - yposROI];
%         end
    else
        [row,col] = ind2sub([yMax-1,xMax-1], currCoord);
%         disp('raw coord');
%         [row, col]
        cellCoordinates{end+1} = [col - xposROI,row - yposROI]; % store coordinate as coordinate wrt the rectangle ROI
%         cellCoordinates{end}
    end
%     xMax
%     yMax
%     yposROI
%     xposROI
end


if CC.NumObjects == 1
    [row,col] = ind2sub([xMax, yMax], CC.PixelIdxList{1});
end

if (CC.NumObjects == 1) && (row(1) == round(1+(xMax-1)/2)) && (col(1) == round((yMax-1)/2)) % conditions when no cells detected, counts as 1 cell with center coordinates
    cellCoordinates = {}; % reset
end

end %end of function

