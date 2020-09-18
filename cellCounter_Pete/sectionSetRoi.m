function [roiPatchArray, xMin, xMax, yMin, yMax] = sectionSetRoi(IchooseRoi)
% sectionSetRoiThresh: function to read image and pick ROIs and set thresholds in order to quantify fluorescence
% Jason Keller, 5th Feb 2016
%
% inputs:
%       [I] image including your region of interest; for better machine-readability, consider doing the following:
%           (1) Use soma/nuclear-expressing line or virus - then you could even use the control/GFP channel to compute cellular ROIs, no counterstain required
%           (2) Use NeuN or Nissl counterstain to look at only neurons
%           (3) Always normalize image exposure, for example, by saturating a certain % of pixels in each channel
%           (4) Image stacks are beter to get entire cell outlines, but can make cell overlap ambiguous (which is possible to clean up with watershed, but messy)
%           (5) use a scope with an auto-focus feature and advanced image stitching capabilities to avoid bleaching lines, uneven section focus, etc. for larger images
%
% outputs:
% (A)   manual ROI masks for area of interest
% (B)   ROI limits in X & Y
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

% ROI selection
% manually perform selection of ROIs using selected channel:
numRoi = 1;
roiPhrases = {'Select ROI #1:'};
figure; hIm = imagesc(IchooseRoi); %create image handle and show image
set(gcf, 'Name','Zoom in first:');
hAxes = gca;  % get handle to axes for 'waitfor' funtion below
colormap(gray); %set to grayscale
% set(gcf, 'Position', get(0,'Screensize')); % Maximize figure - commented since this can stretch figure if not at same aspect ratio as monitor
zoom on
waitfor(hAxes,'CameraPosition');  % wait until zoom changes - note that you only have one chance to change ;-)
zoom off
roiXLim = get(hAxes, 'XLim'); %save axis limits for later plotting
roiYLim = get(hAxes, 'YLim'); %save axis limits for later plotting
xMin = round(uint16(roiXLim(1)));
xMax = round(uint16(roiXLim(2)));
yMin = round(uint16(roiYLim(1)));
yMax = round(uint16(roiYLim(2)));

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
    hRoi = imfreehand(gca);
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
close gcf;
end %end of function

