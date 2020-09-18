function [roiPatchArray, Gfp, Tdt, GfpThresh, TdtThresh] = vmhQuant(filepathIn, mouseNum, sectionNum, roiPatchArrayInput)
%% vmhQuant: function to read image and quantify viral expression in the VMH
% Jason Keller, 5th Dec 2013
% 
% TO DO: 
%       MAKE a separate script to do all manual steps first: choosing ROIs and thresholds!!!!!
%       find out why watershed sometimes does not return obvious cells
%       reinstate using "regionprops" to get rid of regions < 50 pixel area
%       ***make "threshTool" into "adaptiveThreshTool" for GFP/TDT
%       try a couple different methods (one oversegment, one undersegment) and average the results
%       try using GFP and cellular ROIs (simpler, but prob. needs nuclear expression)
%
% inputs:
%       2D image including your region of interest; for better machine-readability, consider doing the following:
%           (1) Use Ai6 line, or other Nissl/nuclear-expressing line or virus - then you could even use the control/GFP channel to compute cellular ROIs, no counterstain required
%           (2) Use Nissl counterstain so that you will only look at neurons
%           (3) Always normalize image exposure, for example, by saturating a certain % of pixels in each channel
%           (4) Always take an image stack (ex. 3 sections spanning 15um) so that you get entire cell outline
%           (5) use a scope with an auto-focus feature and advanced image stitching capabilities to avoid bleaching lines, uneven section focus, etc.
% 
%       optional: roiPatchArrayInput, input ROIs from previous segmentation
%       (in format set below)
%
% outputs:
% (A)   manual ROI masks for (1) total VMH, (2) VMHdm, (3) VMHvl, and
%       (3) VMHc; then for each of these masks, compute the rest of the
%       metrics below:
% (B)   struct of GFP mean intensities in areas overlapping with segmented Hoechst cells
% (C)   struct of cFosCre TDT mean intensities in areas overlapping with segmented Hoechst cells

%% image read & format
% read nd2 directly to cut out the ElementsViewer step (can also access image metadata)
% using imreadBF function from MATLAB Central: http://www.mathworks.com/matlabcentral/fileexchange/32920-imread-for-multiple-life-science-image-file-formats 
% NOTE that Java heap size probably needs to be increased in MATLAB Preferences

% read in optional ROI mask input:
if nargin  > 3
    createRois = false;
    roiPatchArray = roiPatchArrayInput;  %probably should do some error checking here...
else
    createRois = true;
end

% Note that these are actually 12-bit images, so MATLAB will represent with
% uint16; however, the highest value (saturation) will be 2^12 - 1 = 4095,
% so we can just normalize by this to make a standard colormap

% possibly use imadjust/histeq/adapthisteq to fill 16 bits to optimize contrast at the 16-bit
% level, but this can also warp image intesities significantly

max16BitValue = 2^16-1;
maxPixelValue = 2^12-1;

filepathTotal = [filepathIn, mouseNum, '\', sectionNum, '.nd2'];
zplanes = 1; tframes = 1; channel = 1;  %only one channel can be imported at once
[vol] = imreadBF(filepathTotal,zplanes,tframes,channel); %read 1st channel in ND2 file - Hoechst
Ib = uint16(vol);
% Ib = imadjust(uint16(vol));
channel = 2;
[vol] = imreadBF(filepathTotal,zplanes,tframes,channel); %read 2nd channel in ND2 file - GFP
Ig = uint16(vol);
% meanGfp = mean(mean(Ig)); % perhaps use total image mean as threshold later
channel = 3;
[vol] = imreadBF(filepathTotal,zplanes,tframes,channel); %read 3rd channel in ND2 file - TDT
Ir = uint16(vol);
% meanTdt = mean(mean(Ir)); % perhaps use total image mean as threshold later
clear vol; %free up memory

% uncomment below to create and view RGB image to test:
% Irgb = zeros(size(Ib,1),size(Ib,2),3);
% Irgb(:,:,1) = double(Ir)./ double(max(max(Ir))); 
% Irgb(:,:,2) = double(Ig)./ double(max(max(Ig)));
% Irgb(:,:,3) = double(Ib)./ double(max(max(Ib)));
% imshow(Irgb);

%% ROI selection
% manually perform selection of 3 ROIs using the Hoechst channel: (VMHc can
% be computed)
if createRois
    numRoi = 3;
    roiPhrases = {'(1) Select total VMH ROI:'; '(2) Select VMH_dm ROI:'; '(3) Select VMH_vl ROI:'};
    figure; hIm = imagesc(Ib); %create image handle and show image
    hAxes = gca;  % get handle to axes for 'waitfor' funtion below
    colormap(gray); %set to grayscale
    % set(gcf, 'Position', get(0,'Screensize')); % Maximize figure - commented since this can stretch figure if not at same aspect ratio as monitor
    zoom on
    waitfor(hAxes,'CameraPosition');  % wait until zoom changes - note that you only have one chance to change ;-)
    zoom off
    vmhXLim = get(hAxes, 'XLim'); %save axis limits for later plotting
    vmhYLim = get(hAxes, 'YLim'); %save axis limits for later plotting
    xMin = round(uint16(vmhXLim(1)));
    xMax = round(uint16(vmhXLim(2)));
    yMin = round(uint16(vmhYLim(1)));
    yMax = round(uint16(vmhYLim(2)));

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
    %close figure
    % close(gcf);
    set(gcf, 'Name','all ROIs');
    
else %if not creating ROI, just set window around it:
    [xVal, yVal] = find(roiPatchArray(1, 1).patchMask);
    xMin = round(uint16(min(xVal)));
    xMax = round(uint16(max(xVal)));
    yMin = round(uint16(min(yVal)));
    yMax = round(uint16(max(yVal)));
end

% uncomment below to create and view RGB image to test:
% Irgb = zeros(size(Ib,1),size(Ib,2),3);
% Ir(~roiPatchArray(1,1).patchMask) = 0; Irgb(:,:,1) = double(Ir)./  double(max(max(Ir))); 
% Ig(~roiPatchArray(1,1).patchMask) = 0; Irgb(:,:,2) = double(Ig)./  double(max(max(Ig)));
% Ib(~roiPatchArray(1,1).patchMask) = 0; Irgb(:,:,3) = double(Ib)./  double(max(max(Ib)));
% imshow(Irgb);

%% Image segmentation
% now clean up / preprocess and then segment cell bodies, only for nuclear stain:
% FILTERING
%first use a median filter in [m n] neighborhood to remove salt-and-pepper noise:
Ibfilt = medfilt2(Ib, [3 3]);
H = fspecial('gaussian', [5 5], 1);
Ibfilt = imfilter(Ibfilt, H, 'replicate');  %also smooth out a bit to help with watershed oversegmentation
figure; 
middle = rot90(round(mean(roiPatchArray(1,1).patchPosition))); % as a spot check, plot a couple hundred pixels around the patch
subplot(2,2,1); imagesc(Ib(middle(1)-100:middle(1)+100,middle(2)-100:middle(2)+100)); title('original')
subplot(2,2,2); imagesc(Ibfilt(middle(1)-100:middle(1)+100,middle(2)-100:middle(2)+100)); title('median/gaussian filtered')

% THRESHOLD & MORPHOLOGICAL PREPROCESS
% now convert to binary based on nuclear stain; this is not trivial, and there are many methods of thresholding; one option is to use the
% interactive "thresh_tool" from MATLAB Central, or to use a local adaptive method to try and correct for imaging issues such as
% bleaching in a stitched image and uneven focus / illumination; this could probably work well if the imaging process is very controlled, but in the 
% absence of that I will try the interactive "thresh_tool" from MATLAB Central to manually choose a global threshold based on the VMH

% However, first preprocess using morphological operations, adapted from MATLAB 'Marker-Controlled Watershed Segmentation'
% example:
se = strel('disk', 3);  %make disk filter with 3 pixel radius
Ie = imerode(Ibfilt, se); %erode by filter
Iobr = imreconstruct(Ie, Ibfilt); 
Iobrd = imdilate(Iobr, se);
Iobrcbr = imreconstruct(imcomplement(Iobrd), imcomplement(Iobr));
Iobrcbr = imcomplement(Iobrcbr);
% figure; imagesc(Iobrcbr), title('Opening/closing-by-reconstruction (Iobr)')

%       Ibthresh = niblack(Ibfilt, [50 50], 0.1); %window=50 seems to be the sweet spot of ~1 nucleus; need window greater than cell radius, so that average will be pulled down by background
%       thresh = graythresh(Ibfilt);  %use automatic Otsu threshold

% make 8-bit smaller image to feed into threshTool:
subimage = uint8(Iobrcbr(middle(1)-100:middle(1)+100,middle(2)-100:middle(2)+100) .* (255/maxPixelValue));  %this needs to be scaled by max 16-bit value since the threshold will be eventually applied as such
thresh = threshToolMod(subimage) / 255 * (maxPixelValue/max16BitValue);  % thresh must be normalized between 0 & 1 for 16-bit image
% thresh = 31/255; %this was manually chose using tool
Ibthresh = im2bw(Iobrcbr,thresh);

% WATERSHED & MORPHOLOGICAL CLEAN-UP
% finally separate overlapping nuclei using the watershed algorithm, which only allows a single local minumum per watershed
% also, use morphological operations (ex. see 'Marker-Controlled Watershed Segmentation' in HELP of 'bwmorph') to clean up
% oversegmentation:
Ibthresh(~roiPatchArray(1,1).patchMask) = 0; %only count within VMH!!!
Ibthresh = bwareaopen(Ibthresh,50,8); % removes from a binary image all connected components (objects) that have fewer than X pixels; 50 pixels seems to be nice for Nikon 20x stitch
Ibdist = bwdist(~Ibthresh); %compute distance transform (effectively makes small grayscale 'basins' around each cell)
Ibdist = -Ibdist;  %preprocessing from MATLAB 'watershed' example code
Ibdist2 = Ibdist; %create a copy to manipulate basins and combine small ones (presumably where nuclear staining is punctate) to prevent oversegmentation
Ibdist2(~Ibthresh) = -Inf; %from MATLAB 'watershed' example code - set basins & background to trenches

% NOTE: since the watershed algorithm can oversegment and overestimate the number
% of nuclei, perhaps try without (i.e. use bwconncomp to label cells) as well and average results

%clean up bwdist output before watershed to prevent oversegmentation:
bw = imregionalmin(Ibdist);  %regional min predicts watershed - but each nucleus may have several local minima, so we will try to combine them
se2 = strel('disk', 3);  %we will dilate the regional min with a disk to connect entire nuclei
bw2 = imdilate(bw,se2);  %dilate to connect together very-regional minima within nuclei; note that there is a tradeoff here in that some cells that overlap a lot will be combined
bw3 = bwmorph(bw2, 'bridge'); %connect caddy-corner pixels; this should have little impact
bw4 = bwmorph(bw3,'thin', 2); %now thin back down to skeleton before setting the basins;  n=Inf takes a long time
Ibdist2(bw4) = -Inf;
subplot(2,2,3); imagesc(Ibdist2(middle(1)-100:middle(1)+100,middle(2)-100:middle(2)+100)); title('before watershed')

% now finally compute the watershed regions as our cellular ROIs:
Ibwater = watershed(Ibdist2,8);

% see "regionprops" for all options to further process the watershed output
%ex. get rid of all watersheds with area < 50 pixels
% stats = regionprops(Ibwater, 'Area');
% idxTooSmall = find([stats.Area] <= 50);
% % note that doing the above will mean that certain ROIs get NaN means - try "length(find(isnan(Gfp.total)))", for example
% Ibwater(ismember(Ibwater, idxTooSmall)) = 1; %set pixels in regions with small areas = 1 (background)

IbwaterRegions = label2rgb(Ibwater,'jet', 'w', 'shuffle');  %create colored regions surrounded by white
subplot(2,2,4); imshow(IbwaterRegions(middle(1)-100:middle(1)+100,middle(2)-100:middle(2)+100)); title('watershed')
% subplot(2,2,4); imagesc(Ibdist(middle(1)-100:middle(1)+100,middle(2)-100:middle(2)+100)); title('watershed'); % view bwdist image to see how watershed works...
colormap(jet);  %but maybe easier to see peaks/valleys with 'jet'
set(gcf, 'Position', get(0,'Screensize')); % Maximize figure
waitforbuttonpress; %display figure before moving on, to spot check segmentation for major errors

%% Cell counting
% for every label/region/nucleus found by the segmentation above,
% store in an array the mean GFP intensity in that region, and the mean TDT

% Note that the Ibwater region corresponding to '1' is the dead space between nuclei, so subtract 1
% idxLargeEnough = find([stats.Area] > 50);
% numNuclei = length(idxLargeEnough) - 1;
numNuclei = max(max(Ibwater)) - 1;

Gfp = struct(   ...    %structure array containing mean GFP intensity data for each cellular ROI found in the image        
                                'total',uint16(zeros(numNuclei,1)),...          %total VMH
                                'dm',uint16(zeros(numNuclei,1)),...   %VMHdm
                                'vl',uint16(zeros(numNuclei,1)));      %VMHvl                                                                                
                            
Tdt = struct(   ...    %structure array containing mean TDT intensity data for each cellular ROI found in the image        
                                'total',uint16(zeros(numNuclei,1)),...          %total VMH
                                'dm',uint16(zeros(numNuclei,1)),...   %VMHdm
                                'vl',uint16(zeros(numNuclei,1)));      %VMHvl     

IgDm = Ig; IgDm(~roiPatchArray(1,2).patchMask) = 0; %second ROI mask in for DM - set pixels outside to zero
IgVl = Ig; IgVl(~roiPatchArray(1,3).patchMask) = 0; %third ROI mask in for VL - set pixels outside to zero
IrDm = Ir; IrDm(~roiPatchArray(1,2).patchMask) = 0; 
IrVl = Ir; IrVl(~roiPatchArray(1,3).patchMask) = 0;

% Note that the Ibwater region corresponding to '1' is the dead space between nuclei, so only take (2:end)                          
GfpStatsTotal = regionprops(Ibwater, Ig, 'MeanIntensity');
GfpStatsDm = regionprops(Ibwater, IgDm, 'MeanIntensity');
GfpStatsVl = regionprops(Ibwater, IgVl, 'MeanIntensity');
Gfp.total = [GfpStatsTotal(2:end).MeanIntensity];
Gfp.dm = [GfpStatsDm(2:end).MeanIntensity];
Gfp.vl = [GfpStatsVl(2:end).MeanIntensity];

TdtStatsTotal = regionprops(Ibwater, Ir, 'MeanIntensity');
TdtStatsDm = regionprops(Ibwater, IrDm, 'MeanIntensity');
TdtStatsVl = regionprops(Ibwater, IrVl, 'MeanIntensity');
Tdt.total = [TdtStatsTotal(2:end).MeanIntensity];
Tdt.dm = [TdtStatsDm(2:end).MeanIntensity];
Tdt.vl = [TdtStatsVl(2:end).MeanIntensity];

% ex. 'Gfp-positive' cells are those that are one standard deviation above mean intensity:
% numGfpDm = length(find(Gfp.dm > (mean(Gfp.dm)+std(Gfp.dm))));
% numTdtDm = length(find(Tdt.dm > (mean(Tdt.dm)+std(Tdt.dm))));
% numGfpTotal = length(find(Gfp.total > (mean(Gfp.total))));
% efficiencyGfp = numGfpTotal/double(numNuclei);

%% Choose counting thresholds
% the idea here is to apply a mean filter about the size of a nucleus, such
% that when choosing the threshold, the user should select a level that
% will catch entire nuclei, since we will take the mean value of a cellular
% region of interest in deciding if it is GFP/TDT positive

% simpler alternatives include:
%     gfpThresh = (mean(Gfp.total(Gfp.total>0)));
%     tdtThresh = (mean(Tdt.total(Tdt.total>0)));
%     etc.

% for each channel (GFP & TDT), we make a filtered, 8-bit version. GFP first:
hNucleusDisk = fspecial('disk', 5);  % this is about a nucleus size @ 20x
Igfilt = imfilter(Ig, hNucleusDisk, 'replicate'); 
% make 8-bit smaller image to feed into threshTool:
% subimageIg = uint8(Igfilt .* (255/maxPixelValue));
subimageIg = uint8(Igfilt(yMin:yMax, xMin:xMax) .* (255/maxPixelValue));
GfpThresh = threshToolMod(subimageIg) / 255*maxPixelValue; % convert back to 12-bit relevance

% TDT:
Irfilt = imfilter(Ir, hNucleusDisk, 'replicate'); 
% make 8-bit smaller image to feed into threshTool:
% subimageIr = uint8(Irfilt .* (255/maxPixelValue));
subimageIr = uint8(Irfilt(yMin:yMax, xMin:xMax) .* (255/maxPixelValue));
TdtThresh = threshToolMod(subimageIr) / 255*maxPixelValue; % convert back to 12-bit relevance

end %end of function

