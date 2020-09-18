function [roiPatchArray, Gfp, Tdt] = vmhQuant(filepathIn, mouseNum, sectionNum)
%% vmhQuant: function to read image and quantify viral expression in the VMH
% Jason Keller, 19th Nov 2013
% 
% TO DO: make input option for previous ROI masks
%
% inputs:
%       stiched image (from Nikon C2) in native nd2 format; it would
%       probably be best if this was a max projection over a few sections in a ~20um area,
%       since we want to make sure and capture a nice, bright nucleus for
%       segmentation
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

% Note that NikonC2 images are actually 12-bit images, so I represent with
% uint16; however, the highest value (saturation) will be 2^12 - 1 = 4095,
% so we can just normalize by this to make a standard colormap below
% also maybe use 'imadjust' to fill 16 bits
filepathTotal = [filepathIn, mouseNum, '\', sectionNum, '.nd2'];
zplanes = 1; tframes = 1; channel = 1;  %only one channel can be imported at once
[vol] = imreadBF(filepathTotal,zplanes,tframes,channel); %read 1st channel in ND2 file - Hoechst
Ib = uint16(vol);
channel = 2;
[vol] = imreadBF(filepathTotal,zplanes,tframes,channel); %read 2nd channel in ND2 file - GFP
Ig = uint16(vol);
channel = 3;
[vol] = imreadBF(filepathTotal,zplanes,tframes,channel); %read 3rd channel in ND2 file - TDT
Ir = uint16(vol);
clear vol; %free up memory

% uncomment below to create and view RGB image to test:
% bitDepth = 12;
% maxValue = 2^bitDepth - 1;
% Irgb = zeros(size(Ib,1),size(Ib,2),3);
% Irgb(:,:,1) = Ir;
% Irgb(:,:,2) = Ig;
% Irgb(:,:,3) = Ib;
% Irgb = double(Irgb)./maxValue;
% imshow(Irgb);

%% ROI selection
% manually perform selection of 3 ROIs using the Hoechst channel: (VMHc can
% be computed)
numRoi = 3;
roiPhrases = {'(1) Select total VMH ROI:'; '(2) Select VMH_dm ROI:'; '(3) Select VMH_vl ROI:'};
figure; hIm = imagesc(Ib); %create image handle and show image
hAxes = gca;  % get handle to axes for 'waitfor' funtion below
colormap(gray); %set to grayscale
% set(gcf, 'Position', get(0,'Screensize')); % Maximize figure - commented since this can stretch figure if not at same aspect ratio as monitor
zoom on
waitfor(hAxes,'CameraPosition');  % wait until zoom changes - note that you only have one chance to change ;-)
zoom off

% TO DO: add option to zoom in here first

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
set(gcf, 'Position', get(0,'Screensize')); % Maximize figure

% THRESHOLD
% now convert to binary based on nuclear stain; this is not trivial, and there are many methods of thresholding; one option is to use the
% interactive "thresh_tool" from MATLAB Central, or to use a local adaptive method to try and correct for imaging issues such as
% bleaching in a stitched image and uneven focus / illumination; this could probably work well if the imaging process is very controlled, but in the 
% absence of that I will try the interactive "thresh_tool" from MATLAB Central to manually choose a global threshold based on the VMH
% Ibthresh = niblack(Ibfilt, [50 50], 0.1); %window=50 seems to be the sweet spot of ~1 nucleus; need window greater than cell radius, so that average will be pulled down by background
% thresh = graythresh(Ibfilt);  %use automatic Otsu threshold
subimage = Ibfilt(middle(1)-100:middle(1)+100,middle(2)-100:middle(2)+100);
thresh = thresh_tool(subimage);
Ibthresh = im2bw(Ibfilt,thresh);
subplot(2,2,3); imshow(Ibthresh(middle(1)-100:middle(1)+100,middle(2)-100:middle(2)+100)); title('thresholded')

% WATERSHED
% finally separate overlapping nuclei using the watershed algorithm, which only allows a single local minumum per watershed
%try morphological operations (ex. imerode, imdilate, imopen, bwareaopen) to clean up first to
%prevent small watersheds:
Ibthresh(~roiPatchArray(1,1).patchMask) = 0; %only count within VMH!!!
Ibthresh = bwareaopen(Ibthresh,100,8); % removes from a binary image all connected components (objects) that have fewer than X pixels
Ibdist = bwdist(~Ibthresh); %compute distance transform (efectively makes small grayscale 'basins' around each cell)
Ibdist = -Ibdist;  %preprocessing from MATLAB 'watershed' example code
Ibdist(~Ibthresh) = -Inf; %from MATLAB 'watershed' example code
% TO DO: the watershed algorithm can oversegment and overestimate the number
% of nuclei, so perhaps try without (i.e. use bwconncomp to label cells) as well and average results
Ibwater = watershed(Ibdist,8);
IbwaterRegions = label2rgb(Ibwater,'jet', 'w', 'shuffle');  %create colored regions surrounded by white
subplot(2,2,4); imshow(IbwaterRegions(middle(1)-100:middle(1)+100,middle(2)-100:middle(2)+100)); title('watershed')
% subplot(2,2,4); imagesc(Ibdist(middle(1)-100:middle(1)+100,middle(2)-100:middle(2)+100)); title('watershed'); % view bwdist image to see how watershed works...
colormap(jet);  %but maybe easier to see peaks/valleys with 'jet'
waitforbuttonpress; %display figure before moving on
% perhaps use "regionprops" to further process the watershed output, but
% for now...

%% Cell counting
%algorithm: for every label/region/nucleus found by the segmentation above,
% store in an array the mean GFP intensity in that region, and the mean TDT
numNuclei = max(max(Ibwater));
Gfp = struct(   ...    %structure array containing mean GFP intensity data for each cellular ROI found in the image        
                                'total',uint16(zeros(numNuclei,1)),...          %total VMH
                                'dm',uint16(zeros(numNuclei,1)),...   %VMHdm
                                'vl',uint16(zeros(numNuclei,1)));      %VMHvl                                                                                
                            
Tdt = struct(   ...    %structure array containing mean TDT intensity data for each cellular ROI found in the image        
                                'total',uint16(zeros(numNuclei,1)),...          %total VMH
                                'dm',uint16(zeros(numNuclei,1)),...   %VMHdm
                                'vl',uint16(zeros(numNuclei,1)));      %VMHvl     

hWait = waitbar(0,'Finding cells...');  %create waitbar because this is the part that takes a while

% NOTE: start at index 2 since the region corresponding to '1' is the dead space between nuclei, and end at numNuclei-1 since last one seems to be zero:                            
for i = 2:numNuclei-1
    currentIndeces = find(Ibwater==i);
    Gfp.total(i-1) = mean(Ig(currentIndeces));  %extract mean pixels intensity in region in green channel
    Tdt.total(i-1) = mean(Ir(currentIndeces));  %extract mean pixels intensity in region in red channel
    
    % for subregion counts, set values outside the mask (i.e. entire cell needs to be inside ROI to count) to 0 for later
    % parsing; to do this I ask below if the length of the number of indices found for the original watershed matrix is the same as the
    % length of the number of indices found if we 'AND' with the VMHdm mask, for which values outside the mask are zero:
    insideDmRoi = length(currentIndeces)==length(find(Ibwater==i & roiPatchArray(1,2).patchMask)==i);  %second ROI mask in for DM
    if insideDmRoi
        Gfp.dm(i-1) = Gfp.total(i-1);
        Tdt.dm(i-1) = Tdt.total(i-1);
    else
        Gfp.dm(i-1) = 0;  %if not in ROI, set to 0 to parse out in later analysis
        Tdt.dm(i-1) = 0;
    end
    
    insideVLRoi = length(currentIndeces)==length(find(Ibwater==i & roiPatchArray(1,3).patchMask)==i); %third ROI mask in for VL   
    if insideVLRoi
        Gfp.vl(i-1) = Gfp.total(i-1);
        Tdt.vl(i-1) = Tdt.total(i-1);
    else
        Gfp.vl(i-1) = 0;
        Tdt.vl(i-1) = 0;
    end
    
    waitbar(double(i)/double(numNuclei), hWait);
end
delete(hWait); 

end %end of function

