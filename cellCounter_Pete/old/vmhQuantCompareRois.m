% top level script to ccompare ROIs drawn by 2 different people
clear all; close all;

filepathNd2 = 'C:\data\jason\peterVirusData\';
% filepathNd2 = 'C:\Users\keller\School\UCSD\stowersLab\code\cellSegmentation\';
mouseNums = {'p1', 'j1'};  %ex. 2nd entry is mouseNums{1,2}
totalMice = size(mouseNums, 2);
lastSection = 4; % to keep track of how many sections were imaged for each mouse

for j = 1:lastSection 
    mouseNum = mouseNums{1,1};
    %first load original image to compare:
    filepathTotal = [filepathNd2, mouseNum, '\', num2str(j), '.nd2']

    zplanes = 1; tframes = 1; channel = 1;  %only one channel can be imported at once
    [vol] = imreadBF(filepathTotal,zplanes,tframes,channel); %read 1st channel in ND2 file - Hoechst
    Ibz = uint16(vol);;
%     channel = 2;
%     [vol] = imreadBF(filepathTotal,zplanes,tframes,channel); %read 2nd channel in ND2 file - GFP
%     Igz = uint16(vol);
%     % meanGfp = mean(mean(Ig)); % perhaps use total image mean as threshold later
%     channel = 3;
%     [vol] = imreadBF(filepathTotal,zplanes,tframes,channel); %read 3rd channel in ND2 file - TDT
%     Irz = uint16(vol);
    % meanTdt = mean(mean(Ir)); % perhaps use total image mean as threshold later
    clear vol; %free up memory
    Irgb = zeros(size(Ibz,1),size(Ibz,2),3); %set R & G channels to zero
    Irgb(:,:,3) = double(Ibz)./ double(max(max(Ibz)));
    
    filepathTotal = [filepathNd2, mouseNum, '\', mouseNum, '_', num2str(j), '.mat']  %display full path
    % do all manual steps for each section first (i.e. ROIs and thresholds):
    loadStr = ['load(', sprintf( '\''' ), filepathTotal, sprintf( '\''' ), ')'];
    eval(loadStr);
    roiPatchArray_1 = roiPatchArray; %#ok<*SAGROW>
    vmhXLim_1 = vmhXLim;
    vmhYLim_1 = vmhYLim;
    GfpThresh_1 = GfpThresh;
%         GfpThreshImage
    TdtThresh_1 = TdtThresh;
%         TdtThreshImage
    Ig_1 = Ig;
    Ir_1 = Ir;
    Ibwater_1 = Ibwater;
    
    mouseNum = mouseNums{1,2};
    filepathTotal = [filepathNd2, mouseNum, '\', mouseNum, '_', num2str(j), '.mat']  %display full path
    % do all manual steps for each section first (i.e. ROIs and thresholds):
    loadStr = ['load(', sprintf( '\''' ), filepathTotal, sprintf( '\''' ), ')'];
    eval(loadStr);
    roiPatchArray_2 = roiPatchArray; %#ok<*SAGROW>
    vmhXLim_2 = vmhXLim;
    vmhYLim_2 = vmhYLim;
    GfpThresh_2 = GfpThresh;
%         GfpThreshImage
    TdtThresh_2 = TdtThresh;
%         TdtThreshImage
    Ig_2 = Ig;
    Ir_2 = Ir;
    Ibwater_2 = Ibwater;  
    
    pixelOffset = 300;
    vmhLimAltMax = max(roiPatchArray_1(1, 1).patchPosition);
    vmhLimAltMin = min(roiPatchArray_1(1, 1).patchPosition);
    xMinAlt = uint16(vmhLimAltMin(1)-pixelOffset);
    xMaxAlt = uint16(vmhLimAltMax(1)+pixelOffset);
    yMinAlt = uint16(vmhLimAltMin(2)-pixelOffset);
    yMaxAlt = uint16(vmhLimAltMax(2)+pixelOffset);
    
    subimage1 = Irgb(yMinAlt:yMaxAlt, xMinAlt:xMaxAlt, :); %original nuclear thresholded image
    
    subimage2 = zeros(size(subimage1,1),size(subimage1,2),3); %total VMH overlay
    subimage2(:,:,1) = roiPatchArray_1(1,1).patchMask(yMinAlt:yMaxAlt, xMinAlt:xMaxAlt); %red is first person's ROI
    subimage2(:,:,2) = roiPatchArray_2(1,1).patchMask(yMinAlt:yMaxAlt, xMinAlt:xMaxAlt); %green is first person's ROI
    
    subimage3 = zeros(size(subimage1,1),size(subimage1,2),3);%VMH DM overlay
    subimage3(:,:,1) = ((roiPatchArray_1(1,1).patchMask(yMinAlt:yMaxAlt, xMinAlt:xMaxAlt)) & (roiPatchArray_1(1,2).patchMask(yMinAlt:yMaxAlt, xMinAlt:xMaxAlt))); %red is first person's ROI
    subimage3(:,:,2) = ((roiPatchArray_2(1,1).patchMask(yMinAlt:yMaxAlt, xMinAlt:xMaxAlt)) & (roiPatchArray_2(1,2).patchMask(yMinAlt:yMaxAlt, xMinAlt:xMaxAlt))); %green is first person's ROI
    
    subimage4 = zeros(size(subimage1,1),size(subimage1,2),3);%VMH VL overlay
    subimage4(:,:,1) = ((roiPatchArray_1(1,1).patchMask(yMinAlt:yMaxAlt, xMinAlt:xMaxAlt)) & (~roiPatchArray_1(1,2).patchMask(yMinAlt:yMaxAlt, xMinAlt:xMaxAlt))); %red is first person's ROI
    subimage4(:,:,2) = ((roiPatchArray_2(1,1).patchMask(yMinAlt:yMaxAlt, xMinAlt:xMaxAlt)) & (~roiPatchArray_2(1,2).patchMask(yMinAlt:yMaxAlt, xMinAlt:xMaxAlt))); %green is first person's ROI
    
    hFig = figure;
    subplot(2,2,1); imshow(subimage1);
    subplot(2,2,2); imshow(subimage2);
    subplot(2,2,3); imshow(subimage3);
    subplot(2,2,4); imshow(subimage4);
    set(gcf, 'Position', get(0,'Screensize')); % Maximize figure
    filenameFig = [filepathNd2, mouseNum, '\', num2str(j), '.eps'];
    print(hFig, '-depsc', filenameFig);
end



