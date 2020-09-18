% top level script to ccompare ROIs drawn by 2 different people
clear all; close all;

filepathNd2 = 'C:\data\jason\peterVirusData\';
% filepathNd2 = 'C:\Users\keller\School\UCSD\stowersLab\code\cellSegmentation\';
mouseNums = {'p1', 'j1'};  %ex. 2nd entry is mouseNums{1,2}
totalMice = size(mouseNums, 2);
lastSections = {4, 4, 4, 4}; % to keep track of how many sections were imaged for each mouse

figure;

for k = 1:totalMice
    mouseNum = mouseNums{1,k};
    lastSection = lastSections{1,k};
    for j = 1:lastSection %lastSections{1,1};
        filepathTotal = [filepathNd2, mouseNum, '\', mouseNum, '_', num2str(j), '.mat']  %display full path
        % do all manual steps for each section first (i.e. ROIs and thresholds):
        loadStr = ['load(', sprintf( '\''' ), filepathTotal, sprintf( '\''' ), ')'];
        eval(loadStr);
        
        roiPatchArray_comp{j*k} = roiPatchArray; %#ok<*SAGROW>
        vmhXLim_comp{j*k} = vmhXLim;
        vmhYLim_comp{j*k} = vmhYLim;
        GfpThresh_comp{j*k} = GfpThresh;
%         GfpThreshImage
        TdtThresh_comp{j*k} = TdtThresh;
%         TdtThreshImage
        Ig_comp{j*k} = Ig;
        Ir_comp{j*k} = Ir;
        Ibwater_comp{j*k} = Ibwater;
        
        vmhLimAltMax = max(roiPatchArray(1, 1).patchPosition);
        vmhLimAltMin = min(roiPatchArray(1, 1).patchPosition);
        xMinAlt = uint16(vmhLimAltMin(1));
        xMaxAlt = uint16(vmhLimAltMax(1));
        yMinAlt = uint16(vmhLimAltMin(2));
        yMaxAlt = uint16(vmhLimAltMax(2));
        subimageIbwater = Ibwater(yMinAlt:yMaxAlt, xMinAlt:xMaxAlt);
        subimageIbwater(subimageIbwater<=1) = 0;  %watershed 1 is backround, so discard
        subimageIbwater(subimageIbwater>1) = 1;
        
        subimage1(:,:,j*k) = %original nuclear thresholded image
        subimage2(:,:,j*k) = 
        subimage3(:,:,j*k) = 
        subimage4(:,:,j*k) = 
        
        
    end

end