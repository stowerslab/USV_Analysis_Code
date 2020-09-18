% top level script to count all sections and plot results
clear all; close all;
filepathNd2 = 'C:\data\Jason\microscope\2016_01_CrhTdtwithEsrImmuno\male\';

%%% INPUT:
mouseNums = {'crhTdtMale5wHoechst\count' 'crhTdtMale6wHoechst\count' 'crhTdtMale7wHoechst\count' 'crhTdtMale8wHoechst\count' 'crhTdtMale9wHoechst\count' 'crhTdtMale10wHoechst\count'};
sectionsL =  {'1L' '2L' '3L' '4L' '5L' '6L'};
sectionsR = {'1R' '2R' '3R' '4R' '5R' '6R'}; 
numSections = size(sectionsL, 2);
totalMice = size(mouseNums, 2);

redTotal = 0;
greenTotal = 0;
overlapTotal = 0;

totalNisslCount = zeros(totalMice,1);

for k = 1:totalMice
    mouseNum = mouseNums{1,k};
    for j = 1:numSections
        
        thisSectionL = sectionsL{1,j};
        matName = [filepathNd2, mouseNum, '\', thisSectionL, '_Nissl.mat']; %use saved data already counted
        display(['counting ', mouseNum, ', section '  thisSectionL]);
        maxImageBitValue = 2^16-1; %for ND2 format
        maxPixelValue = 2^12-1;
        [Idapi, Iesr, Icrh, Inissl] = readNd2_4ch(filepathNd2, mouseNum, thisSectionL);
        
        Iroi = zeros(size(Icrh,1),size(Icrh,2),3);
        Iroi(:,:,1) = double(Icrh)./ maxPixelValue; 
        Iroi(:,:,2) = double(Iesr)./ maxPixelValue;
        Iroi(:,:,3) = double(Inissl)./ maxPixelValue;
%         imshow(Iroi);

        [roiPatchArray, xMin, xMax, yMin, yMax] = sectionSetRoiRgb(Iroi); % select ROI
        IroiMask = roiPatchArray.patchMask(yMin:yMax,xMin:xMax);
        IsubNissl = Iroi(yMin:yMax,xMin:xMax,3); %reduce & mask
        IsubNissl(~IroiMask) = 0;
        
        [Ioverlay, ~, nisslCC] = sectionSegmentCountRoi(IsubNissl);   
        totalNisslCount(k) = totalNisslCount(k) + nisslCC.NumObjects;
        save(matName, 'totalNisslCount', 'roiPatchArray', 'xMin', 'xMax', 'yMin', 'yMax', 'Ioverlay')
        
        thisSectionR = sectionsR{1,j};
        matName = [filepathNd2, mouseNum, '\', thisSectionR, '_Nissl.mat']; %use saved data already counted
        display(['counting ', mouseNum, ', section '  thisSectionR]);
        [Idapi, Iesr, Icrh, Inissl] = readNd2_4ch(filepathNd2, mouseNum, thisSectionR);
        
        Iroi = zeros(size(Icrh,1),size(Icrh,2),3);
        Iroi(:,:,1) = double(Icrh)./ maxPixelValue; 
        Iroi(:,:,2) = double(Iesr)./ maxPixelValue;
        Iroi(:,:,3) = double(Inissl)./ maxPixelValue;
%         imshow(Iroi);

        [roiPatchArray, xMin, xMax, yMin, yMax] = sectionSetRoiRgb(Iroi); % select ROI
        IroiMask = roiPatchArray.patchMask(yMin:yMax,xMin:xMax);
        IsubNissl = Iroi(yMin:yMax,xMin:xMax,3); %reduce & mask
        IsubNissl(~IroiMask) = 0;
        
        [Ioverlay, ~, nisslCC] = sectionSegmentCountRoi(IsubNissl);   
        totalNisslCount(k) = totalNisslCount(k) + nisslCC.NumObjects;
        save(matName, 'totalNisslCount', 'roiPatchArray', 'xMin', 'xMax', 'yMin', 'yMax', 'Ioverlay')

    end    
end

