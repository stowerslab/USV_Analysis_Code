% top level script to count all VMH sections and plot results
clear all; close all;

filepathNd2 = 'C:\data\Jason\microscope\2016_01_CrhTdtwithEsrImmuno\female\crh_EsrImm_perf\';
viewOverlay = 0;

%%% INPUT:
mouseNums = {'esrF1'};  %ex. 2nd entry is mouseNums{1,2}
lastSections =  {8}; % to keep track of how many sections were imaged for each mouse
firstSections = {1};

totalMice = size(mouseNums, 2);

%% do all manual steps for each section first (i.e. ROIs and thresholds):
for k = 1:totalMice
    mouseNum = mouseNums{1,k};
    lastSection = lastSections{1,k};
    firstSection = firstSections{1,k};
    for j = firstSection:lastSection
        filepathTotal = [filepathNd2, mouseNum, '\', num2str(j), '.nd2']  %#ok<NOPTS> %display full path
        [roiPatchArray, vmhXLim, vmhYLim, GfpThresh, GfpThreshImage, TdtThresh, TdtThreshImage, Ig, Ir, Ibwater] = vmhSetRoiThresh(filepathNd2, mouseNum, num2str(j)); % do all manual steps for each section first (i.e. ROIs and thresholds)
        matName = [filepathNd2, mouseNum, '\', mouseNum, '_', num2str(j), '.mat']; %save data for each individual section:
        save(matName, 'roiPatchArray', 'vmhXLim', 'vmhYLim', 'GfpThresh', 'GfpThreshImage', 'TdtThresh', 'TdtThreshImage', 'Ig', 'Ir', 'Ibwater');
        close all; %prevent large number figures. and also clear for next wait function
    end
end

%% Now do all of the post-manual steps together and pool data:
% for k = 1:totalMice
%     mouseNum = mouseNums{1,k};
%     lastSection = lastSections{1,k};
%     firstSection = firstSections{1,k};
%     for j = firstSection:lastSection
%         filepathTotal = [filepathNd2, mouseNum, '\', mouseNum, '_', num2str(j), '.mat']  %display full path
%         load(filepathTotal);
% 
%         % uncomment below to create and view RGB overlay of VMH cellular ROIs & thresholded R & G channels, to test:
%         if viewOverlay
%             vmhLimAltMax = max(roiPatchArray(1, 1).patchPosition);
%             vmhLimAltMin = min(roiPatchArray(1, 1).patchPosition);
%             xMinAlt = uint16(vmhLimAltMin(1));
%             xMaxAlt = uint16(vmhLimAltMax(1));
%             yMinAlt = uint16(vmhLimAltMin(2));
%             yMaxAlt = uint16(vmhLimAltMax(2));
%             subimageIbwater = Ibwater(yMinAlt:yMaxAlt, xMinAlt:xMaxAlt);
%             subimageIbwater(subimageIbwater<=1) = 0;  %watershed 1 is backround, so discard
%             subimageIbwater(subimageIbwater>1) = 1;
%             Irgb = zeros(size(subimageIbwater,1),size(subimageIbwater,2),3);
%             Irgb(:,:,1) = TdtThreshImage; 
%             Irgb(:,:,2) = GfpThreshImage;
%             Irgb(:,:,3) = subimageIbwater;
%             figure; imshow(Irgb); pause on; pause(1); pause off;
%             imFilename = [filepathNd2, mouseNum, '\', mouseNum, '_', num2str(j), '.jpg'] ;
%             imwrite(Irgb, imFilename, 'jpg');
%         end
% 
%         [Gfp, Tdt, numCells] = vmhQuant(Ig, Ir, Ibwater, roiPatchArray); 
%         %compute percentages of cell in different fluo channels & regions, and build of array to compute statistics later:
%         [numCellsTotal(j), numGfpCellsTotal(j), numGfpCellsInDm(j), numGfpCellsInVl(j), numTdtCellsTotal(j), numTdtCellsInDm(j), numTdtCellsInVl(j), percentGfpCellsOfTotal(j) , percentTdtOfGfpTotal(j), percentTdtOfGfpDm(j), percentTdtOfGfpVl(j)]...
%             = vmhCountCells(Gfp, Tdt, numCells, GfpThresh, TdtThresh); %#ok<*SAGROW>
%     end
%     %save summary data for entire mouse:
%     matName = [filepathNd2, mouseNum, '\', mouseNum, '_all.mat'];
%     save(matName, '-regexp', '^percent', '^num');  %save all variables starting with 'percent', etc.
% end

