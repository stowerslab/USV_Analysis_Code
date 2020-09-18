% top level script to count all VMH sections and plot results
clear all; close all;

filepathNd2 = 'C:\data\jason\peterVirusData\';
% filepathNd2 = 'C:\Users\keller\School\UCSD\stowersLab\code\cellSegmentation\';

% mouseNum = '19249';
% sectionNum = 1:10;
% for j = 1:length(sectionNum)
%     filepathTotal = [filepathNd2, mouseNum, '\', num2str(sectionNum(j)), '.nd2']  %display full path
%     [roiPatchArray, gfp, tdt, gfpThresh, tdtThresh] = vmhQuant(filepathNd2, mouseNum, num2str(sectionNum(j)));
%     %save data for each individual section:
%     matName = [filepathNd2, mouseNum, '\', mouseNum, '_', num2str(sectionNum(j)), '.mat'];
%     save(matName, 'roiPatchArray', 'gfp', 'tdt', 'gfpThresh', 'tdtThresh');
%     %compute percentages of cell in different fluo channels & regions, and build of array to compute statistics later:
%     [numCellsTotal(j), numGfpCellsTotal(j), numGfpCellsInDm(j), numGfpCellsInVl(j), numTdtCellsTotal(j), numTdtCellsInDm(j), numTdtCellsInVl(j), percentGfpCellsOfTotal(j) , percentTdtOfGfpTotal(j), percentTdtOfGfpDm(j), percentTdtOfGfpVl(j)]...
%         = vmhCountCells(gfp, tdt, gfpThresh, tdtThresh); %#ok<*SAGROW>
%     close all; %prevent large number figures. and also clear for next wait function
% end
% %save summary data for entire mouse:
% matName = [filepathNd2, mouseNum, '\', mouseNum, '_all.mat'];
% save(matName, '-regexp', '^percent', '^num');  %save all variables starting with 'percent', etc.

% mouseNum = '19485';
% sectionNum = 1:9;
% for j = 1:length(sectionNum)
%     filepathTotal = [filepathNd2, mouseNum, '\', num2str(sectionNum(j)), '.nd2']  %display full path
%     [roiPatchArray, gfp, tdt, gfpThresh, tdtThresh] = vmhQuant(filepathNd2, mouseNum, num2str(sectionNum(j)));
%     %save data for each individual section:
%     matName = [filepathNd2, mouseNum, '\', mouseNum, '_', num2str(sectionNum(j)), '.mat'];
%     save(matName, 'roiPatchArray', 'gfp', 'tdt', 'gfpThresh', 'tdtThresh');
%     %compute percentages of cell in different fluo channels & regions, and build of array to compute statistics later:
%     [numCellsTotal(j), numGfpCellsTotal(j), numGfpCellsInDm(j), numGfpCellsInVl(j), numTdtCellsTotal(j), numTdtCellsInDm(j), numTdtCellsInVl(j), percentGfpCellsOfTotal(j) , percentTdtOfGfpTotal(j), percentTdtOfGfpDm(j), percentTdtOfGfpVl(j)]...
%         = vmhCountCells(gfp, tdt, gfpThresh, tdtThresh); %#ok<*SAGROW>
%     close all; %prevent large number figures. and also clear for next wait function
% end
% %save summary data for entire mouse:
% matName = [filepathNd2, mouseNum, '\', mouseNum, '_all.mat'];
% save(matName, '-regexp', '^percent', '^num');  %save all variables starting with 'percent', etc.

% mouseNum = '70321';
% sectionNum = 1:9;
% for j = 1:length(sectionNum)
%     filepathTotal = [filepathNd2, mouseNum, '\', num2str(sectionNum(j)), '.nd2']  %display full path
%     [roiPatchArray, gfp, tdt, gfpThresh, tdtThresh] = vmhQuant(filepathNd2, mouseNum, num2str(sectionNum(j)));
%     %save data for each individual section:
%     matName = [filepathNd2, mouseNum, '\', mouseNum, '_', num2str(sectionNum(j)), '.mat'];
%     save(matName, 'roiPatchArray', 'gfp', 'tdt', 'gfpThresh', 'tdtThresh');
%     %compute percentages of cell in different fluo channels & regions, and build of array to compute statistics later:
%     [numCellsTotal(j), numGfpCellsTotal(j), numGfpCellsInDm(j), numGfpCellsInVl(j), numTdtCellsTotal(j), numTdtCellsInDm(j), numTdtCellsInVl(j), percentGfpCellsOfTotal(j) , percentTdtOfGfpTotal(j), percentTdtOfGfpDm(j), percentTdtOfGfpVl(j)]...
%         = vmhCountCells(gfp, tdt, gfpThresh, tdtThresh); %#ok<*SAGROW>
%     close all; %prevent large number figures. and also clear for next wait function
% end
% %save summary data for entire mouse:
% matName = [filepathNd2, mouseNum, '\', mouseNum, '_all.mat'];
% save(matName, '-regexp', '^percent', '^num');  %save all variables starting with 'percent', etc.

% mouseNum = '70324';
% sectionNum = 1:7;
% for j = 1:length(sectionNum)
%     filepathTotal = [filepathNd2, mouseNum, '\', num2str(sectionNum(j)), '.nd2']  %display full path
%     [roiPatchArray, gfp, tdt, gfpThresh, tdtThresh] = vmhQuant(filepathNd2, mouseNum, num2str(sectionNum(j)));
%     %save data for each individual section:
%     matName = [filepathNd2, mouseNum, '\', mouseNum, '_', num2str(sectionNum(j)), '.mat'];
%     save(matName, 'roiPatchArray', 'gfp', 'tdt', 'gfpThresh', 'tdtThresh');
%     %compute percentages of cell in different fluo channels & regions, and build of array to compute statistics later:
%     [numCellsTotal(j), numGfpCellsTotal(j), numGfpCellsInDm(j), numGfpCellsInVl(j), numTdtCellsTotal(j), numTdtCellsInDm(j), numTdtCellsInVl(j), percentGfpCellsOfTotal(j) , percentTdtOfGfpTotal(j), percentTdtOfGfpDm(j), percentTdtOfGfpVl(j)]...
%         = vmhCountCells(gfp, tdt, gfpThresh, tdtThresh); %#ok<*SAGROW>
%     close all; %prevent large number figures. and also clear for next wait function
% end
% %save summary data for entire mouse:
% matName = [filepathNd2, mouseNum, '\', mouseNum, '_all.mat'];
% save(matName, '-regexp', '^percent', '^num');  %save all variables starting with 'percent', etc.

% mouseNum = '70750';
% sectionNum = 1:7;
% for j = 1:length(sectionNum)
%     filepathTotal = [filepathNd2, mouseNum, '\', num2str(sectionNum(j)), '.nd2']  %display full path
%     [roiPatchArray, gfp, tdt, gfpThresh, tdtThresh] = vmhQuant(filepathNd2, mouseNum, num2str(sectionNum(j)));
%     %save data for each individual section:
%     matName = [filepathNd2, mouseNum, '\', mouseNum, '_', num2str(sectionNum(j)), '.mat'];
%     save(matName, 'roiPatchArray', 'gfp', 'tdt', 'gfpThresh', 'tdtThresh');
%     %compute percentages of cell in different fluo channels & regions, and build of array to compute statistics later:
%     [numCellsTotal(j), numGfpCellsTotal(j), numGfpCellsInDm(j), numGfpCellsInVl(j), numTdtCellsTotal(j), numTdtCellsInDm(j), numTdtCellsInVl(j), percentGfpCellsOfTotal(j) , percentTdtOfGfpTotal(j), percentTdtOfGfpDm(j), percentTdtOfGfpVl(j)]...
%         = vmhCountCells(gfp, tdt, gfpThresh, tdtThresh); %#ok<*SAGROW>
%     close all; %prevent large number figures. and also clear for next wait function
% end
% %save summary data for entire mouse:
% matName = [filepathNd2, mouseNum, '\', mouseNum, '_all.mat'];
% save(matName, '-regexp', '^percent', '^num');  %save all variables starting with 'percent', etc.

% mouseNum = '70753';
% sectionNum = 1:8;
% for j = 1:length(sectionNum)
%     filepathTotal = [filepathNd2, mouseNum, '\', num2str(sectionNum(j)), '.nd2']  %display full path
%     [roiPatchArray, gfp, tdt, gfpThresh, tdtThresh] = vmhQuant(filepathNd2, mouseNum, num2str(sectionNum(j)));
%     %save data for each individual section:
%     matName = [filepathNd2, mouseNum, '\', mouseNum, '_', num2str(sectionNum(j)), '.mat'];
%     save(matName, 'roiPatchArray', 'gfp', 'tdt', 'gfpThresh', 'tdtThresh');
%     %compute percentages of cell in different fluo channels & regions, and build of array to compute statistics later:
%     [numCellsTotal(j), numGfpCellsTotal(j), numGfpCellsInDm(j), numGfpCellsInVl(j), numTdtCellsTotal(j), numTdtCellsInDm(j), numTdtCellsInVl(j), percentGfpCellsOfTotal(j) , percentTdtOfGfpTotal(j), percentTdtOfGfpDm(j), percentTdtOfGfpVl(j)]...
%         = vmhCountCells(gfp, tdt, gfpThresh, tdtThresh); %#ok<*SAGROW>
%     close all; %prevent large number figures. and also clear for next wait function
% end
% %save summary data for entire mouse:
% matName = [filepathNd2, mouseNum, '\', mouseNum, '_all.mat'];
% save(matName, '-regexp', '^percent', '^num');  %save all variables starting with 'percent', etc.

% mouseNum = '70755';
% sectionNum = 1:8;
% for j = 1:length(sectionNum)
%     filepathTotal = [filepathNd2, mouseNum, '\', num2str(sectionNum(j)), '.nd2']  %display full path
%     [roiPatchArray, gfp, tdt, gfpThresh, tdtThresh] = vmhQuant(filepathNd2, mouseNum, num2str(sectionNum(j)));
%     %save data for each individual section:
%     matName = [filepathNd2, mouseNum, '\', mouseNum, '_', num2str(sectionNum(j)), '.mat'];
%     save(matName, 'roiPatchArray', 'gfp', 'tdt', 'gfpThresh', 'tdtThresh');
%     %compute percentages of cell in different fluo channels & regions, and build of array to compute statistics later:
%     [numCellsTotal(j), numGfpCellsTotal(j), numGfpCellsInDm(j), numGfpCellsInVl(j), numTdtCellsTotal(j), numTdtCellsInDm(j), numTdtCellsInVl(j), percentGfpCellsOfTotal(j) , percentTdtOfGfpTotal(j), percentTdtOfGfpDm(j), percentTdtOfGfpVl(j)]...
%         = vmhCountCells(gfp, tdt, gfpThresh, tdtThresh); %#ok<*SAGROW>
%     close all; %prevent large number figures. and also clear for next wait function
% end
% %save summary data for entire mouse:
% matName = [filepathNd2, mouseNum, '\', mouseNum, '_all.mat'];
% save(matName, '-regexp', '^percent', '^num');  %save all variables starting with 'percent', etc.
 
mouseNum = '70756';
sectionNum = 1:9;
for j = 1:length(sectionNum)
    filepathTotal = [filepathNd2, mouseNum, '\', num2str(sectionNum(j)), '.nd2']  %display full path
    [roiPatchArray, gfp, tdt, gfpThresh, tdtThresh] = vmhQuant(filepathNd2, mouseNum, num2str(sectionNum(j)));
    %save data for each individual section:
    matName = [filepathNd2, mouseNum, '\', mouseNum, '_', num2str(sectionNum(j)), '.mat'];
    save(matName, 'roiPatchArray', 'gfp', 'tdt', 'gfpThresh', 'tdtThresh');
    %compute percentages of cell in different fluo channels & regions, and build of array to compute statistics later:
    [numCellsTotal(j), numGfpCellsTotal(j), numGfpCellsInDm(j), numGfpCellsInVl(j), numTdtCellsTotal(j), numTdtCellsInDm(j), numTdtCellsInVl(j), percentGfpCellsOfTotal(j) , percentTdtOfGfpTotal(j), percentTdtOfGfpDm(j), percentTdtOfGfpVl(j)]...
        = vmhCountCells(gfp, tdt, gfpThresh, tdtThresh); %#ok<*SAGROW>
    close all; %prevent large number figures. and also clear for next wait function
end
%save summary data for entire mouse:
matName = [filepathNd2, mouseNum, '\', mouseNum, '_all.mat'];
save(matName, '-regexp', '^percent', '^num');  %save all variables starting with 'percent', etc.

% total 8 mice, 68 sections (~8 sections/mouseVMH)