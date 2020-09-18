% top level script to count all VMH sections and plot results
clear all; close all;

gfpThresh = 1;  %this is the threshold of mean intensity above which we consider a positive GFP cell
tdtThresh = 1;  %this is the threshold of mean intensity above which we consider a positive TDT cell

filepathNd2 = 'C:\data\jason\peterVirusData\';
mouseNum = '19249';
sectionNum = 1:10;
for j = 1:length(sectionNum)
    filepathTotal = [filepathNd2, mouseNum, '\', num2str(sectionNum(j)), '.nd2']  %display full path
    [roiPatchArray, Gfp, Tdt] = vmhQuant(filepathNd2, mouseNum, num2str(sectionNum(j)));
    %save data:
    matName = [filepathNd2, mouseNum, '\', mouseNum, '_', num2str(sectionNum(j)), '.mat'];
    save(matName, 'roiPatchArray', 'Gfp', 'Tdt');
    %compute percentages of cell in different fluo channels & regions, and build of array to compute statistics later:
    [percentGfpCellsOfTotal(j) , percentTdtOfGfpTotal(j), percentTdtOfGfpDm(j), percentTdtOfGfpVl(j)]...
        = vmhCountCells(Gfp, Tdt, gfpThresh, tdtThresh); %#ok<*SAGROW>
    close all; %prevent large number figures. and also clear for next wait function
end
matName = [filepathNd2, mouseNum, '\', mouseNum, '_all.mat'];
save(matName, '-regexp', '^percent');  %save all variables starting with 'percent'

% mouseNum = '19485';
% sectionNum = 1:9;
% 
% mouseNum = '70321';
% sectionNum = 1:9;
% 
% mouseNum = '70324';
% sectionNum = 1:8;
% 
% mouseNum = '70750';
% sectionNum = 1:7;
% 
% mouseNum = '70753';
% sectionNum = 1:8;
% 
% mouseNum = '70755';
% sectionNum = 1:8;
% 
% mouseNum = '70756';
% sectionNum = 1:9;


% total 8 mice, 68 sections (~8 sections/mouseVMH)