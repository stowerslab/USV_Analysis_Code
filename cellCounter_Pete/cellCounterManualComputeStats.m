% top level script to count all sections and plot results
clear all; close all;

filepathNd2 = 'C:\data\Jason\microscope\2016_01_CrhTdtwithEsrImmuno\male\';

%%% INPUT:
mouseNumbers = {'esrF1' 'esrF2'};
lastSects =  {8, 5}; % to keep track of how many sections were imaged for each mouse
firstSects = {1, 1};

totalMice = size(mouseNumbers, 2);

totalRed = [];
totalGreen = [];
totalBlue = [];
totalOverlap = [];

for k = 1:totalMice
    mouseNumber = mouseNumbers{1,k};
    lastSect = lastSects{1,k};
    firstSect = firstSects{1,k};
    for j = firstSect:lastSect
        matName = [filepathNd2, mouseNumber, '\', mouseNumber, '_', num2str(j), '.mat']; %use saved data already counted
        load(matName);
        totalRed(end+1) = numRed; %#ok<*SAGROW>
        totalGreen(end+1) = numGreen;
        totalBlue(end+1) = numBlue;
        totalOverlap(end+1) = numOverlap;
    end    
end

