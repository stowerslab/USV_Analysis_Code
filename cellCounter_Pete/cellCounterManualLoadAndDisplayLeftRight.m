% top level script to count all sections and plot results
clear all; close all;
filepathNd2 = 'C:\data\Jason\cellCounting\wtDREADD\';

%%% INPUT:
mouseNums = {'h7'};
sections =  {'1R'}; % '2R' '3R' '1L' '2L' '3L'}; 
numSections = size(sectionsL, 2);
totalMice = size(mouseNums, 2);

redTotal = 0;
greenTotal = 0;
overlapTotal = 0;

for k = 1:totalMice
    mouseNum = mouseNums{1,k};
    for j = 1:numSections
        
        thisSectionL = sectionsL{1,j};
        matName = [filepathNd2, mouseNum, '\', thisSectionL, '.mat']; %use saved data already counted
        load(matName); %careful not to save/load loop variables!!!
        redTotal = numRed;
        greenTotal = numGreen;
        overlapTotal = numOverlap;
        
        thisSectionR = sectionsR{1,j};
        matName = [filepathNd2, mouseNum, '\', thisSectionR, '.mat']; %use saved data already counted
        load(matName); %careful not to save/load loop variables!!!
        redTotal = redTotal + numRed;
        greenTotal = greenTotal + numGreen;
        overlapTotal = overlapTotal+ numOverlap;
        
        display(['counting ', mouseNum, ', section '  thisSectionL]);
        display(['numRed = ', num2str(redTotal), '; numGreen = ', num2str(greenTotal), '; numOverlap = ', num2str(overlapTotal)]);
         
    end    
end

