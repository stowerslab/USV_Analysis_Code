% top level script to count all sections and plot results
clear all; close all;
% filepathNd2 = 'C:\data\Jason\cellCounting\crhDREADD\';
% filepathNd2 = 'C:\data\Jason\cellCounting\esrDREADD\';
% filepathNd2 = 'C:\data\Jason\cellCounting\wtDREADD\';
% filepathNd2 = 'C:\data\Jason\cellCounting\crhChR2\';
filepathNd2 = 'C:\data\Jason\cellCounting\esrChR2\';
% filepathNd2 = 'C:\data\Jason\cellCounting\controlChR2\';
% filepathNd2 = 'C:\data\Jason\cellCounting\crhArchT\';
% filepathNd2 = 'C:\data\Jason\cellCounting\esrArchT\';


%%% INPUT:
% mouseNums = {'crh1' 'crh4' 'crh7' 'crh9' 'crh15' 'crh16' 'mc3' 'mc6'}; %Crh DREADD
% mouseNums = {'i1' 'i3' 'k3' 'm1' 'm4' 'n1' 'n2' 'm3'}; %Esr DREADD
% mouseNums = {'h7' 'h8' 'h9' 'h11' 'h12'}; %WT DREADD
% mouseNums = {'cc2\red' 'cc3\red' 'cc6\red' 'cc7\red' 'cc9\red' 'cc16\red'}; %Crh ChR2, red only
% mouseNums = {'cc2' 'cc3' 'cc6' 'cc7' 'cc9' 'cc12' 'cc13' 'cc16'}; %Crh ChR2, green only
mouseNums = {'ec4' 'ec5' 'ec7' 'ec8' 'ec10' 'ec11' 'ec12' 'n7' 'n15'}; %Esr ChR2, red only
% mouseNums = {'econ1' 'econ2' 'econ3'}; % control ChR2
% mouseNums = {'a2' 'a4' 'ac9' 'ac10' }; % Crh ArchT
% mouseNums = {'a3' 'ae6' 'ae8'}; % Esr ArchT


totalMice = size(mouseNums, 2);


for k = 1:totalMice
    mouseNum = mouseNums{1,k};
    matName = [filepathNd2, mouseNum, '\', mouseNum, '_totalCounts.mat']; %
%     matName = [filepathNd2, mouseNum, '\', mouseNums{1,k}(1:end-4), '_totalCountsRed.mat']; %for Crh ChR2 red only
    load(matName); %careful not to save/load loop variables!!!
    redTotalVec(k) = redTotal;
%     redOutsideVec(k) = redOutside;
    greenTotalVec(k) = greenTotal; %#ok<*SAGROW>
    greenOutsideVec(k) = greenOutside;
    greenVLVec(k) = greenVL;
    greenDMVec(k) = greenDM;
end

meanGreen = mean(greenTotalVec)
stdGreen = std(greenTotalVec)
meanRed = mean(redTotalVec)
stdRed = std(redTotalVec)