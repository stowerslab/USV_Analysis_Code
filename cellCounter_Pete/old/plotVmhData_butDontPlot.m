%script to read in VMH data and plot
clear all; close all;

mice = ['p1'; 'j1'; 'p2'; 'j2'];
aggr = []; %indeces of aggression mice

stat1 = 'ratioTdtCellsVlToDm';
stat2 = 'ratioPercentTdtOfGfpVlToDm';
stat3 = 'percentTdtOfGfpVl';
stat4 = 'percentTdtOfGfpDm';
stat5 = 'numTdtCellsInDm';
stat6 = 'numTdtCellsInVl';
stat7 = 'numGfpCellsInDm';
stat8 = 'numGfpCellsInVl';


aggrStats = struct(stat1,[], stat2,[], stat3,[], stat4,[], stat5,[], stat6,[], stat7,[], stat8,[]);
contStats = struct(stat1,[], stat2,[], stat3,[], stat4,[], stat5,[], stat6,[], stat7,[], stat8,[]);
totalCellsCounted = 0;

filepathNd2 = 'C:\data\jason\peterVirusData\';

for i = 1:size(mice,1)
    mouseNum = mice(i,:);
    filename = [filepathNd2 mouseNum '\' mouseNum '_all.mat'];
    load(filename);

    if ismember(i,aggr)  % if aggression mouse
        aggrStats.ratioTdtCellsVlToDm(end+1) = mean(numGfpCellsInVl)/mean(numGfpCellsInDm);
        aggrStats.ratioPercentTdtOfGfpVlToDm(end+1) = mean(percentTdtOfGfpVl)/mean(percentTdtOfGfpDm);
        aggrStats.percentTdtOfGfpVl(end+1) = mean(percentTdtOfGfpVl);
        aggrStats.percentTdtOfGfpDm(end+1) = mean(percentTdtOfGfpDm);
        
    else                 % if control mouse
        contStats.ratioTdtCellsVlToDm(end+1) = mean(numGfpCellsInVl)/mean(numGfpCellsInDm);
        contStats.ratioPercentTdtOfGfpVlToDm(end+1) = mean(percentTdtOfGfpVl)/mean(percentTdtOfGfpDm);
        contStats.percentTdtOfGfpVl(end+1) = mean(percentTdtOfGfpVl);
        contStats.percentTdtOfGfpDm(end+1) = mean(percentTdtOfGfpDm);
        contStats.numTdtCellsInDm(end+1) = sum(numTdtCellsInDm);
        contStats.numTdtCellsInVl(end+1) = sum(numTdtCellsInVl);
        contStats.numGfpCellsInDm(end+1) = sum(numGfpCellsInDm);
        contStats.numGfpCellsInVl(end+1) = sum(numGfpCellsInVl);   
    end

    totalCellsCounted = totalCellsCounted + sum(numCellsTotal);  % multiply by approx. #sec/cell (w/ overhead for loading files, etc.) and then divide by (3600) to get the approx. # hours it would have take to count manually
end

% [mean1,std1] = grpstats(aggrStats.percentTdtOfGfpVl,[],{'mean','std'});
[mean2,std2] = grpstats(contStats.percentTdtOfGfpVl,[],{'mean','std'});
% [mean3,std3] = grpstats(aggrStats.percentTdtOfGfpDm,[],{'mean','std'});
[mean4,std4] = grpstats(contStats.percentTdtOfGfpDm,[],{'mean','std'});





