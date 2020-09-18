%script to read in VMH data and plot
clear all; close all;

% there are several problems with this data:
%       - axons are brighter at injection locus (esp. bad w/ non-Cre-dep GFP), so more cells likely detected at this location b/c of averaging
%       - about 1/4 to 1/2 of all sections had bad striping problem on Nikon C2, which does not affect DM & Vl equally and makes it hard to find DM & VL

mice = ['p1';'p2';'p3';'p4';'p5';'p6';'p7'];
aggr = [1 3 6]; %indeces of aggression mice

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

filepathNd2 = 'C:\data\jason\peterVirusData\firstAggrTmpDmso\';

for i = 1:size(mice,1)
    mouseNum = mice(i,:);
    filename = [filepathNd2 mouseNum '\' mouseNum '_all.mat'];
    load(filename);

    if ismember(i,aggr)  % if aggression mouse
        aggrStats.ratioTdtCellsVlToDm(end+1) = mean(numGfpCellsInVl)/mean(numGfpCellsInDm);
        aggrStats.ratioPercentTdtOfGfpVlToDm(end+1) = mean(percentTdtOfGfpVl)/mean(percentTdtOfGfpDm);
        aggrStats.percentTdtOfGfpVl(end+1) = mean(percentTdtOfGfpVl);
        aggrStats.percentTdtOfGfpDm(end+1) = mean(percentTdtOfGfpDm);
        aggrStats.numTdtCellsInVl(end+1) = mean(numGfpCellsInVl);
        
    else                 % if control mouse
        contStats.ratioTdtCellsVlToDm(end+1) = mean(numGfpCellsInVl)/mean(numGfpCellsInDm);
        contStats.ratioPercentTdtOfGfpVlToDm(end+1) = mean(percentTdtOfGfpVl)/mean(percentTdtOfGfpDm);
        contStats.percentTdtOfGfpVl(end+1) = mean(percentTdtOfGfpVl);
        contStats.percentTdtOfGfpDm(end+1) = mean(percentTdtOfGfpDm);
        contStats.numTdtCellsInVl(end+1) = mean(numGfpCellsInVl);
    end

    totalCellsCounted = totalCellsCounted + sum(numCellsTotal);  % multiply by approx. #sec/cell (w/ overhead for loading files, etc.) and then divide by (3600) to get the approx. # hours it would have take to count manually
end

[mean1,std1] = grpstats(aggrStats.percentTdtOfGfpVl,[],{'mean','std'});
[mean2,std2] = grpstats(contStats.percentTdtOfGfpVl,[],{'mean','std'});
[mean3,std3] = grpstats(aggrStats.percentTdtOfGfpDm,[],{'mean','std'});
[mean4,std4] = grpstats(contStats.percentTdtOfGfpDm,[],{'mean','std'});

% BOXPLOT: On each box, the central mark is the median, the edges of the box are the 25th and 75th percentiles, 
% the whiskers extend to the most extreme data points not considered
% outliers, and outliers are plotted individually:
% figure;
% boxplot([contStats.percentTdtOfGfpVl' aggrStats.percentTdtOfGfpVl'],'labels',{'control','aggression'}); %'boxstyle','filled'
% title('% TDT of GFP in VMHvl');
% 
% figure;
% boxplot([contStats.percentTdtOfGfpDm' aggrStats.percentTdtOfGfpDm'],'labels',{'control','aggression'});
% title('% TDT of GFP in VMHdm');

figure;
barvalues = [mean2 mean1];
errors = [std2 std1];
width = 1;
groupnames = {};
bw_title = '% TDT of GFP in VMHvl';
bw_xlabel = [];
bw_ylabel = '%';
bw_colormap = jet;
gridstatus = 'none';
bw_legend = {'control','aggression'};
error_sides = 2;
legend_type = 'plot';
handles1 = barweb(barvalues, errors, width, groupnames, bw_title, bw_xlabel, bw_ylabel, bw_colormap, gridstatus, bw_legend, error_sides,legend_type);

figure;
barvalues = [mean4 mean3];
errors = [std4 std3];
bw_title = '% TDT of GFP in VMHdm';
handles2 = barweb(barvalues, errors, width, groupnames, bw_title, bw_xlabel, bw_ylabel, bw_colormap, gridstatus, bw_legend, error_sides,legend_type);








