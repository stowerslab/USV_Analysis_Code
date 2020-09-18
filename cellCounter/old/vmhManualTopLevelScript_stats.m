% top level script to read count data from all VMH sections and plot results
clear all; close all;

filepathRoot = 'C:\data\jason\peterVirusData\vmhIso';  %root for all mice

%%% INPUT folder names for individual mice:
mouseNums = {'v5','v6','v7','v8','v12','v15','v16','v17','v20',...
             'v21','v23','v24','v25','v26','v28','v29','v30'};  %ex. 2nd entry is mouseNums{1,2}
aggr = [3 4 7 10 11]; %indeces of aggression mice
cont = [1 5 6 9 12 13]; %indeces of control mice
leak = [2 8 14 15 16 17]; %indeces of leak mice

%All stats are for VMHvl only; declare fields to be calculated:
stat1 = 'totalRedCells';
stat2 = 'ratioRedToGreenCells';
stat3 = 'totalAvgRed';
stat4 = 'totalAvgRedToGreen';
stat5 = 'stat5';
stat6 = 'stat6';
stat7 = 'stat7';
stat8 = 'stat8';

%initialize total stat arrays:
aggrStats = struct(stat1,[], stat2,[], stat3,[], stat4,[], stat5,[], stat6,[], stat7,[], stat8,[]);
contStats = struct(stat1,[], stat2,[], stat3,[], stat4,[], stat5,[], stat6,[], stat7,[], stat8,[]);
leakStats = struct(stat1,[], stat2,[], stat3,[], stat4,[], stat5,[], stat6,[], stat7,[], stat8,[]);

totalMice = size(mouseNums, 2);

%% Read in cell counts for each channel and save summary data for each mouse:
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
for k = 1:totalMice
    mouseNum = mouseNums{1,k};
    %reset each mouse to zero:
    totalNumGreen = 0; totalNumBlue = 0; totalNumRed = 0;
    totalAvgGreen = zeros(4,1); totalAvgBlue = zeros(4,1); totalAvgRed = zeros(4,1);
    for j = 1:4
        matName = [filepathRoot, '\', mouseNum, '\', mouseNum, '_', num2str(j), '_greenCount.mat']; %GREEN
        load(matName);
        totalNumGreen = totalNumGreen + numGreen;
        totalAvgGreen(j) = avgGreen;
        matName = [filepathRoot, '\', mouseNum, '\', mouseNum, '_', num2str(j), '_blueCount.mat']; %BLUE
        load(matName);
        totalNumBlue = totalNumBlue + numBlue;
        totalAvgBlue(j) = avgBlue;
        matName = [filepathRoot, '\', mouseNum, '\', mouseNum, '_', num2str(j), '_redCount.mat']; %RED
        load(matName);
        totalNumRed = totalNumRed + numRed;
        totalAvgRed(j) = avgRed;
    end
    %save data for each mouse:
    totalMatName = [filepathRoot, '\', mouseNum, '\', mouseNum, '_totalCount.mat'];
    save(totalMatName, '-regexp', '^total');  %save all variables starting with 'total', etc.
end

%% Read in summary data for all mice and compute stats, plot:
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
for k = 1:totalMice
    mouseNum = mouseNums{1,k};
    totalMatName = [filepathRoot, '\', mouseNum, '\', mouseNum, '_totalCount.mat'];
    load(totalMatName);
    
    if ismember(k,aggr)  % if aggression mouse
        aggrStats.totalRedCells(end+1) = totalNumRed;
        aggrStats.ratioRedToGreenCells(end+1) = totalNumRed/totalNumGreen;
        aggrStats.totalAvgRed(end+1) = mean(totalAvgRed); %NOTE: variability across sections is collapsed here, so we're loosing some statistical power
        aggrStats.totalAvgRedToGreen(end+1) = mean(totalAvgRed)/mean(totalAvgGreen);
        %
    elseif ismember(k,cont)                % if control mouse
        contStats.totalRedCells(end+1) = totalNumRed;
        contStats.ratioRedToGreenCells(end+1) = totalNumRed/totalNumGreen;
        contStats.totalAvgRed(end+1) = mean(totalAvgRed); 
        contStats.totalAvgRedToGreen(end+1) = mean(totalAvgRed)/mean(totalAvgGreen);
        %
    else % if leak mouse
        leakStats.totalRedCells(end+1) = totalNumRed;
        leakStats.ratioRedToGreenCells(end+1) = totalNumRed/totalNumGreen;
        leakStats.totalAvgRed(end+1) = mean(totalAvgRed);
        leakStats.totalAvgRedToGreen(end+1) = mean(totalAvgRed)/mean(totalAvgGreen);
        %
    end   
end

%% Compute stats
[mean1,std1] = grpstats(aggrStats.totalRedCells,[],{'mean','std'});
[mean2,std2] = grpstats(contStats.totalRedCells,[],{'mean','std'});
[mean3,std3] = grpstats(leakStats.totalRedCells,[],{'mean','std'});
[s_NumRedCells,p_NumRedCells] = ttest2(aggrStats.totalRedCells,contStats.totalRedCells); %alternative hypothesis is that the data in x and y comes from populations with unequal means

[mean4,std4] = grpstats(aggrStats.ratioRedToGreenCells,[],{'mean','std'});
[mean5,std5] = grpstats(contStats.ratioRedToGreenCells,[],{'mean','std'});
[mean6,std6] = grpstats(leakStats.ratioRedToGreenCells,[],{'mean','std'});
[s_ratioRedToGreenCells,p_ratioRedToGreenCells] = ttest2(aggrStats.ratioRedToGreenCells,contStats.ratioRedToGreenCells);

[mean7,std7] = grpstats(aggrStats.totalAvgRed,[],{'mean','std'});
[mean8,std8] = grpstats(contStats.totalAvgRed,[],{'mean','std'});
[mean9,std9] = grpstats(leakStats.totalAvgRed,[],{'mean','std'});
[s_totalAvgRed,p_totalAvgRed] = ttest2(aggrStats.totalAvgRed,contStats.totalAvgRed);

[mean10,std10] = grpstats(aggrStats.totalAvgRedToGreen,[],{'mean','std'});
[mean11,std11] = grpstats(contStats.totalAvgRedToGreen,[],{'mean','std'});
[mean12,std12] = grpstats(leakStats.totalAvgRedToGreen,[],{'mean','std'});
[s_totalAvgRedToGreen,p_totalAvgRedToGreen] = ttest2(aggrStats.totalAvgRedToGreen,contStats.totalAvgRedToGreen);


%% Make plots
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

%BAR CHARTS: ('barweb' is from MATLAB Central; also see 'sigstar')
figure;
barvalues = [mean1 mean2 mean3];
errors = [std1 std2 std3];
width = 1;
groupnames = {};
bw_title = ['# TDT cells in VMHvl, p = ', num2str(p_NumRedCells)];
bw_xlabel = [];
bw_ylabel = '# cells';
bw_colormap = jet;
gridstatus = 'none';
bw_legend = {'aggression','control','leak'};
error_sides = 2;
legend_type = 'plot';
handles1 = barweb(barvalues, errors, width, groupnames, bw_title, bw_xlabel, bw_ylabel, bw_colormap, gridstatus, bw_legend, error_sides,legend_type);

figure;
barvalues = [mean4 mean5 mean6];
errors = [std4 std5 std6];
bw_title = ['% TDT of GFP cells in VMHvl, p = ', num2str(p_ratioRedToGreenCells)];
bw_ylabel = 'TDT/GFP cell ratio';
handles2 = barweb(barvalues, errors, width, groupnames, bw_title, bw_xlabel, bw_ylabel, bw_colormap, gridstatus, bw_legend, error_sides,legend_type);

figure;
barvalues = [mean7 mean8 mean9];
errors = [std7 std8 std9];
bw_title = ['Avg. Red Fluo. Intensity in VMHvl, p = ', num2str(p_totalAvgRed)];
bw_ylabel = 'avg. fluo. intensity (a.u.)';
handles3 = barweb(barvalues, errors, width, groupnames, bw_title, bw_xlabel, bw_ylabel, bw_colormap, gridstatus, bw_legend, error_sides,legend_type);

figure;
barvalues = [mean10 mean11 mean12];
errors = [std10 std11 std12];
bw_title = ['Avg. Red/Green Fluo. Intensity in VMHvl, p = ', num2str(p_totalAvgRedToGreen)];
bw_ylabel = 'TDT/GFP fluo. intensity ratio (a.u.)';
handles4 = barweb(barvalues, errors, width, groupnames, bw_title, bw_xlabel, bw_ylabel, bw_colormap, gridstatus, bw_legend, error_sides,legend_type);

































