clear all; close all;

filepathNd2Crh = 'C:\data\Jason\microscope\2017_02_CrhSpinalRatioCounts\';
mouseNumsCrh = {'a2', 'a4', 'ac9', 'ac10', 'cc5', 'cc6', 'cc7', 'cc9'};
numSectionsCrh = {[1:1:20], [1:1:18], [1:1:33], [1:1:23], [1:1:28], [1:1:20], [1:1:28], [1:1:30]};
totalMiceCrh = size(mouseNumsCrh, 2);

filepathNd2 = 'C:\data\Jason\microscope\2017_02_EsrSpinalRatioCounts\';
mouseNums = {'a3', 'ae6', 'ae8', 'ec4', 'ec5', 'ec7'}; %nd2 paths
numSections = {[1:1:12], [1:1:31], [1:1:13], [1:1:25], [1:1:27], [1:1:27]};
totalMice = size(mouseNums, 2);

intRatioCrh = [];
leftIntCrh = [];
middleIntCrh = [];
rightIntCrh = [];
intRatio = [];
leftInt = [];
middleInt = [];
rightInt = [];

axonThresh = 1000000; %make sure there is not just 1 axon in image
grayMatterThresh = 500; %light threshold for gray matter autofluorescence
    
for k = 1:totalMice
    mouseNum = mouseNums{1,k};
    tframes = numSections{1,k};
    
    for j = 1:max(tframes) %for all sections
        matName = [filepathNd2, mouseNum, '_Section', num2str(j), '.mat'];
        load(matName); %currentNisslRot, currentAxonRot, roiPatch
        
        xmin = roiPatch.patchPosition(1);
        ymin = roiPatch.patchPosition(2);
        width = roiPatch.patchPosition(3);
        height = roiPatch.patchPosition(4);
        
        %QC for ROIs, comment out normally:
%         currentAxonRot(floor(ymin:ymin+height), floor(xmin:xmin+width/3)) = 2^12 - 1;  %LEFT
%         imagesc(currentAxonRot);
%         currentAxonRot(floor(ymin:ymin+height), floor(xmin+width/3:xmin+2*(width/3))) = 2^12 - 1;  %MID
%         imagesc(currentAxonRot);
%         currentAxonRot(floor(ymin:ymin+height), floor(xmin+2*(width/3):xmin+3*(width/3))) = 2^12 - 1;  %RIGHT
%         imagesc(currentAxonRot);
        currentAxonRot(currentAxonRot<grayMatterThresh) = 0; %threshold gray matter noise, which can skew results
        
        leftInt(end+1) = sum(sum(currentAxonRot(floor(ymin:ymin+height), floor(xmin:xmin+width/3)))); %#ok<*SAGROW>
        middleInt(end+1) = sum(sum(currentAxonRot(floor(ymin:ymin+height), floor(xmin+width/3:xmin+2*(width/3)))));
        rightInt(end+1) = sum(sum(currentAxonRot(floor(ymin:ymin+height), floor(xmin+2*(width/3):xmin+3*(width/3)))));
        
        % could exclude here on the basis of very little labelling, say if
        if leftInt(end)+middleInt(end)+rightInt(end) > axonThresh
            intRatio(end+1) = middleInt(end)/(leftInt(end)+rightInt(end));
        end    
    end
    meanMouseIntRatio(k) = mean(intRatio);
end


for k = 1:totalMiceCrh
    mouseNum = mouseNumsCrh{1,k};
    tframes = numSectionsCrh{1,k};
    
    for j = 1:max(tframes) %for all sections
        matName = [filepathNd2Crh, mouseNum, '_Section', num2str(j), '.mat'];
        load(matName); %currentNisslRot, currentAxonRot, roiPatch
        
        xmin = roiPatch.patchPosition(1);
        ymin = roiPatch.patchPosition(2);
        width = roiPatch.patchPosition(3);
        height = roiPatch.patchPosition(4);
        
        currentAxonRot(currentAxonRot<grayMatterThresh) = 0;
        
        leftIntCrh(end+1) = sum(sum(currentAxonRot(floor(ymin:ymin+height), floor(xmin:xmin+width/3)))); %#ok<*SAGROW>
        middleIntCrh(end+1) = sum(sum(currentAxonRot(floor(ymin:ymin+height), floor(xmin+width/3:xmin+2*(width/3)))));
        rightIntCrh(end+1) = sum(sum(currentAxonRot(floor(ymin:ymin+height), floor(xmin+2*(width/3):xmin+3*(width/3)))));
        
        % could exclude here on the basis of very little labelling, say if
        if leftIntCrh(end)+middleIntCrh(end)+rightIntCrh(end) > axonThresh
            intRatioCrh(end+1) = middleIntCrh(end)/(leftIntCrh(end)+rightIntCrh(end));
        end    
    end
    meanMouseIntRatioCrh(k) = mean(intRatioCrh);
end

cmapG = [0 0.8 0];
cmapM = [0.8 0 0.8];
[meanRatio, stdRatio, semRatio] = grpstats(intRatio,[],{'mean','std','sem'});
[meanRatioCrh, stdRatioCrh, semRatioCrh] = grpstats(intRatioCrh,[],{'mean','std','sem'});

hRatios = figure;
hold on;

[x,y] = vscatter(intRatio, 0.08);
plot_little_circles(x-1, y, 0.04, cmapG, 0.9);
% plot(x, y, 'g*')

[xCrh,yCrh] = vscatter(intRatioCrh, 0.08);
plot_little_circles(xCrh+1, yCrh, 0.04, cmapM, 0.9);
% plot(xCrh+1, yCrh, 'm*')

% p = ranksum(intRatio,intRatioCrh)  %Mann-Whitney U Test
p = ranksum(meanMouseIntRatio,meanMouseIntRatioCrh) %more conservative

% barvalues = [meanRatio];
% errors = [semRatio];
% width = 0.25;
% groupnames = {};
% bw_title = [];
% bw_xlabel = [];
% bw_ylabel = 'intensity M / (R+L)';
% bw_colormap = [cmapG]; %; cmapM];
% gridstatus = 'none';
% bw_legend = {};
% error_sides = 2;
% legend_type = 'plot';
% handles1 = barweb(barvalues, errors, width, groupnames, bw_title, bw_xlabel, bw_ylabel, bw_colormap, gridstatus, bw_legend, error_sides,legend_type);

plot(-1, meanRatio, 'k.', 'MarkerSize', 40)
errorbar(-1, meanRatio, stdRatio, 'Color', 'k', 'LineWidth', 4)
plot(1, meanRatioCrh, 'k.', 'MarkerSize', 40)
errorbar(1, meanRatioCrh, stdRatioCrh, 'Color', 'k', 'LineWidth', 4) 
axis([-2 2.5 0 4.5])
hold off

%plot2svg([filepathNd2, 'EsrVsCrhScRatios.svg'], gcf)