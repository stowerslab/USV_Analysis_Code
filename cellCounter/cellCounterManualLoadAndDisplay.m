% top level script to count all sections and plot results
clear all; close all;
% filepathNd2 = 'C:\data\Jason\cellCounting\crhDREADD\';
% filepathNd2 = 'C:\data\Jason\cellCounting\esrDREADD\';
% filepathNd2 = 'C:\data\Jason\cellCounting\wtDREADD\';
% filepathNd2 = 'C:\data\Jason\cellCounting\crhChR2\';
% filepathNd2 = 'C:\data\Jason\cellCounting\esrChR2\';
% filepathNd2 = 'C:\data\Jason\cellCounting\controlChR2\';
% filepathNd2 = 'C:\data\Jason\cellCounting\crhArchT\';
filepathNd2 = 'E:\Pete\nd2 test img\CGT14-PMC.nd2';

%%% INPUT:
% mouseNums = {'crh1' 'crh4' 'crh7' 'crh9' 'crh15' 'crh16'}; %3 sections  Crh DREADDs
% mouseNums = {'mc3' 'mc6'}; % 6 sections Crh DREADDs
% mouseNums = {'i1' 'i3' 'k3'}; %3 sections Esr DREADDs
% mouseNums = {'m1' 'm4' 'n1' 'n2'}; % 6 sections Esr DREADDs
% mouseNums = {'m3'}; % 5 sections Esr DREADDs
% mouseNums = {'h7' 'h8' 'h9' 'h11' 'h12'}; %3 sections WT DREADDs
% mouseNums = {'cc2\red' 'cc3\red' 'cc6\red' 'cc7\red' 'cc9\red' 'cc15\red' 'cc16\red'}; %6 sections Crh ChR2, red only
% mouseNums = {'cc2' 'cc3' 'cc6' 'cc7' 'cc9' 'cc12' 'cc13' 'cc15' 'cc16'}; %6 sections Crh ChR2, green only
% mouseNums = {'ec4' 'ec5' 'ec8' 'ec11' 'ec12' 'n7' 'n15'}; %6 sections Esr ChR2, red only
% mouseNums = {'ec7'}; %6 sections Esr ChR2 unilateral, red only
% mouseNums = {'ec10'}; %3 sections Esr ChR2, red only
% mouseNums = {'econ1' 'econ2' 'econ3'}; %6 sections control ChR2
% mouseNums = {'a2' 'a4' }; %3 sections unilateral Crh ArchT
% mouseNums = {'ac9' 'ac10' }; %6 sections unilateral Crh ArchT
mouseNums = {'a3' 'ae6' 'ae8' }; %6 sections unilateral Esr ArchT
    
% sections =  {'1R' '2R' '3R' '1L' '2L' '3L'};  %3 sections bilateral
% sections =  {'1R' '2R' '3R' '4R' '5R' '6R' '1L' '2L' '3L' '4L' '5L' '6L'}; %6 sections bilateral
% sections =  {'1R' '2R' '3R' '4R' '5R' '1L' '2L' '3L' '4L' '5L'}; %5 sections bilateral
sections =  {'1' '2' '3' '4' '5' '6'}; %6 sections unilateral
% sections =  {'1' '2' '3'}; %3 sections unilateral

numSections = size(sections, 2);
totalMice = size(mouseNums, 2);
numVertPoints = 100; % # of points as vertices
subPatchColorVl = [1 0 0]; %use red for VL suboutline
subPatchColorDm = [0 0 1]; %use blue for DM suboutline

saveOverlay = true;
countGreen = true;
countRed = true;
countBlue = false;
countRgOverlap = false;

for k = 1:totalMice
    mouseNum = mouseNums{1,k};

    %reset all counts to zero for each mouse
    redTotal = 0;
    redInBar = 0;
    redOnBarBorder = 0;
    greenTotal = 0;
    greenVL = 0;
    greenDM = 0;
    greenInBar = 0;
    greenOnBarBorder = 0;
    greenInBarVl = 0;
    greenOnBarBorderVl = 0;
    greenInBarDm = 0;
    greenOnBarBorderDm = 0;
    greenOutside = 0;
    overlapTotal = 0;

    for j = 1:numSections
        
        thisSection = sections{1,j};
        matName = [filepathNd2, mouseNum, '\', thisSection, '.mat']; %use saved data already counted
        load(matName); %careful not to save/load loop variables!!!
        
        % make DM & VL patches:
        if ellip.left % if on left of image
            vlPoints = 90:180/(numVertPoints/3):270; %+/- 90 degrees from 180
            dmPoints = -90:180/(numVertPoints/3):90;
        else % if on right of image
            vlPoints = -90:180/(numVertPoints/3):90; %+/- 90 degrees
            dmPoints = 90:180/(numVertPoints/3):270;
        end
        
        %create and view RGB overlay of ROIs, to test:
        Ioverlay = rgb2gray(Irgb(yMin:yMax,xMin:xMax,:));
        figure; hOverlay = imshow(Ioverlay); 
        hold on;
        hAxes = gca;
        %VL mask:
        [x, y] = makeEllipse(ellip.a, ellip.b, ellip.center, ellip.theta, vlPoints);
        hVlPatch = patch('Parent',hAxes,'Visible','on',...  %draw colored patch to visualize ROI
                   'XData',x,'YData',y,'EdgeColor',subPatchColorVl,'LineWidth',2,'FaceColor','none');%
        xVertVl = get(hVlPatch, 'XData');
        yVertVl = get(hVlPatch, 'YData');
        patchMaskVl = poly2mask(xVertVl, yVertVl, size(patchMask,1), size(patchMask,2));  %create binary mask from area enclosed by X & Y data

        %DM mask:
        [x, y] = makeEllipse(ellip.a, ellip.b, ellip.center, ellip.theta, dmPoints);
        hDmPatch = patch('Parent',hAxes,'Visible','on',...  %draw colored patch to visualize ROI
                   'XData',x,'YData',y,'EdgeColor',subPatchColorDm,'LineWidth',2,'FaceColor','none');%
        xVertDm = get(hDmPatch, 'XData');
        yVertDm = get(hDmPatch, 'YData');
        patchMaskDm = poly2mask(xVertDm, yVertDm, size(patchMask,1), size(patchMask,2));  %create binary mask from area enclosed by X & Y data

        if countRed
            if ~isempty(manualPointsIr)
                plot(manualPointsIr(:,1), manualPointsIr(:,2), 'r', 'LineStyle', 'none', 'Marker', '+', 'MarkerSize', 4)
            end
        end
        if countGreen
            if ~isempty(manualPointsIg)
                plot(manualPointsIg(:,1), manualPointsIg(:,2), 'g', 'LineStyle', 'none', 'Marker', '+', 'MarkerSize', 4)
            end
        end
        if countBlue
            plot(manualPointsIb(:,1), manualPointsIb(:,2), 'b', 'LineStyle', 'none', 'Marker', '+', 'MarkerSize', 2)
        end
        if countRgOverlap
            plot(colocalPoints(:,1), colocalPoints(:,2), 'y', 'LineStyle', 'none', 'Marker', 'x', 'MarkerSize', 6)
        end
        hold off;
            
        if saveOverlay 
            saveas(hOverlay, [filepathNd2, mouseNum, '\', mouseNum, '_', num2str(j), '_manualOverlay.jpg'], 'jpg');
%             saveas(hOverlay, [filepathNd2, mouseNum, '\', mouseNums{1,k}(1:end-4), '_', num2str(j), '_manualOverlayRed.jpg'], 'jpg');
        end
        
        %calculate points within polygons (Bar, BarVL, BarDM)
        if ~isempty(manualPointsIg)
            [greenInBar, greenOnBarBorder] = inpolygon(manualPointsIg(:,1),manualPointsIg(:,2), ellip.x, ellip.y);
            [greenInBarVl, greenOnBarBorderVl] = inpolygon(manualPointsIg(:,1),manualPointsIg(:,2), xVertVl, yVertVl);
            [greenInBarDm, greenOnBarBorderDm] = inpolygon(manualPointsIg(:,1),manualPointsIg(:,2), xVertDm, yVertDm);
        end
        
        if ~isempty(manualPointsIr)
            [redInBar, redOnBarBorder] = inpolygon(manualPointsIr(:,1),manualPointsIr(:,2), ellip.x, ellip.y);
            [redInBarVl, redOnBarBorderVl] = inpolygon(manualPointsIr(:,1),manualPointsIr(:,2), xVertVl, yVertVl);
            [redInBarDm, redOnBarBorderDm] = inpolygon(manualPointsIr(:,1),manualPointsIr(:,2), xVertDm, yVertDm);
        end
        
        redTotal = redTotal + sum(redInBar + redOnBarBorder);
        greenTotal = greenTotal + sum(greenInBar + greenOnBarBorder);
        greenVL = greenVL + sum(greenInBarVl + greenOnBarBorderVl);
        greenDM = greenDM + sum(greenInBarDm + greenOnBarBorderDm);
        greenOutside = greenOutside + (numGreen - sum(greenInBar + greenOnBarBorder));
        overlapTotal = numOverlap; 
        
        display(['counting ', mouseNum, ', section '  thisSection]);
        display(['numRed = ', num2str(redTotal), '; numGreen = ', num2str(greenTotal), '; numGreenVL = ', num2str(greenVL), '; numGreenDM = ', num2str(greenDM), '; numOverlap = ', num2str(overlapTotal)]);
    end    
    
    matName = [filepathNd2, mouseNum, '\', mouseNum, '_totalCounts.mat']; %save total data for each mouse
%     matName = [filepathNd2, mouseNum, '\', mouseNums{1,k}(1:end-4), '_totalCountsRed.mat']; %save total data for each mouse
    save(matName, 'greenTotal', 'greenOutside', 'greenVL', 'greenDM', 'redTotal');
end

