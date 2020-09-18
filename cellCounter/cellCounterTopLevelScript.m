% top level script to count all VMH sections and plot results 
clear all; close all;

% settings
filepathNd2 =  'E:\Sourish\For Pete 20200728\control3\'; % make sure fp has ending "\"
viewOverlay = 0; 
nd2format = 1;  % NOTE: changing this will require changing "uint8" to "uint16" in some places

%%% INPUT:
maxProj = 1; % boolean for if file is max projection or not

mouseNums = {'20200116_juv_esr1_ai6-smuf-cfos_control_3_slide_d-e-MaxIP.nd2'}; 
ROIlength = 3000;
ROIwidth = 3800;
% intensity threshold here

% % thresh = 0;
% 
% used before -> 0.019
thresh = 0.03; % threshold to detect cells using brightness/intensity
channelNums = [1,2,3]; % put the numbers of the 3 channels you want to use, such as [2,3,4], etc.

firstSections = {1};
lastSections =  {33}; % to keep track of how many sections were imaged for each mouse

roiCoord = []; % if you already have a position to set the roi rectangle, you can put it here [row, col], if not, just manually move when prompted

flipYaxis = 0; % TODO! : some sections are flipped vertically. put 1 if the section is flipped. might need to do sections 1 by 1 in this case

totalMice = size(mouseNums, 2);

% storage vars
cellCounts = {}; % saved as 3 channels' counts at a time. put in only counts of cells in dictionary with section num and channel num. eg: sec_1_ch_1
roisCropped = {}; % stores the cropped ROI of each section and channel in form like sec_1_ch_1
channelCellCoords = {}; % stores coordinates of cells
numOverlappingCells12 = 0; % overlapping cells between 1st and 2nd channel chosen
numOverlappingCells13 = 0;
numOverlappingCells23 = 0;  
nd2FileType = '.nd2';
% do all manual steps for each section first (i.e. ROIs and thresholds):
for k = 1:totalMice
    mouseNum = mouseNums{1,k};
    lastSection = lastSections{1,k};
    firstSection = firstSections{1,k};
    for j = firstSection:lastSection
        boundaries = {};
        
        if nd2format
            maxImageBitValue = 2^16-1;
            maxPixelValue = 2^12-1;
%             [Ir, Ig, Ib] = readNd2([filepathNd2 mouseNum], 1);%, num2str(j));
            channelArray = readNd2([filepathNd2 mouseNum], j, channelNums);
        else %assume JPG
            maxImageBitValue = 2^8-1;
            maxPixelValue = 2^8-1;
            [Ir, Ig, Ib, Irgb] = readJpg(filepathNd2, mouseNum, num2str(j));
        end
        
        channelCellCoords = {}; % reset coordinates of cells
        cellCounts = {}; % reset counts
        roisCropped = {}; % reset
        
        channelCounterForDisp = 1;
        for m = 1:length(channelArray) % m = 1:2:length(channelArray), step-size might have to be 2 for jp2 files
            Ich = channelArray{m};
            
%             if ~contains(mouseNum, nd2FileType)
%                 IchLower = channelArray{m+1};
%                 [roiPatchArrayR, xMin, xMax, yMin, yMax, thresh, I_bw, channelCellCoord, roiCoord, roiCropped] = sectionSetRoiThresh(maxImageBitValue, maxPixelValue, IchLower, Ich, roiCoord); % do all manual steps for each section first (i.e. ROIs and thresholds)
%             end
            [roiPatchArrayR, xMin, xMax, yMin, yMax, thresh, I_bw, channelCellCoord, roiCoord, roiCropped, boundary] = sectionSetRoiThresh(maxImageBitValue, maxPixelValue, Ich, Ich, roiCoord, ROIlength, ROIwidth, thresh); % do all manual steps for each section first (i.e. ROIs and thresholds)
%             if flipYaxis
%             end
            boundaries{m} = boundary;
            clearvars boundary;
            disp(['(Index ' num2str(channelCounterForDisp) ' of channelNums) # of cells:']);
            disp(length(channelCellCoord));
            cellCounts(end+1,:) = {['sec_' num2str(j) '_ch_' channelCounterForDisp],length(channelCellCoord)}; % store counts
            channelCellCoords{end+1} = channelCellCoord;
            roisCropped(end+1,:) = {['sec_' num2str(j) '_ch_' channelCounterForDisp],roiCropped}; % store cropped image
            channelCounterForDisp = channelCounterForDisp + 1;
        end
        
        
%         
%         [roiPatchArrayR, xMin, xMax, yMin, yMax, threshR, Ir_bw, channel1CellCoord, roiCoord, roiCropped1] = sectionSetRoiThresh(maxImageBitValue, maxPixelValue, Ir, Ir, roiCoord); % do all manual steps for each section first (i.e. ROIs and thresholds)
%         disp("ch 1 # of cells");
%         disp(length(channel1CellCoord));
%         [roiPatchArrayG, xMin, xMax, yMin, yMax, threshG, Ig_bw, channel2CellCoord, roiCoord, roiCropped2] = sectionSetRoiThresh(maxImageBitValue, maxPixelValue, Ig, Ig, roiCoord); 
%         disp("ch 2 # of cells");
%         disp(length(channel2CellCoord));
%         %         [roiPatchArrayB, xMin, xMax, yMin, yMax, threshB, Ib_bw] = sectionSetRoiThresh(maxImageBitValue, maxPixelValue, Ir, Ib); 
%         [roiPatchArrayB, xMin, xMax, yMin, yMax, threshB, Ib_bw,channel3CellCoord, roiCoord, roiCropped3] = sectionSetRoiThresh(maxImageBitValue, maxPixelValue, Ib, Ib, roiCoord); 
%         disp("ch 3 # of cells");
%         disp(length(channel3CellCoord));
        % if there are more channels, have to change readNd2.m to read more in
        
        % overlay roi from channels
        % 1st and 2nd and 3rd
%         C = imfuse(roisCropped{1},roisCropped{2},'falsecolor','Scaling','independent','ColorChannels',[1 2 0]);
%         figure;
%         C = imresize(C,0.5);
%         image(C)
        
        %%
        % pete: disregard overlapping cells by checking coords for each
        % channel and delete. just check distance between each coordinate 
        % from other array and count once.
        
        % when indexing into channelCellCoords, order depends on number of
        % channelNums used. TODO: make it so you don't have to manually
        % change
        if ismember(1, channelNums) 
            channel1CellCoord = channelCellCoords{1}; 
        end
        if ismember(2, channelNums) 
            if length(channelNums) == 1
                channel1CellCoord = channelCellCoords{1};
            end
            if length(channelNums) >= 2
                channel2CellCoord = channelCellCoords{2};
            end
        end
        if ismember(3, channelNums) 
            if length(channelNums) == 1
                channel1CellCoord = channelCellCoords{1};
            end
            if length(channelNums) == 2
                channel2CellCoord = channelCellCoords{2};
            end
            if length(channelNums) == 3
                channel3CellCoord = channelCellCoords{3};
            end
%             channel3CellCoord = channelCellCoords{3}; % fixed this comment -- usually 3. the index to access depends on the number of channels being used
        end

%         numOverlappingCells12 = 0; % reset counts
%         numOverlappingCells13 = 0;
%         numOverlappingCells23 = 0;
% 
%         overlapThresh = 1; % just the distance which counts as overlap
%         
%         for m = 1:length(channel1CellCoord) 
%             for n = 1:length(channel2CellCoord)
%                 [row1, col1] = deal(channel1CellCoord{m}(1), channel1CellCoord{m}(2));
%                 [row2, col2] = deal(channel2CellCoord{n}(1), channel2CellCoord{n}(2));
%                 rowDiffSquared = (row1 - row2).^2;
%                 colDiffSquared = (col1 - col2).^2;
%                 if sqrt(rowDiffSquared + colDiffSquared) <= overlapThresh % ch1 check ch2
%                     numOverlappingCells12 = numOverlappingCells12 + 1;
% %                     disp(['overlap in ' num2str(channelNums(1)) ' ' num2str(channelNums(2))]);
% %                     disp([row1, col1]);
%                 end
%             end
%             
%             for p = 1:length(channel3CellCoord) 
%                 [row1, col1] = deal(channel3CellCoord{p}(1), channel3CellCoord{p}(2));
%                 [row2, col2] = deal(channel1CellCoord{m}(1), channel1CellCoord{m}(2));
%                 rowDiffSquared = (row1 - row2).^2;
%                 colDiffSquared = (col1 - col2).^2;
%                 if sqrt(rowDiffSquared + colDiffSquared) <= overlapThresh % ch1 check ch3
%                     numOverlappingCells13 = numOverlappingCells13 + 1;
% %                     disp(['overlap in ', num2str(channelNums(1)), ' ', num2str(channelNums(3))]);
% %                     disp([row1, col1]);
%                 end
%             end
%         end
% 
%         for o = 1:length(channel3CellCoord) 
%             for p = 1:length(channel2CellCoord) 
%                 [row1, col1] = deal(channel3CellCoord{o}(1), channel3CellCoord{o}(2));
%                 [row2, col2] = deal(channel2CellCoord{p}(1), channel2CellCoord{p}(2));
%                 rowDiffSquared = (row1 - row2).^2;
%                 colDiffSquared = (col1 - col2).^2;
%                 if sqrt(rowDiffSquared + colDiffSquared) <= overlapThresh % ch2 check ch3
%                     numOverlappingCells23 = numOverlappingCells23 + 1;
% %                     disp(['overlap in ', num2str(channelNums(2)), ' ', num2str(channelNums(3))]);
% %                     disp([row2, col2]);
%                 end
%             end
%         end
%             
%         disp(['num overlapping cells channels ', num2str(channelNums(1)), ' ', num2str(channelNums(3)), ': ', num2str(numOverlappingCells13)]);
%         disp(['num overlapping cells channels ',  num2str(channelNums(1)), ' ', num2str(channelNums(2)),  ': ', num2str(numOverlappingCells12)]);
%         disp(['num overlapping cells channels ',  num2str(channelNums(2)), ' ', num2str(channelNums(3)),  ': ', num2str(numOverlappingCells23)]);
%%
%         matName = [filepathNd2, mouseNum, '\', mouseNum, '_', num2str(j), '.mat']; %save data for each individual section:
        matName = [filepathNd2, mouseNum, '_', num2str(j), '.mat']; % pete
        save(matName);
%         close all; %prevent large number figures
        
        if viewOverlay % to create and view RGB overlay of ROIs, to test:
%             figure; imshow(Irgb(yMin:yMax,xMin:xMax,:));
%             hold on;
%             Ioverlay = rgb2gray(Irgb(yMin:yMax,xMin:xMax,:)); 
%             h = imshow(Ioverlay); 
%             hold off;
%             set(h, 'AlphaData', ~Ir_bw(yMin:yMax,xMin:xMax,:));
            [Xr,Yr] = find(Ir_bw(yMin:yMax,xMin:xMax));
            [Xg,Yg] = find(Ig_bw(yMin:yMax,xMin:xMax));
            
            Ioverlay = rgb2gray(Irgb(yMin:yMax,xMin:xMax,:));
            figure; hOverlay = imshow(Ioverlay); 
            hold on;
            plot(Yr, Xr, 'r', 'LineStyle', 'none', 'Marker', '+', 'MarkerSize', 4);
            plot(Yg, Xg, 'g', 'LineStyle', 'none', 'Marker', '+', 'MarkerSize', 4);
%             plot(manualPointsIb(:,1), manualPointsIb(:,2), 'b', 'LineStyle', 'none', 'Marker', '+')
%             plot(colocalPoints(:,1), colocalPoints(:,2), 'y', 'LineStyle', 'none', 'Marker', 'x')
            hold off;
            saveas(hOverlay, [filepathNd2, mouseNum, '\', mouseNum, '_', num2str(j), '_manualOverlay.jpg'], 'jpg');
        end  
%    clear all;
    end
end


