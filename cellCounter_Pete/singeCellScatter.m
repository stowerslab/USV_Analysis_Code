% FOR CELL COUNTER/ HEATMAP. get cell coordinates and scatter
clear all;
close all;
% Input:
sectionFilePath = 'E:\Sourish\FOR PETE\'; % select directory with the .mat files
ROIwidth = 1800;
ROIlength = 1800;
% ----------
sectionDirectory = dir(fullfile(sectionFilePath,'*.mat'));
ch3_cell_coords = {};
ch2_cell_coords = {};
ch1_cell_coords = {};

ch2_names = {}; % (look at img for channel color order) put all names of files you want the (2nd) green channel for

filenameff = '7528_PAG_1-16'; % the part before "-MaxIP.nd2_<num>" etc
% numberTrial = '8150';

ch3_names = {};
ch1_names = {};
startN = 1
endN = 16

% TODO: not have to manually change the end index between "startN" and
% "endN" when plotting single vs multiple

for i = startN:endN
    if i < 10
        num = ['0' num2str(i)];
    else
        num = num2str(i);
    end
    names = dir(fullfile(sectionFilePath,[filenameff '-MaxIP.nd2_' num '*.mat']))
    names
    ch3_names{end+1} = names.name;
end

ch3_names
% ch1_names = {'20200203_vgl_pag_gi_8034-MaxIP.nd2_1.mat';
%     '20200203_vgl_pag_gi_8034-MaxIP.nd2_2.mat'
%     '20200203_vgl_pag_gi_8034-MaxIP.nd2_3.mat'
%     '20200203_vgl_pag_gi_8034-MaxIP.nd2_4.mat'
%     '20200203_vgl_pag_gi_8034-MaxIP.nd2_5.mat'
%     '20200203_vgl_pag_gi_8034-MaxIP.nd2_6.mat'
%     '20200203_vgl_pag_gi_8034-MaxIP.nd2_7.mat'
%     '20200203_vgl_pag_gi_8034-MaxIP.nd2_8.mat'
%     '20200203_vgl_pag_gi_8034-MaxIP.nd2_9.mat'
%     '20200203_vgl_pag_gi_8034-MaxIP.nd2_10.mat'
%     '20200203_vgl_pag_gi_8034-MaxIP.nd2_11.mat'
%     '20200203_vgl_pag_gi_8034-MaxIP.nd2_12.mat'
%     '20200203_vgl_pag_gi_8034-MaxIP.nd2_13.mat'
%     '20200203_vgl_pag_gi_8034-MaxIP.nd2_14.mat'
%     '20200203_vgl_pag_gi_8034-MaxIP.nd2_15.mat'
% };% (look at img for channel color order)

% ch3_names = {};

% 2nd chan
for sectionData = sectionDirectory'
    if ~any(strcmp(sectionData.name, ch2_names)) % if these are the same, continue
        continue; % SKIP
    end
    
    load(fullfile(sectionData.folder, sectionData.name)); 
    disp(sectionData.name);
    
    if contains(sectionData.name, 'flip')
        for cellNum = 1:length(channel2CellCoord)
            if length(channel2CellCoord{cellNum}) > 2 
                    % NOTE !!! put in the height of the rectangle ROI here V <-pointer 
                ch2_cell_coords{end + 1} = [ROIwidth -channel2CellCoord{cellNum}(1,1), ROIlength-channel2CellCoord{cellNum}(1,2)]; % eg, height of ROI is 800, replace with your own
            else
                ch2_cell_coords{end + 1} = [ROIwidth-channel2CellCoord{cellNum}(1), ROIlength-channel2CellCoord{cellNum}(2)]; % channel2CellCoord{cellNum};
            end
        end
    else
        for cellNum = 1:length(channel2CellCoord)
            if length(channel2CellCoord{cellNum}) > 2 
                    % NOTE !!! put in the height of the rectangle ROI here V <-pointer 
                ch2_cell_coords{end + 1} = [channel2CellCoord{cellNum}(1,1), ROIlength-channel2CellCoord{cellNum}(1,2)]; % eg, height of ROI is 800, replace with your own
            else
                ch2_cell_coords{end + 1} = [channel2CellCoord{cellNum}(1), ROIlength-channel2CellCoord{cellNum}(2)]; % channel2CellCoord{cellNum};
            end
        end
    end
        
        
   

    % THIS ONLY ADDS CHANNEL 2, you can copy past and change to channel 1
    
end



figure;
mat = reshape(cell2matPETE(ch2_cell_coords),2,size(ch2_cell_coords,2));
scatter(mat(1,:),mat(2,:),8,'red', 'filled','MarkerFaceAlpha',.4,'MarkerEdgeAlpha',.4);
% ALSO PUT ROI DIMENSIONS HERE
xlim([0, ROIwidth]);
ylim([0, ROIlength]); 

title("Cell Coordinate Scatter");

hold on;


%% 1st chan

for sectionData = sectionDirectory'
    if ~any(strcmp(sectionData.name, ch1_names)) % if these are the same, continue
        continue;
    end
    load(fullfile(sectionData.folder, sectionData.name)); 
    disp(sectionData.name);

    for cellNum = 1:length(channel1CellCoord)
        if length(channel1CellCoord{cellNum}) > 2 
                % NOTE !!! put in the height of the rectangle ROI here V <-pointer 
            ch1_cell_coords{end + 1} = [channel1CellCoord{cellNum}(1,1), ROIlength-channel1CellCoord{cellNum}(1,2)]; % eg, height of ROI is 800, replace with your own
        else
            ch1_cell_coords{end + 1} = [channel1CellCoord{cellNum}(1), ROIlength-channel1CellCoord{cellNum}(2)]; % channel2CellCoord{cellNum};
        end
    end
end


mat = reshape(cell2matPETE(ch1_cell_coords),2,size(ch1_cell_coords,2));
scatter(mat(1,:),mat(2,:),8,'red', 'filled','MarkerFaceAlpha',.4,'MarkerEdgeAlpha',.4);

%% 3rd chan

for sectionData = sectionDirectory'
    if ~any(strcmp(sectionData.name, ch3_names)) % if these are the same, continue
        continue; % SKIP
    end
    
    load(fullfile(sectionData.folder, sectionData.name)); 
    disp(sectionData.name);
    
    if contains(sectionData.name, 'flip')
        for cellNum = 1:length(channel3CellCoord)
            if length(channel3CellCoord{cellNum}) > 2 
                    % NOTE !!! put in the height of the rectangle ROI here V <-pointer 
                ch3_cell_coords{end + 1} = [ROIwidth -channel3CellCoord{cellNum}(1,1), ROIlength-channel3CellCoord{cellNum}(1,2)]; % eg, height of ROI is 800, replace with your own
            else
                ch3_cell_coords{end + 1} = [ROIwidth-channel3CellCoord{cellNum}(1), ROIlength-channel3CellCoord{cellNum}(2)]; % channel2CellCoord{cellNum};
            end
        end
    else
        for cellNum = 1:length(channel3CellCoord)
            if length(channel3CellCoord{cellNum}) > 2 
                    % NOTE !!! put in the height of the rectangle ROI here V <-pointer 
                ch3_cell_coords{end + 1} = [channel3CellCoord{cellNum}(1,1), ROIlength-channel3CellCoord{cellNum}(1,2)]; % eg, height of ROI is 800, replace with your own
            else
                ch3_cell_coords{end + 1} = [channel3CellCoord{cellNum}(1), ROIlength-channel3CellCoord{cellNum}(2)]; % channel2CellCoord{cellNum};
            end
        end
    end
        
end

mat = reshape(cell2matPETE(ch3_cell_coords),2,size(ch3_cell_coords,2));
scatter(mat(1,:),mat(2,:),8,'red', 'filled','MarkerFaceAlpha',.4,'MarkerEdgeAlpha',.4);

daspect([1 1 1])
saveas(gcf, fullfile(sectionFilePath, [ch3_names{1} '.fig']) ) % sectionData.name 

% ALSO PUT ROI DIMENSIONS HERE
% xlim([0, 1600]);
% ylim([0, 1600]); 

% title("Ch 1: Cell Coordinate Scatter");