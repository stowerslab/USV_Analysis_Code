function [ ] = Auto_Heat_Mapper(image,data,file)
%% Ulloa, Andres    Auto Heat Mapper    11/23/15

% The purpose of this script is to take eye traking data from multiple
% participants looking at a single image, smooth it using kernel density 
% estimation assuming a normal distribution, and overlay it on the image

% Image and eye tracking data source 
% [http://www.csc.kth.se/~kootstra/index.php?item=215&menu=200]

% clear all
close all

%% Image and Eye Tracking Data Import 
% eyeTrackData = data.eyeTrackData;
% C            = file;
%% Define as X and Y tracking coordinates
% NA  = file(65:66)
% Xstructname = strcat('animals_',cellstr(NA));
% d1  = getfield(eyeTrackData,'animals',Xstructname{1},'subject_01','fixX')
% d2  = getfield(eyeTrackData,'animals',Xstructname{1},'subject_01','fixY')


C = 'E:\Pete\white.jpeg'; % image background

% Input:
sectionFilePath = 'E:\Jingyi\Dropbox (Scripps Research)\Jingyi Data backup\Imaging\2018-5 WT DREADD infection check\GROUP 1\sec3\';

% ------
sectionDirectory = dir(fullfile(sectionFilePath,'*.mat'));

all_cell_coords = {};

FLIPPED_SECTIONS = {'LEFT_INVERTED_3L 2R (Inverted)-MaxIP.jp2_1.mat';
    'RIGHT_INVERTED_3L 2R (Inverted)-MaxIP.jp2_1.mat'
    };

for sectionData = sectionDirectory'
    load(fullfile(sectionData.folder, sectionData.name)); 
    disp(length(channel2CellCoord));
%     figure;
%     imagesc(roisCropped{2,2}.CData);

    for cellNum = 1:length(channel2CellCoord)
        disp(sectionData.name);
        disp(channel2CellCoord{cellNum});
        
        %ch 2
        if xor(any(strcmp(FLIPPED_SECTIONS,sectionData.name)),  startsWith(sectionData.name,'RIGHT'))
            if length(channel2CellCoord{cellNum}) > 2 
                all_cell_coords{end + 1} = [600-channel2CellCoord{cellNum}(1,1), 600 - channel2CellCoord{cellNum}(1,2)];
            else
                all_cell_coords{end + 1} = [600-channel2CellCoord{cellNum}(1), 600 - channel2CellCoord{cellNum}(2)];
            end
        else
            if length(channel2CellCoord{cellNum}) > 2 
                all_cell_coords{end + 1} = [channel2CellCoord{cellNum}(1,1), 600-channel2CellCoord{cellNum}(1,2)];
            else
                all_cell_coords{end + 1} = [channel2CellCoord{cellNum}(1), 600 - channel2CellCoord{cellNum}(2)];
            end
        end
    end
        
        % ch3
    for cellNum = 1:length(channel3CellCoord)
        disp(sectionData.name);
        disp(channel3CellCoord{cellNum});
        
        if xor(any(strcmp(FLIPPED_SECTIONS,sectionData.name)),  startsWith(sectionData.name,'RIGHT'))
            if length(channel3CellCoord{cellNum}) > 2 
                all_cell_coords{end + 1} = [600-channel3CellCoord{cellNum}(1,1), 600 - channel3CellCoord{cellNum}(1,2)];
            else
                all_cell_coords{end + 1} = [600-channel3CellCoord{cellNum}(1), 600 - channel3CellCoord{cellNum}(2)];
            end
        else
            if length(channel3CellCoord{cellNum}) > 2 
                all_cell_coords{end + 1} = [channel3CellCoord{cellNum}(1,1),600- channel3CellCoord{cellNum}(1,2)];
            else
                all_cell_coords{end + 1} = [channel3CellCoord{cellNum}(1), 600 - channel3CellCoord{cellNum}(2)];
            end
        end
    end

end

figure;
mat = reshape(cell2matPETE(all_cell_coords),2,size(all_cell_coords,2));
scatter(mat(1,:),mat(2,:));
ylim([0, 600]);
xlim([0, 600]);

% d   = [d1;d2];
p   = gkde2(transpose(mat));

figure(2)
surf(p.x,p.y,p.pdf*1000000,'FaceAlpha','interp',...
    'AlphaDataMapping','scaled',...
    'AlphaData',p.pdf*1000000,...
    'FaceColor','interp','EdgeAlpha',0);
colormap(jet)
hold on

img    = imread(C);     %# Load a sample image
x      = size(p.x)
y      = size(p.y)
xImage = [0;600]%[p.x(1,1)+125 p.x(1,x(1))+125; p.x(x(1),1)+125 p.x(x(1),x(1))+125];   %# The x data for the image corners
yImage = [0;600]%[p.y(1,y(1))-10 p.y(1,1)-10; p.y(y(1),y(1))-10 p.y(y(1),1)-10];             %# The y data for the image corners
zImage = [0 0; 0 0];   %# The z data for the image corners
surf(xImage,yImage,zImage,...    %# Plot the surface
     'CData',img,...
     'FaceColor','texturemap'); %'CData',img,...
% set(gca,'Ydir','reverse')
view(2)
axis([p.x(1,1)+125 p.x(1,x(1))+125 p.y(1,y(1))-10 p.y(y(1),y(1))-10])
xlim([0,600]);
ylim([0,600]);

 %% Scatter Check Implementation
% figure(3)
% h1 = imread(C);
% image(h1)
% hold on
% scatter(eyeTrackData.animals.animals_01.subject_01.fixX,eyeTrackData.animals.animals_01.subject_01.fixY);
end
