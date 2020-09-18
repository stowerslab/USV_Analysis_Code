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


C = 'E:\Pete\transparent.png'; % image background

% uncomment the code for the channel you are heatmapping (in the for loop)

% Input:
% sectionFilePath = 'E:\Jingyi\Dropbox (Scripps Research)\Jingyi Data backup\Imaging\2019-4_ChR2_cFOS\BNST\COMBINED_female_male\sec4\'; % select directory with the .mat files
sectionFilePath = 'E:\Sourish\For Pete 20200728\control3\d-e\'
side = ''; % for title and saving purposes
% ----------
sectionDirectory = dir(fullfile(sectionFilePath,'*.mat'));


individual_plots = 1;
combined_plots = 0;

draw_boundary = 1; % yes or no, 1 or 0
boundary_num = 1; % fmu1a-c= 3... fmu1=1... fmu3d-e=1... 1, 2, or 3, get the boundary made for specific channel, just use the best one

ch1_ = 0;
ch2_ = 1;
ch3_ = 1;
% if sections are vertical, then instead of the default horizontal flip, we
% need to use a vertical flip.
vertical_flip = 1; % 1 for true if sections need to be flipped vertically

ch1_cell_coords = {};
ch2_cell_coords = {};
ch3_cell_coords = {};

FLIPPED_SECTIONS = {
    };

%    

% CHANNEL 1!
% counter = 1;
if (ch1_ == 1)
    for sectionData = sectionDirectory'
%         if (counter < startSlice) || (counter > endSlice)
%             disp(counter)
%             counter=counter+1;
%             continue;
%         end
        load(fullfile(sectionData.folder, sectionData.name)); 
        ROI_width=ROIwidth;
        ROI_length=ROIlength;
        disp(sectionData.name);

    %     length(channel2CellCoord)
        for cellNum = 1:length(channel1CellCoord)
    %        temporary for test
    %         disp(sectionData.name);
    %         disp(channel1CellCoord{cellNum});

    %         if ~contains(sectionData.name, 'RIGHT')
    %             continue;
    %         end

            if (vertical_flip == 0) && xor(any(strcmp(FLIPPED_SECTIONS,sectionData.name)),  startsWith(sectionData.name,'RIGHT'))
    %             THIS HANDLES HORIZONTAL FLIP
    
                if length(channel1CellCoord{cellNum}) > 2 
                    ch1_cell_coords{end + 1} = [ROI_width-channel1CellCoord{cellNum}(1,1), ROI_length-channel1CellCoord{cellNum}(1,2)];
                else
                    ch1_cell_coords{end + 1} = [ROI_width-channel1CellCoord{cellNum}(1), ROI_length-channel1CellCoord{cellNum}(2)];
                end

            elseif (vertical_flip == 1) % THIS HANDLES VERTICAL FLIP

                if contains(sectionData.name, 'X_FLIP') % if a horizontal flip is also necessary

                    if length(channel1CellCoord{cellNum}) > 2 && startsWith(sectionData.name,'RIGHT')
                        ch1_cell_coords{end + 1} = [ROI_width-channel1CellCoord{cellNum}(1,1), channel1CellCoord{cellNum}(1,2)];
                    elseif startsWith(sectionData.name,'RIGHT')
                        ch1_cell_coords{end + 1} = [ROI_width-channel1CellCoord{cellNum}(1), channel1CellCoord{cellNum}(2)];
                    else
                        ch1_cell_coords{end + 1} = [ROI_width-channel1CellCoord{cellNum}(1), ROI_length-channel1CellCoord{cellNum}(2)];
                    end

                else
                    if length(channel1CellCoord{cellNum}) > 2 && startsWith(sectionData.name,'RIGHT')
                        ch1_cell_coords{end + 1} = [channel1CellCoord{cellNum}(1,1), channel1CellCoord{cellNum}(1,2)];
                    elseif startsWith(sectionData.name,'RIGHT')
                        
                        ch1_cell_coords{end + 1} = [channel1CellCoord{cellNum}(1), channel1CellCoord{cellNum}(2)];
                    else
                        
                        ch1_cell_coords{end + 1} = [channel1CellCoord{cellNum}(1), ROI_length-channel1CellCoord{cellNum}(2)];
                    end

                end

            else % normal plot
                if length(channel1CellCoord{cellNum}) > 2 
                    ch1_cell_coords{end + 1} = [channel1CellCoord{cellNum}(1,1), ROI_length-channel1CellCoord{cellNum}(1,2)];
                else
                    ch1_cell_coords{end + 1} = [channel1CellCoord{cellNum}(1), ROI_length-channel1CellCoord{cellNum}(2)]; % channel1CellCoord{cellNum};
                end
            end
        end
        
        if (individual_plots == 1) && (length(ch1_cell_coords) > 0)
            
            INDIV_CELL_COORD= ch1_cell_coords(end-length(channel1CellCoord)+1:end);

            figure;
            mat = reshape(cell2matPETE(INDIV_CELL_COORD),2,size(INDIV_CELL_COORD,2));
            scatter(mat(1,:),mat(2,:),'MarkerEdgeAlpha', 0.2);
            ylim([0, ROI_length]);
            xlim([0, ROI_width]);
            daspect([1 1 1])

            % d   = [d1;d2];
            p   = gkde2(transpose(mat));
            make_title=[side, ' ch1 Scatter ', slice_num{1}]
            title(make_title)
            set(gca, 'Color', 'none'); % Sets axes background
            
            savefig([sectionFilePath, make_title]);
            figure;
%             figure(2)

            surf(p.x, p.y, p.pdf*1000000,'FaceAlpha','interp',...
                'AlphaDataMapping','scaled',...
                'AlphaData',p.pdf*1000000,...
                'FaceColor','interp','EdgeAlpha',0);

            % colormap(jet)
            hold on

            img    = imread(C);     %# Load a sample image
            x      = size(p.x);
            y      = size(p.y);
            xImage = [0;ROI_width];%[p.x(1,1)+125 p.x(1,x(1))+125; p.x(x(1),1)+125 p.x(x(1),x(1))+125];   %# The x data for the image corners
            yImage = [0;ROI_length];%[p.y(1,y(1))-10 p.y(1,1)-10; p.y(y(1),y(1))-10 p.y(y(1),1)-10];             %# The y data for the image corners
            zImage = [0 0; 0 0];   %# The z data for the image corners
            surf(xImage,yImage,zImage,...    %# Plot the surface
                 'CData',img,...
                 'FaceColor','texturemap'); %'CData',img,...
            % set(gca,'Ydir','reverse')
            view(2)
            axis([p.x(1,1)+125 p.x(1,x(1))+125 p.y(1,y(1))-10 p.y(y(1),y(1))-10])
            xlim([0,ROI_width]);
            ylim([0,ROI_length]);
            daspect([1 1 1])
            hold off
            make_title=[side, ' ch1 Heatmap ', slice_num{1}];
            title(make_title)
            savefig([sectionFilePath, make_title]);
%             counter = counter +1;
            close all;
        end
    end
    if combined_plots == 1


        figure;
        mat = reshape(cell2matPETE(ch1_cell_coords),2,size(ch1_cell_coords,2));
        scatter(mat(1,:),mat(2,:),'MarkerEdgeAlpha', 0.2);
        ylim([0, ROI_length]);
        xlim([0, ROI_width]);
        daspect([1 1 1])

        % d   = [d1;d2];
        p   = gkde2(transpose(mat));
        title('ch 1 combined scatter');
        set(gca, 'Color', 'none'); % Sets axes background
        savefig([sectionFilePath, 'ch1 combined']);

        % figure(2)
        figure;
        surf(p.x, p.y, p.pdf*1000000,'FaceAlpha','interp',... %1000000
            'AlphaDataMapping','scaled',...
            'AlphaData',p.pdf*1000000,...%1000000
            'FaceColor','interp','EdgeAlpha',0);

        % colormap(jet)
        hold on

        img    = imread(C);     %# Load a sample image
        x      = size(p.x)
        y      = size(p.y)
        xImage = [0;ROI_width]%[p.x(1,1)+125 p.x(1,x(1))+125; p.x(x(1),1)+125 p.x(x(1),x(1))+125];   %# The x data for the image corners
        yImage = [0;ROI_length]%[p.y(1,y(1))-10 p.y(1,1)-10; p.y(y(1),y(1))-10 p.y(y(1),1)-10];             %# The y data for the image corners
        zImage = [0 0; 0 0];   %# The z data for the image corners
        surf(xImage,yImage,zImage,...    %# Plot the surface
             'CData',img,...
             'FaceColor','texturemap'); %'CData',img,...
        % set(gca,'Ydir','reverse')
        view(2)
        axis([p.x(1,1)+125 p.x(1,x(1))+125 p.y(1,y(1))-10 p.y(y(1),y(1))-10])
        xlim([0,ROI_width]);
        ylim([0,ROI_length]);
        daspect([1 1 1])
        title('ch1 combined heatmap')
        set(gca, 'Color', 'none'); % Sets axes background
        savefig([sectionFilePath, 'ch1 heatmapcombined']);
        close all;
    end
end


if (ch2_==1)
    % CHANNEL 2!
%     counter = 1;
    for sectionData = sectionDirectory'
        load(fullfile(sectionData.folder, sectionData.name));
        
        disp(sectionData.name);
        slice_num = split(sectionData.name, '_');
        slice_num = split(slice_num(end), '.mat');
        slice_num = slice_num(1);
        
        ROI_width=ROIwidth;
        ROI_length=ROIlength;
    %     length(channel2CellCoord)
        
        for cellNum = 1:length(channel2CellCoord)
    %        temporary for test
    %         disp(sectionData.name);
    %         disp(channel2CellCoord{cellNum});

    %         if ~contains(sectionData.name, 'RIGHT')
    %             continue;
    %         end

            if (vertical_flip == 0) && xor(any(strcmp(FLIPPED_SECTIONS,sectionData.name)),  startsWith(sectionData.name,'RIGHT'))
    %             THIS HANDLES HORIZONTAL FLIP
                if length(channel2CellCoord{cellNum}) > 2 
                    ch2_cell_coords{end + 1} = [ROI_width-channel2CellCoord{cellNum}(1,1), ROI_length-channel2CellCoord{cellNum}(1,2)];
                else
                    ch2_cell_coords{end + 1} = [ROI_width-channel2CellCoord{cellNum}(1), ROI_length-channel2CellCoord{cellNum}(2)];
                end

            elseif (vertical_flip == 1) % THIS HANDLES VERTICAL FLIP

                if contains(sectionData.name, 'X_FLIP') % if a horizontal flip is also necessary

                    if length(channel2CellCoord{cellNum}) > 2 && startsWith(sectionData.name,'RIGHT')
                        ch2_cell_coords{end + 1} = [ROI_width-channel2CellCoord{cellNum}(1,1), channel2CellCoord{cellNum}(1,2)];
                    elseif startsWith(sectionData.name,'RIGHT')
                        ch2_cell_coords{end + 1} = [ROI_width-channel2CellCoord{cellNum}(1), channel2CellCoord{cellNum}(2)];
                    else
                        ch2_cell_coords{end + 1} = [ROI_width-channel2CellCoord{cellNum}(1), ROI_length-channel2CellCoord{cellNum}(2)];
                    end

                else

                    if length(channel2CellCoord{cellNum}) > 2 && startsWith(sectionData.name,'RIGHT')
                        ch2_cell_coords{end + 1} = [channel2CellCoord{cellNum}(1,1), channel2CellCoord{cellNum}(1,2)];
                    elseif startsWith(sectionData.name,'RIGHT')
                        ch2_cell_coords{end + 1} = [channel2CellCoord{cellNum}(1), channel2CellCoord{cellNum}(2)];
                    else
                        ch2_cell_coords{end + 1} = [channel2CellCoord{cellNum}(1), ROI_length-channel2CellCoord{cellNum}(2)];
                    end

                end

            else % normal plot
                if length(channel2CellCoord{cellNum}) > 2 
                    ch2_cell_coords{end + 1} = [channel2CellCoord{cellNum}(1,1), ROI_length-channel2CellCoord{cellNum}(1,2)];
                else
                    ch2_cell_coords{end + 1} = [channel2CellCoord{cellNum}(1), ROI_length-channel2CellCoord{cellNum}(2)]; % channel2CellCoord{cellNum};
                end
            end
        end
        
        if (individual_plots == 1) && (length(ch2_cell_coords) > 0)
            scatter_title_for_tif = ['scatter combined ch +boundary ', slice_num{1}, '.tif'];
            heatmap_title_for_tif = ['heatmap combined ch +boundary ', slice_num{1}, '.tif'];
            
            INDIV_CELL_COORD= ch2_cell_coords(end-length(channel2CellCoord)+1:end);

            figure;
            mat = reshape(cell2matPETE(INDIV_CELL_COORD),2,size(INDIV_CELL_COORD,2));
            scatter(mat(1,:),mat(2,:),'MarkerEdgeAlpha', 0.2);
            ylim([0, ROI_length]);
            xlim([0, ROI_width]);
            daspect([1 1 1])

            % d   = [d1;d2];
            p   = gkde2(transpose(mat));
            make_title=[side, ' ch2 Scatter ', slice_num{1}]
            title(make_title)
%             savefig([sectionFilePath, make_title]);
            export_fig([sectionFilePath, scatter_title_for_tif], '-transparent', '-append');
            figure;
    %             figure(2)

            surf(p.x, p.y, p.pdf*1000000,'FaceAlpha','interp',...
                'AlphaDataMapping','scaled',...
                'AlphaData',p.pdf*1000000,...
                'FaceColor','interp','EdgeAlpha',0);

            % colormap(jet)
            hold on

            img    = imread(C);     %# Load a sample image
            x      = size(p.x)
            y      = size(p.y)
            xImage = [0;ROI_width]%[p.x(1,1)+125 p.x(1,x(1))+125; p.x(x(1),1)+125 p.x(x(1),x(1))+125];   %# The x data for the image corners
            yImage = [0;ROI_length]%[p.y(1,y(1))-10 p.y(1,1)-10; p.y(y(1),y(1))-10 p.y(y(1),1)-10];             %# The y data for the image corners
            zImage = [0 0; 0 0];   %# The z data for the image corners
            surf(xImage,yImage,zImage,...    %# Plot the surface
                 'CData',img,...
                 'FaceColor','texturemap'); %'CData',img,...
            % set(gca,'Ydir','reverse')
            view(2)
            axis([p.x(1,1)+125 p.x(1,x(1))+125 p.y(1,y(1))-10 p.y(y(1),y(1))-10])
            xlim([0,ROI_width]);
            ylim([0,ROI_length]);
            daspect([1 1 1])
            hold off
            make_title=[side, ' ch2 Heatmap ', slice_num{1}]
            title(make_title)
%             savefig([sectionFilePath, make_title]);
            export_fig([sectionFilePath, heatmap_title_for_tif], '-transparent', '-append');
%             counter = counter +1;
            
%             if draw_boundary == 1
%                 btest = boundaries{boundary_num};
%                 plot(btest(:,2) - roiCoord(1), ROIlength - btest(:,1) + roiCoord(2), 'b', 'LineWidth', 1)
%                 
%                 export_fig([sectionFilePath, scatter_title_for_tif], '-transparent', '-append');
%                 export_fig([sectionFilePath, heatmap_title_for_tif], '-transparent', '-append');
%             end
            close all;
        end

    end
    if combined_plots == 1
        figure;
        mat = reshape(cell2matPETE(ch2_cell_coords),2,size(ch2_cell_coords,2));
        scatter(mat(1,:),mat(2,:),'MarkerEdgeAlpha', 0.2);
        ylim([0, ROI_length]);
        xlim([0, ROI_width]);
        daspect([1 1 1])

        % d   = [d1;d2];
        p   = gkde2(transpose(mat));
        title('ch2 combined scatter');
        set(gca, 'Color', 'none'); % Sets axes background
        savefig([sectionFilePath, 'ch2 combined']);

        % figure(2)
        figure;


        surf(p.x,p.y,p.pdf*1000000,'FaceAlpha','interp',...
            'AlphaDataMapping','scaled',...
            'AlphaData',p.pdf*1000000,...
            'FaceColor','interp','EdgeAlpha',0);
    %     surf(p.x, p.y, p.pdf*1000000,'FaceAlpha','interp',... %1000000
    %         'AlphaDataMapping','scaled',...
    %         'AlphaData',p.pdf*1000000,... % 1000000
    %         'FaceColor','interp','EdgeAlpha',0);

    %     colormap(jet)
        hold on

        img    = imread(C);     %# Load a sample image
        x      = size(p.x)
        y      = size(p.y)
        xImage = [0;ROI_width]%[p.x(1,1)+125 p.x(1,x(1))+125; p.x(x(1),1)+125 p.x(x(1),x(1))+125];   %# The x data for the image corners
        yImage = [0;ROI_length]%[p.y(1,y(1))-10 p.y(1,1)-10; p.y(y(1),y(1))-10 p.y(y(1),1)-10];             %# The y data for the image corners
        zImage = [0 0; 0 0];   %# The z data for the image corners

        % comment this surf() out to have color range
        surf(xImage,yImage,zImage,...    %# Plot the surface
             'CData',img,...
             'FaceColor','texturemap'); %'CData',img,...
    %     set(gca,'Ydir','reverse')

        view(2)
        axis([p.x(1,1)+125 p.x(1,x(1))+125 p.y(1,y(1))-10 p.y(y(1),y(1))-10])
        xlim([0,ROI_width]);
        ylim([0,ROI_length]);
        daspect([1 1 1])
        title('ch2 combined heatmap')
        set(gca, 'Color', 'none'); % Sets axes background
        savefig([sectionFilePath, 'ch2 heatmapcombined']);
        close all;
    end
end

if (ch3_==1)
    % CHANNEL 3!
%     counter = 1;
    for sectionData = sectionDirectory'
        load(fullfile(sectionData.folder, sectionData.name)); 
        disp(sectionData.name);
        
        slice_num = split(sectionData.name, '_');
        slice_num = split(slice_num(end), '.mat');
        slice_num = slice_num(1);
        
        ROI_width=ROIwidth;
        ROI_length=ROIlength;
    %     length(channel2CellCoord)
        
        for cellNum = 1:length(channel3CellCoord)
    %        temporary for test
    %         disp(sectionData.name);
    %         disp(channel2CellCoord{cellNum});

    %         if ~contains(sectionData.name, 'RIGHT')
    %             continue;
    %         end

            if (vertical_flip == 0) && xor(any(strcmp(FLIPPED_SECTIONS,sectionData.name)),  startsWith(sectionData.name,'RIGHT'))
    %             THIS HANDLES HORIZONTAL FLIP
                if length(channel3CellCoord{cellNum}) > 2 
                    ch3_cell_coords{end + 1} = [ROI_width-channel3CellCoord{cellNum}(1,1), ROI_length-channel3CellCoord{cellNum}(1,2)];
                else
                    ch3_cell_coords{end + 1} = [ROI_width-channel3CellCoord{cellNum}(1), ROI_length-channel3CellCoord{cellNum}(2)];
                end

            elseif (vertical_flip == 1) % THIS HANDLES VERTICAL FLIP

                if contains(sectionData.name, 'X_FLIP') % if a horizontal flip is also necessary

                    if length(channel3CellCoord{cellNum}) > 2 && startsWith(sectionData.name,'RIGHT')
                        ch3_cell_coords{end + 1} = [ROI_width-channel3CellCoord{cellNum}(1,1), channel3CellCoord{cellNum}(1,2)];
                    elseif startsWith(sectionData.name,'RIGHT')
                        ch3_cell_coords{end + 1} = [ROI_width-channel3CellCoord{cellNum}(1), channel3CellCoord{cellNum}(2)];
                    else
                        ch3_cell_coords{end + 1} = [ROI_width-channel3CellCoord{cellNum}(1), ROI_length-channel3CellCoord{cellNum}(2)];
                    end

                else

                    if length(channel3CellCoord{cellNum}) > 2 && startsWith(sectionData.name,'RIGHT')
                        ch3_cell_coords{end + 1} = [channel3CellCoord{cellNum}(1,1), channel3CellCoord{cellNum}(1,2)];
                    elseif startsWith(sectionData.name,'RIGHT')
                        ch3_cell_coords{end + 1} = [channel3CellCoord{cellNum}(1), channel3CellCoord{cellNum}(2)];
                    else
                        ch3_cell_coords{end + 1} = [channel3CellCoord{cellNum}(1), ROI_length-channel3CellCoord{cellNum}(2)];
                    end

                end

            else % normal plot
                if length(channel3CellCoord{cellNum}) > 2 
                    ch3_cell_coords{end + 1} = [channel3CellCoord{cellNum}(1,1), ROI_length-channel3CellCoord{cellNum}(1,2)];
                else
                    ch3_cell_coords{end + 1} = [channel3CellCoord{cellNum}(1), ROI_length-channel3CellCoord{cellNum}(2)]; % channel3CellCoord{cellNum};
                end
            end
        end
        
        if (individual_plots == 1) && (length(ch3_cell_coords) > 0)
            
            scatter_title_for_tif = ['scatter combined ch +boundary ', slice_num{1}, '.tif'];
            heatmap_title_for_tif = ['heatmap combined ch +boundary ', slice_num{1}, '.tif'];
            
            INDIV_CELL_COORD= ch3_cell_coords(end-length(channel3CellCoord)+1:end);

            figure;
            mat = reshape(cell2matPETE(INDIV_CELL_COORD),2,size(INDIV_CELL_COORD,2));
            scatter(mat(1,:),mat(2,:),'MarkerEdgeAlpha', 0.2);
            ylim([0, ROI_length]);
            xlim([0, ROI_width]);
            daspect([1 1 1])

            % d   = [d1;d2];
            p   = gkde2(transpose(mat));
            make_title=[side, ' ch3 Scatter ', slice_num{1}]
            title(make_title)
%             savefig([sectionFilePath, make_title]);
            export_fig([sectionFilePath, scatter_title_for_tif], '-transparent', '-append');
            
            figure;
    %             figure(2)

            surf(p.x, p.y, p.pdf*1000000,'FaceAlpha','interp',...
                'AlphaDataMapping','scaled',...
                'AlphaData',p.pdf*1000000,...
                'FaceColor','interp','EdgeAlpha',0);

            % colormap(jet)
            hold on

            img    = imread(C);     %# Load a sample image
            x      = size(p.x);
            y      = size(p.y);
            xImage = [0;ROI_width]%[p.x(1,1)+125 p.x(1,x(1))+125; p.x(x(1),1)+125 p.x(x(1),x(1))+125];   %# The x data for the image corners
            yImage = [0;ROI_length]%[p.y(1,y(1))-10 p.y(1,1)-10; p.y(y(1),y(1))-10 p.y(y(1),1)-10];             %# The y data for the image corners
            zImage = [0 0; 0 0];   %# The z data for the image corners
            surf(xImage,yImage,zImage,...    %# Plot the surface
                 'CData',img,...
                 'FaceColor','texturemap'); %'CData',img,...
            % set(gca,'Ydir','reverse')
            view(2)
            axis([p.x(1,1)+125 p.x(1,x(1))+125 p.y(1,y(1))-10 p.y(y(1),y(1))-10])
            xlim([0,ROI_width]);
            ylim([0,ROI_length]);
            daspect([1 1 1])
            hold off
            make_title=[side, ' ch3 Heatmap ', slice_num{1}]
            title(make_title)
%             savefig([sectionFilePath, make_title]);

            export_fig([sectionFilePath, heatmap_title_for_tif], '-transparent', '-append');
%             counter = counter +1;
            
            if draw_boundary == 1
                bounds_1ch = boundaries{boundary_num};
                clearvars boundaries;
                figure;
                hold on;
                
                for k = 1:length(bounds_1ch)
                   if iscell(bounds_1ch)
                        btest = bounds_1ch{k};
                   else
                       btest = bounds_1ch;
                       plot(btest(:,2) - roiCoord(1), ROIlength - btest(:,1) + roiCoord(2), 'b', 'LineWidth', 1)
                       break;
                   end
                   
                   plot(btest(:,2) - roiCoord(1), ROIlength - btest(:,1) + roiCoord(2), 'b', 'LineWidth', 1)
                   clearvars btest;
                end
                
                title('brain segment rough boundary');
%                 plot(btest(:,2) - roiCoord(1), ROIlength - btest(:,1) + roiCoord(2), 'b', 'LineWidth', 1)
                xlim([0,ROI_width]);
                ylim([0,ROI_length]);
                export_fig([sectionFilePath, scatter_title_for_tif], '-transparent', '-append');
                export_fig([sectionFilePath, heatmap_title_for_tif], '-transparent', '-append');
                
                figure;
                imagesc(channelArray{2});
                export_fig([sectionFilePath, scatter_title_for_tif], '-transparent', '-append');
                export_fig([sectionFilePath, heatmap_title_for_tif], '-transparent', '-append');
            end
            close all;
        end

    end
    if combined_plots == 1
        figure;
        mat = reshape(cell2matPETE(ch3_cell_coords),2,size(ch3_cell_coords,2));
        scatter(mat(1,:),mat(2,:),'MarkerEdgeAlpha', 0.2);
        ylim([0, ROI_length]);
        xlim([0, ROI_width]);
        daspect([1 1 1])

        % d   = [d1;d2];
        p   = gkde2(transpose(mat));

        hold on;

        title('ch3 combined scatter');
        set(gca, 'Color', 'none'); % Sets axes background
        savefig([sectionFilePath, 'ch3 combined']);

        % figure(2)
        figure;


        surf(p.x,p.y,p.pdf*1000000,'FaceAlpha','interp',...
            'AlphaDataMapping','scaled',...
            'AlphaData',p.pdf*1000000,...
            'FaceColor','interp','EdgeAlpha',0);
    %     surf(p.x, p.y, p.pdf*1000000,'FaceAlpha','interp',... %1000000
    %         'AlphaDataMapping','scaled',...
    %         'AlphaData',p.pdf*1000000,... % 1000000
    %         'FaceColor','interp','EdgeAlpha',0);

    %     colormap(jet)
        hold on

        img    = imread(C);     %# Load a sample image
        x      = size(p.x)
        y      = size(p.y)
        xImage = [0;ROI_width]%[p.x(1,1)+125 p.x(1,x(1))+125; p.x(x(1),1)+125 p.x(x(1),x(1))+125];   %# The x data for the image corners
        yImage = [0;ROI_length]%[p.y(1,y(1))-10 p.y(1,1)-10; p.y(y(1),y(1))-10 p.y(y(1),1)-10];             %# The y data for the image corners
        zImage = [0 0; 0 0];   %# The z data for the image corners

        % comment this surf() out to have color range
        surf(xImage,yImage,zImage,...    %# Plot the surface
             'CData',img,...
             'FaceColor','texturemap'); %'CData',img,...
    %     set(gca,'Ydir','reverse')

        view(2)
        axis([p.x(1,1)+125 p.x(1,x(1))+125 p.y(1,y(1))-10 p.y(y(1),y(1))-10])
        xlim([0,ROI_width]);
        ylim([0,ROI_length]);
        daspect([1 1 1])



        title('ch3 combined heatmap')
        set(gca, 'Color', 'none'); % Sets axes background
        savefig([sectionFilePath, 'ch3 heatmapcombined']);
        close all;
        
    end
end
