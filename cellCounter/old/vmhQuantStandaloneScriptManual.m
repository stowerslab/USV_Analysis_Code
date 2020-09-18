%% vmhQuantStandaloneScriptManual: top-level script to read image and assist counting viral expression in the VMH
%   Jason Keller, 12th Mar 2013

% inputs:
%       stiched image (from Nikon C2) in ND2 format
%       manual ellipse-based VMHvl ROI patch

% outputs:
%       in VMHvl, total cell counts for each image color channel

clear all; close all;

%% image read & format
% alternatively read nd2 directly to cut out the ElementsViewer step (and also to access
% image metadata)
% this may be possible from MATLAB Central:http://www.mathworks.com/matlabcentral/fileexchange/32920-imread-for-multiple-life-science-image-file-formats 
% BUT, Java heap size needs to be increased in MATLAB Preferences

filepath = 'C:\data\jason\peterVirusData\';
% filepath = 'C:\Users\keller\School\UCSD\stowersLab\code\cellSegmentation\';
% mouseNum = '19249';
% sectionNum = '8';
mouseNum = '70324';
sectionNum = '2';  % this one has distortion and is very dim
% mouseNum = '19485';
% sectionNum = '4';  % this one is very bright
filepathTotal = [filepath, mouseNum, '\', sectionNum, '.tif']  %display full path

% Note that these are actually 12-bit images, so MATLAB will represent with
% uint16; however, the highest value (saturation) will be 2^12 - 1 = 4095,
% so we can normalize by this to make a standard colormap

% possibly use imadjust/histeq/adapthisteq to fill 16 bits to optimize contrast at the 16-bit
% level, but this can also warp image intesities significantly
max16bitValue = 2^16-1;
% Ib = imadjust(imread(filepathTotal, 1)); %read 1st image in TIF file - Hoechst
% Ig = imadjust(imread(filepathTotal, 2)); %read 2nd image in TIF file - GFP
% Ir = imadjust(imread(filepathTotal, 3)); %read 3rd image in TIF file - TDT
Ib = imread(filepathTotal, 1); %read 1st image in TIF file - Hoechst
Ig = imread(filepathTotal, 2); %read 2nd image in TIF file - GFP
Ir = imread(filepathTotal, 3); %read 3rd image in TIF file - TDT


%% ROI selection
% manually perform selection of 3 ROIs using the Hoechst channel: (VMHc can
% be computed from these if necessary)
numRoi = 3;
roiPhrases = {'(1) Select total VMH ROI:'; '(2) Select VMH_dm ROI:'; '(3) Select VMH_vl ROI:'};
figure; hIm = imagesc(Ib); %create image handle and show image
hAxes = gca;  % get handle to axes for 'waitfor' funtion below
colormap(gray); %set to grayscale
% set(gcf, 'Position', get(0,'Screensize')); % Maximize figure - commented since this can stretch figure if not at same aspect ratio as monitor
zoom on
waitfor(hAxes,'CameraPosition');  % wait until zoom changes - note that you only have one chance to change ;-)
zoom off

% TO DO: add option to zoom in here first

cmap =  [0         0    1.0000;...
         0    0.5000         0;...
    1.0000         0         0;...
         0    0.7500    0.7500]; %colormap to ID different ROIs, first four lines from colormap(lines)
     
roiPatchArray = struct(   ...    %structure array containing data for each ROI selected in the image        
                                'hPatch',{},...          %Handles to ROI patch objects
                                'patchPosition',{},...   %Positions of patches, specified in pixel coordinates
                                'patchMask',{},...       %binary masks for each ROI patch (pixels inside patch = 1)                                                                                
                                'patchColor',{});        %color value of patch as index into 'lines' colomap
                            
for i=1:numRoi
    set(gcf, 'Name',roiPhrases{i});
    % User selects ROI with mouse:
    hRoi = imfreehand(gca);
    hMask = createMask(hRoi,hIm);
    pos = hRoi.getPosition;
    color = cmap(i,:);
    hPatch = patch('Parent',gca,'Visible','on',...  %draw colored patch to visualize ROI
                   'XData',pos(:,1),'YData',pos(:,2),...
                   'EdgeColor',color,'LineWidth',2,'FaceColor','none');%,...
%                  'FaceAlpha',0.15,'FaceColor',color);  %NOTE: something about specifying FaceAlpha really slows down computer, graphics rendering                

    %save ROIs:   
    roiPatchArray(end+1).hPatch = hPatch;
    roiPatchArray(end).patchPosition = pos;
    roiPatchArray(end).patchMask = hMask;  %access ex: imshow(roiPatchArray(1,1).patchMask)
    roiPatchArray(end).patchColor = i; %store color for later use
%     delete(hRoi); %delete imfreehand object
end
%close figure
% close(gcf);
set(gcf, 'Name','all ROIs');

%create and view RGB image to test:
% Irgb = zeros(size(Ib,1),size(Ib,2),3);
% Ir(~roiPatchArray(1,1).patchMask) = 0; Irgb(:,:,1) = double(Ir)./  double(max(max(Ir))); 
% Ig(~roiPatchArray(1,1).patchMask) = 0; Irgb(:,:,2) = double(Ig)./  double(max(max(Ig)));
% Ib(~roiPatchArray(1,1).patchMask) = 0; Irgb(:,:,3) = double(Ib)./  double(max(max(Ib)));
% imshow(Irgb);

%% Manual cell marking:
% first create image with only nuclei in ROI:
Ibz = Ib; Ibz(~roiPatchArray(1,1).patchMask) = 0;
hf = figure; hImIbz = imagesc(Ibz); colormap(gray);
set(hf, 'Name','zoom in, mouse click on cells, then press f when finished');
hAxesIbz = gca;
zoom on
waitfor(hAxesIbz,'CameraPosition');  % wait until zoom changes - note that you only have one chance to change ;-)
zoom off

exit = false; %flag to wait for button press
i = 1; %index to cells counted

while(~exit)
    hPb(i) = impoint(gca); %#ok<*SAGROW>  %perhaps better to use ginput for this
%     setString(hPb(i),num2str(i))
    manualPointsIb(i,:) = getPosition(hPb(i));  % to plot later use:  plot(manualPointsIb(:,1), manualPointsIb(:,2), 'y', 'LineStyle', 'none', 'Marker', '+')
    
    % set DM points:
    if roiPatchArray(1,2).patchMask(round(manualPointsIb(i,2)),round(manualPointsIb(i,1))) % check if point is within DM (must round off)
        manualPointsIbDm(i,:) = manualPointsIb(i,:);
    else
        manualPointsIbDm(i,:) = [0 0];  %set to zero if not in DM
    end
    % set VL points:
    if roiPatchArray(1,3).patchMask(round(manualPointsIb(i,2)),round(manualPointsIb(i,1))) % check if point is within VL (must round off)
        manualPointsIbVl(i,:) = manualPointsIb(i,:);
    else
        manualPointsIbVl(i,:) = [0 0];
    end

    setColor(hPb(i),'y')
    i = i + 1;
    currButton = get(hf,'CurrentCharacter');  %press f key while window active when done
    if isempty(currButton)
        exit = false;
    else
        exit = currButton == 'f';
    end
end

% NOW count GFP cells:
Igz = Ig; Igz(~roiPatchArray(1,1).patchMask) = 0;
hfg = figure; hImIgz = imagesc(Igz); colormap(gray);
set(hfg, 'Name','Gfp: zoom in, mouse click on cells, then press f when finished');
hAxesIgz = gca;
zoom on
waitfor(hAxesIgz,'CameraPosition');  % wait until zoom changes - note that you only have one chance to change ;-)
zoom off

exit = false; %flag to wait for button press
i = 1; %index to cells counted

while(~exit)
    hPg(i) = impoint(gca); %#ok<*SAGROW>  %perhaps better to use ginput for this
%     setString(hPg(i),num2str(i))
    manualPointsIg(i,:) = getPosition(hPg(i));  % to plot later use:  plot(manualPointsIg(:,1), manualPointsIg(:,2), 'y', 'LineStyle', 'none', 'Marker', '+')
    
    % set DM points:
    if roiPatchArray(1,2).patchMask(round(manualPointsIg(i,2)),round(manualPointsIg(i,1))) % check if point is within DM (must round off)
        manualPointsIgDm(i,:) = manualPointsIg(i,:);
    else
        manualPointsIgDm(i,:) = [0 0];  %set to zero if not in DM
    end
    % set VL points:
    if roiPatchArray(1,3).patchMask(round(manualPointsIg(i,2)),round(manualPointsIg(i,1))) % check if point is within VL (must round off)
        manualPointsIgVl(i,:) = manualPointsIg(i,:);
    else
        manualPointsIgVl(i,:) = [0 0];
    end

    setColor(hPg(i),'y')
    i = i + 1;
    currButton = get(hfg,'CurrentCharacter');  %press f key while window active when done
    if isempty(currButton)
        exit = false;
    else
        exit = currButton == 'f';
    end
end

% NOW count TDT cells:
Irz = Ir; Irz(~roiPatchArray(1,1).patchMask) = 0;
hfr = figure; hImIrz = imagesc(Irz); colormap(gray);
set(hfr, 'Name','Tdt: zoom in, mouse click on cells, then press f when finished');
hAxesIrz = gca;
zoom on
waitfor(hAxesIrz,'CameraPosition');  % wait until zoom changes - note that you only have one chance to change ;-)
zoom off

exit = false; %flag to wait for button press
i = 1; %index to cells counted

while(~exit)
    hPr(i) = impoint(gca); %#ok<*SAGROW>  %perhaps better to use ginput for this
%     setString(hPg(i),num2str(i))
    manualPointsIr(i,:) = getPosition(hPr(i));  % to plot later use:  plot(manualPointsIr(:,1), manualPointsIr(:,2), 'y', 'LineStyle', 'none', 'Marker', '+')
    
    % set DM points:
    if roiPatchArray(1,2).patchMask(round(manualPointsIr(i,2)),round(manualPointsIr(i,1))) % check if point is within DM (must round off)
        manualPointsIrDm(i,:) = manualPointsIr(i,:);
    else
        manualPointsIrDm(i,:) = [0 0];  %set to zero if not in DM
    end
    % set VL points:
    if roiPatchArray(1,3).patchMask(round(manualPointsIr(i,2)),round(manualPointsIr(i,1))) % check if point is within VL (must round off)
        manualPointsIrVl(i,:) = manualPointsIr(i,:);
    else
        manualPointsIrVl(i,:) = [0 0];
    end

    setColor(hPr(i),'y')
    i = i + 1;
    currButton = get(hfr,'CurrentCharacter');  %press f key while window active when done
    if isempty(currButton)
        exit = false;
    else
        exit = currButton == 'f';
    end
end



%% Cell counting
manualRoiPatchArray = roiPatchArray;
% manualPointsIr = hPr;
% manualPointsIg = hPg;

manualNumNucleiTotal = length(hPb);
manualNumNucleiDm = length(find(manualPointsIbDm(:,1) ~= 0 & manualPointsIbDm(:,2) ~= 0));
manualNumNucleiVl = length(find(manualPointsIbVl(:,1) ~= 0 & manualPointsIbVl(:,2) ~= 0));

manualNumGfpTotal = length(hPg);
manualNumGfpDm = length(find(manualPointsIgDm(:,1) ~= 0 & manualPointsIgDm(:,2) ~= 0));
manualNumGfpVl = length(find(manualPointsIgVl(:,1) ~= 0 & manualPointsIgVl(:,2) ~= 0));

manualNumTdtTotal = length(hPr);
manualNumTdtDm = length(find(manualPointsIrDm(:,1) ~= 0 & manualPointsIrDm(:,2) ~= 0));
manualNumTdtVl = length(find(manualPointsIrVl(:,1) ~= 0 & manualPointsIrVl(:,2) ~= 0));

%% save variables
matName = [filepath, mouseNum, '\', mouseNum, '_', sectionNum, '_manualCount.mat'];
save(matName, '-regexp', '^manual');  %save all variables starting with 'manual'

    