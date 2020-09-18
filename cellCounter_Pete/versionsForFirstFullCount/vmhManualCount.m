function [numRed, numGreen, numBlue] = vmhManualCount(IrSub, IgSub, IbSub)
%% vmhManualCount: script to read image and assist counting viral expression in the VMH
%  Cell marking process: 
%
%   Jason Keller, 12th Mar 2013
% inputs:
%       filepath to MAT file with stiched image (from Nikon C2) 3-matrix format (one each for R/G/B channels) - these are already cropped
%       manual ellipse-based VMHvl ROI patch
% outputs:
%       in VMHvl, total cell counts for each image color channel

% TO DO: 
%       make done buttons!
%       auto full screen
%       make separate programs with independent save files for each channel
%
%       really want to use imagesc & jet colormap?
%       want to threshold first to make counting easier?

%% Manual cell marking; 
% first create image with only nuclei in ROI:
hf = figure; hImIbz = imagesc(IbSub); colormap(jet);  %
set(gcf, 'Position', get(0,'Screensize')); % maximize figure
set(hf, 'Name','Left-mouse-click on Nissl cells, then press DONE when finished');
exit = false; %flag to wait for button press
i = 1; %index to cells counted

hDoneControl = uicontrol(hFig, 'Style', 'pushbutton', 'String', 'Done Rotating',...  %make button to press when done manipulating
                          'Callback', @doneCallback, 'Position', [20 20 60 20]);      

%% first mark Nissl bodies
while(~exit)
    hPb(i) = impoint(gca); %#ok<*SAGROW>  %perhaps better to use ginput for this
%     setString(hPb(i),num2str(i))
    manualPointsIb(i,:) = getPosition(hPb(i));  % to plot later use:  plot(manualPointsIb(:,1), manualPointsIb(:,2), 'y', 'LineStyle', 'none', 'Marker', '+')
    
    setColor(hPb(i),'w')
    if exit && i==1 %if DONE is pressed before counting any cells, this means no cells counted
        hPb = [];
    else
        i = i + 1;
    end
end

%% NOW mark GFP cells:
hfg = figure; hImIgz = imagesc(IgSub); colormap(jet);
set(hfg, 'Name','Gfp: mouse click on cells, then press f when finished');
exit = false; %flag to wait for button press
i = 1; %index to cells counted

while(~exit)
    hPg(i) = impoint(gca); %#ok<*SAGROW>  %perhaps better to use ginput for this
%     setString(hPg(i),num2str(i))
    manualPointsIg(i,:) = getPosition(hPg(i));  % to plot later use:  plot(manualPointsIg(:,1), manualPointsIg(:,2), 'y', 'LineStyle', 'none', 'Marker', '+')

    setColor(hPg(i),'w')
    i = i + 1;
    currButton = get(hfg,'CurrentCharacter');  %press f key while window active when done
    if isempty(currButton)
        exit = false;
    else
        exit = (currButton == 'f');
    end
end

%% NOW mark TDT cells:
hfr = figure; hImIrz = imagesc(IrSub); colormap(jet);
set(hfr, 'Name','Tdt: mouse click on cells, then press f when finished');
exit = false; %flag to wait for button press
i = 1; %index to cells counted

while(~exit)
    hPr(i) = impoint(gca); %#ok<*SAGROW>  %perhaps better to use ginput for this
%     setString(hPg(i),num2str(i))
    manualPointsIr(i,:) = getPosition(hPr(i));  % to plot later use:  plot(manualPointsIr(:,1), manualPointsIr(:,2), 'y', 'LineStyle', 'none', 'Marker', '+')
    
    setColor(hPr(i),'w')
    i = i + 1;
    currButton = get(hfr,'CurrentCharacter');  %press f key while window active when done
    if isempty(currButton)
        exit = false;
    else
        exit = (currButton == 'f');
    end
end

%% Cell counting
numBlue = length(hPb);
numGreen = length(hPg);
numRed = length(hPr);


%% Callbacks %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
                           
function doneCallback(hObject, ~)
    delete(hObject);
    exit = true;
end

end %end of function

