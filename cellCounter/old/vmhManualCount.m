function [numCells, avgFluo, manualPoints] = vmhManualCount(ISub, hVmhVlPatchMask)
%% vmhManualCount: script to read image and assist counting viral expression in the VMH
%  Cell marking process: click with left mouse button
%   Jason Keller, 12th Mar 2013

% inputs:
%       image matrix saved to MAT file, already masked & cropped for the VMHvl
% outputs:
%       in VMHvl, total cell count and avg. fluo. intensity

% TO DO: 
%       threshold or smooth first to make counting easier?

%% Manual cell marking; 
% first create image with only nuclei in ROI:
hf = figure; imagesc(ISub); colormap(gray);  %use 'imagesc' to spread out image histogram possibly 
set(gcf, 'Position', get(0,'Screensize')); % maximize figure
axis off;
set(hf, 'Name','Left-mouse-click on cells, then press DONE when finished');
exit = false; %flag to wait for button press
i = 1; %index to cells counted

hDoneControl = uicontrol(hf, 'Style', 'pushbutton', 'String', 'Done Counting',...  %make button to press when done manipulating
                          'Callback', @doneCallback, 'Position', [50 50 100 20]);       %#ok<NASGU>

%% first mark Nissl bodies
while(~exit)
    hPoints(i) = impoint(gca); %#ok<AGROW,*SAGROW>  %perhaps better to use ginput for this
%     setString(hPoints(i),num2str(i))
    manualPoints(i,:) = getPosition(hPoints(i));  %#ok<AGROW> % to plot later use:  plot(manualPointsIb(:,1), manualPointsIb(:,2), 'y', 'LineStyle', 'none', 'Marker', '+')
    
    setColor(hPoints(i),'r')  %use red on gray background
    if exit %if DONE is pressed, delete last point after pressing DONE
        hPoints(end) = [];
        manualPoints(end,:) = [];
    else
        i = i + 1;
    end
end

% count cells:
numCells = length(hPoints);

% calculate average fluo intensity within the mask:
avgFluo = mean(ISub(hVmhVlPatchMask));

%% Callbacks %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
                           
function doneCallback(hObject, ~)
    delete(hObject);
    exit = true;
end

end %end of function

