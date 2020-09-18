function [numCells, manualPoints] = manualCount(ISub)
%% manualCount: script to read image and assist counting viral expression within a ROI
%  Cell marking process: click with left mouse button
%   Jason Keller, 8th Feb 2016

% inputs:
%       image alread masked with ROI
% outputs:
%       in ROI, total cell count and markers

% TO DO: 
%       threshold or smooth first to make counting easier?

%% Manual cell marking; 
% first create image with only selected channel ROI:
hf = figure; imagesc(ISub); colormap(gray);  %use 'imagesc' to spread out image histogram 
% set(gcf, 'Position', get(0,'Screensize')); % maximize figure - careful as this can distort cell shapes
axis off;
set(hf, 'Name','Left-mouse-click on cells, then press DONE when finished, then click image');
exit = false; %flag to wait for button press
i = 1; %index to cells counted

hDoneControl = uicontrol(hf, 'Style', 'pushbutton', 'String', 'Done Counting',...  %make button to press when done manipulating
                          'Callback', @doneCallback, 'Position', [50 50 100 20]);       %#ok<NASGU>
while(~exit)
    hPoints(i) = impoint(gca); %#ok<AGROW,*SAGROW>  %perhaps better to use ginput for this
%     setString(hPoints(i),num2str(i))
    manualPoints(i,:) = getPosition(hPoints(i));  %#ok<AGROW> % to plot later use:  plot(manualPoints(:,1), manualPoints(:,2), 'y', 'LineStyle', 'none', 'Marker', '+')
    
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

%% Callbacks %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
                           
function doneCallback(hObject, ~)
    delete(hObject);
    exit = true;
end

end %end of function

