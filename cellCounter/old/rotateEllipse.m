function [hVmhVlPatchMask, maskCoord, vlVertPoints, ellip] = rotateEllipse(xCenter, yCenter, hAxes, hFig, hIm) %mainfunction
%rotateEllipse - Interactively build and manipulate ellipse ROI, return a
%structure of all relevant ellipse parameters

totalPatchColor = [1 0 0]; % use red for patch outline
subPatchColor = [0 1 0]; %use green for VL suboutline
hold on;

% first build the ellipse (this should probably be converted to OOP):
numVertPoints = 100; % # of points as vertices
ellip.vertPoints = 0:360/numVertPoints:360;  % do everything in degrees, not radians
ellip.theta = 0;  % initial angle from x-axis in degrees
ellip.a = 600;  % initial major semiaxis length 
ellip.b = 300; % initial minor semiaxis length 
ellip.center = [xCenter yCenter];
[x, y] = makeEllipse(ellip.a, ellip.b, ellip.center, ellip.theta, ellip.vertPoints);
verts = [x' y'];  % putting into column format for rotation matrix; these are original to be read by slider callback for rotation
ellip.x = verts(:,1);  %initialize
ellip.y = verts(:,2);
ellip.left = true;

% Decide if we are using the right or left side of the section:
choice = questdlg('Using right or left side of image?', ...
	'R or L?', ...
	'Left','Right','Right');
% Handle response
switch choice
    case 'Left'
        ellip.left = true;
    case 'Right'
        ellip.left = false;
end

% make a patch to visualize ellipse now:
hPatch = patch('Parent',hAxes,'Visible','on',...  %draw colored patch to visualize ROI
           'XData',verts(:,1),'YData',verts(:,2),'EdgeColor',totalPatchColor,'LineWidth',2,'FaceColor','none');%
                           
hDoneControl = uicontrol(hFig, 'Style', 'pushbutton', 'String', 'Done Rotating',...  %make button to press when done manipulating
                          'Callback', @doneCallback, 'Position', [50 50 150 20]);      

set(hFig,'KeyPressFcn',@keyCallback);  %set a callback function for key pressed while the figure is in focus
                                            
waitfor(hDoneControl);  %wait until hDoneControl is deleted, which will happen upon button press

% Once manipulation is finished, draw another patch for VMHvl boundary (+/- 60 degrees in parametric equation):
if ellip.left % if on left of image
    vlVertPoints = 105:180/(numVertPoints/3):255; %+/- 75 degrees, between 1/2 and 1/3
else % if on right of image
    vlVertPoints = -75:180/(numVertPoints/3):75;
end
[x, y] = makeEllipse(ellip.a, ellip.b, ellip.center, ellip.theta, vlVertPoints);
hVlPatch = patch('Parent',hAxes,'Visible','on',...  %draw colored patch to visualize ROI
           'XData',x,'YData',y,'EdgeColor',subPatchColor,'LineWidth',2,'FaceColor','none');%

xVert = get(hVlPatch, 'XData');
yVert = get(hVlPatch, 'YData');
hVmhVlPatchMask = poly2mask(xVert, yVert, size(get(hIm,'CData'),1), size(get(hIm,'CData'),2));  %create binary mask from area enclosed by X & Y data

%also return min-max of mask for later plotting:
maskCoord.xMin = uint16(min(xVert)-5);
maskCoord.xMax = uint16(max(xVert)+5);
maskCoord.yMin = uint16(min(yVert)-5);
maskCoord.yMax = uint16(max(yVert)+5);

hold off;

%%%%%%%%%%%%%%%%%%%% Callbacks %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
                           
function doneCallback(hObject, ~)
    delete(hObject);
end

function keyCallback(~, ~)
    val = double(get(hFig,'CurrentCharacter')); %get value of last get pressed
    switch val
        case 28 %leftArrow - move center left
            moveCenter(-2,0);
        case 29 %rightArrow - move center
            moveCenter(2,0);      
        case 30 %downArrow - move center
            moveCenter(0,-2);
        case 31 %upArrow - move center
            moveCenter(0,2);
        case 97 %'a' - increase major axis
            chgAxes(5,0);
        case 115 %'s' - decrease major axis
            chgAxes(-5,0);
        case 98 %'b' - increase minor axis
            chgAxes(0,5);
        case 110 %'n' - decrease minor axis
            chgAxes(0,-5);
        case 114 %'r' - rotate clockwise
            chgRotation(1);
        case 116 %'t' - rotate counterclockwise
            chgRotation(-1);
        otherwise
            %do nothing
    end   
end

%%%%%%%%%%%%%%%%%%%% Other Functions %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

function moveCenter(xChg, yChg)
    % subroutine to move the center of the ellipse x/y pixels and update all
    % associated values (x & y s/b positive or negative integers)
    ellip.center = [ellip.center(1) + xChg, ellip.center(2) + yChg];
    [x, y] = makeEllipse(ellip.a, ellip.b, ellip.center, ellip.theta, ellip.vertPoints);
    verts = [x' y'];  % putting into column format for rotation matrix;
    ellip.x = verts(:,1); 
    ellip.y = verts(:,2);
    set(hPatch, 'XData',ellip.x,'YData',ellip.y);
end

function chgAxes(aChg, bChg)
    % subroutine to update major (by aChg) and minor (by bChg) axes of the ellipse and update all
    % associated values (a & b s/b positive or negative integers)
    ellip.a = ellip.a + aChg;
    ellip.b = ellip.b + bChg;
    [x, y] = makeEllipse(ellip.a, ellip.b, ellip.center, ellip.theta, ellip.vertPoints);
    verts = [x' y'];  % putting into column format for rotation matrix;
    ellip.x = verts(:,1); 
    ellip.y = verts(:,2);
    set(hPatch, 'XData',ellip.x,'YData',ellip.y);
end

function chgRotation(thetaChg)
    % subroutine to update theta (by thetaChg) and update all
    % associated values (thetaChg s/b positive or negative integer)
    ellip.theta = ellip.theta + thetaChg;
    [x, y] = makeEllipse(ellip.a, ellip.b, ellip.center, ellip.theta, ellip.vertPoints);
    verts = [x' y'];  % putting into column format for rotation matrix;
    ellip.x = verts(:,1); 
    ellip.y = verts(:,2);
    set(hPatch, 'XData',ellip.x,'YData',ellip.y);
end

function [xNew, yNew] = makeEllipse(a, b, center, theta, vertPoints)
    % subroutine to update x & y vertices values of the ellipse from
    % trig parametric equation (see Wikipedia for equation details)
    xNew = center(1) + a*cosd(theta)*cosd(vertPoints) - b*sind(theta)*sind(vertPoints); % x-coords of vertices
    yNew = center(2) + a*sind(theta)*cosd(vertPoints) + b*cosd(theta)*sind(vertPoints); % y-coords
end

end
