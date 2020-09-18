function [IrSub, IgSub, IbSub, hVmhVlPatchMaskCrop, ellip] = vmhSetRoiEllipse(filepathIn, mouseNum, sectionNum, useBlue)
%% vmhSetRoiEllipse: function to read image and pick elliptical ROI in order to quantify viral expression in the VMH
% Jason Keller, 12th Mar 2013
%
% This script does some manual steps first (picking ROIs and
% loading images), so that counting set up.
%
% Input is 3-channel TIF, output is separated channels for separate
% counting, each as an image matrix saved to MAT file, already masked & cropped for the VMHvl. The
% VMHvl mask will also be saved on top of GFP channel as JPG so that viral coverage can
% be determined easily.
%
% The counting scripts will then load the separate channel files and calculate average intensities, and 
% save them along with handles for all of the points chosen as cells.
%
%
%% image read & format
% read TIF - Note that these are actually 12-bit images, so MATLAB will represent with
% uint16; however, the highest value (saturation) will be 2^12 - 1 = 4095,
% so we can just normalize by this to make a standard colormap

% possibly use imadjust/histeq/adapthisteq to fill 16 bits to optimize contrast at the 16-bit
% level, but this can also warp image intesities significantly

max16BitValue = 2^16-1;
maxPixelValue = 2^12-1;

filepathTotal = [filepathIn, '\' mouseNum, '\', sectionNum, '.tif'] %#ok<NOPRT>
Ib = imread(filepathTotal, 1); %read 1st image in TIF file - Nissl
% Ig = imread(filepathTotal, 2); %read 2nd image in TIF file - GFP
% Ir = imread(filepathTotal, 3); %read 3rd image in TIF file - TDT
Ir = imread(filepathTotal, 2); %read 3rd image in TIF file - TDT
Ig = Ib;  %temp

% uncomment below to create and view RGB image to test:
% Irgb = zeros(size(Ib,1),size(Ib,2),3);
% Irgb(:,:,1) = double(Ir)./ double(max(max(Ir))); 
% Irgb(:,:,2) = double(Ig)./ double(max(max(Ig)));
% Irgb(:,:,3) = double(Ib)./ double(max(max(Ib)));
% imshow(Irgb);

%% ROI selection
%use Nissl channel only to make ROI (can substitute GFP on case-by-case basis)
hFig = figure;
if(useBlue)
    hIm = imagesc(Ib);  %autoscale to see borders better, and fill 16 bits
else
    hIm = imagesc(Ig); %alternatively, use green channel to draw ROI
end
set(hFig, 'Name', ['Section#', mouseNum, '_', sectionNum]);
hAxes = gca;
axis off;  %turn off axis labels, etc.
colormap(gray);
% set(gcf, 'Position', get(0,'Screensize')); % maximize figure - CAUTION- may distort ellipse

hCenter = impoint(hAxes);  %first draw ellipse center point manually
set(hCenter, 'Visible', 'off');
pos = getPosition(hCenter);
xCenter = pos(1);
yCenter = pos(2);

% now call function to manipulate (using R/T, A/S, B/N, and arrow keys) ellipse and save:
[hVmhVlPatchMask, maskCoord, vlVertPoints, ellip] = rotateEllipse(xCenter, yCenter, hAxes, hFig, hIm);
% figure; imshow(hVmhVlPatchMask)

%read subimage X & Y limits to plot only ROI:
xMinAlt = maskCoord.xMin;
xMaxAlt = maskCoord.xMax;
yMinAlt = maskCoord.yMin;
yMaxAlt = maskCoord.yMax;

% create cropped images, & set to 0 outside ROI
hVmhVlPatchMaskCrop = hVmhVlPatchMask(yMinAlt:yMaxAlt, xMinAlt:xMaxAlt);
IrSub = Ir(yMinAlt:yMaxAlt, xMinAlt:xMaxAlt);  IrSub(~hVmhVlPatchMaskCrop) = 0;
IgSub = Ig(yMinAlt:yMaxAlt, xMinAlt:xMaxAlt);  IgSub(~hVmhVlPatchMaskCrop) = 0;
IbSub = Ib(yMinAlt:yMaxAlt, xMinAlt:xMaxAlt);  IbSub(~hVmhVlPatchMaskCrop) = 0;

% save a copy of GFP channel with white patch outline, to make sure that the VL has been fully filled:
figure; imagesc(Ig); colormap(gray);
[x, y] = makeEllipse(ellip.a, ellip.b, ellip.center, ellip.theta, vlVertPoints);
hVlPatch = patch('Parent',gca,'Visible','on',...  %draw colored patch to visualize ROI
           'XData',x,'YData',y,'EdgeColor',[1 0 0],'LineWidth',2,'FaceColor','none');%
axis off;
saveas(gcf,[filepathIn, '\', mouseNum, '\', sectionNum, '_vlGfpFill.jpg'], 'jpg');

function [xNew, yNew] = makeEllipse(a, b, center, theta, vertPoints)
    % subroutine to update x & y vertices values of the ellipse from
    % trig parametric equation (see Wikipedia for equation details)
    xNew = center(1) + a*cosd(theta)*cosd(vertPoints) - b*sind(theta)*sind(vertPoints); % x-coords of vertices
    yNew = center(2) + a*sind(theta)*cosd(vertPoints) + b*cosd(theta)*sind(vertPoints); % y-coords
end

end  %end function vmhSetRoiEllipse