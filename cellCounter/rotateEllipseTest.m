clear all; close all;

% I = imread('test.jpg');
I = imread('C:\data\Jason\microscope\2016_12_EsrChR2Pmc\EC4\slide_pmc_xy04.tif');
I = I(:,:,3);  %Nissl channel only
% + OPTION TO USE GREEN CHANNEL
hFig = figure;
hIm = imagesc(I);
hAxes = gca;
axis off;  %turn off axis labels, etc.

hCenter = impoint(hAxes);  %first draw center point manually
set(hCenter, 'Visible', 'off');
pos = getPosition(hCenter);
xCenter = pos(1);
yCenter = pos(2);

[hVmhVlPatchMask, maskCoord, ellip] = rotateEllipse(xCenter, yCenter, hAxes, hFig, hIm);