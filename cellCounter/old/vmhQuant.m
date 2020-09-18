function [Gfp, Tdt, numCells] = vmhQuant(Ig, Ir, Ibwater, roiPatchArray)
%% vmhQuant: function to quantify viral expression in the VMH
% Jason Keller, 5th Dec 2013
% 
% This function mainly just does the region property calculations now...
%
% TO DO / IDEAS: 
%       prob. rename function...
%
% inputs:
%       roiPatchArrayInput, input ROIs from previous segmentation (in format set below)
%       Ig, Ir, and Ibwater (from vmhSetRoiThresh.m) - note that Ibwater only contains watershed regions within the VMH mask, others have already been set to zero
%
% outputs:
% (A)   struct of GFP mean intensities in areas overlapping with segmented Hoechst cells
% (B)   struct of cFosCre TDT mean intensities in areas overlapping with segmented Hoechst cells

%% Calculate cell properties
% for every label/region/nucleus found by the segmentation above,
% store in an array the mean GFP intensity in that region, and the mean TDT

% Note that the Ibwater region corresponding to '1' is the dead space between nuclei, so subtract 1
% idxLargeEnough = find([stats.Area] > 50);
% numCells.total = length(idxLargeEnough) - 1; %
numCells.total = max(max(Ibwater)) - 1; 

IbwaterTemp = Ibwater; %make temp version of Ibwater to set to zero outside ROI
IbwaterTemp(roiPatchArray(1,2).patchMask == 0) = 1;  %set to background
IbwaterTemp(IbwaterTemp ~= 1) = 0;    % set others to zero (same as region outlines)
% figure; imagesc(IbwaterTemp)
[L, NUM] =  bwlabeln(~IbwaterTemp); % now relabel connected components & recount
numCells.dm = NUM; 

IbwaterTemp = Ibwater; %make temp version of Ibwater to set to zero outside ROI
IbwaterTemp(((roiPatchArray(1,1).patchMask) & (~roiPatchArray(1,2).patchMask)) == 0) = 1;  %set to background
IbwaterTemp(IbwaterTemp ~= 1) = 0;    % set others to zero (same as region outlines)
% figure; imagesc(IbwaterTemp)
[L, NUM] =  bwlabeln(~IbwaterTemp); % now relabel connected components & recount
numCells.vl = NUM; 

Gfp = struct(   ...    %structure array containing mean GFP intensity data for each cellular ROI found in the image        
                                'total',uint16(zeros(numCells.total,1)),...          %total VMH
                                'dm',uint16(zeros(numCells.total,1)),...   %VMHdm
                                'vl',uint16(zeros(numCells.total,1)));      %VMHvl                                                                                
                            
Tdt = struct(   ...    %structure array containing mean TDT intensity data for each cellular ROI found in the image        
                                'total',uint16(zeros(numCells.total,1)),...          %total VMH
                                'dm',uint16(zeros(numCells.total,1)),...   %VMHdm
                                'vl',uint16(zeros(numCells.total,1)));      %VMHvl     

IgDm = Ig; IgDm(~roiPatchArray(1,2).patchMask) = 0; %second ROI mask in for DM - set pixels outside to zero
IgVl = Ig; IgVl(~((roiPatchArray(1,1).patchMask) & (~roiPatchArray(1,2).patchMask))) = 0; %ROI mask in for VL is ~DM - set pixels outside to zero
IrDm = Ir; IrDm(~roiPatchArray(1,2).patchMask) = 0; 
IrVl = Ir; IrVl(~((roiPatchArray(1,1).patchMask) & (~roiPatchArray(1,2).patchMask))) = 0;

% Note that the Ibwater region corresponding to '1' is the dead space between nuclei, so only take (2:end)                          
GfpStatsTotal = regionprops(Ibwater, Ig, 'MeanIntensity');
GfpStatsDm = regionprops(Ibwater, IgDm, 'MeanIntensity');
GfpStatsVl = regionprops(Ibwater, IgVl, 'MeanIntensity');
Gfp.total = [GfpStatsTotal(2:end).MeanIntensity];
Gfp.dm = [GfpStatsDm(2:end).MeanIntensity];
Gfp.vl = [GfpStatsVl(2:end).MeanIntensity];

TdtStatsTotal = regionprops(Ibwater, Ir, 'MeanIntensity');
TdtStatsDm = regionprops(Ibwater, IrDm, 'MeanIntensity');
TdtStatsVl = regionprops(Ibwater, IrVl, 'MeanIntensity');
Tdt.total = [TdtStatsTotal(2:end).MeanIntensity];
Tdt.dm = [TdtStatsDm(2:end).MeanIntensity];
Tdt.vl = [TdtStatsVl(2:end).MeanIntensity];

% ex. 'Gfp-positive' cells are those that are one standard deviation above mean intensity:
% numGfpDm = length(find(Gfp.dm > (mean(Gfp.dm)+std(Gfp.dm))));
% numTdtDm = length(find(Tdt.dm > (mean(Tdt.dm)+std(Tdt.dm))));
% numGfpTotal = length(find(Gfp.total > (mean(Gfp.total))));
% efficiencyGfp = numGfpTotal/double(numCells.total);

% % uncomment below to create and view RGB overlay of ROIs, to test (not working yet):
% xMin = round(uint16(vmhXLim(1)));
% xMax = round(uint16(vmhXLim(2)));
% yMin = round(uint16(vmhYLim(1)));
% yMax = round(uint16(vmhYLim(2)));
% temp = Ibwater(yMin:yMax, xMin:xMax);
% Irgb = zeros(size(temp,1),size(temp,2),3);
% Irgb(:,:,1) = roiPatchArray(1,1).patchMask(yMin:yMax, xMin:xMax); %red is full VMH ROI
% Irgb(:,:,2) = (roiPatchArray(1,1).patchMask(yMin:yMax, xMin:xMax)) & (roiPatchArray(1,2).patchMask(yMin:yMax, xMin:xMax)); %green is DM
% Irgb(:,:,3) = ~((roiPatchArray(1,1).patchMask(yMin:yMax, xMin:xMax)) & (~roiPatchArray(1,2).patchMask(yMin:yMax, xMin:xMax))); %blue is VL
% imshow(Irgb);
end %end of function

