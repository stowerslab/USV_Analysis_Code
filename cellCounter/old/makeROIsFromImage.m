function [rois, regions]=makeROIsFromImage(ccimage,disp)

winX = 10;% window determined emperically
winY = 20;
temp = ccimage;
temp(isnan(temp)) = []; %discard NaNs to get noise estimate from correlation image
noiseVarEst = evar(temp);
threshOffset = -(noiseVarEst*500); %this is the (negative of the) pixel intensity greater 
                                    %than the local avg that a pixel must be to be included
                                    %Use a value nominally larger than the
                                    %noise - not sure yet if this does any good
% threshOffset = -0.05;
tm = 0; %0 for mean, 1 for median
bw=adaptivethreshold(ccimage,[winX winY],threshOffset,tm);
%minArea=(winX*winY)/4;
minArea=12;

if nargin<2
    disp=0;
end

regions=bwlabeln(bw);
s=regionprops(regions);

roiNum=1;
for r=1:max(regions(:))
   if s(r,1).Area<minArea
       regions(regions==r)=0;
   else
       rois(:,:,roiNum)=double(regions==r);
       roiNum=roiNum+1;
   end
end

regions_th=regions>0;
regions=bwlabeln(regions_th);

if disp==1
    figure
    subplot(1,3,1)
    imagesc(ccimage)
    subplot(1,3,2)
    imagesc(bw)
    subplot(1,3,3)
    imagesc(regions)
end
