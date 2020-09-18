% filepath = 'C:\data\jason\peterVirusData\';
% % filepath = 'C:\Users\keller\School\UCSD\stowersLab\code\cellSegmentation\';
% % mouseNum = '19249';
% % sectionNum = '8';
% % mouseNum = '70324';
% % sectionNum = '2';  % this one has distortion and is very dim
% mouseNum = '19485';
% sectionNum = '4';  % this one is very bright

% Note that these are actually 12-bit images, so MATLAB will represent with
% uint16; however, the highest value (saturation) will be 2^12 - 1 = 4095,
% so we can just normalize by this to make a standard colormap

% possibly use imadjust/histeq/adapthisteq to fill 16 bits to optimize contrast at the 16-bit
% level, but this can also warp image intesities significantly

max16BitValue = 2^16-1;
maxPixelValue = 2^12-1;

% filepathTotal = [filepath, mouseNum, '\', sectionNum, '.nd2'];
filepathTotal = 'C:\data\peter\fromJason\2014_Ai9FearAggr\aggrHomecage\aggr1C_20xVmh_section6_corectShading592crop20overlap.nd2';

zplanes = 1; tframes = 1; channel = 1;  %only one channel can be imported at once
[vol] = imreadBF(filepathTotal,zplanes,tframes,channel); %read 1st channel in ND2 file - Hoechst
Ib = uint16(vol);
% Ib = imadjust(uint16(vol));
channel = 2;
[vol] = imreadBF(filepathTotal,zplanes,tframes,channel); %read 2nd channel in ND2 file - GFP
Ir = uint16(vol);
% meanTdt = mean(mean(Ir)); % perhaps use total image mean as threshold later
clear vol; %free up memory

% uncomment below to create and view RGB image to test:
Irgb = zeros(size(Ib,1),size(Ib,2),3);
Irgb(:,:,1) = double(Ir)./ double(max(max(Ir))); 
Irgb(:,:,3) = double(Ib)./ double(max(max(Ib)));
figure; imshow(Irgb);