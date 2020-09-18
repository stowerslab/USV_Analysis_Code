%% image read & format
% read TIF - Note that these are actually 12-bit images, so MATLAB will represent with
% uint16; however, the highest value (saturation) will be 2^12 - 1 = 4095,
% so we can just normalize by this to make a standard colormap

% possibly use imadjust/histeq/adapthisteq to fill 16 bits to optimize contrast at the 16-bit
% level, but this can also warp image intesities significantly

max16BitValue = 2^16-1;
maxPixelValue = 2^12-1;
filepathIn = 'C:\data\jason\peterAi9data\2014_Ai9FearAggr\aggrCastrateWithUrine\';
filepathTotal = [filepathIn, 'aggr3f_20xvmh_section04_512crop.tif'];
Ib = imread(filepathTotal, 1); %read 1st image in TIF file - Nissl
Ir = imread(filepathTotal, 2); %read 2nd image in TIF file - TDT

% uncomment below to create and view RGB image to test:
% Irgb = zeros(size(Ib,1),size(Ib,2),3);
% Irgb(:,:,1) = double(Ir)./ double(max(max(Ir))); 
% Irgb(:,:,2) = double(Ig)./ double(max(max(Ig)));
% Irgb(:,:,3) = double(Ib)./ double(max(max(Ib)));
% imshow(Irgb);

