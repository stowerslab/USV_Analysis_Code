function [I1, I2, I3, I4] = readNd2_4ch(filepathIn, mouseNum, sectionNum)
% readNd2

% read nd2 directly to cut out the ElementsViewer step (can also access image metadata)
% using imreadBF function from MATLAB Central: http://www.mathworks.com/matlabcentral/fileexchange/32920-imread-for-multiple-life-science-image-file-formats 
% NOTE that Java heap size probably needs to be increased in MATLAB Preferences

% Note that these are actually 12-bit images, so MATLAB will represent with
% uint16; however, the highest value (saturation) will be 2^12 - 1 = 4095,
% so we can just normalize by this to make a standard colormap

% possibly use imadjust/histeq/adapthisteq to fill 16 bits to optimize contrast at the 16-bit
% level, but this can also warp image intesities significantly

filepathTotal = [filepathIn, mouseNum, '\', sectionNum, '.nd2'];
zplanes = 1; tframes = 1; channel = 1;  %only one channel can be imported at once
[vol] = imreadBF(filepathTotal,zplanes,tframes,channel); %read 1st channel in ND2 file
I1 = uint16(vol);
% Ib = imadjust(uint16(vol));
channel = 2;
[vol] = imreadBF(filepathTotal,zplanes,tframes,channel); %read 2nd channel in ND2 file
I2 = uint16(vol);
channel = 3;
[vol] = imreadBF(filepathTotal,zplanes,tframes,channel); %read 3rd channel in ND2 file
I3 = uint16(vol);
channel = 4;
[vol] = imreadBF(filepathTotal,zplanes,tframes,channel); %read 4th channel in ND2 file
I4 = uint16(vol);

% uncomment below to create and view RGB image to test:
% Irgb = zeros(size(Ib,1),size(Ib,2),3);
% maxPixelValue = 2^12-1;
% Irgb(:,:,1) = double(Ir)./ maxPixelValue; 
% Irgb(:,:,2) = double(Ig)./ maxPixelValue;
% Irgb(:,:,3) = double(Ib)./ maxPixelValue;
% imshow(Irgb);
% keyboard
end %end function
