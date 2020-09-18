subimage = Ibfilt(middle(1)-100:middle(1)+100,middle(2)-100:middle(2)+100);
suborig = Ib(middle(1)-100:middle(1)+100,middle(2)-100:middle(2)+100);
figure; 
subplot(2,3,1);
imagesc(suborig); title('Original')
% hy = fspecial('sobel');
% hx = hy';
% Iy = imfilter(double(subimage), hy, 'replicate');
% Ix = imfilter(double(subimage), hx, 'replicate');
% gradmag = sqrt(Ix.^2 + Iy.^2);
% figure, imshow(gradmag,[]), title('Gradient magnitude (gradmag)')

se = strel('disk', 3);  %disk3 seems to smooth out single nuclei

% Io = imopen(subimage, se);
% Ioc = imclose(Io, se);
% subplot(1,3,2);
% % imagesc(Io); title('Opening (Io)');
% imagesc(Ioc); title('Opening/closing (Ioc)');

Ie = imerode(subimage, se);
Iobr = imreconstruct(Ie, subimage);
Iobrd = imdilate(Iobr, se);
Iobrcbr = imreconstruct(imcomplement(Iobrd), imcomplement(Iobr));
Iobrcbr = imcomplement(Iobrcbr);
subplot(2,3,2);
% imagesc(Iobr), title('Opening-by-reconstruction (Iobr)')
imagesc(Iobrcbr), title('Opening/closing-by-reconstruction (Iobrcbr)')

subimageThresh = uint8(Iobrcbr .* (255/maxPixelValue));
thresh = threshToolMod(subimageThresh) / 255 * (maxPixelValue/max16BitValue); % thresh must be normalized between 0 & 1 for 16-bit image;
% thresh =  0.0132; %this was manually chose using tool
Ibthresh = im2bw(Iobrcbr,thresh);


% fgm = imregionalmax(Iobrcbr,8);
subplot(2,3,3);
% imagesc(fgm), title('Regional maxima of opening-closing by reconstruction (fgm)')
imagesc(Ibthresh); title('threshd')

Ibthresh = bwareaopen(Ibthresh,50,8); % removes from a binary image all connected components (objects) that have fewer than X pixels
Ibdist = bwdist(~Ibthresh); %compute distance transform (effectively makes small grayscale 'basins' around each cell)
Ibdist = -Ibdist;  %preprocessing from MATLAB 'watershed' example code
Ibdist2 = Ibdist;
Ibdist2(~Ibthresh) = -Inf; %from MATLAB 'watershed' example code

%clean up bwdist output before watershed to prevent oversegmentation:
bw = imregionalmin(Ibdist);  %regional min predicts watershed - but each nucleus may have several, so we will try to combine them
se2 = strel('disk', 3);  %we will dilate the regional min with a disk to connect entire nuclei
bw2 = imdilate(bw,se2);  %dilate to connect together very-regional minima within nuclei; note that there is a tradeoff here in that some cells that overlap a lot will be combined
bw3 = bwmorph(bw2, 'bridge'); %connect caddy-corner pixels; this should have little impact
bw4 = bwmorph(bw3,'thin', 2); %now thin back down to skeleton before setting the basins; n=Inf takes a long time
Ibdist2(bw4) = -Inf;
subplot(2,3,4);
imagesc(bw4); title('regional min, cleaned up with morph')
% imagesc(Ibdist2); title('before watershed')

% TO DO: the watershed algorithm can oversegment and overestimate the number
% of nuclei, so perhaps try without (i.e. use bwconncomp to label cells) as well and average results
Ibwater = watershed(Ibdist2,8);
IbwaterRegions = label2rgb(Ibwater,'jet', 'w', 'shuffle'); 
subplot(2,3,5);
imagesc(suborig); title('Original')
subplot(2,3,6);
imshow(IbwaterRegions)

% bgm = Ibwater == 0;
% figure, imshow(bgm), title('Watershed ridge lines (bgm)')

















