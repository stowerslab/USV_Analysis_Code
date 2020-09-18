Ibfilt = medfilt2(Ig, [3 3]); %GREEN!
H = fspecial('gaussian', [5 5], 1);
Ibfilt2 = imfilter(Ibfilt, H, 'replicate');

Ibthresh = bradley(Ibfilt2, [100 100], 5);  %I like this one for its simplicity
% Ibthresh = niblack(Ibfilt2, [300 300], 0.1); %OK as well, but bradley is more simple
% Ibthresh = wolf(Ibfilt2, [300 300], 0.5);  %too complicated, similar, result to niblack
% Ibthresh = nick(Ibfilt2, [300 300], 0.1); %too complicated, similar, result to niblack
% Ibthresh = sauvola(Ibfilt2, [300 300], 0.34); %too complicated, similar, result to niblack

Ibthresh = bwareaopen(Ibthresh,50,8); % removes from a binary image all connected components (objects) that have fewer than X pixels; 50 pixels seems to be nice for Nikon 20x stitch
Ibdist = bwdist(~Ibthresh); %compute distance transform (effectively makes small grayscale 'basins' around each cell)
Ibdist = -Ibdist;  %preprocessing from MATLAB 'watershed' example code
Ibdist2 = Ibdist; %create a copy to manipulate basins and combine small ones (presumably where nuclear staining is punctate) to prevent oversegmentation
Ibdist2(~Ibthresh) = -Inf; %from MATLAB 'watershed' example code - set basins & background to trenches

%clean up bwdist output before watershed to prevent oversegmentation:
bw = imregionalmin(Ibdist);  %regional min predicts watershed - but each nucleus may have several local minima, so we will try to combine them
se2 = strel('disk', 1);  %we will dilate the regional min with a disk to connect entire nuclei
bw2 = imdilate(bw,se2);  %dilate to connect together very-regional minima within nuclei; note that there is a tradeoff here in that some cells that overlap a lot will be combined
bw3 = bwmorph(bw2, 'bridge'); %connect caddy-corner pixels; this should have little impact
bw4 = bwmorph(bw3,'thin', 2); %now thin back down to skeleton before setting the basins;  n=Inf takes a long time
Ibdist2(bw4) = -Inf;

Ibwater = watershed(Ibdist2, 8);
IbwaterRegions = label2rgb(Ibwater,'jet', 'w', 'shuffle'); 

figure; 
subplot(1,3,1); 
imagesc(Ibdist2(middle(1)-200:middle(1)+200,middle(2)-200:middle(2)+200));
subplot(1,3,2); 
imagesc(IbwaterRegions(middle(1)-200:middle(1)+200,middle(2)-200:middle(2)+200));
subplot(1,3,3); 
imagesc(Ibfilt(middle(1)-200:middle(1)+200,middle(2)-200:middle(2)+200));
colormap(jet)

% NOTES for local threshold:
% need window greater than cell radius, so that average will be pulled down by background
% 