function channelArray = readNd2(filepathIn, tframes, channelsToRead)
% pete NOTE: sometimes have to change from uint16() functions to uint8()
% depending on the image type, 8 or 16bit

% readNd2

% read nd2 directly to cut out the ElementsViewer step (can also access image metadata)
% using imreadBF function from MATLAB Central: http://www.mathworks.com/matlabcentral/fileexchange/32920-imread-for-multiple-life-science-image-file-formats 
% NOTE that Java heap size probably needs to be increased in MATLAB Preferences

% Note that these are actually 12-bit images, so MATLAB will represent with
% uint16; however, the highest value (saturation) will be 2^12 - 1 = 4095,
% so we can just normalize by this to make a standard colormap

% possibly use imadjust/histeq/adapthisteq to fill 16 bits to optimize contrast at the 16-bit
% level, but this can also warp image intesities significantly

% nd2 files have color channels, as well as frames, each with their own
% color channel.
% channel order - green, red, blue

% numZplanes = ;
% zplanes = 2; %only deal with maxIP images % PETE THIS CONTROLS WHICH FRAME IN THE SEQUENCE
% zplanes [1 2 3 4 5 6 7];
% channel = 1;  %only one channel can be imported at once

% [vol] = imreadBF(filepathIn,zplanes,tframes,channel); %read 1st channel in ND2 file -
% Iaxon = uint16(vol);
% 
% Ich1 = imadjust(uint16(vol));

channelArray = {};

nd2FileType = '.nd2';
jp2FileType = '.jp2';

disp('displaying channels');
channelsToRead
for i = 1:length(channelsToRead)
%     disp('fusing z-stack if exists');
    channel = channelsToRead(i);
    meta = imreadBFmeta(filepathIn);
%     filepathIn
%     disp(meta.channels)
%     
%     figure;
%     imagesc(uint8(vol));
    
    if ~contains(filepathIn, nd2FileType)
        disp(['loading jp2 ', filepathIn , ' slice ', num2str(tframes), ' ch  ',num2str(channel)]);
        [vol] = imreadBF(filepathIn, 1 ,tframes, channel);
%         vol2 = uint8(vol); % lower half of intensity
%         vol; % will contain upper half of intensity vals
        vol= uint16(vol);
    else
        disp(['loading nd2 ', filepathIn , ' slice ', num2str(tframes), ' ch  ',num2str(channel)]);
        [vol] = imreadBF(filepathIn, 1, tframes, channel); % sometimes
%         invalid T index error... 
%         [vol] = imreadBF(filepathIn, tframes,1, channel);
        vol = uint16(vol);
    end
    
    fusedImCh = vol;
    % sometimes if the image looks weird, might have to check the image
    % matrix, and if there are many negatives, try to fix it as below?

    % COMMENT THESE 2 BLOCKS OUT IF 16-bit, and KEEP IF 8-bit
%     if ~contains(filepathIn, nd2FileType)
%         vol = vol*-1;
%         vol(vol > 0) = 256 - vol(vol > 0);
%     end
% 
%     if ~contains(filepathIn, nd2FileType)
%         fusedImCh= uint8(vol); 
%     else
%         fusedImCh= uint16(vol); 
%     end
    
    figure; 
    imagesc(fusedImCh);
    title(['CH: ' num2str(channel)]);

    %% FILTERING BLOCK
%     fusedImCh = medfilt2(fusedImCh, [3 3]); % original: [4 4], use [m n] of ~cell size, can change values. 
%     vol2 = medfilt2(vol2, [3 3]);
    %     gaussian filter
%     H1 = fspecial('gaussian', [4 4], 1);
%     fusedImCh = imfilter(fusedImCh, H1, 'replicate');
%     vol2 = imfilter(vol2, H1, 'replicate');
    
%% BLENDING Z-STACK BLOCK
%     for zNum = 2:meta.zsize % blending together sections
%         disp('blending z');
%         meta.zsize
%         [vol] = imreadBF(filepathIn, zNum ,tframes,channel); % read next z-slice
%         if ~contains(filepathIn, nd2FileType)
%             nextImCh = uint8(vol);  
%         else
%             nextImCh = uint16(vol);
%         end
% 
%         nextImCh = medfilt2(nextImCh, [4 4]); % use [m n] of ~cell size [4 4], PETE: sometimes use [5 5] if lots of veins to filter out
%         H2 = fspecial('gaussian', [4 4], 1);
%         nextImCh = imfilter(nextImCh, H2, 'replicate');
%         %% image intensity scaling done here
%         if mean2(nextImCh) < 20.0 % pete: value found by testing. red channel usually very dark. intensity scaling is really bad if the image is dark, or not a lot of contrast. it will make contrast out of nothing
%             fusedImCh = imfuse(fusedImCh, nextImCh, 'blend', 'Scaling', 'none');  % SOMETIMES HAVE TO CHANGE THIS MANUALLY from 'none' to 'joint' depending on contrast
%         else
%             fusedImCh = imfuse(fusedImCh, nextImCh, 'blend', 'Scaling', 'joint');
%         end
%     end
%     if ~contains(filepathIn, nd2FileType)
%         Ich = imadjust(uint8(fusedImCh)); %(uint16)
%         vol2 = imadjust(uint8(vol2));
%     else
        Ich = (uint16(fusedImCh)); %(uint16)
%     end
        
    
    % show image as guide to select area - PETE
%     figure;
%     Rfilt = medfilt2(Ich, [4 4]); % use [m n] of ~cell size
%     H = fspecial('gaussian', [4 4], 1);
%     Rfilt = imfilter(Rfilt, H, 'replicate');  %also smooth out a bit to help with oversegmentation w/ regionalmax
%     imagesc(Ich);
%     title(['Filtered Channel ' num2str(channel)]);
    if ~contains(filepathIn, nd2FileType)
        channelArray{end+1} = Ich;
%         channelArray{end+1} = vol2;
    else
        channelArray{end+1} = Ich;
    end
end


%%
% show filtered image as guide to select area - PETE
% figure;
% Rfilt = medfilt2(Ich1, [4 4]); % use [m n] of ~cell size
% H = fspecial('gaussian', [4 4], 1);
% Rfilt = imfilter(Rfilt, H, 'replicate');  %also smooth out a bit to help with oversegmentation w/ regionalmax
% RIm = imagesc(Rfilt);
% title("1st Channel Filtered")
% 
% %green
% figure;
% Gfilt = medfilt2(Ich2, [4 4]); % use [m n] of ~cell size
% H = fspecial('gaussian', [4 4], 1);
% Gfilt = imfilter(Gfilt, H, 'replicate');  %also smooth out a bit to help with oversegmentation w/ regionalmax
% GIm = imagesc(Gfilt);
% title("2nd Channel Filtered")
% 
% figure;
% Bfilt = medfilt2(Ich3, [4 4]); % use [m n] of ~cell size
% H = fspecial('gaussian', [4 4], 1);
% Bfilt = imfilter(Bfilt, H, 'replicate');  %also smooth out a bit to help with oversegmentation w/ regionalmax
% BIm = imagesc(Bfilt);
% title("3nd Channel Filtered")
% 

end %end function
