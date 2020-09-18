% Pete Sheurpukdi - 7/31/2020, helper to grab boundary coordinates from
% image for brain slices.

function [boundary] = detectBoundary(image, num_boundaries_to_return)
    
%     [counts,x] = imhist(image,16);
%     stem(x,counts)
%     T = otsuthresh(counts);
%     BW = imbinarize(image,T);

    image = imgaussfilt(image,2);

%     [counts,x] = imhist(image,256);
    [counts,x] = imhist(image, 1000);

    %     stem(x,counts)
    T = otsuthresh(counts);

    BW = imbinarize(image, T);

    % BW = imcomplement(BW);
    % BW = bwareaopen(BW,20);

    figure;
    imshow(BW);
    % % 
    % 
    [B,L] = bwboundaries(BW,'noholes');
    imshow(label2rgb(L, @jet, [.5 .5 .5]))
    hold on

    % sort index ascending length
    [~,index] = sort(cellfun(@length,B));

    index = index(end:-1:1); % reverses index

    boundary = B(index);

    if length(boundary) < num_boundaries_to_return
        boundary = boundary(1:length(boundary)); % save only the largest 5 boundaries
    else
        boundary = boundary(1:num_boundaries_to_return);
    end


%     for k = 1:1
%        btest = boundary{1};
%        plot(btest(:,2), btest(:,1), 'r', 'LineWidth', 2)
%     end