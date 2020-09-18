function [numColocalized, percent1Colocalized, percent2Colocalized, colocalPoints] = countColocalize( manualPoints1, manualPoints2)
%countColocalize: function to count colocalized cells based on manual X&Y coordinates' euclidean distances
pixelDistanceThresh = 3; %for now use 3 for 512x512 10x image

size1 = size(manualPoints1, 1);
size2 = size(manualPoints2, 1);
% concat = cat(1,manualPoints1,manualPoints2); %concatenate all points
% 
% distances = pdist2(concat(:,1),concat(:,2),'euclidean'); %computes distances between ALL concatenated points
% keyboard

if size1 > size2
%     distances = distances(1:size1,size1+1:end);
    d = zeros(size1,size2);
    for i = 1:size1
        for j = 1:size2
            xDist = (manualPoints1(i,1) - manualPoints2(j,1))^2;
            yDist = (manualPoints1(i,2) - manualPoints2(j,2))^2;
            d(i,j) = sqrt(xDist + yDist);
        end
    end
else %size2 >= size1
%     distances = distances(size1+1:end,1:size1); 
    d = zeros(size2,size1);
    for i = 1:size2
        for j = 1:size1
            xDist = (manualPoints1(i,1) - manualPoints2(j,1))^2;
            yDist = (manualPoints1(i,2) - manualPoints2(j,2))^2;
            d(i,j) = sqrt(xDist + yDist);
        end
    end
end

dHalf = tril(d); % take lower triangular part to avoid duplicates
[X,Y] = find(dHalf < pixelDistanceThresh & dHalf ~= 0); 
numColocalized = length(X);
if size1 > size2
    colocalPoints = manualPoints1(X,:);
else
    colocalPoints = manualPoints2(X,:);
end
percent1Colocalized = numColocalized/size1;
percent2Colocalized = numColocalized/size2;

end  % end function