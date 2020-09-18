% top level script to count all sections and plot results
clear all; close all;
filepathNd2 = 'C:\data\Jason\cellCounting\esrChR2\';

viewOverlay = false;
loadDataInstead = false;
nd2format = false;  %else assume jpeg; NOTE: changing this may require changing "uint8" to "uint16" in some places
computeStats = false;
somaColor = 'all';
countGreen = false;
countRed = true;
countBlue = false;
countRgOverlap = false;
useRoiMask = false;

%%% INPUT:
mouseNums = {'ec9'};
% sections =  {'1' '2' '3'};
% sections =  {'1' '2' '3' '4' '5' '6'}; 
% sections =  {'1R' '2R' '3R' '1L' '2L' '3L'}; 
sections =  {'1R' '2R' '3R' '4R' '5R' '6R' '1L' '2L' '3L' '4L' '5L' '6L'}; 
numSections = size(sections, 2);
totalMice = size(mouseNums, 2);

numRedTotal = 0;
numGreenTotal = 0;
numBlueTotal = 0;
numOverlapTotal = 0;

for k = 1:totalMice
    mouseNum = mouseNums{1,k};
    for j = 1:numSections
        thisSection = sections{1,j};
        if loadDataInstead
            matName = [filepathNd2, mouseNum, '\', thisSection, '.mat']; %use saved data already counted
            load(matName); %careful not to save/load loop variables!!!
            loadDataInstead = true;
        else %count cells manually:
            display(['counting ', mouseNum, ', section '  thisSection]);
            if nd2format
                maxImageBitValue = 2^16-1;
                maxPixelValue = 2^12-1;
                [Idapi, Ig, Ir, Ib] = readNd2_4ch(filepathNd2, mouseNum, thisSection);
            else %assume JPG
                maxImageBitValue = 2^8-1;
                maxPixelValue = 2^8-1;
                [Ir, Ig, Ib, Irgb] = readJpg(filepathNd2, mouseNum, thisSection);
            end
            
            % first set color for cell body/soma and select ROI:
            if strcmp(somaColor,'blue')
                Iroi = Ib;
            elseif strcmp(somaColor,'red') %for magenta or using Crh-tdT-defined area
                Iroi = Ir;
            elseif strcmp(somaColor,'all') %for magenta or using Crh-tdT-defined area
                Iroi = Irgb;
            end
            [patchMask, maskCoord, ellip, xMin, xMax, yMin, yMax] = sectionSetRoiEllipseJpeg(Iroi); % select ROI
            
            
            IroiMask = patchMask(yMin:yMax,xMin:xMax);
            
            % BLUE (first):
            if countBlue
                IsubB = Ib(yMin:yMax,xMin:xMax);
                if useRoiMask
                    IsubB(~IroiMask) = 0;
                end
                [numBlue, manualPointsIb] = manualCount(IsubB);
            else
                numBlue = 0;
                manualPointsIb = [];
            end
            
            % RED:
            if countRed
                IsubR = Ir(yMin:yMax,xMin:xMax);
                if useRoiMask
                    IsubR(~IroiMask) = 0;
                end
                [numRed, manualPointsIr] = manualCount(IsubR);
            else
                numRed = 0;
                manualPointsIr = [];
            end
            
            % GREEN:
            if countGreen
                IsubG = Ig(yMin:yMax,xMin:xMax);
                if useRoiMask
                    IsubG(~IroiMask) = 0;
                end
                [numGreen, manualPointsIg] = manualCount(IsubG);
            else
                numGreen = 0;
                manualPointsIg = [];
            end
            
            % RG OVERLAP:
            if countRgOverlap
                pixelShift = 1;
                IsubGshift = circshift(IsubG,pixelShift,1); %shift G channel pixelShift pixels in first dimension to aid overlap count
                IsubRG = zeros(size(IsubR,1),size(IsubR,2),3); %maximize R & G channels & set B=0
                IsubRG(:,:,1) = double(IsubR)./ (double(max(max(Ir)))*2); % divide red by extra factor of 2 since it is much brighter
                IsubRG(:,:,2) = double(IsubGshift)./ double(max(max(Ig)));
                [numOverlap, colocalPoints] = manualCount2Ch(IsubRG)
                percentRColocalized = numOverlap/numRed
                percentGColocalized = numOverlap/numGreen
            else
                numOverlap = 0;
                colocalPoints = [];
                percentRColocalized = 0;
                percentGColocalized = 0;
            end
            
            matName = [filepathNd2, mouseNum, '\', thisSection, '.mat']; %save data for each individual section:
            save(matName, 'patchMask', 'maskCoord', 'ellip', 'xMin', 'xMax', 'yMin', 'yMax', 'Irgb',...
                          'numBlue', 'manualPointsIb',...
                          'numGreen', 'manualPointsIg',...
                          'numRed', 'manualPointsIr',...
                          'numOverlap', 'colocalPoints', 'percentRColocalized', 'percentGColocalized');  % save relevant variables
            close all; %prevent large number figures
        
            if viewOverlay % to create and view RGB overlay of ROIs, to test:
                Ioverlay = rgb2gray(Irgb(yMin:yMax,xMin:xMax,:));
                figure; hOverlay = imshow(Ioverlay); 
                hold on;
                if countRed
                    plot(manualPointsIr(:,1), manualPointsIr(:,2), 'r', 'LineStyle', 'none', 'Marker', '+', 'MarkerSize', 4)
                end
                if countGreen
                    plot(manualPointsIg(:,1), manualPointsIg(:,2), 'g', 'LineStyle', 'none', 'Marker', '+', 'MarkerSize', 4)
                end
                if countBlue
                    plot(manualPointsIb(:,1), manualPointsIb(:,2), 'b', 'LineStyle', 'none', 'Marker', '+', 'MarkerSize', 2)
                end
                if countRgOverlap
                    plot(colocalPoints(:,1), colocalPoints(:,2), 'y', 'LineStyle', 'none', 'Marker', 'x', 'MarkerSize', 6)
                end
                hold off;
                saveas(hOverlay, [filepathNd2, mouseNum, '\', mouseNum, '_', num2str(j), '_manualOverlay.jpg'], 'jpg');
            end
        
        end
        
        if computeStats
            numRedTotal = numRedTotal + numRed; %#ok<*UNRCH>
            numBlueTotal = numBlueTotal + numBlue;
            numGreenTotal = numGreenTotal + numGreen;
            numOverlapTotal = numOverlapTotal + numOverlap;
        end
        
    end    
end

