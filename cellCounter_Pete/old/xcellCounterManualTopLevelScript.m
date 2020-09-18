%need to reset loop
%variables (& do not save them)!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!

% top level script to count all sections and plot results
clear all; close all;

filepathNd2 = 'C:\data\Jason\microscope\2016_01_CrhTdtwithEsrImmuno\female\';
viewOverlay = false;
loadDataInstead = true;
nd2format = false;  % NOTE: changing this may require changing "uint8" to "uint16" in some places
computeStats = true;

%%% INPUT:
mouseNums = {'esrF1' 'esrF2'};
lastSections =  {8, 5}; % to keep track of how many sections were imaged for each mouse
firstSections = {1, 1};

totalMice = size(mouseNums, 2);

for k = 1:totalMice
    mouseNum = mouseNums{1,k};
    lastSection = lastSections{1,k};
    firstSection = firstSections{1,k};
    for j = firstSection:lastSection
        if loadDataInstead
            matName = [filepathNd2, mouseNum, '\', mouseNum, '_', num2str(j), '.mat']; %use saved data already counted
            load(matName); %need to reset loop variables!!!
        else %count cells manually:
            if nd2format
                maxImageBitValue = 2^16-1;
                maxPixelValue = 2^12-1;
                [Ir, Ig, Ib] = readNd2(filepathNd2, mouseNum, num2str(j));
            else %assume JPG
                maxImageBitValue = 2^8-1;
                maxPixelValue = 2^8-1;
                [Ir, Ig, Ib, Irgb] = readJpg(filepathNd2, mouseNum, num2str(j));
            end
            % RED:
            [roiPatchArray, xMin, xMax, yMin, yMax] = sectionSetRoi(Ir); % use Tdt channel in this case to select ROI
            IsubR = Ir(yMin:yMax,xMin:xMax);
            IsubRMask = roiPatchArray.patchMask(yMin:yMax,xMin:xMax);
            IsubR(~IsubRMask) = 0;
            [numRed, manualPointsIr] = manualCount(IsubR);
            % GREEN:
            IsubG = Ig(yMin:yMax,xMin:xMax);
            IsubG(~IsubRMask) = 0;
            [numGreen, manualPointsIg] = manualCount(IsubG);
            % BLUE:
            IsubB = Ib(yMin:yMax,xMin:xMax);
            IsubB(~IsubRMask) = 0;
            [numBlue, manualPointsIb] = manualCount(IsubB);
            % OVERLAP:
            IsubRG = zeros(size(IsubR,1),size(IsubR,2),3); %maximize R & G channels & set B=0
            IsubRG(:,:,1) = double(IsubR)./ (double(max(max(Ir)))*2); % divide red by extra factor of 2 since it is much brighter
            IsubRG(:,:,2) = double(IsubG)./ double(max(max(Ig)));
            [numOverlap, colocalPoints] = manualCount2Ch(IsubRG);
            
            percentRColocalized = numOverlap/numRed;
            percentGColocalized = numOverlap/numGreen;
            
            matName = [filepathNd2, mouseNum, '\', mouseNum, '_', num2str(j), '.mat']; %save data for each individual section:
            save(matName);  % save all variables
            close all; %prevent large number figures
        
            if viewOverlay % to create and view RGB overlay of ROIs, to test:
                Ioverlay = rgb2gray(Irgb(yMin:yMax,xMin:xMax,:));
                figure; hOverlay = imshow(Ioverlay); 
                hold on;
                plot(manualPointsIr(:,1), manualPointsIr(:,2), 'r', 'LineStyle', 'none', 'Marker', '+', 'MarkerSize', 4)
                plot(manualPointsIg(:,1), manualPointsIg(:,2), 'g', 'LineStyle', 'none', 'Marker', '+', 'MarkerSize', 4)
                plot(manualPointsIb(:,1), manualPointsIb(:,2), 'b', 'LineStyle', 'none', 'Marker', '+', 'MarkerSize', 2)
                plot(colocalPoints(:,1), colocalPoints(:,2), 'y', 'LineStyle', 'none', 'Marker', 'x', 'MarkerSize', 6)
                hold off;
                saveas(hOverlay, [filepathNd2, mouseNum, '\', mouseNum, '_', num2str(j), '_manualOverlay.jpg'], 'jpg');
            end
        
        end
    end    
end

