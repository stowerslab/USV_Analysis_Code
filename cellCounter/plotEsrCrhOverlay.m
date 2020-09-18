% top level script to count all sections and plot results
clear all; close all;
filepathNd2 = 'C:\Users\keller\School\UCSD\stowersLab\data\microscope\2016_01_CrhTdtwithEsrImmuno\male\';

%%% INPUT:
mouseNums = {'crhTdtMale5wHoechst\count' 'crhTdtMale6wHoechst\count' 'crhTdtMale7wHoechst\count' 'crhTdtMale8wHoechst\count' 'crhTdtMale9wHoechst\count' 'crhTdtMale10wHoechst\count'};
sectionsR = {'2R' '3R' '4R' '5R'};
sectionsL = {'2L' '3L' '4L' '5L'};  
numSectionsR = size(sectionsR, 2);
numSectionsL = size(sectionsL, 2);
totalMice = size(mouseNums, 2);

% numOverlapTotal = 0;

for k = 1:totalMice
    mouseNum = mouseNums{1,k};
    for j = 1:numSectionsR
        thisSectionR = sectionsR{1,j};
        
        matName = [filepathNd2, mouseNum, '\', thisSectionR, '.mat']; %use saved data already counted
        load(matName); %careful not to save/load loop variables!!!
        
        % Calculate centroid and ranges
        crhCentroidX = mean(manualPointsIr(:,1)); %defined as average of Crh points
        crhCentroidY = mean(manualPointsIr(:,2));
        
        % Register oval centroid to middle of much larger area
        crhPointsNormX = manualPointsIr(:,1) - crhCentroidX;
        crhPointsNormY = manualPointsIr(:,2) - crhCentroidY;
        esrPointsNormX = manualPointsIg(:,1) - crhCentroidX;
        esrPointsNormY = manualPointsIg(:,2) - crhCentroidX;
        
        % Plot all cell locations with transparent open cicles
        hold on;
        for cellC = 1:size(manualPointsIr, 1);
            plot_little_circles(crhPointsNormX(cellC), crhPointsNormY(cellC), 5, 'r', 0.5);
        end
        for cellE = 1:size(manualPointsIg, 1);
            plot_little_circles(esrPointsNormX(cellE), esrPointsNormY(cellE), 5, 'g', 0.5);
        end
        
    end 
    
    for j = 1:numSectionsL
        thisSectionL = sectionsL{1,j};
        
        matName = [filepathNd2, mouseNum, '\', thisSectionL, '.mat']; %use saved data already counted
        load(matName); %careful not to save/load loop variables!!!
        
        % Calculate centroid and ranges
        crhCentroidX = mean(manualPointsIr(:,1)); %defined as average of Crh points
        crhCentroidY = mean(manualPointsIr(:,2));
        
        % Register oval centroid to middle of much larger area, and FLIP x-points around (0,0) to match right side
        crhPointsNormX = -(manualPointsIr(:,1) - crhCentroidX);
        crhPointsNormY = manualPointsIr(:,2) - crhCentroidY;
        esrPointsNormX = -(manualPointsIg(:,1) - crhCentroidX);
        esrPointsNormY = manualPointsIg(:,2) - crhCentroidX;
        
        % Plot all cell locations with transparent open cicles
        for cellC = 1:size(manualPointsIr, 1);
            plot_little_circles(crhPointsNormX(cellC), crhPointsNormY(cellC), 5, 'r', 0.5);
        end
        for cellE = 1:size(manualPointsIg, 1);
            plot_little_circles(esrPointsNormX(cellE), esrPointsNormY(cellE), 5, 'g', 0.5);
        end
        
    end 
    
    hold off;
    set(gca, 'XTick', []);
    set(gca, 'YTick', []);
    axis image;
    %plot2svg([filepathNd2, 'CrhEsrOverlay.svg'], gcf)
end

