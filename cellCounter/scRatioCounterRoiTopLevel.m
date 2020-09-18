clear all; close all;

filepathNd2 = 'C:\data\Jason\microscope\2017_02_EsrSpinalRatioCounts\';
mouseNums = {'a3', 'ae6', 'ae8', 'ec4', 'ec5', 'ec7'}; %nd2 paths
numSections = {[1:1:12], [1:1:31], [1:1:13], [1:1:25], [1:1:27], [1:1:27]}; %#ok<*NBRAK> %how many sections in each nd2 file
% filepathNd2 = 'C:\data\Jason\microscope\2017_02_CrhSpinalRatioCounts\';
% mouseNums = {'a2', 'a4', 'ac9', 'ac10', 'cc5', 'cc6', 'cc7', 'cc9'}; %nd2 paths
% numSections = {[1:1:20], [1:1:18], [1:1:33], [1:1:23], [1:1:28], [1:1:20], [1:1:28], [1:1:30]}; %how many sections in each nd2 file

totalMice = size(mouseNums, 2);
zplanes = 1; %only deal with maxIP images
    
for k = 1:totalMice
    mouseNum = mouseNums{1,k} %#ok<NOPTS>
    filepathTotal = [filepathNd2, mouseNum, '.nd2'];
    tframes = numSections{1,k};
    [Iaxon, Inissl] = readNd2(filepathTotal, tframes);
    
    for j = min(tframes):max(tframes) %for all sections
        j %#ok<NOPTS>
        currentNissl = Inissl(:,:,j);
        currentAxon = Iaxon(:,:,j);
        [currentNisslRot, angle] = rotateImage(currentNissl);
        currentAxonRot = imrotate(currentAxon,angle,'bilinear','loose');
        roiPatch = sectionSetRoiRectangle(currentNisslRot);
        
        matName = [filepathNd2, mouseNum, '_Section', num2str(j), '.mat'];
        save(matName, 'currentNisslRot', 'currentAxonRot', 'roiPatch');% save rotated image, ROI rectangle
    end

end

