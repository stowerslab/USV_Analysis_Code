% top level script to count all VMH sections and plot results
clear all; close all;

filepathRoot = 'C:\data\jason\peterVirusData\vmhIso';  %root for all mice

%%% INPUT folder names for individual mice:
% mouseNums = {'v12','v15','v16','v17','v20',...
%              'v21','v23','v24','v25','v26','v28','v29','v30'};  %ex. 2nd entry is mouseNums{1,2}
 mouseNums = {'v30'}; 

totalMice = size(mouseNums, 2);

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% do all manual steps for each section first (i.e. ROIs and thresholds):
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% for k = 1:totalMice
%     mouseNum = mouseNums{1,k};
%     for j = 1:4  %loop over sections for each mouse
%         % do image loading and ROIs for each section first:
%         useBlue = true;  %true if using blue
%         [Ir, Ig, Ib, hVmhVlPatchMask, ellip] = vmhSetRoiEllipse(filepathRoot, mouseNum, num2str(j), useBlue); 
%         %save data for each individual section:
%         matName = [filepathRoot, '\', mouseNum, '\', mouseNum, '_', num2str(j), '_roiImages.mat'];
%         save(matName, 'Ir', 'Ig', 'Ib', 'hVmhVlPatchMask', 'ellip');
%         close all; %prevent large number figures
%         
%         % Use the green channel to draw ROI alternatively:
%         choice = questdlg('Using green to draw another ROI?', ...
%             'a query?', ...
%             'Yes','No','No');
%         % Handle response
%         switch choice
%             case 'Yes'
%                 useBlue = false;
%             case 'No'
%                 useBlue = true;
%         end
%         
%         if ~useBlue %make an ROI using GFP channel and overwrite that using blue
%             [Ir, Ig, Ib, hVmhVlPatchMask, ellip] = vmhSetRoiEllipse(filepathRoot, mouseNum, num2str(j), useBlue); 
%             %save data for each individual section:
%             matName = [filepathRoot, '\', mouseNum, '\', mouseNum, '_', num2str(j), '_roiImages.mat'];
%             save(matName, 'Ir', 'Ig', 'Ib', 'hVmhVlPatchMask', 'ellip');
%             close all; %prevent large number figures
%         end
%         
%     end
% end

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% Now uncomment below to do all of the counting:
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
for k = 1:totalMice
    mouseNum = mouseNums{1,k};
    for j = 1:4
        matName = [filepathRoot, '\', mouseNum, '\', mouseNum, '_', num2str(j), '_roiImages.mat'] %#ok<NOPTS>
        load(matName);
        [numBlue, avgBlue, manualPointsIb] = vmhManualCount(Ib,hVmhVlPatchMask); %blue only
        %save data for each section:
        matName = [filepathRoot, '\', mouseNum, '\', mouseNum, '_', num2str(j), '_blueCount2.mat'];
        save(matName, 'numBlue', 'avgBlue', 'manualPointsIb');
        close all; %prevent large number figures
    end
 
end




