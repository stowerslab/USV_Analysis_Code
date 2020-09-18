%plot PMC Esr vs. Crh pie chart
% totalNissl = 1000;
% esr = 250;
% crh = 500;
% overlap = 50;
% esrOnlyPercent = (esr - overlap)/totalNissl;
% crhOnlyPrecent = (crh - overlap)/totalNissl;
% overlapPercent = overlap/totalNissl;
% otherPercent = 1 - esrOnlyPercent - crhOnlyPrecent - overlapPercent;
esrOnlyPercent = .192-.052;
crhOnlyPrecent = .599-.052;
overlapPercent = .052;
otherPercent = 1 - .192 - .599 + .052;

labels = {'Esr1','overlap','Crh', 'other'};
slices = [esrOnlyPercent overlapPercent crhOnlyPrecent otherPercent];

fontSz = 20;

hPmcPieChart = figure;
pie(slices,labels)
colormap([0 1 0; .5 .5 .5; 1 0 1; 0 0 1])      %green, grey, magenta, blue



% axis equal, axis off
% A = [totalNissl esr crh]; I = [esr crh overlap 0];
% %For three-circle venn diagrams, A is a three element vector [c1 c2 c3], 
% %and I is a four element vector [i12 i13 i23 i123], specifiying the 
% %two-circle intersection areas i12, i13, i23, and the three-circle
% %intersection i123.
% venn(A,I,'FaceColor',{'b','g','m'},'FaceAlpha',{0.5,0.5,0.5},'EdgeColor','black')
% 
% 


