% script to analyze & plot total results from all brains/sections

gfpThresh = 1;
tdtThresh = 1;

% for each section, read in summary data
[percentGfpCellsOfTotal, percentTdtOfGfpTotal, percentTdtOfGfpDm, percentTdtOfGfpVl] = vmhCountCells(Gfp, Tdt, gfpThresh, tdtThresh);
