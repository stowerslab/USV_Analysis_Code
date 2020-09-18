function [ percentGfpCellsOfTotal, percentTdtOfGfpTotal, percentTdtOfGfpDm, percentTdtOfGfpVl ] = vmhCountCells( Gfp, Tdt, gfpThresh, tdtThresh )
%vmhCountCells: function to count cells and return percentages in different fluo channels & regions

numCellsTotal = length(Gfp.total);
numGfpCellsTotal = length(find(Gfp.total > gfpThresh));
numGfpCellsInDm = length(find(Gfp.dm > gfpThresh));
numGfpCellsInVl = length(find(Gfp.vl > gfpThresh));
numTdtCellsTotal = length(find(Tdt.total > tdtThresh));
numTdtCellsInDm = length(find(Tdt.dm > tdtThresh));
numTdtCellsInVl = length(find(Tdt.vl > tdtThresh));

% note that I'm NOT checking yet if each TDT cell has corresponding GFP expression (this would require a loop over the number of cells) 

percentGfpCellsOfTotal = numGfpCellsTotal / numCellsTotal;  %efficiency of GFP infection in VMH, compared to total nuclei (not total neurons)
percentTdtOfGfpTotal = numTdtCellsTotal / numGfpCellsTotal;
percentTdtOfGfpDm = numTdtCellsInDm / numGfpCellsInDm;
percentTdtOfGfpVl = numTdtCellsInVl / numGfpCellsInVl;

end %end of function

