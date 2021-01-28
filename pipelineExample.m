expressMat = importdata('AdjustCount.mat'); % lod expression matrix
geneTable = readtable('filterGeneList.csv'); % load the gene list
gene = table2cell(geneTable(:, 2));

[coordinate, idx, all_lPGC, all_lPGC_pair] = lPGCCoor_simplePipeline(expressMat, gene); % call the simplest pipeline from the example