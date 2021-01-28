function [lPGCScore, allPairIntersect, allPairUnion] = lPGCCoor(expressMat, coordinate, threshold)
% + lPGCCoor compare all lPGC (locally-scaled polar gini curve) correlations
% between two genes in a single-cell expression matrix. It is highly
% recommended that the matrix is preprocessed such that the cell are within
% a cluster and the genes pass some certain expressing threshold. Only
% choosing genes that express in at least 40% of the cells is recommended.
% + overview workflow: (expressMat, coordinate) --> lPGC -->
% (allPairInterest, allPairUnion) --> lPGCScore
% + usage: [lPGCScore, allPairIntersect, allPairUnion] = lPGCCoor(expressMat)
% - expressMat: single-cell expression matrix. The matrix could be in both
% sparse and full format. Each row represents a gene, each column
% represents a cell.
% - coordinate: the cell 2D embedding (tSNE, umap, and others) of expressMat.
% For each gene, the expression and its expressing cell embedding are used
% to create its l-PGC
% - threshold (optional): the list of thresholds to choode highly
% expressing (foreground) 
% - lPGCScore: the l-PGC correlation, between 0 to 1, for all gene pairs.
% This is a square matrix; the dimension is the number of genes in the
% input matrix. Higher lPGCScore implies more positive correlation.  The matrix is upper-triangle
% formula: lPGCScore = allPairIntersect ./ allPairUnion (pairwise division)
% - allPairIntersect: the intersecting area between two l-PGCs
% corresponding to a pair of genes.
% - allPairUnion: the union area between two l-PGCs
% corresponding to a pair of genes.

assert(nargin >= 2, 'We need at least the expression matrix (expressMat) and the embedding (coordinate)');
if nargin < 3 %if the threshol is not supplied, set all gene threshold to 0
    threshold = zeros( size(expressMat, 1), 1 );
end

numGene = size(expressMat, 1); % number of genes
numCell = size(expressMat, 2); % number of cells

allPairIntersect = zeros(numGene);
allPairUnion = zeros(numGene);

for i = 1 : numGene
    
    markerExpression = full( expressMat(i, :) )';
    index1 = find(markerExpression > threshold(i));
    [angleList, allGini1] = make2DGiniMinMax([coordinate; coordinate(index1, :)],...
        [ones(numCell, 1); 2*ones(length(index1), 1)], {'background', 'foreground'}, [1, 1]);
    allGini1 = allGini1{2}; % l-PGC for gene i
    
    for j = (i+1) : numGene % may change to parfor j = (i+1) : numGene if having enough memory for parallel processing
        markerExpression = full( expressMat(j, :) )';
        index1 = find(markerExpression > threshold(j));
        [angleList, allGini2] = make2DGiniMinMax([coordinate; coordinate(index1, :)],...
            [ones(numCell, 1); 2*ones(length(index1), 1)], {'background', 'foreground'}, [1, 1]);
        allGini2 = allGini1{2}; % l-PGC for gene j
        
        % compute the union area between l-PGCs of genes i and j in polar coordinate
        allJoinPt = max(allGini1, allGini2);
        angleFactor = 2*pi/1000;
        allJoinPtArea = 0;
        for k = 1 : length(allJoinPt)-1
            allJoinPtArea = allJoinPtArea + allJoinPt(k)^2 * angleFactor / 2;
        end
        allPairUnion(i, j) = allJoinPtArea;
        
        % compute the intersect area between l-PGCs of genes i and j in polar coordinate
        allIntersectPt = min(allGini1, allGini2);
        allIntersectPtArea = 0;
        for k = 1 : length(allIntersectPt)-1
            allIntersectPtArea = allIntersectPtArea + allIntersectPt(k)^2 * angleFactor / 2;
        end
        allPairIntersect(i, j) = allIntersectPtArea;
    end
end

lPGCScore = allPairIntersect ./ allPairUnion;