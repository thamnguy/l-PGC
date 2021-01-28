function [coordinate, idx, all_lPGC, all_lPGC_pair] = lPGCCoor_simplePipeline(expressMat, gene)
% this is the 'baseline' pipeline of locally-scaled Polar Gini Curve
% (l-PGC) for a single-cell dataset.
% usage : [coordinate, idx, all_lPGC, all_lPGC_pair] = lPGCCoor_simplePipeline(expressMat, gene)
% Input: 
% - expressMat: single-cell expression matrix. The matrix could be in both
% sparse and full format.
% - gene: list of genes in the expression. The order of gene is according
% to the order of rows in expressMat.
% Output:
% - coordinate: embedding cordinate. Here, I choose Matlab tSNE
% - idx: clusterID for each cell. Here I use the dbscan algorithm to find
% cluster ID.
% - all_lPGC: l-PGC for each gene in each cluster. To plot the l-PGC, run:
% polarplot(pi*(0:360/resolution:360) / 180, all_lPGC{i, j}). This plots
% l-PGC of gene i in cluster j.
% - all_lPGC_pair: all matrix of pairwise l-PGC correlations in all
% clusters. For cluster i, all_lPGC_pair{i, 1} is the matrix,
% all_lPGC_pair{i, 2} is the list of genes associating with the matrix.

coordinate = tsne(expressMat'); % compute tSNE embedding
idx = dbscan(coordinate, 0.6, 50); % carefully chose the right parameters to run dbscan!

% for each gene in each cluster, compare the % of cell expressing the gene
% in the cluster
percenExp = zeros(length(gene), max(idx));
for i = 1 : length(gene)
    for j = 1 : max(idx)
        percenExp(i, j) = length( find( expressMat(i, clusterIndex) > 0 ) ) / length (clusterIndex);      
    end
end

%% all l-PGC for all genes in the cluster
angleFactor = 2*pi/1000;
threshold = 0;
all_lPGC = cell(length(gene), max(idx));

for clusterID = 1 : max(idx)
    countCell = length(find(idx==clusterID));
    if countCell > 500 % only plot l-PGC if the cluster has at least 500 cells
        
        for j = 1 : length(gene) % may use  parfor j = 1 : length(gene)
            index = j;
            markerExpression = full( expressMat(index, :) )';
            
            if percenExp(j, clusterID) > 0.05 && length( find(markerExpression>0) ) > 3
                index1 = find(markerExpression > threshold & idx==clusterID);
                [angleList, allGini2] = make2DGiniMinMax([coordinate(clusterIdx, :); coordinate(index1, :)],...
                    [ones(length(clusterIdx), 1); 2*ones(length(index1), 1)], {'background', 'foreground'}, [1, 1]);
                allGini2 = allGini2{2};
                
                for k = 1 : length(allGini1)-1
                    all_lPGC(j, clusterID) = all_lPGC(j, clusterID) + allGini2(k)^2 * angleFactor / 2;
                end
            end
        end
    end
end

%% all l-PGC correlation matrix for all clusters
all_lPGC_pair = cell(max(idx), 2);
for clusterID = 1 : max(idx)
    countCell = length(find(idx==clusterID));
    if countCell > 1000
        
        expressingList = gene (find(percenExp(:, clusterID) > 0.4) );
        % only compute l-PGC correlations with gene expressing in at least 40% of the cluster cell
        all_lPGC_pair{clusterID, 2} = expressingList;
        
        allPairIntersect = zeros(length(expressingList));
        allPairUnion = zeros(length(expressingList));
        
        threshold = 0;
        
        for i = 1 : length(expressingList)
            marker1 = expressingList{i} ;
            [~, index] = ismember(marker1, gene);
            markerExpression = full( expressMat(index, :) )';
            index1 = find(markerExpression > threshold & idx==clusterID);
            clusterIdx = find(idx==clusterID);
            [angleList, allGini1] = make2DGiniMinMax([coordinate(clusterIdx, :); coordinate(index1, :)],...
                [ones(length(clusterIdx), 1); 2*ones(length(index1), 1)], {'background', 'foreground'}, [1, 1]);
            allGini1 = allGini1{2};
            
            for j = (i+1) : length(expressingList) % may use par for j = (i+1) : length(expressingList)
                marker2 = expressingList{j} ;
                [~, index] = ismember(marker2, gene);
                markerExpression = full( expressMat(index, :) )';
                index1 = find(markerExpression > threshold & idx==clusterID);
                [angleList, allGini2] = make2DGiniMinMax([coordinate(clusterIdx, :); coordinate(index1, :)],...
                    [ones(length(clusterIdx), 1); 2*ones(length(index1), 1)], {'background', 'foreground'}, [1, 1]);
                allGini2 = allGini2{2};
                
                allJoinPt = max(allGini1, allGini2);
                angleFactor = 2*pi/1000;
                allJoinPtArea = 0;
                for k = 1 : length(allJoinPt)-1
                    allJoinPtArea = allJoinPtArea + allJoinPt(k)^2 * angleFactor / 2;
                end
                allPairUnion(i, j) = allJoinPtArea;
                
                allIntersectPt = min(allGini1, allGini2);
                allIntersectPtArea = 0;
                for k = 1 : length(allIntersectPt)-1
                    allIntersectPtArea = allIntersectPtArea + allIntersectPt(k)^2 * angleFactor / 2;
                end
                allPairIntersect(i, j) = allIntersectPtArea;
            end
            
        end
        
        all_lPGC_pair{i, 1} = allPairIntersect ./ allPairUnion;
    end
end
