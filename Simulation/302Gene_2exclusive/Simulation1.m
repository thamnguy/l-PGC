%% setup the correlation matrix for the first 300 genes, as showed in supplemental note 4
strongCorr = 0.8 + 0.2*rand(100, 100) + 1; % this indicates positive correlation
antiStrongCorr = -( 0.8 + 0.2*rand(100, 100) ) -1; % this indicate negative correlation
 
corrMat = zeros(300);
 
corrMat(1:100, 101:200) = strongCorr; % block-1 positively correlate to block-2
corrMat(101:200, 1:100) = strongCorr';
 
corrMat(1:100, 201:300) = antiStrongCorr;% block-1 negatively correlates to block-3
corrMat(201:300, 1:100) = antiStrongCorr';
 
corrMat(101:200, 201:300) = 0.1*antiStrongCorr; % block-2 negatively correlates to block-3
corrMat(201:300, 101:200) = 0.1*antiStrongCorr';
 
corrMat(1:100, 1:100) = strongCorr; % genes within the block positively correlate.
corrMat(101:200, 101:200) = strongCorr;
corrMat(201:300, 201:300) = strongCorr;
 
det(corrMat)
 
corrMat = corrMat  * corrMat';
corrMat = corrMat / max(max(corrMat));
 
%% synthesize the expression of the first 300 genes
meanExpress = zeros(300, 1); % the log-scale mean is 0 for these 300 genes
Expression = mvnrnd(meanExpress, corrMat,100000);
Expression = 2.^ Expression;
figure, heatmap(Expression); %heatmap for the expression of these first 300 genes
[coordinate, umap, clusterIdentifiers]=run_umap(Expression); % run umap when the expression matrix only has just these 300 genes. Umap library is from https://www.mathworks.com/matlabcentral/fileexchange/71902-uniform-manifold-approximation-and-projection-umap

figure, gscatter (coordinate(:, 1), coordinate(:, 2)); % plot the umap for these 300 genes
title ('Umap if only using the first 300 genes');
xlabel('umap 1'); ylabel('umap 2');

%% randomly choose 5000 cell to express gene 301; the other 5000 cell expresses gene 302
randVec = randperm(100000);
extraGene = 2.^mvnrnd([0 0], eye(2),100000) - 1; % the log-scale mean is 0 for genes 301 and 302
extraGene(randVec(1:5000), 1) = 0;
extraGene(randVec(5001:100000), 2) = 0;
 
Expression = [Expression, extraGene];
[coordinate, umap, clusterIdentifiers]=run_umap(Expression);

figure, gscatter (coordinate(:, 1), coordinate(:, 2)); % plot the umap for these 302 genes
title ('Umap if only using all 302 genes');
xlabel('umap 1'); ylabel('umap 2');

%% plot gene 301
markerExpression = full( Expression(:, 301) );
threshold = 0;

index1 = find(markerExpression <= threshold );
figure, scatter(coordinate(index1, 1), coordinate(index1, 2), 15, [192, 192, 192]/255, '.');
hold on
index1 = find(markerExpression > threshold );
scatter(coordinate(index1, 1), coordinate(index1, 2), 15, markerExpression(index1), '.');
colorbar;
colormap('jet');
title('Gene 301 expression');
xlabel('umap 1'); ylabel('umap 2');
set(gca,'FontSize',16)
hold off

%% plot gene 302
markerExpression = full( Expression(:, 302) );
threshold = 0;

index1 = find(markerExpression <= threshold );
figure, scatter(coordinate(index1, 1), coordinate(index1, 2), 15, [192, 192, 192]/255, '.');
hold on
index1 = find(markerExpression > threshold );
scatter(coordinate(index1, 1), coordinate(index1, 2), 15, markerExpression(index1), '.');
colorbar;
colormap('jet');
title('Gene 301 expression');
xlabel('umap 1'); ylabel('umap 2');
set(gca,'FontSize',16)
hold off




