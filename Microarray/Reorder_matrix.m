%% Reorder coexpression matrix according to distances between assigned samples
% 
% cd ('/gpfs/M2Scratch/Monash076/aurina/Gen_Cog/code/Microarray/S01-S06_combined');
% load('82_DistThresh5_aparcaseg_combined_ExpressionProbePCA.mat');
% 
% Coord = DataCoordinatesCombinedSorted(:,1:3);
% Dist = pdist2(Coord, Coord);
% % reorder distances according to similarity.
% d_row = pdist(Dist,'corr');
% links_row = linkage(d_row);
% [~,~,ord_row] = dendrogram(links_row,0);
% 
% %calculate coexpression
% Exp = DataExpressionCombinedSorted(:,2:19344);
% ExpN = BF_NormalizeMatrix(Exp); % normalise expression
% R = corr(ExpN');
% 
% % plot reordered distances and coexpression (reordered accodring to distances)
% figure; imagesc(R(ord_row, ord_row)); title('Coexpression reordered according to distances between samples'); 
% figure; imagesc(Dist(ord_row,ord_row)); colormap(flipud(colormap)); title('Distances between samples reordered'); 
% 
% %

%dataMatrixNorm = BF_NormalizeMatrix(DataExpressionCombinedSorted(:,3:19345));

% for i=1:6
%     
% A = DataExpression{i,1};
% data = A(:,3:19345);
% DataNORM = BF_NormalizeMatrix(data); 
% 
% Inf = A(:,1:2); 
% DataExpNORM{i} = cat(2, Inf, DataNORM);
% end
% 
% Combined = cat(1,DataExpNORM{1}, DataExpNORM{2}, DataExpNORM{3}, DataExpNORM{4}, DataExpNORM{5}, DataExpNORM{6});


% dataMatrixNorm = Combined(:,3:19345); 

Data = DataExpressionCombinedSorted(:,3:size(ExpSubj,2));
dataMatrixNorm = BF_NormalizeMatrix(Data); 

Dist = pdist2(dataMatrixNorm, dataMatrixNorm);
figure; imagesc(Dist);
d_row = pdist(Dist,'corr');
links_row = linkage(d_row);
[~,~,ord_row] = dendrogram(links_row,0); 

% Subj = Combined(:,2); 
% ROI = Combined(:,1); 
Subj = DataExpressionCombinedSorted(:,1);
ROI = DataExpressionCombinedSorted(:,2);
figure; subplot(1,25,1); imagesc(Subj(ord_row)); 
        subplot(1,25,2); imagesc(ROI(ord_row)); 
        subplot(1,25, [3 25]);  imagesc(dataMatrixNorm(ord_row,:)); caxis([-2 2]);
% 
%         
% 

%dataMatrixNorm = BF_NormalizeMatrix(DataExpressionCombinedSorted(:,3:19345));

% DataExpNORM = cell(6,1);
% for i=1:6
%     
% A = DataTable.Expression{i,1}; 
% A =A';
% S = zeros(size(A,1),1);
% S(:,1)=i;
% 
% DataExpNORM{i} = cat(2, S, A);
% end
% 
% Combined = cat(1,DataExpNORM{1}, DataExpNORM{2}, DataExpNORM{3}, DataExpNORM{4}, DataExpNORM{5}, DataExpNORM{6});
% 
% 
% dataMatrixNorm = Combined(:,3:19345); 
% 
% Data = DataExpressionCombinedSorted(:,3:19345);
% dataMatrixNorm = BF_NormalizeMatrix(Data); 
% 
% Dist = pdist2(dataMatrixNorm, dataMatrixNorm);
% figure; imagesc(Dist);
% d_row = pdist(Dist,'corr');
% links_row = linkage(d_row);
% [~,~,ord_row] = dendrogram(links_row,0); 
% 
% Subj = Combined(:,2); 
% ROI = Combined(:,1); 
% %Subj = DataExpressionCombined(:,2);
% %ROI = DataExpressionCombined(:,1);
% figure; subplot(1,25,1); imagesc(Subj(ord_row)); 
%         subplot(1,25,2); imagesc(ROI(ord_row)); 
%         subplot(1,25, [3 25]);  imagesc(dataMatrixNorm(ord_row,:));