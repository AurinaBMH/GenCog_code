% connection length and weight relationship in a group connectome
clear all; close all;

DS = 5;

cd('/Users/Aurina/GoogleDrive/Genetics_connectome/Gen_Cog/Data/Microarray/')
load(sprintf('%dDSSpearmanExpressionCust100CorrectedSigmoid.mat', DS));
load('FACT_custom200ANDaseg20.mat');

Thresholds = [0.3 0.6 0.8];
weightMeasure = {'count'};
WEIGHT = count;
NumWeights = 1000;


ConUncon = zeros(5, 2);
ConUnconSD = zeros(5, 2, length(DS));
% get connection distances
for i=1:88
    L(i,:,:) = coords{1,i};
end

Len = mean(L,1);
Len = squeeze(Len);
l = Len(1:100,:)';
l = dist(l);

for i=1:length(Thresholds)
    
    grpThr = Thresholds(i);
    adjGrp = connectomeGroupThreshold(WEIGHT, grpThr, 2);
    Adj = adjGrp(1:100, 1:100);
    %lmask = l.*logical(Adj);
    
    
    if strcmp(weightMeasure{1}, 'density') || strcmp(weightMeasure{1}, 'count')
        
        wmin = min(log(nonzeros(Adj(:)))); wmax = max(log(Adj(:)));
        weights = linspace(wmin,wmax,NumWeights);
        k=1;
        FigHandle = figure;
        set(FigHandle, 'Position', [100, 100, 700, 1400]);
        for j=[100 200 400 650 800]
            
            w = weights(j);
            mask = log(Adj)>w;
            lmaskCon = l.*mask;
            maskUncon = 1-mask;
            lmaskUncon = l.*maskUncon;
            lmaskUncon(logical(eye(size( lmaskUncon)))) = 0.001;
            [y,e, dataCell,t, p] = CoexpConnectVSUnconnect(Exp,mask);
          
            
            
            B = [nonzeros(lmaskCon(:)) dataCell{1,1}(:)];
            B = sortrows(B);
            C = [nonzeros(lmaskUncon(:)) dataCell{2,1}(:)];
            C = sortrows(C);
            
            subplot(5,1,k); scatter(C(:,1), C(:,2)); hold on; scatter(B(:,1), B(:,2), '.r');
            title(sprintf('%d group threshold\nWeight threshold %d\nCoexpression for connected versus unconnected links as a function of distance', grpThr*100, j));
            xlabel('Euclidean distance between ROIs'); ylabel('Coexpression'); 
            k=k+1;
        end
        
        
        
    end
end
    
    
    
    