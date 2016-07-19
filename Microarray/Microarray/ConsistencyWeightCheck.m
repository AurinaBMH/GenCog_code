% this script compares coexpression differences varying consistency when
% forming a group matrix and weights. Purpose - to set up a threshold for
% binarising the matrix. Coexpression between distance Uncorrected ROIs
% should be higher for connected versus Connected nodes.
%clear all; close all;

cd('/Users/Aurina/GoogleDrive/Genetics_connectome/Gen_Cog/Data/Microarray/');
load ('FACT_cust100_HCP_good.mat');


DSs = [2 5 10 50 100];
Thresholds = [0 0.1 0.2 0.3 0.4 0.5 0.6 0.7 0.8 0.9 1];
NumWeights = 1000;
weightMeasure = {'density'};
WEIGHT = density; 


for l=1:length(DSs)
    
    DS = DSs(l);
    cd('/Users/Aurina/GoogleDrive/Genetics_connectome/Gen_Cog/Data/Microarray/')
    load(sprintf('%dDSSpearmanExpressionCust100CorrectedSigmoid.mat', DS));

    Tvals = zeros(length(Thresholds), NumWeights);
    Pvals = zeros(length(Thresholds), NumWeights);
    
    for i=1:length(Thresholds)
        
        grpThr = Thresholds(i);
        adjGrp = connectomeGroupThreshold(WEIGHT, grpThr, 2);
        Adj = adjGrp(1:100, 1:100);
        
        if strcmp(weightMeasure{1}, 'density') || strcmp(weightMeasure{1}, 'count')
            
            wmin = min(log(nonzeros(Adj(:)))); wmax = max(log(Adj(:)));
            weights = linspace(wmin,wmax,NumWeights);
            
            for j=1:length(weights)
                
                w = weights(j);
                mask = log(Adj)>w;
                [y,e, Con, Uncon, t, p] = CoexpConnectVSUnconnect(Exp,mask);
                Tvals(i,j) = t;
                Pvals(i,j) = p;
                
                
            end
            
        elseif strcmp(weightMeasure{1}, 'FA')
            
            wmin = min(nonzeros(Adj(:))); wmax = max(Adj(:));
            weights = linspace(wmin,wmax,NumWeights);
            
            for j=1:length(weights)
                
                w = weights(j);
                mask = Adj>w;
                [y,e, Con, Uncon, t, p] = CoexpConnectVSUnconnect(Exp,mask);
                Tvals(i,j) = t;
                Pvals(i,j) = p;
                
                
            end
            
        end
        
        
        
    end
    
    figure; imagesc(Tvals); colormap([flipud(BF_getcmap('blues',20));1,1,1;BF_getcmap('reds',20)]);
    caxis([-7 7]);
    title(sprintf('%dDS Corrected T values', DS)); xlabel('Weight threshold'); ylabel('Consistency between subjects, %'); colorbar;
    
    
end

% [maxValue, linearIndexesOfMaxes] = max(Tvals(:));
% [rowsOfMaxes colsOfMaxes] = find(Tvals == maxValue);
%
% consist = 0.3;
% weightID = 63;
%
% % check distribution for max value
% adjGrp = connectomeGroupThreshold(FA, consist, 2);
% Adj = adjGrp(1:100, 1:100);
% % for each group threshold scan through weights
% wmin = min(log(nonzeros(Adj(:)))); wmax = max(log(Adj(:)));
% weights = linspace(wmin,wmax,NumWeights);
% w = weights(weightID);
% mask = log(Adj)>w;
% [y,e, Con, Uncon, t, p] = CoexpConnectVSUnconnect(Exp,mask);
%
%
% [row col] = max(Tvals);
