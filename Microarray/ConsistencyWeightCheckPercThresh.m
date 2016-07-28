% this script compares coexpression differences varying consistency when
% forming a group matrix and weights. Purpose - to set up a threshold for
% binarising the matrix. Coexpression between distance Uncorrected ROIs
% should be higher for connected versus Connected nodes.

%% this script uses [AdjThr] = KeepPercWeights(Adj, perc); function to threshold weights

%clear all; close all;

cd('/Users/Aurina/GoogleDrive/Genetics_connectome/Gen_Cog/Data/Microarray/');
load ('FACT_custom200ANDaseg_withoutdenoise.mat');


DSs = [5];
Thresholds = [0 0.1 0.2 0.3 0.4 0.5 0.6 0.7 0.8 0.9 1];
NumWeights = 100;
weights = linspace(5,100,NumWeights);
weights = flip(weights);
weightMeasure = {'count'};
WEIGHT = count; 


for l=1:length(DSs)
    
    DS = DSs(l);
    cd('/Users/Aurina/GoogleDrive/Genetics_connectome/Gen_Cog/Data/Microarray/')
    %load(sprintf('%dDSSpearmanExpressionCust100CorrectedSigmoid.mat', DS));

    Tvals = zeros(length(Thresholds), NumWeights);
    Pvals = zeros(length(Thresholds), NumWeights);
    ConUncon = zeros(NumWeights, 2);
    ConUnconSD = zeros(NumWeights, 2, length(DSs));
    
    for i=1:length(Thresholds)
        
        grpThr = Thresholds(i);
        adjGrp = connectomeGroupThreshold(WEIGHT, grpThr, 2);
        Adj = adjGrp(1:100, 1:100);
        Adj(isnan(Adj))=0;
        

            for j=1:length(weights)
                
                % threshold matrix heeping a set proportion of nonzero weights
                perc = weights(j);
                [AdjThr] = KeepPercWeights(Adj, perc);
                mask = logical(AdjThr);
                [y,e, dataCell,t, p] = CoexpConnectVSUnconnect(Exp,mask);
                Tvals(i,j) = t;
                Pvals(i,j) = p;
                ConUncon(j,:,l) = y;
                ConUnconSD(j,:,l) = e;
                
                
            end

    end
    
    figure; imagesc(Tvals); colormap([flipud(BF_getcmap('blues',20));1,1,1;BF_getcmap('reds',20)]);
    caxis([-7 7]);colorbar; %colorbar('southoutside');
    title(sprintf('%dDS Corrected T values', DS)); xlabel('Weight threshold'); ylabel('Consistency between subjects, %');
    
    
end
