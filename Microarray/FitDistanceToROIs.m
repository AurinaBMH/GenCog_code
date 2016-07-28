% this scripts performes distance coexpression-distance fitting on average
% ROI values (ROI-ROI coexpression and average euclidean distances between
% ROIs - ROI coordinate for each subject is calculated according to the
% mean sample coordinates (assigned to the same ROI), then the average of
% those (between subjects) is taken. 
clear all; close all; 

load('TestUncorrectedROI.mat'); 

LeftCortex = 100;
A = vertcat(CoordinatesSubjROI{1}, CoordinatesSubjROI{2}, CoordinatesSubjROI{3}, ...
    CoordinatesSubjROI{4}, CoordinatesSubjROI{5}, CoordinatesSubjROI{6});


B = sortrows(A,1);
B = B(:,[1,3:5]);
ROI = B(:,1);
W = unique(ROI);

for i=1:length(W)
    
    C = find(ROI == W(i));
    if size(C,1)~=1
    MeanROI(i,:) = mean(B(C,:));
    else
    MeanROI(i,:) = B(C,:);
    end
        
    
end

ROIcoor = MeanROI(:,2:4);
ROIdist = pdist2(ROIcoor, ROIcoor); 

V(:,1) = ROIdist(:);
V(:,2) = ParcelCoexpression(:);

%V( ~any(V,2), : ) = [];  %rows

[curve, goodness, output] = fit(V(:,1),V(:,2),'poly2');
figure; scatter(V(:,1), V(:,2)); hold on; plot(curve); 

O = reshape(output.residuals,[99 99]);

p = 10;
Exp = O;
N1 = nan(LeftCortex,1);
N2 = nan(1, LeftCortex-1);


B = vertcat(Exp(1:p-1,:), N2, Exp(p:end,:));
Exp = horzcat(B(:,1:p-1), N1, B(:,p:end));


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
    caxis([0 12]);colorbar; %colorbar('southoutside');
    title(sprintf('%dDS Corrected T values', DS)); xlabel('Weight threshold'); ylabel('Consistency between subjects, %');
    
    
end
