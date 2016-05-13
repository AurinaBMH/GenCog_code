

clear all; close all; 

DataExpressionS = cell(6,1);
DataCoordinatesS = cell(6,1);
Parcellation = {'cust250'};
Threshold = 2; 
Thr = 0;
NormMethod = {'robustSigmoid'};
LEFTcortex = 1; 
load (sprintf('82_DistThresh2_aparcaseg_combined_ExpressionProbePCA_GeneThr%d.mat', Thr));

% choose 1 if want to normalise samples assigned to left cortex separately; 
% choose 2 if want to normalise LEFT cortex + left subcortex together
% choose 3 if you want to normalise the whole brain. 
% choose 4 if you want to normalise left cortex + right cortex. 

if strcmp(Parcellation, 'aparcaseg')
     
                NumNodes = 82;
                LeftCortex = 34;
                LeftSubcortex = 41; 
                RightCortex = 75;
                RightSubcortex = NumNodes;
    elseif strcmp(Parcellation, 'cust100')
                NumNodes = 220;
                LeftCortex = 100;
                LeftSubcortex = 110; 
                RightCortex = 210;
                RightSubcortex = NumNodes;
    elseif strcmp(Parcellation, 'cust250')
                NumNodes = 530;
                LeftCortex = 250;
                LeftSubcortex = 265; 
                RightCortex = 515;
                RightSubcortex = NumNodes;
end

for i=1:6
  % normalise data for each subject separately using ROIs
    ExpSubj1 = DataExpression{i,1};
    Coord1 = DataCoordinates{i,1};
    
    
    if LEFTcortex == 1
        ExpSubj = ExpSubj1((ExpSubj1(:,2)<=LeftCortex),:);
        Coord = Coord1((Coord1(:,5)<=LeftCortex),2:4);
    elseif LEFTcortex == 2
        ExpSubj = ExpSubj1((ExpSubj1(:,2)<=LeftSubcortex),:);
        Coord = Coord1((Coord1(:,5)<=LeftSubcortex),2:4);
    elseif LEFTcortex == 3
        ExpSubj = ExpSubj1;
        Coord = Coord1(:,2:4);
    elseif LEFTcortex == 4
        ExpSubj3 = ExpSubj1((ExpSubj1(:,2)>=LeftSubcortex & ExpSubj1(:,2)<=RightCortex),:);
        Coord3 = Coord1((Coord1(:,5)>=LeftSubcortex & Coord1(:,5)<=RightCortex),2:4);
        
        ExpSubj2 = ExpSubj1((ExpSubj1(:,2)<=LeftCortex),:);
        Coord2 = Coord1((Coord1(:,5)<=LeftCortex),2:4);
        
        ExpSubj = cat(1, ExpSubj2, ExpSubj3);
        Coord = cat(1, Coord2, Coord3); 
    end
    Inf = ExpSubj(:,1:2);
    DataExpressionNorm = BF_NormalizeMatrix(ExpSubj(:,3:size(ExpSubj,2)),NormMethod{1}) ;
    DataExpressionS{i} = [Inf DataExpressionNorm];
    DataCoordinatesS{i} = [Inf Coord];
    

end


CombinedExp = cat(1,DataExpressionS{1}, DataExpressionS{2}, DataExpressionS{3}, ...
    DataExpressionS{4}, DataExpressionS{5}, DataExpressionS{6});
CombinedCoord = cat(1,DataCoordinatesS{1}, DataCoordinatesS{2}, DataCoordinatesS{3},...
    DataCoordinatesS{4}, DataCoordinatesS{5}, DataCoordinatesS{6});




DATA = CombinedExp(:,3:size(ExpSubj,2)); 
Dist = pdist2(DATA, DATA);
d_row = pdist(Dist,'corr');
links_row = linkage(d_row);
[~,~,ord_row] = dendrogram(links_row,0); 

ROI = CombinedExp(:,2); 
Subj = CombinedExp(:,1); 

figure; subplot(1,25,1); imagesc(Subj(ord_row)); 
        subplot(1,25,2); imagesc(ROI(ord_row)); 
        subplot(1,25, [3 25]);  imagesc(DATA(ord_row,:)); caxis([0 1]); title(sprintf('LeftCortex %s', NormMethod{1}));

%%
ExpressionNorm = corr(DATA');
Coordinates = CombinedCoord(:,3:5);
MRIvoxCoordinates = pdist2(Coordinates, Coordinates);
 Dvect = triu(MRIvoxCoordinates,1);
 Rvect = triu(ExpressionNorm,1);
L = logical(Dvect);
Rvect = Rvect.*L;
 Dvect = nonzeros(Dvect(:));
 Rvect = nonzeros(Rvect(:));
 
%make a vector for coexpression and distances


%Plot quantiles
 BF_PlotQuantiles(Dvect,Rvect,100,0,1)
 figure; imagesc(ExpressionNorm); caxis([-1,1])
colormap([flipud(BF_getcmap('blues',9));BF_getcmap('reds',9)]);
        

Fit = {'exp_1_0'};
[f_handle,Stats,c] = GiveMeFit(Dvect,Rvect,Fit{1});

figure; plot(c, Dvect, Rvect); 
figure; plot(c,Dvect, Rvect,'Residuals')

switch Fit{1}
case 'linear'
    FitCurve = c.p1*Dvect + c.p2;
case 'exp'
    FitCurve = c.A*exp(-c.n*Dvect) + c.B;
case 'exp_1_0'
    FitCurve = exp(-c.n*Dvect);
case 'decay'
    FitCurve = c.A/Dvect + c.B;
    Residuals = Rvect' - FitCurve;
case 'exp0'
    FitCurve = c.A.*exp(-c.n*Dvect);
case 'exp1'
    FitCurve = exp(-c.n*Dvect) + c.B;
end
    
Residuals = Rvect - FitCurve;
