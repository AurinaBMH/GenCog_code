clear all; close all;

Parcellation = {'cust100'}; 
Threshold = 2; 
%NormMethod = {'zscore'};
LEFTcortex = 4; 
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
cd ('/Users/Aurina/Documents/Genetics_connectome/Gen_Cog/Data/Microarray/S01-S06_combined/');
load(sprintf('%d_DistThresh%d_%s_combined_ExpressionProbePCA_GeneThr3.mat', NumNodes, Threshold, Parcellation{1}));

 ExpressionSubj = cell(6,1);
 CoordinatesSubj = cell(6,1);
 
for i=1:2
  % normalise data for each subject separately using samples  
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
    
    data = ExpSubj(:,3:size(ExpSubj,2));
    
    
      ROI = ExpSubj(:,2);
      S = ExpSubj(:,1); 
      
      ExpressionSubj{i} = [S, ROI, data];
      CoordinatesSubj{i} = [S, ROI, Coord];
end

CombinedExp = cat(1,ExpressionSubj{1}, ExpressionSubj{2}, ExpressionSubj{3}, ExpressionSubj{4}, ExpressionSubj{5}, ExpressionSubj{6});
CombinedCoord = cat(1,CoordinatesSubj{1}, CoordinatesSubj{2}, CoordinatesSubj{3}, CoordinatesSubj{4}, CoordinatesSubj{5}, CoordinatesSubj{6});

% normalise all samples together
DataNORM = BF_NormalizeMatrix(CombinedExp(:,3:end), 'zscore');
%% plot normalised combined data clustered according to provile similarity. 

DATA = DataNORM;
Dist = pdist2(DATA, DATA);
d_row = pdist(Dist,'corr');
links_row = linkage(d_row);
[~,~,ord_row] = dendrogram(links_row,0); 

Subj = CombinedExp(:,1); 
ROI = CombinedExp(:,2); 


figure; subplot(1,25,1); imagesc(Subj(ord_row)); 
        subplot(1,25,2); imagesc(ROI(ord_row)); 
        subplot(1,25, [3 25]);  imagesc(DATA(ord_row,:)); caxis([-3 3]); title(sprintf('LeftCortex %s', 'z-score'));
        colormap([flipud(BF_getcmap('blues',9));BF_getcmap('reds',9)]);