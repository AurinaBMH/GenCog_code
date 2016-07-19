
ExpSample = cell(6,1);
CoordSample = cell(6,1); 
Parcellation = {'cust100'};
NormMethod = {'hampel'};

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
    
    ExpSubj1 = DataExpression{i,1};
    Coord1 = DataCoordinates{i,1};
    
    ExpSubj = ExpSubj1((ExpSubj1(:,2)<=LeftCortex),:);
    Coord = Coord1((Coord1(:,5)<=LeftCortex),2:4);
    
    if strcmp(NormMethod, 'hampel')
            ExpSample{i,1} = Norm_hampel(ExpSubj(:,3:size(ExpSubj,2)));
    else
            ExpSample{i,1} = BF_NormalizeMatrix(ExpSubj(:,3:size(ExpSubj,2)),NormMethod{1});
    end
 
   CoordSample{i,1} = Coord;
   
   
end

Samples = cat(1, ExpSample{1},  ExpSample{2},  ExpSample{3},  ExpSample{4},  ExpSample{5},  ExpSample{6});
CoordSample = cat(1, CoordSample{1},  CoordSample{2},  CoordSample{3},  CoordSample{4},  CoordSample{5},  CoordSample{6});


SamplesCorr = corr(Samples');
MRIvoxCoordinates = pdist2(CoordSample, CoordSample);
 Dvect = triu(MRIvoxCoordinates,1);
 Rvect = triu(SamplesCorr,1);
  L = logical(Dvect);
  Rvect = Rvect.*L;
 Dvect = nonzeros(Dvect(:));
 Rvect = nonzeros(Rvect(:));
 
%make a vector for coexpression and distances


%Plot quantiles
 BF_PlotQuantiles(Dvect,Rvect,100,0,1)
 figure; imagesc(SamplesCorr); caxis([-1,1])
colormap([flipud(BF_getcmap('blues',9));BF_getcmap('reds',9)]);
        
% Data = BF_NormalizeMatrix(SamplesRaw, 'scaledRobustSigmoid');
% figure; imagesc(Data);
% Dist = pdist2(Data, Data);
% figure; imagesc(Dist);
% d_row = pdist(Dist,'corr');
% links_row = linkage(d_row);
% [~,~,ord_row] = dendrogram(links_row,0); 
% 
% % Subj = Combined(:,2); 
% % ROI = Combined(:,1); 
% 
% figure; subplot(1,25,1); imagesc(Subjects(ord_row)); 
%         subplot(1,25,2); imagesc(ROIS(ord_row)); 
%         subplot(1,25, [3 25]);  imagesc(Data(ord_row,:));