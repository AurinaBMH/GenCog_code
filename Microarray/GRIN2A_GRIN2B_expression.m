% script to get GRIN2A and GRIN2B gene expression for Nigel
% Aurina Arnatkeviciute 14/02/2017

%% run a part of DS code to get left cortex normalised data assigned to anatomical parcellation
clear all; close all;

Parcellation = {'aparcaseg'};
Threshold = 2;
NormMethod = {'zscore'};
LEFTcortex = 1;
% choose 1 if want to normalise samples assigned to left cortex separately;
% choose 2 if want to normalise LEFT cortex + left subcortex together
% choose 3 if you want to normalise the whole brain.
% choose 4 if you want to normalise left cortex + right cortex.
Thr = 0;

% choose 1 if want to normalise samples assigned to left cortex separately; choose 0 if want to normalise all samples together.
if LEFTcortex == 1 || LEFTcortex == 2
    NumSubjects = 6;
elseif LEFTcortex == 3 || LEFTcortex == 4
    NumSubjects = 2;
end

if LEFTcortex == 1
    numR = 34; NAME = 'NLC';
elseif LEFTcortex==2
    numR = 41; NAME = 'NLCS';
end

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
cd ('/Users/Aurina/GoogleDrive/Genetics_connectome/Gen_Cog/Data/Microarray/S01-S06_combined/');
load(sprintf('%d_DistThresh%d_%s_combined_ExpressionProbePCA_GeneThr%d.mat', NumNodes, Threshold, Parcellation{1}, round(Thr)));


ExpressionSubjROI = cell(6,1);
CoordinatesSubjROI = cell(6,1);
CoordSample = cell(6,1);
ExpSampNorm = cell(6,1);

for i=1:NumSubjects
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
    CoordSample{i} = Coord;
    % normalise sample x gene data for each subject separately
    %% commented - not normalise
    if strcmp(NormMethod, 'hampel')
        ExpressionSampleNorm = Norm_hampel(data);
    else
        ExpressionSampleNorm = BF_NormalizeMatrix(data, NormMethod{1});
    end
    
    ROI = ExpSubj(:,2);
    
    ExpSampNorm{i} = [ROI, ExpressionSampleNorm];
    ROIs = unique(ExpSubj(:,2));
    
    % combine noramlised data for all subjects.
end
ExpSampNormalisedAll = vertcat(ExpSampNorm{1}, ExpSampNorm{2},ExpSampNorm{3},ExpSampNorm{4},ExpSampNorm{5},ExpSampNorm{6});
CombinedCoord = cat(1,CoordSample{1}, CoordSample{2}, CoordSample{3},...
    CoordSample{4}, CoordSample{5}, CoordSample{6});

% get probeIDs for selected genes
cd('/Users/Aurina/GoogleDrive/Genetics_connectome/Gen_Cog/Data/Microarray/');
load('MicroarrayDataPCAS01.mat');
Probes = ProbeInformation.GeneSymbol;

%% get data for GRIN2A and GRIN2B genes
[~, GRIN2Aind] = intersect(Probes, 'GRIN2A');
[~, GRIN2Bind] = intersect(Probes, 'GRIN2B');

%% get values for those two genes for each brain separately
% create a cell to store data
expROISUBJECTS = cell(6,2); 
for pp=1:6
    expression = ExpSampNorm{pp,1}(:,2:end);
    ROIS = ExpSampNorm{pp,1}(:,1);
    expression2genes = expression(:,[GRIN2Aind GRIN2Bind]);
    
    %% average samples according to ROIs
    numROIs = length(unique(ROIS));
    uROIS = unique(ROIS);
    expROI = zeros(numROIs,2);
    aa=1;
    for kk=1:numROIs
        R = uROIS(kk);
        ROIind = find(ROIS==R);
        expROI(aa,:) = mean(expression2genes(ROIind,:));
        aa=aa+1;
    end
    
    cd('/Users/Aurina/GoogleDrive/Genetics_connectome/Gen_Cog/Data/Microarray/')
    % get ROInames
    fid = fopen('ROInames.txt');
    ROInames = textscan(fid, '%s%s','Delimiter','\t');
    fclose(fid);
    ROIsLEFT = ROInames{1,1}(1:numR,:);
    ROIsLEFT = ROIsLEFT(uROIS);
    
    expROISUBJECTS{pp,1} = expROI;
    expROISUBJECTS{pp,2} = ROIsLEFT; 
    
    %% plot expression for those 2 genes
    h = figure; imagesc(expROI); caxis([-2,2]);title(sprintf('GRIN2A, GRIN2B gene expression in subject %d', pp));
    colormap([flipud(BF_getcmap('blues',9));[1 1 1]; BF_getcmap('reds',9)]); colorbar;
    gene1 = Probes(GRIN2Aind);
    gene2 = Probes(GRIN2Bind);
    genes = [gene1,gene2];
    set(gca,'xtick',[1 2], 'ytick', 1:length(uROIS), 'xticklabel',genes, 'yticklabel',ROIsLEFT, 'FontSize',15);
    cd('/Users/Aurina/GoogleDrive/Genetics_connectome/Gen_Cog/Data/Microarray/NigelEXPRESSION')
    print(h,sprintf('%s_GRIN2A_GRIN2B_expression_subject_%d.jpg', NAME,pp),'-dpng');
    savefig(h,sprintf('%s_GRIN2A_GRIN2B_expression_subject_%d.fig', NAME,pp)); 
end

%% get expression for those two genes average among all subjects
expROIAVERAGE = cell(1,2); 
expressionALL = ExpSampNormalisedAll(:,2:end);
ROIS = ExpSampNormalisedAll(:,1);
expression2genes = expressionALL(:,[GRIN2Aind GRIN2Bind]);

%% average samples according to ROIs
numROIs = length(unique(ROIS));
uROIS = unique(ROIS);
expROI = zeros(numROIs,2);
aa=1;
for kk=1:numROIs
    R = uROIS(kk);
    ROIind = find(ROIS==R);
    expROI(aa,:) = mean(expression2genes(ROIind,:));
    aa=aa+1;
end
cd('/Users/Aurina/GoogleDrive/Genetics_connectome/Gen_Cog/Data/Microarray/')
% get ROInames
fid = fopen('ROInames.txt');
ROInames = textscan(fid, '%s%s','Delimiter','\t');
fclose(fid);
ROIsLEFT = ROInames{1,1}(1:numR,:);
ROIsLEFT = ROIsLEFT(uROIS);
expROIAVERAGE{1,1} = expROI; 
expROIAVERAGE{1,2} = ROIsLEFT; 
%% plot expression for those 2 genes
h = figure; imagesc(expROI); caxis([-1.5,1.5]);title('AVERAGE GRIN2A, GRIN2B gene expression');
colormap([flipud(BF_getcmap('blues',9));[1 1 1]; BF_getcmap('reds',9)]); colorbar;
genes = {'GRIN2A', 'GRIN2B'};
set(gca,'xtick',[1 2], 'ytick', 1:length(uROIS), 'xticklabel',genes, 'yticklabel',ROIsLEFT, 'FontSize',15);
cd('/Users/Aurina/GoogleDrive/Genetics_connectome/Gen_Cog/Data/Microarray/NigelEXPRESSION')
print(h, sprintf('%s_GRIN2A_GRIN2B_expression_AVERAGE.jpg',NAME)   ,'-dpng');
savefig(h, sprintf('%s_GRIN2A_GRIN2B_expression_AVERAGE.fig',NAME)); 

save(sprintf('%s_GRIN2A_GRIN2B_data',NAME),'expROISUBJECTS','expROIAVERAGE','genes'); 