% differences between 10 and 50 mln streamlines
% clear all; close all; 
cd ('/Users/Aurina/GoogleDrive/Genetics_connectome/Gen_Cog/Data/Connectomes/10_50_mlnStreamlines/');
load('FACT_custom500ANDaseg30_10_mil.mat');
density10 = density; 
len10 = len;
count10 = count;
deg10 = zeros(length(density10), size(density10{1},1));
NumNodes = size(density{1},1);
for i=1:length(density10)
    deg10(i,:) = degrees_und(density10{i}); 
end
load('FACT_custom500ANDaseg30_50_mil.mat');
density50 = density; 
len50 = len;
count50 = count;
deg50 = zeros(length(density), size(density{1},1));
for i=1:length(density)
    deg50(i,:) = degrees_und(density50{i}); 
end
dif = deg50-deg10;
for i=1:length(density)
    FigHandle = figure('Position', [100, 100, 1200, 1200]);
    subplot(3,2,5);
    scatter(deg10(i,:), deg50(i,:));
    xlabel('10 mil streamlines degrees'); ylabel('50 mil streamlines degrees');
    hold on; 
    x = linspace(1,size(density{1},1), size(density{1},1));
    y = x;
    plot(x,y);
     
    subplot(3,2,6); bar(x,dif(i,:)); title('Difference between degrees (50mln - 10mln streamlines)');
    xlabel('Node');
    ylabel('Difference in degrees'); 
    Matrix10 = fRandom2AnatomicalParcellations_annot('custom500', 'ANDaseg30', density10{i});
    Matrix50 = fRandom2AnatomicalParcellations_annot('custom500', 'ANDaseg30', density50{i});
    subplot(3,2,1); imagesc(log(Matrix10)); axis square; title ('Connectome 10 mil streamlines');
    subplot(3,2,2); imagesc(log(Matrix50)); axis square; title ('Connectome 50 mil streamlines');
   
    subplot(3,2,3); histogram(deg10(i,:), 20); title('Degree distribution 10 mil'); xlabel('Degree');
    subplot(3,2,4); histogram(deg50(i,:), 20); title('Degree distribution 50 mil'); xlabel('Degree');

end
A = mean(dif,1); 
figure; bar(x,A);

Ranks10 = zeros(length(density),NumNodes); 
Ranks50 = zeros(length(density),NumNodes); 

for i=1:length(density)
    [ignore, idx10] = sort( deg10(i,:) , 'descend');
    [ignore, idx50] = sort( deg50(i,:) , 'descend');
    Ranks10(i,idx10 ) = 1:numel(idx10);
    Ranks50(i,idx50 ) = 1:numel(idx50);
end

R = zeros(length(density),NumNodes); 
for i=1:length(density)
R10 = Ranks10(i,:); 
R50 = Ranks50(i,:);
Ra = R50-R10;
R(i,:) = Ra;
figure; 
%FigHandle = figure('Position', [100, 100, 1200, 1200]);
subplot(1,2,1); bar(Ra); title(sprintf('%d subj Difference between ranks 50mln - 10mln\ncustom500 parcellation', sub_id(i))) ;xlabel('Node'); ylabel('Rank difference'); axis square;
subplot(1,2,2); scatter(R10, R50); title(sprintf('%d subj 50mln and 10mln rank correlation\ncustom500 parcellation', sub_id(i))); xlabel('10 mln rank'); ylabel('50 mln rank');axis square;
hold on; 
x = linspace(1,size(density{1},1), size(density{1},1));
y = x;
plot(x,y);
end
%mean rank difference with SD for each ROI across subjects
MeanR = mean(R,1);
SDR = std(R,1);

% weight and distance correlation
for i=1:length(density)
    Weights10 = density10{i}(:);
    Length10 = len10{i}(:);
    Weights50 = density50{i}(:);
    Length50 = len50{i}(:);
    A10 = [Length10 Weights10];
    A50 = [Length50 Weights50];
    A10( ~any(A10,2), : ) = [];  %rows
    A50( ~any(A50,2), : ) = [];  %rows
    BF_PlotQuantiles(A10(:,1),log10(A10(:,2)),500,0,1)
    title (sprintf('Subject %d 10 mln custom500ANDaseg30', sub_id(i))); xlabel('Connection length'); ylabel('log10 (Connection weight (count))');
    BF_PlotQuantiles(A10(:,1),log10(A10(:,2)),500,0,1)
    title (sprintf('Subject %d 50 mln custom500ANDaseg30', sub_id(i))); xlabel('Connection length'); ylabel('log10 (Connection weight (count))');
end