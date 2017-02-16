cd('/Users/Aurina/GoogleDrive/Genetics_connectome/Gen_Cog/Data/Connectomes/Denoise')
load('FACT_custom200ANDaseg_withdenoise.mat');
With = density; 
load('FACT_custom200ANDaseg_withoutdenoise.mat');
Without = density;

Adjall = cell(2,1);
AdjCort = cell(2,1);

% make group connectomes
[Adjall{1,1}] = connectomeGroupThreshold(With, 0.6, 2) ;
[Adjall{2,1}] = connectomeGroupThreshold(Without, 0.6, 2) ;

%take cortex
AdjCort{1,1} = Adjall{1,1}(1:100,1:100);
AdjCort{2,1} = Adjall{2,1}(1:100,1:100);


Inds = cell(2,1);
klevel = max(Degree); 
Hub = 20; 


for i=1:2

Degree = degrees_und(AdjCort{i,1}); 
[~, ind] = find(Degree > Hub); 
Inds{i,1} = ind;

end

Overlap = intersect(Inds{1,1}, Inds{2,1}); 
% 
% for k = Huball:klevel
%         
%         isHub = Degree > k;
%         
%         
% %         mask(isHub,isHub) = 1; % rich
% %         mask(isHub,~isHub) = 2; % feeder
% %         mask(~isHub,isHub) = 2;
% %         
% end
% figure; imagesc(mask); title('