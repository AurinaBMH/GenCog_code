
clear all; 
close all; 

cd ('/Users/Aurina/Documents/Genetics_connectome/Gen_Cog/Data/Connectomes/FinalCorrections/');
Parcellation = {'cust100'}; % choose parcellation
Tract = {'FACT_'}; % choose tract method

Name = strcat(Tract{1}, Parcellation{1});
load (sprintf('%s.mat', Name));

Density = zeros(size(density,2),1);

for i=1:size(density,2)
    Mask = count{i};
    density{i}(Mask<10) = 0;
    figure; imagesc(log(density{i})); title(sprintf('Subject %d', sub_id{i}));
    Deg = degrees_und(density{i});
    Density(i) = density_und(density{i});
    Strength = sum(density{i});
    Len = len{i};
     figure; 
     subplot(2,1,1);histogram(Deg, 20); title(sprintf('Degree Subject %d', sub_id{i}));
     subplot(2,1,2);histogram(nonzeros(Len(:)), 20); title(sprintf('length Subject %d', sub_id{i}));
end


% cd ('/Users/Aurina/Documents/Genetics_connectome/Gen_Cog/Data/Connectomes/WithACT/');
% Parcellation = {'aparcaseg'};
% Tract = {'FACT_'};
% 
% Name = strcat(Tract{1}, Parcellation{1});
% load (sprintf('%s.mat', Name));
% 
% for i=1:size(density,2)
%     figure; imagesc(log(density{i})); title(sprintf('Subject %d', sub_id{i}));
%     Deg = degrees_und(density{i});
%     Len = length{i};
%     figure; 
%     subplot(2,1,1);histogram(Deg, 20); title(sprintf('degrees Subject %d', sub_id{i}));
%     subplot(2,1,2);histogram(nonzeros(Len(:)), 20); title(sprintf('length Subject %d', sub_id{i}));
% end
% 














