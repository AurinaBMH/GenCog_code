
clear all; 
close all; 

cd ('/Users/Aurina/Documents/Genetics_connectome/Gen_Cog/Data/Connectomes/WithoutACT/'); %choose the folder
Parcellation = {'aparcaseg'}; % choose parcellation
Tract = {'iFOD2_'}; % choose tract method

Name = strcat(Tract{1}, Parcellation{1});
load (sprintf('%s.mat', Name));

Density = zeros(size(density,2),1);

for i=1:size(density,2)
    figure; imagesc(log(density{i})); title(sprintf('Subject %d', sub_id{i}));
   
    Strength = sum(density{i});
    
    figure; histogram(Strength, 20); title(sprintf('Strength Subject %d', sub_id{i}));
     
end
% clear all; 
% close all; 
% 
% cd ('/Users/Aurina/Documents/Genetics_connectome/Gen_Cog/Data/Connectomes/WithACT/'); %choose the folder
% Parcellation = {'aparcaseg'}; % choose parcellation
% Tract = {'iFOD2_'}; % choose tract method
% 
% Name = strcat(Tract{1}, Parcellation{1});
% load (sprintf('%s.mat', Name));
% 
% Density = zeros(size(density,2),1);
% 
% for i=1:size(density,2)
%     figure; imagesc(log(density{i})); title(sprintf('Subject %d', sub_id{i}));
%    
%     Strength = sum(density{i});
%     
%     figure; histogram(Strength, 20); title(sprintf('Strength Subject %d', sub_id{i}));
%      
% end