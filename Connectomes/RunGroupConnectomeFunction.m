%% this script creates group connecromes for different thresholds and saves them in this format:

% Adj_all(threshold,NumNodes, NumNodes)
% first dimention in a matrix corresponds to the threshold (connections present in % of subjects): 
% 0 (all connections), 30, 50, 60, 75, 90 % respectively. 
% note: fGroupMatrix_commonConnections function removes low density matrices before calculating mst common connections. 

clear all; close all;
% go to where the matrices are
cd ('/Users/Aurina/Documents/Genetics_connectome/Gen_Cog/Data/Connectomes/Final88/');
% choose the type of tractography (FACT or iFOD2)
Tract = {'FACT_'};
% choose the parcellations
Parcellation = {'aparcANDaseg', 'custom200ANDaseg20'};

for parc = Parcellation

    Name = strcat(Tract{1}, parc{1});
    % load a file
    load (sprintf('%s.mat', Name));
    % generate group connectomme using count and density as weights separately
    Weights = {count, density};
    
    % make a group connectome for different thresholds: 0 (all connections), 30, 50, 60, 75, 90% 
    % First dimention in Adj_all corresponds to those thresholds in this order. 
    Adj_all = fGroupMatrix_commonConnections(Weights,Tract, parc);
    
end