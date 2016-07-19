
cd ('/Users/Aurina/Documents/Genetics_connectome/Gen_Cog/Data/Connectomes/Final88/');

% load anatomical parcellation
Tract = {'FACT_'};
parc = {'aparcANDaseg'};

Name = strcat(Tract{1}, parc{1});
load (sprintf('%s.mat', Name));

NumNodesAnat = size(density{1},1);
NumSubj = size(density,2);

nodesAnat= zeros(NumSubj,NumNodesAnat,3);

for i=1:NumSubj 
    for j=1:NumNodesAnat
        nodesAnat(i,j,:) = coords{1,i}(j,:); 
    end
end

meanCoordsAnat = squeeze(mean(nodesAnat,1));

% load random parcellation

parc = {'custom200ANDaseg20'};
Name = strcat(Tract{1}, parc{1});

load (sprintf('%s.mat', Name));

NumNodesRand = size(density{1},1);
nodesRand = zeros(NumSubj,NumNodesRand,3);

for i=1:NumSubj 
    for j=1:NumNodesRand 
        nodesRand(i,j,:) = coords{1,i}(j,:); 
    end
end

meanCoordsRand = squeeze(mean(nodesRand,1));

% for each node in random parcellation calculate the losest point in
% anatomical parcellation according to the average coordinates
k = dsearchn(meanCoordsAnat,meanCoordsRand);
numNodes = length(k);
Q = 1:1:numNodes;
list = [Q',k];

% reorder list of random parcellation according to anatomica parcellation
p = sortrows(list,2);

%reorder a random matrix according to new node order
densityReordered = cell(1, NumSubj);
countReordered = cell(1, NumSubj);

for l=1:NumSubj
    densityReordered{l} = density{l}(p(:,1), p(:,1));
    countReordered{l} = count{l}(p(:,1), p(:,1));
end
