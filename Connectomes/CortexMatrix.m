
function Adj = CortexMatrix(Parcel, Adj)
Parcel = {'aparcaseg'};

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


%Adj = squeeze(Adj_all(4,:,:));
L = Adj(1:LeftCortex, 1:LeftCortex);
R = Adj(LeftSubcortex+1:RightCortex, LeftSubcortex+1:RightCortex);
RL1 = Adj(1:LeftCortex, LeftSubcortex+1:RightCortex); RL2 = Adj(LeftSubcortex+1:RightCortex, 1:LeftCortex);
U = horzcat(L, RL1); L = horzcat(RL2, R); M = vertcat(U,L);
Adj = M; 



[r,c, v] = find(Adj);
edges = [r,c, v];

