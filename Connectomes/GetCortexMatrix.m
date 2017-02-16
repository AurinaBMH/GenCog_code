
function Adj = GetCortexMatrix(Adj, Parcel)


if strcmp(Parcel, 'aparcANDaseg_FA') || strcmp(Parcel, 'aparcANDaseg')
    
    NumNodes = 82;
    LeftCortex = 34;
    LeftSubcortex = 41;
    RightCortex = 75;
    RightSubcortex = NumNodes;
    
elseif strcmp(Parcel, 'custom200ANDaseg20_FA') || strcmp(Parcel, 'custom200ANDaseg20') 
    NumNodes = 220;
    LeftCortex = 100;
    LeftSubcortex = 110;
    RightCortex = 210;
    RightSubcortex = NumNodes;
    
elseif strcmp(Parcel, 'custom200ANDaseg_FA') || strcmp(Parcel, 'custom200ANDaseg')
    NumNodes = 214;
    LeftCortex = 100;
    LeftSubcortex = 107;
    RightCortex = 207;
    RightSubcortex = NumNodes;
    
elseif strcmp(Parcel, 'custom200ANDfslatlas_FA') || strcmp(Parcel, 'custom200ANDfslatlas')
    NumNodes = 228;
    LeftCortex = 100;
    LeftSubcortex = 114;
    RightCortex = 214;
    RightSubcortex = NumNodes;
    
elseif strcmp(Parcel, 'cust250')
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
end

