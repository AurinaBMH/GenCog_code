function [y,e, dataCell,t, p] = CoexpConnectVSUnconnect(Exp,Adj)

% connected vs unconnected
% Rexp - coexpression matrix
% Adj - adjecency matrix

%% output
% y - mean values for connected and unconnected links
% e - SD for connected and unconnected links
% dataCell{1,1} - distribution of connected values
% dataCell{2,1} - - distribution of unconnected values
% t - ttest2 t value
% p - ttest2 p value


mask = logical(Adj);

        
       Con = nonzeros(Exp(mask==1));
       Connected = nanmean(Con);
       ConnectedSD = nanstd(Con);
       
       Uncon = nonzeros(Exp(mask==0));
       Unconnected = nanmean(Uncon);
       UnconnectedSD = nanstd(Uncon);
       
       y = [Connected Unconnected];
       e = [ConnectedSD UnconnectedSD];
        
       dataCell = cell(2,1);
       dataCell{1,1} = Con; dataCell{2,1} = Uncon; 
       
%JitteredParallelScatter(dataCell)
%ylabel('Average coexpression'); legend('Connecteced', 'Unconnected');
[h,p,ci,stats] = ttest2(Con,Uncon, 'Vartype','unequal');
t = stats.tstat;
end
