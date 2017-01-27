
% run Differential_stability.m script first

% this script fit coexpression-distance relationship separately for subsets of data

IND1 = find(DistExpVect(:,1)>0); 
IND2 = find(DistExpVect(:,1)<=120); 
IND = intersect(IND1,IND2);

i=1; 
S = struct;
P=struct; 
R= struct; 

while length(IND)>1
    
if length(IND)>2000
INDrun = randsample(IND,2000); 
else
    INDrun=IND; 
end

%[param,stat] = sigm_fit(DistExpVect(INDrun,1), DistExpVect(INDrun,2)); 
[curve, goodness, output] = fit(DistExpVect(INDrun,1), DistExpVect(INDrun,2), 'exp1');
name = sprintf('run%d',i);  
%S.(name) = stat; 
%P.(name) = param; 
%Residuals = DistExpVect(INDrun,2) - stat.ypred;
coords = DistExpVect(INDrun,1); 
R.(name) = [coords, DistExpVect(INDrun,2), curve, output.residuals]; 

IND = setdiff(IND, INDrun);
i=i+1; 
end


for i=20:30
name = sprintf('run%d',i);
figure; scatter(R.(name)(:,1), R.(name)(:,3),'MarkerEdgeColor',[0 .5 .5], 'MarkerFaceColor',[0 .7 .7],'LineWidth',1.5);
xlabel('Distance between samples, mm'); 
ylabel('Coexpression'); 
end