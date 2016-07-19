
clear all; close all; 


Lcortex = cell(6,1);
Thr = 0;

for i=1:6
    FileName = sprintf('MicroarrayDataPCAS0%d_GeneThr%d_leftCortex.mat', i, Thr);
    load(FileName);
    Lcortex{i} = Expression;
end


%t = zeros(6,6,size(Lcortex{1,1},2));
t=zeros(6,6);
Score2gene = zeros(size(Lcortex{1,1},2),1);

for j=1:size(Lcortex{1,1},2)
    
    for i=1:6
        
        for k=1:6
        
        GeneS1 = Lcortex{i,1}(:,j);
        GeneS2 = Lcortex{k,1}(:,j);
        [~,~,ks2stat] = kstest2(GeneS1,GeneS2);
        t(i,k) = ks2stat;
 
        end
    end
        tu = triu(t,1);
        tu = (nonzeros(tu)); 
        Score2gene(j)= mean(tu);

end

