

%data = CombinedExp(:,3:end);
data = DataNORM; 
%data = BF_NormalizeMatrix(data,'zscore' );
%data = Combined(:,2:30281);
%data = B; 
% data = DATA;
Subj = CombinedExp(:,1);
[W, pc] = pca(data); pc=pc'; W=W';


figure; 
for i=1:length(Subj)
    if Subj(i) == 1; 
        scatter(pc(1,i),pc(2,i),'.', 'k');hold on; 
    elseif Subj(i) == 2; 
        scatter(pc(1,i),pc(2,i),'.', 'r');hold on; 
    elseif Subj(i) == 3; 
        scatter(pc(1,i),pc(2,i),'.', 'g');hold on; 
    elseif Subj(i) == 4; 
        scatter(pc(1,i),pc(2,i),'.', 'm');hold on; 
    elseif Subj(i) == 5; 
        scatter(pc(1,i),pc(2,i),'.', 'b'); hold on; 
    elseif Subj(i) == 6; 
        scatter(pc(1,i),pc(2,i),'.', 'c'); hold on; 
    end
end