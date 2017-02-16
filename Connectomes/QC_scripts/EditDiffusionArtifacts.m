
%clear all; 
cd ('/gpfs/M2Scratch/Monash076/aurina/Gen_Cog/code/Diffusion_artefacts/Correction_scripts');

%load('CorrectionList.mat');

subjects = correct(:,1);
SubjIDs = unique(subjects);

for k = 1:length(SubjIDs)
    subj = SubjIDs(k);
    ind = find(subjects == subj);
    [d,h] = OpenDiffImage(subj);
    for l=1:length(ind)
        j = ind(l);
            % % average adjacent for fault slices input z value, gradients value
            % % choose what gradient value
            
            gradient = correct(j,2);
            z =  correct(j,3);
            
            [d] = Averageslices(gradient, z, d);
       
    end

   SaveDiffImg (subj, h,d );
end








