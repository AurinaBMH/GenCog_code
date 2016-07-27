% Run MakeCopy

cd('/gpfs/M2Scratch/Monash076/aurina/Gen_Cog/code/Diffusion_artefacts/Correction_scripts');
% make copy for subjects that needed slice corrections but no grad
% corrections. 
%load('Subjects.mat');
%run separately for NoCorrections variable and for OnlySliceCorrections variable)
% make copy for subjects that didn't need corrections


for i=1:length(Correct)
    subj = Correct(i);
    
    %open original image and save it with _qc
    cd(['/gpfs/M2Scratch/Monash076/simon/GenCog/subjects/1008.2.57.' num2str(subj) '/diffusion/']);
%     grad = dlmread('dwscheme_orig.bvec');
%     bval = dlmread('dwscheme_orig.bval');
%     
%     dlmwrite('dwscheme_qc.bvec', grad);
%     dlmwrite('dwscheme_qc.bval', bval);
% 
%     
      [h, d] = read('forward.nii');
      write(h,d,'forward_qcraw.nii');
    
end