
function [d,h] = OpenDiffImage(subj)

cd(['/gpfs/M2Scratch/Monash076/simon/GenCog/subjects/1008.2.57.' num2str(subj) '/diffusion/']);
% Read the diffusion images file
%cd(['/gpfs/M2Scratch/Monash076/aurina/Gen_Cog/dwidenoise/1008.2.57.' num2str(subj) ]);
[h, d] = read('forward.nii');
% save original diffusion file with name forward_beforeqc.nii
%write(h,d,'forward_beforeqc.nii');
% bvals = dlmread('dwscheme_corrected.bval')'; 
% grad = dlmread('dwscheme_corrected.txt')'; 
% inds = linspace(1, length(grad),length(grad));
end