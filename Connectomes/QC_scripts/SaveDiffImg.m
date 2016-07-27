
function SaveDiffImg (subj, h,d)
%cd(['/gpfs/M2Scratch/Monash076/aurina/Gen_Cog/dwidenoise/1008.2.57.' num2str(subj) ]);

cd(['/gpfs/M2Scratch/Monash076/simon/GenCog/subjects/1008.2.57.' num2str(subj) '/diffusion/']);
% cd ('/gpfs/M2Scratch/Monash076/aurina/Gen_Cog/code/Diffusion_artefacts/AfterArtifactRemoval/');
% if exist(sprintf('%d', subj), 'dir')==0
% mkdir(sprintf('%d', subj));
% end
% cd (sprintf('/gpfs/M2Scratch/Monash076/aurina/Gen_Cog/code/Diffusion_artefacts/AfterArtifactRemoval/%d', subj));
%write(h,d,'forward_qc.nii');

write(h,d,'forward_qcraw.nii');
end