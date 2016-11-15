
function [grad] = Removebvecs(subj, ind)
ind = ind+1;
cd(['/gpfs/M2Scratch/Monash076/simon/GenCog/subjects/1008.2.57.' num2str(subj) '/diffusion/']);
%read original files
grad = dlmread('dwscheme_orig.bvec');
bval = dlmread('dwscheme_orig.bval'); 

grad = grad';
bval = bval';
for i=ind
grad(ind,:)=NaN;
bval(ind,:) = NaN;
end
grad(~any(~isnan(grad),2),:)=[];
bval(~any(~isnan(bval),2),:)=[];

cd(['/gpfs/M2Scratch/Monash076/simon/GenCog/subjects/1008.2.57.' num2str(subj) '/diffusion/']);
%cd (sprintf('/gpfs/M2Scratch/Monash076/aurina/Gen_Cog/code/Diffusion_artefacts/AfterArtifactRemovalNEW43/%d', subj));
%save corrected files
dlmwrite('dwscheme_qc.bvec', grad');
dlmwrite('dwscheme_qc.bval', bval');




end

%% edit bvecs

