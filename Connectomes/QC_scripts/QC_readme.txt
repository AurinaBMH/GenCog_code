Steps for diffusion image QC:

1. Manual inspection of all images and artifact registration (all corrections saved in CorrectionList.mat).
1.1. for subjects that don't need any corrections (Subjects(NoCorrections)) create a copy of original file named forward_qcraw.nii and make a copy of bvec and bval files and name them with "_qc" at the end.
1.2. for subjects that need slice dropout correction (but no gradient correction), make a copy of bvec and bval files with "_qc" at the end (Subjects(OnlySliceCorrections)).
Use MakeCopyforOKsubjects.m script for it. 

2. Averaging adjecent slices for slice drop-outs using EditDiffusionArtifacts.m
	it loads correction list (format in each row: subject, volume, slice)
	opens forward.nii file and averages slices
	saves a forward_qcraw.nii file

3. for subjects that had corrections before, bad gradients are removed using forward_qcraw.nii file
   for subjects that had no corrections before (Subjects(GradCorrections) - 38,82,166,191,204), bad gradients are removed using forward.nii file. (run separately for those two cases)
   in both cases save them as forward_qcraw.nii.
	- change files in RemoveSingleVolume script.

	RunRemoveSingleVolume - removes single volume.
	RemoveMultipleVolumes - removes multiple volumes for 5 subjects (44,51,116,144,280).

4. Update bvec and bval files according to removed volumes using RunRemovebvecs.m
	
