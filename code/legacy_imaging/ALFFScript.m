% HISTORICAL REFERENCE ONLY. Requires original inputs/toolboxes and path configuration.
% Use code/RUN_ALL.m for the portable release analyses.
clear

% rp_Cohe_ReHo_Brain('F:\Mito_Rest\FunImgNormalized\WITHHELD_SUBJECT\', 27, 'LOCAL_SOURCE\Documents\GitHub\CanlabCore\CanlabCore\canlab_canonical_brains\Canonical_brains_surfaces\brainmask_canlab.nii', 'F:\Mito_Rest\Result_Reho',0.46,0.01,0.1,0,782,1)


% data=fmri_data('F:\Mito_Rest\FunImgNormalized\WITHHELD_SUBJECT\sub-WITHHELD_SUBJECT_ses-PicardMiSBIE_task-rest_run-1_space-MNI152NLin2009cAsym_desc-preproc_bold.nii')
% Mask=fmri_data('LOCAL_SOURCE\Documents\GitHub\CanlabCore\CanlabCore\canlab_canonical_brains\Canonical_brains_surfaces\brainmask_canlab.nii')
% data_resample=resample_space(Mask,data);
% 
% fname=strcat('Resampled_CanlabMask.nii');
% write(data_resample, 'fname', fname, 'overwrite');

% reho('F:\Mito_Rest\FunImgNormalized\WITHHELD_SUBJECT\', 27, 'LOCAL_SOURCE\Documents\GitHub\CanlabCore\CanlabCore\canlab_canonical_brains\Canonical_brains_surfaces\brainmask_canlab.nii', 'F:\Mito_Rest\Result_Reho')
AllVolume='F:\Mito_Rest\FunImgNormalized\WITHHELD_SUBJECT\';
% AllVolume='F:\Mito_Preprocessed_Onemodel\Rest\sub-WITHHELD_SUBJECT_ses-PicardMiSBIE_task-rest_NuissanceRegress.nii';
ASamplePeriod=0.46;
HighCutoff=0.1; LowCutoff=0.01;
AMaskFilename='E:\Mito_DICOM\Code\Resampled_CanlabMask.nii'; 
AResultFilename_Reho='F:\Mito_Rest\Result\Reho';
AResultFilename_ALFF='F:\Mito_Rest\Result\ALFF';
NVoxel=27;
TemporalMask=''; 
% ScrubbingMethod, 
% Header;
% CUTNUMBER;

 [ALFFBrain, fALFFBrain, Header] = y_alff_falff(AllVolume,ASamplePeriod, HighCutoff, LowCutoff, AMaskFilename, AResultFilename_ALFF, TemporalMask)


 Band=[0.01 0.08];
 IsNeedDetrend=1;
 TR=0.46;
[ReHoBrain, Header] = y_reho(AllVolume, NVoxel, AMaskFilename, AResultFilename_Reho, IsNeedDetrend, Band, TR, TemporalMask)

