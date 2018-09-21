function OutFiles = Aparcmgz_To_aparcnii(FreesurferDir);
%FreesurferDir = '/media/Data/PROCESSING_RESULTS/PEPS/5-freesurfer_processing';
a = dir([FreesurferDir filesep '*']);
[namet{1:size(a,1)}] = deal(a.name);namet = char(namet);
[direc{1:size(a,1)}] = deal(a.isdir);direc = cell2mat(direc);
SubjIds = namet(direc,:);
SubjIds = SubjIds(3:end,:);
%% =============== Verifying if they are freesurfer results ============ %%
a = zeros(size(SubjIds,1),1);
T1Files = '';
for i = 1:size(SubjIds,1)
    if exist([FreesurferDir filesep deblank(SubjIds(i,:)) filesep 'mri' ],'dir');
        if exist([FreesurferDir filesep deblank(SubjIds(i,:)) filesep 'mri' filesep 'aparc+aseg.mgz' ],'file')
            T1Files = strvcat(T1Files,[FreesurferDir filesep deblank(SubjIds(i,:)) filesep 'mri' filesep 'aparc+aseg.mgz' ]);
            a(i) =1;
        end
    end
end
SubjIds = SubjIds(find(a),:);
Ns = size(SubjIds,1);
for i = 1:1
    disp(['Processing Subject ' num2str(i) ' of ' num2str(Ns)]);
    cad = ['mri_convert -i ' FreesurferDir filesep deblank(SubjIds(i,:)) filesep 'mri' filesep 'aparc+aseg.mgz' ' -o ' FreesurferDir filesep deblank(SubjIds(i,:)) filesep 'tmp' filesep 'aparc+aseg.nii'];
    system(cad);
    outputs = flip_to_axial_nii([FreesurferDir filesep deblank(SubjIds(i,:)) filesep 'tmp' filesep 'aparc+aseg.nii'],0);
    delete([FreesurferDir filesep deblank(SubjIds(i,:)) filesep 'tmp' filesep 'aparc+aseg.nii']);
    movefile(outputs,[FreesurferDir filesep deblank(SubjIds(i,:)) filesep 'tmp' filesep 'aparc+aseg.nii']);
end
return;