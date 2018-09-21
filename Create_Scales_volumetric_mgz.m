function NNames = Create_Scales_volumetric_mgz(FreeSDir, IdsFile, scales);
%
% Syntax :
% Create_Scales_volumetric_mgz(FreeSDir, IdsFile, scales);
%
% This function creates volumetric mgz files from different scales parcellation
%
% Input Parameters:
%      FreeSDir         : FreeSurfer Directory
%      IdsFile          : Ids File
%      scales           : Scales identification
%
% Output Parameters:
%
%
% Related references:
%
%
% See also:
%__________________________________________________
% Authors: Yasser Aleman Gomez
% LIM
% September 12th 2012
% Version $1.0

scales = strvcat('aparc_scale03_Nzones0270','aparc_scale04_Nzones0364','aparc_scale05_Nzones0410');
Ns = size(scales,1); 
IdsFile = '/media/COSAS/scripts/Pipeline_Manage/Ids.txt';
opts.pipe.freesdir = '/media/Data/PROCESSING_RESULTS/PEPS/5-freesurfer_processing';
Ids = textread(IdsFile,'%s');
%Ids{1} = '0002-20090515';
Nsub = size(Ids,1);
NNames = '';
for i = 1:Ns
    aparcId = deblank(scales(i,:));
    for j = 1:Nsub
        opts.pipe.subjId = char(Ids{j});
        disp(['Computing Subject ' num2str(j) ' of ' num2str(Nsub) ': Scales ' aparcId]);
        cad = ['mri_convert -i ' opts.pipe.freesdir filesep opts.pipe.subjId filesep 'mri' filesep 'ribbon.mgz' ' -o ' opts.pipe.freesdir filesep opts.pipe.subjId filesep 'tmp' filesep 'ribbon.nii'];
        system(cad);
        V = spm_vol([opts.pipe.freesdir filesep opts.pipe.subjId filesep 'tmp' filesep 'ribbon.nii']);
        I = spm_read_vols(V);
        It = I*0;
        ind = find(I == 3);
        It(ind) = 1;
        ind = find(I == 42);
        It(ind) = 1;
        VA = spm_vol([opts.pipe.freesdir filesep opts.pipe.subjId filesep 'tmp' filesep aparcId '+aseg+1mm.nii']);
        IA = spm_read_vols(VA);
        IA(IA<1000) = 0;
        It = It.*IA;
        
        Nname = [opts.pipe.freesdir filesep opts.pipe.subjId filesep 'tmp' filesep aparcId '+aseg.nii'];
        cad = ['mri_convert -i ' opts.pipe.freesdir filesep opts.pipe.subjId filesep 'mri' filesep 'aparc+aseg.mgz' ' -o ' opts.pipe.freesdir filesep opts.pipe.subjId filesep 'tmp' filesep aparcId '+aseg.nii'];
        system(cad);
        VA = spm_vol(Nname);
        IA = spm_read_vols(VA);
        ind = find(It);
        IA(ind) = It(ind);
        spm_write_vol(VA,IA);
        cad = ['mri_convert -i ' opts.pipe.freesdir filesep opts.pipe.subjId filesep 'tmp' filesep aparcId '+aseg.nii' ' -o ' opts.pipe.freesdir filesep opts.pipe.subjId filesep 'mri' filesep aparcId '+aseg.mgz'];
        system(cad);
        delete(Nname);
        delete([opts.pipe.freesdir filesep opts.pipe.subjId filesep 'tmp' filesep 'ribbon.nii']);
        NNames = strvcat(NNames,Nname);
     end
end
return;