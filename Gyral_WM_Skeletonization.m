function [OutFile] = Gyral_WM_Skeletonization(Id);
%
% Syntax :
% [OutFile, cads] = Gyral_WM_Skeletonization(FreeSDir, Id, OutAtlasFile);
%
% This script skeletonize the fa image and label the FA inside each gyri
% using the FreeSurfer cortical parcellation.
%
% Input Parameters:
%   OutputDir         : Pipeline Output directory 
%   Id                : Subject Id
%  
%
% Output Parameters:
%
%     OutFile         : Labeled skeleton
%
% Related references:
%
%
% See also: 
% 
%
%__________________________________________________
% Authors: Yasser Aleman Gomez
% LIM, HUGGM
% July 27th 2013
% Version $1.0

% OutputDir = '/media/Data/PROCESSING_RESULTS/PEPS';
% OutAtlasFile = '/media/Data/Joost/gspan/volumes/gyri/0019x-20091120.wmparc.gyri.nii';

FreeSDir = '/media/Data/PROCESSING_RESULTS/5-freesurfer_processing';
AtlasId = 'wmparc';
 
cad = ['mri_convert ' FreeSDir filesep Id filesep 'mri' filesep  AtlasId '.mgz ' FreeSDir filesep Id filesep 'tmp' filesep AtlasId '.nii'];
system(cad);

AtlasFile = [FreeSDir filesep Id filesep 'tmp' filesep AtlasId '.nii'];

V = spm_vol(AtlasFile);
I = spm_read_vols(V);
ind = find(I <3000);
I(ind) = 0;
V.fname =  [FreeSDir filesep Id filesep 'tmp' filesep Id '_wm.nii'];
spm_write_vol(V,logical(I));

 cad = ['fslmaths ' V.fname ' -s 2 ' FreeSDir filesep Id filesep 'tmp' filesep Id '_swm.nii'];
 system(cad);
try
    cad = ['gunzip -df ' FreeSDir filesep Id filesep 'tmp' filesep Id '_swm.nii.gz'];
    system(cad);
end

cad = ['tbss_skeleton -i ' FreeSDir filesep Id filesep 'tmp' filesep Id '_swm.nii' ' -o ' FreeSDir filesep Id filesep 'tmp' filesep Id '_swm_skel.nii'];
system(cad);
try
    cad = ['gunzip -df ' FreeSDir filesep Id filesep 'tmp' filesep Id '_swm_skel.nii.gz'];
    system(cad);
end

V = spm_vol([ FreeSDir filesep Id filesep 'tmp' filesep Id '_swm_skel.nii']);
I = spm_read_vols(V);
VA = spm_vol(AtlasFile);
IA = spm_read_vols(VA);
IA = logical(I).*IA;
ind = find(IA<3000);
IA(ind) = 0;
ind = find(IA>5000);
IA(ind) = 0;
VA.fname =  [FreeSDir filesep Id filesep 'mri' filesep Id '_skelwmatlas.nii'];
spm_write_vol(VA,IA);
delete([FreeSDir filesep Id filesep 'tmp' filesep Id '_wm.nii']);
delete([FreeSDir filesep Id filesep 'tmp' filesep Id '_swm.nii.gz']);
delete([ FreeSDir filesep Id filesep 'tmp' filesep Id '_swm_skel.nii']);
delete(AtlasFile);

OutFile = [FreeSDir filesep Id filesep 'mri' filesep Id '_skelwmatlas.nii'];
cad = ['mri_convert ' OutFile ' ' FreeSDir filesep Id filesep 'mri' filesep Id '_skelwmatlas.mgz'];
system(cad);
delete(OutFile);
OutFile = [FreeSDir filesep Id filesep 'mri' filesep Id '_skelwmatlas.mgz'];

return;
