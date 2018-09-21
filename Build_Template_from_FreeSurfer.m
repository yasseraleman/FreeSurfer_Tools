function OutTempFile = Build_Template_from_FreeSurfer(FreeSDir,Idlist, TemplateFilename);
%
% Syntax :
% OutTempFile = Build_Template_from_FreeSurfer(Idlist, TemplateFilename);
%
% This script uses ANTS subroutines to create anatomical templates for 
% preprocessed freesurfer data.
%
% Input Parameters:
%   Images             : Native T1 Images
%   TempFilename       : Template Filename
%  
%
% Output Parameters:
%
%     OutTempFile      : Output Template Filename
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
Images = '';
Ids = textread(Idlist,'%s');
Ids = char(Ids);
Ns = size(Ids,1);
for i = 1:Ns
    Id = deblank(Ids(i,:));
    cad = ['mri_convert -i ' FreeSDir filesep Id filesep 'mri' filesep 'T1.mgz' ' -o ' FreeSDir filesep Id filesep 'tmp' filesep Id '_T1.nii'];
    system(cad);
    outputs = flip_to_axial_nii([FreeSDir filesep Id filesep 'tmp' filesep Id '_T1.nii'],1);
    Images = strvcat(Images,outputs);
end

OutTempFile = Build_Template(Images, TemplateFilename);