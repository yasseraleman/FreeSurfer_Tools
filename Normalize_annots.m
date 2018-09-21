function varargout = Normalize_annots(varargin);
%
% Syntax :
%     OutAnnotFiles = Normalize_annots(FreeSDir, Idfile, trgId);
%
% This function moves individual annot files to a target sapace.
% See
%
% Input Parameters:
%       FreeSDir                : FreeSurfer Subjects Directory.
%       Idfile                  : Ids File
%       trgId                   : Target Id
%
% Output Parameters:
%        OutAnnotFile           : Saving the sulci parcellation in an
%                                 annotation file. If nargin <3
%                                 OutAnnotFile is a vector file containing
%                                 sulci parcellation.
%
% See also: save_annotfiles
%__________________________________________________
% Authors: Yasser Aleman Gomez 
% LIM, HUGGM
% November 13th 2014
% Version $1.0

FreeSDir = '/media/Data/PROCESSING_RESULTS/5-freesurfer_processing';
%Idfile = '/media/MyDisk/PROCESSING_RESULTS/5-freesurfer_processing/HCPrest_Ids.txt';
Idfile = '/media/MyDisk/PROCESSING_RESULTS/5-freesurfer_processing/HCPwmPR1_Ids.txt';
trgId = 'fsaverage';
if exist(Idfile,'file')
    Ids = char(textread(Idfile,'%s'));
else
    Ids = Idfile;
end

Ns = size(Ids,1);
failed = '';
cad = '';
Ido = '';
for  i = 1:Ns
    
    Id = deblank(Ids(i,:));
    
    lhannot = [FreeSDir filesep Id filesep 'label' filesep 'lh.sulcbasins.aparc2sulci.annot'];
    outannotfile = [FreeSDir filesep Id filesep 'label' filesep 'lh.fsaverage.aparc2sulci.annot'];
    cad = ['mri_surf2surf --srcsubject ' Id ' --trgsubject ' trgId ' --hemi lh --sval-annot ' lhannot ' --tval ' outannotfile];
    system(cad);
    
    lhannot = [FreeSDir filesep Id filesep 'label' filesep 'rh.sulcbasins.aparc2sulci.annot'];
    outannotfile = [FreeSDir filesep Id filesep 'label' filesep 'rh.fsaverage.aparc2sulci.annot'];
    cad = ['mri_surf2surf --srcsubject ' Id ' --trgsubject ' trgId ' --hemi rh --sval-annot ' lhannot ' --tval ' outannotfile];
    system(cad);
    disp(' ')
end
return;

