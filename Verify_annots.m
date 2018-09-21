function VolAtlasfiles = Verify_annots;
%
% Syntax :
% VolAtlasfiles = Verify_annots;
%
% This script creates different parcellations from a known annot file.
% The parcellations can be organized in multiple scales. For each scale the
% parcellation obtained in the previos scale is taken as initial
% parcellation for the next scale.
%
%
% Input Parameters:
%    optst                      : Structure variable containing the options
%                                 for freesurfer average subject
%    optst.pipe.freesdir        : FreeSurfer Directory
%    optst.pipe.subjId          : Subject Id
%    optst.anat.atlas.nzones    : Number of zones per hemisphere in the fine
%                                 scale
%    optst.anat.atlas.grow      : Gray matter growing option. Distance in mm (default = 1 mm)
%
%    opts                       : Structure variable containing the options
%                                 for individual subjects
%    opts.pipe.freesdir        : FreeSurfer Directory
%    opts.pipe.subjId          : Subject Id
%    opts.anat.atlas.grow      : Gray matter growing option. Distance in mm (default = 1 mm)
%     Nscales                   : Number of scales (default = 1)
%
% Output Parameters:
%     VolAtlasfiles            : Volumetric parcellations
%
% See also:
%__________________________________________________
% Authors: Yasser Aleman Gomez
% LIM, HUGGM
% March 22th 2013
% Version $1.0

clc;


opts.pipe.outdir = '/media/Data/PROCESSING_RESULTS/PEPS';
opts.pipe.freesdir = [opts.pipe.outdir filesep '5-freesurfer_processing'];
setenv('SUBJECTS_DIR',opts.pipe.freesdir);
opts.anat.atlas.corrdiff = 1;
annotIds = strvcat('aparc_scale01_Nzones0097','aparc_scale02_Nzones0191','aparc_scale03_Nzones0270','aparc_scale04_Nzones0364','aparc_scale05_Nzones0410');


%% ============ End of Setting Enviroment for Freesurfer Average ======= %%

%% ================== Subjects Processing ============================== %%
IdsFile = '/media/COSAS/scripts/Pipeline_Manage/Ids.txt';
IdsFile = '/media/Data/PROCESSING_RESULTS/REDES/5-freesurfer_processing/redesids.txt';
opts.pipe.freesdir  = '/media/Data/PROCESSING_RESULTS/REDES/5-freesurfer_processing/';
Ids = char(textread(IdsFile,'%s'));
FailedFiles = '';
FailedFilesV = '';
cont = 0;
Ns = size(annotIds,1);
Nsubj = size(Ids,1);
 a = 0;
for i =1:Ns
    annotId = deblank(annotIds(i,:));
    for j= 1:Nsubj
        Id = deblank(Ids(j,:));
        tic;
        disp(['Processing =======>  Scale: ' num2str(i) ' of ' num2str(Ns) ' . =====> Subject ID: ' Id ' . ---  ' num2str(j) ' of ' num2str(Nsubj)]);
        %% =========== Moving Surfaces to individual space ============= %%
                    cadl = [opts.pipe.freesdir  filesep Id filesep 'label' filesep 'lh.' annotId '.annot'];
                    cadr = [opts.pipe.freesdir  filesep Id filesep 'label' filesep 'rh.' annotId '.annot'];
                    [txtl,ctabl] = read_cfiles(cadl);
                    [txtr,ctabr] = read_cfiles(cadr);
                    try
                        
                    if (sum(unique(txtl) - sort(ctabl.table(:,5))) == 0) & (sum(unique(txtr) - sort(ctabr.table(:,5))) == 0) & (sum(unique(txtl) - unique(txtr)== 0))
        
                    else
                        FailedFiles = strvcat(FailedFiles, [Id ' ====== > ' annotId]);
                    end
                    catch
                         FailedFiles = strvcat(FailedFiles, [Id ' ====== > ' annotId]);
                    end
%         atlas = ['/media/Data/PROCESSING_RESULTS/PEPS/7-connectome/' Id '/preproc/' Id '_DiffparcKmeans_' annotId '+1mm.nii'];
%         % atlas = ['/media/Data/PROCESSING_RESULTS/PEPS/5-freesurfer_processing/' Id '/tmp/' annotId '+aseg+1mm.nii'];
%         V = spm_vol(atlas);
%         I = spm_read_vols(V);
%         ind = find(I == 5001);
%         I(ind) = 0;
%         ind = find(I == 5002);
%         I(ind) = 0;
%         indl = find((I>1000)&(I<2000));
%         stl = unique(I(indl));
%         indr = find((I>2000)&(I<3000));
%         str = unique(I(indr));
%         if length(stl) ~= length(str)
%             FailedFilesV = strvcat(FailedFilesV, [Id ' ====== > ' annotId]);
%             a = 1;
%         else
%             if (sum((unique(str) - unique(stl))-1000) == 0)
%                 
%             else
%                 FailedFilesV = strvcat(FailedFilesV, [Id ' ====== > ' annotId]);
%                 a = 1;
%             end
%         end
%         if a == 1
%             opts.pipe.outdir = '/media/Data/PROCESSING_RESULTS/PEPS';
%             Conf_File = [opts.pipe.outdir filesep '7-connectome' filesep Id filesep 'logs' filesep Id '_Out_ConfFile.txt'];
%             opts = Reading_Pipe_Configuration_File(Conf_File);
%             Outatlas = [opts.pipe.outdir filesep '7-connectome' filesep opts.pipe.subjId filesep 'preproc' filesep opts.pipe.subjId '_DiffparcKmeans_' annotId '+' num2str(opts.anat.atlas.grow) 'mm.nii'];
%             try
%                 cad = ['WarpImageMultiTransform 3 ' ufilename ' ' Outatlas ' -R ' opts.diff.pdata.b0 ' ' opts.transforms.t1_2_diff.affine ' --use-NN'];
%                 system(cad);
%             catch
%                 FailedFiles = strvcat(FailedFiles,Outatlas);
%             end
%         end
%         a = 0;
        
    end
end
disp(FailedFiles)
return;
