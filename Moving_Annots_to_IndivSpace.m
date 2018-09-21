function varargout = Moving_Annots_to_IndivSpace(FreeSDir, Idfile, varargin);
%
% Syntax :
%     OutAnnotFiles = Moving_Annots_to_IndivSpace(FreeSDir, Idfile, varargin);
%
% This function moves individual annot files to individual space.
%
% Input Parameters:
%       FreeSDir                : FreeSurfer Subjects Directory.
%       Idfile                  : Ids File
%       parcType                : Parcellation type (default = lausanne).
%                                 This script do not work for parcellations
%                                 already included in freesurfer.
%        refId                  : Reference Id (default = fsaverage).
%
% Output Parameters:
%        OutAnnotFile           : Saving the sulci parcellation in an
%                                 annotation file. If nargin <3
%                                 OutAnnotFile is a vector file containing
%                                 sulci parcellation.
%
% Examples of usage: 
%   OutAnnotFiles = Moving_Annots_to_IndivSpace('/home/user/myFreesurfersubjects', '/home/user/myFreesurfersubjects/Ids.txt');
%   OutAnnotFiles = Moving_Annots_to_IndivSpace('/home/user/myFreesurfersubjects', '/home/user/myFreesurfersubjects/Ids.txt','refId','fsaverage');
%   OutAnnotFiles = Moving_Annots_to_IndivSpace('/home/yaleman/PROCESSING_RESULTS/5-freesurfer_processing', '/home/yaleman/PROCESSING_RESULTS/IDs.txt','refId','fsaverage','parcType','myparc' );
%
%   Files lh.myparc.annot and rh.myparc.annot must be stored in the label
%   folder inside the freesurfer folder for the reference Id.
% 
% See also: save_annotfiles
%__________________________________________________
% Authors: Yasser Aleman Gomez
% CIBM, MIAL, Radiology and Psychiatry Departments, CHUV
% January 31th 2018
% Version $1.0


%% ====================== Checking input parameters ===================== %
if nargin<2 % the indispensable input arguments are not provided
    error('Two inputs are mandatory');
else
    if ~exist(FreeSDir,'dir')
        error('Wrong Freesurfer Subjects  Directory');
        return;
    end
    if length(strfind(Idfile,filesep))
        if ~exist(Idfile,'file')
            error('Wrong Freesurfer Subjects  Directory');
            return;
        end
    end
    
    % Default parcellation
    parcType = 'lausanne';
    
    % Reference Id
    refId = 'fsaverage';
end


% deal with the input arguments
if nargin<1 % the indispensable input arguments are not provided
    error('One input is mandatory');
else
    if numel(varargin)>0 % optional input arguments are provided
        while ~isempty(varargin)
            if numel(varargin)<2
                error('You need to provide optional input arguments as ''ParameterName''-''ParameterValue'' pairs.');
            end
            switch varargin{1}
                case 'parcType' % Parcellation type 
                    parcType=varargin{2};
                case 'refId' % Reference subject
                    refId=varargin{2};
                otherwise
                    error('Unexpected ''ParameterName'' input: %s\n',varargin{1});
            end
            varargin(1:2)=[]; % this pair of optional input arguments has been dealt with -- remove...
        end
    end
end
%% ===================== End of Checking input parameters =============== %

% FreeSDir = '/home/yaleman/PROCESSING_RESULTS/5-freesurfer_processing';
% Idfile = '/home/yaleman/PROCESSING_RESULTS/5-freesurfer_processing/Giedre_Ids.txt';

% Selecting the cortical parcellation. Lausanne parcellation includes the
% parcellation of both cortical hemispheres in 5 different scales
switch lower(parcType)
    case 'lausanne'
        annotIds = strvcat('lausanne2008.scale1','lausanne2008.scale2','lausanne2008.scale3','lausanne2008.scale4','lausanne2008.scale5');
    case {'aparc','aparc.a2009s','aparc.a2005s'} % These parcellations are already included in the freesurfer output folder.
        annotIds = '';
    otherwise
        annotIds = parcType;
end

% Reading the Ids file
if exist(Idfile,'file')
    Ids = char(textread(Idfile,'%s'));
else
    Ids = Idfile;
end


%% ======================= Main Program ================================ %%
% Freesurfer IDs for subcortical structures and brain stem
left_subcIds = [10 11 12 13 26 17 18];
right_subcIds = [49 50 51 52 58 53 54];
brain_stem = 16;

if ~isempty(annotIds)
    OutAnnotFiles = '';
    Ns = size(Ids,1); % Number of subjects
    Nannots = size(annotIds,1); % Number of annots
    parfor  i = 1:Ns % Loop over subjects
        Id = deblank(Ids(i,:));
        for j = 1:Nannots % Loop over annots
            annotId = deblank(annotIds(j,:));
            
            % Moving Left annot to individual space
            lhannot = [FreeSDir filesep refId filesep 'label' filesep 'lh.' annotId '.annot'];   % Annot in reference space (commonly fsaverage)
            outannotfile = [FreeSDir filesep Id filesep 'label' filesep 'lh.' annotId '.annot']; % Annot in individual space
            cad = ['mri_surf2surf --srcsubject ' refId  ' --trgsubject ' Id ' --hemi lh --sval-annot ' lhannot ' --tval ' outannotfile];
            system(cad);
            
            % Moving Right annot to individual space
            rhannot = [FreeSDir filesep refId filesep 'label' filesep 'rh.' annotId '.annot'];   % Annot in reference space (commonly fsaverage)
            outannotfile = [FreeSDir filesep Id filesep 'label' filesep 'rh.' annotId '.annot']; % Annot in individual space
            cad = ['mri_surf2surf --srcsubject ' refId ' --trgsubject ' Id ' --hemi rh --sval-annot ' rhannot ' --tval ' outannotfile];
            system(cad);
            disp(' ')
            
            OutAnnotFiles = strvcat(OutAnnotFiles,[FreeSDir filesep Id filesep 'label' filesep 'lh.' annotId '.annot'],[FreeSDir filesep Id filesep 'label' filesep 'rh.' annotId '.annot']);
            
            % Creating the volume in conform space (freesurfer)
            switch parcType
                case 'lausanne'
                    cad = ['mri_aparc2aseg --s ' Id ' --annot ' annotId ' --wmparc-dmax 0 --labelwm --hypo-as-wm --o ' [FreeSDir filesep Id filesep 'tmp' filesep annotId '+aseg.nii']];
                    system(cad);
                    
                    V = spm_vol([FreeSDir filesep Id filesep 'tmp' filesep annotId '+aseg.nii']);
                    I = spm_read_vols(V);
                    It = I*0;
                    
                    % Relabelling Right hemisphere
                    ind = find(I > 2000 & I <3000);
                    It(ind) = I(ind) - 2000;
                    nlabel = max(It(:));
                    
                    % Relabelling Subcortical Right hemisphere
                    newLabels = [nlabel+1:nlabel+length(right_subcIds)];
                    [bool,ord] = ismember(I,right_subcIds);
                    ind = find(bool);
                    It(ind) = newLabels(ord(ind));
                    nlabel = max(It(:));
                    
                    % Relabelling Left hemisphere
                    ind = find(I > 1000 & I <2000);
                    It(ind) = I(ind) - 1000 + nlabel;
                    nlabel = max(It(:));
                    
                    % Relabelling Subcortical left hemisphere
                    newLabels = [nlabel+1:nlabel+length(left_subcIds)];
                    [bool,ord] = ismember(I,left_subcIds);
                    ind = find(bool);
                    It(ind) = newLabels(ord(ind));
                    nlabel = max(It(:));
                    
                    % Relabelling Brain Stem
                    ind = find(I == brain_stem);
                    It(ind) = nlabel+1;
                    
                    % Saving the new parcellation
                    Vout = V;
                    Vout.fname = [FreeSDir filesep Id filesep 'tmp' filesep annotId '+aseg.nii'];
                    spm_write_vol(Vout,It);
                    
                    % Converting the new parcellation to mgz
                    cad = ['mri_convert -i '  [FreeSDir filesep Id filesep 'tmp' filesep annotId '+aseg.nii'] ' -o ' [FreeSDir filesep Id filesep 'mri' filesep annotId '+aseg.mgz']];
                    system(cad)
                    
                    % Removing temporary files
                    delete([FreeSDir filesep Id filesep 'tmp' filesep annotId '+aseg.nii']);
                otherwise
                    cad = ['mri_aparc2aseg --s ' Id ' --annot ' annotId ' --wmparc-dmax 0 --labelwm --hypo-as-wm --o ' [FreeSDir filesep Id filesep 'mri' filesep annotId '+aseg.mgz']];
                    system(cad);
            end
        end
    end
end
%% ==================== End of Main Program ============================ %%
% Outputs
varargout{1} = OutAnnotFiles;
return;