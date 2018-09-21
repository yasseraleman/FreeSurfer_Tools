function varargout = Correct_aparc_and_lobar_parcellation(varargin);
%
% Syntax :
%     OutLobesAnnots = Correct_aparc_and_lobar_parcellation(FreeSDir, Idfile);
%
% This function does an approximate mapping of individual 'Desikan-Killiany' ROIs 
% (found in ?h.aparc.annot) to the lobes and corrects aparc parcellation
% files.
%
% Input Parameters:
%        FreeSDir               : Freesurfer subjects directory.
%        Idfile                 : Subjects Ids text file. It can be also a
%                                 char matlab variable (vertical concatenation 
%                                 of subjects Ids.
%
% Output Parameters:
%        OutAnnotFiles          : Saving the lobar parcellation in an
%                                 annotation file.
%                                 OutAnnotFiles is a vector file containing
%                                 lobar parcellation.
%
% See also: save_annotfiles
%__________________________________________________
% Authors: Yasser Aleman Gomez 
% LIM, HUGGM
% November 13th 2014
% Version $1.0

%% ============================= Checking Inputs ======================= %%
% % % % % % % Idfile = '/media/Data/PROCESSING_RESULTS/HCP/5-freesurfer_processing/temp1.txt';
% % % % % % % FreeSDir = '/media/Data/PROCESSING_RESULTS/HCP/5-freesurfer_processing';
FreeSDir = varargin{1};
Idfile   = varargin{2};
if nargin > 2
    errordlg('To Many Input Parameters');
    return;
end
if nargout > 1
    errordlg('To Many Output Parameters');
    return;
end
%% ========================= End of Checking Inputs ==================== %%

%% ============================= Reading Ids  ========================== %%

if exist(Idfile,'file')
    Ids = char(textread(Idfile,'%s'));
else
    Ids = Idfile;
end
Ns = size(Ids,1);
%% ======================== End of Reading Ids  ======================== %%

%% =========================== Main Program  =========================== %%
OutAnnotFiles = '';
for  i = 1:Ns
    Id = deblank(Ids(i,:));
    disp(['Processing Subject: ' Id '. No ' num2str(i) ' of ' num2str(Ns)] );
    
    %% ======================== Left Hemisphere ========================= %
    % Reading aparc annotation file
    AnnotFile = [FreeSDir filesep Id filesep filesep 'label' filesep 'lh.aparc.annot'];
    [txt,colortable] = read_cfiles(AnnotFile);
    
    % Reading Surface
    Surftemp = Read_Surface([FreeSDir filesep Id filesep 'surf' filesep 'lh.white']);
    
    % Detecting 0s in label values
    ind = find(txt == 0);
    txt(ind) = 1;
    Surftemp.Is = txt;
    
    % Removing isolated clusters and refilling holes
    Surf = Surf_Corr(Surftemp);
    ind = find(Surf.Is == 1);
    txt(ind) = 0;
    Surf.Is(ind) = 0;
    
    % Saving the correct annot file
    OutAnnotFile = save_annotfiles(Surf.Is,AnnotFile,colortable);
    
    % Creating Lobar parcellation the correct annot file
    OutAnnotFile = [FreeSDir filesep Id filesep filesep 'label' filesep 'lh.aparc.lobes.annot'];
    OutAnnotFile = Aparc2Lobes(AnnotFile, 0, OutAnnotFile);
    OutAnnotFiles = strvcat(OutAnnotFiles, OutAnnotFile);
    %% =================== End of Left Hemisphere ======================= %
    
    %% ======================== Right Hemisphere ======================== %
    % Reading aparc annotation file
    AnnotFile = [FreeSDir filesep Id filesep filesep 'label' filesep 'rh.aparc.annot'];
    [txt,colortable] = read_cfiles(AnnotFile);
    
    % Reading Surface
    Surftemp = Read_Surface([FreeSDir filesep Id filesep 'surf' filesep 'rh.white']);
    
    % Detecting 0s in label values
    ind = find(txt == 0);
    txt(ind) = 1;
    Surftemp.Is = txt;
    
    % Removing isolated clusters and refilling holes
    Surf = Surf_Corr(Surftemp);
    ind = find(Surf.Is == 1);
    txt(ind) = 0;
    Surf.Is(ind) = 0;
    
    % Saving the correct annot file
    OutAnnotFile = save_annotfiles(Surf.Is,AnnotFile,colortable);
    
    % Creating Lobar parcellation the correct annot file
    OutAnnotFile = [FreeSDir filesep Id filesep filesep 'label' filesep 'rh.aparc.lobes.annot'];
    OutAnnotFile = Aparc2Lobes(AnnotFile, 0, OutAnnotFile);
    OutAnnotFiles = strvcat(OutAnnotFiles, OutAnnotFile);
    %% =================== End of Left Hemisphere ======================= %
    
end
varargout{1} = OutAnnotFiles;
end