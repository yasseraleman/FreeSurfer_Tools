function varargout = Aparc2Lobes(varargin);
%
% Syntax :
%     OutAnnotFile = Aparc2Lobes(AnnotFile, boolcingulate, OutAnnotFile);
%
% This function does an approximate mapping of individual 'Desikan-Killiany' ROIs 
% (found in ?h.aparc.annot) to the lobes
%
% Input Parameters:
%        AnnotFile              : Annotation File
%        boolcingulate          : Boolean variable to include or not Cingulate Lobe 
%                                 (0 = do not include, 1 = include)
%        OutAnnotFile           : Saving the lobar parcellation in an
%                                 annotation file
%
% Output Parameters:
%        OutAnnotFile           : Saving the lobar parcellation in an
%                                 annotation file. If nargin <3
%                                 OutAnnotFile is a vector file containing
%                                 lobar parcellation.
%
% See also: save_annotfiles
%__________________________________________________
% Authors: Yasser Aleman Gomez 
% LIM, HUGGM
% November 13th 2014
% Version $1.0

%% ============================= Checking Inputs ======================= %%
AnnotFile = varargin{1};
if ~exist(AnnotFile,'file')
    errordlg('The Annotation File Does not exist');
    return;
else
    try
        % Reading Parcellation
        [txt,colortable] = read_cfiles(AnnotFile);
        if ~isfield(colortable,'table')
            errordlg('The File is not an annotation file');
            return;
        end
    catch
        errordlg('Error reading Annotation File');
        return;
    end
end

if nargin == 1
    boolcingulate = 0;
end
if nargin == 2
    boolcingulate = logical(varargin{2});
end
if nargin == 3
    OutAnnotFile = varargin{3};
    boolcingulate = logical(varargin{2});
end
if nargin > 3
    errordlg('To Many Input Parameters');
    return;
end
if nargout > 1
    errordlg('To Many Output Parameters');
    return;
end
%% ========================= End of Checking Inputs ==================== %%

%% ============================= Defining Lobes ======================== %%
if boolcingulate
    % Frontal Lobe
    Fron = strvcat('caudalmiddlefrontal',...
        'lateralorbitofrontal',...
        'medialorbitofrontal',...
        'paracentral',...
        'parsopercularis',...
        'parsorbitalis',...
        'parstriangularis',...
        'precentral',...
        'rostralmiddlefrontal',...
        'superiorfrontal',...
        'frontalpole');
    
    % Parietal Lobe
    Par = strvcat('inferiorparietal',...
        'isthmuscingulate',...
        'postcentral',...
        'precuneus',...
        'superiorparietal',...
        'supramarginal');
    
    % Cingulate Lobe
    Cing = strvcat('rostralanteriorcingulate',...
                   'caudalanteriorcingulate',...
                   'posteriorcingulate',...
                   'isthmuscingulate');
    
else
    % Frontal Lobe
    Fron = strvcat('caudalanteriorcingulate',...
        'caudalmiddlefrontal',...
        'lateralorbitofrontal',...
        'medialorbitofrontal',...
        'paracentral',...
        'parsopercularis',...
        'parsorbitalis',...
        'parstriangularis',...
        'precentral',...
        'rostralanteriorcingulate',...
        'rostralmiddlefrontal',...
        'superiorfrontal',...
        'frontalpole');
    
    % Parietal Lobe
    Par = strvcat('inferiorparietal',...
        'isthmuscingulate',...
        'postcentral',...
        'posteriorcingulate',...
        'precuneus',...
        'superiorparietal',...
        'supramarginal');
end

% Temporal Lobe
Temp = strvcat('bankssts',...
    'entorhinal',...
    'fusiform',...
    'inferiortemporal',...
    'middletemporal',...
    'parahippocampal',...
    'superiortemporal',...
    'temporalpole',...
    'transversetemporal');

% Occipital Lobe
Occ = strvcat('cuneus',...
    'lateraloccipital',...
    'lingual',...
    'pericalcarine');

% Insula Lobe
Ins = 'insula';
%% =======================End of Defining Lobes ======================== %%


%% ================= Removing white spaces from Structure Names ======== %%
totnames1 = char(colortable.struct_names);
totnames = '';
for i = 1:size(totnames1,1)
    totnames = strvcat(totnames,deblank(totnames1(i,:)));
end
ind = find(sum(isspace(totnames))==size(colortable.struct_names,1));
totnames(:,ind) = [];
%% ========== End of Removing white spaces from Structure Names ======== %%

%% ================= Removing white spaces from Frontal Lobe =========== %%
ind = find(sum(isspace(Fron))==size(Fron,1));
Fron(:,ind) = [];
ind = ismember(totnames(:,1:size(Fron,2)),Fron,'rows');
FronIds = colortable.table(ind,5);
indfron = ismember(txt,FronIds);
%% ========== End of Removing white spaces from Frontal Lobe =========== %%

%% ================= Removing white spaces from Parietal Lobe ========== %%
ind = find(sum(isspace(Par))==size(Par,1));
Par(:,ind) = [];
ind = ismember(totnames(:,1:size(Par,2)),Par,'rows');
ParIds = colortable.table(ind,5);
indpar = ismember(txt,ParIds);
%% ========== End of Removing white spaces from Parietal Lobe ========== %%

%% ================= Removing white spaces from Temporal Lobe ========== %%
ind = find(sum(isspace(Temp))==size(Temp,1));
Temp(:,ind) = [];
ind = ismember(totnames(:,1:size(Temp,2)),Temp,'rows');
TempIds = colortable.table(ind,5);
indtemp = ismember(txt,TempIds);
%% ========== End of Removing white spaces from Temporal Lobe ========== %%

%% ================= Removing white spaces from Occipital Lobe ========= %%
ind = find(sum(isspace(Occ))==size(Occ,1));
Occ(:,ind) = [];
ind = ismember(totnames(:,1:size(Occ,2)),Occ,'rows');
OccIds = colortable.table(ind,5);
indocc = ismember(txt,OccIds);
%% ========== End of Removing white spaces from Occipital Lobe ========= %%

if boolcingulate
    %% ================= Removing white spaces from Occipital Lobe ========= %%
    ind = find(sum(isspace(Cing))==size(Cing,1));
    Cing(:,ind) = [];
    ind = ismember(totnames(:,1:size(Cing,2)),Cing,'rows');
    CingIds = colortable.table(ind,5);
    indcing = ismember(txt,CingIds);
    %% ========== End of Removing white spaces from Occipital Lobe ========= %%
end

%% ============= Creating New Label Id vector and ColorTable =========== %%
lobparc = txt*0;

lobparc(indfron) = 1;  % Frontal Lobe
stnames{2,1} = 'frontallobe ';
lobparc(indpar) = 2;   % Parietal Lobe
stnames{3,1} = 'parietallobe ';
lobparc(indtemp) = 3;  % Temporal Lobe
stnames{4,1} = 'temporallobe ';
lobparc(indocc) = 4;   % Occipital Lobe
if boolcingulate
    lobparc(indcing) = 5;   % Cingulate Lobe
    stnames{5,1} = 'linguallobe ';
end

if find(unique(lobparc) == 0)
    stnames{1,1} = 'unknown';
    stnames{2,1} = 'frontallobe';
    lobparc(indpar) = 2;   % Parietal Lobe
    stnames{3,1} = 'parietallobe';
    lobparc(indtemp) = 3;  % Temporal Lobe
    stnames{4,1} = 'temporallobe';
    lobparc(indocc) = 4;   % Occipital Lobe
    stnames{5,1} = 'occipitallobe';
    if boolcingulate
        lobparc(indcing) = 5;   % Cingulate Lobe
        stnames{6,1} = 'linguallobe';
    end
else
    stnames{1,1} = 'frontallobe';
    lobparc(indpar) = 2;   % Parietal Lobe
    stnames{2,1} = 'parietallobe';
    lobparc(indtemp) = 3;  % Temporal Lobe
    stnames{3,1} = 'temporallobe';
    lobparc(indocc) = 4;   % Occipital Lobe
    stnames{4,1} = 'occipitallobe';
    if boolcingulate
        lobparc(indcing) = 5;   % Cingulate Lobe
        stnames{5,1} = 'linguallobe';
    end
end
    
sts = unique(lobparc);
if find(sts ==0)
    colors = [255 255 255; 12 73 153; 7 96 32; 210 113 26; 126 18 33; 187 205 31];
else
    colors = [12 73 153; 7 96 32; 210 113 26; 126 18 33; 187 205 31];
end
ctab = [colors colors(:,1)*0 colors(:,1)+colors(:,2)*2^8+colors(:,3)*2^16];
Nstruct = length(sts);
for i = 1:Nstruct
    ind = find(lobparc==sts(i));
    if sts(i) == 0
        ctab(i,:)  = [255 255 255 0 16777215];
        lobparc(ind) = 16777215;
    else
        lobparc(ind) = ones(size(ind,1),1)*ctab(i,5);
    end
end
colortable.numEntries = length(sts);
colortable.orig_tab = 'Custom Colortable';
colortable.struct_names = stnames;
colortable.table = ctab;

%% ====================== End of Defining Lobes ======================== %%

if nargin <= 2;
    [a, lobparc] = ismember(lobparc,colortable.table(:,5));
    varargout{1} = lobparc;
else
    OutAnnotFile = save_annotfiles(lobparc,OutAnnotFile,colortable);
    varargout{1} = OutAnnotFile;
end
return;