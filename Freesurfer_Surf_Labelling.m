function Surfa = Freesurfer_Surf_Labelling(IndAtlas,FreeDir, SubjID, OutFile);
%
% Syntax :
% Surfa = Freesurfer_Surf_Labelling(IndAtlas,FreeDir, SubjID, OutFile);
%
% This function says which points of the surface belongs to each anatomical
% structures taking into an account a given reference atlas.  The inter hemispheres
% points are labeled 0. It can also compute the mean cortical thickness
% using FreeSurfer Results.
%
% Input Parameters:
%   IndAtlas  : Individual Atlas.
%   FreeDir   : Freesurfer Directory
%   SubjID    : Subject ID
%   OutFile   : Output Thickness File
%
% Output Parameters:
%   Surfa     : Matlab variable containing surface information.
%
% Related references:
%
%
% See also: Red_Surf Smooth_Surf Plot_Surf Surf_Comp Plot_oversurf
% Exp_Surf
%
%__________________________________________________
% Authors: Yasser Aleman Gomez
% LIM
% December 17th 2012
% Version $1.0

IndAtlas = '/media/COSAS/Test/freesurfer/ch2/tmp/Atlased/c1sT1_atlas116.nii';
FreeDir = '/media/COSAS/Test/freesurfer';
SubjID = 'ch2';

opts.pipe.freesdir = FreeDir;
opts.pipe.subjId = SubjID;
if nargin < 4 
    atype = 'aal';
end
V = spm_vol(IndAtlas);
It = spm_read_vols(V);
switch atype
    case 'aal'
        CodeFile = which('atlas116.cod');
        fid = fopen(CodeFile);
        cont = 0;
        while 1
            cont = cont + 1;
            line = fgetl(fid);
            if ~ischar(line),   break,   end
            ind = strfind(line,'=');
            StructNames{1, cont} = line(ind+1:end);
            StructCodes(1,cont) = str2num(deblank(line(1:ind-1)));
        end
        fclose(fid);
        LIds = [1:2:90]; % Structures IDs, Left Hemisphere 
        ind = ismember(LIds,[37 38 41 42 71:78]);
        LIds(find(ind)) = [];
        RIds = [2:2:90]; % Structures IDs, Right Hemisphere
        ind = ismember(RIds,[37 38 41 42 71:78]);
        RIds(find(ind)) = [];
    otherwise
       StructCodes =  unique(It(:));
       StructCodes(StructCodes ==0) =[];
       Ns = length(StructCodes);
       for i =1:Ns
           StructNames{1, i} = sprintf('%.6d',StructCodes(i));
       end
end

%% ================== Mandatory input Files =============================%%
talfile = [opts.pipe.freesdir filesep opts.pipe.subjId filesep 'mri' filesep 'transforms' filesep 'talairach.lta'];
ctfiles = strvcat([opts.pipe.freesdir filesep opts.pipe.subjId  filesep 'surf' filesep 'lh.thickness'],[opts.pipe.freesdir filesep opts.pipe.subjId filesep 'surf' filesep 'rh.thickness']);
atfiles = strvcat([opts.pipe.freesdir filesep opts.pipe.subjId  filesep 'label' filesep 'lh.aparc.annot'],[opts.pipe.freesdir filesep opts.pipe.subjId filesep 'label' filesep 'rh.aparc.annot']);
surffiles = strvcat([opts.pipe.freesdir filesep opts.pipe.subjId  filesep 'surf' filesep 'lh.pial'],[opts.pipe.freesdir filesep opts.pipe.subjId filesep 'surf' filesep 'rh.pial']);
%% ================== End of Mandatory input Files ======================%%

%% ================== Reading Talairach Information =====================%%
cras1 = textread(talfile,'%s',5,'headerlines',20);
cras = char(cras1);
cras = [str2num(cras(3,:))  str2num(cras(4,:)) str2num(cras(5,:))];
%% ============== End of Reading Talairach Information ==================%%

%% =========%%%%%%%%%%%%%%% Left Hemisphere %%%%%%%%%%%%%%%% ============%%
% ----------- Reading Surface
[OutFiles, SurfF] = Exp_Surf(deblank(surffiles(1,:)), '0', '','', 'imp','n');
Surfl= SurfF{1};
Surfl.Name = 'LH.Pial';
Surfl.SurfData.vertices = Surfl.SurfData.vertices + repmat(cras,[size(Surfl.SurfData.vertices,1) 1]);
vertvoxl = (inv(V.mat)*[Surfl.SurfData.vertices ones(size(Surfl.SurfData.vertices),1)]')';

% ----------- Labelling

Itaux = It*0;
ind = ismember(It,[LIds]);
Itaux(ind) = It(ind);
Itaux = imdilate(imdilate(Itaux,strel([1 1 1; 1 1 1; 1 1 1])),strel([1 1 1; 1 1 1; 1 1 1]));
Is = spm_sample_vol(Itaux,double(vertvoxl(:,1)),double(vertvoxl(:,2)),double(vertvoxl(:,3)),0);

% ---------- Correcting Surface and Labelling
Surfl.Is = Is;
Is = Recur(Is, Surfl,3);
ind = find(Is == 0);
if ~isempty(ind)
    Is = Recur(Is, Surfl,2);
end
ind = find(Is == 0);
if ~isempty(ind)
    Is = Recur(Is, Surfl,1);
end
Surfl.Is = int16(Is);
[Surfl] = Surf_Corr(Surfl);

% --------- Removing Medial Wall
[txt,colortable] = read_cfiles(deblank(atfiles(1,:)));
names = char(colortable.struct_names);
ind = ismember(names(:,1:7),'unknown','rows');
MwId = colortable.table(find(ind),5);
ind = find(txt == MwId);
Surfl.Is(ind) = 0;
ind = find(txt == 0);
Surfl.Is(ind) = 0;

% --------- Computing Mean Cortical Thickness
[txt] = read_cfiles(deblank(ctfiles(1,:)));
sts = unique(Surfl.Is);
sts(sts==0) =[];
Ns = length(sts);
namesl = '';
for i = 1:Ns
    ind = find(Surfl.Is == sts(i));
    ctl(i) = mean(txt(ind));
    ind = find(StructCodes == sts(i));
    codesl(i) = StructCodes(ind);
    namesl = strvcat(namesl,char(StructNames{1,ind}));
end

%% =========%%%%%%%%%%% End of Left Hemisphere %%%%%%%%%%%%% ============%%



%% =========%%%%%%%%%%%%%%% Right Hemisphere %%%%%%%%%%%%%%%% ============%%
% ----------- Reading Surface
[OutFiles, SurfF] = Exp_Surf(deblank(surffiles(2,:)), '0', '','', 'imp','n');
Surfr= SurfF{1};
Surfr.Name = 'LH.Pial';
Surfr.SurfData.vertices = Surfr.SurfData.vertices + repmat(cras,[size(Surfr.SurfData.vertices,1) 1]);
vertvoxr = (inv(V.mat)*[Surfr.SurfData.vertices ones(size(Surfr.SurfData.vertices),1)]')';

% ----------- Labelling
Itaux = It*0;
ind = ismember(It,[RIds]);
Itaux(ind) = It(ind);
Itaux = imdilate(imdilate(Itaux,strel([1 1 1; 1 1 1; 1 1 1])),strel([1 1 1; 1 1 1; 1 1 1]));
Is = spm_sample_vol(Itaux,double(vertvoxr(:,1)),double(vertvoxr(:,2)),double(vertvoxr(:,3)),0);

% ---------- Correcting Surface and Labelling
Surfr.Is = Is;
Is = Recur(Is, Surfr,3);
ind = find(Is == 0);
if ~isempty(ind)
    Is = Recur(Is, Surfr,2);
end
ind = find(Is == 0);
if ~isempty(ind)
    Is = Recur(Is, Surfr,1);
end
Surfr.Is = int16(Is);
[Surfr] = Surf_Corr(Surfr);

% --------- Removing Medial Wall
[txt,colortable] = read_cfiles(deblank(atfiles(2,:)));
names = char(colortable.struct_names);
ind = ismember(names(:,1:7),'unknown','rows');
MwId = colortable.table(find(ind),5);
ind = find(txt == MwId);
Surfr.Is(ind) = 0;
ind = find(txt == 0);
Surfr.Is(ind) = 0;

% --------- Computing Mean Cortical Thickness
[txt] = read_cfiles(deblank(ctfiles(2,:)));
sts = unique(Surfr.Is);
sts(sts==0) =[];
Ns = length(sts);
namesr = '';
for i = 1:Ns
    ind = find(Surfr.Is == sts(i));
    ctr(i) = mean(txt(ind));
    ind = find(StructCodes == sts(i));
    codesr(i) = StructCodes(ind);
    namesr = strvcat(namesr,char(StructNames{1,ind}));
end
%% =========%%%%%%%%%%% End of Right Hemisphere %%%%%%%%%%%%% ============%%


%% ===================== Saving Results =================================%%
% -----------Reordering
codes = [codesl(:);codesr(:)];
names = strvcat(namesl,namesr);
ct = [ctl(:);ctr(:)];
[ncodes,order] = sort(codes);
nnames  = names(order,:);
nct = ct(order);

% --------- Saving 
if nargin < 4
    [pth, nm, ext] = fileparts(IndAtlas);
    OutFile = [pth filesep nm '_thickness.txt'];
end
fid = fopen(OutFile,'wt');
Ns = length(ncodes);
fprintf(fid,'%s             %s\n','Structure','C_Thickness(mm)');
for i = 1: Ns
    fprintf(fid,'%s     %s\n',nnames(i,:),num2str(nct(i)));
end
fclose(fid);

%% ================ End of Saving Results ===============================%%
Surfa{1,1} = Surfl;
Surfa{1,2} = Surfr;
return;
