function varargout = Create_VolumetricAtlas_from_Annot(varargin);

fsDirectory = '/media/Data/PROCESSING_RESULTS/HCP/5-freesurfer_processing/';
Id = 'fsaverage';
outAtlas = [fsDirectory filesep Id filesep 'mri' filesep 'aparc_scale05_Nzones0410+aseg.nii'];

talairachFile = [fsDirectory filesep Id filesep 'mri' filesep 'transforms' filesep 'talairach.lta'];

lhwhite =  [fsDirectory filesep Id filesep 'surf' filesep 'lh.white' ];
lhpial  =  [fsDirectory filesep Id filesep 'surf' filesep 'lh.pial' ];
lhannot =  [fsDirectory filesep Id filesep 'label' filesep 'lh.aparc_scale05_Nzones0410.annot' ];

rhwhite =  [fsDirectory filesep Id filesep 'surf' filesep 'rh.white' ];
rhpial  =  [fsDirectory filesep Id filesep 'surf' filesep 'rh.pial' ];
rhannot =  [fsDirectory filesep Id filesep 'label' filesep 'rh.aparc_scale05_Nzones0410.annot' ];

% Converting Ribbon
cad = ['mri_convert -i ' fsDirectory filesep Id filesep 'mri' filesep 'ribbon.mgz -o ' fsDirectory filesep Id filesep 'tmp' filesep 'ribbon.nii'];  
% system(cad);
% outputs = flip_to_axial_nii([fsDirectory filesep Id filesep 'tmp' filesep 'ribbon.nii'],0);
% delete([fsDirectory filesep Id filesep 'tmp' filesep 'ribbon.nii']);
% movefile(outputs,[fsDirectory filesep Id filesep 'tmp' filesep 'ribbon.nii']);
    
Vr = spm_vol([fsDirectory filesep Id filesep 'tmp' filesep 'ribbon.nii']);

Ir =spm_read_vols(Vr);
Itrtemp = Ir;
ind = find(Ir == 41);
Ir(ind) =0;
ind = find(Ir == 2);
Ir(ind) =0;
ind = find(Ir>0);
Ir(ind) =1;

% Converting Aparc+Aseg
cad = ['mri_convert -i ' fsDirectory filesep Id filesep 'mri' filesep 'aparc+aseg.mgz -o ' fsDirectory filesep Id filesep 'tmp' filesep 'aparc+aseg.nii'];  
system(cad);
outputs = flip_to_axial_nii([fsDirectory filesep Id filesep 'tmp' filesep 'aparc+aseg.nii'],0);
delete([fsDirectory filesep Id filesep 'tmp' filesep 'aparc+aseg.nii']);
movefile(outputs,[fsDirectory filesep Id filesep 'tmp' filesep 'aparc+aseg.nii']);

Vaparc = spm_vol([fsDirectory filesep Id filesep 'tmp' filesep 'aparc+aseg.nii']);

%%   Creciendo el RIBBON  Left Hemisphere
[IAL,stCodes,stNames] = refill_thicknessSpace(lhwhite,lhpial, lhannot, talairachFile, Vr);
lhNames = [repmat('ctx-lh-',[size(stNames,1) 1]) char(stNames)];
lhCodes = stCodes;


[IAR,stCodes,stNames] = refill_thicknessSpace(rhwhite,rhpial, rhannot, talairachFile, Vr);
rhNames = [repmat('ctx-rh-',[size(stNames,1) 1]) char(stNames)];
rhCodes = stCodes + 1000;


IAR(IAR~=0) = IAR(IAR~=0)+1000;
IA = zeros(Vr.dim(1),Vr.dim(2),Vr.dim(3));
 
ind = find((IAL~=0)&(IAR~=0));
if unique(IAL(ind)) == unique(IAR(ind))-1000
    IA = IAR + IAL;
    IA(ind) = 1000;
end

%% Correcting Ribbon

[IA] = Atlas_Corr(Ir,Ir.*IA);


% Subcortical Structures
SubCL = [10:13 17:18 26];
SubCR = [ 49:54 58];

subcLHNames = strvcat('Left-Thalamus-Proper','Left-Caudate',...
'Left-Putamen','Left-Pallidum','Left-Hippocampus','Left-Amygdala',...                  
'Left-Accumbens-area');

subcRHNames = strvcat('Right-Thalamus-Proper','Right-Caudate',...
'Right-Putamen','Right-Pallidum','Right-Hippocampus','Right-Amygdala',...                  
'Right-Accumbens-area');

Iaparc = spm_read_vols(Vaparc);

Iaparc(ismember(Iaparc,[SubCL SubCR]) == 0) = 0;

ind = find(Iaparc);
IA(ind) = Iaparc(ind);
IA(IA == 1000) = 0;
IA(IA == 2000) = 0;
Vout = Vaparc;
Vout.fname = outAtlas;
spm_write_vol(Vout,IA);

return;



 
 
function [IA,stCodes,stNames] = refill_thicknessSpace(whiteSurf,pialSurf, annotFile, talairachFile,V);
Surfw = Read_Surface(whiteSurf);
Surfp = Read_Surface(pialSurf);

[txt, ctab, colors] = read_cfiles(annotFile);

N = 100;
if exist(talairachFile,'file');
    cras1 = textread(talairachFile,'%s',5,'headerlines',20);
    cras = char(cras1);cras = [str2num(cras(3,:))  str2num(cras(4,:)) str2num(cras(5,:))];
else
    cras = [0 0 0];
end
IA = zeros(V.dim(1),V.dim(2),V.dim(3));
for i = 1:size(ctab.table,1)
    ind = find(txt == ctab.table(i,5));
    xw = Surfw.SurfData.vertices(ind,1); yw = Surfw.SurfData.vertices(ind,2); zw = Surfw.SurfData.vertices(ind,3);
    xp = Surfp.SurfData.vertices(ind,1); yp = Surfp.SurfData.vertices(ind,2); zp = Surfp.SurfData.vertices(ind,3);
    norms = sqrt((xp-xw).^2+(yp-yw).^2+(zp-zw).^2);ind = find(norms==0);xp(ind) = [];xw(ind) = [];yw(ind) = [];yp(ind) = [];zp(ind) = [];zw(ind) = [];norms(ind) = [];
    u = (xp-xw)./norms; v = (yp-yw)./norms; w = (zp-zw)./norms;
    Cp = [xp yp zp];
    Cw = [xw yw zw];
    dista = (sqrt((Cp(:,1)-Cw(:,1)).^2+(Cp(:,2)-Cw(:,2)).^2+(Cp(:,3)-Cw(:,3)).^2));
    steps = [0:1/N:1];dista = repmat(steps,[size(dista,1) 1]).*repmat(dista,[1 N+1]);
    CintX = [repmat(Cw(:,1),[1 size(dista,2)])+repmat(u,[1 size(dista,2)]).*dista ]';%siz = size(CintX);%CintX = CintX(:);
    CintY = [repmat(Cw(:,2),[1 size(dista,2)])+repmat(v,[1 size(dista,2)]).*dista ]';%CintY = CintY(:);
    CintZ = [repmat(Cw(:,3),[1 size(dista,2)])+repmat(w,[1 size(dista,2)]).*dista ]';%CintZ = CintZ(:);
    voxs = (inv(V.mat)*[CintX(:)+cras(1) CintY(:)+cras(2) CintZ(:)+cras(3) ones(size(CintX(:),1),1)]')';
    %voxX = reshape(voxs(:,1),[siz]);voxY = reshape(voxs(:,2),[siz]);voxZ = reshape(voxs(:,3),siz);
    a = unique(round(voxs),'rows');ind = sub2ind(V.dim(1:3),a(:,1),a(:,2),a(:,3));
    IA(ind) = 1000+(i-1);
end
stCodes = [1000:1000+size(ctab.table,1)-1]';
stNames = ctab.struct_names;
return;