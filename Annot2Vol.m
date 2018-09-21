function OutVol = Annot2Vol(SubjId,OutVol);
%
% Syntax :
% OutAnnot = Masking_Annot(AnnotFile,CurvFile,thr,OutAnnot);
%
% This script masks annotation files using the curvature information
% 
%
% Input Parameters:
%   AnnotFile            : Annotation File
%   CurvFile             : Curvature file
%   thr                  : Curvature threshold for binarization
%   OutAnnot             : New Annotation filename
%
% Output Parameters:
%   OutAnnot             : New Annotation filename
%
% See also: 
%__________________________________________________
% Authors: Yasser Aleman Gomez 
% LIM, HUGGM
% March 22th 2012
% Version $1.0

%=====================Checking Input Parameters===========================%

%=========================================================================%
%=====================     Main Program    ===============================%
SubjId = 'ch2';
Annottype = 'aparc.a2009s';
distance = 2; %mm
names = strvcat('Medial_wall','corpuscallosum','Unknown');
% [a,temp] = system('echo $SUBJECTS_DIR');
% freesurferdir = [deblank(temp)']';
freesurferdir = '/media/COSAS/Test/freesurfer'
maindir =[freesurferdir filesep SubjId];
%% Left Hemisphere
[vertices, faces] = freesurfer_read_surf([maindir filesep 'surf' filesep 'lh.white' ]);Surfw.SurfData.vertices = vertices;Surfw.SurfData.faces = faces;Surfw.Name = 'White';
[vertices, faces] = freesurfer_read_surf([maindir filesep 'surf' filesep 'lh.pial' ]);Surfp.SurfData.vertices = vertices;Surfp.SurfData.faces = faces;Surfp.Name = 'Pial';
[label1, txt2, colortable] = read_annotation([maindir filesep 'label' filesep 'lh.' Annottype '.annot' ]); %
normals = Surfp.SurfData.vertices-Surfw.SurfData.vertices;
norma = sqrt(sum((normals').^2))+eps;
norms = normals./repmat(norma',[1 3]);
namestemp = char(colortable.struct_names);
stid = 0;
for i = 1:size(names,1)
    ind = ismember(namestemp(:,1:length(deblank(names(i,:)))),deblank(names(i,:)),'rows');
    stid = [stid colortable.table(find(ind),5);]
end
stid(1) = [];
indexs = ismember(txt2,stid(:));
xw = Surfw.SurfData.vertices(indexs,1); yw = Surfw.SurfData.vertices(indexs,2); zw = Surfw.SurfData.vertices(indexs,3);distw=sqrt(xw.^2+yw.^2+zw.^2);norms = norms(indexs,:);
Cp = [xw xw xw]+norms*distance;
Cw = [xw yw zw]-norms*distance;
N=5;
dista = (sqrt((Cp(:,1)-Cw(:,1)).^2+(Cp(:,2)-Cw(:,2)).^2+(Cp(:,3)-Cw(:,3)).^2));
steps = [0:1/N:1];dista = repmat(steps,[size(dista,1) 1]).*repmat(dista,[1 N+1]);
CintX = [repmat(Cw(:,1),[1 size(dista,2)])+repmat(norms(:,1),[1 size(dista,2)]).*dista ]';siz = size(CintX);%CintX = CintX(:);
CintY = [repmat(Cw(:,2),[1 size(dista,2)])+repmat(norms(:,2),[1 size(dista,2)]).*dista ]';%CintY = CintY(:);
CintZ = [repmat(Cw(:,3),[1 size(dista,2)])+repmat(norms(:,3),[1 size(dista,2)]).*dista ]';%CintZ = CintZ(:);
%% Right Hemisphere

%% Creating Volume
%---------------------------------

%%%%%%%%% Leyendo el centro de coordenadas RAS
if isunix
    cras = textread([dira 'ch2/mri/transforms/talairach.lta'],'%s',5,'headerlines',20);
else
    cras = textread([dira 'ch2\mri\transforms\talairach.lta'],'%s',5,'headerlines',20);
end
cras = char(cras);cras = str2num(cras(3:end,:))';

V = spm_vol([dira 'ch2.nii']); It =spm_read_vols(V);       I =0*It;

voxs = (inv(V.mat)*[CintX(:)+cras(1) CintY(:)+cras(2) CintZ(:)+cras(3) ones(size(CintX(:),1),1)]')';voxX = reshape(voxs(:,1),[siz]);voxY = reshape(voxs(:,2),[siz]);voxZ = reshape(voxs(:,3),[siz]);