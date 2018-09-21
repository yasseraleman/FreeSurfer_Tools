function hf = Plot_FreeS_results(SubjId,hemi,Sfile,CFile,cl,thmin,thmax);
% Syntax :
% Plot_FreeS_results(SubjId,hemi,Sfile,CFiles,cl);
%
% This function plots Freesurfer Statistics results over a surface
%
% Input Parameters:
%   SubjId      : Subject Id.
%   hemi        : Hemisphere (lh or rh)
%   Sfile       : Surface type
%   CFiles      : Characteristic file
%   cl          : Colormap
%   plotvar     : Plot variant. 1: No transparency. ~=1 Transparency for
%   different levels
%
% Output Parameters:
%    hf         : Figure handle.
%
% Related references:
%
%
% See also: Smooth_Surf Surf_Comp Red_Surf Plot_Surf Atlas_Surf Exp_Surf
%__________________________________________________
% Authors: Yasser Aleman Gomez
% Neuroimaging Department
% Cuban Neuroscience Center
% September 29th 2007
% Version $1.0

% SubjId = 'ch2';
% hemi = 'lh';
% Sfile = 'inflated';
% cl = 'hot';
 plotvar = 1;


[a,temp] = system('echo $SUBJECTS_DIR');
temp = temp';temp = deblank(temp(:)');
optst.pipe.freesdir = temp;

SurfFile = [ optst.pipe.freesdir filesep SubjId filesep 'surf' filesep hemi '.' Sfile];
CurvFile = [ optst.pipe.freesdir filesep SubjId filesep 'surf' filesep hemi '.curv'];
Annotfile = [ optst.pipe.freesdir filesep SubjId filesep 'label' filesep hemi '.aparc.annot'];

CFiles = [ optst.pipe.freesdir filesep SubjId filesep 'surf' filesep hemi CFile];
%CFiles = strvcat('/media/COSAS/Test/freesurfer/ch2/label/lh.aparc+Yeo2011_7Networks_N1000.annot');
%=====================Checking Input Parameters===========================%
% if nargin < 7
%     cl = 'jet';
% end
% 
% thmin = 0;
% thmax = 10;

[OutFiles, SurfF] = Exp_Surf(deblank(SurfFile), '0', '','', 'imp','n');
Surf = SurfF{1};
if ~exist('CFiles','var')|(isempty(CFiles))
    [CFiles,sts] = spm_select([1 5],'any','Selecting Thickness/Curvature Files','',cd);
end
tic;[Trip] = Vert_Neibp(double(Surf.SurfData.faces),size(Surf.SurfData.vertices,1),size(Surf.SurfData.faces,1));toc;
Temp = sum(Trip);
Trip(:,Temp==0) = [];
Trip = int32(Trip);
temp = Trip(:,3:end);
indz = find(temp == 0);
temp(indz) = 1;

Npoints = size(Surf.SurfData.vertices,1);
if isfield(Surf.SurfData,'FaceVertexCData')
    rmfield(Surf.SurfData,'FaceVertexCData');
end
%========================= Reading and Coloring Curvature ================%
[txtcurv] = read_cfiles(CurvFile);
at = txtcurv;
indm= find(txtcurv<0);
inda= find(txtcurv>=0);
at(indm)= 1;
at(inda)= 2;
Surf.Is = at;
col = [0.9 0.9 0.9;0.51 0.51 0.51];
Ncolor = size(col,1);
ut = sort(unique(Surf.Is));
re = floor(length(ut)/Ncolor); col = repmat(col,[re+1 1]);
for j = 1:size(ut,1)
    indpos = find(Surf.Is == ut(j)); Colors(indpos,:) = repmat(col(j,:),[size(indpos,1) 1]);
end
Surf.SurfData.FaceVertexCData = Colors;
Surfcurv = Surf;
%========================= Reading Characteristic  =======================%
[txtc,ctab] = read_cfiles(CFiles);
if  ~exist('thmin','var')
    thmin = min(txtc);
end
if  ~exist('thmax','var')
    thmax = max(txtc);
end
txtc((txtc<thmin)&(txtc>thmax)) = 0;
Surf.Is = txtc;
[Colors] = Surf_Color(Surf,cl);
ind = find(txtc~=0);
Surf.SurfData.FaceVertexCData(ind,:) = Colors(ind,:);
Surfchar = Surf;
%========================= Reading and Annot File  =======================%
[txt,ctab] = read_cfiles(Annotfile);
Names = char(ctab.struct_names);
indu = find(ismember(Names(:,1:7),'unknown','rows') == 1);
indcc = find(ismember(Names(:,1:14),'corpuscallosum','rows') == 1);
Names([indu indcc],:) = [];
rids = ctab.table([indu indcc],5);
ctab.table([indu indcc],:) = [];
indcc = find(ismember(txt,rids));
txt(indcc) = 0;

sts = unique(txt);
sts(sts ==0) =[];
Nst = length(sts);
ColorsV = Surf.SurfData.FaceVertexCData*0;
for i = 1:Nst
    ind = find(txt==sts(i));
    indc = find(ctab.table(:,5)==sts(i));
    if isempty(indc)
        Matname = 'unknown_structure';
        Color =   [1 1 1];      % Color
    else
        Matname = Names(indc,:);
        Color =  ctab.table(indc,1:3)/255;       % Color
    end
    ColorsV(ind,:) = repmat(Color,[length(ind) 1 ]);
end
temp1 = txt(temp);
temp1(indz) =  max(temp1(:))+1;
NewMat = temp1-repmat(min(temp1')',[1 size(temp1,2)]);
NewMat(indz) = 0;
a = logical(sum(logical(NewMat)')');
indc = find(a);
Surf.SurfData.FaceVertexCData(indc,:) = ColorsV(indc,:)*0;
indfa = find(sum(ismember(Surf.SurfData.faces,indc)')<2);
Surfbound = Surf;
Surfbound.SurfData.faces(indfa,:) = [];
%=========================================================================%

% --------------- Plotting Surface  --------------------------------------%
colordef black;
hf = figure('numbertitle','off','Color','white','Position',[0 0 1200 900]);
if plotvar == 1
    custom_plotsurf(Surf);
else
    custom_plotsurf(Surfcurv);
    hold on
    Surfchar.SurfData.vertices = Surfchar.SurfData.vertices*1.1;
    strsurf=patch(Surfchar.SurfData,'edgecolor','none', 'tag','patch','facelighting','gouraud');
    set(strsurf,'FaceAlpha',.7);
    Surfbound.SurfData.vertices = Surfbound.SurfData.vertices*1.2;
    strsurf=patch(Surfbound.SurfData,'edgecolor','none', 'tag','patch','facelighting','gouraud');
end
if strcmp(hemi,'lh')
    view([270 0]);axis off;axis tight;axis equal; %h=title(['Left Hemisphere. Lateral View']);set(h,'FontSize',15,'FontName','Arial');
    camlight;
elseif strcmp(hemi,'rh')
    view([90 0]);axis off;axis tight;axis equal;% h=title(['Right Hemisphere. Lateral View']);set(h,'FontSize',15,'FontName','Arial');
    camlight;
end

% ----------- Creating Legend ---------------------------------------%
if sum(txtc-floor(txtc)) ~=0
    h = colorbar;
    if ~strcmp(cl,'spectral')
        colormap(cl);
    else
        colormap(spectral(64));
    end
    range = thmax-thmin;values = thmin:range/10:thmax;
    set(h,'YTickLabel',num2str(values'));
end
%========================End of main program==============================%
return;


function custom_plotsurf(Surf);
Surf.SurfData.FaceColor = 'interp';
if isunix
    strsurf=patch(Surf.SurfData,'edgecolor','none', 'tag','patch','facelighting','gouraud');
else
    strsurf=patch(Surf.SurfData,'edgecolor','black','tag', 'patch','facelighting','gouraud');
end
return