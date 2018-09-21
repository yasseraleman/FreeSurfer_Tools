function Multi_Create_Lobar_ConvexHull(FreesDir,Ids);
opts.pipe.freesdir = FreesDir;

setenv('SUBJECTS_DIR',opts.pipe.freesdir);
%Ids = textread(Ids,'%s');
Ns = size(Ids,1);
for i = 1:Ns
    %Cadini = 'Subject ID      CHull_L_Frontal     CHull_L_Parietal     CHull_L_Temporal     CHull_L_Occipital     CHull_R_Frontal     CHull_R_Parietal     CHull_R_Temporal     CHull_R_Occipital';
    Cadini = 'Subject ID      CHull_Left     CHull_Right';
    disp(['Processing Subject ' num2str(i) ' of ' num2str(Ns)]);
    opts.pipe.subjId = char(deblank(Ids(i,:)));
    if ~exist([opts.pipe.freesdir filesep opts.pipe.subjId filesep 'tmp' filesep 'aparc+aseg+3mm.nii'],'file')
        cad = ['mri_aparc2aseg --s ' opts.pipe.subjId ' --annot aparc --wmparc-dmax 3 --labelwm --hypo-as-wm --o ' [opts.pipe.freesdir filesep opts.pipe.subjId filesep 'tmp' filesep 'aparc+aseg+3mm.nii']];
        system(cad);
        [ufilename] = Unified_ctx([opts.pipe.freesdir filesep opts.pipe.subjId filesep 'tmp' filesep 'aparc+aseg+3mm.nii'],'aparc+aseg',1);
    else
        [ufilename] = Unified_ctx([opts.pipe.freesdir filesep opts.pipe.subjId filesep 'tmp' filesep 'aparc+aseg+3mm.nii'],'aparc+aseg',1);
    end
    [Surfl,Surfr] = Create_Lobar_ConvexHull(ufilename);
    
    hdf = figure('Color',[0 0 0]);
    subplot(1,2,1);
    custom_plotsurf(Surfl); h=title(['Left Hemisphere: Lateral View']);set(h,'Color',[1 1 1]);
    view([270 0]); 
    axis off;
    axis tight;axis equal;
    camlight;
    subplot(1,2,2);
    custom_plotsurf(Surfr);h=title(['Right Hemisphere: Lateral View']);set(h,'Color',[1 1 1]);
    view([90 0]);
    axis off;
    axis tight;axis equal; 
    camlight;
%     imname =  [opts.pipe.freesdir filesep opts.pipe.subjId '-CHull_Lobar_Parc.jpg'];
%     F = getframe(gcf);
%     imwrite(F.cdata,imname,'jpg');
%     close(hdf);
%%  ====================== Just for Lobes ============================== %%
    ind = find(Surfl.StructS(:,1) == 1); Fl = Surfl.StructS(ind,2);
    ind = find(Surfl.StructS(:,1) == 2); Pl = Surfl.StructS(ind,2);
    ind = find(Surfl.StructS(:,1) == 3); Tl = Surfl.StructS(ind,2);
    ind = find(Surfl.StructS(:,1) == 4); Ol = Surfl.StructS(ind,2);
    
    ind = find(Surfr.StructS(:,1) == 1); Fr = Surfr.StructS(ind,2);
    ind = find(Surfr.StructS(:,1) == 2); Pr = Surfr.StructS(ind,2);
    ind = find(Surfr.StructS(:,1) == 3); Tr = Surfr.StructS(ind,2);
    ind = find(Surfr.StructS(:,1) == 4); Or = Surfr.StructS(ind,2);
    cadt = [char(deblank(Ids(i,:))) '     ' num2str(Fl) '            ' num2str(Pl) '              ' num2str(Tl) '              ' num2str(Ol) '              ' num2str(Fr) '              ' num2str(Pr) '              ' num2str(Tr) '             ' num2str(Or)];
    Cadini = strvcat(Cadini,cadt);
    Outfile = [opts.pipe.freesdir filesep opts.pipe.subjId 'Results_CHull_Lobar.txt'];
    fid = fopen(Outfile,'wt');
    for j = 1:size(Cadini,1)
        fprintf(fid,'%s\n',Cadini(j,:));
    end
    fclose(fid);
%%  ================= End of Just for Lobes ============================ %%
%%  ================ Just for Hemispheres ============================== %%
%      ind = find(Surfl.StructS(:,1) == 1); Fl = Surfl.StructS(ind,2);
%      ind = find(Surfl.StructS(:,1) == 2); Pl = Surfl.StructS(ind,2);
%      ind = find(Surfl.StructS(:,1) == 3); Tl = Surfl.StructS(ind,2);
%      ind = find(Surfl.StructS(:,1) == 4); Ol = Surfl.StructS(ind,2);
%      ind = find(Surfl.StructS(:,1) == 5); Cl = Surfl.StructS(ind,2);
%      
%      ind = find(Surfr.StructS(:,1) == 1); Fr = Surfr.StructS(ind,2);
%      ind = find(Surfr.StructS(:,1) == 2); Pr = Surfr.StructS(ind,2);
%      ind = find(Surfr.StructS(:,1) == 3); Tr = Surfr.StructS(ind,2);
%      ind = find(Surfr.StructS(:,1) == 4); Or = Surfr.StructS(ind,2);
%      ind = find(Surfr.StructS(:,1) == 5); Cr = Surfr.StructS(ind,2);
%      
%      cadt = [char(deblank(Ids(i,:))) '     ' num2str(Fl+Pl+Tl+Ol+Cl) '            ' num2str(Fr+Pr+Tr+Or+Cr)];
%      Cadini = strvcat(Cadini,cadt);
%      Outfile = [opts.pipe.freesdir filesep opts.pipe.subjId '-Results_CHull.txt'];
%      fid = fopen(Outfile,'wt');
%      for j = 1:size(Cadini,1)
%          fprintf(fid,'%s\n',Cadini(j,:));
%      end
%      fclose(fid);
     %%  ================ End of Just for Hemispheres ================== %%
end

return

function custom_plotsurf(Surf,cl);
if nargin ==1
    cl = 'jet';
end
[Colors] = Surf_Color(Surf,cl);
Surf.SurfData.FaceVertexCData = Colors;
Surf.SurfData.FaceColor = 'interp';
strsurf=patch(Surf.SurfData,'edgecolor','none', 'tag','patch','facelighting','gouraud');
return
function [Colors] = Surf_Color(Surf,cl);
% Syntax :
% [Colors] = Surf_Color(Surf,cl);
%
% This function plots the surfaces contained in the cell file Surfa. 
%
% Input Parameters:
%   Surf        : Surface variable.
%   cl          : Colormap used to see the results.
%
% Output Parameters:
%  Colors       : Output colormap matrix.
%
% Related references:
% 
%
% See also: Smooth_Surf Surf_Comp Plot_Surf Plot_oversurf Exp_Surf
%__________________________________________________
% Authors: Yasser Aleman Gomez
% Neuroimaging Department
% Cuban Neuroscience Center
% September 29th 2007
% Version $1.0

%=========================Main program====================================%
txt=Surf.Is;
if sum(txt-floor(txt))
    [ut,i1t,it2] = unique(txt);
    switch cl
        case 'hot'
            col = hot(size(ut,1));
        case 'hsv'
            col = hsv(size(ut,1));
        case 'jet'
            col = jet(size(ut,1));
        case 'cool'
            col = cool(size(ut,1));
        case 'bone'
            col = bone(size(ut,1));
        case 'pink'
            col = pink(size(ut,1));
        case 'winter'
            col = winter(size(ut,1));
        case 'autumn'
            col = autumn(size(ut,1));
        case 'summer'
            col = summer(size(ut,1));
        case 'spring'
            col = spring(size(ut,1));
        case 'copper'
            col = copper(size(ut,1));
        case 'spectral'
            col = spectral(size(ut,1));
    end
    Colors = col(it2,:);
elseif sum(txt-floor(txt)) ==0
    %col = [213 221 227;175 206 227;149 196 228; 120 186 230;87 172 231;24 146 232;6 73 119;244 207 154;244 192 117;244 179 86;244 161 43; 212 133 20;158 101 19; 113 71 12]/255;
    col = [1 0 0;0 1 0; 0 0 1; 1 1 0;0 1 1;1 0 1;1 0.5 0;0 0.5 1;1 0 0.5;0.5 1 0; 0.25 0.5 1;0.6 0.3 0.52;0.7 0.5 0.9];
    Ncolor = size(col,1);
    ut = sort(unique(txt));
    re = floor(length(ut)/Ncolor); col = repmat(col,[re+1 1]);
    for j = 1:size(ut,1)
        indpos = find(txt == ut(j)); Colors(indpos,:) = repmat(col(j,:),[size(indpos,1) 1]);
    end
    
    %Creating boundary lines
%     Is=txt;Is(Is==0)=max(Is)+1;
%     St=unique(Is);
%     Ns =size(St,1);
%     neibb = size(Surf.SurfData.vertices,1)+1;
%     Ist = 0*Is;
%     for i =1:Ns
%         ind = find(Is == St(i));
%         for j=1:size(ind,1)
%             Neigh=Surf.SurfData.faces(nonzeros(Surf.Tri(ind(j),3:end)),:);Neigh=unique(Neigh(:));Neigh(Neigh==ind(j))=[];
%             c=accumarray(double(Is(Neigh)),ones(size(Is(Neigh),1),1));
%             if size(nonzeros(c),1)~=1&(~logical(sum(ismember(Neigh,neibb))))
%                 Ist(ind(j))=Is(ind(j));
%             end
%         end
%         neibb = unique([neibb;find(Ist~=0)]);
%     end
%     indl =find(Ist~=0);
%     Colors(indl,:) =repmat([0    0    0],[size(indl,1) 1]);
end
%========================End of main program==============================%
return

function s = spectral(m)
%SPECTRAL Black-purple-blue-green-yellow-red-white color map.
%
%         map = spectral(num_colors)
%
% SPECTRAL(M) returns an M-by-3 matrix containing a "spectral" colormap.
% SPECTRAL, by itself, is the same length as the current colormap.
%
% For example, to reset the colormap of the current figure:
%
%           colormap(spectral)
%
% See also HSV, GRAY, PINK, HOT, COOL, BONE, COPPER, FLAG,
%          COLORMAP, RGBPLOT.

if nargin < 1, m = size(get(gcf,'colormap'),1); end
base = [
  0.0000 0.0000 0.0000
  0.4667 0.0000 0.5333
  0.5333 0.0000 0.6000
  0.0000 0.0000 0.6667
  0.0000 0.0000 0.8667
  0.0000 0.4667 0.8667
  0.0000 0.6000 0.8667
  0.0000 0.6667 0.6667
  0.0000 0.6667 0.5333
  0.0000 0.6000 0.0000
  0.0000 0.7333 0.0000
  0.0000 0.8667 0.0000
  0.0000 1.0000 0.0000
  0.7333 1.0000 0.0000
  0.9333 0.9333 0.0000
  1.0000 0.8000 0.0000
  1.0000 0.6000 0.0000
  1.0000 0.0000 0.0000
  0.8667 0.0000 0.0000
  0.8000 0.0000 0.0000
  0.8000 0.8000 0.8000
];
n = length(base);
X0 = linspace (1, n, m);
s = interp1(1:n,base,X0);

return

function [outfilename] = Unified_ctx(infilename,afile,delopt);
%
% Syntax :
% [outfilename] = Unified_ctx(infilename,delopt);
%
% This script was developed to join freesurfer ctx-GM and ctx-WM regions. 
% The ctx-WM regions are generated using the command mri_aparc2aseg from
% freesurfer. This command allows the growing of gray matter regions into
% the white matter
%
% Input Parameters:
%   infilename     : Atlas filename
%   afile          : Freesurfer segmentation file
%   delopt         : Delete the original atlas file
%
% Output Parameters:
%   outfilename    : Output atlas filename.
%
% See also: 
%__________________________________________________
% Authors:  Yasser Aleman Gomez 
% LIM
% February 27th 2012
% Version $1.0
if nargin<3 
    delopt = 0;
end
if nargin<2 
    afile = 'aparc+aseg';
end
switch afile
    case 'aparc+aseg'
        V = spm_vol(infilename);
        IA = spm_read_vols(V);
        ind = find((IA(:)>=3000)&(IA(:)<3100));  %Left Hemisphere
        ind2 = find((IA(:)>=4000)&(IA(:)<4100)); %Right Hemisphere
        ind = [ind;ind2];
        It =IA;It(ind) = It(ind)-2000;       
        [pth,nm,ext] = fileparts(infilename);
        if delopt == 0
            V.fname = [pth filesep nm '_gwj.nii'];
        end
        spm_write_vol(V,It);
        outfilename = V.fname;
    case 'a2005s+aseg'
        V = spm_vol(infilename);
        IA = spm_read_vols(V);
        ind = find((IA(:)>=3100)&(IA(:)<3200));  %Left Hemisphere
        ind2 = find((IA(:)>=4100)&(IA(:)<4200)); %Right Hemisphere
        ind = [ind;ind2];
        It =IA;It(ind) = It(ind)-2000;       
        [pth,nm,ext] = fileparts(infilename);
        if delopt == 0
            V.fname = [pth filesep nm '_gwj.nii'];
        end
        spm_write_vol(V,It);
        outfilename = V.fname;
    case 'a2009s+aseg'
        V = spm_vol(infilename);
        IA = spm_read_vols(V);
        ind = find((IA(:)>=13100)&(IA(:)<13200));  %Left Hemisphere
        ind2 = find((IA(:)>=14100)&(IA(:)<14200)); %Right Hemisphere
        ind = [ind;ind2];
        It =IA;It(ind) = It(ind)-2000;       
        [pth,nm,ext] = fileparts(infilename);
        if delopt == 0
            V.fname = [pth filesep nm '_gwj.nii'];
        end
        spm_write_vol(V,It);
        outfilename = V.fname;
end
return;
