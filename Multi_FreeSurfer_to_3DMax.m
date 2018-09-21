function Multi_FreeSurfer_to_3DMax
SurfFiles = '/usr/local/freesurfer/subjects/fsaverage/surf/lh.pial';
CharacFiles = '/media/COSAS/Test/freesurfer/fsaverage/label/lh.aparc.annot';
talfile = ['/media/COSAS/Test/freesurfer/fsaverage' filesep 'mri' filesep 'transforms' filesep 'talairach.lta'];

customplot = 0;
Ns = size(SurfFiles,1);
for i =1:Ns
    SurfFile = deblank(SurfFiles(i,:));
    CharacFile = deblank(CharacFiles(i,:));
    [pth,nm,ext] = fileparts(SurfFile);
    nm = [deblank(nm) deblank(ext(2:end))];
    ind = strfind(nm,'.');nm(ind) = [];
%     OBJFile = [pth filesep nm '_Stats.obj'];
    OBJFile = ['/media/UserBackup/Erika_Left.obj'];
    [txt,colortable] = read_cfiles(CharacFile);
    %[vertices, faces] = freesurfer_read_surf(SurfFile);
    
%     
%     cras1 = textread(talfile,'%s',5,'headerlines',20);
%     cras = char(cras1);
%     cras = [str2num(cras(3,:))  str2num(cras(4,:)) str2num(cras(5,:))];
    [OutFiles, SurfF] = Exp_Surf(SurfFile, '0', '','', 'imp','n');
    Surf= SurfF{1};
    Surf.Name = 'LH.Pial';
%     Surf.SurfData.vertices = Surf.SurfData.vertices + repmat(cras,[size(Surf.SurfData.vertices,1) 1]);
    %
    %
    %
%     
%     
%     Surf.SurfData.vertices = vertices;
%     Surf.SurfData.faces = faces;
%     Surf.Name = 'Pial';
    Surf.Is = txt;
    if colortable.table == 0
        col = [1 0 0;0 1 0; 0 0 1; 1 1 0;0 1 1;1 0 1;1 0.5 0;0 0.5 1;1 0 0.5;0.5 1 0; 0.25 0.5 1;0.6 0.3 0.52;0.7 0.5 0.9];
        Ncolor = size(col,1);
        structs = unique(Surf.Is);
        Nst = length(structs);
        re = floor(Nst/Ncolor); 
        col = repmat(col,[re+1 1])*255;
        col = floor(col(1:Nst,:));
        colortable.table = [col repmat(0,[Nst 1]) structs];
        for i = 1:Nst
            colortable.struct_names{i,1} = sprintf('Region_%.5d',structs(i));
        end
        ind = find(structs == 0);
        if ~isempty(ind)
            colortable.table(ind,:) = repmat([0 0 0 0 0 ],[length(ind) 1]);
        end
    end
    colortable.table(:,6) = ones(size(colortable.table,1),1);
    if customplot == 1
        TextFile = '/media/Data/yasser/ch_brainvisa/Results/aparc_perc_WM_thickness_change.txt';
        [Num, Snames, Vals] = textread(TextFile,'%u%s%f','delimiter',' ');
        Nnames = size(Snames,1);
        Scolor = [0.5 0.5 0.5];
        Colors = repmat([Scolor],[size(colortable.table,1) 1]);
        colortable.table(:,1:3) = round(Colors*255);
        [Colors] = Surf_Color_custom(Vals,'hsv');
        names = char(colortable.struct_names);
        colortable.table
        for i = 1:Nnames
            ind = find(ismember(names(:,1:length(char(Snames{i,1}))),char(Snames{i,1}),'rows'));
            colortable.table(ind,1:3) = round(Colors(i,:)*255) ;
            colortable.table
        end
    end
    OBJFile = FreeSurfer_to_3DMax(Surf, colortable, OBJFile);
end
return

function [Colors] = Surf_Color_custom(txt,cl);
% Syntax :
% [Colors] = Surf_Color(txt,cl);
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