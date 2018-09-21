function Plot_Stat_results(Sfile,Afile,Txtfile,cl)
%% Para que Joost pueda plotear sobre una superficie de freesurfer los valores de p para cada region
%%
[vertices, faces] = freesurfer_read_surf(Sfile);Surf.SurfData.vertices = vertices;Surf.SurfData.faces = faces;Surf.Name = 'White';
if strfind(Sfile,'lh.')
    hemi = 'lh';
elseif strfind(Sfile,'rh.')
    hemi = 'rh';
end
%%   Leyendo Annot para dividir %%%
[label, labels, colortable] = read_annotation(Afile); % [txt] = read_cfiles(Oldannot);

%% Leyendo Text File %%%
txt = textread(Txtfile,'%f');
ind = find(isnan(txt) ==0);
Colors = repmat([0.51 0.51 0.51],[size(txt,1) 1]);
Surf.SurfData.FaceVertexCData = repmat([0.51 0.51 0.51],[size(labels,1) 1]);
txt2 = txt(ind);color2 = repmat([0.51 0.51 0.51],[size(txt2,1) 1]);
 [pth,nm,ext] = fileparts(cl);
if ~isempty(ext)
    colors = textread(cl,'%f');
    valcolors = reshape(colors,[size(colors,1)/4 4])';
    
    for i = 1:size(valcolors,1)
        ind = find(txt ==valcolors(i,1));
        ind = find(labels == colortable.table(ind,5));
        %Surf.Is = txt(i);
        Surf.SurfData.FaceVertexCData(ind,:) = repmat(valcolors(i,2:end),[size(ind,1) 1]);
    end
else
    %% Asignando colores
    [ut,i1t,it2] = unique(txt2);
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
    Colors(ind,:) = col(it2,:);
    for i = 1:size(colortable.table,1)
        ind = find(labels ==colortable.table(i,5));
        %Surf.Is = txt(i);
        Surf.SurfData.FaceVertexCData(ind,:) = repmat(Colors(i,:),[size(ind,1) 1]);
    end
end
%%    Subdividiendo Regiones

Surf.SurfData.FaceColor = 'interp';
nm = Surf.Name;
%% Plotting results
%== Lateral==%
h = figure('numbertitle','off','Color', 'white','name',['IBASPM Surface Ploting...:  Plot ' upper(hemi) ' _ Hemisphere']);
subplot(3,2,1);
strsurf=patch(Surf.SurfData,'edgecolor','none', 'tag','patch','facelighting','gouraud');
h=title(['Lateral View']);set(h,'FontSize',8,'FontName','Arial','Color',[0 0 0]);

if hemi == 'lh'
    view([270 0]);
    %     set(get(gca,'ZLabel'),'String','Anterior');
    %     set(gca,'ZTickLabel','');
    %     set(get(gca,'YLabel'),'String','Bottom');
    %     set(gca,'YTickLabel','');
elseif hemi == 'rh'
    view([90 0]);
    %     set(get(gca,'ZLabel'),'String','Posterior');
    %     set(gca,'ZTickLabel','');
    %    % set(gca,'ZColor',[0 0 0]);
    %     set(get(gca,'YLabel'),'String','Bottom');
    %     set(gca,'YTickLabel','');
    %set(gca,'YColor',[0 0 0]);
end
axis off;
% set(gca,'ZColor',[0 0 0]);
% set(gca,'YColor',[0 0 0]);
axis tight;axis equal;
camlight;
%== Medial==%
subplot(3,2,2);
strsurf=patch(Surf.SurfData,'edgecolor','none', 'tag','patch','facelighting','gouraud');
h=title(['Medial View']);set(h,'FontSize',8,'FontName','Arial','Color',[0 0 0]);
if hemi == 'lh'
    view([90 0]);
    %     set(get(gca,'ZLabel'),'String','Posterior');
    %     set(gca,'ZTickLabel','');
    %     set(get(gca,'YLabel'),'String','Bottom');
    %     set(gca,'YTickLabel','');
elseif hemi == 'rh'
    view([270 0]);
    %     set(get(gca,'ZLabel'),'String','Anterior');
    %     set(gca,'ZTickLabel','');
    %     set(get(gca,'YLabel'),'String','Bottom');
    %     set(gca,'YTickLabel','');
end
axis off;
axis tight;axis equal;
camlight;
%== Bottom==%
subplot(3,2,3);
strsurf=patch(Surf.SurfData,'edgecolor','none', 'tag','patch','facelighting','gouraud');
h=title(['Bottom View']);set(h,'FontSize',8,'FontName','Arial','Color',[0 0 0]);
if hemi == 'lh'
    view([180 270]);
    %     set(get(gca,'XLabel'),'String','Posterior');
    %     set(gca,'XTickLabel','');
    %     set(get(gca,'YLabel'),'String','Left');
    %     set(gca,'YTickLabel','');
elseif hemi == 'rh'
    view([180 270]);
    %     set(get(gca,'XLabel'),'String','Posterior');
    %     set(gca,'XTickLabel','');
    %     set(get(gca,'YLabel'),'String','Right');
    %     set(gca,'YTickLabel','');
end
axis off;
axis tight;axis equal; %h=title(['Lateral View']);set(h,'FontSize',15,'FontName','Arial');
camlight;
%== Top==%
subplot(3,2,4);
strsurf=patch(Surf.SurfData,'edgecolor','none', 'tag','patch','facelighting','gouraud');
h=title(['Top View']);set(h,'FontSize',8,'FontName','Arial','Color',[0 0 0]);
if hemi == 'lh'
    view([0 90]);
    %     set(get(gca,'XLabel'),'String','Posterior');
    %     set(gca,'XTickLabel','');
    %     set(get(gca,'YLabel'),'String','Left');
    %     set(gca,'YTickLabel','');
elseif hemi == 'rh'
    view([0 90]);
    %     set(get(gca,'XLabel'),'String','Posterior');
    %     set(gca,'XTickLabel','');
    %     set(get(gca,'YLabel'),'String','Left');
    %     set(gca,'YTickLabel','');
end
axis off;
axis tight;axis equal;
camlight;

%== Front==%
subplot(3,2,5);
strsurf=patch(Surf.SurfData,'edgecolor','none', 'tag','patch','facelighting','gouraud');
h=title(['Frontal View']);set(h,'FontSize',8,'FontName','Arial','Color',[0 0 0]);
if hemi == 'lh'
    view([180 0]);
elseif hemi == 'rh'
    view([180 0]);
end
axis off;
axis tight;axis equal; %h=title(['Lateral View']);set(h,'FontSize',15,'FontName','Arial');
camlight;
%== Back==%
subplot(3,2,6);
strsurf=patch(Surf.SurfData,'edgecolor','none', 'tag','patch','facelighting','gouraud');
h=title(['Posterior View']);set(h,'FontSize',8,'FontName','Arial','Color',[0 0 0]);
if hemi == 'lh'
    view([0 0]);
elseif hemi == 'rh'
    view([0 0]);
end
axis off;
axis tight;axis equal;
camlight;
if ~isempty(ext)
else
    h = colorbar;
    if cl ~= 'spectral'
        colormap(cl);
    else
        colormap(spectral(64));
    end
    range = max(txt2)-min(txt2);values =  min(txt2):range/5:max(txt2);
    set(h,'YTickLabel',num2str(values'));
end
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










