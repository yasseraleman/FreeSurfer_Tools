function Plot_SPAMs

aparcId = 'fsaverage.aparc2sulci';
plotthr = 0.1;
Surfl = Read_Surface('/media/Data/PROCESSING_RESULTS/HCP/5-freesurfer_processing/fsaverage/surf/lh.inflated','curv');
curvmapl = Surfl.Is;
Surfr = Read_Surface('/media/Data/PROCESSING_RESULTS/HCP/5-freesurfer_processing/fsaverage/surf/rh.inflated','curv');
curvmapr = Surfr.Is;
sPamsDir = '/media/UserBackup/Publicaciones/Articulo_SulcalLines/Figures/SPAMs/'
sulcIds = [660701;3988703;717055;11842623;11822181;6609981;14423183;13143061;15733781;1367261];
%sulcIds = [660701;3988703;717055;11842623;11822181;6609981;14423183;13143061];
%names = {'centralsulcus';'superiortemporalsulcus';'cingulatesulcus';'calcarinefissure';'parietooccipitalfissure';'superiorfrontalsulcus';'intraparietalsulcus';'postcentralsulcus'};
names = {'centralsulcus';'superiortemporalsulcus';'cingulatesulcus';'calcarinefissure';'parietooccipitalfissure';'superiorfrontalsulcus';'intraparietalsulcus';'postcentralsulcus';'precentralsulcus';'inferiorfrontalsulcus'};

lhFiles = [repmat([sPamsDir filesep 'lh.' aparcId '.' ],[length(names) 1]) char(names)];
rhFiles = [repmat([sPamsDir filesep 'rh.' aparcId '.' ],[length(names) 1]) char(names)];

rawColors = [255 0 0;255 255 0;255 150 0;0 0 255;0 255 0;0 255 255;255 0 255;50 127 92;116    10   255; 255   204   153];

% Left Hemisphere
TempColors = zeros(size(Surfl.SurfData.vertices,1),1);
sCales= zeros(size(Surfl.SurfData.vertices,1),1);
%% =================== Computing Neighoborhoods ======================== %%
[Trip] = Vert_Neibp(double(Surfl.SurfData.faces),size(Surfl.SurfData.vertices,1),size(Surfl.SurfData.faces,1));
Temp = sum(Trip);
Trip(:,Temp==0) = [];
temp = Trip(:,3:end);
indz = find(temp == 0);
temp(indz) = 1;

for i = 1:length(sulcIds)
    T = read_cfiles(deblank(lhFiles(i,:)));    
%     temp1 = T(temp);temp1(indz) = 0;
%     T = sum(temp1,2)./(sum(logical(temp1),2)+eps);
    
  
%     ind = find(T>plotthr);
%     Surfl.Is = T*0;
%     Surfl.Is(ind) = 1;
%     indzeros = find(Surfl.Is == 0);
%     Surfl.Is(indzeros) = max(Surfl.Is) + 1;
%     [Surfl] = Surf_Corr(Surfl);
%     indzeros = find(Surfl.Is == max(Surfl.Is));
%     Surfl.Is(indzeros) = 0;
%     T = T.*Surfl.Is;
    ind = find(T);
    
    [ut,i1t,it2] = unique(T(ind));
    col = colorGradient([0 0 0]*.5,rawColors(i,:),size(ut,1));
    Colortemp = zeros(size(Surfl.SurfData.vertices,1),3) ;
    Colors = col(it2,:);
    Colortemp(ind,:) = Colors;
    
    
    TempColors = [TempColors Colortemp.*[logical(T) logical(T) logical(T)]];
    sCales = [sCales T];

end
if isfield(Surfl.SurfData, 'FaceVertexCData')
    Surfl.SurfData = rmfield(Surfl.SurfData, 'FaceVertexCData');
end
if isfield(Surfl,'Is')
    Surfl = rmfield(Surfl, 'Is');
end

TempColors(:,1) = [];
sCales(:,1) = [];

r= TempColors(:,1:3:end);
g= TempColors(:,2:3:end);
b= TempColors(:,3:3:end);

r = sum(sCales.*r,2)./(sum(sCales,2)+eps);
g = sum(sCales.*g,2)./(sum(sCales,2)+eps);
b = sum(sCales.*b,2)./(sum(sCales,2)+eps);

Surfl.SurfData.FaceVertexCData = [r g b];
temp = sum(sCales,2)./(sum(logical(sCales),2)+eps);
ind = find(temp < plotthr);
[ut,i1t,it2] = unique(temp(ind));[ut,i1t,it2] = unique(temp(ind));
col = colorGradient([1 1 1]*.8,[0 0 0],size(ut,1));
Surfl.SurfData.FaceVertexCData(ind,:) = col(it2,:);



ind = find(temp == 0);
indg = find(curvmapl <= 0); 
inds = find(curvmapl > 0);
curvColors =  zeros(size(Surfl.SurfData.vertices,1),3);
curvColors(indg,:) = repmat([1 1 1]*.8,[length(indg) 1]);
curvColors(inds,:) = repmat([1 1 1]*.5,[length(inds) 1]);
Surfl.SurfData.FaceVertexCData(ind,:) = curvColors(ind,:);


h = Plot_Surf(Surfl);
set(h(1),'SpecularExponent',1000000);
set(gcf,'name','New Figure','Color',[1 1 1],'units','centimeters','Position',[20 10 10 10],'InvertHardcopy','off','Visible','on');
axis off;
colorbar off;
set(h(2),'Color',[1 1 1]);
h2 = camlight;
set(h2,'Color',[1 1 1]);

% Lateral
view([270 0]);lightangle(h2,215, 0);lightangle(h(2),305, 0);
 
Figurename = [sPamsDir filesep 'Left_Hemisphere-Sulci_SPAMs-' aparcId '-Lateral.tiff'];
[pth, nm, ext] = fileparts(Figurename);
warning off;mkdir(pth);
export_fig(Figurename,'-tiff', gcf, '-nocrop','-transparent','-opengl','-r300' );

 %Medial
 view([90 0]);lightangle(h(2),45, 0);lightangle(h2,135, 0);set(h(2),'Color',[1 1 1]*.8); set(h2,'Color',[1 1 1]*.8);
 
 Figurename = [sPamsDir filesep 'Left_Hemisphere-Sulci_SPAMs-' aparcId '-Medial.tiff'];
[pth, nm, ext] = fileparts(Figurename);
warning off;mkdir(pth);
export_fig(Figurename,'-tiff', gcf, '-nocrop','-transparent','-opengl','-r300' );
close(gcf);


%% Right Hemisphere
TempColors = zeros(size(Surfr.SurfData.vertices,1),1);
sCales= zeros(size(Surfr.SurfData.vertices,1),1);

%% =================== Computing Neighoborhoods ======================== %%
[Trip] = Vert_Neibp(double(Surfr.SurfData.faces),size(Surfr.SurfData.vertices,1),size(Surfr.SurfData.faces,1));
Temp = sum(Trip);
Trip(:,Temp==0) = [];
temp = Trip(:,3:end);
indz = find(temp == 0);
temp(indz) = 1;

for i = 1:length(sulcIds)
    T = read_cfiles(deblank(rhFiles(i,:)));    
%     temp1 = T(temp);temp1(indz) = 0;
%     T = sum(temp1,2)./(sum(logical(temp1),2)+eps);
    
  
%     ind = find(T>plotthr);
%     Surfr.Is = T*0;
%     Surfr.Is(ind) = 1;
%     indzeros = find(Surfr.Is == 0);
%     Surfr.Is(indzeros) = max(Surfr.Is) + 1;
%     [Surfr] = Surf_Corr(Surfr);
%     indzeros = find(Surfr.Is == max(Surfr.Is));
%     Surfr.Is(indzeros) = 0;
%     T = T.*Surfr.Is;
    ind = find(T);
    
    [ut,i1t,it2] = unique(T(ind));
    col = colorGradient([0 0 0]*.5,rawColors(i,:),size(ut,1));
    Colortemp = zeros(size(Surfr.SurfData.vertices,1),3) ;
    Colors = col(it2,:);
    Colortemp(ind,:) = Colors;
    
    
    TempColors = [TempColors Colortemp.*[logical(T) logical(T) logical(T)]];
    sCales = [sCales T];

end
if isfield(Surfr.SurfData, 'FaceVertexCData')
    Surfr.SurfData = rmfield(Surfr.SurfData, 'FaceVertexCData');
end
if isfield(Surfr,'Is')
    Surfr = rmfield(Surfr, 'Is');
end

TempColors(:,1) = [];
sCales(:,1) = [];

r= TempColors(:,1:3:end);
g= TempColors(:,2:3:end);
b= TempColors(:,3:3:end);

r = sum(sCales.*r,2)./(sum(sCales,2)+eps);
g = sum(sCales.*g,2)./(sum(sCales,2)+eps);
b = sum(sCales.*b,2)./(sum(sCales,2)+eps);

Surfr.SurfData.FaceVertexCData = [r g b];
temp = sum(sCales,2)./(sum(logical(sCales),2)+eps);
ind = find(temp < plotthr);
[ut,i1t,it2] = unique(temp(ind));[ut,i1t,it2] = unique(temp(ind));
col = colorGradient([1 1 1]*.8,[0 0 0],size(ut,1));
Surfr.SurfData.FaceVertexCData(ind,:) = col(it2,:);



ind = find(temp == 0);
indg = find(curvmapr <= 0); 
inds = find(curvmapr > 0);
curvColors =  zeros(size(Surfr.SurfData.vertices,1),3);
curvColors(indg,:) = repmat([1 1 1]*.8,[length(indg) 1]);
curvColors(inds,:) = repmat([1 1 1]*.5,[length(inds) 1]);
Surfr.SurfData.FaceVertexCData(ind,:) = curvColors(ind,:);


h = Plot_Surf(Surfr);
set(h(1),'SpecularExponent',1000000);
set(gcf,'name','New Figure','Color',[1 1 1],'units','centimeters','Position',[20 10 10 10],'InvertHardcopy','off','Visible','off');
axis off;
colorbar off;
set(h(2),'Color',[1 1 1]);
h2 = camlight;
set(h2,'Color',[1 1 1]);

% Lateral
 view([90 0]);lightangle(h(2),45, 0);lightangle(h2,135, 0);
 
Figurename = [sPamsDir filesep 'Right_Hemisphere-Sulci_SPAMs-' aparcId '-Lateral.tiff'];
[pth, nm, ext] = fileparts(Figurename);
warning off;mkdir(pth);
export_fig(Figurename,'-tiff', gcf, '-nocrop','-transparent','-opengl','-r300' );

 %Medial
view([270 0]);lightangle(h(2),215, 0);lightangle(h2,305, 0); set(h(2),'Color',[1 1 1]*.8); set(h2,'Color',[1 1 1]*.8);

 
 Figurename = [sPamsDir filesep 'Right_Hemisphere-Sulci_SPAMs-' aparcId '-Medial.tiff'];
[pth, nm, ext] = fileparts(Figurename);
warning off;mkdir(pth);
export_fig(Figurename,'-tiff', gcf, '-nocrop','-transparent','-opengl','-r300' );
close(gcf);
return;