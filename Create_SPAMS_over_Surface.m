function Outfile = Create_SPAMS_over_Surface(FreeSDir,Idfile, Outfile)

FreeSDir = '/media/Data/PROCESSING_RESULTS/5-freesurfer_processing';
%Idfile = '/media/MyDisk/PROCESSING_RESULTS/5-freesurfer_processing/HCPrest_Ids.txt';
Idfile = '/media/MyDisk/PROCESSING_RESULTS/5-freesurfer_processing/HCPwmPR1_Ids.txt';
if exist(Idfile,'file')
    Ids = char(textread(Idfile,'%s'));
else
    Ids = Idfile;
end
Ns = size(Ids,1);
failed = '';
cad = '';
Ido = '';
aparcId = 'fsaverage.aparc2sulci';
Surfl = Read_Surface('/media/Data/PROCESSING_RESULTS/HCP/5-freesurfer_processing/fsaverage/surf/lh.inflated');
Surfr = Read_Surface('/media/Data/PROCESSING_RESULTS/HCP/5-freesurfer_processing/fsaverage/surf/rh.inflated');
sPamsDir = '/media/UserBackup/Publicaciones/Articulo_SulcalLines/Figures/SPAMs/';
% % % % for  i = 1:Ns
% % % %     
% % % %     Id = deblank(Ids(i,:));
% % % %     Id
% % % %     lhannot = [FreeSDir filesep Id filesep 'label' filesep 'lh.' aparcId '.annot'];
% % % %     [txt,ctab] = read_cfiles(lhannot);
% % % %     txt(txt<0) = 0;
% % % % %     opts.fill = 0;
% % % % %     opts.iterations = 3;
% % % % %     [dilLabels] = Dilate_Surface_Label(Surf, txt, opts);
% % % %     
% % % %     allAnnotL(:,i) = txt;
% % % %    % Surf.Is = txt;
% % % % %     
% % % %     rhannot = [FreeSDir filesep Id filesep 'label' filesep 'rh.' aparcId '.annot'];
% % % %     [txt,ctab] = read_cfiles(rhannot);
% % % %     txt(txt<0) = 0;
% % % %     allAnnotR(:,i) = txt;
% % % % end
load('/media/COSAS/scripts/matlab.mat');

sulcIds = [660701;3988703;717055;11842623;11822181;6609981;14423183;13143061;15733781;1367261];
%sulcIds = [660701;3988703;717055;11842623;11822181;6609981;14423183;13143061];
%names = {'centralsulcus';'superiortemporalsulcus';'cingulatesulcus';'calcarinefissure';'parietooccipitalfissure';'superiorfrontalsulcus';'intraparietalsulcus';'postcentralsulcus'};
names = {'centralsulcus';'superiortemporalsulcus';'cingulatesulcus';'calcarinefissure';'parietooccipitalfissure';'superiorfrontalsulcus';'intraparietalsulcus';'postcentralsulcus';'precentralsulcus';'inferiorfrontalsulcus'};


% ind2rem = find(ismember(allAnnotL, sulcIds) == 0);
% allAnnotL(ind2rem) = 0;
% load('/media/COSAS/scripts/matlab.mat');

% ind2rem = find(ismember(allAnnotR, sulcIds) == 0);
% allAnnotL(ind2rem) = 0;
thr = .35;
%Surf.SurfData.FaceVertexCData = repmat([1 1 1]*0.8,[size(Surf.SurfData.vertices,1) 1]);
[Tripl] = Vert_Neibp(double(Surfl.SurfData.faces),size(Surfl.SurfData.vertices,1),size(Surfl.SurfData.faces,1));
Temp = sum(Tripl);
Tripl(:,Temp==0) = [];
templ = Tripl(:,3:end);
indzl = find(templ == 0);
templ(indzl) = 1;

[Tripr] = Vert_Neibp(double(Surfr.SurfData.faces),size(Surfr.SurfData.vertices,1),size(Surfr.SurfData.faces,1));
Temp = sum(Tripr);
Tripr(:,Temp==0) = [];
tempr = Tripr(:,3:end);
indzr = find(tempr == 0);
tempr(indzr) = 1;

for i = 4:length(sulcIds)
    
    i
    
    
    
    %% Left Hemisphere
    T = allAnnotL*0;
    ind = find(allAnnotL == sulcIds(i));
    T(ind) = 1;
    T =sum(T,2)/Ns;

    
%      temp1 = T(templ);temp1(indzl) = 0;
%      T = sum(temp1,2)./(sum(logical(Tripl(:,3:end)),2)+eps);
     Torig = T;
     ind2rem = find(T < thr);
     T(ind2rem) = 0;
    
    ind = find(T);
    
    Surfl.Is = T*0;
    Surfl.Is(ind) = 1;
    indzeros = find(Surfl.Is == 0);
    Surfl.Is(indzeros) = max(Surfl.Is) + 1;
    [Surfl] = Surf_Corr(Surfl);
    indzeros = find(Surfl.Is == max(Surfl.Is));
    Surfl.Is(indzeros) = 0;   
    
    opts.iterations = 2;
    opts.fill = 0;
    Surfl.Is = Dilate_Surface_Label(Surfl, Surfl.Is, opts);
    T = Torig.*Surfl.Is;
    
%      temp1 = T(templ);temp1(indzl) = 0;
%      T = sum(temp1,2)./(sum(logical(Tripl(:,3:end)),2)+eps);
     
    
    fname = [sPamsDir 'lh.' aparcId '.' char(names{i})];
    save_char('/tmp/ltmp.map', T, size(Surfl.SurfData.faces,1));
    try
        copyfile('/tmp/ltmp.map',fname);
    end
    delete('/tmp/ltmp.map');
    
    smoothlhchar = [sPamsDir 'lh.' aparcId '.smooth.' char(names{i})];
%     cad = ['mri_surf2surf --srcsubject ' trgId ' --trgsubject ' trgId ' --hemi lh --sval ' fname ' --tval ' smoothlhchar   ' --nsmooth-in 2' ];
%     system(cad);
    
    %% Right Hemisphere
    T = allAnnotR*0;
    ind = find(allAnnotR == sulcIds(i));
    T(ind) = 1;
    T =sum(T,2)/Ns;

    
%      temp1 = T(tempr);temp1(indzr) = 0;
%      T = sum(temp1,2)./(sum(logical(Tripr(:,3:end)),2)+eps);
     Torig = T;
     ind2rem = find(T < thr);
     T(ind2rem) = 0;
    
    ind = find(T);
    
    Surfr.Is = T*0;
    Surfr.Is(ind) = 1;
    indzeros = find(Surfr.Is == 0);
    Surfr.Is(indzeros) = max(Surfr.Is) + 1;
    [Surfr] = Surf_Corr(Surfr);
    indzeros = find(Surfr.Is == max(Surfr.Is));
    Surfr.Is(indzeros) = 0;   
    
    opts.iterations = 2;
    opts.fill = 0;
    Surfr.Is = Dilate_Surface_Label(Surfr, Surfr.Is, opts);
    T = Torig.*Surfr.Is;
    
%      temp1 = T(tempr);temp1(indzr) = 0;
%      T = sum(temp1,2)./(sum(logical(Tripr(:,3:end)),2)+eps);

    
    
    
    fname = [sPamsDir 'rh.' aparcId '.' char(names{i})];
    save_char('/tmp/rtmp.map', T, size(Surfr.SurfData.faces,1));
    try
        copyfile('/tmp/rtmp.map',fname);
    end
    delete('/tmp/rtmp.map');
    
    
    smoothlhchar = [sPamsDir 'rh.' aparcId '.smooth.' char(names{i})];
%     cad = ['mri_surf2surf --srcsubject ' trgId ' --trgsubject ' trgId ' --hemi rh --sval ' fname ' --tval ' smoothlhchar   ' --nsmooth-in 2' ];
%     system(cad);
end
return;



