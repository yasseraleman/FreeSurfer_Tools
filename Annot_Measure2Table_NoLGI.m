function OutStatFile = Annot_Measure2Table_NoLGI(FreeSDir, Idfile, OutStatFile, annotIds, thtype);
%
% Syntax :
% OutStatFile = Annot_Measure2Table(FreeSDir, Idfile, OutStatFile, annotIds, thtype);
%
% This script obtains different morphometric measures for a specified
% annotation file.
%
%
% Input Parameters:
%       FreeSDir                : FreeSurfer Directory
%       Idfile                  : Text file containing the Ids List
%     OutStatFile               : Output statistics file
%       annotIds                : Annotation Files
%       thtype                  : Cortical thickness type
%
% Output Parameters:
%      OutFiles                 : Output stat files
%
% See also:
%__________________________________________________
% Authors: Yasser Aleman Gomez
% LIM, HUGGM
% March 22th 2013
% Version $1.0


%% =======================  FreeSurfer IDs  ============================ %%
%Idfile = '/media/Data/PROCESSING_RESULTS/PCAMPO/5-freesurfer_processing/Ids.txt';
if exist(Idfile,'file')
    Ids = char(textread(Idfile,'%s'));
else
    Ids = Idfile;
end
 %FreeSDir = '/media/Data/PROCESSING_RESULTS/PCAMPO/5-freesurfer_processing';
 %OutStatFile = '/media/Data/PROCESSING_RESULTS/PCAMPO/5-freesurfer_processing//AnnotStat-21-7-2014.txt';
%  Ids = '0822x-20110430';
%% =============== End of Detecting FreeSurfer IDs  ==================== %%
if nargin <4
    annotIds = 'aparc';
end
if nargin <5
    thtype = 'thickness';
end
switch thtype
    case 'thickness'
        Ident = 'cth_freesurfer';      % Cortical thickness computed via freesurfer
    case 'segload.thickness'
        Ident = 'cth_direct_segload';  % Cortical thickness computed via Direct and segmentations obtained from segload
    case 'fslfast.thickness'
        Ident = 'cth_direct_fast';     % Cortical thickness computed via Direct and segmentations obtained from fast
    case 'spmvbm8.thickness'
        Ident = 'cth_direct_vbm8';     % Cortical thickness computed via Direct and segmentations obtained from spm
end

%% ==================== Lobes Parcellation ============================== %
% Selecting Frontal regions
FroIds = sort([1028 1003 1027 1018 1019 1020 1012 1014 1024 1017 1032 1026 1002]);
FroIdsR = FroIds + 1000;

% Selecting Temporal regions
TempIds = sort([1009 1015 1030 1001 1007 1034 1006 1033 1016]);
TempIdsR = TempIds + 1000;

% Selecting Parietal regions
ParIds = sort([1029 1008 1031 1022 1025 1010 1023]);
ParIdsR = ParIds + 1000;

% Selecting Occipital regions
OccIds = sort([1011 1013 1005 1021]);
OccIdsR = OccIds + 1000;

% Selecting Insula regions
InsIds = [1035];
InsIdsR = [2035];

%% ================ End of Lobes Parcellation =========================== %


%% ================== Subjects Processing ============================== %%
tempd = '';
Nc = size(annotIds,1);
OutFiles = '';
for i = 1:Nc
    aparcId = deblank(annotIds(i,:));
    
    Nsubj = size(Ids,1);
    for z= 1:Nsubj
        subjId = deblank(Ids(z,:));
        disp(strvcat(' ',' '));
        disp(['Processing =======>  Annot Files: ' num2str(i) ' of ' num2str(Nc) ' . =====> Subject ID: ' subjId ' . ---  ' num2str(z) ' of ' num2str(Nsubj)]);
        
        %% =================== Reading files ================================ %
        % Parcellation Maps
        lhannot = [FreeSDir filesep subjId filesep 'label' filesep 'lh.' aparcId '.annot'];
        rhannot = [FreeSDir filesep subjId filesep 'label' filesep 'rh.' aparcId '.annot'];
        Aparclhannot = [FreeSDir filesep subjId filesep 'label' filesep 'lh.aparc.annot'];
        Aparcrhannot = [FreeSDir filesep subjId filesep 'label' filesep 'rh.aparc.annot'];
        
        
        % Surfaces
        lhsurf = [FreeSDir filesep subjId filesep 'surf' filesep 'lh.pial'];
        rhsurf = [FreeSDir filesep subjId filesep 'surf' filesep 'rh.pial'];
        
        % Curvature Maps
        lhcurvf = [FreeSDir filesep subjId filesep 'surf' filesep 'lh.curv'];
        rhcurvf = [FreeSDir filesep subjId filesep 'surf' filesep 'rh.curv'];
        
        % Thickness Maps
        lhcthf = [FreeSDir filesep subjId filesep 'surf' filesep 'lh.' thtype];
        rhcthf = [FreeSDir filesep subjId filesep 'surf' filesep 'rh.' thtype];
        
        % Volumetric Atlas
        atlasfile = [FreeSDir filesep subjId filesep 'mri' filesep aparcId '+aseg.mgz'];
        
        
        %% =================== End of Reading files ========================= %
        
        
        
        %% ================= Computing Things =================================== %
        %----- Reading Surfaces
        [OutFiles, SurfF] = Exp_Surf(lhsurf, '0', '','', 'imp','n');
        Surfwl= SurfF{1};
        Surfwl.Name = 'LH.WHITE';
        
        [OutFiles, SurfF] = Exp_Surf(rhsurf, '0', '','', 'imp','n');
        Surfwr= SurfF{1};
        Surfwr.Name = 'RH.White';
        
        %----- Reading Annot Files
        [txtl,ctabl] = read_cfiles(lhannot);
        ctabl.table = [ctabl.table 1000+[0:size(ctabl.table,1)-1]' ];
        tempname = char(ctabl.struct_names);
        indu = find(ismember(tempname(:,1:7),'unknown','rows') == 1);
        indcc = find(ismember(tempname(:,1:14),'corpuscallosum','rows') == 1);
        ctabl.table([indu indcc],:) = [];
        ctabl.struct_names([indu indcc]) = [];
        [txtr,ctabr] = read_cfiles(rhannot);
        ctabr.table = [ctabr.table 2000+[0:size(ctabr.table,1)-1]' ];
        tempname = char(ctabr.struct_names);
        indu = find(ismember(tempname(:,1:7),'unknown','rows') == 1);
        indcc = find(ismember(tempname(:,1:14),'corpuscallosum','rows') == 1);
        ctabr.table([indu indcc],:) = [];
        ctabr.struct_names([indu indcc]) = [];
        
        %----- Reading Curvature
        [left_curv,ctab] = read_cfiles(lhcurvf);
        left_curv = abs(left_curv);
        [right_curv,ctab] = read_cfiles(rhcurvf);
        right_curv = abs(right_curv);
        
        %----- Reading Thickness Map
        [left_cth,ctab] = read_cfiles(lhcthf);
        [right_cth,ctab] = read_cfiles(rhcthf);
        
        
        % --- Converting Volumetric Image
        if ~exist([FreeSDir filesep subjId filesep 'tmp' filesep aparcId '+aseg.nii'],'file')
            cad = ['mri_convert -i ' atlasfile ' -o ' FreeSDir filesep subjId filesep 'tmp' filesep aparcId '+aseg.nii'];
            system(cad);
        end
        VA =  spm_vol([FreeSDir filesep subjId filesep 'tmp' filesep aparcId '+aseg.nii']);
        IA = spm_read_vols(VA);
        temp = nonzeros(IA(:));
        c = accumarray(temp,ones(length(temp),1));
        voxsize = sqrt(sum(VA.mat(1:3,1:3).^2));
        
        
        % ======================== Left Hemisphere ===============================%
        
        % ------ Gray Matter regions
        Nr = size(ctabl.table,1);
        fv.vertices = Surfwl.SurfData.vertices;
        for j = 1:Nr
            
            ind = find(txtl == ctabl.table(j,5));
            lnvert(j,1) = length(ind); % Number of vertices
            lcthm(j,1) = mean(left_cth(ind)); % Mean Cortical Thickness
            lcths(j,1) = std(left_cth(ind)); % Std Cortical Thickness
            lcurv(j,1) = mean(left_curv(ind)); % Mean curvature
            
            %----- Computing Area
            indf = find(sum(ismember( Surfwl.SurfData.faces,ind)')==3);
            At = 0;
            fv.faces = Surfwl.SurfData.faces(indf,:);
            
            At = Comp_Surf_Area_custom(fv);
            larea(j,1) = real(At);% mm^2
        end
        %----- Computing Volume
        lvols = c(ctabl.table(:,6))*prod(voxsize);
        
        % ------ Brain and Hemisphere
        ind = find(ismember(ctabl.table(:,6),[FroIds(:);ParIds(:);TempIds(:);OccIds(:);InsIds(:)]));
        lharea(1,1) = sum(larea(ind,1)); % Area
        Sids = ctabl.table(ind,5);
        ind = find(ismember(txtl,Sids));
        lhnvert(1,1) = length(ind);
        lhcthm(1,1) = mean(left_cth(ind)); % Mean Cortical Thickness
        lhcths(1,1) = std(left_cth(ind)); % Std Cortical Thickness
        lhcurv(1,1) = mean(left_curv(ind)); % Mean curvature
        lhvols(1,1) = sum(c([FroIds(:);ParIds(:);TempIds(:);OccIds(:);InsIds(:)])*prod(voxsize)); % Volume
        
        % ------ Lobes ----- Frontal
        ind = find(ismember(ctabl.table(:,6),FroIds(:)));
        llobearea(1,1) = sum(larea(ind,1)); % Area
        Sids = ctabl.table(ind,5);
        ind = find(ismember(txtl,Sids));
        llobenvert(1,1) = length(ind);
        llobecthm(1,1) = mean(left_cth(ind)); % Mean Cortical Thickness
        llobecths(1,1) = std(left_cth(ind)); % Std Cortical Thickness
        llobecurv(1,1) = mean(left_curv(ind)); % Mean curvature
        llobevols(1,1) = sum(c(FroIds(:))*prod(voxsize)); % Volume
        
        % ------ Lobes ----- Parietal
        ind = find(ismember(ctabl.table(:,6),ParIds(:)));
        llobearea(2,1) = sum(larea(ind,1)); % Area
        Sids = ctabl.table(ind,5);
        ind = find(ismember(txtl,Sids));
        llobenvert(2,1) = length(ind);
        llobecthm(2,1) = mean(left_cth(ind)); % Mean Cortical Thickness
        llobecths(2,1) = std(left_cth(ind)); % Std Cortical Thickness
        llobecurv(2,1) = mean(left_curv(ind)); % Mean curvature
        llobevols(2,1) = sum(c(ParIds(:))*prod(voxsize)); % Volume
        
        % ------ Lobes ----- Temporal
        ind = find(ismember(ctabl.table(:,6),TempIds(:)));
        llobearea(3,1) = sum(larea(ind,1)); % Area
        Sids = ctabl.table(ind,5);
        ind = find(ismember(txtl,Sids));
        llobenvert(3,1) = length(ind);
        llobecthm(3,1) = mean(left_cth(ind)); % Mean Cortical Thickness
        llobecths(3,1) = std(left_cth(ind)); % Std Cortical Thickness
        llobecurv(3,1) = mean(left_curv(ind)); % Mean curvature
        llobevols(3,1) = sum(c(TempIds(:))*prod(voxsize)); % Volume
        
        % ------ Lobes ----- Occipital
        ind = find(ismember(ctabl.table(:,6),OccIds(:)));
        llobearea(4,1) = sum(larea(ind,1)); % Area
        Sids = ctabl.table(ind,5);
        ind = find(ismember(txtl,Sids));
        llobenvert(4,1) = length(ind);
        llobecthm(4,1) = mean(left_cth(ind)); % Mean Cortical Thickness
        llobecths(4,1) = std(left_cth(ind)); % Std Cortical Thickness
        llobecurv(4,1) = mean(left_curv(ind)); % Mean curvature
        llobevols(4,1) = sum(c(OccIds(:))*prod(voxsize)); % Volume
        
        % ======================End of Left Hemisphere ==========================%
        
        
        % ======================== Right Hemisphere ===============================%
        Nr = size(ctabr.table,1);
        fv.vertices = Surfwr.SurfData.vertices;
        for j = 1:Nr
            ind = find(txtr == ctabr.table(j,5));
            rnvert(j,1) = length(ind); % Number of vertices
            rcthm(j,1) = mean(right_cth(ind)); % Mean Cortical Thickness
            rcths(j,1) = std(right_cth(ind)); % Std Cortical Thickness
            rcurv(j,1) = mean(right_curv(ind)); % Mean curvature

            %----- Computing Area
            indf = find(sum(ismember( Surfwr.SurfData.faces,ind)')==3);
            At = 0;
            fv.faces = Surfwr.SurfData.faces(indf,:);
            At = Comp_Surf_Area_custom(fv);
            rarea(j,1) = real(At);% mm^2
        end
        %----- Computing Volume
        rvols = c(ctabr.table(:,6))*prod(voxsize);
        
        % ------ Brain and Hemisphere
        ind = find(ismember(ctabr.table(:,6),[FroIdsR(:);ParIdsR(:);TempIdsR(:);OccIdsR(:);InsIdsR(:)]));
        rharea(1,1) = sum(rarea(ind,1)); % Area
        Sids = ctabr.table(ind,5);
        ind = find(ismember(txtr,Sids));
        rhnvert(1,1) = length(ind);
        rhcthm(1,1) = mean(right_cth(ind)); % Mean Cortical Thickness
        rhcths(1,1) = std(right_cth(ind)); % Std Cortical Thickness
        rhcurv(1,1) = mean(right_curv(ind)); % Mean curvature
        rhvols(1,1) = sum(c([FroIdsR(:);ParIdsR(:);TempIdsR(:);OccIdsR(:);InsIdsR(:)])*prod(voxsize)); % Volume
        
        % ------ Lobes ----- Frontal
        ind = find(ismember(ctabr.table(:,6),FroIdsR(:)));
        rlobearea(1,1) = sum(rarea(ind,1)); % Area
        Sids = ctabr.table(ind,5);
        ind = find(ismember(txtr,Sids));
        rlobenvert(1,1) = length(ind);
        rlobecthm(1,1) = mean(right_cth(ind)); % Mean Cortical Thickness
        rlobecths(1,1) = std(right_cth(ind)); % Std Cortical Thickness
        rlobecurv(1,1) = mean(right_curv(ind)); % Mean curvature
        rlobevols(1,1) = sum(c(FroIdsR(:))*prod(voxsize)); % Volume
        
        % ------ Lobes ----- Parietal
        ind = find(ismember(ctabr.table(:,6),ParIdsR(:)));
        rlobearea(2,1) = sum(rarea(ind,1)); % Area
        Sids = ctabr.table(ind,5);
        ind = find(ismember(txtr,Sids));
        rlobenvert(2,1) = length(ind);
        rlobecthm(2,1) = mean(right_cth(ind)); % Mean Cortical Thickness
        rlobecths(2,1) = std(right_cth(ind)); % Std Cortical Thickness
        rlobecurv(2,1) = mean(right_curv(ind)); % Mean curvature
        rlobevols(2,1) = sum(c(ParIdsR(:))*prod(voxsize)); % Volume
        
        % ------ Lobes ----- Temporal
        ind = find(ismember(ctabr.table(:,6),TempIdsR(:)));
        rlobearea(3,1) = sum(rarea(ind,1)); % Area
        Sids = ctabr.table(ind,5);
        ind = find(ismember(txtr,Sids));
        rlobenvert(3,1) = length(ind);
        rlobecthm(3,1) = mean(right_cth(ind)); % Mean Cortical Thickness
        rlobecths(3,1) = std(right_cth(ind)); % Std Cortical Thickness
        rlobecurv(3,1) = mean(right_curv(ind)); % Mean curvature
        rlobevols(3,1) = sum(c(TempIdsR(:))*prod(voxsize)); % Volume
        
        % ------ Lobes ----- Occipital
        ind = find(ismember(ctabr.table(:,6),OccIdsR(:)));
        rlobearea(4,1) = sum(rarea(ind,1)); % Area
        Sids = ctabr.table(ind,5);
        ind = find(ismember(txtr,Sids));
        rlobenvert(4,1) = length(ind);
        rlobecthm(4,1) = mean(right_cth(ind)); % Mean Cortical Thickness
        rlobecths(4,1) = std(right_cth(ind)); % Std Cortical Thickness
        rlobecurv(4,1) = mean(right_curv(ind)); % Mean curvature
        rlobevols(4,1) = sum(c(OccIdsR(:))*prod(voxsize)); % Volume
        % ======================End of Right Hemisphere ==========================%
        
        
        % Computing Hull Surface
        [Surfl,Surfr] = Create_Region_ConvexHull(VA.fname);
        Surfa(1) = Surfl;
        Surfa(2) = Surfr;
        try
            save(outhsurf,'Surfa');
        end
        
        % Hemispheres
        indhsl = ismember(Surfl.StructS(:,1),[FroIds(:);ParIds(:);TempIds(:);OccIds(:);InsIds(:)]);
        lhhulls = sum(Surfl.StructS(indhsl,2));
        indhsr = ismember(Surfr.StructS(:,1),[FroIdsR(:);ParIdsR(:);TempIdsR(:);OccIdsR(:);InsIdsR(:)]);
        rhhulls = sum(Surfr.StructS(indhsr,2));
        hullsurf = rhhulls + lhhulls;
        
        % ------ Left Lobes ----- Frontal
        indhtemp = ismember(Surfl.StructS(:,1),[FroIds(:)]);
        llobehulls(1,1) = sum(Surfl.StructS(indhtemp,2));
        % ------ Left Lobes ----- Parietal
        indhtemp = ismember(Surfl.StructS(:,1),[ParIds(:)]);
        llobehulls(2,1) = sum(Surfl.StructS(indhtemp,2));
        % ------ Left Lobes ----- Temporal
        indhtemp = ismember(Surfl.StructS(:,1),[TempIds(:)]);
        llobehulls(3,1) = sum(Surfl.StructS(indhtemp,2));
        % ------ Left Lobes ----- Occipital
        indhtemp = ismember(Surfl.StructS(:,1),[OccIds(:)]);
        llobehulls(4,1) = sum(Surfl.StructS(indhtemp,2));
        
        % ------ Right Lobes ----- Frontal
        indhtemp = ismember(Surfr.StructS(:,1),[FroIdsR(:)]);
        rlobehulls(1,1) = sum(Surfr.StructS(indhtemp,2));
        % ------ Right Lobes ----- Parietal
        indhtemp = ismember(Surfr.StructS(:,1),[ParIdsR(:)]);
        rlobehulls(2,1) = sum(Surfr.StructS(indhtemp,2));
        % ------ Right Lobes ----- Temporal
        indhtemp = ismember(Surfr.StructS(:,1),[TempIdsR(:)]);
        rlobehulls(3,1) = sum(Surfr.StructS(indhtemp,2));
        % ------ Right Lobes ----- Occipital
        indhtemp = ismember(Surfr.StructS(:,1),[OccIdsR(:)]);
        rlobehulls(4,1) = sum(Surfr.StructS(indhtemp,2));
        
        % Left ROIs
        [a,b] = ismember(ctabl.table(:,6),Surfl.StructS(:,1));
        lhhsurf = zeros(size(ctabl.table(:,6),1),1);
        lhhsurf(a) = Surfl.StructS(nonzeros(b),2);
        
        % Right ROIs
        [a,b] = ismember(ctabr.table(:,6),Surfr.StructS(:,1));
        rhhsurf = zeros(size(ctabr.table(:,6),1),1);
        rhhsurf(a) = Surfr.StructS(nonzeros(b),2);
        
              
        
        LIds = [FroIds(:);ParIds(:);TempIds(:);OccIds(:);InsIds(:)];
        indl = find(ismember(ctabl.table(:,6),LIds(:)));
        Sids = ctabl.table(indl,5);
        indl = find(ismember(txtl,Sids));
        
        RIds = [FroIdsR(:);ParIdsR(:);TempIdsR(:);OccIdsR(:);InsIdsR(:)];
        indr = find(ismember(ctabr.table(:,6),RIds(:)));
        Sids = ctabr.table(indr,5);
        indr = find(ismember(txtr,Sids));
        
        cortnvert = lhnvert + rhnvert;
        cortcthm = mean([left_cth(indl);right_cth(indr)]); % Mean Cortical Thickness
        cortcths = std([left_cth(indl);right_cth(indr)]); % Std Cortical Thickness
        cortcurv = mean([left_curv(indl);right_curv(indr)]); % Cortical Mean curvature
        cortarea = rharea+lharea; % Cortical Area
        
        BIds = [FroIds(:);ParIds(:);TempIds(:);OccIds(:);InsIds(:);FroIdsR(:);ParIdsR(:);TempIdsR(:);OccIdsR(:);InsIdsR(:)];
        cortvol = sum(c(BIds)*prod(voxsize)); % Cortical Volume
        
        [Codes,Names] = Brain_GM_codes('aparc+aseg');
        ind = find(ismember(IA,[Codes(:);2;41]));
        TBV = length(ind).*prod(voxsize); % Total Brain Volume
        ind = find(IA);
        ICV = length(ind).*prod(voxsize); % Intracranial Volume
        
        % Number of vertex
        NvertNames = {'Total_NumVert';'LH_NumVert';'RH_NumVert'};
        Temp = [cellstr([repmat('LH_NumVert-',[4 1]) strvcat('Frontal','Parietal','Temporal','Occipital')]) cellstr([repmat('RH_NumVert-',[4 1]) strvcat('Frontal','Parietal','Temporal','Occipital')])]';
        Temp= Temp(:);
        Temp1 = [cellstr([repmat('LH_NumVert-',[size(ctabl.struct_names,1) 1]) char(ctabl.struct_names)]) cellstr([repmat('RH_NumVert-',[size(ctabr.struct_names,1) 1]) char(ctabr.struct_names)])]';
        Temp1= Temp1(:);
        NvertNames = [NvertNames;Temp;Temp1];
        lobesnvert = [llobenvert rlobenvert]';lobesnvert = lobesnvert(:);
        roisnvert = [lnvert rnvert]';roisnvert = roisnvert(:);
        Numvert = [cortnvert;lhnvert;rhnvert;lobesnvert;roisnvert];
        
        % Cortical Volume
        VCortNames = {'Total_CortVol';'LH_CortVol';'RH_CortVol'};
        Temp = [cellstr([repmat('LH_CortVol-',[4 1]) strvcat('Frontal','Parietal','Temporal','Occipital')]) cellstr([repmat('RH_CortVol-',[4 1]) strvcat('Frontal','Parietal','Temporal','Occipital')])]';
        Temp= Temp(:);
        Temp1 = [cellstr([repmat('LH_CortVol-',[size(ctabl.struct_names,1) 1]) char(ctabl.struct_names)]) cellstr([repmat('RH_CortVol-',[size(ctabr.struct_names,1) 1]) char(ctabr.struct_names)])]';
        Temp1= Temp1(:);
        VCortNames = [VCortNames;Temp;Temp1];
        lobesvcort = [llobevols rlobevols]';lobesvcort = lobesvcort(:);
        roisvcort = [lvols rvols]';roisvcort = roisvcort(:);
        CortVol = [cortvol;lhvols;rhvols;lobesvcort;roisvcort];
        
        % Cortical Area
        ACortNames = {'Total_CortArea';'LH_CortArea';'RH_CortArea'};
        Temp = [cellstr([repmat('LH_CortArea-',[4 1]) strvcat('Frontal','Parietal','Temporal','Occipital')]) cellstr([repmat('RH_CortArea-',[4 1]) strvcat('Frontal','Parietal','Temporal','Occipital')])]';
        Temp= Temp(:);
        Temp1 = [cellstr([repmat('LH_CortArea-',[size(ctabl.struct_names,1) 1]) char(ctabl.struct_names)]) cellstr([repmat('RH_CortArea-',[size(ctabr.struct_names,1) 1]) char(ctabr.struct_names)])]';
        Temp1= Temp1(:);
        ACortNames = [ACortNames;Temp;Temp1];
        lobesarea = [llobearea rlobearea]';lobesarea = lobesarea(:);
        roisacort = [larea rarea]';roisacort = roisacort(:);
        CortArea = [cortarea;lharea;rharea;lobesarea;roisacort];
        
        % Mean Cortical Thickness
        CortThickNames = {'Total_MeanCortThick';'LH_MeanCortThick';'RH_MeanCortThick'};
        Temp = [cellstr([repmat('LH_MeanCortThick-',[4 1]) strvcat('Frontal','Parietal','Temporal','Occipital')]) cellstr([repmat('RH_MeanCortThick-',[4 1]) strvcat('Frontal','Parietal','Temporal','Occipital')])]';
        Temp= Temp(:);
        Temp1 = [cellstr([repmat('LH_MeanCortThick-',[size(ctabl.struct_names,1) 1]) char(ctabl.struct_names)]) cellstr([repmat('RH_MeanCortThick-',[size(ctabr.struct_names,1) 1]) char(ctabr.struct_names)])]';
        Temp1= Temp1(:);
        CortThickNames = [CortThickNames;Temp;Temp1];
        lobescortt = [llobecthm rlobecthm]';lobescortt = lobescortt(:);
        roiscortt = [lcthm rcthm]';roiscortt = roiscortt(:);
        CortCthM = [cortcthm;lhcthm;rhcthm;lobescortt;roiscortt];
        
        % STD Cortical Thickness
        CortThickSNames = {'Total_StdCortThick';'LH_StdCortThick';'RH_StdCortThick'};
        Temp = [cellstr([repmat('LH_StdCortThick-',[4 1]) strvcat('Frontal','Parietal','Temporal','Occipital')]) cellstr([repmat('RH_StdCortThick-',[4 1]) strvcat('Frontal','Parietal','Temporal','Occipital')])]';
        Temp= Temp(:);
        Temp1 = [cellstr([repmat('LH_StdCortThick-',[size(ctabl.struct_names,1) 1]) char(ctabl.struct_names)]) cellstr([repmat('RH_StdCortThick-',[size(ctabr.struct_names,1) 1]) char(ctabr.struct_names)])]';
        Temp1= Temp1(:);
        CortThickSNames = [CortThickSNames;Temp;Temp1];
        lobescortts = [llobecths rlobecths]';lobescortts = lobescortts(:);
        roiscortts = [lcths rcths]';roiscortts = roiscortts(:);
        CortCthS = [cortcths;lhcths;rhcths;lobescortts;roiscortts];
        
        % Mean Curvature
        CurvNames = {'Total_MeanCurv';'LH_MeanCurv';'RH_MeanCurv'};
        Temp = [cellstr([repmat('LH_MeanCurv-',[4 1]) strvcat('Frontal','Parietal','Temporal','Occipital')]) cellstr([repmat('RH_MeanCurv-',[4 1]) strvcat('Frontal','Parietal','Temporal','Occipital')])]';
        Temp= Temp(:);
        Temp1 = [cellstr([repmat('LH_MeanCurv-',[size(ctabl.struct_names,1) 1]) char(ctabl.struct_names)]) cellstr([repmat('RH_MeanCurv-',[size(ctabr.struct_names,1) 1]) char(ctabr.struct_names)])]';
        Temp1= Temp1(:);
        CurvNames = [CurvNames;Temp;Temp1];
        lobescurvs = [llobecurv rlobecurv]';lobescurvs = lobescurvs(:);
        roiscurvs = [lcurv rcurv]';roiscurvs = roiscurvs(:);
        CortCurv = [cortcurv;lhcurv;rhcurv;lobescurvs;roiscurvs];
        
        % Hull Surface
        HullNames = {'Total_HullArea';'LH_HullArea';'RH_HullArea'};
        Temp = [cellstr([repmat('LH_HullArea-',[4 1]) strvcat('Frontal','Parietal','Temporal','Occipital')]) cellstr([repmat('RH_HullArea-',[4 1]) strvcat('Frontal','Parietal','Temporal','Occipital')])]';
        Temp= Temp(:);
        Temp1 = [cellstr([repmat('LH_HullArea-',[size(ctabl.struct_names,1) 1]) char(ctabl.struct_names)]) cellstr([repmat('RH_HullArea-',[size(ctabr.struct_names,1) 1]) char(ctabr.struct_names)])]';
        Temp1= Temp1(:);
        HullNames = [HullNames;Temp;Temp1];
        lobeshulls = [llobehulls rlobehulls]';lobeshulls = lobeshulls(:);
        roishulls = [lhhsurf rhhsurf]';roishulls = roishulls(:);
        HullSurf = [hullsurf;lhhulls;rhhulls;lobeshulls;roishulls];
        
        % Regional Gyrification Index
        RGyrifNames = {'Total_RegGyrifIndex';'LH_RegGyrifIndex';'RH_RegGyrifIndex'};
        Temp = [cellstr([repmat('LH_RegGyrifIndex-',[4 1]) strvcat('Frontal','Parietal','Temporal','Occipital')]) cellstr([repmat('RH_RegGyrifIndex-',[4 1]) strvcat('Frontal','Parietal','Temporal','Occipital')])]';
        Temp= Temp(:);
        Temp1 = [cellstr([repmat('LH_RegGyrifIndex-',[size(ctabl.struct_names,1) 1]) char(ctabl.struct_names)]) cellstr([repmat('RH_RegGyrifIndex-',[size(ctabr.struct_names,1) 1]) char(ctabr.struct_names)])]';
        Temp1= Temp1(:);
        RGyrifNames = [RGyrifNames;Temp;Temp1];
        indh = find(HullSurf);
        inds = find(CortArea);
        indt = find(ismember(indh,inds));
        RGyrifIndex = HullSurf*0;
        RGyrifIndex(indh(indt)) = CortArea(indh(indt))./HullSurf(indh(indt));
        
        for k=1:size(CurvNames,1);
            temp = [deblank(char(VCortNames(k))) '(mm^3)'];
            VCortNames{k} = temp;
            temp = [deblank(char(ACortNames(k))) '(mm^2)'];
            ACortNames{k} = temp;
            temp = [deblank(char(CortThickNames(k))) '(mm)'];
            CortThickNames{k} = temp;
            temp = [deblank(char(CortThickSNames(k))) '(mm)'];
            CortThickSNames{k} = temp;
            temp = [deblank(char(CurvNames(k))) ];
            CurvNames{k} = temp;
        end
        
        
        cadvarname = 'Subject_ID;Intracranial_Volume(mm^3);Total_Brain_Volume(mm^3)';
        cadvars = [subjId ';' num2str(ICV) ';' num2str(TBV)];
        for j = 1:length(CortCurv)
            cadvarname = [cadvarname ';' char(NvertNames(j)) ';' char(VCortNames(j)) ';' char(ACortNames(j)) ';' char(CortThickNames(j)) ';' char(CortThickSNames(j)) ';' char(CurvNames(j)) ';' char(HullNames(j))  ';' char(RGyrifNames(j))];
            cadvars = [cadvars ';' num2str(Numvert(j)) ';' num2str(CortVol(j)) ';' num2str(CortArea(j)) ';' num2str(CortCthM(j)) ';' num2str(CortCthS(j)) ';' num2str(CortCurv(j))  ';' num2str(HullSurf(j)) ';' num2str(RGyrifIndex(j))];
        end
        if z == 1
            CadTotal = strvcat(cadvarname,cadvars);
        else
            CadTotal = strvcat(CadTotal,cadvars);
        end
        OStatFile = [FreeSDir filesep subjId filesep 'stats' filesep subjId '_whole_freesurfer_stats_' aparcId '-NoLGI.txt'];
        fid = fopen(OStatFile,'wt');
        fprintf(fid,'%s\n',cadvarname);
        fprintf(fid,'%s\n',cadvars);
        fclose(fid);
    end
    
    % ------------------- Saving Stat File ------------------------------------
    if nargin <3
        OutStatFile = CadTotal;
    else
        fid = fopen(OutStatFile,'wt');
        for z = 1:size(CadTotal,1)
            fprintf(fid,'%s\n',CadTotal(z,:));
        end
        fclose(fid);
    end
end
return;

function At = Comp_Surf_Area_custom(fv);
%=========================Main program====================================%
d12 = sqrt(sum((fv.vertices(fv.faces(:,1),:) - fv.vertices(fv.faces(:,2),:)).^2,2));
d23 = sqrt(sum((fv.vertices(fv.faces(:,2),:) - fv.vertices(fv.faces(:,3),:)).^2,2));
d13 = sqrt(sum((fv.vertices(fv.faces(:,1),:) - fv.vertices(fv.faces(:,3),:)).^2,2));
per = (d12+d23+d13)/2;
Af = sqrt(per.*(per-d12).*(per-d23).*(per-d13))/100;% cm^2
At = sum(Af); % Area per each face
%========================End of main program==============================%
return

function [curv, fnum] = read_char(fname);

fid = fopen(fname, 'rb', 'b') ;
if (fid < 0)
    str = sprintf('could not open file %s.', fname) ;
    error(str) ;
end
% vnum = fread3(fid) ;
b1 = fread(fid, 1, 'uchar') ;
b2 = fread(fid, 1, 'uchar') ;
b3 = fread(fid, 1, 'uchar') ;
vnum = bitshift(b1, 16) + bitshift(b2,8) + b3 ;

NEW_VERSION_MAGIC_NUMBER = 16777215;
if (vnum == NEW_VERSION_MAGIC_NUMBER)
    vnum = fread(fid, 1, 'int32') ;
    fnum = fread(fid, 1, 'int32') ;
    vals_per_vertex = fread(fid, 1, 'int32') ;
    curv = fread(fid, vnum, 'float') ;
    
    fclose(fid) ;
else
    b1 = fread(fid, 1, 'uchar') ;
    b2 = fread(fid, 1, 'uchar') ;
    b3 = fread(fid, 1, 'uchar') ;
    vnum = bitshift(b1, 16) + bitshift(b2,8) + b3 ;
    curv = fread(fid, vnum, 'int16') ./ 100 ;
    fclose(fid) ;
end

function [retval] = fread3(fid)

b1 = fread(fid, 1, 'uchar') ;
b2 = fread(fid, 1, 'uchar') ;
b3 = fread(fid, 1, 'uchar') ;
retval = bitshift(b1, 16) + bitshift(b2,8) + b3 ;

function [vol, M, mr_parms, volsz] = load_mgh(fname,slices,frames,headeronly)
% [vol, M, mr_parms, Mdc, volsz] = load_mgh(fname,<slices>,<frames>,<headeronly>)
%
% fname - path of the mgh file
%
% slices - list of one-based slice numbers to load. All
%   slices are loaded if slices is not specified, or
%   if slices is empty, or if slices(1) <= 0.
%
% frames - list of one-based frame numbers to load. All
%   frames are loaded if frames is not specified, or
%   if frames is empty, or if frames(1) <= 0.
%
% M is the 4x4 vox2ras transform such that
% y(i1,i2,i3), xyz1 = M*[i1 i2 i3 1] where the
% indices are 0-based. If the input has multiple frames,
% only the first frame is read.
%
% mr_parms = [tr flipangle te ti fov]
%
% volsz = size(vol). Helpful when using headeronly as vol is [].
%
% See also: save_mgh, vox2ras_0to1
%


%
% load_mgh.m
%
% Original Author: Bruce Fischl
% CVS Revision Info:
%    $Author: greve $
%    $Date: 2009/07/01 17:13:08 $
%    $Revision: 1.16.2.2 $
%
% Copyright (C) 2002-2007,
% The General Hospital Corporation (Boston, MA).
% All rights reserved.
%
% Distribution, usage and copying of this software is covered under the
% terms found in the License Agreement file named 'COPYING' found in the
% FreeSurfer source code root directory, and duplicated here:
% https://surfer.nmr.mgh.harvard.edu/fswiki/FreeSurferOpenSourceLicense
%
% General inquiries: freesurfer@nmr.mgh.harvard.edu
% Bug reports: analysis-bugs@nmr.mgh.harvard.edu
%

vol = [];
M = [];
mr_parms = [];
volsz = [];

if(nargin < 1 | nargin > 4)
    msg = 'USAGE: [vol M] = load_mgh(fname,<slices>,<frames>,<headeronly>)';
    fprintf('%s',msg);
    return;
end

% unzip if it is compressed
if (strcmpi(fname((strlen(fname)-3):strlen(fname)), '.MGZ') | ...
        strcmpi(fname((strlen(fname)-3):strlen(fname)), '.GZ'))
    rand('state', sum(100*clock));
    gzipped =  round(rand(1)*10000000 + ...
        sum(int16(fname))) + round(cputime);
    ind = findstr(fname, '.');
    new_fname = sprintf('/tmp/tmp%d.mgh', gzipped);
    if(strcmp(computer,'MAC') || strcmp(computer,'MACI') || ismac)
        unix(sprintf('gunzip -c %s > %s', fname, new_fname)) ;
    else
        unix(sprintf('zcat %s > %s', fname, new_fname)) ;
    end
    fname = new_fname ;
else
    gzipped = -1 ;
end


if(exist('slices')~=1) slices = []; end
if(isempty(slices)) slices = 0; end
if(slices(1) <= 0) slices = 0; end

if(exist('frames')~=1) frames = []; end
if(isempty(frames)) frames = 0; end
if(frames(1) <= 0) frames = 0; end

if(exist('headeronly')~=1) headeronly = 0; end

fid    = fopen(fname, 'rb', 'b') ;
if(fid == -1)
    fprintf('ERROR: could not open %s for reading\n',fname);
    return;
end
v       = fread(fid, 1, 'int') ;
if(isempty(v))
    fprintf('ERROR: problem reading fname\n');
    if(gzipped >=0) unix(sprintf('rm %s', fname)); end
end
ndim1   = fread(fid, 1, 'int') ;
ndim2   = fread(fid, 1, 'int') ;
ndim3   = fread(fid, 1, 'int') ;
nframes = fread(fid, 1, 'int') ;
type    = fread(fid, 1, 'int') ;
dof     = fread(fid, 1, 'int') ;

if(slices(1) > 0)
    ind = find(slices > ndim3);
    if(~isempty(ind))
        fprintf('ERROR: load_mgh: some slices exceed nslices\n');
        return;
    end
end

if(frames(1) > 0)
    ind = find(frames > nframes);
    if(~isempty(ind))
        fprintf('ERROR: load_mgh: some frames exceed nframes\n');
        return;
    end
end

UNUSED_SPACE_SIZE= 256;
USED_SPACE_SIZE = (3*4+4*3*4);  % space for ras transform

unused_space_size = UNUSED_SPACE_SIZE-2 ;
ras_good_flag = fread(fid, 1, 'short') ;
if (ras_good_flag)
    delta  = fread(fid, 3, 'float32') ;
    Mdc    = fread(fid, 9, 'float32') ;
    Mdc    = reshape(Mdc,[3 3]);
    Pxyz_c = fread(fid, 3, 'float32') ;
    
    D = diag(delta);
    
    Pcrs_c = [ndim1/2 ndim2/2 ndim3/2]'; % Should this be kept?
    
    Pxyz_0 = Pxyz_c - Mdc*D*Pcrs_c;
    
    M = [Mdc*D Pxyz_0;  ...
        0 0 0 1];
    ras_xform = [Mdc Pxyz_c; ...
        0 0 0 1];
    unused_space_size = unused_space_size - USED_SPACE_SIZE ;
end

fseek(fid, unused_space_size, 'cof') ;
nv = ndim1 * ndim2 * ndim3 * nframes;
volsz = [ndim1 ndim2 ndim3 nframes];

MRI_UCHAR =  0 ;
MRI_INT =    1 ;
MRI_LONG =   2 ;
MRI_FLOAT =  3 ;
MRI_SHORT =  4 ;
MRI_BITMAP = 5 ;

% Determine number of bytes per voxel
switch type
    case MRI_FLOAT,
        nbytespervox = 4;
    case MRI_UCHAR,
        nbytespervox = 1;
    case MRI_SHORT,
        nbytespervox = 2;
    case MRI_INT,
        nbytespervox = 4;
end

if(headeronly)
    fseek(fid,nv*nbytespervox,'cof');
    if(~feof(fid))
        [mr_parms count] = fread(fid,4,'float32');
        if(count ~= 4)
            fprintf('WARNING: error reading MR params\n');
        end
    end
    fclose(fid);
    if(gzipped >=0)  unix(sprintf('rm %s', fname));  end
    return;
end


%------------------ Read in the entire volume ----------------%
if(slices(1) <= 0 & frames(1) <= 0)
    switch type
        case MRI_FLOAT,
            vol = fread(fid, nv, 'float32') ;
        case MRI_UCHAR,
            vol = fread(fid, nv, 'uchar') ;
        case MRI_SHORT,
            vol = fread(fid, nv, 'short') ;
        case MRI_INT,
            vol = fread(fid, nv, 'int') ;
    end
    
    if(~feof(fid))
        [mr_parms count] = fread(fid,4,'float32');
        if(count ~= 4)
            fprintf('WARNING: error reading MR params\n');
        end
    end
    fclose(fid) ;
    if(gzipped >=0)  unix(sprintf('rm %s', fname));  end
    
    nread = prod(size(vol));
    if(nread ~= nv)
        fprintf('ERROR: tried to read %d, actually read %d\n',nv,nread);
        vol = [];
        return;
    end
    vol = reshape(vol,[ndim1 ndim2 ndim3 nframes]);
    
    return;
end

%----- only gets here if a subest of slices/frames are to be loaded ---------%


if(frames(1) <= 0) frames = [1:nframes]; end
if(slices(1) <= 0) slices = [1:ndim3]; end

nvslice = ndim1 * ndim2;
nvvol   = ndim1 * ndim2 * ndim3;
filepos0 = ftell(fid);
vol = zeros(ndim1,ndim2,length(slices),length(frames));
nthframe = 1;
for frame = frames
    
    nthslice = 1;
    for slice = slices
        filepos = ((frame-1)*nvvol + (slice-1)*nvslice)*nbytespervox + filepos0;
        fseek(fid,filepos,'bof');
        
        switch type
            case MRI_FLOAT,
                [tmpslice nread]  = fread(fid, nvslice, 'float32') ;
            case MRI_UCHAR,
                [tmpslice nread]  = fread(fid, nvslice, 'uchar') ;
            case MRI_SHORT,
                [tmpslice nread]  = fread(fid, nvslice, 'short') ;
            case MRI_INT,
                [tmpslice nread]  = fread(fid, nvslice, 'int') ;
        end
        
        if(nread ~= nvslice)
            fprintf('ERROR: load_mgh: reading slice %d, frame %d\n',slice,frame);
            fprintf('  tried to read %d, actually read %d\n',nvslice,nread);
            fclose(fid);
            if(gzipped >=0) unix(sprintf('rm %s', fname)); end
            return;
        end
        
        vol(:,:,nthslice,nthframe) = reshape(tmpslice,[ndim1 ndim2]);
        nthslice = nthslice + 1;
    end
    
    nthframe = nthframe + 1;
end

% seek to just beyond the last slice/frame %
filepos = (nframes*nvvol)*nbytespervox + filepos0;
fseek(fid,filepos,'bof');

if(~feof(fid))
    [mr_parms count] = fread(fid,5,'float32');
    if(count < 4)
        fprintf('WARNING: error reading MR params\n');
    end
end

fclose(fid) ;
if(gzipped >=0) unix(sprintf('rm %s', fname)); end

return;

function [vertices, label, colortable] = Read_Brain_Annotation(filename)
% [vertices, label, colortable] = Read_Brain_Annotation(annotfilename.annot)
%
% vertices expected to be simply from 0 to number of vertices - 1;
% label is the vector of annotation
%
% colortable is empty struct if not embedded in .annot. Else, it will be
% a struct.
% colortable.numEntries = number of Entries
% colortable.orig_tab = name of original colortable
% colortable.struct_names = list of structure names (e.g. central sulcus and so on)
% colortable.table = n x 5 matrix. 1st column is r, 2nd column is g, 3rd column
% is b, 4th column is flag, 5th column is resultant integer values
% calculated from r + g*2^8 + b*2^16 + flag*2^24. flag expected to be all 0.


%
% read_annotation.m
%
% Original Author: Bruce Fischl
% CVS Revision Info:
%    $Author: nicks $
%    $Date: 2007/01/10 22:55:09 $
%    $Revision: 1.4 $
%
% Copyright (C) 2002-2007,
% The General Hospital Corporation (Boston, MA).
% All rights reserved.
%
% Distribution, usage and copying of this software is covered under the
% terms found in the License Agreement file named 'COPYING' found in the
% FreeSurfer source code root directory, and duplicated here:
% https://surfer.nmr.mgh.harvard.edu/fswiki/FreeSurferOpenSourceLicense
%
% General inquiries: freesurfer@nmr.mgh.harvard.edu
% Bug reports: analysis-bugs@nmr.mgh.harvard.edu
%

fp = fopen(filename, 'r', 'b');

if(fp < 0)
    disp('Annotation file cannot be opened');
    return;
end

A = fread(fp, 1, 'int');

tmp = fread(fp, 2*A, 'int');
vertices = tmp(1:2:end);
label = tmp(2:2:end);

bool = fread(fp, 1, 'int');
if(isempty(bool)) %means no colortable
    disp('No Colortable found.');
    colortable = struct([]);
    fclose(fp);
    return;
end

if(bool)
    
    %Read colortable
    numEntries = fread(fp, 1, 'int');
    
    if(numEntries > 0)
        
        disp(['Reading from Original Version']);
        colortable.numEntries = numEntries;
        len = fread(fp, 1, 'int');
        colortable.orig_tab = fread(fp, len, '*char')';
        colortable.orig_tab = colortable.orig_tab(1:end-1);
        
        colortable.struct_names = cell(numEntries,1);
        colortable.table = zeros(numEntries,5);
        for i = 1:numEntries
            len = fread(fp, 1, 'int');
            colortable.struct_names{i} = fread(fp, len, '*char')';
            colortable.struct_names{i} = colortable.struct_names{i}(1:end-1);
            colortable.table(i,1) = fread(fp, 1, 'int');
            colortable.table(i,2) = fread(fp, 1, 'int');
            colortable.table(i,3) = fread(fp, 1, 'int');
            colortable.table(i,4) = fread(fp, 1, 'int');
            colortable.table(i,5) = colortable.table(i,1) + colortable.table(i,2)*2^8 + colortable.table(i,3)*2^16 + colortable.table(i,4)*2^24;
        end
        disp(['colortable with ' num2str(colortable.numEntries) ' entries read (originally ' colortable.orig_tab ')']);
        
    else
        version = -numEntries;
        if(version~=2)
            disp(['Error! Does not handle version ' num2str(version)]);
        else
            disp(['Reading from version ' num2str(version)]);
        end
        numEntries = fread(fp, 1, 'int');
        colortable.numEntries = numEntries;
        len = fread(fp, 1, 'int');
        colortable.orig_tab = fread(fp, len, '*char')';
        colortable.orig_tab = colortable.orig_tab(1:end-1);
        
        colortable.struct_names = cell(numEntries,1);
        colortable.table = zeros(numEntries,5);
        
        numEntriesToRead = fread(fp, 1, 'int');
        for i = 1:numEntriesToRead
            structure = fread(fp, 1, 'int')+1;
            if (structure < 0)
                disp(['Error! Read entry, index ' num2str(structure)]);
            end
            if(~isempty(colortable.struct_names{structure}))
                disp(['Error! Duplicate Structure ' num2str(structure)]);
            end
            len = fread(fp, 1, 'int');
            colortable.struct_names{structure} = fread(fp, len, '*char')';
            colortable.struct_names{structure} = colortable.struct_names{structure}(1:end-1);
            colortable.table(structure,1) = fread(fp, 1, 'int');
            colortable.table(structure,2) = fread(fp, 1, 'int');
            colortable.table(structure,3) = fread(fp, 1, 'int');
            colortable.table(structure,4) = fread(fp, 1, 'int');
            colortable.table(structure,5) = colortable.table(structure,1) + colortable.table(structure,2)*2^8 + colortable.table(structure,3)*2^16 + colortable.table(structure,4)*2^24;
        end
        disp(['colortable with ' num2str(colortable.numEntries) ' entries read (originally ' colortable.orig_tab ')']);
    end
else
    disp('Error! Should not be expecting bool = 0');
end

fclose(fp);
return;


