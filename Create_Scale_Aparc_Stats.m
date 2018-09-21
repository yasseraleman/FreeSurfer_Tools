function OutStats = Create_Scale_Aparc_Stats(FreeSDir, subjId, aparcId, thtype);
%
% Syntax :
% Create_Scale_Aparc(FreeSDir, subjId, aparcId);
%
% This function creates new stat files for left and right hemisphere
% according to a new cortical parcellation.
%
% Input Parameters:
%      FreeSDir         : FreeSurfer Directory
%      subjId           : Subject Id 
%      aparcId          : Parcellation identification
%
% Output Parameters:
%      OutStats         : New Stat Files
%
% Related references:
%
%
% See also:
%__________________________________________________
% Authors: Yasser Aleman Gomez
% LIM
% September 12th 2012
% Version $1.0

tic
% aparcId = 'aparc_scale01_Nzones0097';
% FreeSDir = '/media/Data/PROCESSING_RESULTS/PEPS/5-freesurfer_processing';
% subjId = '0002-20090515';
if nargin <4
    thtype = '.thickness';
end
switch thtype
    case '.thickness'
        Ident = 'cth_freesurfer';  % Cortical thickness computed via freesurfer
    case '.segload.thickness'
        Ident = 'cth_direct_segload';  % Cortical thickness computed via Direct and segmentations obtained from segload
    case '.fslfast.thickness'
        Ident = 'cth_direct_fast'; % Cortical thickness computed via Direct and segmentations obtained from fast
    case '.spmvbm8.thickness'
        Ident = 'cth_direct_vbm8'; % Cortical thickness computed via Direct and segmentations obtained from fast
end

%% ================= Reading Necessary files ============================ %
% Statistic Files
lhaparc = [FreeSDir filesep subjId filesep 'stats' filesep 'lh.aparc.stats'];
rhaparc = [FreeSDir filesep subjId filesep 'stats' filesep 'rh.aparc.stats'];
NLaparc = [FreeSDir filesep subjId filesep 'stats' filesep 'lh.' Ident '.' aparcId '.stats'];
NRaparc = [FreeSDir filesep subjId filesep 'stats' filesep 'rh.' Ident '.' aparcId '.stats'];

% Parcellation Maps
lhannot = [FreeSDir filesep subjId filesep 'label' filesep 'lh.' aparcId '.annot'];
rhannot = [FreeSDir filesep subjId filesep 'label' filesep 'rh.' aparcId '.annot'];

% Surfaces
lhsurf = [FreeSDir filesep subjId filesep 'surf' filesep 'lh.white'];
rhsurf = [FreeSDir filesep subjId filesep 'surf' filesep 'rh.white'];

% Curvature Maps
lhcurv = [FreeSDir filesep subjId filesep 'surf' filesep 'lh.curv'];
rhcurv = [FreeSDir filesep subjId filesep 'surf' filesep 'rh.curv'];

% Thickness Maps
lhcth = [FreeSDir filesep subjId filesep 'surf' filesep 'lh' thtype];
rhcth = [FreeSDir filesep subjId filesep 'surf' filesep 'rh' thtype];

% Volumetric Atlas
atlasfile = [FreeSDir filesep subjId filesep 'mri' filesep aparcId '+aseg.mgz'];

%% ================= End of Reading Necessary files ===================== %



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
[left_curv,ctab] = read_cfiles(lhcurv);
[right_curv,ctab] = read_cfiles(rhcurv);

%----- Reading Thickness Map
[left_cth,ctab] = read_cfiles(lhcth);
[right_cth,ctab] = read_cfiles(rhcth);

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
Nr = size(ctabl.table,1);
fv.vertices = Surfwl.SurfData.vertices;
for j = 1:Nr
    ind = find(txtl == ctabl.table(j,5));
    nvert(j,1) = length(ind); % Number of vertices
    lcthm(j,1) = mean(left_cth(ind)); % Mean Cortical Thickness
    lcths(j,1) = std(left_cth(ind)); % Std Cortical Thickness
    lcurv(j,1) = mean(left_curv(ind)); % Mean curvature
    
    %----- Computing Area
    indf = find(sum(ismember( Surfwl.SurfData.faces,ind)')==3);
    At = 0;
    fv.faces = Surfwl.SurfData.faces(indf,:);
    N = size(fv.faces,1);
    for i = 1:N;
        di = dist(fv.vertices(fv.faces(i,:),:)');
        A = abs(di(1,2));
        B = abs(di(1,3));
        C = abs(di(2,3));
        p = (A+B+C)/2;
        Ar = sqrt(p*(p-A)*(p-B)*(p-C));
        %Ar = (A*B/2)*(sqrt(1-((A^2+B^2-C^2)^2)/((2*A*B)^2)));
        At = At+Ar;
    end
    larea(j,1) = real(At);% cm^2
end
%----- Computing Volume
lvols = c(ctabl.table(:,6))*prod(voxsize);

% ---- Reading old aparc
fid = fopen(lhaparc);
cont = 0;line = '#';
lines = '';
while strcmp(line(1),'#')
    cont = cont + 1;
    line = fgetl(fid);
    if ~ischar(line),   break,   end
    lines = strvcat(lines,line);
end
lines = lines(1:end-1,:);
fclose(fid);

% ---- Creating New aparc
ids = char(ctabl.struct_names);
Total = [ids  repmat('               ',[size(ids,1) 1]) num2str(nvert) repmat('  ',[size(ids,1) 1]) ...
        num2str(larea) repmat('  ',[size(ids,1) 1]) num2str(lvols) repmat('  ',[size(ids,1) 1])...
        num2str(lcthm) repmat('  ',[size(ids,1) 1]) num2str(lcths) repmat('  ',[size(ids,1) 1])...
        num2str(abs(lcurv)) repmat('  ',[size(ids,1) 1]) num2str(zeros(size(ids,1),1)) repmat('  ',[size(ids,1) 1])...
        num2str(zeros(size(ids,1),1)) repmat('  ',[size(ids,1) 1]) num2str(zeros(size(ids,1),1))];
lines = strvcat(lines,Total);

% ==================== Saving the New Color File ======================= %
fid = fopen(NLaparc,'wt');
Ns = size(lines,1);
for i = 1:Ns
    fprintf(fid, '%s\n',deblank(lines(i,:)));
end
fclose(fid);
% ============== End of Saving the New Color File ====================== %
% ======================End of Left Hemisphere ==========================%

% ======================== Right Hemisphere ===============================%
Nr = size(ctabr.table,1);
fv.vertices = Surfwr.SurfData.vertices;
for j = 1:Nr
    ind = find(txtr == ctabr.table(j,5));
    nvert(j,1) = length(ind); % Number of vertices
    rcthm(j,1) = mean(right_cth(ind)); % Mean Cortical Thickness
    rcths(j,1) = std(right_cth(ind)); % Std Cortical Thickness
    rcurv(j,1) = mean(right_curv(ind)); % Mean curvature
    
    %----- Computing Area
    indf = find(sum(ismember( Surfwr.SurfData.faces,ind)')==3);
    At = 0;
    fv.faces = Surfwr.SurfData.faces(indf,:);
    N = size(fv.faces,1);
    for i = 1:N;
        di = dist(fv.vertices(fv.faces(i,:),:)');
        A = abs(di(1,2));
        B = abs(di(1,3));
        C = abs(di(2,3));
        p = (A+B+C)/2;
        Ar = sqrt(p*(p-A)*(p-B)*(p-C));
        %Ar = (A*B/2)*(sqrt(1-((A^2+B^2-C^2)^2)/((2*A*B)^2)));
        At = At+Ar;
    end
    rarea(j,1) = real(At);% cm^2
end
%----- Computing Volume
rvols = c(ctabr.table(:,6))*prod(voxsize);

% ---- Reading old aparc
fid = fopen(rhaparc);
cont = 0;line = '#';
lines = '';
while strcmp(line(1),'#')
    cont = cont + 1;
    line = fgetl(fid);
    if ~ischar(line),   break,   end
    lines = strvcat(lines,line);
end
lines = lines(1:end-1,:);
fclose(fid);

% ---- Creating New aparc
ids = char(ctabr.struct_names);
Total = [ids  repmat('               ',[size(ids,1) 1]) num2str(nvert) repmat('  ',[size(ids,1) 1]) ...
        num2str(rarea) repmat('  ',[size(ids,1) 1]) num2str(rvols) repmat('  ',[size(ids,1) 1])...
        num2str(rcthm) repmat('  ',[size(ids,1) 1]) num2str(rcths) repmat('  ',[size(ids,1) 1])...
        num2str(abs(rcurv)) repmat('  ',[size(ids,1) 1]) num2str(zeros(size(ids,1),1)) repmat('  ',[size(ids,1) 1])...
        num2str(zeros(size(ids,1),1)) repmat('  ',[size(ids,1) 1]) num2str(zeros(size(ids,1),1))];
lines = strvcat(lines,Total);

% ==================== Saving the New Color File ======================= %
fid = fopen(NRaparc,'wt');
Ns = size(lines,1);
for i = 1:Ns
    fprintf(fid, '%s\n',deblank(lines(i,:)));
end
fclose(fid);
% ============== End of Saving the New Color File ====================== %
% ======================End of Left Hemisphere ==========================%
OutStats = strvcat(NLaparc,NRaparc);

return;

