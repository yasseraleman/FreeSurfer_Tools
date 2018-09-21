function OutFigures = FreeSurfer_Qcontrol(FreeDir,Id,plotopts);
%
% Syntax :
% OutFigures = FreeSurfer_Qcontrol(FreeDir,Id,plotopts);
%
% This function generates quality control images(*.jpg) using freesurfer outputs
%
% Input Parameters:
%     FreeDir       : Freesurfer Subjects Directory
%       Id          : Subject Id
%      plotopts     :  Plot options
%
% Output Parameters:
%     OutFigures          : Structure Variable containing all figures
%                           filenames
%   --- QcontrolDir ---   : Quality Control Output Directory
%   --- BrainMask ---     : Brain Mask Extraction Images
%    --- APARC ---        : Automatic Structures Segmentation Images (aparc)
%    --- APARCs ---       : Automatic Structures Segmentation Images (aparc.a2009s)
%     --- VPS ---         : Surface Extraction Images
%     --- LHP ---         : Left Hemisphere Pial Surface overlayed wtih
%                           thickness, curv, aparc.annot and
%                           aparc.a2009s.annot
%     --- LHW ---         : Left Hemisphere White Surface overlayed wtih
%                           thickness, curv, aparc.annot and
%                           aparc.a2009s.annot
%     --- LHI ---         : Left Hemisphere Inflated Surface overlayed wtih
%                           thickness, curv, aparc.annot and
%                           aparc.a2009s.annot
%     --- RHP ---         : Right Hemisphere Pial Surface overlayed wtih
%                           thickness, curv, aparc.annot and
%                           aparc.a2009s.annot
%     --- RHW ---         : Right Hemisphere White Surface overlayed wtih
%                           thickness, curv, aparc.annot and
%                           aparc.a2009s.annot
%     --- RHI ---         : Right Hemisphere Inflated Surface overlayed wtih
%                           thickness, curv, aparc.annot and
%                           aparc.a2009s.annot
%
% Related references:
%
% See also:
%__________________________________________________
% Authors:  Yasser Aleman Gomez
% Unidad de Medicina y Cirugía Experimental, Hospital General Universitario Gregorio Marañón, Madrid
%
% February 17th 2012
% Version $1.0

%  FreeDir = '/media/Data/PROCESSING_RESULTS/ASPERGER/5-freesurfer_processing';
%  Id = 'ASPER_00001__101-20060510';

opts.pipe.freesdir = FreeDir;
opts.pipe.subjId = Id;
if nargin < 3
    plotopts.type = 'pial+white+color';
    plotopts.subc = 0;
    plotopts.nslices = 3;
end
color = [1 0 0;0 1 0; 0 0 1; 1 1 0;0 1 1;1 0 1;1 0.5 0;0 0.5 1;1 0 0.5;0.5 1 0; 0.25 0.5 1;0.6 0.3 0.52;0.7 0.5 0.9];
warning off;

% Inicialization as Empty Sctruct
OutFigures.BrainMask = '';
OutFigures.APARC = '';
OutFigures.APARCs = '';
OutFigures.VPS = '';
OutFigures.LHP = '';
OutFigures.LHW = '';
OutFigures.LHI = '';
OutFigures.RHP = '';
OutFigures.RHW = '';
OutFigures.RHI = '';
OutFigures.Notes = '';
cont = 0;

%% ================== Cheking Input Parameters ========================= %%
% TempVexist: Variable to check the existence of needed Volume Files
TempVexist = logical([exist([opts.pipe.freesdir filesep opts.pipe.subjId filesep 'mri' filesep 'T1.mgz'],'file')...
    exist([opts.pipe.freesdir filesep opts.pipe.subjId filesep 'mri' filesep 'brainmask.mgz'],'file')...
    exist([opts.pipe.freesdir filesep opts.pipe.subjId filesep 'mri' filesep 'aparc+aseg.mgz'],'file')...
    exist([opts.pipe.freesdir filesep opts.pipe.subjId filesep 'mri' filesep 'aparc.a2009s+aseg.mgz'],'file')]);
if sum(TempVexist) == 0
    error(['Quality control can not be done for subject ' opts.pipe.subjId]);
    return;
end
[Outboolean] =FreeSurfer_verif(opts.pipe.freesdir,opts.pipe.subjId);
if ~Outboolean
    warndlg('FreeSurfer Processing was not finished. Please run FreeSurfer to guarantee that all output files exist');
    OutFigures.Notes = strvcat(OutFigures.Notes,'Note 0: There are some missing files in FreeSurfer Output Directory');
else
    OutFigures.Notes = strvcat(OutFigures.Notes,'Note 0: All FreeSurfer output files exist');
end
%% ============ End of Cheking Input Parameters ======================== %%

%% ======================= Main Program ================================ %%

%% ======================== Reading Images ============================= %%
% -------------------------- Converting T1 -------------------------------%
if TempVexist(1) == 1
    cad = ['mri_convert -i ' opts.pipe.freesdir filesep opts.pipe.subjId filesep 'mri' filesep 'T1.mgz' ' -o ' opts.pipe.freesdir filesep opts.pipe.subjId filesep 'tmp' filesep 'T1t.nii'];
    system(cad);
    Image = flip_to_axial_nii([opts.pipe.freesdir filesep opts.pipe.subjId filesep 'tmp' filesep 'T1t.nii'],1);
else
    cont = cont + 1;
    OutFigures.Notes = strvcat(OutFigures.Notes,['Note ' num2str(cont) ': mri/T1.mgz file does not exist']);
end

% ---------------------- Converting Brain Mask ---------------------------%
if TempVexist(2) == 1
    cad = ['mri_convert -i ' opts.pipe.freesdir filesep opts.pipe.subjId filesep 'mri' filesep 'brainmask.mgz' ' -o ' opts.pipe.freesdir filesep opts.pipe.subjId filesep 'tmp' filesep 'brainmaskt.nii'];
    system(cad);
    ImageM = flip_to_axial_nii([opts.pipe.freesdir filesep opts.pipe.subjId filesep 'tmp' filesep 'brainmaskt.nii'],1);
else
    cont = cont + 1;
    OutFigures.Notes = strvcat(OutFigures.Notes,['Note ' num2str(cont) ': mri/brainmask.mgz file does not exist']);
end

% ---------------------- Converting Aparc+Aseg Mask ----------------------%
if TempVexist(3) == 1
    cad = ['mri_convert -i ' opts.pipe.freesdir filesep opts.pipe.subjId filesep 'mri' filesep 'aparc+aseg.mgz' ' -o ' opts.pipe.freesdir filesep opts.pipe.subjId filesep 'tmp' filesep 'aparc+asegt.nii'];
    system(cad);
    ImageA = flip_to_axial_nii([opts.pipe.freesdir filesep opts.pipe.subjId filesep 'tmp' filesep 'aparc+asegt.nii'],1);
else
    cont = cont + 1;
    OutFigures.Notes = strvcat(OutFigures.Notes,['Note ' num2str(cont) ': mri/aparc+aseg.mgz file does not exist']);
end

% ---------------------- Converting Aparc2009s+Aseg Mask -----------------%
if TempVexist(4) == 1
    cad = ['mri_convert -i ' opts.pipe.freesdir filesep opts.pipe.subjId filesep 'mri' filesep 'aparc.a2009s+aseg.mgz' ' -o ' opts.pipe.freesdir filesep opts.pipe.subjId filesep 'tmp' filesep 'aparc.a2009s+asegt.nii'];
    system(cad);
    ImageAs = flip_to_axial_nii([opts.pipe.freesdir filesep opts.pipe.subjId filesep 'tmp' filesep 'aparc.a2009s+asegt.nii'],1);
else
    cont = cont + 1;
    OutFigures.Notes = strvcat(OutFigures.Notes,['Note ' num2str(cont) ': mri/aparc.a2009s+aseg.mgz file does not exist']);
end
%% =================== End of Reading Images =========================== %%

%% ================== Extracting Subcortical Structures ================ %%
if sum(TempVexist(3:4)) == 0;
    plotopts.subc = 0; % No Surfaces for Subcortical Structures
elseif (TempVexist(3) == 1)&(TempVexist(4) == 0);
    ImageAt = ImageA;
elseif (TempVexist(3) == 0)&(TempVexist(4) == 1);
    ImageAt = ImageAs;
else (TempVexist(3) == 1)&(TempVexist(4) == 1);
    ImageAt = ImageA;
end
if plotopts.subc == 1
    [SurfSub] = Surf_Extraction(ImageAt,[10:13 17:18 26 49:54 58],'vox');
    switch plotopts.type
        case 'pial+white'
            for i = 1:length(SurfSub)
                SurfSub(i).color = [0.6 0.3 0.52];
            end
        case 'pial+white+color'
            SurfSub(1).color = color(1,:); SurfSub(8).color = color(1,:);
            SurfSub(2).color = color(2,:); SurfSub(9).color = color(2,:);
            SurfSub(3).color = color(3,:); SurfSub(10).color = color(3,:);
            SurfSub(4).color = color(4,:); SurfSub(11).color = color(4,:);
            SurfSub(5).color = color(5,:); SurfSub(12).color = color(5,:);
            SurfSub(6).color = color(6,:); SurfSub(13).color = color(6,:);
            SurfSub(7).color = color(7,:); SurfSub(14).color = color(7,:);
    end
end
%% ================ End of Extracting Subcortical Structures =========== %%
%% ========================== Creating Slices ========================== %%
% ---------------------- Verifying Brain Mask ----------------------------%
if TempVexist(2) == 1
    ImageMt = ImageM;
elseif (TempVexist(2) == 0)&(TempVexist(3) == 1)
    ImageMt = ImageA;
elseif (TempVexist(2) == 0)&(TempVexist(3) == 0)&(TempVexist(4) == 1)
    ImageMt = ImageAs;
else 
    ImageMt = Image;
    Nomask = 1; % Boolean Variable to plot (0) or not (1) brain mask
end
% ---------------------- Detecting Slices --------------------------------%
VM = spm_vol(ImageMt);
IM = logical(spm_read_vols(VM));
% Sagital Slice
temp = sum(sum(IM,3),2);
ind = find(temp == max(temp));
xslice = [ind(1) 0 0];
indx = find(temp);
Xmin = min(indx);
Xmax = max(indx);

% Coronal Slice
temp = sum(sum(IM,1),3);
ind = find(temp == max(temp));
yslice = [0 ind(1)  0];
indy = find(temp);
Ymin = min(indy);
Ymax = max(indy);

% Axial Slice
temp = sum(sum(IM,1),2);
ind = find(temp == max(temp));
zslice = [0 0 ind(1)];
indz = find(temp);
Zmin = min(indz);
Zmax = max(indz);

% All Slices
if plotopts.nslices == 1
    slices = [xslice;yslice;zslice]; % Maximun Number of points in each slice
else
    % Equidistant Sagital Slices
    if plotopts.nslices < (Xmax-Xmin-2)
        inc = (Xmax-Xmin)/(plotopts.nslices+1);
        Xslices = round([Xmin:inc:Xmax]);
    else
        Xslices = round([Xmin:Xmax]);
    end
    Xslices = Xslices(2:end-1);
    
    % Equidistant Coronal Slices
    if plotopts.nslices < (Ymax-Ymin-2)
        inc = (Ymax-Ymin)/(plotopts.nslices+1);
        Yslices = round([Ymin:inc:Ymax]);
    else
        Yslices = round([Ymin:Ymax]);
    end
    Yslices = Yslices(2:end-1);
    
    % Equidistant Axial Slices
    if plotopts.nslices < (Ymax-Ymin-2)
        inc = (Zmax-Zmin)/(plotopts.nslices+1);
        Zslices = round([Zmin:inc:Zmax]);
    else
        Zslices = round([Zmin:Zmax]);
    end
    Zslices = Zslices(2:end-1);
    
    % All Slices
    slices = [Xslices(:) Yslices(:) Zslices(:)];
end
%% =================== End of Creating Slices ========================== %%

%% ================== Reading Talairach Information ==================== %%
talfile = [opts.pipe.freesdir filesep opts.pipe.subjId filesep 'mri' filesep 'transforms' filesep 'talairach.lta'];
if exist(talfile,'file')
    cras1 = textread(talfile,'%s',5,'headerlines',20);
    cras = char(cras1);
    cras = [str2num(cras(3,:))  str2num(cras(4,:)) str2num(cras(5,:))];
else
    cont = cont + 1;
    cras = [0 0 0];
    OutFigures.Notes = strvcat(OutFigures.Notes,['Note ' num2str(cont) ': mri/transforms/talairach.lta file does not exist. Cras  = 0 0 0. This can lead to mismatch errors between surfaces and volumes']);
end
%% ============= End of Reading Talairach Information ================== %%

%% ====================== Reading Surfaces ============================= %%
% TempSexist: Variable to check the existence of needed Surface Files
TempSexist = logical([exist([opts.pipe.freesdir filesep opts.pipe.subjId filesep 'surf' filesep 'lh.pial'],'file')...
                      exist([opts.pipe.freesdir filesep opts.pipe.subjId filesep 'surf' filesep 'lh.white'],'file')...
                      exist([opts.pipe.freesdir filesep opts.pipe.subjId filesep 'surf' filesep 'rh.pial'],'file')...
                      exist([opts.pipe.freesdir filesep opts.pipe.subjId filesep 'surf' filesep 'rh.white'],'file')]);
ASurfcad = '';
% ---------- Reading Left Hemisphere Pial Surface ------------------------%
ASurf = '';
if TempSexist(1) == 1
    [OutFiles, SurfF] = Exp_Surf([opts.pipe.freesdir filesep opts.pipe.subjId filesep 'surf' filesep 'lh.pial'], '0', '','', 'imp','n');
    Surftl= SurfF{1};
    Surftl.Name = 'LH.Pial';
    Surftl.SurfData.vertices = Surftl.SurfData.vertices + repmat(cras,[size(Surftl.SurfData.vertices,1) 1]);
    Surftl.color = color(1,:);
    ASurf = [ASurf Surftl];
else
    cont = cont + 1;
    OutFigures.Notes = strvcat(OutFigures.Notes,['Note ' num2str(cont) ': surf/lh.pial file does not exist']);
end

% ---------- Reading Left Hemisphere White Surface -----------------------%
if TempSexist(2) == 1
    [OutFiles, SurfF] = Exp_Surf([opts.pipe.freesdir filesep opts.pipe.subjId filesep 'surf' filesep 'lh.white'], '0', '','', 'imp','n');
    Surftlw= SurfF{1};
    Surftlw.Name = 'LH.WHITE';
    Surftlw.SurfData.vertices = Surftlw.SurfData.vertices + repmat(cras,[size(Surftlw.SurfData.vertices,1) 1]);
    Surftlw.color = color(2,:);
    ASurf = [ASurf Surftlw];
else
    cont = cont + 1;
    OutFigures.Notes = strvcat(OutFigures.Notes,['Note ' num2str(cont) ': surf/lh.white file does not exist']);
end

% ---------- Reading Right Hemisphere Pial Surface -----------------------%
if TempSexist(3) == 1
    [OutFiles, SurfF] = Exp_Surf([opts.pipe.freesdir filesep opts.pipe.subjId filesep 'surf' filesep 'rh.pial'], '0', '','', 'imp','n');
    Surftr= SurfF{1};
    Surftr.Name = 'RH.Pial';
    Surftr.SurfData.vertices = Surftr.SurfData.vertices + repmat(cras,[size(Surftr.SurfData.vertices,1) 1]);
    Surftr.color = color(3,:);
    ASurf = [ASurf Surftr];
else
    cont = cont + 1;
    OutFigures.Notes = strvcat(OutFigures.Notes,['Note ' num2str(cont) ': surf/rh.pial file does not exist']);
end
% ---------- Reading Right Hemisphere White Surface ----------------------%
if TempSexist(4) == 1
    [OutFiles, SurfF] = Exp_Surf([opts.pipe.freesdir filesep opts.pipe.subjId filesep 'surf' filesep 'rh.white'], '0', '','', 'imp','n');
    Surftrw= SurfF{1};
    Surftrw.Name = 'RH.White';
    Surftrw.SurfData.vertices = Surftrw.SurfData.vertices + repmat(cras,[size(Surftrw.SurfData.vertices,1) 1]);
    Surftrw.color = color(4,:);
    ASurf = [ASurf Surftrw];
else
    cont = cont + 1;
    OutFigures.Notes = strvcat(OutFigures.Notes,['Note ' num2str(cont) ': surf/rh.white file does not exist']);
end
%% ================== End of Reading Surfaces ========================== %%

%% ====================== Selecting Surfaces to Plot =================== %%
Surfa = '';
switch plotopts.type
    case 'pial+white'
        % Joining Surfaces
        if plotopts.subc == 1
            Surfa = [ASurf SurfSub]';
        else
            Surfa = [ASurf]';
        end
    case 'pial+white+color'
        % -------------------- Dividing Left Hemisphere ------------------%
        if exist([opts.pipe.freesdir filesep opts.pipe.subjId filesep 'label' filesep 'lh.aparc.annot'],'file')
            [txt,colortable] = read_cfiles([opts.pipe.freesdir filesep opts.pipe.subjId filesep 'label' filesep 'lh.aparc.annot']);
            if TempSexist(1) == 1
                % ------------ Dividing Pial Left Hemisphere -------------%
                Surftl.Is = txt;
                [Surfoutl] = Divide_Surf(Surftl, colortable);
            else
                Surfoutl = '';
            end
            if TempSexist(2) == 1
                % ------------ Dividing White Left Hemisphere ------------%
                Surftlw.Is = txt;
                [Surfoutlw] = Divide_Surf(Surftlw, colortable);
            else
                Surfoutlw = '';
            end
        else
            cont = cont + 1;
            OutFigures.Notes = strvcat(OutFigures.Notes,['Note ' num2str(cont) ': label/lh.aparc.white file does not exist']);
            Surfoutl = '';
            Surfoutlw = '';
        end
        % -------------------- Dividing Right Hemisphere -----------------%
        if exist([opts.pipe.freesdir filesep opts.pipe.subjId filesep 'label' filesep 'rh.aparc.annot'],'file')
            [txt,colortable] = read_cfiles([opts.pipe.freesdir filesep opts.pipe.subjId filesep 'label' filesep 'rh.aparc.annot']);
            if TempSexist(3) == 1
                % ------------ Dividing Pial Left Hemisphere -------------%
                Surftr.Is = txt;
                [Surfoutr] = Divide_Surf(Surftr, colortable);
            else
                Surfoutr = '';
            end
            if TempSexist(4) == 1
                % ------------ Dividing White Left Hemisphere ------------%
                Surftrw.Is = txt;
                [Surfoutrw] = Divide_Surf(Surftrw, colortable);
            else
                Surfoutrw = '';
            end
        else
            cont = cont + 1;
            OutFigures.Notes = strvcat(OutFigures.Notes,['Note ' num2str(cont) ': label/rh.aparc.white file does not exist']);
            Surfourl = '';
            Surfoutrw = '';
        end
       
        % Joining Surfaces
        if plotopts.subc == 1
            for k = 1:length(SurfSub)
                Surft = SurfSub(k);
                Surft = rmfield(Surft,'Area');
                Surft = rmfield(Surft,'Tri');
                Surft = rmfield(Surft,'Orig');
                Surft = rmfield(Surft,'Dim');
                Surft = rmfield(Surft,'VoxSize');
                Surft = rmfield(Surft,'Code');
                Surft.Is = k*ones(size(Surft.SurfData.vertices,1),1);
                SurfSubt(k) = Surft;
            end
            Surfa = [Surfoutl Surfoutlw Surfoutr Surfoutrw SurfSubt]'; clear SurfSub;
        else
            Surfa = [Surfoutl Surfoutlw Surfoutr Surfoutrw]';
        end
end
%% ======================== Creating Figures =========================== %%
mkdir([opts.pipe.freesdir filesep opts.pipe.subjId],'qcontrol');
OutFigures.QcontrolDir = [opts.pipe.freesdir filesep opts.pipe.subjId filesep 'qcontrol'];
system(['rm -r ' opts.pipe.freesdir filesep opts.pipe.subjId filesep 'qcontrol' filesep '*']);
% 1. Brain Mask Extraction
plotopts.transp = 0.7;
plotopts.colormap = 'jet';
plotopts.itype = 'mask';
plotopts.atype = 'aparc+aseg';
plotopts.dtype.type = 'vectorial';
plotopts.dtype.modul = 0;
plotopts.linelength = 1/2;
if TempVexist(1)&TempVexist(2)
    FiguresBM = Image_plus_overlay(Image, ImageM, slices, plotopts, [opts.pipe.freesdir filesep opts.pipe.subjId filesep 'qcontrol']);
elseif ~TempVexist(1)&TempVexist(2)
    FiguresBM = Image_plus_overlay(ImageM, ImageM, slices, plotopts, [opts.pipe.freesdir filesep opts.pipe.subjId filesep 'qcontrol']);
else
    FiguresBM = '';
end

% 2. Aparc Segmentation
plotopts.transp = 0.7;
plotopts.colormap = 'jet';
plotopts.itype = 'atlas';
plotopts.atype = 'aparc+aseg';
plotopts.dtype.type = 'vectorial';
plotopts.dtype.modul = 0;
plotopts.linelength = 1/2;
if TempVexist(1)&TempVexist(3)
    FiguresAPARC = Image_plus_overlay(Image, ImageA, slices, plotopts, [opts.pipe.freesdir filesep opts.pipe.subjId filesep 'qcontrol']);
elseif ~TempVexist(1)&TempVexist(3)
    plotopts.transp = 1;
    FiguresAPARC = Image_plus_overlay(ImageA, ImageA, slices, plotopts, [opts.pipe.freesdir filesep opts.pipe.subjId filesep 'qcontrol']);
else
    FiguresAPARC = '';
end

% 3. Aparc.a2009s Segmentation
plotopts.transp = 0.7;
plotopts.colormap = 'jet';
plotopts.itype = 'atlas';
plotopts.atype = 'a2009s+aseg';
plotopts.dtype.type = 'vectorial';
plotopts.dtype.modul = 0;
plotopts.linelength = 1/2;
if TempVexist(1)&TempVexist(4)
    FiguresAPARCs = Image_plus_overlay(Image, ImageAs, slices, plotopts, [opts.pipe.freesdir filesep opts.pipe.subjId filesep 'qcontrol']);
elseif ~TempVexist(1)&TempVexist(3)
    plotopts.transp = 1;
    FiguresAPARCs = Image_plus_overlay(ImageAs, ImageAs, slices, plotopts, [opts.pipe.freesdir filesep opts.pipe.subjId filesep 'qcontrol']);
else
    FiguresAPARCs = '';
end
% ----------------------------------------------------------------------- %

% 4. Surface Extraction
if ~isempty(Surfa)&TempVexist(1)
    FiguresVPS = Image_plus_surface(Image,Surfa, slices, [opts.pipe.freesdir filesep opts.pipe.subjId filesep 'qcontrol']);
else
    FiguresVPS = '';
end

% 5. Characteristics
%% ====================== Left Hemisphere ============================== %%
Characs = strvcat([opts.pipe.freesdir filesep opts.pipe.subjId filesep 'surf' filesep 'lh.curv'],...
                  [opts.pipe.freesdir filesep opts.pipe.subjId filesep 'label' filesep 'lh.aparc.annot'],...
                  [opts.pipe.freesdir filesep opts.pipe.subjId filesep 'label' filesep 'lh.aparc.a2009s.annot'],...
                  [opts.pipe.freesdir filesep opts.pipe.subjId filesep 'surf' filesep 'lh.thickness']);
              
% TempCexist: Variable to check the existence of needed Characteristic Files
TempCexist = logical([exist([opts.pipe.freesdir filesep opts.pipe.subjId filesep 'surf' filesep 'lh.curv'],'file')...
                      exist([opts.pipe.freesdir filesep opts.pipe.subjId filesep 'label' filesep 'lh.aparc.annot'],'file')...
                      exist([opts.pipe.freesdir filesep opts.pipe.subjId filesep 'label' filesep 'lh.aparc.a2009s.annot'],'file')...
                      exist([opts.pipe.freesdir filesep opts.pipe.subjId filesep 'surf' filesep 'lh.thickness'],'file')]);
                  
% ----------- Verifying the existence of characteristics files -----------%
if ~TempCexist(1)
    cont = cont + 1;
    OutFigures.Notes = strvcat(OutFigures.Notes,['Note ' num2str(cont) ': surf/lh.curv file does not exist']);
end
if ~TempCexist(2)
    cont = cont + 1;
    OutFigures.Notes = strvcat(OutFigures.Notes,['Note ' num2str(cont) ': label/lh.aparc.annot file does not exist']);
end
if ~TempCexist(3)
    cont = cont + 1;
    OutFigures.Notes = strvcat(OutFigures.Notes,['Note ' num2str(cont) ': label/lh.aparc.a2009s.annot file does not exist']);
end
if ~TempCexist(4)
    cont = cont + 1;
    OutFigures.Notes = strvcat(OutFigures.Notes,['Note ' num2str(cont) ': surf/lh.thickness file does not exist']);
end

%   a) Pial Left Hemisphere
if  exist([opts.pipe.freesdir filesep opts.pipe.subjId filesep 'surf' filesep 'lh.pial'],'file')
    if sum(TempCexist)>0
        Characs = Characs(TempCexist,:);
        FiguresLHP = Surface_plus_overlay([opts.pipe.freesdir filesep opts.pipe.subjId filesep 'surf' filesep 'lh.pial'],Characs , [opts.pipe.freesdir filesep opts.pipe.subjId filesep 'qcontrol']);
    else
        FiguresLHP = '';
    end
else
    FiguresLHP = '';
end

%  b) White Left Hemisphere
if  exist([opts.pipe.freesdir filesep opts.pipe.subjId filesep 'surf' filesep 'lh.white'],'file')
    if sum(TempCexist)>0
        Characs = Characs(TempCexist,:);
        FiguresLHW = Surface_plus_overlay([opts.pipe.freesdir filesep opts.pipe.subjId filesep 'surf' filesep 'lh.white'],Characs , [opts.pipe.freesdir filesep opts.pipe.subjId filesep 'qcontrol']);
    else
        FiguresLHW = '';
    end
else
    FiguresLHW = '';
end

%  c) Inflated Left Hemisphere
if  exist([opts.pipe.freesdir filesep opts.pipe.subjId filesep 'surf' filesep 'lh.inflated'],'file')
    if sum(TempCexist)>0
        Characs = Characs(TempCexist,:);
        FiguresLHI = Surface_plus_overlay([opts.pipe.freesdir filesep opts.pipe.subjId filesep 'surf' filesep 'lh.inflated'],Characs , [opts.pipe.freesdir filesep opts.pipe.subjId filesep 'qcontrol']);
    else
        FiguresLHI = '';
    end
else
    FiguresLHI = '';
end
%% =================== End of Left Hemisphere ========================== %%

%% ===================== Right Hemisphere ============================== %%
Characs = strvcat([opts.pipe.freesdir filesep opts.pipe.subjId filesep 'surf' filesep 'rh.curv'],...
    [opts.pipe.freesdir filesep opts.pipe.subjId filesep 'label' filesep 'rh.aparc.annot'],...
    [opts.pipe.freesdir filesep opts.pipe.subjId filesep 'label' filesep 'rh.aparc.a2009s.annot'],...
    [opts.pipe.freesdir filesep opts.pipe.subjId filesep 'surf' filesep 'rh.thickness']);

% TempCexist: Variable to check the existence of needed Characteristic Files
TempCexist = logical([exist([opts.pipe.freesdir filesep opts.pipe.subjId filesep 'surf' filesep 'rh.curv'],'file')...
                      exist([opts.pipe.freesdir filesep opts.pipe.subjId filesep 'label' filesep 'rh.aparc.annot'],'file')...
                      exist([opts.pipe.freesdir filesep opts.pipe.subjId filesep 'label' filesep 'rh.aparc.a2009s.annot'],'file')...
                      exist([opts.pipe.freesdir filesep opts.pipe.subjId filesep 'surf' filesep 'rh.thickness'],'file')]);
                  
% ----------- Verifying the existence of characteristics files -----------%
if ~TempCexist(1)
    cont = cont + 1;
    OutFigures.Notes = strvcat(OutFigures.Notes,['Note ' num2str(cont) ': surf/rh.curv file does not exist']);
end
if ~TempCexist(2)
    cont = cont + 1;
    OutFigures.Notes = strvcat(OutFigures.Notes,['Note ' num2str(cont) ': label/rh.aparc.annot file does not exist']);
end
if ~TempCexist(3)
    cont = cont + 1;
    OutFigures.Notes = strvcat(OutFigures.Notes,['Note ' num2str(cont) ': label/rh.aparc.a2009s.annot file does not exist']);
end
if ~TempCexist(4)
    cont = cont + 1;
    OutFigures.Notes = strvcat(OutFigures.Notes,['Note ' num2str(cont) ': surf/rh.thickness file does not exist']);
end

%   a) Pial Right Hemisphere
if  exist([opts.pipe.freesdir filesep opts.pipe.subjId filesep 'surf' filesep 'rh.pial'],'file')
    if sum(TempCexist)>0
        Characs = Characs(TempCexist,:);
        FiguresRHP = Surface_plus_overlay([opts.pipe.freesdir filesep opts.pipe.subjId filesep 'surf' filesep 'rh.pial'],Characs , [opts.pipe.freesdir filesep opts.pipe.subjId filesep 'qcontrol']);
    else
        FiguresRHP = '';
    end
else
    FiguresRHP = '';
end

%  b) White Right Hemisphere
if  exist([opts.pipe.freesdir filesep opts.pipe.subjId filesep 'surf' filesep 'rh.white'],'file')
    if sum(TempCexist)>0
        Characs = Characs(TempCexist,:);
        FiguresRHW = Surface_plus_overlay([opts.pipe.freesdir filesep opts.pipe.subjId filesep 'surf' filesep 'rh.white'],Characs , [opts.pipe.freesdir filesep opts.pipe.subjId filesep 'qcontrol']);
    else
        FiguresRHW = '';
    end
else
    FiguresRHW = '';
end

%  c) Inflated Right Hemisphere
if  exist([opts.pipe.freesdir filesep opts.pipe.subjId filesep 'surf' filesep 'rh.inflated'],'file')
    if sum(TempCexist)>0
        Characs = Characs(TempCexist,:);
        FiguresRHI = Surface_plus_overlay([opts.pipe.freesdir filesep opts.pipe.subjId filesep 'surf' filesep 'rh.inflated'],Characs , [opts.pipe.freesdir filesep opts.pipe.subjId filesep 'qcontrol']);
    else
        FiguresRHI = '';
    end
else
    FiguresRHI = '';
end
%% ===================== End of Right Hemisphere ======================= %%

OutFigures.BrainMask = FiguresBM;
OutFigures.APARC = FiguresAPARC;
OutFigures.APARCs = FiguresAPARCs;
OutFigures.VPS = FiguresVPS;
OutFigures.LHP = FiguresLHP;
OutFigures.LHW = FiguresLHW;
OutFigures.LHI = FiguresLHI;
OutFigures.RHP = FiguresRHP;
OutFigures.RHW = FiguresRHW;
OutFigures.RHI = FiguresRHI;

%% =================== Removing Temporary Files ======================== %%
if exist(Image,'file')
    delete(Image);
end
if exist(ImageM,'file')
    delete(ImageM);
end
if exist(ImageA,'file')
    delete(ImageA);
end
if exist(ImageAs,'file')
    delete(ImageAs);
end
%% =================== End of Removing Temporary Files ================= %%
return

%% ==================== End of Main Program ============================ %%
%=======================Internal functions================================%

function [Xline,Yline,Zline] = cut_surf(plane,fv)

% This function cuts a surface (defined by vertices and faces as represented
% by MATLAB patches) by a plane, leading to a curve in 3D space. The
% resulting curve is represented by a set of contigous lines in the space
%
% Syntax:
% [Xline,Yline,Zline] = cut_surf(plane,SurfData)
%
% INPUTS:
% plane : A 4-length vector with the parameters of the plane. If plane = [A
% B C D] then every 3D point P = (x,y,z) belonging to the plane satisfies
% plane*[P; 1]' = A*x + B*y + C*z + D = 0
% SurfData : surface structure as represented in MATLAB by patches:
%            SurfData.vertices
%            SurfData.faces
%
% OUTPUTS:
% Xline,Yline,Zline: Matrices with the line coordinates.
% The entire curve can be plotted by simply typing:
% line(Xline,Yline,Zline,'Properties1',Value1,...);
%
% Pedro Antonio Vald�s Hern�ndez
%
% October 29 2008

warning off %#ok
Xline = []; Yline = []; Zline = [];
oo = ones(size(fv.vertices,1),1);
maxdist = sqrt(max(sum((fv.vertices(fv.faces(:,1),:)-fv.vertices(fv.faces(:,2),:)).^2,2)));
maxdist = max(maxdist,sqrt(max(sum((fv.vertices(fv.faces(:,1),:)-fv.vertices(fv.faces(:,3),:)).^2,2))));
maxdist = max(maxdist,sqrt(max(sum((fv.vertices(fv.faces(:,2),:)-fv.vertices(fv.faces(:,3),:)).^2,2))));
vertx = find(abs(dot([fv.vertices oo],plane(oo,:),2)/norm(plane(1:3)))<maxdist);
indf = ismember(fv.faces,vertx);
[rindf,cindf] = find(indf); %#ok
rindf = unique(rindf);
Nf = length(rindf);
% h = waitbar(0,'cutting surface...');
for i = 1:Nf
    verts = fv.vertices(fv.faces(rindf(i),:),:);
    verts(:,4) = 1;
    diffv(1,:) = diff(verts([1 2],:));
    diffv(2,:) = diff(verts([2 3],:));
    diffv(3,:) = diff(verts([3 1],:));
    alpha = -verts*plane'./(diffv*plane');
    % NaN   : contains
    % < 0   : not contains down
    % -Inf  : parallel down
    % > 1   : not contains up
    % +Inf  : parallel up
    ind = find((alpha<1 & alpha >=0) | (alpha<=1 & alpha >0))  ;
    if ~isempty(ind) && length(ind)==2
        points = verts(ind,1:3) + alpha(ind,[1 1 1]).*diffv(ind,1:3);
        Xline = [Xline points(:,1)]; %#ok
        Yline = [Yline points(:,2)]; %#ok
        Zline = [Zline points(:,3)]; %#ok
    end
    %    waitbar(i/Nf,h);
end
% close(h);
warning on %#ok
return

function R = getrot(A)

M = diag(sqrt(sum(A.^2))); A = A*inv(M); [s,ind] = max(abs(A)); %#ok
R(3,3) = 0; for i = 1:3, R(ind(i),i) = sign(A(ind(i),i)); end
return

function [Surfout] = Divide_Surf(Surf, colortable);

txt = Surf.Is;
sts = unique(txt);
Ns = length(sts);
Nfaces = 0;
for i = 1:Ns
    ind = find(txt==sts(i));
    indc = find(colortable.table(:,5)==sts(i));
    if isempty(indc)
        Matname = 'unknown_structure';
        Color =   [1 1 1];      % Color
    else
        Matname = char(colortable.struct_names{indc});
        Color =  colortable.table(indc,1:3)/255;       % Color
    end
    Surft = Surf;
    Faces = Surf.SurfData.faces;
    a = ismember(Surf.SurfData.faces,ind);
    ind2 = find(sum(a')==0);
    Surft.SurfData.faces(ind2,:) = [];
    [Surft] = Reorg_Surf(Surft);
    Surft.color = Color;
    Surfout(i) = Surft;
end
return
%=========================================================================%
