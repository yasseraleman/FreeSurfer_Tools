function HTML_File = Create_FreeSurfer_Subject_HTML(FreeDir,Id,OutFigures);
%
% Syntax :
%   HTML_File = Create_FreeSurfer_HTML(FreeDir,Id,HTML_File,OutFigures);
%
% This script creates a HTML file as quality control file. The HTML shows
% results for different processing steps.
%
% Input Parameters:
%      FreeDir            : FreeSurfer Subjects Directory
%       Id                : Subject Id
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
% Output Parameters:
%     HTML_File          :  Output HTML Filename
%
% Related references:
%
%
% See also: Graph_Fiber_Tracking Probtrack_Fiber_Tracking
% Trackvis_Connectivity Diffusion_HARDI_ODF_Estimation
% Diffusion_Tensor_Estimation DWI_Correction_Pipeline Atlas_To_Diff_Space
% Graph_Fiber_Tracking Save_Connectivity_Matrix Read_Connectivity_Matrix Plot_Matrix
%
%
%__________________________________________________
% Authors: Yasser Aleman Gomez
% LIM, HUGGM
% August 14th 2012
% Version $1.0

opts.pipe.freesdir = FreeDir;
opts.pipe.subjId = Id;
if nargin < 3
    OutFigures = FreeSurfer_Qcontrol(FreeDir,Id);
end
mkdir(OutFigures.QcontrolDir);
HTML_File = [OutFigures.QcontrolDir filesep opts.pipe.subjId '-FS_index.html'];
fid  = fopen(HTML_File,'wt');
% End of creating Folders and files

%% ===================== Creating Header ================================%%
fprintf(fid,'%s\n','<!DOCTYPE html PUBLIC "-//W3C//DTD XHTML 1.0 Transitional//EN" "http://www.w3.org/TR/xhtml1/DTD/xhtml1-transitional.dtd">');
fprintf(fid,'%s\n','<html xmlns="http://www.w3.org/1999/xhtml">');
fprintf(fid,'%s\n','<head>');
fprintf(fid,'%s\n','<meta http-equiv="Content-Type" content="text/html; charset=iso-8859-1" />');
fprintf(fid,'%s\n',['<title>HTML File for Subject ' opts.pipe.subjId ' </title>']);
fprintf(fid,'%s\n','<link rel="stylesheet" type="text/css" href="main_style.css" />');
fprintf(fid,'%s\n','</head>');

% General information
fprintf(fid,'%s\n','<body bgcolor="#FFFFCC"> </body>'); % Background Color
fprintf(fid,'%s\n',['<h1>' opts.pipe.subjId ' </h1>']);
fprintf(fid,'%s\n',['<p><strong>FreeSurfer Output Directory:</strong> ' [opts.pipe.freesdir filesep opts.pipe.subjId] ' </p>']);

% Stats Links
% Aseg Stats
asegstat = [opts.pipe.freesdir filesep opts.pipe.subjId filesep 'stats' filesep 'aseg.stats'];
if exist(asegstat,'file')
    fprintf(fid,'%s\n',['<p><strong>Automatic Segmentation Stats File:</strong> <a href="../stats/aseg.stats"target="_blank"> aseg.stats </a> </p>']);
end

% LH Stats
asegstat = [opts.pipe.freesdir filesep opts.pipe.subjId filesep 'stats' filesep 'lh.aparc.stats'];
if exist(asegstat,'file')
    fprintf(fid,'%s\n',['<p><strong>Left Hemisphere APARC Stats File:</strong> <a href="../stats/lh.aparc.stats"target="_blank">  lh.aparc.stats </a> </p>']);
end
% RH Stats
asegstat = [opts.pipe.freesdir filesep opts.pipe.subjId filesep 'stats' filesep 'rh.aparc.stats'];
if exist(asegstat,'file')
    fprintf(fid,'%s\n',['<p><strong>Right Hemisphere APARC Stats File:</strong> <a href="../stats/rh.aparc.stats"target="_blank"> rh.aparc.stats </a> </p>']);
end

% LH A2009s Stats
asegstat = [opts.pipe.freesdir filesep opts.pipe.subjId filesep 'stats' filesep 'lh.aparc.a2009s.stats'];
if exist(asegstat,'file')
    fprintf(fid,'%s\n',['<p><strong>Left Hemisphere APARC.a2009s Stats File:</strong> <a href="../stats/lh.aparc.a2009s.stats"target="_blank"> lh.aparc.a2009s.stats </a> </p>']);
end
% RH A2009s Stats
asegstat = [opts.pipe.freesdir filesep opts.pipe.subjId filesep 'stats' filesep 'rh.aparc.a2009s.stats'];
if exist(asegstat,'file')
    fprintf(fid,'%s\n',['<p><strong>Right Hemisphere APARC.a2009s Stats File:</strong> <a href="../stats/rh.aparc.a2009s.stats"target="_blank"> rh.aparc.a2009s.stats </a> </p>']);
end

%%  =================== Inserting Notes ================================ %%
for i = 0:size(OutFigures.Notes,1)-1
    Note = OutFigures.Notes(i+1,:);
    Ct = strread(Note,'%s','delimiter',':');
    Corr = char(Ct{end});
    fprintf(fid,'%s\n',['<p><strong> Note ' num2str(i) ': </strong>' Corr ' </p>']);
end
fprintf(fid,'%s\n',['<h1>                         </h1>']);
%% ==================== End of Inserting Notes ==========================%%

%% ============= Creating T1+ BrainMask Table ========================== %%
if ~isempty(OutFigures.BrainMask)
    fprintf(fid,'%s\n',['<h2> T1 + Brain Mask </h2>']);
    fprintf(fid,'%s\n',['<table width="1600" height="400" border="1">']);
    fprintf(fid,'%s\n',['<tr height="20" align="center">']);
    fprintf(fid,'%s\n',['<td><b> Axial Slices </b></td>']);
    fprintf(fid,'%s\n',['<td><b> Sagital Slices </b></td>']);
    fprintf(fid,'%s\n',['<td><b> Coronal Slices </b></td>']);
    fprintf(fid,'%s\n',['</tr>']);
    % --------------------- Detecting Slices ---------------------------------%
    SagFiles = '';conts = 0;
    CorFiles= '';contc = 0;
    AxFiles = '';conta = 0;
    for i = 1:size(OutFigures.BrainMask,1)
        [pth,nm,ext] = fileparts(deblank(OutFigures.BrainMask(i,:)));
        [IMtype] = strread(nm,'%s','delimiter','-');
        SliceN = IMtype{end};
        [SliceName,SliceNumber] = strread(char(SliceN),'%s%s','delimiter','_');
        SliceName = lower(deblank(char(SliceName)));
        switch SliceName
            case 'sagital'
                conts = conts+1;
                SagFiles = strvcat(SagFiles,deblank(OutFigures.BrainMask(i,:)));
                Numbers(conts) = str2num(char(SliceNumber));
            case 'coronal'
                contc = contc+1;
                CorFiles = strvcat(CorFiles,deblank(OutFigures.BrainMask(i,:)));
                Numberc(conts) = str2num(char(SliceNumber));
            case 'axial'
                conta = conta+1;
                AxFiles = strvcat(AxFiles,deblank(OutFigures.BrainMask(i,:)));
                Numbera(conts) = str2num(char(SliceNumber));
        end
    end
    % --------------------- Sorting Slices -----------------------------------%
    [a,b] = sort(Numbers);
    SagFiles = SagFiles(b,:);
    Numbers = Numbers(b);
    [a,b] = sort(Numberc);
    CorFiles = CorFiles(b,:);
    Numberc = Numberc(b);
    [a,b] = sort(Numbera);
    AxFiles = AxFiles(b,:);
    Numbera = Numbera(b);
    % --------------------- Table Organization -------------------------------%
    for i = 1:size(SagFiles,1);
        % ----------- Creating Column Titles Figures -------------------------%
        fprintf(fid,'%s\n',['<tr >']);
        fprintf(fid,'%s\n',['<td height="10" align="center"><I> Axial Slice  ' num2str(Numbera(i)) ' vox </I></td>']);
        fprintf(fid,'%s\n',['<td height="10" align="center"><I> Sagital Slice  ' num2str(Numbers(i)) ' vox </I></td>']);
        fprintf(fid,'%s\n',['<td height="10" align="center"><I> Coronal Slice  ' num2str(Numberc(i)) ' vox </I></td>']);
        fprintf(fid,'%s\n',['</tr>']);
        % ----------- Inserting Figures --------------------------------------%
        fprintf(fid,'%s\n',['<tr >']);
        % Axial
        [pth, nm,ext] = fileparts(AxFiles(i,:));
        T = imread(deblank(AxFiles(i,:)));
        resx = 400;
        resy = round(resx/size(T,1)*size(T,2));
        Filename = [nm ext];
        fprintf(fid,'%s\n',['<td><div align="center"><a href="' Filename '" target="_blank"><img src="' Filename ...
            '" alt="Axial Slice Not Found" width="' num2str(resy) '" height=' num2str(resx) ' border="0" /></a></div></td>']);
        % Sagital
        [pth, nm,ext] = fileparts(SagFiles(i,:));
        T = imread(deblank(SagFiles(i,:)));
        resx = 400;
        resy = round(resx/size(T,1)*size(T,2));
        Filename = [nm ext];
        fprintf(fid,'%s\n',['<td><div align="center"><a href="' Filename '" target="_blank"><img src="' Filename ...
            '" alt="Sagital Slice Not Found" width="' num2str(resy) '" height=' num2str(resx) ' border="0" /></a></div></td>']);
        % Coronal
        [pth, nm,ext] = fileparts(CorFiles(i,:));
        Filename = [nm ext];
        T = imread(deblank(CorFiles(i,:)));
        resx = 400;
        resy = round(resx/size(T,1)*size(T,2));
        fprintf(fid,'%s\n',['<td><div align="center"><a href="' Filename '" target="_blank"><img src="' Filename ...
            '" alt="Coronal Slice Not Found" width="' num2str(resy) '" height=' num2str(resx) ' border="0" /></a></div></td>']);
        fprintf(fid,'%s\n',['</tr>']);
    end
    fprintf(fid,'%s\n',['</table>']);
    %% ============= End of Creating T1 + BrainMask Table ================== %%
end

if ~isempty(OutFigures.APARC)
    %% ============= Creating T1 + Aparc+Aseg Atlas  Table ================= %%
    fprintf(fid,'%s\n',['<h2> T1 + Aparc+Aseg Atlas </h2>']);
    fprintf(fid,'%s\n',['<table width="1600" height="400" border="1">']);
    fprintf(fid,'%s\n',['<tr height="20" align="center">']);
    fprintf(fid,'%s\n',['<td><b> Axial Slices </b></td>']);
    fprintf(fid,'%s\n',['<td><b> Sagital Slices </b></td>']);
    fprintf(fid,'%s\n',['<td><b> Coronal Slices </b></td>']);
    fprintf(fid,'%s\n',['</tr>']);
    % --------------------- Detecting Slices ---------------------------------%
    SagFiles = '';conts = 0;
    CorFiles= '';contc = 0;
    AxFiles = '';conta = 0;
    for i = 1:size(OutFigures.APARC,1)
        [pth,nm,ext] = fileparts(deblank(OutFigures.APARC(i,:)));
        [IMtype] = strread(nm,'%s','delimiter','-');
        SliceN = IMtype{end};
        [SliceName,SliceNumber] = strread(char(SliceN),'%s%s','delimiter','_');
        SliceName = lower(deblank(char(SliceName)));
        switch SliceName
            case 'sagital'
                conts = conts+1;
                SagFiles = strvcat(SagFiles,deblank(OutFigures.APARC(i,:)));
                Numbers(conts) = str2num(char(SliceNumber));
            case 'coronal'
                contc = contc+1;
                CorFiles = strvcat(CorFiles,deblank(OutFigures.APARC(i,:)));
                Numberc(conts) = str2num(char(SliceNumber));
            case 'axial'
                conta = conta+1;
                AxFiles = strvcat(AxFiles,deblank(OutFigures.APARC(i,:)));
                Numbera(conts) = str2num(char(SliceNumber));
        end
    end
    % --------------------- Sorting Slices -----------------------------------%
    [a,b] = sort(Numbers);
    SagFiles = SagFiles(b,:);
    Numbers = Numbers(b);
    [a,b] = sort(Numberc);
    CorFiles = CorFiles(b,:);
    Numberc = Numberc(b);
    [a,b] = sort(Numbera);
    AxFiles = AxFiles(b,:);
    Numbera = Numbera(b);
    % --------------------- Table Organization -------------------------------%
    for i = 1:size(SagFiles,1);
        % ----------- Creating Column Titles Figures -------------------------%
        fprintf(fid,'%s\n',['<tr >']);
        fprintf(fid,'%s\n',['<td height="10" align="center"><I> Axial Slice  ' num2str(Numbera(i)) ' vox </I></td>']);
        fprintf(fid,'%s\n',['<td height="10" align="center"><I> Sagital Slice  ' num2str(Numbers(i)) ' vox </I></td>']);
        fprintf(fid,'%s\n',['<td height="10" align="center"><I> Coronal Slice  ' num2str(Numberc(i)) ' vox </I></td>']);
        fprintf(fid,'%s\n',['</tr>']);
        % ----------- Inserting Figures --------------------------------------%
        fprintf(fid,'%s\n',['<tr >']);
        % Axial
        [pth, nm,ext] = fileparts(AxFiles(i,:));
        T = imread(deblank(AxFiles(i,:)));
        resx = 400;
        resy = round(resx/size(T,1)*size(T,2));
        Filename = [nm ext];
        fprintf(fid,'%s\n',['<td><div align="center"><a href="' Filename '" target="_blank"><img src="' Filename ...
            '" alt="Sagital Slice Not Found" width="' num2str(resy) '" height=' num2str(resx) ' border="0" /></a></div></td>']);
        % Sagital
        [pth, nm,ext] = fileparts(SagFiles(i,:));
        T = imread(deblank(SagFiles(i,:)));
        resx = 400;
        resy = round(resx/size(T,1)*size(T,2));
        Filename = [nm ext];
        fprintf(fid,'%s\n',['<td><div align="center"><a href="' Filename '" target="_blank"><img src="' Filename ...
            '" alt="Axial Slice Not Found" width="' num2str(resy) '" height=' num2str(resx) ' border="0" /></a></div></td>']);
        % Coronal
        [pth, nm,ext] = fileparts(CorFiles(i,:));
        Filename = [nm ext];
        T = imread(deblank(CorFiles(i,:)));
        resx = 400;
        resy = round(resx/size(T,1)*size(T,2));
        fprintf(fid,'%s\n',['<td><div align="center"><a href="' Filename '" target="_blank"><img src="' Filename ...
            '" alt="Coronal Slice Not Found" width="' num2str(resy) '" height=' num2str(resx) ' border="0" /></a></div></td>']);
        fprintf(fid,'%s\n',['</tr>']);
    end
    fprintf(fid,'%s\n',['</table>']);
    %% ============= End of Creating T1 + Aparc+Aseg Atlas  Table ========== %%
end
if ~isempty(OutFigures.APARCs)
    %% ========= Creating T1 + Creating T1 + Aparc+A2009s Atlas Table ====== %%
    fprintf(fid,'%s\n',['<h2> T1 + Aparc2009s+Aseg Atlas </h2>']);
    fprintf(fid,'%s\n',['<table width="1600" height="400" border="1">']);
    fprintf(fid,'%s\n',['<tr height="20" align="center">']);
    fprintf(fid,'%s\n',['<td><b> Axial Slices </b></td>']);
    fprintf(fid,'%s\n',['<td><b> Sagital Slices </b></td>']);
    fprintf(fid,'%s\n',['<td><b> Coronal Slices </b></td>']);
    fprintf(fid,'%s\n',['</tr>']);
    % --------------------- Detecting Slices ---------------------------------%
    SagFiles = '';conts = 0;
    CorFiles= '';contc = 0;
    AxFiles = '';conta = 0;
    for i = 1:size(OutFigures.APARCs,1)
        [pth,nm,ext] = fileparts(deblank(OutFigures.APARCs(i,:)));
        [IMtype] = strread(nm,'%s','delimiter','-');
        SliceN = IMtype{end};
        [SliceName,SliceNumber] = strread(char(SliceN),'%s%s','delimiter','_');
        SliceName = lower(deblank(char(SliceName)));
        switch SliceName
            case 'sagital'
                conts = conts+1;
                SagFiles = strvcat(SagFiles,deblank(OutFigures.APARCs(i,:)));
                Numbers(conts) = str2num(char(SliceNumber));
            case 'coronal'
                contc = contc+1;
                CorFiles = strvcat(CorFiles,deblank(OutFigures.APARCs(i,:)));
                Numberc(conts) = str2num(char(SliceNumber));
            case 'axial'
                conta = conta+1;
                AxFiles = strvcat(AxFiles,deblank(OutFigures.APARCs(i,:)));
                Numbera(conts) = str2num(char(SliceNumber));
        end
    end
    % --------------------- Sorting Slices -----------------------------------%
    [a,b] = sort(Numbers);
    SagFiles = SagFiles(b,:);
    Numbers = Numbers(b);
    [a,b] = sort(Numberc);
    CorFiles = CorFiles(b,:);
    Numberc = Numberc(b);
    [a,b] = sort(Numbera);
    AxFiles = AxFiles(b,:);
    Numbera = Numbera(b);
    % --------------------- Table Organization -------------------------------%
    for i = 1:size(SagFiles,1);
        % ----------- Creating Column Titles Figures -------------------------%
        fprintf(fid,'%s\n',['<tr >']);
        fprintf(fid,'%s\n',['<td height="10" align="center"><I> Axial Slice  ' num2str(Numbera(i)) ' vox </I></td>']);
        fprintf(fid,'%s\n',['<td height="10" align="center"><I> Sagital Slice  ' num2str(Numbers(i)) ' vox </I></td>']);
        fprintf(fid,'%s\n',['<td height="10" align="center"><I> Coronal Slice  ' num2str(Numberc(i)) ' vox </I></td>']);
        fprintf(fid,'%s\n',['</tr>']);
        % ----------- Inserting Figures --------------------------------------%
        fprintf(fid,'%s\n',['<tr >']);
        % Axial
        [pth, nm,ext] = fileparts(AxFiles(i,:));
        T = imread(deblank(AxFiles(i,:)));
        resx = 400;
        resy = round(resx/size(T,1)*size(T,2));
        Filename = [nm ext];
        fprintf(fid,'%s\n',['<td><div align="center"><a href="' Filename '" target="_blank"><img src="' Filename ...
            '" alt="Axial Slice Not Found" width="' num2str(resy) '" height=' num2str(resx) ' border="0" /></a></div></td>']);
        % Sagital
        [pth, nm,ext] = fileparts(SagFiles(i,:));
        T = imread(deblank(SagFiles(i,:)));
        resx = 400;
        resy = round(resx/size(T,1)*size(T,2));
        Filename = [nm ext];
        fprintf(fid,'%s\n',['<td><div align="center"><a href="' Filename '" target="_blank"><img src="' Filename ...
            '" alt="Sagital Slice Not Found" width="' num2str(resy) '" height=' num2str(resx) ' border="0" /></a></div></td>']);
        % Coronal
        [pth, nm,ext] = fileparts(CorFiles(i,:));
        Filename = [nm ext];
        T = imread(deblank(CorFiles(i,:)));
        resx = 400;
        resy = round(resx/size(T,1)*size(T,2));
        fprintf(fid,'%s\n',['<td><div align="center"><a href="' Filename '" target="_blank"><img src="' Filename ...
            '" alt="Coronal Slice Not Found" width="' num2str(resy) '" height=' num2str(resx) ' border="0" /></a></div></td>']);
        fprintf(fid,'%s\n',['</tr>']);
    end
    fprintf(fid,'%s\n',['</table>']);
    %% ============ End of Creating T1 + Aparc+A2009s Atlas Table ========== %%
end

if ~isempty(OutFigures.VPS)
    %% ============== Creating T1 + Surfaces Table ========================= %%
    fprintf(fid,'%s\n',['<h2> T1 + White and Pial Surfaces </h2>']);
    fprintf(fid,'%s\n',['<table width="1600" height="400" border="1">']);
    fprintf(fid,'%s\n',['<tr height="20" align="center">']);
    fprintf(fid,'%s\n',['<td><b> Axial Slices </b></td>']);
    fprintf(fid,'%s\n',['<td><b> Sagital Slices </b></td>']);
    fprintf(fid,'%s\n',['<td><b> Coronal Slices </b></td>']);
    fprintf(fid,'%s\n',['</tr>']);
    % --------------------- Detecting Slices ---------------------------------%
    SagFiles = '';conts = 0;
    CorFiles= '';contc = 0;
    AxFiles = '';conta = 0;
    for i = 1:size(OutFigures.VPS,1)
        [pth,nm,ext] = fileparts(deblank(OutFigures.VPS(i,:)));
        [IMtype] = strread(nm,'%s','delimiter','-');
        SliceN = IMtype{end};
        [SliceName,SliceNumber] = strread(char(SliceN),'%s%s','delimiter','_');
        SliceName = lower(deblank(char(SliceName)));
        switch SliceName
            case 'sagital'
                conts = conts+1;
                SagFiles = strvcat(SagFiles,deblank(OutFigures.VPS(i,:)));
                Numbers(conts) = str2num(char(SliceNumber));
            case 'coronal'
                contc = contc+1;
                CorFiles = strvcat(CorFiles,deblank(OutFigures.VPS(i,:)));
                Numberc(conts) = str2num(char(SliceNumber));
            case 'axial'
                conta = conta+1;
                AxFiles = strvcat(AxFiles,deblank(OutFigures.VPS(i,:)));
                Numbera(conts) = str2num(char(SliceNumber));
        end
    end
    % --------------------- Sorting Slices -----------------------------------%
    [a,b] = sort(Numbers);
    SagFiles = SagFiles(b,:);
    Numbers = Numbers(b);
    [a,b] = sort(Numberc);
    CorFiles = CorFiles(b,:);
    Numberc = Numberc(b);
    [a,b] = sort(Numbera);
    AxFiles = AxFiles(b,:);
    Numbera = Numbera(b);
    % --------------------- Table Organization -------------------------------%
    for i = 1:size(SagFiles,1);
        % ----------- Creating Column Titles Figures -------------------------%
        fprintf(fid,'%s\n',['<tr >']);
        fprintf(fid,'%s\n',['<td height="10" align="center"><I> Axial Slice  ' num2str(Numbera(i)) ' vox </I></td>']);
        fprintf(fid,'%s\n',['<td height="10" align="center"><I> Sagital Slice  ' num2str(Numbers(i)) ' vox </I></td>']);
        fprintf(fid,'%s\n',['<td height="10" align="center"><I> Coronal Slice  ' num2str(Numberc(i)) ' vox </I></td>']);
        fprintf(fid,'%s\n',['</tr>']);
        % ----------- Inserting Figures --------------------------------------%
        fprintf(fid,'%s\n',['<tr >']);
        % Axial
        [pth, nm,ext] = fileparts(AxFiles(i,:));
        T = imread(deblank(AxFiles(i,:)));
        resx = 400;
        resy = round(resx/size(T,1)*size(T,2));
        Filename = [nm ext];
        fprintf(fid,'%s\n',['<td><div align="center"><a href="' Filename '" target="_blank"><img src="' Filename ...
            '" alt="Axial Slice Not Found" width="' num2str(resy) '" height=' num2str(resx) ' border="0" /></a></div></td>']);
        % Sagital
        [pth, nm,ext] = fileparts(SagFiles(i,:));
        T = imread(deblank(SagFiles(i,:)));
        resx = 400;
        resy = round(resx/size(T,1)*size(T,2));
        Filename = [nm ext];
        fprintf(fid,'%s\n',['<td><div align="center"><a href="' Filename '" target="_blank"><img src="' Filename ...
            '" alt="Sagital Slice Not Found" width="' num2str(resy) '" height=' num2str(resx) ' border="0" /></a></div></td>']);
        % Coronal
        [pth, nm,ext] = fileparts(CorFiles(i,:));
        Filename = [nm ext];
        T = imread(deblank(CorFiles(i,:)));
        resx = 400;
        resy = round(resx/size(T,1)*size(T,2));
        fprintf(fid,'%s\n',['<td><div align="center"><a href="' Filename '" target="_blank"><img src="' Filename ...
            '" alt="Coronal Slice Not Found" width="' num2str(resy) '" height=' num2str(resx) ' border="0" /></a></div></td>']);
        fprintf(fid,'%s\n',['</tr>']);
    end
    fprintf(fid,'%s\n',['</table>']);
    %% ============ End of Creating T1 + Aparc+A2009s Atlas Table ========== %%
end

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

%% ============== Creating Surface + Cortical Thickness Table ============== %%
fprintf(fid,'%s\n',['<h2> Surfaces + Cortical Thickness</h2>']);
fprintf(fid,'%s\n',['<table width="1600" height="400" border="1">']);

%% ========================= Left Hemisphere =========================== %%
fprintf(fid,'%s\n',['<tr height="20" align="center">']);
fprintf(fid,'%s\n',['<td><b>  Left Pial Surface + Cortical Thickness </b></td>']);
fprintf(fid,'%s\n',['<td><b>  Left White Surface + Cortical Thickness </b></td>']);
fprintf(fid,'%s\n',['<td><b>  Left Inflated Surface + Cortical Thickness </b></td>']);
fprintf(fid,'%s\n',['</tr>']);
%% ========================== Lateral View ============================= %%
% --------------------- Table Organization -------------------------------%
fprintf(fid,'%s\n',['<tr >']);
fprintf(fid,'%s\n',['<td height="10" align="center"><I> Lateral View  </I></td>']);
fprintf(fid,'%s\n',['<td height="10" align="center"><I> Lateral View  </I></td>']);
fprintf(fid,'%s\n',['<td height="10" align="center"><I> Lateral View  </I></td>']);
fprintf(fid,'%s\n',['</tr>']);
% ----------- Inserting Figures --------------------------------------%
fprintf(fid,'%s\n',['<tr >']);

% 1. Inserting LH Pial + Cortical Thickness  Lateral View
SFile = [OutFigures.QcontrolDir filesep 'LH_pial_plus_thickness-Lateral_View.jpg'];
if exist(SFile,'file')
    [pth, nm,ext] = fileparts(SFile);
    T = imread(deblank(SFile));
    resx = 400;
    resy = round(resx/size(T,1)*size(T,2));
    Filename = [nm ext];
else
    Filename ='';
    resx = 400;
    resy = 400;
end
fprintf(fid,'%s\n',['<td><div align="center"><a href="' Filename '" target="_blank"><img src="' Filename ...
    '" alt="LH_pial_plus_thickness-Lateral_View Surface Not Found" width="' num2str(resy) '" height=' num2str(resx) ' border="0" /></a></div></td>']);

% 2. Inserting LH White + Cortical Thickness  Lateral View
SFile = [OutFigures.QcontrolDir filesep 'LH_white_plus_thickness-Lateral_View.jpg'];
if exist(SFile,'file')
    [pth, nm,ext] = fileparts(SFile);
    T = imread(deblank(SFile));
    resx = 400;
    resy = round(resx/size(T,1)*size(T,2));
    Filename = [nm ext];
else
    Filename ='';
    resx = 400;
    resy = 400;
end
fprintf(fid,'%s\n',['<td><div align="center"><a href="' Filename '" target="_blank"><img src="' Filename ...
    '" alt="LH_white_plus_thickness-Lateral_View Surface Not Found" width="' num2str(resy) '" height=' num2str(resx) ' border="0" /></a></div></td>']);

% 3. Inserting LH Inflated + Cortical Thickness  Lateral View
SFile = [OutFigures.QcontrolDir filesep 'LH_inflated_plus_thickness-Lateral_View.jpg'];
if exist(SFile,'file')
    [pth, nm,ext] = fileparts(SFile);
    T = imread(deblank(SFile));
    resx = 400;
    resy = round(resx/size(T,1)*size(T,2));
    Filename = [nm ext];
else
    Filename ='';
    resx = 400;
    resy = 400;
end
fprintf(fid,'%s\n',['<td><div align="center"><a href="' Filename '" target="_blank"><img src="' Filename ...
    '" alt="LH_inflated_plus_thickness-Lateral_View Surface Not Found" width="' num2str(resy) '" height=' num2str(resx) ' border="0" /></a></div></td>']);
%% ====================== End of Lateral View ========================== %%

%% =========================== Medial View ============================= %%
%--------------------- Table Organization -------------------------------%
fprintf(fid,'%s\n',['<tr >']);
fprintf(fid,'%s\n',['<td height="10" align="center"><I> Medial View  </I></td>']);
fprintf(fid,'%s\n',['<td height="10" align="center"><I> Medial View  </I></td>']);
fprintf(fid,'%s\n',['<td height="10" align="center"><I> Medial View  </I></td>']);
fprintf(fid,'%s\n',['</tr>']);
% ----------- Inserting Figures --------------------------------------%
fprintf(fid,'%s\n',['<tr >']);

% 1. Inserting LH Pial + Cortical Thickness  Medial View
SFile = [OutFigures.QcontrolDir filesep 'LH_pial_plus_thickness-Medial_View.jpg'];
if exist(SFile,'file')
    [pth, nm,ext] = fileparts(SFile);
    T = imread(deblank(SFile));
    resx = 400;
    resy = round(resx/size(T,1)*size(T,2));
    Filename = [nm ext];
else
    Filename ='';
    resx = 400;
    resy = 400;
end
fprintf(fid,'%s\n',['<td><div align="center"><a href="' Filename '" target="_blank"><img src="' Filename ...
    '" alt="LH_pial_plus_thickness-Medial_View Surface Not Found" width="' num2str(resy) '" height=' num2str(resx) ' border="0" /></a></div></td>']);

% 2. Inserting LH White + Cortical Thickness  Medial View
SFile = [OutFigures.QcontrolDir filesep 'LH_white_plus_thickness-Medial_View.jpg'];
if exist(SFile,'file')
    [pth, nm,ext] = fileparts(SFile);
    T = imread(deblank(SFile));
    resx = 400;
    resy = round(resx/size(T,1)*size(T,2));
    Filename = [nm ext];
else
    Filename ='';
    resx = 400;
    resy = 400;
end
fprintf(fid,'%s\n',['<td><div align="center"><a href="' Filename '" target="_blank"><img src="' Filename ...
    '" alt="LH_white_plus_thickness-Medial_View Surface Not Found" width="' num2str(resy) '" height=' num2str(resx) ' border="0" /></a></div></td>']);

% 3. Inserting LH Inflated + Cortical Thickness  Medial View
SFile = [OutFigures.QcontrolDir filesep 'LH_inflated_plus_thickness-Medial_View.jpg'];
if exist(SFile,'file')
    [pth, nm,ext] = fileparts(SFile);
    T = imread(deblank(SFile));
    resx = 400;
    resy = round(resx/size(T,1)*size(T,2));
    Filename = [nm ext];
else
    Filename ='';
    resx = 400;
    resy = 400;
end
fprintf(fid,'%s\n',['<td><div align="center"><a href="' Filename '" target="_blank"><img src="' Filename ...
    '" alt="LH_inflated_plus_thickness-Medial_View Surface Not Found" width="' num2str(resy) '" height=' num2str(resx) ' border="0" /></a></div></td>']);
%% ======================== End of Medial View ========================= %%
%% ====================== End of Left Hemisphere ======================= %%

%% ========================= Right Hemisphere =========================== %%
fprintf(fid,'%s\n',['<tr height="20" align="center">']);
fprintf(fid,'%s\n',['<td><b>  Right Pial Surface + Cortical Thickness </b></td>']);
fprintf(fid,'%s\n',['<td><b>  Right White Surface + Cortical Thickness </b></td>']);
fprintf(fid,'%s\n',['<td><b>  Right Inflated Surface + Cortical Thickness </b></td>']);
fprintf(fid,'%s\n',['</tr>']);
%% ========================== Lateral View ============================= %%
% --------------------- Table Organization -------------------------------%
fprintf(fid,'%s\n',['<tr >']);
fprintf(fid,'%s\n',['<td height="10" align="center"><I> Lateral View  </I></td>']);
fprintf(fid,'%s\n',['<td height="10" align="center"><I> Lateral View  </I></td>']);
fprintf(fid,'%s\n',['<td height="10" align="center"><I> Lateral View  </I></td>']);
fprintf(fid,'%s\n',['</tr>']);
% ----------- Inserting Figures --------------------------------------%
fprintf(fid,'%s\n',['<tr >']);

% 1. Inserting RH Pial + Cortical Thickness  Lateral View
SFile = [OutFigures.QcontrolDir filesep 'RH_pial_plus_thickness-Lateral_View.jpg'];
if exist(SFile,'file')
    [pth, nm,ext] = fileparts(SFile);
    T = imread(deblank(SFile));
    resx = 400;
    resy = round(resx/size(T,1)*size(T,2));
    Filename = [nm ext];
else
    Filename ='';
    resx = 400;
    resy = 400;
end
fprintf(fid,'%s\n',['<td><div align="center"><a href="' Filename '" target="_blank"><img src="' Filename ...
    '" alt="RH_pial_plus_thickness-Lateral_View Surface Not Found" width="' num2str(resy) '" height=' num2str(resx) ' border="0" /></a></div></td>']);

% 2. Inserting RH White + Cortical Thickness  Lateral View
SFile = [OutFigures.QcontrolDir filesep 'RH_white_plus_thickness-Lateral_View.jpg'];
if exist(SFile,'file')
    [pth, nm,ext] = fileparts(SFile);
    T = imread(deblank(SFile));
    resx = 400;
    resy = round(resx/size(T,1)*size(T,2));
    Filename = [nm ext];
else
    Filename ='';
    resx = 400;
    resy = 400;
end
fprintf(fid,'%s\n',['<td><div align="center"><a href="' Filename '" target="_blank"><img src="' Filename ...
    '" alt="RH_white_plus_thickness-Lateral_View Surface Not Found" width="' num2str(resy) '" height=' num2str(resx) ' border="0" /></a></div></td>']);

% 3. Inserting RH Inflated + Cortical Thickness  Lateral View
SFile = [OutFigures.QcontrolDir filesep 'RH_inflated_plus_thickness-Lateral_View.jpg'];
if exist(SFile,'file')
    [pth, nm,ext] = fileparts(SFile);
    T = imread(deblank(SFile));
    resx = 400;
    resy = round(resx/size(T,1)*size(T,2));
    Filename = [nm ext];
else
    Filename ='';
    resx = 400;
    resy = 400;
end
fprintf(fid,'%s\n',['<td><div align="center"><a href="' Filename '" target="_blank"><img src="' Filename ...
    '" alt="RH_inflated_plus_thickness-Lateral_View Surface Not Found" width="' num2str(resy) '" height=' num2str(resx) ' border="0" /></a></div></td>']);
%% ====================== End of Lateral View ========================== %%

%% =========================== Medial View ============================= %%
%--------------------- Table Organization -------------------------------%
fprintf(fid,'%s\n',['<tr >']);
fprintf(fid,'%s\n',['<td height="10" align="center"><I> Medial View  </I></td>']);
fprintf(fid,'%s\n',['<td height="10" align="center"><I> Medial View  </I></td>']);
fprintf(fid,'%s\n',['<td height="10" align="center"><I> Medial View  </I></td>']);
fprintf(fid,'%s\n',['</tr>']);
% ----------- Inserting Figures --------------------------------------%
fprintf(fid,'%s\n',['<tr >']);

% 1. Inserting RH Pial + Cortical Thickness  Medial View
SFile = [OutFigures.QcontrolDir filesep 'RH_pial_plus_thickness-Medial_View.jpg'];
if exist(SFile,'file')
    [pth, nm,ext] = fileparts(SFile);
    T = imread(deblank(SFile));
    resx = 400;
    resy = round(resx/size(T,1)*size(T,2));
    Filename = [nm ext];
else
    Filename ='';
    resx = 400;
    resy = 400;
end
fprintf(fid,'%s\n',['<td><div align="center"><a href="' Filename '" target="_blank"><img src="' Filename ...
    '" alt="RH_pial_plus_thickness-Medial_View Surface Not Found" width="' num2str(resy) '" height=' num2str(resx) ' border="0" /></a></div></td>']);

% 2. Inserting RH White + Cortical Thickness  Medial View
SFile = [OutFigures.QcontrolDir filesep 'RH_white_plus_thickness-Medial_View.jpg'];
if exist(SFile,'file')
    [pth, nm,ext] = fileparts(SFile);
    T = imread(deblank(SFile));
    resx = 400;
    resy = round(resx/size(T,1)*size(T,2));
    Filename = [nm ext];
else
    Filename ='';
    resx = 400;
    resy = 400;
end
fprintf(fid,'%s\n',['<td><div align="center"><a href="' Filename '" target="_blank"><img src="' Filename ...
    '" alt="RH_white_plus_thickness-Medial_View Surface Not Found" width="' num2str(resy) '" height=' num2str(resx) ' border="0" /></a></div></td>']);

% 3. Inserting RH Inflated + Cortical Thickness  Medial View
SFile = [OutFigures.QcontrolDir filesep 'RH_inflated_plus_thickness-Medial_View.jpg'];
if exist(SFile,'file')
    [pth, nm,ext] = fileparts(SFile);
    T = imread(deblank(SFile));
    resx = 400;
    resy = round(resx/size(T,1)*size(T,2));
    Filename = [nm ext];
else
    Filename ='';
    resx = 400;
    resy = 400;
end
fprintf(fid,'%s\n',['<td><div align="center"><a href="' Filename '" target="_blank"><img src="' Filename ...
    '" alt="RH_inflated_plus_thickness-Medial_View Surface Not Found" width="' num2str(resy) '" height=' num2str(resx) ' border="0" /></a></div></td>']);
fprintf(fid,'%s\n',['</tr>']);
fprintf(fid,'%s\n',['</table>']);
%% ======================== End of Medial View ========================= %%
%% ====================== End of Right Hemisphere ====================== %%

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%% ============== Creating Surface + Curvature Maps Table ============== %%
fprintf(fid,'%s\n',['<h2> Surfaces + Curvature Map</h2>']);
fprintf(fid,'%s\n',['<table width="1600" height="400" border="1">']);

%% ========================= Left Hemisphere =========================== %%
fprintf(fid,'%s\n',['<tr height="20" align="center">']);
fprintf(fid,'%s\n',['<td><b>  Left Pial Surface + Curvature Map </b></td>']);
fprintf(fid,'%s\n',['<td><b>  Left White Surface + Curvature Map </b></td>']);
fprintf(fid,'%s\n',['<td><b>  Left Inflated Surface + Curvature Map </b></td>']);
fprintf(fid,'%s\n',['</tr>']);
%% ========================== Lateral View ============================= %%
% --------------------- Table Organization -------------------------------%
fprintf(fid,'%s\n',['<tr >']);
fprintf(fid,'%s\n',['<td height="10" align="center"><I> Lateral View  </I></td>']);
fprintf(fid,'%s\n',['<td height="10" align="center"><I> Lateral View  </I></td>']);
fprintf(fid,'%s\n',['<td height="10" align="center"><I> Lateral View  </I></td>']);
fprintf(fid,'%s\n',['</tr>']);
% ----------- Inserting Figures --------------------------------------%
fprintf(fid,'%s\n',['<tr >']);

% 1. Inserting LH Pial + Curvature  Lateral View
SFile = [OutFigures.QcontrolDir filesep 'LH_pial_plus_curv-Lateral_View.jpg'];
if exist(SFile,'file')
    [pth, nm,ext] = fileparts(SFile);
    T = imread(deblank(SFile));
    resx = 400;
    resy = round(resx/size(T,1)*size(T,2));
    Filename = [nm ext];
else
    Filename ='';
    resx = 400;
    resy = 400;
end
fprintf(fid,'%s\n',['<td><div align="center"><a href="' Filename '" target="_blank"><img src="' Filename ...
    '" alt="LH_pial_plus_curv-Lateral_View Surface Not Found" width="' num2str(resy) '" height=' num2str(resx) ' border="0" /></a></div></td>']);

% 2. Inserting LH White + Curvature  Lateral View
SFile = [OutFigures.QcontrolDir filesep 'LH_white_plus_curv-Lateral_View.jpg'];
if exist(SFile,'file')
    [pth, nm,ext] = fileparts(SFile);
    T = imread(deblank(SFile));
    resx = 400;
    resy = round(resx/size(T,1)*size(T,2));
    Filename = [nm ext];
else
    Filename ='';
    resx = 400;
    resy = 400;
end
fprintf(fid,'%s\n',['<td><div align="center"><a href="' Filename '" target="_blank"><img src="' Filename ...
    '" alt="LH_white_plus_curv-Lateral_View Surface Not Found" width="' num2str(resy) '" height=' num2str(resx) ' border="0" /></a></div></td>']);

% 3. Inserting LH Inflated + Curvature  Lateral View
SFile = [OutFigures.QcontrolDir filesep 'LH_inflated_plus_curv-Lateral_View.jpg'];
if exist(SFile,'file')
    [pth, nm,ext] = fileparts(SFile);
    T = imread(deblank(SFile));
    resx = 400;
    resy = round(resx/size(T,1)*size(T,2));
    Filename = [nm ext];
else
    Filename ='';
    resx = 400;
    resy = 400;
end
fprintf(fid,'%s\n',['<td><div align="center"><a href="' Filename '" target="_blank"><img src="' Filename ...
    '" alt="LH_inflated_plus_curv-Lateral_View Surface Not Found" width="' num2str(resy) '" height=' num2str(resx) ' border="0" /></a></div></td>']);
%% ====================== End of Lateral View ========================== %%

%% =========================== Medial View ============================= %%
%--------------------- Table Organization -------------------------------%
fprintf(fid,'%s\n',['<tr >']);
fprintf(fid,'%s\n',['<td height="10" align="center"><I> Medial View  </I></td>']);
fprintf(fid,'%s\n',['<td height="10" align="center"><I> Medial View  </I></td>']);
fprintf(fid,'%s\n',['<td height="10" align="center"><I> Medial View  </I></td>']);
fprintf(fid,'%s\n',['</tr>']);
% ----------- Inserting Figures --------------------------------------%
fprintf(fid,'%s\n',['<tr >']);

% 1. Inserting LH Pial + Curvature  Medial View
SFile = [OutFigures.QcontrolDir filesep 'LH_pial_plus_curv-Medial_View.jpg'];
if exist(SFile,'file')
    [pth, nm,ext] = fileparts(SFile);
    T = imread(deblank(SFile));
    resx = 400;
    resy = round(resx/size(T,1)*size(T,2));
    Filename = [nm ext];
else
    Filename ='';
    resx = 400;
    resy = 400;
end
fprintf(fid,'%s\n',['<td><div align="center"><a href="' Filename '" target="_blank"><img src="' Filename ...
    '" alt="LH_pial_plus_curv-Medial_View Surface Not Found" width="' num2str(resy) '" height=' num2str(resx) ' border="0" /></a></div></td>']);

% 2. Inserting LH White + Curvature  Medial View
SFile = [OutFigures.QcontrolDir filesep 'LH_white_plus_curv-Medial_View.jpg'];
if exist(SFile,'file')
    [pth, nm,ext] = fileparts(SFile);
    T = imread(deblank(SFile));
    resx = 400;
    resy = round(resx/size(T,1)*size(T,2));
    Filename = [nm ext];
else
    Filename ='';
    resx = 400;
    resy = 400;
end
fprintf(fid,'%s\n',['<td><div align="center"><a href="' Filename '" target="_blank"><img src="' Filename ...
    '" alt="LH_white_plus_curv-Medial_View Surface Not Found" width="' num2str(resy) '" height=' num2str(resx) ' border="0" /></a></div></td>']);

% 3. Inserting LH Inflated + Curvature  Medial View
SFile = [OutFigures.QcontrolDir filesep 'LH_inflated_plus_curv-Medial_View.jpg'];
if exist(SFile,'file')
    [pth, nm,ext] = fileparts(SFile);
    T = imread(deblank(SFile));
    resx = 400;
    resy = round(resx/size(T,1)*size(T,2));
    Filename = [nm ext];
else
    Filename ='';
    resx = 400;
    resy = 400;
end
fprintf(fid,'%s\n',['<td><div align="center"><a href="' Filename '" target="_blank"><img src="' Filename ...
    '" alt="LH_inflated_plus_curv-Medial_View Surface Not Found" width="' num2str(resy) '" height=' num2str(resx) ' border="0" /></a></div></td>']);

%% ======================== End of Medial View ========================= %%
%% ====================== End of Left Hemisphere ======================= %%

%% ========================= Right Hemisphere =========================== %%
fprintf(fid,'%s\n',['<tr height="20" align="center">']);
fprintf(fid,'%s\n',['<td><b>  Right Pial Surface + Curvature Map </b></td>']);
fprintf(fid,'%s\n',['<td><b>  Right White Surface + Curvature Map </b></td>']);
fprintf(fid,'%s\n',['<td><b>  Right Inflated Surface + Curvature Map </b></td>']);
fprintf(fid,'%s\n',['</tr>']);
%% ========================== Lateral View ============================= %%
% --------------------- Table Organization -------------------------------%
fprintf(fid,'%s\n',['<tr >']);
fprintf(fid,'%s\n',['<td height="10" align="center"><I> Lateral View  </I></td>']);
fprintf(fid,'%s\n',['<td height="10" align="center"><I> Lateral View  </I></td>']);
fprintf(fid,'%s\n',['<td height="10" align="center"><I> Lateral View  </I></td>']);
fprintf(fid,'%s\n',['</tr>']);
% ----------- Inserting Figures --------------------------------------%
fprintf(fid,'%s\n',['<tr >']);

% 1. Inserting RH Pial + Curvature  Lateral View
SFile = [OutFigures.QcontrolDir filesep 'RH_pial_plus_curv-Lateral_View.jpg'];
if exist(SFile,'file')
    [pth, nm,ext] = fileparts(SFile);
    T = imread(deblank(SFile));
    resx = 400;
    resy = round(resx/size(T,1)*size(T,2));
    Filename = [nm ext];
else
    Filename ='';
    resx = 400;
    resy = 400;
end
fprintf(fid,'%s\n',['<td><div align="center"><a href="' Filename '" target="_blank"><img src="' Filename ...
    '" alt="RH_pial_plus_curv-Lateral_View Surface Not Found" width="' num2str(resy) '" height=' num2str(resx) ' border="0" /></a></div></td>']);

% 2. Inserting RH White + Curvature  Latral View
SFile = [OutFigures.QcontrolDir filesep 'RH_white_plus_curv-Lateral_View.jpg'];
if exist(SFile,'file')
    [pth, nm,ext] = fileparts(SFile);
    T = imread(deblank(SFile));
    resx = 400;
    resy = round(resx/size(T,1)*size(T,2));
    Filename = [nm ext];
else
    Filename ='';
    resx = 400;
    resy = 400;
end
fprintf(fid,'%s\n',['<td><div align="center"><a href="' Filename '" target="_blank"><img src="' Filename ...
    '" alt="RH_white_plus_curv-Lateral_View Surface Not Found" width="' num2str(resy) '" height=' num2str(resx) ' border="0" /></a></div></td>']);

% 3. Inserting RH Inflated + Curvature  Lateral View
SFile = [OutFigures.QcontrolDir filesep 'RH_inflated_plus_curv-Lateral_View.jpg'];
if exist(SFile,'file')
    [pth, nm,ext] = fileparts(SFile);
    T = imread(deblank(SFile));
    resx = 400;
    resy = round(resx/size(T,1)*size(T,2));
    Filename = [nm ext];
else
    Filename ='';
    resx = 400;
    resy = 400;
end
fprintf(fid,'%s\n',['<td><div align="center"><a href="' Filename '" target="_blank"><img src="' Filename ...
    '" alt="RH_inflated_plus_curv-Lateral_View Surface Not Found" width="' num2str(resy) '" height=' num2str(resx) ' border="0" /></a></div></td>']);
%% ====================== End of Lateral View ========================== %%

%% =========================== Medial View ============================= %%
%--------------------- Table Organization -------------------------------%
fprintf(fid,'%s\n',['<tr >']);
fprintf(fid,'%s\n',['<td height="10" align="center"><I> Medial View  </I></td>']);
fprintf(fid,'%s\n',['<td height="10" align="center"><I> Medial View  </I></td>']);
fprintf(fid,'%s\n',['<td height="10" align="center"><I> Medial View  </I></td>']);
fprintf(fid,'%s\n',['</tr>']);
% ----------- Inserting Figures --------------------------------------%
fprintf(fid,'%s\n',['<tr >']);

% 1. Inserting RH Pial + Curvature  Medial View
SFile = [OutFigures.QcontrolDir filesep 'RH_pial_plus_curv-Medial_View.jpg'];
if exist(SFile,'file')
    [pth, nm,ext] = fileparts(SFile);
    T = imread(deblank(SFile));
    resx = 400;
    resy = round(resx/size(T,1)*size(T,2));
    Filename = [nm ext];
else
    Filename ='';
    resx = 400;
    resy = 400;
end
fprintf(fid,'%s\n',['<td><div align="center"><a href="' Filename '" target="_blank"><img src="' Filename ...
    '" alt="RH_pial_plus_curv-Medial_View Surface Not Found" width="' num2str(resy) '" height=' num2str(resx) ' border="0" /></a></div></td>']);

% 2. Inserting RH White + Curvature  Medial View
SFile = [OutFigures.QcontrolDir filesep 'RH_white_plus_curv-Medial_View.jpg'];
if exist(SFile,'file')
    [pth, nm,ext] = fileparts(SFile);
    T = imread(deblank(SFile));
    resx = 400;
    resy = round(resx/size(T,1)*size(T,2));
    Filename = [nm ext];
else
    Filename ='';
    resx = 400;
    resy = 400;
end
fprintf(fid,'%s\n',['<td><div align="center"><a href="' Filename '" target="_blank"><img src="' Filename ...
    '" alt="RH_white_plus_curv-Medial_View Surface Not Found" width="' num2str(resy) '" height=' num2str(resx) ' border="0" /></a></div></td>']);

% 3. Inserting RH Inflated + Curvature  Medial View
SFile = [OutFigures.QcontrolDir filesep 'RH_inflated_plus_curv-Medial_View.jpg'];
if exist(SFile,'file')
    [pth, nm,ext] = fileparts(SFile);
    T = imread(deblank(SFile));
    resx = 400;
    resy = round(resx/size(T,1)*size(T,2));
    Filename = [nm ext];
else
    Filename ='';
    resx = 400;
    resy = 400;
end
fprintf(fid,'%s\n',['<td><div align="center"><a href="' Filename '" target="_blank"><img src="' Filename ...
    '" alt="RH_inflated_plus_curv-Medial_View Surface Not Found" width="' num2str(resy) '" height=' num2str(resx) ' border="0" /></a></div></td>']);
fprintf(fid,'%s\n',['</tr>']);
fprintf(fid,'%s\n',['</table>']);
%% ======================== End of Medial View ========================= %%
%% ====================== End of Right Hemisphere ====================== %%

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

%% ============== Creating Surface + APARC Atlas Table ============== %%
fprintf(fid,'%s\n',['<h2> Surfaces + Aparc Atlas</h2>']);
fprintf(fid,'%s\n',['<table width="1600" height="400" border="1">']);

%% ========================= Left Hemisphere =========================== %%
fprintf(fid,'%s\n',['<tr height="20" align="center">']);
fprintf(fid,'%s\n',['<td><b>  Left Pial Surface + Aparc Atlas </b></td>']);
fprintf(fid,'%s\n',['<td><b>  Left White Surface + Aparc Atlas </b></td>']);
fprintf(fid,'%s\n',['<td><b>  Left Inflated Surface + Aparc Atlas </b></td>']);
fprintf(fid,'%s\n',['</tr>']);
%% ========================== Lateral View ============================= %%
% --------------------- Table Organization -------------------------------%
fprintf(fid,'%s\n',['<tr >']);
fprintf(fid,'%s\n',['<td height="10" align="center"><I> Lateral View  </I></td>']);
fprintf(fid,'%s\n',['<td height="10" align="center"><I> Lateral View  </I></td>']);
fprintf(fid,'%s\n',['<td height="10" align="center"><I> Lateral View  </I></td>']);
fprintf(fid,'%s\n',['</tr>']);
% ----------- Inserting Figures --------------------------------------%
fprintf(fid,'%s\n',['<tr >']);

% 1. Inserting LH Pial + APARC Atlas  Lateral View
SFile = [OutFigures.QcontrolDir filesep 'LH_pial_plus_aparc-Lateral_View.jpg'];
if exist(SFile,'file')
    [pth, nm,ext] = fileparts(SFile);
    T = imread(deblank(SFile));
    resx = 400;
    resy = round(resx/size(T,1)*size(T,2));
    Filename = [nm ext];
else
    Filename ='';
    resx = 400;
    resy = 400;
end
fprintf(fid,'%s\n',['<td><div align="center"><a href="' Filename '" target="_blank"><img src="' Filename ...
    '" alt="LH_pial_plus_aparc-Lateral_View Surface Not Found" width="' num2str(resy) '" height=' num2str(resx) ' border="0" /></a></div></td>']);

% 2. Inserting LH White + APARC Atlas  Lateral View
SFile = [OutFigures.QcontrolDir filesep 'LH_white_plus_aparc-Lateral_View.jpg'];
if exist(SFile,'file')
    [pth, nm,ext] = fileparts(SFile);
    T = imread(deblank(SFile));
    resx = 400;
    resy = round(resx/size(T,1)*size(T,2));
    Filename = [nm ext];
else
    Filename ='';
    resx = 400;
    resy = 400;
end
fprintf(fid,'%s\n',['<td><div align="center"><a href="' Filename '" target="_blank"><img src="' Filename ...
    '" alt="LH_white_plus_aparc-Lateral_View Surface Not Found" width="' num2str(resy) '" height=' num2str(resx) ' border="0" /></a></div></td>']);

% 3. Inserting LH Inflated + APARC Atlas  Lateral View
SFile = [OutFigures.QcontrolDir filesep 'LH_inflated_plus_aparc-Lateral_View.jpg'];
if exist(SFile,'file')
    [pth, nm,ext] = fileparts(SFile);
    T = imread(deblank(SFile));
    resx = 400;
    resy = round(resx/size(T,1)*size(T,2));
    Filename = [nm ext];
else
    Filename ='';
    resx = 400;
    resy = 400;
end
fprintf(fid,'%s\n',['<td><div align="center"><a href="' Filename '" target="_blank"><img src="' Filename ...
    '" alt="LH_inflated_plus_aparc-Lateral_View Surface Not Found" width="' num2str(resy) '" height=' num2str(resx) ' border="0" /></a></div></td>']);
%% ====================== End of Lateral View ========================== %%

%% =========================== Medial View ============================= %%
%--------------------- Table Organization -------------------------------%
fprintf(fid,'%s\n',['<tr >']);
fprintf(fid,'%s\n',['<td height="10" align="center"><I> Medial View  </I></td>']);
fprintf(fid,'%s\n',['<td height="10" align="center"><I> Medial View  </I></td>']);
fprintf(fid,'%s\n',['<td height="10" align="center"><I> Medial View  </I></td>']);
fprintf(fid,'%s\n',['</tr>']);
% ----------- Inserting Figures --------------------------------------%
fprintf(fid,'%s\n',['<tr >']);

% 1. Inserting LH Pial + APARC Atlas  Medial View
SFile = [OutFigures.QcontrolDir filesep 'LH_pial_plus_aparc-Medial_View.jpg'];
if exist(SFile,'file')
    [pth, nm,ext] = fileparts(SFile);
    T = imread(deblank(SFile));
    resx = 400;
    resy = round(resx/size(T,1)*size(T,2));
    Filename = [nm ext];
else
    Filename ='';
    resx = 400;
    resy = 400;
end
fprintf(fid,'%s\n',['<td><div align="center"><a href="' Filename '" target="_blank"><img src="' Filename ...
    '" alt="LH_pial_plus_aparc-Medial_View Surface Not Found" width="' num2str(resy) '" height=' num2str(resx) ' border="0" /></a></div></td>']);

% 2. Inserting LH White + APARC Atlas  Medial View
SFile = [OutFigures.QcontrolDir filesep 'LH_white_plus_aparc-Medial_View.jpg'];
if exist(SFile,'file')
    [pth, nm,ext] = fileparts(SFile);
    T = imread(deblank(SFile));
    resx = 400;
    resy = round(resx/size(T,1)*size(T,2));
    Filename = [nm ext];
else
    Filename ='';
    resx = 400;
    resy = 400;
end
fprintf(fid,'%s\n',['<td><div align="center"><a href="' Filename '" target="_blank"><img src="' Filename ...
    '" alt="LH_white_plus_aparc-Medial_View Surface Not Found" width="' num2str(resy) '" height=' num2str(resx) ' border="0" /></a></div></td>']);

% 3. Inserting LH Inflated + APARC Atlas  Medial View
SFile = [OutFigures.QcontrolDir filesep 'LH_inflated_plus_aparc-Medial_View.jpg'];
if exist(SFile,'file')
    [pth, nm,ext] = fileparts(SFile);
    T = imread(deblank(SFile));
    resx = 400;
    resy = round(resx/size(T,1)*size(T,2));
    Filename = [nm ext];
else
    Filename ='';
    resx = 400;
    resy = 400;
end
fprintf(fid,'%s\n',['<td><div align="center"><a href="' Filename '" target="_blank"><img src="' Filename ...
    '" alt="LH_inflated_plus_aparc-Medial_View Surface Not Found" width="' num2str(resy) '" height=' num2str(resx) ' border="0" /></a></div></td>']);
%% ======================== End of Medial View ========================= %%
%% ====================== End of Left Hemisphere ======================= %%

%% ========================= Right Hemisphere =========================== %%
fprintf(fid,'%s\n',['<tr height="20" align="center">']);
fprintf(fid,'%s\n',['<td><b>  Right Pial Surface + Aparc Atlas </b></td>']);
fprintf(fid,'%s\n',['<td><b>  Right White Surface + Aparc Atlas </b></td>']);
fprintf(fid,'%s\n',['<td><b>  Right Inflated Surface + Aparc Atlas </b></td>']);
fprintf(fid,'%s\n',['</tr>']);
%% ========================== Lateral View ============================= %%
% --------------------- Table Organization -------------------------------%
fprintf(fid,'%s\n',['<tr >']);
fprintf(fid,'%s\n',['<td height="10" align="center"><I> Lateral View  </I></td>']);
fprintf(fid,'%s\n',['<td height="10" align="center"><I> Lateral View  </I></td>']);
fprintf(fid,'%s\n',['<td height="10" align="center"><I> Lateral View  </I></td>']);
fprintf(fid,'%s\n',['</tr>']);
% ----------- Inserting Figures --------------------------------------%
fprintf(fid,'%s\n',['<tr >']);

% 1. Inserting RH Pial + APARC Atlas   Lateral View
SFile = [OutFigures.QcontrolDir filesep 'RH_pial_plus_aparc-Lateral_View.jpg'];
if exist(SFile,'file')
    [pth, nm,ext] = fileparts(SFile);
    T = imread(deblank(SFile));
    resx = 400;
    resy = round(resx/size(T,1)*size(T,2));
    Filename = [nm ext];
else
    Filename ='';
    resx = 400;
    resy = 400;
end
fprintf(fid,'%s\n',['<td><div align="center"><a href="' Filename '" target="_blank"><img src="' Filename ...
    '" alt="RH_pial_plus_aparc-Lateral_View Surface Not Found" width="' num2str(resy) '" height=' num2str(resx) ' border="0" /></a></div></td>']);

% 2. Inserting RH White + APARC Atlas   Lateral View
SFile = [OutFigures.QcontrolDir filesep 'RH_white_plus_aparc-Lateral_View.jpg'];
if exist(SFile,'file')
    [pth, nm,ext] = fileparts(SFile);
    T = imread(deblank(SFile));
    resx = 400;
    resy = round(resx/size(T,1)*size(T,2));
    Filename = [nm ext];
else
    Filename ='';
    resx = 400;
    resy = 400;
end
fprintf(fid,'%s\n',['<td><div align="center"><a href="' Filename '" target="_blank"><img src="' Filename ...
    '" alt="RH_white_plus_aparc-Lateral_View Surface Not Found" width="' num2str(resy) '" height=' num2str(resx) ' border="0" /></a></div></td>']);

% 3. Inserting RH Inflated + APARC Atlas  Lateral View
SFile = [OutFigures.QcontrolDir filesep 'RH_inflated_plus_aparc-Lateral_View.jpg'];
if exist(SFile,'file')
    [pth, nm,ext] = fileparts(SFile);
    T = imread(deblank(SFile));
    resx = 400;
    resy = round(resx/size(T,1)*size(T,2));
    Filename = [nm ext];
else
    Filename ='';
    resx = 400;
    resy = 400;
end
fprintf(fid,'%s\n',['<td><div align="center"><a href="' Filename '" target="_blank"><img src="' Filename ...
    '" alt="RH_inflated_plus_aparc-Lateral_View Surface Not Found" width="' num2str(resy) '" height=' num2str(resx) ' border="0" /></a></div></td>']);
%% ====================== End of Lateral View ========================== %%

%% =========================== Medial View ============================= %%
%--------------------- Table Organization -------------------------------%
fprintf(fid,'%s\n',['<tr >']);
fprintf(fid,'%s\n',['<td height="10" align="center"><I> Medial View  </I></td>']);
fprintf(fid,'%s\n',['<td height="10" align="center"><I> Medial View  </I></td>']);
fprintf(fid,'%s\n',['<td height="10" align="center"><I> Medial View  </I></td>']);
fprintf(fid,'%s\n',['</tr>']);
% ----------- Inserting Figures --------------------------------------%
fprintf(fid,'%s\n',['<tr >']);

% 1. Inserting RH Pial + APARC Atlas  Medial View
SFile = [OutFigures.QcontrolDir filesep 'RH_pial_plus_aparc-Medial_View.jpg'];
if exist(SFile,'file')
    [pth, nm,ext] = fileparts(SFile);
    T = imread(deblank(SFile));
    resx = 400;
    resy = round(resx/size(T,1)*size(T,2));
    Filename = [nm ext];
else
    Filename ='';
    resx = 400;
    resy = 400;
end
fprintf(fid,'%s\n',['<td><div align="center"><a href="' Filename '" target="_blank"><img src="' Filename ...
    '" alt="RH_pial_plus_aparc-Medial_View Surface Not Found" width="' num2str(resy) '" height=' num2str(resx) ' border="0" /></a></div></td>']);

% 2. Inserting RH White + APARC Atlas  Medial View
SFile = [OutFigures.QcontrolDir filesep 'RH_white_plus_aparc-Medial_View.jpg'];
if exist(SFile,'file')
    [pth, nm,ext] = fileparts(SFile);
    T = imread(deblank(SFile));
    resx = 400;
    resy = round(resx/size(T,1)*size(T,2));
    Filename = [nm ext];
else
    Filename ='';
    resx = 400;
    resy = 400;
end
fprintf(fid,'%s\n',['<td><div align="center"><a href="' Filename '" target="_blank"><img src="' Filename ...
    '" alt="RH_white_plus_aparc-Medial_View Surface Not Found" width="' num2str(resy) '" height=' num2str(resx) ' border="0" /></a></div></td>']);

% 3. Inserting RH Inflated + APARC Atlas  Medial View
SFile = [OutFigures.QcontrolDir filesep 'RH_inflated_plus_aparc-Medial_View.jpg'];
if exist(SFile,'file')
    [pth, nm,ext] = fileparts(SFile);
    T = imread(deblank(SFile));
    resx = 400;
    resy = round(resx/size(T,1)*size(T,2));
    Filename = [nm ext];
else
    Filename ='';
    resx = 400;
    resy = 400;
end
fprintf(fid,'%s\n',['<td><div align="center"><a href="' Filename '" target="_blank"><img src="' Filename ...
    '" alt="RH_inflated_plus_aparc-Medial_View Surface Not Found" width="' num2str(resy) '" height=' num2str(resx) ' border="0" /></a></div></td>']);
fprintf(fid,'%s\n',['</tr>']);
fprintf(fid,'%s\n',['</table>']);
%% ======================== End of Medial View ========================= %%
%% ====================== End of Right Hemisphere ====================== %%

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

%% ============== Creating Surface + APARC2009s Atlas Table ============== %%
fprintf(fid,'%s\n',['<h2> Surfaces + APARC2009s Atlas</h2>']);
fprintf(fid,'%s\n',['<table width="1600" height="400" border="1">']);

%% ========================= Left Hemisphere =========================== %%
fprintf(fid,'%s\n',['<tr height="20" align="center">']);
fprintf(fid,'%s\n',['<td><b>  Left Pial Surface + APARC2009s Atlas </b></td>']);
fprintf(fid,'%s\n',['<td><b>  Left White Surface + APARC2009s Atlas </b></td>']);
fprintf(fid,'%s\n',['<td><b>  Left Inflated Surface + APARC2009s Atlas </b></td>']);
fprintf(fid,'%s\n',['</tr>']);
%% ========================== Lateral View ============================= %%
% --------------------- Table Organization -------------------------------%
fprintf(fid,'%s\n',['<tr >']);
fprintf(fid,'%s\n',['<td height="10" align="center"><I> Lateral View  </I></td>']);
fprintf(fid,'%s\n',['<td height="10" align="center"><I> Lateral View  </I></td>']);
fprintf(fid,'%s\n',['<td height="10" align="center"><I> Lateral View  </I></td>']);
fprintf(fid,'%s\n',['</tr>']);
% ----------- Inserting Figures --------------------------------------%
fprintf(fid,'%s\n',['<tr >']);

% 1. Inserting LH Pial + APARC2009s Atlas  Lateral View
SFile = [OutFigures.QcontrolDir filesep 'LH_pial_plus_a2009s-Lateral_View.jpg'];
if exist(SFile,'file')
    [pth, nm,ext] = fileparts(SFile);
    T = imread(deblank(SFile));
    resx = 400;
    resy = round(resx/size(T,1)*size(T,2));
    Filename = [nm ext];
else
    Filename ='';
    resx = 400;
    resy = 400;
end
fprintf(fid,'%s\n',['<td><div align="center"><a href="' Filename '" target="_blank"><img src="' Filename ...
    '" alt="LH_pial_plus_a2009s-Lateral_View Surface Not Found" width="' num2str(resy) '" height=' num2str(resx) ' border="0" /></a></div></td>']);

% 2. Inserting LH White + APARC2009s Atlas  Lateral View
SFile = [OutFigures.QcontrolDir filesep 'LH_white_plus_a2009s-Lateral_View.jpg'];
if exist(SFile,'file')
    [pth, nm,ext] = fileparts(SFile);
    T = imread(deblank(SFile));
    resx = 400;
    resy = round(resx/size(T,1)*size(T,2));
    Filename = [nm ext];
else
    Filename ='';
    resx = 400;
    resy = 400;
end
fprintf(fid,'%s\n',['<td><div align="center"><a href="' Filename '" target="_blank"><img src="' Filename ...
    '" alt="LH_white_plus_a2009s-Lateral_View Surface Not Found" width="' num2str(resy) '" height=' num2str(resx) ' border="0" /></a></div></td>']);

% 3. Inserting LH Inflated + APARC2009s Atlas  Lateral View
SFile = [OutFigures.QcontrolDir filesep 'LH_inflated_plus_a2009s-Lateral_View.jpg'];
if exist(SFile,'file')
    [pth, nm,ext] = fileparts(SFile);
    T = imread(deblank(SFile));
    resx = 400;
    resy = round(resx/size(T,1)*size(T,2));
    Filename = [nm ext];
else
    Filename ='';
    resx = 400;
    resy = 400;
end
fprintf(fid,'%s\n',['<td><div align="center"><a href="' Filename '" target="_blank"><img src="' Filename ...
    '" alt="LH_inflated_plus_a2009s-Lateral_View Surface Not Found" width="' num2str(resy) '" height=' num2str(resx) ' border="0" /></a></div></td>']);
%% ====================== End of Lateral View ========================== %%

%% =========================== Medial View ============================= %%
%--------------------- Table Organization -------------------------------%
fprintf(fid,'%s\n',['<tr >']);
fprintf(fid,'%s\n',['<td height="10" align="center"><I> Medial View  </I></td>']);
fprintf(fid,'%s\n',['<td height="10" align="center"><I> Medial View  </I></td>']);
fprintf(fid,'%s\n',['<td height="10" align="center"><I> Medial View  </I></td>']);
fprintf(fid,'%s\n',['</tr>']);
% ----------- Inserting Figures --------------------------------------%
fprintf(fid,'%s\n',['<tr >']);

% 1. Inserting LH Pial + APARC2009s Atlas  Medial View
SFile = [OutFigures.QcontrolDir filesep 'LH_pial_plus_a2009s-Medial_View.jpg'];
if exist(SFile,'file')
    [pth, nm,ext] = fileparts(SFile);
    T = imread(deblank(SFile));
    resx = 400;
    resy = round(resx/size(T,1)*size(T,2));
    Filename = [nm ext];
else
    Filename ='';
    resx = 400;
    resy = 400;
end
fprintf(fid,'%s\n',['<td><div align="center"><a href="' Filename '" target="_blank"><img src="' Filename ...
    '" alt="LH_pial_plus_a2009s-Medial_View Surface Not Found" width="' num2str(resy) '" height=' num2str(resx) ' border="0" /></a></div></td>']);

% 2. Inserting LH White + APARC2009s Atlas  Medial View
SFile = [OutFigures.QcontrolDir filesep 'LH_white_plus_a2009s-Medial_View.jpg'];
if exist(SFile,'file')
    [pth, nm,ext] = fileparts(SFile);
    T = imread(deblank(SFile));
    resx = 400;
    resy = round(resx/size(T,1)*size(T,2));
    Filename = [nm ext];
else
    Filename ='';
    resx = 400;
    resy = 400;
end
fprintf(fid,'%s\n',['<td><div align="center"><a href="' Filename '" target="_blank"><img src="' Filename ...
    '" alt="LH_white_plus_a2009s-Medial_View Surface Not Found" width="' num2str(resy) '" height=' num2str(resx) ' border="0" /></a></div></td>']);

% 3. Inserting LH Inflated + APARC2009s Atlas  Medial View
SFile = [OutFigures.QcontrolDir filesep 'LH_inflated_plus_a2009s-Medial_View.jpg'];
if exist(SFile,'file')
    [pth, nm,ext] = fileparts(SFile);
    T = imread(deblank(SFile));
    resx = 400;
    resy = round(resx/size(T,1)*size(T,2));
    Filename = [nm ext];
else
    Filename ='';
    resx = 400;
    resy = 400;
end
fprintf(fid,'%s\n',['<td><div align="center"><a href="' Filename '" target="_blank"><img src="' Filename ...
    '" alt="LH_inflated_plus_a2009s-Medial_View Surface Not Found" width="' num2str(resy) '" height=' num2str(resx) ' border="0" /></a></div></td>']);
%% ======================== End of Medial View ========================= %%
%% ====================== End of Left Hemisphere ======================= %%

%% ========================= Right Hemisphere =========================== %%
fprintf(fid,'%s\n',['<tr height="20" align="center">']);
fprintf(fid,'%s\n',['<td><b>  Right Pial Surface + APARC2009s Atlas </b></td>']);
fprintf(fid,'%s\n',['<td><b>  Right White Surface + APARC2009s Atlas </b></td>']);
fprintf(fid,'%s\n',['<td><b>  Right Inflated Surface + APARC2009s Atlas </b></td>']);
fprintf(fid,'%s\n',['</tr>']);
%% ========================== Lateral View ============================= %%
% --------------------- Table Organization -------------------------------%
fprintf(fid,'%s\n',['<tr >']);
fprintf(fid,'%s\n',['<td height="10" align="center"><I> Lateral View  </I></td>']);
fprintf(fid,'%s\n',['<td height="10" align="center"><I> Lateral View  </I></td>']);
fprintf(fid,'%s\n',['<td height="10" align="center"><I> Lateral View  </I></td>']);
fprintf(fid,'%s\n',['</tr>']);
% ----------- Inserting Figures --------------------------------------%
fprintf(fid,'%s\n',['<tr >']);

% 1. Inserting RH Pial + APARC2009s Atlas  Medial View
SFile = [OutFigures.QcontrolDir filesep 'RH_pial_plus_a2009s-Lateral_View.jpg'];
if exist(SFile,'file')
    [pth, nm,ext] = fileparts(SFile);
    T = imread(deblank(SFile));
    resx = 400;
    resy = round(resx/size(T,1)*size(T,2));
    Filename = [nm ext];
else
    Filename ='';
    resx = 400;
    resy = 400;
end
fprintf(fid,'%s\n',['<td><div align="center"><a href="' Filename '" target="_blank"><img src="' Filename ...
    '" alt="RH_pial_plus_a2009s-Lateral_View Surface Not Found" width="' num2str(resy) '" height=' num2str(resx) ' border="0" /></a></div></td>']);

% 2. Inserting RH White + APARC2009s Atlas  Medial View
SFile = [OutFigures.QcontrolDir filesep 'RH_white_plus_a2009s-Lateral_View.jpg'];
if exist(SFile,'file')
    [pth, nm,ext] = fileparts(SFile);
    T = imread(deblank(SFile));
    resx = 400;
    resy = round(resx/size(T,1)*size(T,2));
    Filename = [nm ext];
else
    Filename ='';
    resx = 400;
    resy = 400;
end
fprintf(fid,'%s\n',['<td><div align="center"><a href="' Filename '" target="_blank"><img src="' Filename ...
    '" alt="RH_white_plus_a2009s-Lateral_View Surface Not Found" width="' num2str(resy) '" height=' num2str(resx) ' border="0" /></a></div></td>']);

% 3. Inserting RH Inflated + APARC2009s Atlas  Medial View
SFile = [OutFigures.QcontrolDir filesep 'RH_inflated_plus_a2009s-Lateral_View.jpg'];
if exist(SFile,'file')
    [pth, nm,ext] = fileparts(SFile);
    T = imread(deblank(SFile));
    resx = 400;
    resy = round(resx/size(T,1)*size(T,2));
    Filename = [nm ext];
else
    Filename ='';
    resx = 400;
    resy = 400;
end
fprintf(fid,'%s\n',['<td><div align="center"><a href="' Filename '" target="_blank"><img src="' Filename ...
    '" alt="RH_inflated_plus_a2009s-Lateral_View Surface Not Found" width="' num2str(resy) '" height=' num2str(resx) ' border="0" /></a></div></td>']);
%% ====================== End of Lateral View ========================== %%

%% =========================== Medial View ============================= %%
%--------------------- Table Organization -------------------------------%
fprintf(fid,'%s\n',['<tr >']);
fprintf(fid,'%s\n',['<td height="10" align="center"><I> Medial View  </I></td>']);
fprintf(fid,'%s\n',['<td height="10" align="center"><I> Medial View  </I></td>']);
fprintf(fid,'%s\n',['<td height="10" align="center"><I> Medial View  </I></td>']);
fprintf(fid,'%s\n',['</tr>']);
% ----------- Inserting Figures --------------------------------------%
fprintf(fid,'%s\n',['<tr >']);

% 1. Inserting RH Pial + APARC2009s Atlas  Medial View
SFile = [OutFigures.QcontrolDir filesep 'RH_pial_plus_a2009s-Medial_View.jpg'];
if exist(SFile,'file')
    [pth, nm,ext] = fileparts(SFile);
    T = imread(deblank(SFile));
    resx = 400;
    resy = round(resx/size(T,1)*size(T,2));
    Filename = [nm ext];
else
    Filename ='';
    resx = 400;
    resy = 400;
end
fprintf(fid,'%s\n',['<td><div align="center"><a href="' Filename '" target="_blank"><img src="' Filename ...
    '" alt="RH_pial_plus_a2009s-Medial_View Surface Not Found" width="' num2str(resy) '" height=' num2str(resx) ' border="0" /></a></div></td>']);

% 2. Inserting RH White + APARC2009s Atlas  Medial View
SFile = [OutFigures.QcontrolDir filesep 'RH_white_plus_a2009s-Medial_View.jpg'];
if exist(SFile,'file')
    [pth, nm,ext] = fileparts(SFile);
    T = imread(deblank(SFile));
    resx = 400;
    resy = round(resx/size(T,1)*size(T,2));
    Filename = [nm ext];
else
    Filename ='';
    resx = 400;
    resy = 400;
end
fprintf(fid,'%s\n',['<td><div align="center"><a href="' Filename '" target="_blank"><img src="' Filename ...
    '" alt="RH_white_plus_a2009s-Medial_View Surface Not Found" width="' num2str(resy) '" height=' num2str(resx) ' border="0" /></a></div></td>']);

% 3. Inserting RH Inflated + APARC2009s Atlas  Medial View
SFile = [OutFigures.QcontrolDir filesep 'RH_inflated_plus_a2009s-Medial_View.jpg'];
if exist(SFile,'file')
    [pth, nm,ext] = fileparts(SFile);
    T = imread(deblank(SFile));
    resx = 400;
    resy = round(resx/size(T,1)*size(T,2));
    Filename = [nm ext];
else
    Filename ='';
    resx = 400;
    resy = 400;
end
fprintf(fid,'%s\n',['<td><div align="center"><a href="' Filename '" target="_blank"><img src="' Filename ...
    '" alt="RH_inflated_plus_a2009s-Medial_View Surface Not Found" width="' num2str(resy) '" height=' num2str(resx) ' border="0" /></a></div></td>']);
fprintf(fid,'%s\n',['</tr>']);
fprintf(fid,'%s\n',['</table>']);
%% ======================== End of Medial View ========================= %%
%% ====================== End of Right Hemisphere ====================== %%
fclose(fid);
return

