function HTML_File = Create_Project_FS_HTML(ProjectDir,IdFile);
%
% Syntax :
%   HTML_File = Create_Project_FS_HTML(FreeDir);
%
% This script creates a HTML file as quality control file. The HTML shows
% results for different processing steps.
%
% Input Parameters:
%      FreeDir            : FreeSurfer Subjects Directory
%       Id                : Subject Id
%
%
% Output Parameters:
%     HTML_File          :  Output HTML Filename
%
% Related references:
%
%
% See also: Create_FreeSurfer_Subject_HTML FreeSurfer_Qcontrol
% Surface_plus_overlay Image_plus_surface Image_plus_overlay
%
%
%__________________________________________________
% Authors: Yasser Aleman Gomez
% LIM, HUGGM
% August 14th 2012
% Version $1.0
% ProjectDir = '/media/Data/PROCESSING_RESULTS/PEPS';
FreeDir = [ProjectDir filesep '5-freesurfer_processing'];
if nargin <2
    a = dir([FreeDir filesep '*']);
    [namet{1:size(a,1)}] = deal(a.name);namet = char(namet);
    [direc{1:size(a,1)}] = deal(a.isdir);direc = cell2mat(direc);
    SubjIds = namet(direc,:);
    SubjIds = SubjIds(3:end,:);
    %% =============== Verifying if they are freesurfer results ======== %%
    a = zeros(size(SubjIds,1),1);
    for i = 1:size(SubjIds,1)
        if exist([FreeDir filesep deblank(SubjIds(i,:)) filesep 'mri' ],'dir');
            a(i) = 1;
        end
    end
    SubjIds = SubjIds(find(a),:);
else
    SubjIds = char(textread(IdFile,'%s'));
end
%% ============End of Verifying if they are freesurfer results ========= %%

%% ====================== Creating Individual HTML ===================== %%
% SubjIds = repmat(SubjIds(1,:),[30 1]);
indd = 0;
Ns = size(SubjIds,1);
cont = 0;
for i = 1:Ns
    Id = deblank(SubjIds(i,:));
    if ~exist([FreeDir filesep Id filesep 'qcontrol' filesep Id '-FS_index.html'],'file')
        try
            HTML_File = Create_FreeSurfer_Subject_HTML(FreeDir,Id);
        catch
            cont = cont + 1;
            indd(cont) = i;
        end
    end
end
%% =================== End of Creating Individual HTML ================= %%
if indd~=0
    SubjIds(indd,:) = [];
end
%% =================== Creating Project HTML =========================== %%
[pth,PName,ext] = fileparts(ProjectDir);
HTML_Filep = [FreeDir filesep PName '-FS_index.html'];
fido  = fopen(HTML_Filep,'wt');
fprintf(fido,'%s\n','<!DOCTYPE html PUBLIC "-//W3C//DTD XHTML 1.0 Transitional//EN" "http://www.w3.org/TR/xhtml1/DTD/xhtml1-transitional.dtd">');
fprintf(fido,'%s\n','<html xmlns="http://www.w3.org/1999/xhtml">');
fprintf(fido,'%s\n','<head>');
fprintf(fido,'%s\n','<meta http-equiv="Content-Type" content="text/html; charset=iso-8859-1" />');
fprintf(fido,'%s\n',['<title>HTML File for Project ' PName ' </title>']);
fprintf(fido,'%s\n','<link rel="stylesheet" type="text/css" href="main_style.css" />');
fprintf(fido,'%s\n','</head>');

% General information
fprintf(fido,'%s\n','<body bgcolor="#FFFFCC"> </body>'); % Background Color
fprintf(fido,'%s\n',['<h1> Project Name: ' PName ' </h1>']);
fprintf(fido,'%s\n',['<h2> Subjects Ids </h2>']);
Ns = size(SubjIds,1);
fprintf(fido,'%s\n',['<table width="800" height="40" border="0">']);
for i = 1:Ns
    Id = deblank(SubjIds(i,:));
    HTML_File = [FreeDir filesep Id filesep 'qcontrol' filesep Id '-FS_index.html'];
    fprintf(fido,'%s\n','<tr >');
    fprintf(fido,'%s\n',['<td><div align="center"><strong> <a href="./' Id filesep 'qcontrol' filesep Id '-FS_index.html " target="_blank" class="inlined">' Id ' </a> </strong></p></div></td>']);
%     fprintf(fido,'%s\n',['<td><div align="center"><img src="./' Id filesep 'qcontrol' filesep  'LH_inflated_plus_aparc-Lateral_View.jpg" alt="No" width="' num2str(53) '" height=' num2str(40) ' border="0" /></div></td>']);
%     fprintf(fido,'%s\n',['<td><div align="center"><img src="./' Id filesep 'qcontrol' filesep  'RH_inflated_plus_aparc-Lateral_View.jpg" alt="No" width="' num2str(53) '" height=' num2str(40) ' border="0" /></div></td>']);
%     fprintf(fido,'%s\n',['<td><div align="center"><img src="./' Id filesep 'qcontrol' filesep  'LH_inflated_plus_a2009s-Lateral_View.jpg" alt="No" width="' num2str(53) '" height=' num2str(40) ' border="0" /></div></td>']);
%     fprintf(fido,'%s\n',['<td><div align="center"><img src="./' Id filesep 'qcontrol' filesep  'RH_inflated_plus_a2009s-Lateral_View.jpg" alt="No" width="' num2str(53) '" height=' num2str(40) ' border="0" /></div></td>']);
%     fprintf(fido,'%s\n',['<td><div align="center"><img src="./' Id filesep 'qcontrol' filesep  'LH_inflated_plus_thickness-Lateral_View.jpg" alt="No" width="' num2str(53) '" height=' num2str(40) ' border="0" /></div></td>']);
%     fprintf(fido,'%s\n',['<td><div align="center"><img src="./' Id filesep 'qcontrol' filesep  'RH_inflated_plus_thickness-Lateral_View.jpg" alt="No" width="' num2str(53) '" height=' num2str(40) ' border="0" /></div></td>']);
%     fprintf(fido,'%s\n',['<td><div align="center"><img src="./' Id filesep 'qcontrol' filesep  'LH_inflated_plus_curv-Lateral_View.jpg" alt="No" width="' num2str(53) '" height=' num2str(40) ' border="0" /></div></td>']);
%     fprintf(fido,'%s\n',['<td><div align="center"><img src="./' Id filesep 'qcontrol' filesep  'RH_inflated_plus_curv-Lateral_View.jpg" alt="No" width="' num2str(53) '" height=' num2str(40) ' border="0" /></div></td>']);
    fprintf(fido,'%s\n',['</tr>']);
    %         fprintf(fido,'%s\n',['<p><strong> <a href="./' Id filesep 'qcontrol' filesep Id '-FS_index.html " target="_blank" class="inlined">' Id ' </a> </strong></p>']);
end
fprintf(fido,'%s\n',['</table>']);
%% =================== End of Creating Project HTML ==================== %%
fclose(fido);
return
