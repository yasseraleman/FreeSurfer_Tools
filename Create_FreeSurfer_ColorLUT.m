function OutFile = Create_FreeSurfer_ColorLUT(aparcId,OutDir);
%
% Syntax :
%  OutFile = Create_FreeSurfer_ColorLUT(aparcId,OutDir);
%
% This function creates a ColorLUT table for new freesurfer parcellations.
%
% Input Parameters:
%   aparcId     : Atlas type.
%   OutDir      : Output directory
%
% Output Parameters:
%   OutFile       : Output colorLUT
%
%
% Related references:
%
%
% See also: Create_Scales_volumetric_mgz Create_Scale_Aparc
%
%
%__________________________________________________
% Authors: Yasser Aleman Gomez
% LIM, HUGGM
% March 20th 2013
% Version $1.0

%% =============== Reading Colors for Subcortical Parcellation ========== %
aparcId = 'aparc_scale05_Nzones0410';
txtfile = '/media/COSAS/scripts/FreeSurfer_myTools/Subcort_FreeSurferColorLUT.txt';
OutDir = '/media/COSAS/Test/freesurfer/fsaverage/label';
fid = fopen(txtfile);
cont = 0;lines = '';
while 1
    cont = cont + 1;
    line = fgetl(fid);
    if ~ischar(line),   break,   end
    lines = strvcat(lines,line);
end
fclose(fid);
%% ========== End of Reading Colors for Subcortical Parcellation ======== %

%% =============== Creating Colors for Cortical Parcellation ============ %
ctabfile = ['/media/COSAS/Test/freesurfer/fsaverage/label' filesep aparcId '_ColorFile.txt'];
[stid,stname,r,g,b,t] = textread(ctabfile,'%u %s %u %u %u %u');
ctabl = [1000+stid r g b r*0 ];
lNames = [repmat('ctx-lh-',[size(ctabl,1) 1]) char(stname)];
ctabr = [2000+stid r g b r*0 ];
rNames = [repmat('ctx-rh-',[size(ctabr,1) 1]) char(stname)];
r = num2str(r);
g = num2str(g);
b = num2str(b);
ids = num2str(ctabl(:,1));
Total = [ids repmat(' ',[size(ids,1) 1]) lNames repmat('             ',[size(ids,1) 1]) r repmat(' ',[size(ids,1) 1]) g repmat(' ',[size(ids,1) 1]) b repmat(' ',[size(ids,1) 1]) num2str(zeros(size(b,1),1))];
lines = strvcat(lines,'# Left Hemisphere. Cortical Structures',Total);
ids = num2str(ctabr(:,1));
Total = [ids repmat(' ',[size(ids,1) 1]) rNames repmat('             ',[size(ids,1) 1]) r repmat(' ',[size(ids,1) 1]) g repmat(' ',[size(ids,1) 1]) b repmat(' ',[size(ids,1) 1]) num2str(zeros(size(b,1),1))];
lines = strvcat(lines,'# Right Hemisphere. Cortical Structures',Total);

%% =============== End of Creating Colors for Cortical Parcellationn ==== %
lines = strvcat(lines,'# White Matter','5001    Left-UnsegmentedWhiteMatter                20  30  40  0','5002    Right-UnsegmentedWhiteMatter               20  30  40  0');

%% ==================== Saving the New Color File ======================= %
OutFile = [OutDir filesep aparcId '_FreeSurferColorLUT.txt'];
fid = fopen(OutFile,'wt');
Ns = size(lines,1);
for i = 1:Ns
    fprintf(fid, '%s\n',deblank(lines(i,:)));
end
fclose(fid);
%% ============== End of Saving the New Color File ====================== %
return