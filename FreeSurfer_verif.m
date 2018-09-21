function [Outboolean] = FreeSurfer_verif(InputDir,Id,step);
%
% Syntax :
% [Outboolean] =FreeSurfer_verif(InputDir,Id,step);
%
% This function verifies the existence of all freesurfer outputs. It can be
% limited to predefined freesurfer stages.
%
% Input Parameters:
%  InputDir       : Freesurfer output directory
%  Id             : Subject ID
%  step            : Freesurfer stage (ie. autorecon1, autorecon2,
%  autorecon3);
%
% Output Parameter:
%   Outboolean    : Boolean variable to express if the freesurfer stage is
%                   completed and all the results are saved inside the freesurfer
%                   directory.
%
% See also: 
%__________________________________________________
% Authors: Yasser Aleman Gomez
%LIM
% February 1th 2012
% Version $1.0

if nargin < 3
    step = 'all';
end
%Autorecon1 No skull
filenamesa1 = strvcat([ 'mri' filesep 'orig' filesep '001.mgz'],[ 'mri' filesep 'rawavg.mgz'],[ 'mri' filesep 'orig.mgz'],[ 'mri' filesep 'nu.mgz'],[ 'mri' filesep 'T1.mgz'],[ 'mri' filesep 'transforms' filesep 'talairach.auto.xfm']);
%Autorecon1
filenamesa1 = strvcat([ 'mri' filesep 'orig' filesep '001.mgz'],[ 'mri' filesep 'rawavg.mgz'],[ 'mri' filesep 'orig.mgz'],[ 'mri' filesep 'nu.mgz'],[ 'mri' filesep 'T1.mgz'],[ 'mri' filesep 'brainmask.mgz'],[ 'mri' filesep 'brainmask.auto.mgz'],[ 'mri' filesep 'transforms' filesep 'talairach.auto.xfm']);

%Autorecon2
filenamesa2 = '';
filenamesa2 = strvcat(filenamesa2,[ 'mri' filesep 'transforms' filesep 'talairach.lta'],[ 'mri' filesep 'norm.mgz'],[ 'mri' filesep 'transforms' filesep 'talairach.m3z']);
filenamesa2 = strvcat(filenamesa2,[ 'mri' filesep 'transforms' filesep 'talairach.m3z.inv.x.mgz'],[ 'mri' filesep 'transforms' filesep 'talairach.m3z.inv.y.mgz'],[ 'mri' filesep 'transforms' filesep 'talairach.m3z.inv.z.mgz']);
filenamesa2 = strvcat(filenamesa2,[ 'mri' filesep 'nu_noneck.mgz'],[ 'mri' filesep 'transforms' filesep 'talairach_with_skull.lta'],[ 'mri' filesep 'aseg.auto_noCCseg.mgz'],[ 'mri' filesep 'aseg.auto.mgz'],[ 'mri' filesep 'aseg.mgz']);
filenamesa2 = strvcat(filenamesa2,[ 'mri' filesep 'brain.mgz'],[ 'mri' filesep 'brain.finalsurfs.mgz'],[ 'mri' filesep 'wm.seg.mgz'],[ 'mri' filesep 'wm.asegedit.mgz'],[ 'mri' filesep 'wm.mgz'],[ 'mri' filesep 'filled.mgz'],['scripts' filesep 'ponscc.cut.log']);
filenamesa2 = strvcat(filenamesa2,['surf' filesep 'lh.orig.nofix'],['surf' filesep 'lh.orig.nofix'],['surf' filesep 'lh.smoothwm.nofix']);
filenamesa2 = strvcat(filenamesa2,['surf' filesep 'rh.smoothwm.nofix'],['surf' filesep 'lh.inflated.nofix'],['surf' filesep 'rh.inflated.nofix'],['surf' filesep 'lh.qsphere.nofix'],['surf' filesep 'rh.qsphere.nofix']);
filenamesa2 = strvcat(filenamesa2,['surf' filesep 'lh.orig'],['surf' filesep 'rh.orig'],['surf' filesep 'lh.inflated'],['surf' filesep 'rh.inflated'],['surf' filesep 'lh.white'],['surf' filesep 'rh.white'],['surf' filesep 'lh.pial'],['surf' filesep 'rh.pial']);
filenamesa2 = strvcat(filenamesa2,['surf' filesep 'lh.thickness'],['surf' filesep 'rh.thickness'],['surf' filesep 'lh.curv'],['surf' filesep 'rh.curv'],['surf' filesep 'lh.area'],['surf' filesep 'rh.area'],['label' filesep 'lh.cortex.label'],['label' filesep 'rh.cortex.label']);
filenamesa2 = strvcat(filenamesa2,['surf' filesep 'lh.area.mid'],['surf' filesep 'rh.area.mid'],['surf' filesep 'lh.volume'],['surf' filesep 'rh.volume'],['surf' filesep 'lh.smoothwm'],['surf' filesep 'rh.smoothwm'],['surf' filesep 'lh.sulc'],['surf' filesep 'rh.sulc']);
filenamesa2 = strvcat(filenamesa2,['surf' filesep 'lh.inflated.H'],['surf' filesep 'rh.inflated.H'],['surf' filesep 'lh.inflated.K'],['surf' filesep 'rh.inflated.K']);


%Autorecon3
filenamesa3 = '';
filenamesa3 = strvcat(filenamesa3,['surf' filesep 'lh.sphere'],['surf' filesep 'rh.sphere'],['surf' filesep 'lh.sphere.reg'],['surf' filesep 'rh.sphere.reg'],['surf' filesep 'lh.jacobian_white'],['surf' filesep 'rh.jacobian_white'],['surf' filesep 'lh.avg_curv'],['surf' filesep 'rh.avg_curv']);
filenamesa3 = strvcat(filenamesa3,['label' filesep 'lh.aparc.annot'],['label' filesep 'rh.aparc.annot'],['stats' filesep 'lh.aparc.stats'],['stats' filesep 'rh.aparc.stats'],['label' filesep 'aparc.annot.ctab'],['label' filesep 'lh.aparc.a2009s.annot'],['label' filesep 'rh.aparc.a2009s.annot']);
filenamesa3 = strvcat(filenamesa3,['stats' filesep 'lh.aparc.a2009s.stats'],['stats' filesep 'rh.aparc.a2009s.stats'],['label' filesep 'aparc.annot.a2009s.ctab'],['stats' filesep 'aseg.stats'],[ 'mri' filesep 'lh.ribbon.mgz'],[ 'mri' filesep 'rh.ribbon.mgz'],[ 'mri' filesep 'aparc+aseg.mgz'],[ 'mri' filesep 'aparc.a2009s+aseg.mgz']);
filenamesa3 = strvcat(filenamesa3,[ 'mri' filesep 'wmparc.mgz'],['stats' filesep 'wmparc.stats']);

%TRACULA
%filentrabedpostx

switch step
    case 'autorecon1'
        filenames = filenamesa1;
    case 'autorecon2'
        filenames = strvcat(filenamesa1,filenamesa2);
    case 'autorecon2-cp'
        filenames = strvcat(filenamesa1,filenamesa2(1:13,:));
    case 'autorecon23-cp'
        filenames = strvcat(filenamesa2,filenamesa3);
    case 'autorecon23'
        filenames = strvcat(filenamesa2,filenamesa3);
    case 'autorecon2wm3'
        filenames = strvcat(filenamesa2,filenamesa3);
    case 'autorecon3'
        filenames = strvcat(filenamesa2,filenamesa3);
    case 'all'
        filenames = strvcat(filenamesa1,filenamesa2,filenamesa3);
    otherwise
        filenames = strvcat(filenamesa1,filenamesa2,filenamesa3);
end

A = zeros(1,size(filenames,1));
names = [repmat([InputDir filesep Id filesep],[size(filenames,1) 1]) filenames];
for j = 1:size(filenames,1)
    
    ind=  exist(deblank(names(j,:)));
    if ind>0
        A(j) = 1;
    end
end
Outboolean = prod(A);
return;