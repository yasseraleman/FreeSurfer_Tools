function [Surfml, Surfmr, MidSurfFiles] = Create_Mid_FreeSurfer_Surface(FreeSDir, SubjId, sa);
%
% Syntax :
%  MidSurfFiles = Create_Mid_FreeSurfer_Surface(FreeSDir, SubjId);
%
% This function generates mid surfaces from freesurfer pial and white
% surfaces.
%
% Input Parameters:
%     FreeSDir            : Freesurfer Subjects Directory
%       SubjId            : Subject Id
%        sa               : Boolean Variable to save or not the Mid Surface
%
% Output Parameters:
%     MidSurfFiles          : FreeSurfer Mid Surfaces
%
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
MidSurfFiles = '';
Surfml = '';
Surfmr = '';
if nargin < 3
    sa = 0;
end
%% ================ Reading Talairach File ============================== %
talfile = [FreeSDir filesep SubjId filesep 'mri' filesep 'transforms' filesep 'talairach.lta'];
cras1 = textread(talfile,'%s',5,'headerlines',20);
cras = char(cras1);
cras = [str2num(cras(3,:))  str2num(cras(4,:)) str2num(cras(5,:))];
%% ================ End of Reading Talairach File ======================= %


%% ================ Reading Surfaces ==================================== %
[OutFiles, SurfF] = Exp_Surf([FreeSDir filesep SubjId filesep 'surf' filesep 'lh.pial'], '0', '','', 'imp','n');
Surfpl= SurfF{1};
Surfpl.Name = 'LH.Pial';
Surfpl.SurfData.vertices = Surfpl.SurfData.vertices + repmat(cras,[size(Surfpl.SurfData.vertices,1) 1]);


[OutFiles, SurfF] = Exp_Surf([FreeSDir filesep SubjId filesep 'surf' filesep 'lh.white'], '0', '','', 'imp','n');
Surfwl= SurfF{1};
Surfwl.Name = 'LH.WHITE';
Surfwl.SurfData.vertices = Surfwl.SurfData.vertices + repmat(cras,[size(Surfwl.SurfData.vertices,1) 1]);


[OutFiles, SurfF] = Exp_Surf([FreeSDir filesep SubjId filesep 'surf' filesep 'rh.pial'], '0', '','', 'imp','n');
Surfpr= SurfF{1};
Surfpr.Name = 'RH.Pial';
Surfpr.SurfData.vertices = Surfpr.SurfData.vertices + repmat(cras,[size(Surfpr.SurfData.vertices,1) 1]);


[OutFiles, SurfF] = Exp_Surf([FreeSDir filesep SubjId filesep 'surf' filesep 'rh.white'], '0', '','', 'imp','n');
Surfwr= SurfF{1};
Surfwr.Name = 'RH.White';
Surfwr.SurfData.vertices = Surfwr.SurfData.vertices + repmat(cras,[size(Surfwr.SurfData.vertices,1) 1]);
%% ================ End of Reading Surfaces ============================= %

%% ==================== Creating Mid Surface ============================ %
% ----------------------- Left Hemisphere ------------------------------- %
xw = Surfwl.SurfData.vertices(:,1); yw = Surfwl.SurfData.vertices(:,2); zw = Surfwl.SurfData.vertices(:,3);
xp = Surfpl.SurfData.vertices(:,1); yp = Surfpl.SurfData.vertices(:,2); zp = Surfpl.SurfData.vertices(:,3);
norms = sqrt((xp-xw).^2+(yp-yw).^2+(zp-zw).^2);
u = xp*0;
v = yp*0;
w = zp*0;
ind = find(norms ~= 0);
u(ind) = (xp(ind)-xw(ind))./norms(ind); v(ind) = (yp(ind)-yw(ind))./norms(ind); w(ind) = (zp(ind)-zw(ind))./norms(ind);
vf =ones(size(norms));
[Ca] = [xw yw zw]+[repmat(norms.*vf/2,[ 1 3]).*[u v w]];
Surfml = Surfpl;
Surfml.SurfData.vertices = Ca;
OutFile = [FreeSDir filesep SubjId filesep 'surf' filesep 'lh.mid'];
if sa
    OutFile = save_free_surf(Surfml,OutFile);
end

% ----------------------- Right Hemisphere ------------------------------- %
xw = Surfwr.SurfData.vertices(:,1); yw = Surfwr.SurfData.vertices(:,2); zw = Surfwr.SurfData.vertices(:,3);
xp = Surfpr.SurfData.vertices(:,1); yp = Surfpr.SurfData.vertices(:,2); zp = Surfpr.SurfData.vertices(:,3);
norms = sqrt((xp-xw).^2+(yp-yw).^2+(zp-zw).^2);
u = xp*0;
v = yp*0;
w = zp*0;
ind = find(norms ~= 0);
u(ind) = (xp(ind)-xw(ind))./norms(ind); v(ind) = (yp(ind)-yw(ind))./norms(ind); w(ind) = (zp(ind)-zw(ind))./norms(ind);
vf =ones(size(norms));
[Ca] = [xw yw zw]+[repmat(norms.*vf/2,[ 1 3]).*[u v w]];
Surfmr = Surfpr;
Surfmr.SurfData.vertices = Ca;
OutFile = [FreeSDir filesep SubjId filesep 'surf' filesep 'rh.mid'];
if sa
    OutFile = save_free_surf(Surfmr,OutFile);
    MidSurfFiles = strvcat([FreeSDir filesep SubjId filesep 'surf' filesep 'lh.mid'],[FreeSDir filesep SubjId filesep 'surf' filesep 'rh.mid'])
end

%% ==================== End of Creating Mid Surface ===================== %
return;