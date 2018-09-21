function [Surf] = Surf_Sep(Surf);
%
% Syntax :
% [Surf] = Surf_Corr(Surf);
%
% This function removes small isolated points ( or group of points).     
%
% Input Parameters:
%   Surf      : Individual Surfaces.
%
% Output Parameters:
%   Surf  : Corrected Surfaces.
%
% Related references:
%
%
% See also: Red_Surf Smooth_Surf Plot_Surf Surf_Comp Plot_oversurf 
% Exp_Surf
%
%__________________________________________________
% Authors: Yasser Aleman Gomez
% Neuroimaging Department
% Cuban Neuroscience Center
% December 1st 2007
% Version $1.0

Surf.Is = ones(size(Surf.SurfData.vertices,1),1);
strl = unique(Surf.Is);strl(strl==0) = [];
if ~isfield(Surf,'Tri')
    Npoints = size(Surf.SurfData.vertices,1);
    Nfaces = size(Surf.SurfData.faces,1);
    [Tri] = Vert_Neib(double(Surf.SurfData.faces),Npoints,Nfaces);
    Temp = sum(Tri);
    Tri(:,Temp==0) = [];
    Surf.Tri = Tri; clear Tri;
end
lab = Surf.Is;
set(0,'RecursionLimit',1000);
for i = 1:size(strl,1)
    lab = strl(i);
    [labid] = Recur_Corr(Surf,lab,zeros(size(Surf.Is)),1);
    Surf.Is = labid;
end
return;