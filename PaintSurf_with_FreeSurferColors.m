function Surf = PaintSurf_with_FreeSurferColors(Surf);
%
% Syntax :
% Surf = PaintSurf_with_FreeSurferColors(Surf);
%
% This scripts colors the surface using the freesurfer color map
%
% Input Parameters:
%   Surf      : Surface Struct
%
% Output Parameters:
%   Surf      : Surface Struct
%
% Related references:
%
%
% See also: Pipeline_Manage_perl
%__________________________________________________
% Authors: Yasser Aleman Gomez
% LIM
% September 18th 2012
% Version $1.0

sts = unique(Surf.Is);
sts(sts == 0) = [];
[GMcodes,Names,Colors] = Brain_GM_codes('aparc+aseg');
ctab = Colors(:,1)+Colors(:,2)*2^8+Colors(:,3)*2^16;
Surf.SurfData.FaceVertexCData = ones(length(Surf.Is),3);
for i = 1:length(sts)
    indl = find(GMcodes(:) == sts(i)|ctab == sts(i));
    if ~isempty(indl)
        inds = find(Surf.Is == sts(i));
        Surf.SurfData.FaceVertexCData(inds,:) = repmat(Colors(indl(1),:)/255, [length(inds) 1]);
    end
end
return