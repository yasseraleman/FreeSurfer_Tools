function [Surfout] = Divide_Surf(Surf, colortable);
%
% Syntax :
%   [Surfout] = Divide_Surf(Surf, colortable);
%
% This script divides a Surface matlab variable into individual surfaces
% according to a colortable.
%
% Input Parameters:
%        Surf             :  Matlab Surface Variable
%       colortable        :  Colortable
%
% Output Parameters:
%       Surfout           :  Diivided Surfaces
%
% Related references:
%
%
% See also: 
%
%
%__________________________________________________
% Authors: Yasser Aleman Gomez
% LIM, HUGGM
% August 14th 2012
% Version $1.0

txt = Surf.Is;
sts = unique(txt);
Ns = length(sts);
Nfaces = 0;
col = [1 0 0;0 1 0; 0 0 1; 1 1 0;0 1 1;1 0 1;1 0.5 0;0 0.5 1;1 0 0.5;0.5 1 0; 0.25 0.5 1;0.6 0.3 0.52;0.7 0.5 0.9];
Ncolor = size(col,1);
re = floor(length(sts)/Ncolor); col = repmat(col,[re+1 1]);
for i = 1:Ns
    ind = find(txt==sts(i));
    if nargin >1
        indc = find(colortable.table(:,5)==sts(i));
        if isempty(indc)
            Matname = 'unknown_structure';
            Color =   [1 1 1];      % Color
        else
            Matname = char(colortable.struct_names{indc});
            Color =  colortable.table(indc,1:3)/255;       % Color
        end
    else
        Matname = sprintf('%.5d',i);
        Color = col(i,:);
    end
    Surft = Surf;
    Faces = Surf.SurfData.faces;
    a = ismember(Surf.SurfData.faces,ind);
    ind2 = find(sum(a')==0);
    Surft.SurfData.faces(ind2,:) = [];
    [Surft] = Reorg_Surf(Surft);
    Surft.color = Color;
    Surft.SurfData.FaceVertexCData = repmat(Surft.color, [size(Surft.SurfData.vertices,1) 1]);
    Surfout(i) = Surft;
end
return