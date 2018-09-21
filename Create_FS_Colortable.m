function [colortable,labels, outColors] = Create_FS_Colortable(txt,Newnames);
%
% Syntax :
% [colortable,labels, outColors] = Create_FS_Colortable(txt);
%
% This script creates a Freesurfer colortable for a specified labels vector
% 
%
% Input Parameters:
%      txt             : Input Structures Labels
%     Newnames         : Structures Names
%
% Output Parameters:
%   colortable         : Freesurfer colortable in Matlab structure
%    labels            : Structures Labels 
%    outColors         : Generated Surface Colors
%
% See also: 
%__________________________________________________
% Authors: Yasser Aleman Gomez 
% LIM, HUGGM
% March 22th 2012
% Version $1.0


%============================ Main Program ===============================%
sts = unique(txt);
Nstruct = length(sts);
t = 1;
while t ~=0
    colors = randi([0 255],[Nstruct 3]);
    t = size(unique(colors,'rows'),1)-size(colors,1);
end
ctab = [colors colors(:,1)*0 colors(:,1)+colors(:,2)*2^8+colors(:,3)*2^16];
outColors = ones(length(txt),3)*255;
labels = txt*0;
for i = 1:Nstruct
    ind = find(txt==sts(i));
    if sts(i) == 0
        strnames{i,1} = 'unknown';
         ctab(i,:)  = [255 255 255 0 16777215];
         labels(ind) = 16777215;
%          ctab(i,:)  = [0 0 0 0 0];
%          labels(ind) = 0;
    else
        labels(ind) = ones(size(ind,1),1)*ctab(i,5);
        outColors(ind,:) = repmat(ctab(i,1:3),[length(ind) 1]);
        if nargin == 1
            strnames{i,1} = ['Region_' sprintf('%.5d',i)];
        else
            strnames{i,1} = deblank(Newnames(i,:));
        end
    end
end
colortable.numEntries = length(sts);
colortable.orig_tab = 'Custom Colortable';
colortable.struct_names = strnames;
colortable.table = ctab;
%============================ Main Program ===============================%
return