function OutAnnot = save_annotfiles(labelsd,OutAnnot,colortable);
%
% Syntax :
% OutAnnot = save_annotfiles(labelsd,OutAnnot,colortable);
%
% This script creates an annotation file for cortical parcelation.
% 
%
% Input Parameters:
%   labelsd            : Structures Labels
%   OutAnnot           : New Annotation filename
%   colortable         : Freesurfer colortable in Matlab structure
%
% Output Parameters:
%   OutAnnot             : New Annotation filename
%
% See also: 
%__________________________________________________
% Authors: Yasser Aleman Gomez 
% LIM, HUGGM
% March 22th 2012
% Version $1.0

%=====================Checking Input Parameters===========================%
sts = unique(labelsd);
Nstruct = length(sts);
t = 1;
if nargin<3
    [colortable,labels] = Create_FS_Colortable(labelsd);
    ctab = colortable.table;
else
    labels = labelsd;
    ctab = colortable.table;
end
%=========================================================================%
%=====================     Main Program    ===============================%
Np = length(labels);
fid = fopen(OutAnnot,'wb','b');
fwrite(fid, Np, 'int');
Avalues(2:2:Np*2) = labels;
Avalues(1:2:Np*2) = [0:Np-1]';
fwrite(fid, Avalues, 'int');
fwrite(fid, 1, 'int');
fwrite(fid, -2, 'int');
fwrite(fid,length(unique(labels)),'int');
fwrite(fid, length(OutAnnot), 'int');
fwrite(fid, OutAnnot, '*char');
fwrite(fid,length(unique(labels)),'int');
for i = 1:Nstruct
    fwrite(fid, i-1, 'int');
    name = deblank(colortable.struct_names{i});
    len = length(name);
    fwrite(fid, len+1, 'int');
    fwrite(fid, [name ' '], 'char');
    fwrite(fid, ctab(i,1), 'int');
    fwrite(fid, ctab(i,2), 'int');
    fwrite(fid, ctab(i,3), 'int');
    fwrite(fid, ctab(i,4), 'int');
end
fclose(fid);
%=========================================================================%
return;


