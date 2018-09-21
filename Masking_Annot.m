function OutAnnot = Masking_Annot(AnnotFile,CurvFile,thr,OutAnnot);
%
% Syntax :
% OutAnnot = Masking_Annot(AnnotFile,CurvFile,thr,OutAnnot);
%
% This script masks annotation files using the curvature information
% 
%
% Input Parameters:
%   AnnotFile            : Annotation File
%   CurvFile             : Curvature file
%   thr                  : Curvature threshold for binarization
%   OutAnnot             : New Annotation filename
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
if nargin <4
    [pth,nm,ext] = fileparts(AnnotFile);
    OutAnnot  = [pth filesep  nm '_sulci' ext];
end
if nargin <3
    thr = 0;
    [pth,nm,ext] = fileparts(AnnotFile);
    OutAnnot  = [pth filesep  nm '_sulci' ext];
end
%=========================================================================%
%=====================     Main Program    ===============================%
[label1, txt2, colortable] = read_annotation(AnnotFile); %
[txt] = read_cfiles(CurvFile);
inds = find(txt>thr);
indg = find(txt<=thr);
txt(inds) = 1;
txt(indg) = 0;
% Find unknown
namestemp = char(colortable.struct_names);
ind = ismember(namestemp(:,1:7),'unknown','rows');
stid = colortable.table(find(ind),5);
txt2(indg) = stid;
%=========================================================================%
%===================== Saving new annot File =============================%
save_annotfiles(txt2,OutAnnot, colortable);
%=========================================================================%
return

function OutAnnot = save_annotfiles(labelsd,OutAnnot,colortable);
%
% Syntax :
% OutAnnot = Masking_Annot(AnnotFile,CurvFile,thr,OutAnnot);
%
% This script creates an annotation file for cortical parcelation values.
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
Nstruct = length(unique(labelsd));
if nargin<3
    ctabaux =randi([10 255+255*2^8+255*2^16],[1 Nstruct])';[X,Y,Z] = ind2sub([255 255 255],ctabaux);
    ctabaux(1) = 0;
    ctab = [X(:) Y(:) Z(:) zeros(size(Y,1),1)];
    ctab(1,:) = [130 130 130 0];
    labels = labelsd.*0;
    for i = 1:size(ctabaux,1)
        ind = find(labelsd==i);
        labels(ind) =ones(size(ind,1),1)*ctabaux(i);
    end
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
Avalues(1:2:Np*2-1) = [0:Np-1]';
fwrite(fid, Avalues, 'int');
fwrite(fid, 1, 'int');
fwrite(fid, -2, 'int');
fwrite(fid,length(unique(labels)),'int');
fwrite(fid, length(OutAnnot), 'int');
fwrite(fid, OutAnnot, '*char');
fwrite(fid,length(unique(labels)),'int');
for i = 1:Nstruct
    fwrite(fid, i-1, 'int');
    if nargin<3
        len = length(sprintf('%.5d',i-1));
        name = sprintf('%.5d',i);
    else
        name = deblank(colortable.struct_names{i});
        len = length(name);
    end
    fwrite(fid, len, 'int');
    fwrite(fid, name, 'char');
    fwrite(fid, ctab(i,1), 'int');
    fwrite(fid, ctab(i,2), 'int');
    fwrite(fid, ctab(i,3), 'int');
    fwrite(fid, ctab(i,4), 'int');
end
fclose(fid);
%=========================================================================%
return;

