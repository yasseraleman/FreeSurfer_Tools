function OutAnnotFile = Merge_Annots(AnnotFiles,OutAnnotFile);
%
% Syntax :
% OutAnnotFile = Merge_Annots(AnnotFiles,OutAnnotFile);
%
% This script merges different binary annot files and creates a single
% annot file with different structures labels.
% 
%
% Input Parameters:
%      AnnotFiles               : Annotation files
%     OutAnnotFile              : Merged Annotation Filename
%
% Output Parameters:
%     OutAnnotFile              : Merged Annotation Filename
%
% See also: 
%__________________________________________________
% Authors: Yasser Aleman Gomez 
% LIM, HUGGM
% March 22th 2013
% Version $1.0

AnnotFiles = strvcat('/media/COSAS/Test/Joost/SCAN2013ANNOT/rh.atc.annot','/media/COSAS/Test/Joost/SCAN2013ANNOT/rh.BA10medial.annot',...
    '/media/COSAS/Test/Joost/SCAN2013ANNOT/rh.pSTS.annot','/media/COSAS/Test/Joost/SCAN2013ANNOT/rh.TPJp.annot');
OutAnnotFile ='/media/COSAS/Test/Joost/SCAN2013ANNOT/rh.SOC.annot';

Ns = size(AnnotFiles,1);
ctabtot.numEntries = Ns+1;
ctabtot.orig_tab = 'Merged_AnnotFile';
ctabtot.struct_names{1,1} = 'unknown';
ctabtot.table = [ 0 0 0 0 0];
col = [0 0 0;1 0 0;0 1 0; 0 0 1; 1 1 0;0 1 1;1 0 1;1 0.5 0;0 0.5 1;1 0 0.5;0.5 1 0; 0.25 0.5 1;0.6 0.3 0.52;0.7 0.5 0.9];
Ncolor = size(col,1);
re = floor(Ns/Ncolor); col = repmat(col,[re+1 1]);
for i = 1:Ns
    [txt,ctab] = read_cfiles(deblank(AnnotFiles(i,:)));
    tempname = lower(char(ctab.struct_names));
    indu = find(ismember(tempname(:,1:7),'unknown','rows') == 1);
    indcc = find(ismember(tempname(:,1:14),'corpuscallosum','rows') == 1);
    ctab.table([indu indcc],:) = [];
    ctab.struct_names([indu indcc]) = [];
    ctab.struct_names
    if i == 1
        Ntxt = txt;
    end
    ind = find(txt);
    Ntxt(ind) = i;
    ctabtot.struct_names = [ctabtot.struct_names;ctab.struct_names];
    ctabtot.table = [ctabtot.table;ctab.table];
end
Ns = size(ctabtot.table,1);
colors = floor(col*255);
ctab = [colors colors(:,1)*0 colors(:,1)+colors(:,2)*2^8+colors(:,3)*2^16];
Ntxttotal = Ntxt*0;
for i = 2:Ns
    ind = find(Ntxt == i-1);%ctabtot.table(i,5));
    Ntxttotal(ind) = ctab(i,5);
    ctabtot.table(i,:) = ctab(i,:);
end

OutAnnot = save_annotfiles(Ntxttotal,OutAnnotFile,ctabtot);
return;