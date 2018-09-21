function NewAnnot = Intersecting_AnnotFiles(AnnotFiles,NewAnnot);
% Syntax :
% NewAnnot = Intersecting_AnnotFiles(AnnotFiles,NewAnnot);
%
% This function creates an annot file using the information contained in
% two different annot files.
%
% Input Parameters:
%   AnnotFiles  : AnnotFiles strvcat('Annot1','Annot2')
%   NewAnnot    : New Annot File
%
%
% Output Parameters:
%   NewAnnot    : New Annot File
%
% Related references:
%
%
% See also: Smooth_Surf Surf_Comp Red_Surf Plot_Surf Atlas_Surf Exp_Surf
%__________________________________________________
% Authors: Yasser Aleman Gomez
% Neuroimaging Department
% Cuban Neuroscience Center
% September 29th 2007
% Version $1.0

AnnotFiles = strvcat('/media/COSAS/Test/freesurfer/ch2/label/lh.Yeo2011_7Networks_N1000.annot','/media/COSAS/Test/freesurfer/ch2/label/lh.aparc.annot');
SurfFile = '/media/COSAS/Test/freesurfer/ch2/surf/lh.inflated';
NewAnnot = '/media/COSAS/Test/freesurfer/ch2/label/lh.aparc+Yeo2011_7Networks_N1000.annot';
%=========================Main Program====================================%
% Removing Medial Wall
[txt1,ctab1] = read_cfiles(deblank(AnnotFiles(1,:)));
[txt2,ctab2] = read_cfiles(deblank(AnnotFiles(2,:)));
nt = char(ctab1.struct_names);
indd = find(ismember(lower(nt(:,1:7)),'unknown','rows'));
indd = ctab1.table(indd,5);
txt1(txt1== indd) = 0;
indd = find(ismember(txt1,ctab1.table(:,5)) == 0);
txt1(indd) = 0;
nt = char(ctab2.struct_names);
indd = find(ismember(lower(nt(:,1:7)),'unknown','rows'));
indd = ctab2.table(indd,5);
txt2(txt2 == indd) = 0;
indd = find(ismember(txt2,ctab2.table(:,5)) == 0);
txt2(indd) = 0;

ind1 = find(txt1 == 0);
ind2 = find(txt2 == 0);
indw = unique([ind1;ind2]);
[OutFiles, SurfF] = Exp_Surf(deblank(SurfFile), '0', '','', 'imp','n');
Surf = SurfF{1};
% 
sts = unique(txt1);
sts(sts == 0) = [];
newtxt = txt2*0;
cont = 0;
Newnames = 'unknown';
for i = 1:length(sts)
    ind = find(txt1 == sts(i));
    indn = find(ctab1.table(:,5) == sts(i));
    if isempty(indn)
        name = 'unknown'
    else
    name = ctab1.struct_names{indn};
    end
    
    temp = txt2(ind);
    sts2 = unique(temp);
    sts2(sts2 == 0) = [];
    for j = 1:length(sts2)
        ind2 = find(temp == sts2(j));
        if ~isempty(ind2)
            cont = cont+ 1;
            indn = find(ctab2.table(:,5) == sts2(j));
            name2 = ctab2.struct_names{indn};
            newtxt(ind(ind2)) = cont;
            Newnames = strvcat(Newnames,[name '_' name2 ]);
        end
    end
end
newtxt(indw) = 0;
[colortable,labels] = Create_FS_Colortable(newtxt,Newnames);
nnewtxt = labels;
NewAnnot = save_annotfiles(nnewtxt,NewAnnot,colortable);
%========================End of main program==============================%
return;


function [curv, fnum] = read_char(fname);

fid = fopen(fname, 'rb', 'b') ;
if (fid < 0)
    str = sprintf('could not open file %s.', fname) ;
    error(str) ;
end
% vnum = fread3(fid) ;
b1 = fread(fid, 1, 'uchar') ;
b2 = fread(fid, 1, 'uchar') ;
b3 = fread(fid, 1, 'uchar') ;
vnum = bitshift(b1, 16) + bitshift(b2,8) + b3 ;

NEW_VERSION_MAGIC_NUMBER = 16777215;
if (vnum == NEW_VERSION_MAGIC_NUMBER)
    vnum = fread(fid, 1, 'int32') ;
    fnum = fread(fid, 1, 'int32') ;
    vals_per_vertex = fread(fid, 1, 'int32') ;
    curv = fread(fid, vnum, 'float') ;

    fclose(fid) ;
else
    b1 = fread(fid, 1, 'uchar') ;
    b2 = fread(fid, 1, 'uchar') ;
    b3 = fread(fid, 1, 'uchar') ;
    vnum = bitshift(b1, 16) + bitshift(b2,8) + b3 ;
    curv = fread(fid, vnum, 'int16') ./ 100 ;
    fclose(fid) ;
end

function [retval] = fread3(fid)

b1 = fread(fid, 1, 'uchar') ;
b2 = fread(fid, 1, 'uchar') ;
b3 = fread(fid, 1, 'uchar') ;
retval = bitshift(b1, 16) + bitshift(b2,8) + b3 ;




