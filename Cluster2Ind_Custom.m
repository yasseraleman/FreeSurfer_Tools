function OutFiles = Cluster2Ind_Custom(IndSubjecsFreesDir,Idfile, Clusters, hemi);
%
% Syntax :
% OutFile = Cluster2Ind_Custom(IndSubjecsFreesDir,Idfile, Clusters,hemi);
%
% Example:  OutFiles = Cluster2Ind_Custom('/media/COSAS/Test/freesurfer',...
%                      '/media/COSAS/Test/Alejandro/Ids.txt',...
%         strvcat('/media/COSAS/Test/Alejandro/C/mc-z.abs.sig.ocn.annot',...
%           '/media/COSAS/Test/Alejandro/C/mc-z.abs.sig.test.annot'),'lh');
%
% This script obtains different morphometric measures for a specified 
% cluster annotation file.
%
%
% Input Parameters:
%    IndSubjecsFreesDir         : FreeSurfer Directory
%       Idfile                  : Text file containing the Ids List 
%       Clusters                : Clusters Annotation Files
%       hemi                    : Hemisphere (lh or rh)
%
% Output Parameters:
%      OutFiles                 : Output stat files 
%
% See also:
%__________________________________________________
% Authors: Yasser Aleman Gomez
% LIM, HUGGM
% March 22th 2013
% Version $1.0


%% =======================  FreeSurfer IDs  ============================ %%
Ids = char(textread(Idfile,'%s'));
%Ids = strvcat('ch2','fsaverage');
%% =============== End of Detecting FreeSurfer IDs  ==================== %%

%% ============ Creating variables to process the average subject ====== %%
% ----------- Setting Enviroment for Freesurfer Average ------------- %
[a,temp] = system('echo $FREESURFER_HOME');
temp = temp';temp = deblank(temp(:)');
optst.pipe.freesdir = [temp filesep 'subjects'];
%
% --------------------- Other variables --------------------------------- %
optst.pipe.subjId = 'fsaverage'; % Subject ID

%% ============ End of Setting Enviroment for Freesurfer Average ======= %%


%% ======== Creating variables to process the Individual subjects ====== %%
% ----------- Setting Enviroment for Freesurfer Directory ------------ %
opts.pipe.freesdir = IndSubjecsFreesDir;
setenv('SUBJECTS_DIR',opts.pipe.freesdir);
if ~exist([opts.pipe.freesdir filesep 'fsaverage'],'dir')
    cad = ['cp -r ' optst.pipe.freesdir filesep 'fsaverage ' opts.pipe.freesdir];
    system(cad);
end

%% ============ End of Setting Enviroment for Freesurfer Average ======= %%


%% ================== Subjects Processing ============================== %%
tempd = '';
Nc = size(Clusters,1);
OutFiles = '';
for i = 1:Nc
    Sfile = deblank(Clusters(i,:));
    [pth,nm,ext] = fileparts(Sfile);
    annotId = nm;
    %% ================= Creating Color Codes File ===================== %%
    ColorFile = [opts.pipe.freesdir filesep 'fsaverage' filesep 'label' filesep 'Temp_ColorFile.txt'];
    [txt,ctab] = read_cfiles(Sfile);
    Nclust = unique(txt);
    Nclust(Nclust == 0) = [];
    indo = find(ismember(ctab.table(:,5),Nclust)==0);
    colors = [125 125 125];
    ctab.table(indo,:) =  [colors colors(:,1)*0 colors(:,1)+colors(:,2)*2^8+colors(:,3)*2^16];
    ind = find(txt == 0);
    txt(ind) = ctab.table(indo,5);
    NAnnotFile = save_annotfiles(txt,[opts.pipe.freesdir filesep 'fsaverage' filesep 'label' filesep hemi '.temp.annot'],ctab);
    
    [txtlab,ctablab] = read_cfiles([opts.pipe.freesdir filesep 'fsaverage' filesep 'label' filesep hemi '.aparc.annot']); % Anatomical labels
    
    names = 'unknown';
    for j = 1:length(Nclust)
        ind = find(ctab.table(:,5) == Nclust(j));
        indl = find(txt == Nclust(j));
        stids = accumarray(txtlab(indl),txtlab(indl)*0+1);
        indm = find(stids == max(stids));
        indaparc = find(ctablab.table(:,5) == indm(1)); % APARC INdex
        r(j) = ctab.table(ind,1);
        g(j) = ctab.table(ind,2);
        b(j) = ctab.table(ind,3);
        names = strvcat(names,['Cluster-No_'  num2str(j) '-' deblank(char(ctablab.struct_names(indaparc)))]);
    end
    r = num2str([125;r(:)]);
    g = num2str([125;g(:)]);
    b = num2str([125;b(:)]);
    ids = num2str([0:length(Nclust)]');
    Total = [ids repmat('   ',[size(ids,1) 1]) names repmat('             ',[size(ids,1) 1]) r repmat(' ',[size(ids,1) 1]) g repmat(' ',[size(ids,1) 1]) b repmat(' ',[size(ids,1) 1]) num2str(zeros(size(b,1),1))];
    fid = fopen(ColorFile,'wt');
    for j = 1:size(Total,1)
       fprintf(fid,'%s\n',Total(j,:));
    end
    fclose(fid);
    %% ============== End of Creating Color Codes File ================= %%
    
    %% ==== Creating Multiscale Atlas for template subject ============= %%
    setenv('SUBJECTS_DIR',opts.pipe.freesdir);
    if hemi == 'lh'
        cadl = ['mris_ca_train -n 2 -t ' ColorFile ' lh sphere.reg temp fsaverage ' opts.pipe.freesdir filesep 'fsaverage' filesep 'label' filesep 'lh.temp.gcs'];
        system(cadl);
    elseif hemi == 'rh'
        cadr = ['mris_ca_train -n 2 -t ' ColorFile ' rh sphere.reg temp fsaverage ' opts.pipe.freesdir filesep 'fsaverage' filesep 'label' filesep 'rh.temp.gcs'];
        system(cadr);
    end
    
    %% ==== End of Creating Multiscale Atlas for template subject ====== %%
    Nsubj = size(Ids,1);
    for z= 1:Nsubj
        Id = deblank(Ids(z,:));
        disp(strvcat(' ',' '));
        disp(['Processing =======>  Cluster Files: ' num2str(i) ' of ' num2str(Nc) ' . =====> Subject ID: ' Id ' . ---  ' num2str(z) ' of ' num2str(Nsubj)]);
        %% =========== Moving Surfaces to individual space ============= %%
        if hemi == 'lh'
            cadl = ['mris_ca_label -orig white -novar -t ' ColorFile ' ' Id ' lh sphere.reg ' opts.pipe.freesdir filesep 'fsaverage' filesep 'label' filesep 'lh.temp.gcs '  opts.pipe.freesdir  filesep Id filesep 'label' filesep 'lh.' annotId '.annot'];
            system(cadl);
        elseif hemi == 'rh'
            cadr = ['mris_ca_label -orig white -novar -t ' ColorFile ' ' Id ' rh sphere.reg ' opts.pipe.freesdir filesep 'fsaverage' filesep 'label' filesep 'rh.temp.gcs '  opts.pipe.freesdir  filesep Id filesep 'label' filesep 'rh.' annotId '.annot'];
            system(cadr);
        end
        %% =========== End of Moving Surfaces to individual space ====== %%
        
        %% =================== Reading files ================================ %
        % Parcellation Map
        Annotfile = [opts.pipe.freesdir  filesep Id filesep 'label' filesep hemi '.' annotId '.annot'];
        
        % Surface
        SurfFile = [opts.pipe.freesdir  filesep Id filesep 'surf' filesep hemi '.pial'];
        
        % Curvature Map
        CurvFile = [opts.pipe.freesdir  filesep Id filesep 'surf' filesep hemi '.curv'];
        
        % Thickness Maps
        ThickFile = [opts.pipe.freesdir  filesep Id filesep 'surf' filesep hemi '.thickness'];
        
        %% =================== End of Reading files ========================= %
        
        
        %% ================= Computing Things =================================== %
        %----- Reading Surfaces
        [OutFiles, SurfF] = Exp_Surf(SurfFile, '0', '','', 'imp','n');
        Surf= SurfF{1};
        Surf.Name = 'Pial_Surface';
        
        %----- Reading Annot Files
        [txt,ctabr] = read_cfiles(Annotfile);
        ctabr.table = [ctabr.table 1000+[0:size(ctabr.table,1)-1]' ];
        tempname = char(ctabr.struct_names);
        indu = find(ismember(tempname(:,1:7),'unknown','rows') == 1);
        ctabr.table([indu],:) = [];
        ctabr.struct_names([indu]) = [];
        
        %----- Reading Curvature
        [curv,ctab] = read_cfiles(CurvFile);
        
        
        %----- Reading Thickness Map
        [cth,ctab] = read_cfiles(ThickFile);
        
        
        % ================ Computing Morphometric measures ===================%
        Nr = size(ctabr.table,1);
        fv.vertices = Surf.SurfData.vertices;
        for j = 1:Nr
            ind = find(txt == ctabr.table(j,5));
            nvert(z,j) = length(ind); % Number of vertices
            cthm(z,j) = mean(cth(ind)); % Mean Cortical Thickness
            cths(z,j) = std(cth(ind)); % Std Cortical Thickness
            curva(z,j) = mean(curv(ind)); % Mean curvature
            
            %----- Computing Area
            indf = find(sum(ismember( Surf.SurfData.faces,ind)')==3);
            At = 0;
            fv.faces = Surf.SurfData.faces(indf,:);
            N = size(fv.faces,1);
            for k = 1:N;
                di = dist(fv.vertices(fv.faces(k,:),:)');
                A = abs(di(1,2));
                B = abs(di(1,3));
                C = abs(di(2,3));
                p = (A+B+C)/2;
                Ar = sqrt(p*(p-A)*(p-B)*(p-C));
                %Ar = (A*B/2)*(sqrt(1-((A^2+B^2-C^2)^2)/((2*A*B)^2)));
                At = At+Ar;
            end
            area(z,j) = real(At);% mm^2
            
        end
        vol(z,:) =  area(z,:).*cthm(z,:);
    end
    TotalVals = zeros(size(area,1),1);
    for j = 1:Nr
        TotalVals = [TotalVals nvert(:,j) cthm(:,j) cths(:,j)   area(:,j) vol(:,j) curva(:,j)];
    end
    TotalVals(:,1) = [];
    cads = strvcat('SubjIds',Ids);
    a = char(ctabr.struct_names);[Nrow,Ncol] = size(a);
    a = repmat(char(ctabr.struct_names),[1 6]);
    a = a';
    a = a(:);
    Temp = reshape(a,[Ncol length(a)/Ncol])';
    Varnames = repmat(strvcat('NumVert','Mean_Cthickness','Std_Cthickness','Area(mm^2)','Volume(mm^3)','Mean_Curvature'),[Nr 1]);
    Varnames = [Varnames Temp];
    NVarnames = '';
    for j = 1:size(Varnames,1)
        Temp = deblank(Varnames(j,:));
        ind = find(isspace(Temp));
        Temp(ind) = [];
        NVarnames = strvcat(NVarnames,Temp);
    end
    for j = 1:size(TotalVals,2)
        cads = [cads repmat('    ',[Nsubj+1 1]) strvcat([deblank(NVarnames(j,:))],num2str(TotalVals(:,j))) repmat('    ',[Nsubj+1 1])];
    end
    %% ================= Saving Cluster Values ========================== %
    OutFile = [pth filesep nm '_Clust_Stats.txt'];
    fid = fopen(OutFile,'wt');
    for j = 1:size(cads,1)
        fprintf(fid, '%s\n', cads(j,:));
    end
    fclose all;
    disp(['Output file: ==>>  ' OutFile]);
    OutFiles = strvcat(OutFiles,OutFile);
end
return;


function [txt,ctab] = read_cfiles(CFile);
%
% Syntax :
% [txt,ctab] = read_cfiles(CFile);
%
% This function reads surface textures. It can accept text files,
% freesurfer annotation files and Brainvisa texture files.
%
% Input Parameters:
%   CFile         : Texture File
%
% Output Parameters:
%      txt        : Vertices textures.
%
%   colortable    : Struct variable similar to freesurfer's colortable. It
%                  includes 2 fields:
%                  1. struct_names: Cellarray containing
%                  structures names.
%                  2. table. Nstructures X 6 Matrix containing some
%                  material parameters (Colors, Transparency ( Tr), Ids)
%                  Table order (R G B 0 Id Tr)
%
% Related references:
%
%
% See also:
%__________________________________________________
% Authors: Yasser Aleman Gomez
% LIM
% September 12th 2012
% Version $1.0
[pth, nm, ext]= fileparts(deblank(CFile));
try
    txt = single((textread(CFile,'%f',Npoints)));
    ctab.table = 0;
    ctab.struct_names{1,1} = '';
catch
    try
        [txtt,Format] = read_texBrainvisa(CFile);
        txt = txtt.Values;
        ctab.table = 0;
        ctab.struct_names{1,1} = '';
    catch
        try
               [txt, temp] = read_char(CFile);
                ctab.table = 0;
                ctab.struct_names{1,1} = '';
        catch
            try
                if ext(1:4) ~= '.mgh'
                    [vertices, txt, ctab] = read_annotation(CFile);
                    if isempty(ctab)
                        ctab(1).table= 0;
                    end
                else
                    error;
                end
            catch
                [txt, M, mr_parms, volsz] = load_mgh(CFile);
                ctab.table= 0;
                if isempty(txt)
                    fid = fopen(deblank(CFile),'rb','b');
                    Np = fread(fid, 1, 'int32');
                    Avalues = fread(fid, Np*2, 'int');
                    vb = fread(fid, 1, 'int');
                    if(isempty(vb))
                        disp('No Colortable found.');
                        fclose(fid);
                        return;
                    end
                    tt = fread(fid, 1, 'int');
                    if(tt > 0)
                        temp0 = fread(fid, 1, 'int');
                        temp1 = fread(fid, temp0, 'char')';
                        for k = 1:tt
                            temp2 = fread(fid, 1, 'int');
                            temp3 = fread(fid, temp2, 'char')';
                            ctab.struct_names{k,1} = char(temp3);
                            temp4 = fread(fid, 1, 'int');
                            temp5 = fread(fid, 1, 'int');
                            temp6 = fread(fid, 1, 'int');
                            temp7 = fread(fid, 1, 'int')
                            ord(k) = temp4 + temp5*2^8 + temp6*2^16 + temp7*2^24;
                            ctab.table(k,5) = ord(k);
                        end
                    else
                        tt1 = fread(fid, 1, 'int');
                        temp0 =fread(fid, 1, 'int');
                        temp1 = fread(fid, temp0, 'char')';
                        Nstructs = fread(fid, 1, 'int');
                        for k = 1:Nstructs
                            st = fread(fid, 1, 'int')+1;
                            len = fread(fid, 1, 'int');
                            temp2 = fread(fid, len, 'char')';
                            ctab.struct_names{k,1} = char(temp2);
                            temp3 = fread(fid, 1, 'int');
                            ctab.table(k,1)=temp3;
                            temp4 = fread(fid, 1, 'int');
                            ctab.table(k,2)=temp4;
                            temp5 = fread(fid, 1, 'int');
                            ctab.table(k,3)=temp5;
                            temp6 = fread(fid, 1, 'int');
                            ord(k) = temp3 + temp4*2^8 + temp5*2^16 + temp6*2^24;
                            ctab.table(k,5) = ord(k);
                        end
                    end
                    %ctab(35,:) =[140 220 220];
                    txt = Avalues(2:2:Np*2);
                end
            end
        end
    end
end
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

function [vol, M, mr_parms, volsz] = load_mgh(fname,slices,frames,headeronly)
% [vol, M, mr_parms, Mdc, volsz] = load_mgh(fname,<slices>,<frames>,<headeronly>)
%
% fname - path of the mgh file
%
% slices - list of one-based slice numbers to load. All
%   slices are loaded if slices is not specified, or
%   if slices is empty, or if slices(1) <= 0.
%
% frames - list of one-based frame numbers to load. All
%   frames are loaded if frames is not specified, or
%   if frames is empty, or if frames(1) <= 0.
%
% M is the 4x4 vox2ras transform such that
% y(i1,i2,i3), xyz1 = M*[i1 i2 i3 1] where the
% indices are 0-based. If the input has multiple frames,
% only the first frame is read.
%
% mr_parms = [tr flipangle te ti fov]
%
% volsz = size(vol). Helpful when using headeronly as vol is [].
%
% See also: save_mgh, vox2ras_0to1
%


%
% load_mgh.m
%
% Original Author: Bruce Fischl
% CVS Revision Info:
%    $Author: greve $
%    $Date: 2009/07/01 17:13:08 $
%    $Revision: 1.16.2.2 $
%
% Copyright (C) 2002-2007,
% The General Hospital Corporation (Boston, MA).
% All rights reserved.
%
% Distribution, usage and copying of this software is covered under the
% terms found in the License Agreement file named 'COPYING' found in the
% FreeSurfer source code root directory, and duplicated here:
% https://surfer.nmr.mgh.harvard.edu/fswiki/FreeSurferOpenSourceLicense
%
% General inquiries: freesurfer@nmr.mgh.harvard.edu
% Bug reports: analysis-bugs@nmr.mgh.harvard.edu
%

vol = [];
M = [];
mr_parms = [];
volsz = [];

if(nargin < 1 | nargin > 4)
    msg = 'USAGE: [vol M] = load_mgh(fname,<slices>,<frames>,<headeronly>)';
    fprintf('%s',msg);
    return;
end

% unzip if it is compressed
if (strcmpi(fname((strlen(fname)-3):strlen(fname)), '.MGZ') | ...
        strcmpi(fname((strlen(fname)-3):strlen(fname)), '.GZ'))
    rand('state', sum(100*clock));
    gzipped =  round(rand(1)*10000000 + ...
        sum(int16(fname))) + round(cputime);
    ind = findstr(fname, '.');
    new_fname = sprintf('/tmp/tmp%d.mgh', gzipped);
    if(strcmp(computer,'MAC') || strcmp(computer,'MACI') || ismac)
        unix(sprintf('gunzip -c %s > %s', fname, new_fname)) ;
    else
        unix(sprintf('zcat %s > %s', fname, new_fname)) ;
    end
    fname = new_fname ;
else
    gzipped = -1 ;
end


if(exist('slices')~=1) slices = []; end
if(isempty(slices)) slices = 0; end
if(slices(1) <= 0) slices = 0; end

if(exist('frames')~=1) frames = []; end
if(isempty(frames)) frames = 0; end
if(frames(1) <= 0) frames = 0; end

if(exist('headeronly')~=1) headeronly = 0; end

fid    = fopen(fname, 'rb', 'b') ;
if(fid == -1)
    fprintf('ERROR: could not open %s for reading\n',fname);
    return;
end
v       = fread(fid, 1, 'int') ;
if(isempty(v))
    fprintf('ERROR: problem reading fname\n');
    if(gzipped >=0) unix(sprintf('rm %s', fname)); end
end
ndim1   = fread(fid, 1, 'int') ;
ndim2   = fread(fid, 1, 'int') ;
ndim3   = fread(fid, 1, 'int') ;
nframes = fread(fid, 1, 'int') ;
type    = fread(fid, 1, 'int') ;
dof     = fread(fid, 1, 'int') ;

if(slices(1) > 0)
    ind = find(slices > ndim3);
    if(~isempty(ind))
        fprintf('ERROR: load_mgh: some slices exceed nslices\n');
        return;
    end
end

if(frames(1) > 0)
    ind = find(frames > nframes);
    if(~isempty(ind))
        fprintf('ERROR: load_mgh: some frames exceed nframes\n');
        return;
    end
end

UNUSED_SPACE_SIZE= 256;
USED_SPACE_SIZE = (3*4+4*3*4);  % space for ras transform

unused_space_size = UNUSED_SPACE_SIZE-2 ;
ras_good_flag = fread(fid, 1, 'short') ;
if (ras_good_flag)
    delta  = fread(fid, 3, 'float32') ;
    Mdc    = fread(fid, 9, 'float32') ;
    Mdc    = reshape(Mdc,[3 3]);
    Pxyz_c = fread(fid, 3, 'float32') ;
    
    D = diag(delta);
    
    Pcrs_c = [ndim1/2 ndim2/2 ndim3/2]'; % Should this be kept?
    
    Pxyz_0 = Pxyz_c - Mdc*D*Pcrs_c;
    
    M = [Mdc*D Pxyz_0;  ...
        0 0 0 1];
    ras_xform = [Mdc Pxyz_c; ...
        0 0 0 1];
    unused_space_size = unused_space_size - USED_SPACE_SIZE ;
end

fseek(fid, unused_space_size, 'cof') ;
nv = ndim1 * ndim2 * ndim3 * nframes;
volsz = [ndim1 ndim2 ndim3 nframes];

MRI_UCHAR =  0 ;
MRI_INT =    1 ;
MRI_LONG =   2 ;
MRI_FLOAT =  3 ;
MRI_SHORT =  4 ;
MRI_BITMAP = 5 ;

% Determine number of bytes per voxel
switch type
    case MRI_FLOAT,
        nbytespervox = 4;
    case MRI_UCHAR,
        nbytespervox = 1;
    case MRI_SHORT,
        nbytespervox = 2;
    case MRI_INT,
        nbytespervox = 4;
end

if(headeronly)
    fseek(fid,nv*nbytespervox,'cof');
    if(~feof(fid))
        [mr_parms count] = fread(fid,4,'float32');
        if(count ~= 4)
            fprintf('WARNING: error reading MR params\n');
        end
    end
    fclose(fid);
    if(gzipped >=0)  unix(sprintf('rm %s', fname));  end
    return;
end


%------------------ Read in the entire volume ----------------%
if(slices(1) <= 0 & frames(1) <= 0)
    switch type
        case MRI_FLOAT,
            vol = fread(fid, nv, 'float32') ;
        case MRI_UCHAR,
            vol = fread(fid, nv, 'uchar') ;
        case MRI_SHORT,
            vol = fread(fid, nv, 'short') ;
        case MRI_INT,
            vol = fread(fid, nv, 'int') ;
    end
    
    if(~feof(fid))
        [mr_parms count] = fread(fid,4,'float32');
        if(count ~= 4)
            fprintf('WARNING: error reading MR params\n');
        end
    end
    fclose(fid) ;
    if(gzipped >=0)  unix(sprintf('rm %s', fname));  end
    
    nread = prod(size(vol));
    if(nread ~= nv)
        fprintf('ERROR: tried to read %d, actually read %d\n',nv,nread);
        vol = [];
        return;
    end
    vol = reshape(vol,[ndim1 ndim2 ndim3 nframes]);
    
    return;
end

%----- only gets here if a subest of slices/frames are to be loaded ---------%


if(frames(1) <= 0) frames = [1:nframes]; end
if(slices(1) <= 0) slices = [1:ndim3]; end

nvslice = ndim1 * ndim2;
nvvol   = ndim1 * ndim2 * ndim3;
filepos0 = ftell(fid);
vol = zeros(ndim1,ndim2,length(slices),length(frames));
nthframe = 1;
for frame = frames
    
    nthslice = 1;
    for slice = slices
        filepos = ((frame-1)*nvvol + (slice-1)*nvslice)*nbytespervox + filepos0;
        fseek(fid,filepos,'bof');
        
        switch type
            case MRI_FLOAT,
                [tmpslice nread]  = fread(fid, nvslice, 'float32') ;
            case MRI_UCHAR,
                [tmpslice nread]  = fread(fid, nvslice, 'uchar') ;
            case MRI_SHORT,
                [tmpslice nread]  = fread(fid, nvslice, 'short') ;
            case MRI_INT,
                [tmpslice nread]  = fread(fid, nvslice, 'int') ;
        end
        
        if(nread ~= nvslice)
            fprintf('ERROR: load_mgh: reading slice %d, frame %d\n',slice,frame);
            fprintf('  tried to read %d, actually read %d\n',nvslice,nread);
            fclose(fid);
            if(gzipped >=0) unix(sprintf('rm %s', fname)); end
            return;
        end
        
        vol(:,:,nthslice,nthframe) = reshape(tmpslice,[ndim1 ndim2]);
        nthslice = nthslice + 1;
    end
    
    nthframe = nthframe + 1;
end

% seek to just beyond the last slice/frame %
filepos = (nframes*nvvol)*nbytespervox + filepos0;
fseek(fid,filepos,'bof');

if(~feof(fid))
    [mr_parms count] = fread(fid,5,'float32');
    if(count < 4)
        fprintf('WARNING: error reading MR params\n');
    end
end

fclose(fid) ;
if(gzipped >=0) unix(sprintf('rm %s', fname)); end

return;

function [vertices, label, colortable] = Read_Brain_Annotation(filename)
% [vertices, label, colortable] = Read_Brain_Annotation(annotfilename.annot)
%
% vertices expected to be simply from 0 to number of vertices - 1;
% label is the vector of annotation
%
% colortable is empty struct if not embedded in .annot. Else, it will be
% a struct.
% colortable.numEntries = number of Entries
% colortable.orig_tab = name of original colortable
% colortable.struct_names = list of structure names (e.g. central sulcus and so on)
% colortable.table = n x 5 matrix. 1st column is r, 2nd column is g, 3rd column
% is b, 4th column is flag, 5th column is resultant integer values
% calculated from r + g*2^8 + b*2^16 + flag*2^24. flag expected to be all 0.


%
% read_annotation.m
%
% Original Author: Bruce Fischl
% CVS Revision Info:
%    $Author: nicks $
%    $Date: 2007/01/10 22:55:09 $
%    $Revision: 1.4 $
%
% Copyright (C) 2002-2007,
% The General Hospital Corporation (Boston, MA). 
% All rights reserved.
%
% Distribution, usage and copying of this software is covered under the
% terms found in the License Agreement file named 'COPYING' found in the
% FreeSurfer source code root directory, and duplicated here:
% https://surfer.nmr.mgh.harvard.edu/fswiki/FreeSurferOpenSourceLicense
%
% General inquiries: freesurfer@nmr.mgh.harvard.edu
% Bug reports: analysis-bugs@nmr.mgh.harvard.edu
%

fp = fopen(filename, 'r', 'b');

if(fp < 0)
   disp('Annotation file cannot be opened');
   return;
end

A = fread(fp, 1, 'int');

tmp = fread(fp, 2*A, 'int');
vertices = tmp(1:2:end);
label = tmp(2:2:end);

bool = fread(fp, 1, 'int');
if(isempty(bool)) %means no colortable
   disp('No Colortable found.');
   colortable = struct([]);
   fclose(fp);
   return; 
end

if(bool)
    
    %Read colortable
    numEntries = fread(fp, 1, 'int');

    if(numEntries > 0)
        
        disp(['Reading from Original Version']);
        colortable.numEntries = numEntries;
        len = fread(fp, 1, 'int');
        colortable.orig_tab = fread(fp, len, '*char')';
        colortable.orig_tab = colortable.orig_tab(1:end-1);

        colortable.struct_names = cell(numEntries,1);
        colortable.table = zeros(numEntries,5);
        for i = 1:numEntries
            len = fread(fp, 1, 'int');
            colortable.struct_names{i} = fread(fp, len, '*char')';
            colortable.struct_names{i} = colortable.struct_names{i}(1:end-1);
            colortable.table(i,1) = fread(fp, 1, 'int');
            colortable.table(i,2) = fread(fp, 1, 'int');
            colortable.table(i,3) = fread(fp, 1, 'int');
            colortable.table(i,4) = fread(fp, 1, 'int');
            colortable.table(i,5) = colortable.table(i,1) + colortable.table(i,2)*2^8 + colortable.table(i,3)*2^16 + colortable.table(i,4)*2^24;
        end
        disp(['colortable with ' num2str(colortable.numEntries) ' entries read (originally ' colortable.orig_tab ')']);

    else
        version = -numEntries;
        if(version~=2)    
            disp(['Error! Does not handle version ' num2str(version)]);
        else
            disp(['Reading from version ' num2str(version)]);
        end
        numEntries = fread(fp, 1, 'int');
        colortable.numEntries = numEntries;
        len = fread(fp, 1, 'int');
        colortable.orig_tab = fread(fp, len, '*char')';
        colortable.orig_tab = colortable.orig_tab(1:end-1);
        
        colortable.struct_names = cell(numEntries,1);
        colortable.table = zeros(numEntries,5);
        
        numEntriesToRead = fread(fp, 1, 'int');
        for i = 1:numEntriesToRead
            structure = fread(fp, 1, 'int')+1;
            if (structure < 0)
                disp(['Error! Read entry, index ' num2str(structure)]);
            end
            if(~isempty(colortable.struct_names{structure}))
                disp(['Error! Duplicate Structure ' num2str(structure)]);
            end
            len = fread(fp, 1, 'int');
            colortable.struct_names{structure} = fread(fp, len, '*char')';
            colortable.struct_names{structure} = colortable.struct_names{structure}(1:end-1);
            colortable.table(structure,1) = fread(fp, 1, 'int');
            colortable.table(structure,2) = fread(fp, 1, 'int');
            colortable.table(structure,3) = fread(fp, 1, 'int');
            colortable.table(structure,4) = fread(fp, 1, 'int');
            colortable.table(structure,5) = colortable.table(structure,1) + colortable.table(structure,2)*2^8 + colortable.table(structure,3)*2^16 + colortable.table(structure,4)*2^24;       
        end
        disp(['colortable with ' num2str(colortable.numEntries) ' entries read (originally ' colortable.orig_tab ')']);
    end    
else
    disp('Error! Should not be expecting bool = 0');    
end

fclose(fp);
return;


