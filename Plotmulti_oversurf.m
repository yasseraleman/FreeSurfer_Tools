function Plotmulti_oversurf(SurfFile,CFiles,cl,tr,sw);
% Syntax :
% Plot_oversurf(SurfFile,CFiles,cl,tr,sp,sw);
%
% This function plots the surfaces contained in the cell file Surfa.
%
% Input Parameters:
%   SurfFile    : Surface files.
%   CFiles      : Grosor or atlas text files
%   cl          : Colormap used to see the results.
%   tr          : Transparency vector(values have to be between 0 and 1).
%   sp          : Boolean variable(sp = 1, plot results over the sphere,
%                 sa = 0, do not plot results over the sphere).
%   sw          : Boolean variable(sw = 1, plot surfaces in the same window,
%                 sw = 0, plot surfaces in diferent windows)
%
% Output Parameters:
%
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
SurfFile = '/media/COSAS/Test/freesurfer/fsaverage/surf/lh.inflated';
CFiles = strvcat('/media/COSAS/Test/freesurfer/fsaverage/label/lh.aparc.annot','/media/COSAS/Test/freesurfer/fsaverage/label/lh.aparc_scale02_Nzones0191.annot');



%CFiles = strvcat('/media/COSAS/Test/freesurfer/ch2/label/lh.aparc+Yeo2011_7Networks_N1000.annot');
%=====================Checking Input Parameters===========================%
if ~exist('cl','var')|(isempty(cl))
    cl = 'jet';
end
if ~exist('tr','var')|(isempty(tr))
    tr = 1;
end
if ~exist('sw','var')|(isempty(sw))
    sw = 'y';
end
col = [0 0 0;1 1 1;1 0 0; 0 1 0;0 0 1];
if ~exist('SurfFile','var')|(isempty(SurfFile))
    [SurfFile,sts] = spm_select([1 2],'any','Selecting Surface Files','',cd);
end
wh = whos('SurfFile');
if (strcmp(wh.class,'struct'))|(strcmp(wh.class,'cell'));
    Surfa = SurfFile;
elseif ischar(SurfFile);
    [pth nm ext] = fileparts(deblank(SurfFile));
    if ext(1:4) == '.mat';
        Surf = load('-mat',[pth filesep nm ext(1:4)]);
        if isfield(Surf,'Surf')
            Surf = Surf.Surf;
        end
    else
        [OutFiles, SurfF] = Exp_Surf(deblank(SurfFile), '0', '','', 'imp','n');
        Surf = SurfF{1};
    end
end
if ~exist('CFiles','var')|(isempty(CFiles))
    [CFiles,sts] = spm_select([1 5],'any','Selecting Thickness/Curvature Files','',cd);
end
Nc = size(CFiles,1);
if Nc >5
    CFiles(6:end,:) = [];
end
tic;[Trip] = Vert_Neibp(double(Surf.SurfData.faces),size(Surf.SurfData.vertices,1),size(Surf.SurfData.faces,1));toc;
Temp = sum(Trip);
Trip(:,Temp==0) = [];
Trip = int32(Trip);
temp = Trip(:,3:end);
indz = find(temp == 0);
temp(indz) = 1;

%=========================================================================%
%=========================Main Program====================================%
Npoints = size(Surf.SurfData.vertices,1);
if isfield(Surf.SurfData,'FaceVertexCData')
    rmfield(Surf.SurfData,'FaceVertexCData');
end
indt = 0;
for i = 1:Nc    
    [pth,nm,ext] = fileparts(deblank(CFiles(i,:)));
    CFile = [pth filesep nm deblank(ext)];
    [txt,ctab] = read_cfiles(CFile);
    %% ===================== Painting Surface ============================== %%
    Surf.Is = txt;
    if ctab.table == 0
        [Colors] = Surf_Color(Surf,cl);
        if i == 1
            Surf.SurfData.FaceVertexCData = Colors;
            if Nc >1
                ColorsV = Colors*0;
            end
        else
            ColorsV = Colors;
        end
    else
        sts = unique(txt);
        Nst = length(sts);
        if i == 1
             Surf.SurfData.FaceVertexCData = zeros(size(Surf.SurfData.vertices,1),3);
        end
        for j = 1:Nst
            ind = find(txt==sts(j));
            indc = find(ctab.table(:,5)==sts(j));
            if isempty(indc)
                Matname = 'unknown_structure';
                Color =   [1 1 1];      % Color
            else
                Matname = char(ctab.struct_names{indc});
                Color =  ctab.table(indc,1:3)/255;       % Color
            end
            if i == 1
                Surf.SurfData.FaceVertexCData(ind,:) = repmat(Color,[length(ind) 1 ]);
                if Nc >1
                    ColorsV = Surf.SurfData.FaceVertexCData*0;
                end
            else
                ColorsV(ind,:) = repmat(Color,[length(ind) 1 ]);
            end
        end
    end
    if (i==2)&(Nc == 2)
        temp1 = txt(temp);
        temp1(indz) =  max(temp1(:))+1;
        NewMat = temp1-repmat(min(temp1')',[1 size(temp1,2)]);
        NewMat(indz) = 0;
        a = logical(sum(logical(NewMat)')');
        indc = find(a);
        Surf.SurfData.FaceVertexCData(indc,:) = ColorsV(indc,:);
    elseif (i>1)&(Nc > 2)
        temp1 = txt(temp);
        temp1(indz) =  max(temp1(:))+1;
        NewMat = temp1-repmat(min(temp1')',[1 size(temp1,2)]);
        NewMat(indz) = 0;
        a = logical(sum(logical(NewMat)')');
        indc = find(a);
        inddl = find(ismember(indc,indt));
        indc(inddl) = [];
        indrc = find(sum(logical(Surf.SurfData.FaceVertexCData(indc,:) - repmat(col(i-1,:),[length(indc) 1]))')==0);
        if ~isempty(indrc)
            Surf.SurfData.FaceVertexCData(indc,:) = repmat(col(i-1,:),[length(indc) 1]);
            Surf.SurfData.FaceVertexCData(indc(indrc),:) = repmat([1 1 1] - col(i-1,:),[length(indrc) 1]);
        else
            Surf.SurfData.FaceVertexCData(indc,:) = repmat(col(i-1,:),[length(indc) 1]);
        end
        indt = [indt;indc];
        clear temp1;
    end
    %% ===================== End of Painting Surface ======================= %%
end
%%  Ploting Surfaces
if Nc == 2
    Surft =Surf;
    Surft.SurfData.FaceVertexCData = ColorsV;
    Surfa{1,1}  = Surf;
    Surfa{2,1}  = Surft;
    Plot_Surf(Surfa,[1 0.2],sw,cl);
else
    Plot_Surf(Surf,tr,sw,cl);
end

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




