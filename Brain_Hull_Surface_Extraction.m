function [Surfl,Surfr] = Brain_Hull_Surface_Extraction(Image);
%
% Syntax :
%   [Surfl,Surfr] = Create_Region_ConvexHull(Image);
%
% This script creates hemispheric convex hull from freesurfer parcellation
% atlas. After, each convex hull is labelled according to freesurfer
% cortical segmentation protocol. 
%
% Input Parameters:
%     Image           :  Individual Aparc+Aseg filename
%
% Output Parameters:
%     Surfl, Surfr    :  Labelled Hemispheric Convex Hull Surfaces
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
% May 12th 2012
% Version $1.0
Image = '/media/Data/PROCESSING_RESULTS/PEPS/5-freesurfer_processing/0037-20100326/tmp/aparc+aseg.nii';
%TempI = '/media/Data/PROCESSING_RESULTS/PEPS/5-freesurfer_processing/0002-20090515/tmp/aparc+aseg+3mm_temp.nii'
% SurfFile = '/media/Data/PROCESSING_RESULTS/PEPS/5-freesurfer_processing/0002-20090515/surf/lh.pial';
% ChFiles = '/media/Data/PROCESSING_RESULTS/PEPS/5-freesurfer_processing/0002-20090515/label/lh.aparc.annot';
% [OutFiles, SurfF] = Exp_Surf(SurfFile, '0', '','', 'imp','n'); Surfp= SurfF{1}; % Reading Surface
% 
% talfile = ['/media/Data/PROCESSING_RESULTS/PEPS/5-freesurfer_processing/0002-20090515' filesep 'mri' filesep 'transforms' filesep 'talairach.lta'];
% cras1 = textread(talfile,'%s',5,'headerlines',20);
% cras = char(cras1);
% cras = [str2num(cras(3,:))  str2num(cras(4,:)) str2num(cras(5,:))];
%  Surfp.SurfData.vertices = Surfp.SurfData.vertices + repmat(cras,[size(Surfp.SurfData.vertices,1) 1])

V = spm_vol(Image);
% vertvox = (inv(V.mat)*[Surfp.SurfData.vertices ones(size(Surfp.SurfData.vertices,1),1)]')';
% Surfp.SurfData.vertices = vertvox(:,1:3);
% Surfp.Is = read_cfiles(ChFiles);


%% ==================== ENTRY VERIFICATION ===============================%
atype = 'aparc+aseg';
switch atype
    case 'aparc+aseg'
        [GMcodes,LHMCodes,RHMCodes,Names] = Gray_Matter_codes(atype);
        % Selecting Subcortical regions
        SubC = sort([10:13 17:18 26 49:54 58]);
        SubCL = sort([10:13 17:18 26]);
        SubCR = sort([49:54 58]);
        
        % Selecting Frontal regions
        FroIds = sort([1028 1003 1027 1018 1019 1020 1012 1014 1024 1017 1032]);
        FroIdsR = FroIds+1000;
        
        % Selecting Temporal regions
        TempIds = sort([1009 1015 1030 1001 1007 1034 1006 1033 1016]);
        TempIdsR = TempIds+1000;
        
        % Selecting Parietal regions
        ParIds = sort([1029 1008 1031 1022 1025]);
        ParIdsR = ParIds+1000;
        
        % Selecting Occipital regions
        OccIds = sort([1011 1013 1005 1021]);
        OccIdsR = OccIds+1000;
        
        % Selecting Cingulate regions
        CingIds = sort([1002 1010 1023 1026]);
        CingIdsR = CingIds+1000;
        
        % Selecting Insula regions
        InsIds = [1035];
        InsIdsR = [2035];
end
V = spm_vol(Image);
I = int16(spm_read_vols(V));
ind = find(ismember(I,[15 16 7 8 46 47]'));
I(ind) = 0; 



%% ===================== Creating Hemisphere Masks =======================%
IA = int16(zeros(size(I)));
ind = find(ismember(I,[ 4 28 5 30 31 SubCL FroIds TempIds ParIds OccIds 1000  1002 1010 1023 1026 1035 ]'));
IA(ind) = 1; % Left Hemisphere==1 Gray Matter, Avoiding Freesurfer errors in white Matter Parcelation
% White Matter
ind = find((I == 5001)|(I == 2));
It = IA*0;
It(ind) = 1;
bw = bwlabeln(It);
ind = find(It);
c = accumarray(bw(ind),ones(length(bw(ind)),1));
indpos = find(c == max(c));
It(bw~= indpos) = 0;clear bw;
IA(logical(It)) = 1;
clear It;

%----------------------------------
ind = find(ismember(I,[43 44 60 62 63 SubCR FroIdsR TempIdsR ParIdsR OccIdsR 2000  2002 2010 2023 2026 2035 ]'));
IA(ind) = 2; % Right Hemisphere==2 Gray Matter, Avoiding Freesurfer errors in white Matter Parcelation
% White Matter
ind = find((I == 5002)|(I == 41));
It = IA*0;
It(ind) = 2;
bw = bwlabeln(It);
ind = find(It);
c = accumarray(bw(ind),ones(length(bw(ind)),1));
indpos = find(c == max(c));
It(bw~= indpos) = 0;clear bw;
IA(logical(It)) = 2;
IA = int16(IA);
clear It;
%----------------------------------
%% ===================================================================== %%
%% ===================== Creating Brain Masks ============================%
Iseg = logical(I);
ind = find(ismember(I,[15 16 7 8 46 47]'));
Iseg(ind) = 0; 
%% ===================================================================== %%
%% Mask Refilling 
for z = 1:size(Iseg,3)
    Temp = squeeze(Iseg(:,:,z));
    Temp = imfill(Temp,'holes');
    Iseg(:,:,z) = Temp;
end;
for z = 1:size(Iseg,1)
    Temp = squeeze(Iseg(z,:,:));
    Temp = imfill(Temp,'holes');
    Iseg(z,:,:) = Temp;
end;
for z = 1:size(Iseg,2)
    Temp = squeeze(Iseg(:,z,:));
    Temp = imfill(Temp,'holes');
    Iseg(:,z,:) = Temp;
end;
bw = bwlabeln(Iseg);
ind = find(Iseg);
c = accumarray(bw(ind),ones(length(bw(ind)),1));
indpos = find(c == max(c));
Iseg(bw~= indpos) = 0;clear bw;
%% ===================================================================== %%

%% Closing Operation
 for i =1:1 % Number of closing operations
     Iseg = imerode(imdilate(Iseg,strel(ones(3,3,3))),strel(ones(3,3,3)));
 end
%% ===================================================================== %%
%% Mask Refilling
for z = 1:size(Iseg,3)
    Temp = squeeze(Iseg(:,:,z));
    Temp = imfill(Temp,'holes');
    Iseg(:,:,z) = Temp;
end;
for z = 1:size(Iseg,1)
    Temp = squeeze(Iseg(z,:,:));
    Temp = imfill(Temp,'holes');
    Iseg(z,:,:) = Temp;
end;
for z = 1:size(Iseg,2)
    Temp = squeeze(Iseg(:,z,:));
    Temp = imfill(Temp,'holes');
    Iseg(:,z,:) = Temp;
end;
bw = bwlabeln(Iseg);
ind = find(Iseg);
c = accumarray(bw(ind),ones(length(bw(ind)),1));
indpos = find(c == max(c));
Iseg(bw~= indpos) = 0;clear bw Temp;
%% Closing Operation
%% Closing Operation
 for i =1:10 % Number of closing operations
     Iseg = imerode(imdilate(Iseg,strel(ones(3,3,3))),strel(ones(3,3,3)));
 end
%% ===================================================================== %%
Iseg = Iso_Rem(Iseg,10); % Removing Isolated Points
[IA] = Atlas_Corr_custom(Iseg,IA); % Filling CSF Parts
Iseg = logical(Iseg);
IA = int16(IA);
% %% Removing Isolated points from Hemisphere mask
% IT = IA*0;
% % ---------- Left Hemisphere --------------
% It = IA*0;
% It(IA==1) = 1;
% IT(find(It)) = 1;
% % ---------- Right Hemisphere -------------
% It = IA*0;
% It(IA==2) = 1;
% IT(find(It)) = 2;
% IA = IT; clear IT
% %% ===================================================================== %%
% [IA] = Atlas_Corr_custom(Iseg,IA);
%% Extracting Left Hemisphere Convex Hull
It = IA*0;
ind = find(IA == 1);
It(ind) = 1;

IT = IA*0;
ind = find(ismember(I,[ 4 28 5 30 31 SubCL FroIds TempIds ParIds OccIds 1000  1002 1010 1023 1026 1035 250 251 252 253 254 255]'));
IT(ind) = I(ind);
[IT] = Atlas_Corr_custom(It,IT);
clear Iseg
%% ============== Refilling remaining sulcal space: Z dimension ======== %%
disp(['Left Hemisphere: Refilling Sulcal Space ... Z']);
tic;
[X,Y] = meshgrid(1:size(IT,1),1:size(IT,2));X =X(:);Y =Y(:);
for i = 1:size(IT,3)
    A = IT(:,:,i);
    T = A*0;
    if sum(A(:))
        [B,L] = bwboundaries(logical(A),'noholes');
        for j = 1:length(B)
            boundary = B{j};
            k=LineCurvature2D(boundary);
            ind = find(k >0);
            boundary(ind,:) = [];
            indold = 0;
            while (isequal(ind,indold) ==0)|~isempty(ind)
                indold = ind;
                k=LineCurvature2D(boundary);
                indid = sub2ind(size(A),boundary(:,1),boundary(:,2));
                indr = ismember(A(indid),[1001 1014 1000 1006 1016 1021 7000 1032 1033 1034 1035]);
                ind = find((k >0)&indr==0);
                boundary(ind,:) = [];
            end
            boundary = [boundary;boundary(1,:)];
            indin = inpolygon(Y,X,boundary(:,2),boundary(:,1));
            indin = sub2ind(size(T),X(find(indin)),Y(find(indin)));
            T(indin) = j;
        end
    end
    Temp(:,:,i) = T;
end

% -------------- Convex Hull Lobar Parcelation Left Hemisphere -----------%
IT = uint16(zeros(size(IA)));
ind = find(ismember(I,[ 4 28 5 30 31 SubCL FroIds TempIds ParIds OccIds 1000  1002 1010 1023 1026 1035 250 251 252 253 254 255]'));
IT(ind) = I(ind);
[IT] = int16(Atlas_Corr_custom(logical(Temp),IT));
clear Temp;
toc;
%% =============== End of Refilling remaining sulcal space ============= %%

%% ============== Refilling remaining sulcal space: Y dimension ======== %%
disp(['Left Hemisphere: Refilling Sulcal Space ... Y']);
tic;
[X,Y] = meshgrid(1:size(IT,1),1:size(IT,3));X =X(:);Y =Y(:);
for i = 1:size(IT,2)
    A = squeeze(IT(:,i,:));
    T = A*0;
    if sum(A(:))
        [B,L] = bwboundaries(logical(A),'noholes');
        for j = 1:length(B)
            boundary = B{j};
            k=LineCurvature2D(boundary);
            ind = find(k >0);
            boundary(ind,:) = [];
            indold = 0;
            while (isequal(ind,indold) ==0)|~isempty(ind)
                indold = ind;
                k=LineCurvature2D(boundary);
                indid = sub2ind(size(A),boundary(:,1),boundary(:,2));
                indr = ismember(A(indid),[1001 1014 1000 1006 1016 1021 7000 1032 1033 1034 1035]);
                ind = find((k >0)&indr==0);
                boundary(ind,:) = [];
            end
            boundary = [boundary;boundary(1,:)];
            indin = inpolygon(Y,X,boundary(:,2),boundary(:,1));
            indin = sub2ind(size(T),X(find(indin)),Y(find(indin)));
            T(indin) = j;
        end
    end
    Temp(:,i,:) = T;
end

% -------------- Convex Hull Lobar Parcelation Left Hemisphere -----------%
[IT] = int16(Atlas_Corr_custom(logical(Temp),IT));
clear Temp;
toc
%% =============== End of Refilling remaining sulcal space ============= %%

% -------------- Convex Hull Lobar Parcelation Left Hemisphere -----------%
%% =============== End of Refilling remaining sulcal space ============= %%

ind = find(ismember(I,[ 4 28 5 30 31 SubCL FroIds TempIds ParIds OccIds 1000  1002 1010 1023 1026 1035 250 251 252 253 254 255 5001 2]'));
Itemp = I*0;
Itemp(ind) = I(ind);
IsL = logical(IT*0);
IsL(ind) = 1;
IsL = logical(IT) - IsL;
indr = find(IsL <0);
IT(indr) = 0;
IsL(indr) = 0;
[fv] = Ver_Ext_vox(V,8000,logical(IT));
voxs = (V.mat*[fv.vertices(:,1) fv.vertices(:,2) fv.vertices(:,3) ones(size(fv.vertices(:,1),1),1)]')';
Surfl.SurfData.vertices = [voxs(:,1) voxs(:,2) voxs(:,3)];
Surfl.SurfData.faces = fv.faces;

clear Temp T;
%% ================= Convex Hull Region Parcelation =======================%

%% ===================================================================== %%
%% Extracting Right Hemisphere Convex Hull
It = IA*0;
ind = find(IA == 2);
It(ind) = 1;

IT = IA*0;
T = I(:);
ind = ismember(T,[ 4 28 5 30 31 SubCR FroIdsR TempIdsR ParIdsR OccIdsR 2000  2002 2010 2023 2026 2035 250 251 252 253 254 255]');
IT(ind) = I(ind);
[IT] = Atlas_Corr_custom(It,IT);

%% ============== Refilling remaining sulcal space: Z dimension ======== %%
disp(['Left Hemisphere: Refilling Sulcal Space ... Z']);
tic;
[X,Y] = meshgrid(1:size(IT,1),1:size(IT,2));X =X(:);Y =Y(:);
for i = 1:size(IT,3)
    A = IT(:,:,i);
    T = A*0;
    if sum(A(:))
        [B,L] = bwboundaries(logical(A),'noholes');
        for j = 1:length(B)
            boundary = B{j};
            k=LineCurvature2D(boundary);
            ind = find(k >0);
            boundary(ind,:) = [];
            indold = 0;
            while (isequal(ind,indold) ==0)|~isempty(ind)
                indold = ind;
                k=LineCurvature2D(boundary);
                indid = sub2ind(size(A),boundary(:,1),boundary(:,2));
                indr = ismember(A(indid),[2001 2014 2000 2006 2016 2021 7000 2032 2033 2034 2035]);
                ind = find((k >0)&indr==0);
                boundary(ind,:) = [];
            end
            boundary = [boundary;boundary(1,:)];
            indin = inpolygon(Y,X,boundary(:,2),boundary(:,1));
            indin = sub2ind(size(T),X(find(indin)),Y(find(indin)));
            T(indin) = j;
        end
    end
    Temp(:,:,i) = T;
end

% -------------- Convex Hull Lobar Parcelation Left Hemisphere -----------%
IT = IA*0;
T = I(:);
ind = ismember(T,[ 4 28 5 30 31 SubCR FroIdsR TempIdsR ParIdsR OccIdsR 2000  2002 2010 2023 2026 2035 250 251 252 253 254 255]');
IT(ind) = I(ind);
[IT] = Atlas_Corr_custom(logical(Temp),IT);
toc;
%% =============== End of Refilling remaining sulcal space ============= %%

%% ============== Refilling remaining sulcal space: Y dimension ======== %%
disp(['Left Hemisphere: Refilling Sulcal Space ... Y']);
tic;
[X,Y] = meshgrid(1:size(IT,1),1:size(IT,3));X =X(:);Y =Y(:);
for i = 1:size(IT,2)
    A = squeeze(IT(:,i,:));
    T = A*0;
    if sum(A(:))
        [B,L] = bwboundaries(logical(A),'noholes');
        for j = 1:length(B)
            boundary = B{j};
            k=LineCurvature2D(boundary);
            ind = find(k >0);
            boundary(ind,:) = [];
            indold = 0;
            while (isequal(ind,indold) ==0)|~isempty(ind)
                indold = ind;
                k=LineCurvature2D(boundary);
                indid = sub2ind(size(A),boundary(:,1),boundary(:,2));
                indr = ismember(A(indid),[2001 2014 2000 2006 2016 2021 7000 2032 2033 2034 2035]);
                ind = find((k >0)&indr==0);
                boundary(ind,:) = [];
            end
            boundary = [boundary;boundary(1,:)];
            indin = inpolygon(Y,X,boundary(:,2),boundary(:,1));
            indin = sub2ind(size(T),X(find(indin)),Y(find(indin)));
            T(indin) = j;
        end
    end
    Temp(:,i,:) = T;
end

% -------------- Convex Hull Lobar Parcelation Left Hemisphere -----------%
[IT] = int16(Atlas_Corr_custom(logical(Temp),IT));
toc
%% =============== End of Refilling remaining sulcal space ============= %%

ind = find(ismember(I,[ 4 28 5 30 31 SubCR FroIdsR TempIdsR ParIdsR OccIdsR 2000  2002 2010 2023 2026 2035 250 251 252 253 254 255 41 5002]'));
Itemp = I*0;
Itemp(ind) = I(ind);
IsL = logical(IT*0);
IsL(ind) = 1;
IsL = logical(IT) - IsL;
indr = find(IsL <0);
IT(indr) = 0;
IsL(indr) = 0;
[Itemp] = int16(Atlas_Corr_custom(IsL,Itemp));
ind = find(Itemp);
IT(ind) = Itemp(ind); clear IsL Itemp;
[fv] = Ver_Ext_vox(V,8000,logical(IT));
voxs = (V.mat*[fv.vertices(:,1) fv.vertices(:,2) fv.vertices(:,3) ones(size(fv.vertices(:,1),1),1)]')';
Surfr.SurfData.vertices = [voxs(:,1) voxs(:,2) voxs(:,3)];
Surfr.SurfData.faces = fv.faces;

clear Temp T;
%% ================= Convex Hull Region Parcelation =======================%


return

function [labid] = Recur_Corr(Surf,lab,labid,cont);
%
% Syntax :
% [labid] = Recur_Corr(Surf,lab,labid,cont);
%
% Recursive function for labelling correction.    
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

ind = find(Surf.Is ==lab&(labid ==0));
A = unique(Surf.SurfData.faces(Surf.Tri(ind(1),3:2+Surf.Tri(ind(1),2)),:));
indt = find(Surf.Is(A)~= lab);
A(indt) = [];A(A == ind(1))= [];
T = unique([ind(1); A(:)]);
labid(T) = cont;
An = rand(size(A));
while sum(A)~=sum(An)
    An = A;
    Neib = Surf.Tri(A,3:end); Neib = unique(nonzeros(Neib(:)));
    A = unique(Surf.SurfData.faces(Neib,:));
    indt = find(Surf.Is(A)~= lab);
    A(indt) = [];
    labid(A) = cont;
    T =unique([T;A(:)]);
end
indn = find(Surf.Is ==lab&(labid ==0));
if ~isempty(indn)
    cont = cont+1;
    [labid] = Recur_Corr(Surf,lab,labid,cont);
else
    return;
end
if ~isempty(A)
    for i = 1:size(A,2)
        TA = unique(Surf.SurfData.faces(Surf.Tri(A(i),3:2+Surf.Tri(A(i),2)),:));
        indt = find(Surf.Is(TA(i))~= lab);
        TA(indt) = [] ;
        T = [T; A(:)];
    end
else
    return;
end



%==========================================================================

function I = Iso_Rem(T,Nhood);
%
%This function removes isolated points from mask. 
%
% Input Parameters:
%   T            : Mask
%   Nhood        : Minimun number of neighbors. 
% Output Parameters:
%   I            : Mask without isolated points  
%__________________________________________________________________________
% Authors:  Yasser Alem?n G?mez 
% Neuroimaging Department
% Cuban Neuroscience Center
% Last update: November 15th 2005
% Version $1.0

%=========================Main program====================================%  
warning off
if length(size(T))==3
    I = zeros(size(T)+2);
    I(2:end-1,2:end-1,2:end-1) = T;
    clear T
    ind = find(I>0);
    [x,y,z] = ind2sub(size(I), ind);
    s = size(x,1);
    sROI = zeros(size(I));
    [X, Y, Z] = meshgrid(-1:1,-1:1,-1:1);
    X = X(:);Y = Y(:);Z = Z(:);
    Neib = [X Y Z];clear X Y Z;
    pos = find((Neib(:,1)==0)&(Neib(:,2)==0)&(Neib(:,3)==0));
    Neib(pos,:) = [];
    for i =1:26
        M = Neib(i,:);
        S = [x y z] + M(ones(s,1),:);
        ind2 = sub2ind(size(I),S(:,1),S(:,2),S(:,3));
        sROI(ind) = sROI(ind) + I(ind2);
    end
    ind = find(sROI<Nhood);
    I(ind) =0;
    I = I(2:end-1,2:end-1,2:end-1);
elseif length(size(T))==2
    I = zeros(size(T)+2);
    I(2:end-1,2:end-1) = T;
    clear T
    ind = find(I>0);
    [x,y] = ind2sub(size(I), ind);
    s = size(x,1);
    sROI = zeros(size(I));
    [X, Y] = meshgrid(-1:1,-1:1);
    X = X(:);Y = Y(:);
    Neib = [X Y];clear X Y;
    pos = find((Neib(:,1)==0)&(Neib(:,2)==0));
    Neib(pos,:) = [];
    for i =1:8
        M = Neib(i,:);
        S = [x y] + M(ones(s,1),:);
        ind2 = sub2ind(size(I),S(:,1),S(:,2));
        sROI(ind) = sROI(ind) + I(ind2);
    end
    ind = find(sROI<Nhood);
    I(ind) =0;
    I = I(2:end-1,2:end-1);
end

%========================End of main program==============================%
return;

function [GMcodes,LHMCodes,RHMCodes,Names] = Gray_Matter_codes(atlastype);
%
% Syntax :
% [GMcodes] = Gray_Matter_codes(atlastype);
%
% This function extract gray matter codes from a specified atlas.
%
% Input Parameters:
%   atlastype     : Atlas type.
%
%
% Output Parameters:
%   GMcodes       : Gray matter codes.
%
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
% March 20th 2012
% Version $1.0


switch atlastype
    case 'aparc+aseg'
        [a,temp] = system('echo $FREESURFER_HOME');
        txtfile =  strcat(temp,filesep,'FreeSurferColorLUT.txt');
        fid = fopen(txtfile);
        cont = 0;lines = '';
        while 1
            cont = cont + 1;
            line = fgetl(fid);
            if ~ischar(line),   break,   end
            lines = strvcat(lines,line);
        end
        ind = find(sum(isspace(lines'))==size(lines,2));lines(ind,:) = [];
        Codes = '';
        Names = '';
        for i = 1:size(lines,1)
            if ~strcmp(lower(deblank(lines(i,1))),'#');
                temp = strread(deblank(lines(i,:)),'%s','delimiter',' ');
                code = temp{1};
                name = temp{2};
                Codes = strvcat(Codes,code);
                Names = strvcat(Names,name);
            end
        end
        fclose(fid);
        GMcodes = [10:13 17:18 26 49:54 58 1001:1003 1005:1035 2001:2003 2005:2035];
        LHMCodes = [10:13 17:18 26 1001:1003 1005:1035];
        RHMCodes = [49:54 58 2001:2003 2005:2035];
        index = ismember(str2num(Codes),GMcodes);
        Names = Names(index,:);
    case 'a2009s+aseg'
        [a,temp] = system('echo $FREESURFER_HOME');
        txtfile =  strcat(temp,filesep,'FreeSurferColorLUT.txt');
        fid = fopen(txtfile);
        cont = 0;lines = '';
        while 1
            cont = cont + 1;
            line = fgetl(fid);
            if ~ischar(line),   break,   end
            lines = strvcat(lines,line);
        end
        ind = find(sum(isspace(lines'))==size(lines,2));lines(ind,:) = [];
        Codes = '';
        Names = '';
        for i = 1:size(lines,1)
            if ~strcmp(lower(deblank(lines(i,1))),'#');
                temp = strread(deblank(lines(i,:)),'%s','delimiter',' ');
                code = temp{1};
                name = temp{2};
                Codes = strvcat(Codes,code);
                Names = strvcat(Names,name);
            end
        end
        fclose(fid);
        GMcodes = [9:13 17:18 26 48:54 58 11101:11175 12101:12175];
        index = ismember(str2num(Codes),GMcodes);
        Names = Names(index,:);
    case 'a2005s+aseg'
        [a,temp] = system('echo $FREESURFER_HOME');
        txtfile =  strcat(temp,filesep,'FreeSurferColorLUT.txt');
        fid = fopen(txtfile);
        cont = 0;lines = '';
        while 1
            cont = cont + 1;
            line = fgetl(fid);
            if ~ischar(line),   break,   end
            lines = strvcat(lines,line);
        end
        ind = find(sum(isspace(lines'))==size(lines,2));lines(ind,:) = [];
        Codes = '';
        Names = '';
        for i = 1:size(lines,1)
            if ~strcmp(lower(deblank(lines(i,1))),'#');
                temp = strread(deblank(lines(i,:)),'%s','delimiter',' ');
                code = temp{1};
                name = temp{2};
                Codes = strvcat(Codes,code);
                Names = strvcat(Names,name);
            end
        end
        fclose(fid);
        GMcodes = [9:13 17:18 26 48:54 58 1102:1181 2102:2181];
        index = ismember(str2num(Codes),GMcodes);
        Names = Names(index,:);
    case 'yeoh7networks'
        
    case 'yeoh17networks'
        
    case 'parckmeans'
        
    case 'ibaspm116'
        txtfile = which('atlas116.cod');
        [Codes,b] = textread(txtfile,'%u%s','delimiter','=');
        %         Codes = str2num(char(a));
        Names = char(b);
        GMcodes = [1:90];
        index = ismember(Codes,GMcodes);
        Names = Names(index,:);
    case 'ibaspm71'
        txtfile = which('atlas71.cod');
        [Codes,b] = textread(txtfile,'%u%s','delimiter','=');
        %         Codes = str2num(char(a));
        Names = char(b);
        GMcodes = [1     2     4     5     6     7     9    10    11    12    14    15    16    18    19    20 ...
            23    25    26    27    32    33    36    37    38    39    41    50    52    53    54    56    60 ...
            61    62    63    64    67    69    70    72    74    75    76    80    85    88    90    97    98 ...
            99   101   102   108   110   112   114   119   125   130   132   140   145   154   159   164   165 ...
            175   196   203   251];
        index = ismember(Codes,GMcodes);
        Names = Names(index,:);
end
return

function [IA] = Atlas_Corr_custom(Iseg,IA);
%
% Syntax : 
% [IA] = Atlas_Corr_custom(Iseg,IA);
%
% This function corrects the erros ocurred during the atlasing process. 
% It assings a label to the points of gray matter mask that didn't get any 
% during the labelling step.   
%
% Input Parameters:
%  Iseg           : Gray Matter Mask.
%  IA             : Individual Atlas without correction.
%  
% Output Parameter:
%  IA             : Corrected Individual Gray Matter Atlas.
%__________________________________________________________________________
% Authors:  Yasser Aleman Gomez
% Neuroimaging Department
% Cuban Neuroscience Center
% Last update: November 15th 2006
% Version $1.0
Isegt = zeros(size(Iseg)+2);Isegt(2:end-1,2:end-1,2:end-1) = Iseg;Iseg = Isegt; clear Isegt;
IAt = zeros(size(IA)+2);IAt(2:end-1,2:end-1,2:end-1) = IA;IA = IAt;clear IAt;
ind = find((Iseg==1)&(IA==0));
[X, Y, Z] = meshgrid(-1:1,-1:1,-1:1);
X = int16(X(:));Y = int16(Y(:));Z = int16(Z(:));
Neib = [X Y Z];clear X Y Z;l = 1;indr = 0;
%Neib(find(sum(logical(Neib')) >1),:) = [];

pos = find(ismember(Neib,[0 0 0],'rows'));
Neib(pos,:) = [];
while ~isempty(ind)&(l>0)
    [x,y,z] = ind2sub(size(Iseg), ind);
    pos = find(ismember(Neib,[0 0 0],'rows'));
    Neib(pos,:) = [];
    Temp = zeros(size(ind,1),size(Neib,1)+2);
    Temp(:,1) = ind;
    for j = 1:size(Neib,1)
        S = repmat(double(Neib(j,:)),[size(Temp(:,1),1) 1]) + [x y z];
        indtmp = sub2ind(size(IA),S(:,1),S(:,2),S(:,3));
        Temp(:,j+2) = IA(indtmp);
    end
    Temp(:,2) = sum(logical(Temp(:,3:end))')';
    [xt,yt] = sort(Temp(:,2),'descend');indtemp = find(xt>1);
    Temp1 = Temp(yt(indtemp),:);clear Temp; Temp = Temp1; clear Temp1;
    for j = 1:length(Temp(:,1))
        T = Temp(j,3:end);indp =find(T==0);T(indp) = [];
        c = accumarray(T',ones(length(T),1));
        indpos = find(c == max(c));
        if length(indpos)>1
            [xo,yo,zo] = ind2sub(size(Iseg),Temp(j,1));
            TNeib = Neib;TNeib(indp,:) = [];
            S = repmat(int32([xo yo zo]),[size(TNeib,1) 1]) + int32(TNeib);
            ind2 = sub2ind(size(Iseg),S(:,1),S(:,2),S(:,3));
            temp = ismember(IA(ind2),indpos);
            Vect = [xo yo zo;S(temp,1) S(temp,2) S(temp,3)];
            D = dist(Vect');D = D(1,2:end);
            indt = find(D == min(D));
            indpos = IA(Vect(indt(1)+1,1),Vect(indt(1)+1,2),Vect(indt(1)+1,3));
        end
        IA(Temp(j,1)) = indpos;
    end
    clear x y z X Y Z S xo yo zo ind2 indt indpos Vect T;
    indr= find((Iseg==1)&(IA==0));
    l = size(Temp,1);
    ind = indr; 
end
IAt = IA(2:end-1,2:end-1,2:end-1);
IAt = Iseg(2:end-1,2:end-1,2:end-1).*IAt;
IA = IAt; clear IAt;
return;
