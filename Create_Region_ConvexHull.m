function [Surfl,Surfr] = Create_Region_ConvexHull(Image);
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
%Image = '/media/Data/PROCESSING_RESULTS/PEPS/5-freesurfer_processing/0037-20100326/tmp/aparc+aseg.nii';
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


% Sagital Slice
temp = sum(sum(I,3),2);
ind = find(temp == max(temp));
xslice = [ind(1) 0 0];
indx = find(temp);
Xmin = min(indx);
Xmax = max(indx);

% Coronal Slice
temp = sum(sum(I,1),3);
ind = find(temp == max(temp));
yslice = [0 ind(1)  0];
indy = find(temp);
Ymin = min(indy);
Ymax = max(indy);

% Axial Slice
temp = sum(sum(I,1),2);
ind = find(temp == max(temp));
zslice = [0 0 ind(1)];
indz = find(temp);
Zmin = min(indz);
Zmax = max(indz);

I = I(Xmin-1:Xmax+1,Ymin-1:Ymax+1,Zmin-1:Zmax+1);


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

% %% ============== Refilling remaining sulcal space: X dimension ======== %%
% disp(['Left Hemisphere: Refilling Sulcal Space ... X']);
% tic;
% [X,Y] = meshgrid(1:size(IT,2),1:size(IT,3));X =X(:);Y =Y(:);
% for i = 1:size(IT,1)
%     A = squeeze(IT(i,:,:));
%     T = A*0;
%     if sum(A(:))
%         [B,L] = bwboundaries(logical(A),'noholes');
%         for j = 1:length(B)
%             boundary = B{j};
%             k=LineCurvature2D(boundary);
%             ind = find(k >0);
%             boundary(ind,:) = [];
%             indold = 0;
%             while (isequal(ind,indold) ==0)|~isempty(ind)
%                 indold = ind;
%                 k=LineCurvature2D(boundary);
%                 indid = sub2ind(size(A),boundary(:,1),boundary(:,2));
%                 indr = ismember(A(indid),[1001 1014 1000 1006 1016 1021 7000 1032 1033 1034 1035]);
%                 ind = find((k >0)&indr==0);
%                 boundary(ind,:) = [];
%             end
%             boundary = [boundary;boundary(1,:)];
%            indin = inpolygon(Y,X,boundary(:,2),boundary(:,1));
%            indin = sub2ind(size(T),X(find(indin)),Y(find(indin)));
%            T(indin) = j;
%         end
%     end
%     Temp(i,:,:) = T;
% end
% toc

% -------------- Convex Hull Lobar Parcelation Left Hemisphere -----------%
%[IT] = Atlas_Corr_custom(logical(Temp),IT);
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
[Itemp] = int16(Atlas_Corr_custom(IsL,Itemp));
ind = find(Itemp);
IT(ind) = Itemp(ind); clear IsL Itemp;
[Surf] = Surf_Extraction_custom(V,logical(IT));
Surf.SurfData.vertices = Surf.SurfData.vertices - 2;
clear Temp T;
%[Cmean,Cgaussian,Dir1,Dir2,Lambda1,Lambda2]=patchcurvature(Surf.SurfData,true);
Surfa(1) = Surf;
%% ================= Convex Hull Region Parcelation =======================%

%% Creating Left Regions
ind = find(ismember(IT,[250 251 252 253 254 255 1000]')); IT(ind) = 7000; % Corpus Callosum
ind = find(ismember(IT,[4 28 5 30 31 SubCL ]')); IT(ind) = 7000; % Subcortical Structures, Ventricles and CSF
%%
vert1 = Surfa(1).SurfData.vertices;
Is = spm_sample_vol(IT,double(vert1(:,1)),double(vert1(:,2)),double(vert1(:,3)),0);
while length(find(Is == 0)) > 1000
    [IT] = Atlas_Corr_custom(imdilate(logical(IT),strel(ones(3,3,3))),IT);
    Is = spm_sample_vol(IT,double(vert1(:,1)),double(vert1(:,2)),double(vert1(:,3)),0);
end
Surfa(1).Is = Is;
Surfa(1) = Surf_Corr(Surfa(1));

%Surfa = PaintSurf_with_FreeSurferColors(Surfa);

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

% %% ============== Refilling remaining sulcal space: X dimension ======== %%
% disp(['Left Hemisphere: Refilling Sulcal Space ... X']);
% tic;
% [X,Y] = meshgrid(1:size(IT,2),1:size(IT,3));X =X(:);Y =Y(:);
% for i = 1:size(IT,1)
%     A = squeeze(IT(i,:,:));
%     T = A*0;
%     if sum(A(:))
%         [B,L] = bwboundaries(logical(A),'noholes');
%         for j = 1:length(B)
%             boundary = B{j};
%             k=LineCurvature2D(boundary);
%             ind = find(k >0);
%             boundary(ind,:) = [];
%             indold = 0;
%             while (isequal(ind,indold) ==0)|~isempty(ind)
%                 indold = ind;
%                 k=LineCurvature2D(boundary);
%                 indid = sub2ind(size(A),boundary(:,1),boundary(:,2));
%                 indr = ismember(A(indid),[2001 2014 2000 2006 7000 2035 2033 28]);
%                 ind = find((k >0)&indr==0);
%                 boundary(ind,:) = [];
%             end
%             boundary = [boundary;boundary(1,:)];
%            indin = inpolygon(Y,X,boundary(:,2),boundary(:,1));
%            indin = sub2ind(size(T),X(find(indin)),Y(find(indin)));
%            T(indin) = j;
%         end
%     end
%     Temp(i,:,:) = T;
% end
% toc

% -------------- Convex Hull Lobar Parcelation Left Hemisphere -----------%
%[IT] = Atlas_Corr_custom(logical(Temp),IT);
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
[Surf] = Surf_Extraction_custom(V,logical(IT));
Surf.SurfData.vertices = Surf.SurfData.vertices - 2;
clear Temp T;

%[Cmean,Cgaussian,Dir1,Dir2,Lambda1,Lambda2]=patchcurvature(Surf.SurfData,true);
Surf.Is = zeros(length(Surf.SurfData.vertices),1);
Surfa(2) = Surf;
%% ================= Convex Hull Region Parcelation =======================%

%% Creating Left Regions
ind = find(ismember(IT,[250 251 252 253 254 255 2000]')); IT(ind) = 7000; % Corpus Callosum
ind = find(ismember(IT,[4 28 5 30 31 SubCR ]')); IT(ind) = 7000; % Subcortical Structures, Ventricles and CSF
%%
vert1 = Surfa(2).SurfData.vertices;
Is = spm_sample_vol(IT,double(vert1(:,1)),double(vert1(:,2)),double(vert1(:,3)),0);
while length(find(Is == 0)) > 1000
    [IT] = Atlas_Corr_custom(imdilate(logical(IT),strel(ones(3,3,3))),IT);
    Is = spm_sample_vol(IT,double(vert1(:,1)),double(vert1(:,2)),double(vert1(:,3)),0);
end
Surfa(2).Is = Is;
Surfa(2) = Surf_Corr(Surfa(2));

%% ===================================================================== %%
%% Some Usefull Plots
% close all
% [vertices, faces] = freesurfer_read_surf('/media/COSAS/Test/Joost/LOBE_HULLPARC/lh.pial');
% Surfp.SurfData.vertices = vertices;Surfp.SurfData.faces = faces;Surfp.Name = 'Pial';
% cras1 = textread('/media/COSAS/Test/Joost/LOBE_HULLPARC/talairach.lta','%s',5,'headerlines',20);
% cras = char(cras1);cras = [str2num(cras(3,:))  str2num(cras(4,:)) str2num(cras(5,:))];
% voxs = (inv(V.mat)*[Surfp.SurfData.vertices(:,1)+cras(1) Surfp.SurfData.vertices(:,2)+cras(2) Surfp.SurfData.vertices(:,3)+cras(3) ones(size(Surfp.SurfData.vertices(:,1),1),1)]')';
% showcs3(I);
% hold on
% plot3(voxs(:,1),voxs(:,2),voxs(:,3),'.y');
% [vertices, faces] = freesurfer_read_surf('/media/COSAS/Test/Joost/LOBE_HULLPARC/rh.pial');
% Surfp.SurfData.vertices = vertices;Surfp.SurfData.faces = faces;Surfp.Name = 'Pial';
% cras1 = textread('/media/COSAS/Test/Joost/LOBE_HULLPARC/talairach.lta','%s',5,'headerlines',20);
% cras = char(cras1);cras = [str2num(cras(3,:))  str2num(cras(4,:)) str2num(cras(5,:))];
% voxs = (inv(V.mat)*[Surfp.SurfData.vertices(:,1)+cras(1) Surfp.SurfData.vertices(:,2)+cras(2) Surfp.SurfData.vertices(:,3)+cras(3) ones(size(Surfp.SurfData.vertices(:,1),1),1)]')';
% hold on
% plot3(voxs(:,1),voxs(:,2),voxs(:,3),'.b');
%% ===================================================================== %%

%% ======================= Area Computation ============================ %%


voxs = (V.mat*[Surfa(1).SurfData.vertices(:,1) Surfa(1).SurfData.vertices(:,2) Surfa(1).SurfData.vertices(:,3) ones(size(Surfa(1).SurfData.vertices(:,1),1),1)]')';
Surfa(1).SurfData.vertices = voxs(:,1:3);
[Surfl] = Struct_Area_Comp(Surfa(1));

voxs = (V.mat*[Surfa(2).SurfData.vertices(:,1) Surfa(2).SurfData.vertices(:,2) Surfa(2).SurfData.vertices(:,3) ones(size(Surfa(2).SurfData.vertices(:,1),1),1)]')';
Surfa(2).SurfData.vertices = voxs(:,1:3);
[Surfr] = Struct_Area_Comp(Surfa(2));

%% ===================================================================== %%


return

function [Surf] = Surf_Corr(Surf);
%
% Syntax :
% [Surf] = Surf_Corr(Surf);
%
% This function corrects labeled surfaces. It removes 
% small isolated points( or group of points) that are 
% surrounded by a different label.     
%
% Input Parameters:
%   Surf      : Individual Surfaces.
%
% Output Parameters:
%   Surf  : Corrected Surfaces.
%
% Related references:
%
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

% [fa,intt] =sort(double(Surf.SurfData.faces)');
% [ford,int2] =sortrows(fa');
% temp = ford(1:end-1,:) - ford(2:end,:);
% fnd = find((temp(:,1)== 0)&(temp(:,2)== 0)&(temp(:,3)== 0));
% Surf.SurfData.faces(int2(fnd),:)=[];
% Surf.SurfData.faces =uint32(Surf.SurfData.faces);Nfaces =size(Surf.SurfData.faces,1);
% [Tri] = Vert_Neib(double(Surf.SurfData.faces),size(Surf.SurfData.vertices,1),Nfaces);Surf.Tri=Tri;
% Temp = sum(Tri);
% Tri(:,Temp==0) = [];
strl = unique(Surf.Is);
if ~isfield(Surf,'Tri')
    Npoints = size(Surf.SurfData.vertices,1);
    Nfaces = size(Surf.SurfData.faces,1);
    [Tri] = Vert_Neib(double(Surf.SurfData.faces),Npoints,Nfaces);
    Temp = sum(Tri);
    Tri(:,Temp==0) = [];
    Surf.Tri = Tri; clear Tri;
end
lab = Surf.Is;
set(0,'RecursionLimit',1000);
for i = 1:size(strl,1)
    lab = strl(i);
    [labid] = Recur_Corr(Surf,lab,zeros(size(Surf.Is)),1);
    c = accumarray(nonzeros(labid),ones(length(nonzeros(labid)),1));
    indpos = find(c == max(c));
    if length(indpos)==1
        Surf.Is((labid~= indpos(1))&(labid~=0)) = 0;
    else
        Surf.Is(~(ismember(labid,indpos))&(labid~=0)) = 0;
    end
end
Is = Recur(Surf.Is, Surf,1);
Surf.Is = Is;
ind = find(Surf.Is ==0);
if ~isempty(ind)
    dfac = unique(nonzeros(Surf.Tri(ind,3:end)));
    Surf.SurfData.faces(dfac,:) = [];
    Surf.SurfData.vertices(ind,:) = [];
    if isfield(Surf.SurfData,'VertexNormals')
        Surf.SurfData.VertexNormals(ind,:) = [];
    end
    if isfield(Surf,'Is')
        Surf.Is(ind,:) = [];
    end
    if isfield(Surf.SurfData,'FaceVertexCData')
        Surf.SurfData.FaceVertexCData(ind,:) = [];
    end
    cont = 0;
    Mat = Surf.SurfData.faces;
    for i =1:size(ind,1)
        cont= cont+1;
        dvert = find(Surf.SurfData.faces >ind(i));
        Mat(dvert) = Surf.SurfData.faces(dvert) - cont;
    end
    Surf.SurfData.faces = Mat;
    Npoints = size(Surf.SurfData.vertices,1);
    Nfaces = size(Surf.SurfData.faces,1);
    [Tri] = Vert_Neib(double(Surf.SurfData.faces),Npoints,Nfaces);
    Temp = sum(Tri);
    Tri(:,Temp==0) = [];
    Surf.Tri = Tri;
end
return;

function Is = Recur(Is, Surf,level);
%
% Syntax :
% Is = Recur(Is, Surf,level);
%
% This function refills labeled surfaces that contains non-labeled points. It uses
% a recursive process to use neightboor labeled ponts information for non-labeled ones.
%
% Input Parameters:
%   Is       : Surfaces labels.
%   Surf     : Struct variable that contains surface information.
%   level    : Neightborhood level(How many points do we have to take into 
%              account for labelling correction).
%
% Output Parameters:
%  Is        : Corrected Surface labels
%
% Related references:
%
%
% See also: Red_Surf Smooth_Surf Plot_Surf Surf_Comp Atlas_Surf
% Plot_oversurf
%__________________________________________________
% Authors: Yasser Aleman Gomez
% Neuroimaging Department
% Cuban Neuroscience Center
% November 30th 2006
% Version $1.0

%=========================Main program====================================%
ind = find(Is == 0);
if ind~=0
    for j = 1: length(ind)
        ind1 =  Surf.Tri(ind(j),3:Surf.Tri(ind(j),2)+2);
        Vert = Surf.SurfData.faces(ind1,:);Vert = unique(Vert(:));ind1 = find(Vert ==ind(j));
        Vert(ind1) = [];
        ord(j) = length(find(Is(Vert) ~=0 ));
    end
    [x,y] = sort(ord,'descend');indtemp = find(x);
    indtemp1 = ind(y(indtemp));clear ind; ind = indtemp1;
    for j = 1:length(ind)
        ind1 =  Surf.Tri(ind(j),3:Surf.Tri(ind(j),2)+2);
        Vert = Surf.SurfData.faces(ind1,:);Vert = unique(Vert(:));ind1 = find(Vert == ind(j));
        Vert(ind1) = [];It = Is(Vert); indt = find(It~=0);
        Vert = [ind(j);Vert(indt)];
        D = dist(Surf.SurfData.vertices(Vert,:)'); D = D(1,2:end);
        if size(D,2)>=level
            indt = find(D == min(D));
            Is(ind(j)) = Is(Vert(indt(1)+1));
        else
            Is(ind(j)) = 0;
        end
    end
    indr = find(Is == 0);
    l = ismember(ind,indr);
    if sum(l(:)) ==size(ind,1)
        return;
    end
    if ~isempty(indr)
        Is = Recur(Is, Surf,level);
    end
end
%========================End of main program==============================%
return;


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

function [Surf] = Struct_Area_Comp(Surf);
%
% Syntax :
% [Surf] = Struct_Area_Comp(Surf);
%
% This function computes the parcelated structures surface for a given surface.
%
% Input Parameters:
%   Surf       : Surface.
%
% Output Parameters:
%   Surf      : Output Surfaces. The area values for each structure is stored 
%  in  StructS field inside Surf struct.
%
% Related references:
%
%
% See also: Red_Surf Smooth_Surf Plot_Surf Surf_Comp Atlas_Surf
% Plot_oversurf
%__________________________________________________
% Authors: Yasser Aleman Gomez
% Neuroimaging Department
% Cuban Neuroscience Center
% February 17th 2007
% Version $1.0

%=========================Main program====================================%
uni = sort(unique(Surf.Is));
Temp = zeros(size(Surf.SurfData.faces,1),1);
for k = 1:size(Surf.SurfData.faces,1)
    It = Surf.SurfData.faces(k,:)';
    if sum(Surf.Is(It)) == 0
        Temp(k,1) = 0;
    else
        c = accumarray(double(Surf.Is(It)),ones(length(It),1));
        indpos = find(c == max(c));
        Temp(k,1) =  indpos(1);
    end
end

fv.vertices = Surf.SurfData.vertices;
Nt = length(uni);
Surf.StructS = zeros(Nt,2);
for k =1:Nt
    ind = find(Temp == uni(k));

    if ~isempty(ind)
        fv.faces = Surf.SurfData.faces(ind,:);
        [At] = Area_Comp(fv);
        Surf.StructS(k,2) = At;
    else
        Surf.StructS(k,2) = 0;
    end
    Surf.StructS(k,1) = uni(k);
end
%========================End of main program==============================%
return

function [fv] = Ver_Ext_custom(V,NPoints,dat);
% (subfunction)
% Extract a mesh from a binary mask using matlab script isosurface.
%
% Author: Yasser Aleman Gomez
% Neuroimaging Department
% Cuban Neuroscience Center
% December 1 2006
% Version $1.0

%=========================Main program====================================%
dim = V.dim(1:3)';
siz = sqrt(sum(V.mat(1:3,1:3).^2));;
S = size(dat);
ind = find(dat);
[x,y,z] = ind2sub(S,ind);clear ind;
xc = single([min(x)-1:max(x)+1]');clear x;
yc = single([min(y)-1:max(y)+1]');clear y;
zc = single([min(z)-1:max(z)+1]');clear z;
[meshx,meshy,meshz]=meshgrid(xc,yc,zc);
dat=permute(dat,[2 1 3]);
tic;fv=isosurface(meshx,meshy,meshz,dat(min(yc):max(yc),min(xc):max(xc),min(zc):max(zc)),0,'verbose');toc;
clear meshx meshy meshz;
[fa,indtt] =sort(fv.faces');
[ford,indt2] =sortrows(fa');ford = double(ford);
temp = ford(1:end-1,:) - ford(2:end,:);temp = sum(abs(temp'));
fnd = find(temp == 0);
fv.faces(indt2(fnd),:)=[];
fv.faces = uint32(fv.faces);
fv.vertices = double(fv.vertices);
if (NPoints<=length(fv.vertices))&(NPoints ~=0)
    disp('   ');
    disp('Reducing Surface...');
    factor = 1/(length(fv.vertices)/NPoints);
    tic;fv=reducepatch(fv,factor);toc;
elseif (NPoints>length(fv.vertices))&(NPoints ~=0)
    disp('The number of points is greater than the number of vertex in the original surface. No reduction');
end
[fa,indtt] =sort(fv.faces');
[ford,indt2] =sortrows(fa');ford = double(ford);
temp = ford(1:end-1,:) - ford(2:end,:);temp = sum(abs(temp'));
fnd = find(temp == 0);
fv.faces(indt2(fnd),:)=[];
fv.faces = uint32(fv.faces);
fv.vertices = double(fv.vertices);
xt = -abs(V.mat(1,4))-2*siz(1);
yt = -abs(V.mat(2,4))-2*siz(2);
zt = -abs(V.mat(3,4))-2*siz(3);
% fv.vertices(:,1) =siz(1)*fv.vertices(:,1);
% fv.vertices(:,2) =siz(2)*fv.vertices(:,2);
% fv.vertices(:,3) =siz(3)*fv.vertices(:,3);
% fv.vertices = fv.vertices+repmat([xt yt zt],[size(fv.vertices,1),1]);
fv.normals  = patchnormals(fv);
norma = sqrt(sum((fv.normals').^2));
fv.normals = fv.normals./repmat(norma',[1 3]);

%========================End of main program==============================%
return;

function [Surf] = Surf_Extraction_custom(V,I);
%
% Syntax :
% [Surfa] = Surf_Comp(AImages);
%
% Computes the surface extraction for an atlas file or a mask file. If the
% imput images are atlas files, it extracts the surfaces for the structures
% specified in the structure list(StList).If the imput images are binary or
% mask images it will perform the surface extraction for all this images.
%
% Input Parameters:
%   AImages     : Individual Atlas files.
%
% Output Parameters:
%  Surfa        : Cell Array with surfaces in matlab variables format.
%
% Related references:
%
%
% See also: Red_Surf Smooth_Surf Plot_Surf
%__________________________________________________
% Authors: Yasser Aleman Gomez
% Neuroimaging Department
% Cuban Neuroscience Center
% December 1st 2006
% Version $1.0
warning off
fclose('all');
%=====================Checking Input Parameters===========================%
if nargin ==0
    [AImages,sts] = spm_select([1 Inf],'image','Selecting Atlased Images','',cd);
end

%=========================================================================%
%=========================Main program====================================%
%V = spm_vol(AImages);
Ns  = length(V);
for i = 1:Ns
    [pth nm ext] = fileparts(V(i).fname);
    %I = uint16(spm_read_vols(V(i))); I = imfill(I,'holes');
    I = uint16(I); I = imfill(I,'holes');
    dat = logical(I);
    voxsize = sqrt(sum(V(i).mat(1:3,1:3).^2));
    Surf = struct('Name','','Area','','SurfData','','Tri','','Orig','','Dim','','VoxSize','','Code','');
    bw = bwlabeln(double(dat));
    ind1 = find(dat);
    c = accumarray(bw(ind1),ones(length(bw(ind1)),1));
    indpos = find(c == max(c));
    dat(bw~= indpos) = 0;clear bw;
    It = logical(zeros(size(I,1)+4,size(I,2)+4,size(I,3)+4));It(3:size(I,1)+2,3:size(I,2)+2,3:size(I,3)+2) = dat;
    [fv] = Ver_Ext_custom(V(i),0,It);
    clear dat It;
    disp('  ');
    disp('Computing Neighbor Points .....');
    Surf.SurfData.faces = fv.faces;
    Surf.SurfData.vertices = fv.vertices;
    Surf.SurfData.VertexNormals = fv.normals;clear fv;
    Npoints = size(Surf.SurfData.vertices,1);
    Nfaces = size(Surf.SurfData.faces,1);
    tic;[Tri] = Vert_Neib(double(Surf.SurfData.faces),Npoints,Nfaces);toc;
    Temp = sum(Tri);
    Tri(:,Temp==0) = [];
    Surf.Tri = Tri;
    %try [Surf] = Surf_Ext_Corr(Surf); Surf = Surfa; end
    Surf.Orig = abs(V(i).mat(1:3,4)');
    Surf.Dim = abs(V(i).dim(1:3));
    Surf.VoxSize = voxsize;
    disp('   ');
     disp('Smoothing... ');
     try tic;[OutFiles,sSurf] = Smooth_surf(Surf,'',3,'n','');toc; end
     Surf= sSurf{1};
end
%========================End of main program==============================%
return;

%========================Internal Functions===============================%
function [V, P] = crop_st(P,VA, borders);
S = size(P);
borders = 2*borders;
ind = find(P);
[x,y,z] = ind2sub(S,ind);
P = P(min(x):max(x),min(y):max(y),min(z):max(z));
S = size(P);
P(end + borders(1),end + borders(2),end + borders(3)) = 0;
d = round((size(P)/2) - (S/2));
P = translateImageN0(P,d(1),d(2),d(3));
[pathstr,name,ext,vers] = fileparts(VA.fname);
V = VA;
V.dim(1:3) = size(P);
% vec = [min(x)-1 min(y)-1 min(z)-1] - borders/2; vec = V.mat(1:3,1:3)*vec';
% V.mat(1:3,4) = V.mat(1:3,4) + vec;
return;

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

function Surf = PaintSurf_with_FreeSurferColors(Surf);
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
