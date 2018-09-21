function hf = Plot_FreeSurfer_networks_over_brain(CMatrix, stype, thresh);
%
% Syntax :
% hf = Plot_FreeSurfer_networks_over_brain(CMatrix, stype, thresh);
%
% This script plots freesurfer regions from a connectivity matrix.
%
%
% Input Parameters:
%      CMatrix                :Connectivity Matrix
%       stype                 : Freesurfer surface Type (pial,white,inflated)
%      thresh                 : Connectivity Threshold
%
% Output Parameters:
%        hf                   : Figure Handle
% 
% See also:
%__________________________________________________
% Authors: Yasser Aleman Gomez
% LIM, HUGGM
% March 22th 2013
% Version $1.0
close all;
thresh = 0;
stype = 'pial';
% % ----- No 2 Single connectivity matrix but with a predefined structure organization
% %  CMatrix = zeros(68,68);  % Uncomment this line to see this example
% %  CMatrix(8,42) = 1;      % Uncomment this line to see this example
% %  CMatrix(42,8) = 1;      % Uncomment this line to see this example
% %   CMatrix(9,43) = 1;      % Uncomment this line to see this example
% %   CMatrix(43,9) = 1;      % Uncomment this line to see this example
% %  CMatrix(26,54) = 1;      % Uncomment this line to see this example
% %  CMatrix(54,26) = 1;      % Uncomment this line to see this example
% %  CMatrix(34,68) = 1;      % Uncomment this line to see this example
% %  CMatrix(68,34) = 1;      % Uncomment this line to see this example

%% ====================== Freesurfer parcellation ========================%
% Subcortical Structures
SubCL = [10:13 17:18 26];
SubCR = [ 49:54 58];

% Selecting Frontal regions
FroIds = sort([1028 1003 1027 1018 1019 1020 1012 1014 1024 1017 1032 1026 1002]);
FroIdsR = FroIds + 1000;

% Selecting Temporal regions
TempIds = sort([1009 1015 1030 1001 1007 1034 1006 1033 1016]);
TempIdsR = TempIds + 1000;

% Selecting Parietal regions
ParIds = sort([1029 1008 1031 1022 1025 1010 1023]);
ParIdsR = ParIds + 1000;

% Selecting Occipital regions
OccIds = sort([1011 1013 1005 1021]);
OccIdsR = OccIds + 1000;

% Selecting Insula regions
InsIds = [1035];
InsIdsR = [2035];
OldOrg = sort([FroIds FroIdsR ParIds ParIdsR TempIds TempIdsR OccIds OccIdsR InsIds InsIdsR]);
%% ====================== End of Freesurfer parcellation =================%

%% ===================== Checking input parameters ====================== %

% =================== Detecting Connectivity type ====================== %
if nargin < 3
    thresh = 0;
end
if nargin < 2
    stype = 'inflated';
    thresh = 0;
end
if ischar(CMatrix)
    Connect = Reorder_Connect_Matrix_FreeSurfer(CMatrix);
    C = Connect.Matrix/max(Connect.Matrix(:));
    Names  = Connect.StructNames;
elseif isstruct(CMatrix)
    if isfield(CMatrix,'Matrix')
        C = CMatrix.Matrix;
        if isfield(CMatrix,'StructNames')
            Names  = CMatrix.StructNames;
        else
            if size(CMatrix,1) == 82
                [Names] = Load_Struct_Names;
            elseif size(CMatrix,1) == 68
                [Names] = Load_Struct_Names;
                Names = Names(15:end,:);
            else
                Names = num2str([1:size(C,1)]','%02d');
            end
        end
    else
        error('Connectivity Matrix does not exist');
        return;
    end
else
    C = CMatrix;
    if size(CMatrix,1) == 82
        [Names] = Load_Struct_Names;
        CMatrixtemp.Matrix = C;
        CMatrixtemp.StructNames = Names;
        CMatrixtemp.StructCodes = [SubCL SubCL OldOrg];
        Connect = Reorder_Connect_Matrix_FreeSurfer(CMatrixtemp);
        C = Connect.Matrix;
        Names  = Connect.StructNames;
    elseif size(CMatrix,1) == 68
        [Names] = Load_Struct_Names;
        Names = Names(15:end,:);
        CMatrixtemp.Matrix = C;
        CMatrixtemp.StructNames = Names;
        CMatrixtemp.StructCodes = OldOrg';
        Connect = Reorder_Connect_Matrix_FreeSurfer(CMatrixtemp);
        C = Connect.Matrix;
        Names  = Connect.StructNames;
    else
        if nargin < 3
            Names = num2str([1:size(C,1)]','%02d');
        end
    end
end
% =================== End of Detecting Connectivity type ================%
if sum(C(:))~=0
    LabNet = SubNetworks_Extraction(C,thresh);
    LabNet = triu(LabNet);
    a = accumarray(nonzeros(LabNet(:)),[ones(size(nonzeros(LabNet(:)),1),1)]);
    NNetworks = unique(nonzeros(LabNet(:)));
    NN = length(NNetworks);
    [b,it] = sort(a,'descend');
    Nsize_ord = NNetworks(it);
    col = [1 0 0;0 1 0; 0 0 1; 1 1 0;0 1 1;1 0 1;1 0.5 0;0 0.5 1;1 0 0.5;0.5 1 0; 0.25 0.5 1;0.6 0.3 0.52;0.7 0.5 0.9];
    Ncolor = size(col,1);
    re = floor(NN/Ncolor); col = repmat(col,[re+1 1]);
    %col(it,:) = col(1:length(it),:); % Selecting Colors
    %Colors = Col(LabNet(find(LabNet)),:);
    StructIds = [0 0 0 0];
    for i = 1:NN
        ind = find(LabNet == Nsize_ord(i));
        [X, Y] = ind2sub(size(LabNet), ind);
        Sts = unique([X;Y]);
        StructIds = [StructIds;Connect.StructCodes(Sts) repmat(col(i,:),[length(Sts) 1])];
        NetRep(i).Names = Names(Sts,:);
    end
    StructIds(1,:) = [];
    indl = find(StructIds(:,1) < 2000) ;
    indr = find(StructIds(:,1) >= 2000) ;
    colordef black;
    hf = figure('numbertitle','off','Color','black','Position',[0 0 1200 900]);
    
    if ~isempty(indl)&~isempty(indr)
        TempFilename = which('Plot_FreeSurfer_regions');
        [ptht, nmt, extt] = fileparts(TempFilename);
        switch lower(stype)
            case 'pial'
                if exist([ptht filesep 'AuxFiles' filesep 'fsaverage.lh.pial.mat'],'file');
                    load([ptht filesep 'AuxFiles' filesep 'fsaverage.lh.pial.mat']);
                else
                    load(which('fsaverage.lh.pial.mat'));
                end
                SurfL = Surf;
                if exist([ptht filesep 'AuxFiles' filesep 'fsaverage.rh.pial.mat'],'file');
                    load([ptht filesep 'AuxFiles' filesep 'fsaverage.rh.pial.mat']);
                else
                    load(which('fsaverage.rh.pial.mat'));
                end
                SurfR = Surf;
            case 'white'
                if exist([ptht filesep 'AuxFiles' filesep 'fsaverage.lh.white.mat'],'file');
                    load([ptht filesep 'AuxFiles' filesep 'fsaverage.lh.white.mat']);
                else
                    load(which('fsaverage.lh.pial.mat'));
                end
                SurfL = Surf;
                if exist([ptht filesep 'AuxFiles' filesep 'fsaverage.rh.white.mat'],'file');
                    load([ptht filesep 'AuxFiles' filesep 'fsaverage.rh.white.mat']);
                else
                    load(which('fsaverage.rh.pial.mat'));
                end
                SurfR = Surf;
            case 'inflated'
                if exist([ptht filesep 'AuxFiles' filesep 'fsaverage.lh.inflated.mat'],'file');
                    load([ptht filesep 'AuxFiles' filesep 'fsaverage.lh.inflated.mat']);
                else
                    load(which('fsaverage.lh.pial.mat'));
                end
                SurfL = Surf;
                if exist([ptht filesep 'AuxFiles' filesep 'fsaverage.rh.inflated.mat'],'file');
                    load([ptht filesep 'AuxFiles' filesep 'fsaverage.rh.inflated.mat']);
                else
                    load(which('fsaverage.rh.pial.mat'));
                end
                SurfR = Surf;
        end
        
        
        
%         switch lower(stype)
%             case 'pial'
%                 load('lh.pial.mat');
%                 SurfL = Surf;
%                 load('rh.pial.mat');
%                 SurfR = Surf;
%             case 'white'
%                 load('lh.white.mat');
%                 SurfL = Surf;
%                 load('rh.white.mat');
%                 SurfR = Surf;
%             case 'inflated'
%                 load('lh.inflated.mat');
%                 SurfL = Surf;
%                 load('rh.inflated.mat');
%                 SurfR = Surf;
%         end
        [sts,it] = sort(StructIds(:,1));
        Temp = StructIds(it,:);
        for i = 1:length(Temp(:,1))
            if sts(i) >=2000
                ind = find(SurfR.SCodes(:,6) == Temp(i,1)-1000);
                inds = find(ismember(SurfR.Is,SurfR.SCodes(ind,5)));
                SurfR.SurfData.FaceVertexCData(inds,:) = repmat(Temp(i,2:end),[length(inds) 1]);
            else
                ind = find(SurfL.SCodes(:,6) == Temp(i,1));
                inds = find(ismember(SurfL.Is,SurfL.SCodes(ind,5)));
                SurfL.SurfData.FaceVertexCData(inds,:) = repmat(Temp(i,2:end),[length(inds) 1]);
            end
        end
        subplot(2,2,1);
        h1 = custom_plotsurf(SurfL);
        view([270 0]);axis off;axis tight;axis equal; h=title(['Left Hemisphere. Lateral View']);
        camlight;
        subplot(2,2,3);
        h1 = custom_plotsurf(SurfL);
        view([90 0]);axis off;axis tight;axis equal; h=title(['Left Hemisphere. Medial View']);
        camlight;
        
        
        subplot(2,2,2);
        h1 = custom_plotsurf(SurfR);
        view([90 0]);axis off;axis tight;axis equal; h=title(['Right Hemisphere. Lateral View']);
        camlight;
        subplot(2,2,4);
        h1 = custom_plotsurf(SurfR);
        view([270 0]);axis off;axis tight;axis equal; h=title(['Right Hemisphere. Medial View']);
        camlight;  
        
        
        
    elseif ~isempty(indl)&isempty(indr)
        switch lower(stype)
            case 'pial'
                load('lh.pial.mat');
                SurfL = Surf;
            case 'white'
                load('lh.white.mat');
                SurfL = Surf;
            case 'inflated'
                load('lh.inflated.mat');
                SurfL = Surf;
        end
        [sts,it] = sort(StructIds(:,1));
        Temp = StructIds(it,:);
        for i = 1:length(sts)
                ind = find(SurfL.SCodes(:,6) == Temp(i,1));
                inds = find(ismember(SurfL.Is,SurfL.SCodes(ind,5)));
                SurfL.SurfData.FaceVertexCData(inds,:) = repmat(Temp(i,2:end),[length(inds) 1]);
        end
        subplot(1,2,1);
        h1 = custom_plotsurf(SurfL);
        view([270 0]);axis off;axis tight;axis equal; h=title(['Left Hemisphere. Lateral View']);
        camlight;
        subplot(1,2,2);
        h1 = custom_plotsurf(SurfL);
        view([90 0]);axis off;axis tight;axis equal; h=title(['Left Hemisphere. Medial View']);
        camlight;
    elseif isempty(indl)&~isempty(indr)
        switch lower(stype)
            case 'pial'
                load('rh.pial.mat');
                SurfR = Surf;
            case 'white'
                load('rh.white.mat');
                SurfR = Surf;
            case 'inflated'
                load('rh.inflated.mat');
                SurfR = Surf;
        end
        [sts,it] = sort(StructIds(:,1));
        Temp = StructIds(it,:);
        for i = 1:length(sts)
                ind = find(SurfR.SCodes(:,6) == Temp(i,1)-1000);
                inds = find(ismember(SurfR.Is,SurfR.SCodes(ind,5)));
                SurfR.SurfData.FaceVertexCData(inds,:) = repmat(Temp(i,2:end),[length(inds) 1]);
        end
        
        
        subplot(1,2,1);
        h1 = custom_plotsurf(SurfR);
        view([90 0]);axis off;axis tight;axis equal; h=title(['Right Hemisphere. Lateral View']);
        camlight;
        subplot(1,2,2);
        h1 = custom_plotsurf(SurfR);
        view([270 0]);axis off;axis tight;axis equal; h=title(['Right Hemisphere. Medial View']);
        camlight;  
    end
end
return


function Names = Load_Struct_Names;

Names = strvcat('Left-Thalamus-Proper',...
    'Left-Caudate',...
    'Left-Putamen',...
    'Left-Pallidum',...
    'Left-Hippocampus',...
    'Left-Amygdala',...
    'Left-Accumbens-area',...
    'Right-Thalamus-Proper',...
    'Right-Caudate',...
    'Right-Putamen',...
    'Right-Pallidum',...
    'Right-Hippocampus',...
    'Right-Amygdala',...
    'Right-Accumbens-area',...
    'ctx-lh-bankssts',...
    'ctx-lh-caudalanteriorcingulate',...
    'ctx-lh-caudalmiddlefrontal',...
    'ctx-lh-cuneus',...
    'ctx-lh-entorhinal',...
    'ctx-lh-fusiform',...
    'ctx-lh-inferiorparietal',...
    'ctx-lh-inferiortemporal',...
    'ctx-lh-isthmuscingulate',...
    'ctx-lh-lateraloccipital',...
    'ctx-lh-lateralorbitofrontal',...
    'ctx-lh-lingual',...
    'ctx-lh-medialorbitofrontal',...
    'ctx-lh-middletemporal',...
    'ctx-lh-parahippocampal',...
    'ctx-lh-paracentral',...
    'ctx-lh-parsopercularis',...
    'ctx-lh-parsorbitalis',...
    'ctx-lh-parstriangularis',...
    'ctx-lh-pericalcarine',...
    'ctx-lh-postcentral',...
    'ctx-lh-posteriorcingulate',...
    'ctx-lh-precentral',...
    'ctx-lh-precuneus',...
    'ctx-lh-rostralanteriorcingulate',...
    'ctx-lh-rostralmiddlefrontal',...
    'ctx-lh-superiorfrontal',...
    'ctx-lh-superiorparietal',...
    'ctx-lh-superiortemporal',...
    'ctx-lh-supramarginal',...
    'ctx-lh-frontalpole',...
    'ctx-lh-temporalpole',...
    'ctx-lh-transversetemporal',...
    'ctx-lh-insula',...
    'ctx-rh-bankssts',...
    'ctx-rh-caudalanteriorcingulate',...
    'ctx-rh-caudalmiddlefrontal',...
    'ctx-rh-cuneus',...
    'ctx-rh-entorhinal',...
    'ctx-rh-fusiform',...
    'ctx-rh-inferiorparietal',...
    'ctx-rh-inferiortemporal',...
    'ctx-rh-isthmuscingulate',...
    'ctx-rh-lateraloccipital',...
    'ctx-rh-lateralorbitofrontal',...
    'ctx-rh-lingual',...
    'ctx-rh-medialorbitofrontal',...
    'ctx-rh-middletemporal',...
    'ctx-rh-parahippocampal',...
    'ctx-rh-paracentral',...
    'ctx-rh-parsopercularis',...
    'ctx-rh-parsorbitalis',...
    'ctx-rh-parstriangularis',...
    'ctx-rh-pericalcarine',...
    'ctx-rh-postcentral',...
    'ctx-rh-posteriorcingulate',...
    'ctx-rh-precentral',...
    'ctx-rh-precuneus',...
    'ctx-rh-rostralanteriorcingulate',...
    'ctx-rh-rostralmiddlefrontal',...
    'ctx-rh-superiorfrontal',...
    'ctx-rh-superiorparietal',...
    'ctx-rh-superiortemporal',...
    'ctx-rh-supramarginal',...
    'ctx-rh-frontalpole',...
    'ctx-rh-temporalpole',...
    'ctx-rh-transversetemporal',...
    'ctx-rh-insula');
return

function strsurf=custom_plotsurf(Surf);
Surf.SurfData.FaceColor = 'interp';
if isunix
    strsurf=patch(Surf.SurfData,'edgecolor','none', 'tag','patch','facelighting','gouraud');
else
    strsurf=patch(Surf.SurfData,'edgecolor','black','tag', 'patch','facelighting','gouraud');
end
return