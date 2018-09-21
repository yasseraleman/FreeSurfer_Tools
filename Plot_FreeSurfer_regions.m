function [hl, h1] = Plot_FreeSurfer_regions(stype,hemi,ColorFile);
%
% Syntax :
% Plot_FreeSurfer_regions(stype,hemi,ColorFile);
%
% This script plots freesurfer regions.
%
%
% Input Parameters:
%       stype                 : Freesurfer surface Type (pial,white,inflated)
%       hemi                  : Hemisphere (lh or rh)
%     ColorFile               : Color file
%
%
% Output Parameters:
%
% See also:
%__________________________________________________
% Authors: Yasser Aleman Gomez
% LIM, HUGGM
% March 22th 2013
% Version $1.0



% stype = 'inflated'; % Surface Type (pial, white or inflated)
% hemi = 'lh'; %Hemisphere
% ColorFile = '/media/COSAS/Test/Erika/ColorFile.txt';
% % if nargin < 2
% %     stype = 'pial';
% % end
% % if nargin < 2
% %     hemi = 'lh';
% % end
% % if nargin < 3
% %     ColorFile = '';
% % end
    
if ~exist('ColorFile','var');
    ColorFile = '';
end
%% ==================== Selecting surface to load ========================%
TempFilename = which('Plot_FreeSurfer_regions');
[ptht, nmt, extt] = fileparts(TempFilename);
switch lower(stype)
    case 'pial'
        switch lower(hemi)
            case 'lh'
                if exist(['/media/COSAS/scripts/Share_Data' filesep 'fsaverage.lh.pial.mat'],'file');
                    load(['/media/COSAS/scripts/Share_Data' filesep 'fsaverage.lh.pial.mat']);
                else
                    load(which('fsaverage.lh.pial.mat'));
                end
            case 'rh'
                if exist(['/media/COSAS/scripts/Share_Data' filesep 'fsaverage.rh.pial.mat'],'file');
                    load(['/media/COSAS/scripts/Share_Data' filesep 'fsaverage.rh.pial.mat']);
                else
                    load(which('fsaverage.rh.pial.mat'));
                end
        end
        
    case 'white'
        switch lower(hemi)
            case 'lh'
                if exist(['/media/COSAS/scripts/Share_Data' filesep 'fsaverage.lh.white.mat'],'file');
                    load(['/media/COSAS/scripts/Share_Data' filesep 'fsaverage.lh.white.mat']);
                else
                    load(which('fsaverage.lh.pial.mat'));
                end
            case 'rh'
                if exist(['/media/COSAS/scripts/Share_Data' filesep 'fsaverage.rh.white.mat'],'file');
                    load(['/media/COSAS/scripts/Share_Data' filesep 'fsaverage.rh.white.mat']);
                else
                    load(which('fsaverage.rh.pial.mat'));
                end
        end
    case 'inflated'
        switch lower(hemi)
            case 'lh'
                if exist(['/media/COSAS/scripts/Share_Data' filesep 'fsaverage.lh.inflated.mat'],'file');
                    load(['/media/COSAS/scripts/Share_Data' filesep 'fsaverage.lh.inflated.mat']);
                else
                    load(which('fsaverage.lh.pial.mat'));
                end
            case 'rh'
                if exist(['/media/COSAS/scripts/Share_Data' filesep 'fsaverage.rh.inflated.mat'],'file');
                    load(['/media/COSAS/scripts/Share_Data' filesep 'fsaverage.rh.inflated.mat']);
                else
                    load(which('fsaverage.rh.pial.mat'));
                end
        end
end
ind = find(Surf.Is == 0);
Surf.SurfData.FaceVertexCData(ind,:) = repmat([1 1 1],[length(ind) 1]);
%% ================ End Selecting surface to load ======================= %

%% ==================== Lobar parcellation ============================== %
% Selecting Frontal regions
FroIds = sort([1028 1003 1027 1018 1019 1020 1012 1014 1024 1017 1032 1026 1002]);

% Selecting Temporal regions
TempIds = sort([1009 1015 1030 1001 1007 1034 1006 1033 1016]);
TempIdsR = TempIds + 1000;

% Selecting Parietal regions
ParIds = sort([1029 1008 1031 1022 1025 1010 1023]);

% Selecting Occipital regions
OccIds = sort([1011 1013 1005 1021]);

% Selecting Insula regions
InsIds = [1035];

%% ===================== End of Lobar parcellation =======================%
Names = strvcat(Surf.SNames,'frontallobe','parietallobe','temporallobe','occipitallobe','insula');
r = [Surf.SCodes(:,1);255;0;0;255;255];
g = [Surf.SCodes(:,2);0;255;0;255;192];
b = [Surf.SCodes(:,3);0;0;255;0;32];
if exist(ColorFile,'file')
    % Reading Colorfile
    fid = fopen(ColorFile);
    cont = 0;
    Namesi = '';
    while 1
        cont = cont + 1;
        line = fgetl(fid);
        if ~ischar(line),   break,   end
        Namet = strread(line,'%s','delimiter',' ');
        if length(Namet) < 4
            Namesi = strvcat(Namesi,char(Namet(1)));
        else
            Namesi = strvcat(Namesi,char(Namet(1)));
            nametemp = char(Namet(1));
            indl = find(ismember(Names(:,1:length(nametemp)),nametemp,'rows'));
            r(indl,1) = str2num(char(Namet(2)));
            g(indl,1) = str2num(char(Namet(3)));
            b(indl,1) = str2num(char(Namet(4)));
        end
    end
    fclose(fid);
    clear tline cont;
    Ns = size(Namesi,1);
elseif exist('ColorFile','var')
    Namesi = ColorFile;
    Ns = size(ColorFile,1);
else
    Ns = 39;
    Namesi = Names;
end
r = r/255;g = g/255;b = b/255;
switch hemi
    case 'lh'
        SurfL = Surf;
    case   'rh'
        SurfR = Surf;
    otherwise
        SurfL = Surf;
end
for i = 1:Ns
    Namet = deblank(Namesi(i,:));
    N = length(Namet);
    ind = find(ismember(Names(:,1:N),Namet,'rows'));ind = ind(1);
    if ind < 35
        inds = Surf.SCodes(ind,5);
        indc = find(Surf.Is == inds);
        if (r(ind) + g(ind) +b(ind) == 0)
            color = Surf.SCodes(ind,1:3);
        else
            color = [r(ind)   g(ind)  b(ind)];
        end
        SurfL.SurfData.FaceVertexCData(indc,:) = repmat(color,[length(indc) 1]);
    else
        if ~exist('SurfL','var')
            SurfL = Surf;
        end
        if ind == 35
            inds = find(ismember(Surf.SCodes(:,6),FroIds));
            sids = Surf.SCodes(inds,5);
            indc = find(ismember(Surf.Is,sids));
            if (r(ind) + g(ind) +b(ind) == 0)
                color = [1 0 0];
            else
                color = [r(ind)   g(ind)  b(ind)];
            end
            SurfL.SurfData.FaceVertexCData(indc,:) = repmat([1 0 0],[length(indc) 1]);
        elseif ind == 36
            inds = find(ismember(Surf.SCodes(:,6),ParIds));
            sids = Surf.SCodes(inds,5);
            indc = find(ismember(Surf.Is,sids));
            if (r(ind) + g(ind) +b(ind) == 0)
                color = [0 1 0];
            else
                color = [r(ind)   g(ind)  b(ind)];
            end
            SurfL.SurfData.FaceVertexCData(indc,:) = repmat([0 1 0],[length(indc) 1]);
        elseif ind == 37
                        inds = find(ismember(Surf.SCodes(:,6),TempIds));
            sids = Surf.SCodes(inds,5);
            indc = find(ismember(Surf.Is,sids));
            if (r(ind) + g(ind) +b(ind) == 0)
                color = [0 0 1];
            else
                color = [r(ind)   g(ind)  b(ind)];
            end
            SurfL.SurfData.FaceVertexCData(indc,:) = repmat([0 0 1],[length(indc) 1]);
        elseif ind == 38
            inds = find(ismember(Surf.SCodes(:,6),OccIds));
            sids = Surf.SCodes(inds,5);
            indc = find(ismember(Surf.Is,sids));
            if (r(ind) + g(ind) +b(ind) == 0)
                color = [1 1 0];
            else
                color = [r(ind)   g(ind)  b(ind)];
            end
            SurfL.SurfData.FaceVertexCData(indc,:) = repmat([1 1 0],[length(indc) 1]);
        elseif ind == 39
            inds = find(ismember(Surf.SCodes(:,6),InsIds));
            sids = Surf.SCodes(inds,5);
            indc = find(ismember(Surf.Is,sids));
            if (r(ind) + g(ind) +b(ind) == 0)
                color = [1 1 0];
            else
                color = [r(ind)   g(ind)  b(ind)];
            end
            SurfL.SurfData.FaceVertexCData(indc,:) = repmat([1 1 0],[length(indc) 1]);
        end
    end
end

     
if exist('SurfR','var')|exist('SurfL','var')
    colordef white;
    hf = figure('numbertitle','off','Color','white','Position',[0 0 1200 900]);
    if exist('SurfR','var')&exist('SurfL','var')
        subplot(1,2,1);
        h1 = custom_plotsurf(SurfR);
        if strfind(hemi,'lh')
            view([270 0]);axis off;axis tight;axis equal; %h=title(['Left Hemisphere. Lateral View']);set(h,'FontSize',15,'FontName','Arial');
            hl = camlight;
        elseif strfind(hemi,'rh')
            view([90 0]);axis off;axis tight;axis equal; %h=title(['Right Hemisphere. Lateral View']);set(h,'FontSize',15,'FontName','Arial');
            hl = camlight;
        end
        
        subplot(1,2,2);
        h1 = custom_plotsurf(SurfL);
        if strfind(hemi,'lh')
            view([270 0]);axis off;axis tight;axis equal; %h=title(['Left Hemisphere. Lateral View']);set(h,'FontSize',15,'FontName','Arial');
            hl = camlight;
        elseif strfind(hemi,'rh')
            view([90 0]);axis off;axis tight;axis equal; %h=title(['Right Hemisphere. Lateral View']);set(h,'FontSize',15,'FontName','Arial');
            hl = camlight;
        end
    elseif exist('SurfR','var')&~exist('SurfL','var')
        h1 = custom_plotsurf(SurfR);
        if strfind(hemi,'lh')
            view([270 0]);axis off;axis tight;axis equal; %h=title(['Left Hemisphere. Lateral View']);set(h,'FontSize',15,'FontName','Arial');
            hl = camlight;
        elseif strfind(hemi,'rh')
            view([90 0]);axis off;axis tight;axis equal; %h=title(['Right Hemisphere. Lateral View']);set(h,'FontSize',15,'FontName','Arial');
            hl = camlight;
        end
    elseif ~exist('SurfR','var')&exist('SurfL','var')
        h1 = custom_plotsurf(SurfL);
        if strfind(hemi,'lh')
            view([270 0]);axis off;axis tight;axis equal; %h=title(['Left Hemisphere. Lateral View']);set(h,'FontSize',15,'FontName','Arial');
            hl = camlight;
        elseif strfind(hemi,'rh')
            view([90 0]);axis off;axis tight;axis equal; %h=title(['Right Hemisphere. Lateral View']);set(h,'FontSize',15,'FontName','Arial');
            hl = camlight;
        end
    end
end
return

function strsurf=custom_plotsurf(Surf);
Surf.SurfData.FaceColor = 'interp';
if isunix
    strsurf=patch(Surf.SurfData,'edgecolor','none', 'tag','patch','facelighting','gouraud');
else
    strsurf=patch(Surf.SurfData,'edgecolor','black','tag', 'patch','facelighting','gouraud');
end
return