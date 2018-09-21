function OutFile = Compute_Sulcal_Area(IdsFile,thresh, OutFile);
%
% Syntax :
% OutFile = Compute_Sulcal_Area(IdsFile,thr, OutFile);
%
% This script computes sulcal and medial wall surface area for a specified 
% number of subjects.
%
% Input Parameters:
%    IdsFile          : FreeSurfer IDs File
%     thr             : Curvature Threshold
%   OutFile           : Output File
%
% Output Parameters:
%     OutFile         : OutputFile
%
% See also:
%__________________________________________________
% Authors:  Yasser Aleman Gomez
% LIM
% February 23th 2012
% Version $1.0

Ids = char(textread(IdsFile,'%s'));
% Ids = 'fsaverage';
% thresh = -0.07;
if nargin < 3
   [pth,nm,ext] = fileparts(IdsFile);
   OutFile = [pth nm '_Sulcal_WM_SurfaceArea.txt'];
end
if nargin < 2
    thr = 0;
end
Ns = size(Ids,1);
[a,temp] = system('echo $SUBJECTS_DIR');
temp = temp';FreeSDir = deblank(temp(:)');
disp(['Subject ID    LH-Sucal-Surface    LH-Mwall-Surface    RH-Sucal-Surface    RH-Mwall-Surface ']);
for k = 1:Ns
    subjId = deblank(Ids(k,:));
    
    
    lhannot = [FreeSDir filesep subjId filesep 'label' filesep 'lh.aparc.annot'];
    rhannot = [FreeSDir filesep subjId filesep 'label' filesep 'rh.aparc.annot'];
    
    % Surfaces
    lhsurf = [FreeSDir filesep subjId filesep 'surf' filesep 'lh.pial'];
    rhsurf = [FreeSDir filesep subjId filesep 'surf' filesep 'rh.pial'];
    
    % Curvature Maps
    lhcurv = [FreeSDir filesep subjId filesep 'surf' filesep 'lh.curv'];
    rhcurv = [FreeSDir filesep subjId filesep 'surf' filesep 'rh.curv'];
    
    %----- Reading Surfaces
    [OutFiles, SurfF] = Exp_Surf(lhsurf, '0', '','', 'imp','n');
    Surfwl= SurfF{1};
    Surfwl.Name = 'LH.WHITE';
    
    [OutFiles, SurfF] = Exp_Surf(rhsurf, '0', '','', 'imp','n');
    Surfwr= SurfF{1};
    Surfwr.Name = 'RH.White';
    
    [txtl,ctabl] = read_cfiles(lhannot);
    tempname = char(ctabl.struct_names);
    indu = find(ismember(tempname(:,1:7),'unknown','rows') == 1);
    indcc = find(ismember(tempname(:,1:14),'corpuscallosum','rows') == 1);
    lrids = ctabl.table([indu indcc],5);
    [txtr,ctabr] = read_cfiles(rhannot);
    tempname = char(ctabr.struct_names);
    indu = find(ismember(tempname(:,1:7),'unknown','rows') == 1);
    indcc = find(ismember(tempname(:,1:14),'corpuscallosum','rows') == 1);
    rrids = ctabr.table([indu indcc],5);
    
    % ======================== Left Hemisphere ===========================%
    [left_curv,ctab] = read_cfiles(lhcurv);
    ind = find(left_curv>=thresh);
    ttemp = left_curv*0;
    ttemp(ind) = 1;
    ind = find(ismember(txtl,lrids));
    ttemp(ind) = 2;
    ind = find(txtl==0);
    ttemp(ind) = 2;
    ind = find(ttemp == 0);
    ttemp(ind) = 3;
    [colortable,labels] = Create_FS_Colortable(ttemp,strvcat('sulci','medial_wall','gyri'));
    OutAnnot = [FreeSDir filesep subjId filesep 'label' filesep 'lh.curvthr.annot'];
    OutAnnot = save_annotfiles(labels,OutAnnot,colortable);
    %% Computing Area
    fv.vertices = Surfwl.SurfData.vertices;
    for j = 1:2
        ind = find(ttemp == j);
        
        %----- Computing Area
        indf = find(sum(ismember(Surfwl.SurfData.faces,ind)')==3);
        At = 0;
        fv.faces = Surfwl.SurfData.faces(indf,:);
        N = size(fv.faces,1);
        for i = 1:N;
            di = dist(fv.vertices(fv.faces(i,:),:)');
            A = abs(di(1,2));
            B = abs(di(1,3));
            C = abs(di(2,3));
            p = (A+B+C)/2;
            Ar = sqrt(p*(p-A)*(p-B)*(p-C));
            %Ar = (A*B/2)*(sqrt(1-((A^2+B^2-C^2)^2)/((2*A*B)^2)));
            At = At+Ar;
        end
        larea(k,j) = real(At);% cm^2
    end
    %% End of Computing Area
    
    % ======================== Right Hemisphere ==========================%
    [right_curv,ctab] = read_cfiles(rhcurv);
    ind = find(right_curv>=thresh);
    ttemp = right_curv*0;
    ttemp(ind) = 1;
    ind = find(ismember(txtr,rrids));
    ttemp(ind) = 2;
    ind = find(txtr==0);
    ttemp(ind) = 2;
    ind = find(ttemp == 0);
    ttemp(ind) = 3;
    [colortable,labels] = Create_FS_Colortable(ttemp,strvcat('sulci','medial_wall','gyri'));
    OutAnnot = [FreeSDir filesep subjId filesep 'label' filesep 'rh.curvthr.annot'];
    OutAnnot = save_annotfiles(labels,OutAnnot,colortable);
    %% Computing Area
    fv.vertices = Surfwr.SurfData.vertices;
    for j = 1:2
        ind = find(ttemp == j);
        
        %----- Computing Area
        indf = find(sum(ismember(Surfwr.SurfData.faces,ind)')==3);
        At = 0;
        fv.faces = Surfwr.SurfData.faces(indf,:);
        N = size(fv.faces,1);
        for i = 1:N;
            di = dist(fv.vertices(fv.faces(i,:),:)');
            A = abs(di(1,2));
            B = abs(di(1,3));
            C = abs(di(2,3));
            p = (A+B+C)/2;
            Ar = sqrt(p*(p-A)*(p-B)*(p-C));
            %Ar = (A*B/2)*(sqrt(1-((A^2+B^2-C^2)^2)/((2*A*B)^2)));
            At = At+Ar;
        end
        rarea(k,j) = real(At);% cm^2
    end
    %% End of Computing Area
    disp([subjId '     '  num2str(larea(k,1))  '          ' num2str(larea(k,2)) '           ' num2str(rarea(k,1)) '          '  num2str(rarea(k,2))]);
end
% Saving results
cads = strvcat('Subject ID    LH-Sucal-Surface    LH-Mwall-Surface    RH-Sucal-Surface    RH-Mwall-Surface ',...
    [Ids repmat('     ',[Ns 1]) num2str(larea(:,1)) repmat('          ',[Ns 1]) num2str(larea(:,2)) repmat('           ',[Ns 1]) num2str(rarea(:,1)) repmat('          ',[Ns 1]) num2str(rarea(:,2))]);
%disp(cads);
fid = fopen(OutFile,'wt');
for i = 1:size(cads,1);
    fprintf(fid,'%s\n',cads(i,:))
end
fclose(fid);
return;