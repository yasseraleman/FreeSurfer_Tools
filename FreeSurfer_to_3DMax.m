function OBJFile = FreeSurfer_to_3DMax(Surf, colortable, OBJFile);
%
% Syntax :
% OBJFile = FreeSurfer_to_3DMax(Surf, txt, colortable, OBJFile);
%
% This function exports Matlab Surfaces into WaveFront OBJ format.
% Materials file will be also created.
%
% Input Parameters:
%   Surf         : Surface Variable (Matlab Format)
%   colortable   : Struct variable similar to freesurfer's colortable. It
%                  includes 2 fields:
%                  1. struct_names: Cellarray containing 
%                  structures names.
%                  2. table. Nstructures X 6 Matrix containing some
%                  material parameters (Colors, Transparency ( Tr), Ids)
%                  Table order (R G B 0 Id Tr)
%    OBJFile     : Output Object Filename
%
% Output Parameters:
%  OBJFile      : Output Object Filename
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

fid = fopen(OBJFile,'wt');
fprintf(fid,'%s\n','# Max2Obj Version 4.0 Mar 10th, 2001');
fprintf(fid,'%s\n','#');
[pth,nm,ext] = fileparts(OBJFile);
ind = strfind(nm,'.');nm(ind) = [];
filename = [deblank(pth) filesep deblank(nm) '.mtl' ];
fidmlt = fopen(filename,'wt');
fprintf(fidmlt,'%s\n','# Max2Mtl Version 4.0 Mar 10th, 2001');
fprintf(fidmlt,'%s\n','#');
fprintf(fid, '%s\n',['mtllib ./' deblank(nm) '.mtl']);
txt = Surf.Is;
sts = unique(txt);
d = colortable.table(:,6); % Transparency vector
Ns = length(sts);
Nfaces = 0;
for i = 1:Ns
    ind = find(txt==sts(i));
    indc = find(colortable.table(:,5)==sts(i));
    if isempty(indc)
        Matname = 'unknown_structure';
        Color =   [0 0 0];      % Color
        Trs =   1;              % Transparency
    else
        Matname = char(colortable.struct_names{indc});
        Color =  colortable.table(indc,1:3)/255;       % Color
        Trs =      d(indc)  ;   % Transparency
    end
    disp(['Adding ' Matname ' Surface']);
    Surft = Surf;
    Faces = Surf.SurfData.faces;
    a = ismember(Surf.SurfData.faces,ind);
    ind2 = find(sum(a')==0);
    Surft.SurfData.faces(ind2,:) = [];
    [Surft] = Reorg_Surf(Surft);
    
    %-----------------------------------
%     fv.vertices = double(Surft.SurfData.vertices);
%     fv.faces = Surft.SurfData.faces;
%     factor = 1/(size(fv.vertices,1)/200000);
%     fv=reducepatch(fv,factor);
%     Nv = size(fv.vertices,1);
%     Nf = size(fv.faces,1);
%     Surft.SurfData.vertices = fv.vertices;
%     Surft.SurfData.faces = fv.faces;
% %     [OutFiles,sSurf] = Smooth_surf(Surft,'',3,'n','');
% %     Surft= sSurf{1};
%     %Surft = freesCS2brainvisaCS(Surft,'/media/COSAS/Test/freesurfer/ch2/tmp/T1.nii','b2f');
    %----------------------------------
    
    
    %% Saving Vertices
    fprintf(fid, '\n');
    fprintf(fid, '%s\n','g');
    fprintf(fid, '%s\n',['# object ' Matname ' to come ...']);
    fprintf(fid, '%s\n','#');
    Mat = Surft.SurfData.vertices';
    fprintf(fid,'v %.6f %.6f %.6f\n', Mat(:));
    fprintf(fid, '%s\n',['# ' num2str(size(Surft.SurfData.vertices,1)) ' vertices']);
    fprintf(fid, '\n');
    fprintf(fid, '%s\n',['g ' Matname]);
    %% Saving Structure Faces
    
    Mat = Surft.SurfData.faces'+Nfaces;
    Newfaces = size(Surft.SurfData.vertices,1);
    Nfaces = Nfaces+max(Surft.SurfData.faces(:));
    fprintf(fid, '%s\n',['usemtl ' Matname]);
    fprintf(fid,'%s\n', ['s 2']);
    fprintf(fid,'f %u %u %u\n', Mat(:));
    fprintf(fid, '%s\n',['# ' num2str(size(Surft.SurfData.faces,1)) ' faces']);
    % -----------------------------------------
    %% Creating Material Structure
    fprintf(fidmlt,'%s\n',['newmtl ' Matname]);
    fprintf(fidmlt,'Ka  %.1f %.1f %.1f\n',Color(1),Color(2),Color(3));
    fprintf(fidmlt,'Kd  %.1f %.1f %.1f\n',Color(1),Color(2),Color(3));
    fprintf(fidmlt,'Ks  %.1f %.1f %.1f\n',1,1,1);
    fprintf(fidmlt,'%s\n',['d ' sprintf('%.1f',Trs)]);
    fprintf(fidmlt,'%s\n',['Ns 8.0']);
    fprintf(fidmlt,'%s\n',['illum 2']);
    fprintf(fidmlt,'%s\n',['#']);
    % ----------------------------------------
end
fprintf(fid, '\n');
fprintf(fid, '%s','g');
fprintf(fidmlt,'%s\n','# EOF');
fclose(fid);
fclose(fidmlt);
return

function [Surft] = Reorg_Surf(Surf);
%This scripts reorganize Surfaces in case of deleted faces. 

indf = unique(Surf.SurfData.faces(:));
%%%%%%%%%%%%%%%%%%%%%%%Corriegiendo Superficies para quedarme solo
%%%%%%%%%%%%%%%%%%%%%%%con la parte que me interesa
Surft = Surf;
Surft.SurfData.vertices = zeros(size(indf,1),3);
Surft.Is = zeros(size(indf,1),1);
Surft.SurfData.VertexNormals = zeros(size(indf,1),3);
Surft.SurfData.FaceVertexCData = zeros(size(indf,1),3);
Surft.SurfData.faces = 0*Surft.SurfData.faces;
for i =1:size(indf,1)
    Surft.SurfData.vertices(i,:) = Surf.SurfData.vertices(indf(i),:);
    if isfield(Surf.SurfData,'VertexNormals')
        Surft.SurfData.VertexNormals(i,:) = Surf.SurfData.VertexNormals(indf(i),:);
    end
    if isfield(Surf.SurfData,'FaceVertexCData')
        Surft.SurfData.FaceVertexCData(i,:) = Surf.SurfData.FaceVertexCData(indf(i),:);
    end
    if isfield(Surf.SurfData,'Is')
        Surft.Is(i,:) = Surf.Is(indf(i),:);
    end
    indn = find(Surf.SurfData.faces ==indf(i));Surft.SurfData.faces(indn) = i;
end