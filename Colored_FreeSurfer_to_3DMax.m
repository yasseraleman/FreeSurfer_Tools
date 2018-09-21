function OBJFile = Colored_FreeSurfer_to_3DMax(Surf, OBJFile);
%
% Syntax :
%OBJFile = Colored_FreeSurfer_to_3DMax(Surf, OBJFile);
%
% This function exports Matlab Surfaces into WaveFront OBJ format.
% Materials file will be also created.
%
% Input Parameters:
%   Surf         : Surface Variable (Matlab Format)
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
[Surf] = Reorg_Surf(Surf);
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
Nfaces = size(Surf.SurfData.faces,1);
fprintf(fidmlt,'%s\n',['# Multi/Sub Material__1 (' num2str(Nfaces) ') to come ']); 
fprintf(fidmlt,'%s\n','#');
fprintf(fid, '\n');
fprintf(fid, '%s\n','g');
fprintf(fid, '%s\n',['# object ' nm ' to come ...']);
fprintf(fid, '%s\n','#');
Mat = Surf.SurfData.vertices';
fprintf(fid,'v %.6f %.6f %.6f\n', Mat(:));
fprintf(fid, '%s\n',['# ' num2str(size(Surf.SurfData.vertices,1)) ' vertices']);
fprintf(fid, '\n');
fprintf(fid, '%s\n',['g ' nm]);
d = ones(size(Surf.SurfData.vertices,1),1) % Transparency vector
for i = 1:Nfaces
        Matname = sprintf('Face_%.8d',i);
        Color =   mean(Surf.SurfData.FaceVertexCData(Surf.SurfData.faces(i,:),:));      % Color
        Trs =     mean(d(Surf.SurfData.faces(i,:),:));              % Transparency
       disp(['Adding ' num2str(i) ' of ' num2str(Nfaces)]);

    %% Saving Structure Faces
    
    Mat = Surf.SurfData.faces(i,:)';
    fprintf(fid, '%s\n',['usemtl ' Matname]);
    fprintf(fid,'%s\n', ['s 2']);
    fprintf(fid,'f %u %u %u\n', Mat(:));
    fprintf(fid, '%s\n',['# 1 face']);
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
fprintf(fidmlt,'%s\n','# Multi/Sub Material__1 done');
fprintf(fidmlt,'%s\n','#');
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