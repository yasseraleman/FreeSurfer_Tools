function OutFile = save_free_surf(Surf,OutFile);
%
% Syntax :
% OutFile = save_free_surf(Surf,OutFile);
%
% Script file to save freesurfer surface files
%
% Input Parameters:
%   Surf              :  Surface variable
%   OutFile           :  Output surface filename
%
% Output Parameters:
%   OutFile           :  Output surface filename
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
% April 9th 2012
% Version $1.0

%=================== Main Program ========================================%
if size(Surf.SurfData.faces,2) == 3
    MNumber = 16777214;
elseif size(Surf.SurfData.faces,2) == 4
    MNumber  = 16777215;
    
end
fid = fopen(OutFile, 'wb', 'b') ;
fwrite(fid, bitand(bitshift(MNumber, -16), 255), 'uchar') ;
fwrite(fid, bitand(bitshift(MNumber, -8), 255), 'uchar') ;
fwrite(fid, bitand(MNumber, 255), 'uchar') ;
fwrite(fid, sprintf('created from IBASPM on %s\n',datestr(now)),'char');
fwrite(fid, sprintf('created from IBASPM on %s\n',datestr(now)),'char');
Npoints = size(Surf.SurfData.vertices,1) ;  % number of vertices
Nfaces = size(Surf.SurfData.faces,1) ;  % number of faces
fwrite(fid, Npoints,'int32');
fwrite(fid, Nfaces,'int32');
vertcol = reshape(Surf.SurfData.vertices',size(Surf.SurfData.vertices,1)*size(Surf.SurfData.vertices,2),1);
fwrite(fid, vertcol,'float32');
facecol = reshape(Surf.SurfData.faces',size(Surf.SurfData.faces,1)*size(Surf.SurfData.faces,2),1)-1;
fwrite(fid, facecol,'int32');
fclose(fid) ;
return;