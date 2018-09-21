function OutNames = FreeSurferMaps_for_Statistics(FSDir, Idfile);
%
% Syntax :
%     OutNames = FreeSurferMaps_for_Statistics(FSDir, Idfile);
%
% This function moves freesurfer surface maps to fsaverage space.
%
% Input Parameters:
%       FSDir                    : FreeSurfer Directory
%       Idfile                   : Ids file
%
% Output Parameters:
%       OutNames                 : Surface Maps in fsaverage space
%
% See also: 
%__________________________________________________
% Authors: Yasser Aleman Gomez
% LIM, HUGGM
% November 13th 2014
% Version $1.0

Idfile = '/media/Data/PROCESSING_RESULTS/NEFRO/5-freesurfer_processing/Ids.txt';
FSDir = '/media/Data/PROCESSING_RESULTS/NEFRO/5-freesurfer_processing';

[Ids, NhistClin, TiempoDialisis, VolDiurRes, GanancPeso, DiurRes, RatioUF, Nadir90, Nadir10] = textread('/media/Data/PROCESSING_RESULTS/NEFRO/Stats/Datos_Clínicos.txt','%s%s%u%u%f%u%f%u%u','delimiter',';','headerlines',1);
Ids = char(Ids);


%Ids = char(textread(Idfile,'%s'));
Ns = size(Ids,1);
OutNames = '';
intVariable ='area';
OutDir = '/media/Data/PROCESSING_RESULTS/NEFRO/Stats';
mkdir(OutDir)
contth = 0;
contlgi = 0;
contarea = 0;
contvol = 0;
contsulc = 0;
contcurv = 0;

for i = 1:Ns
    fsID = deblank(Ids(i,:));
    switch intVariable
        case 'thickness'
            if exist([FSDir filesep fsID filesep 'surf'  filesep 'lh.thickness.reg.fsaverage.mgh'],'file')&exist([FSDir filesep fsID filesep 'surf'  filesep 'lh.thickness.reg.fsaverage.mgh'],'file')
                contth = contth + 1;
                Thickness.LH(:,contth) = load_mgh([FSDir filesep fsID filesep 'surf'  filesep 'lh.thickness.reg.fsaverage.mgh']);
                Thickness.RH(:,contth) = load_mgh([FSDir filesep fsID filesep 'surf'  filesep 'rh.thickness.reg.fsaverage.mgh']);
                Thickness.Ids{contth,1} = fsID;
                Thickness.ord(contth) = i;
            end
            
        case 'pial_lgi'
            if exist([FSDir filesep fsID filesep 'surf'  filesep 'lh.pial_lgi.reg.fsaverage.mgh'],'file')&exist([FSDir filesep fsID filesep 'surf'  filesep 'lh.pial_lgi.reg.fsaverage.mgh'],'file')
                contlgi = contlgi + 1;
                Pial_LGI.LH(:,contlgi) = load_mgh([FSDir filesep fsID filesep 'surf'  filesep 'lh.pial_lgi.reg.fsaverage.mgh']);
                Pial_LGI.RH(:,contlgi) = load_mgh([FSDir filesep fsID filesep 'surf'  filesep 'rh.pial_lgi.reg.fsaverage.mgh']);
                Pial_LGI.Ids{contlgi,1} = fsID;
                Pial_LGI.ord(contlgi) = i;
            end
        case 'area'
            if exist([FSDir filesep fsID filesep 'surf'  filesep 'lh.area.reg.fsaverage.mgh'],'file')&exist([FSDir filesep fsID filesep 'surf'  filesep 'lh.area.reg.fsaverage.mgh'],'file')
                contarea = contarea + 1;
                Area.LH(:,contarea) = load_mgh([FSDir filesep fsID filesep 'surf'  filesep 'lh.area.reg.fsaverage.mgh']);
                Area.RH(:,contarea) = load_mgh([FSDir filesep fsID filesep 'surf'  filesep 'rh.area.reg.fsaverage.mgh']);
                Area.Ids{contarea,1} = fsID;
                Area.ord(contarea) = i;
            end
        case 'volume'
            if exist([FSDir filesep fsID filesep 'surf'  filesep 'lh.volume.reg.fsaverage.mgh'],'file')&exist([FSDir filesep fsID filesep 'surf'  filesep 'lh.volume.reg.fsaverage.mgh'],'file')
                contvol = contvol + 1;
                Volume.LH(:,contvol) = load_mgh([FSDir filesep fsID filesep 'surf'  filesep 'lh.volume.reg.fsaverage.mgh']);
                Volume.RH(:,contvol) = load_mgh([FSDir filesep fsID filesep 'surf'  filesep 'rh.volume.reg.fsaverage.mgh']);
                Volume.Ids{contvol,1} = fsID;
                Volume.ord(contvol) = i;
            end
        case 'sulc'
            if exist([FSDir filesep fsID filesep 'surf'  filesep 'lh.sulc.reg.fsaverage.mgh'],'file')&exist([FSDir filesep fsID filesep 'surf'  filesep 'lh.sulc.reg.fsaverage.mgh'],'file')
                contsulc = contsulc + 1;
                Sulc.LH(:,contsulc) = load_mgh([FSDir filesep fsID filesep 'surf'  filesep 'lh.sulc.reg.fsaverage.mgh']);
                Sulc.RH(:,contsulc) = load_mgh([FSDir filesep fsID filesep 'surf'  filesep 'rh.sulc.reg.fsaverage.mgh']);
                Sulc.Ids{contsulc,1} = fsID;
            end
        case 'curv'
            if exist([FSDir filesep fsID filesep 'surf'  filesep 'lh.curv.reg.fsaverage.mgh'],'file')&exist([FSDir filesep fsID filesep 'surf'  filesep 'lh.curv.reg.fsaverage.mgh'],'file')
                contcurv = contcurv + 1;
                Curv.LH(:,i) = load_mgh([FSDir filesep fsID filesep 'surf'  filesep 'lh.curv.reg.fsaverage.mgh']);
                Curv.RH(:,i) = load_mgh([FSDir filesep fsID filesep 'surf'  filesep 'rh.curv.reg.fsaverage.mgh']);
                Curv.Ids{contcurv,1} = fsID;
            end
    end
end

switch intVariable
    case 'thickness'
        save([OutDir filesep 'Thickness.mat'], 'Thickness');
        
    case 'pial_lgi'
        save([OutDir filesep 'Pial_LGI.mat'], 'Pial_LGI');
        
    case 'area'
        save([OutDir filesep 'Area.mat'], 'Area');
        
    case 'volume'
        save([OutDir filesep 'Volume.mat'], 'Volume');
        
    case 'sulc'
        save([OutDir filesep 'Sulc.mat'], 'Sulc');
        
    case 'curv'
        save([OutDir filesep 'Curv.mat'], 'Curv');
end
return
