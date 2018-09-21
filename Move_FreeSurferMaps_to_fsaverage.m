function OutNames = Move_FreeSurferMaps_to_fsaverage(FSDir, Idfile);
%
% Syntax :
%     OutNames = Move_FreeSurferMaps_to_fsaverage(FSDir, Idfile);
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

Idfile = '/media/Data/PROCESSING_RESULTS/PCAMPO/Stats/Final_IDs.txt';
FSDir = '/media/Data/PROCESSING_RESULTS/PCAMPO/5-freesurfer_processing';
setenv('SUBJECTS_DIR',FSDir);


Ids = char(textread(Idfile,'%s'));
Ns = size(Ids,1);
OutNames = '';
for i = 1:Ns
    fsID = deblank(Ids(i,:));
     disp(['Processing Subject ' fsID '. Subject ' num2str(i) ' of ' num2str(Ns)]);
    % Removing processing directories
    files = dir([FSDir filesep fsID filesep 'surf' filesep 'tmp.mris_preproc.*' ]);
    [X{1:size(files,1)}] = deal(files.isdir); a = cell2mat(X); ind = find(a == 0);files(ind) = [];
    [Xn{1:size(files,1)}] = deal(files.name);
    dirnam = [repmat([FSDir filesep fsID filesep 'surf' filesep],[size(Xn,2) 1]) char(Xn)]; clear Xn X;
    for z = 1:size(dirnam,1)
        system(['rm -r ' deblank(dirnam(z,:))]);
    end
    
    % Removing Old files directories
    system(['rm -r ' FSDir filesep fsID filesep 'surf' filesep '*.reg.fsaverage.*' ]);
    % 1. Thickness
    cad = ['mris_preproc --s ' fsID ' --target fsaverage --hemi lh --meas thickness --out ' [FSDir filesep fsID filesep 'surf'  filesep 'lh.thickness.reg.fsaverage.mgh']];
    system(cad);
    % 2. Local Gyrification Index
    cad = ['mris_preproc --s ' fsID ' --target fsaverage --hemi lh --meas pial_lgi --out ' [FSDir filesep fsID filesep 'surf'  filesep 'lh.pial_lgi.reg.fsaverage.mgh']];
    system(cad);
    % 3. Area
    cad = ['mris_preproc --s ' fsID ' --target fsaverage --hemi lh --meas area --out ' [FSDir filesep fsID filesep 'surf'  filesep 'lh.area.reg.fsaverage.mgh']];
    system(cad);
    % 4. Volume
    cad = ['mris_preproc --s ' fsID ' --target fsaverage --hemi lh --meas volume --out ' [FSDir filesep fsID filesep 'surf'  filesep 'lh.volume.reg.fsaverage.mgh']];
    system(cad);
% % % %     % 5. Curvature
% % % %     cad = ['mris_preproc --s ' fsID ' --target fsaverage --hemi lh --meas curv --out ' [FSDir filesep fsID filesep 'surf'  filesep 'lh.curv.reg.fsaverage.mgh']];
% % % %     system(cad);
% % % %     % 6. Volume
% % % %     cad = ['mris_preproc --s ' fsID ' --target fsaverage --hemi lh --meas sulc --out ' [FSDir filesep fsID filesep 'surf'  filesep 'lh.sulc.reg.fsaverage.mgh']];
% % % %     system(cad);
    
    
    % Transforming some freesurfer measures
    % 1. Thickness
    cad = ['mris_preproc --s ' fsID ' --target fsaverage --hemi rh --meas thickness --out ' [FSDir filesep fsID filesep 'surf'  filesep 'rh.thickness.reg.fsaverage.mgh']];
    system(cad);
    % 2. Local Gyrification Index
    cad = ['mris_preproc --s ' fsID ' --target fsaverage --hemi rh --meas pial_lgi --out ' [FSDir filesep fsID filesep 'surf'  filesep 'rh.pial_lgi.reg.fsaverage.mgh']];
    system(cad);
    % 3. Area
    cad = ['mris_preproc --s ' fsID ' --target fsaverage --hemi rh --meas area --out ' [FSDir filesep fsID filesep 'surf'  filesep 'rh.area.reg.fsaverage.mgh']];
    system(cad);
    % 4. Volume
    cad = ['mris_preproc --s ' fsID ' --target fsaverage --hemi rh --meas volume --out ' [FSDir filesep fsID filesep 'surf'  filesep 'rh.volume.reg.fsaverage.mgh']];
    system(cad);
% % % %     % 5. Curvature
% % % %     cad = ['mris_preproc --s ' fsID ' --target fsaverage --hemi rh --meas curv --out ' [FSDir filesep fsID filesep 'surf'  filesep 'rh.curv.reg.fsaverage.mgh']];
% % % %     system(cad);
% % % %     % 6. Volume
% % % %     cad = ['mris_preproc --s ' fsID ' --target fsaverage --hemi rh --meas sulc --out ' [FSDir filesep fsID filesep 'surf'  filesep 'rh.sulc.reg.fsaverage.mgh']];
% % % %     system(cad);
    
    
    OutNames = strvcat(OutNames, [FSDir filesep fsID filesep 'surf'  filesep 'lh.thickness.reg.fsaverage.mgh'],...
                    [FSDir filesep fsID filesep 'surf'  filesep 'rh.thickness.reg.fsaverage.mgh'],...
                    [FSDir filesep fsID filesep 'surf'  filesep 'lh.pial_lgi.reg.fsaverage.mgh'],...
                    [FSDir filesep fsID filesep 'surf'  filesep 'rh.pial_lgi.reg.fsaverage.mgh'],...
                    [FSDir filesep fsID filesep 'surf'  filesep 'lh.area.reg.fsaverage.mgh'],...
                    [FSDir filesep fsID filesep 'surf'  filesep 'rh.area.reg.fsaverage.mgh'],...
                    [FSDir filesep fsID filesep 'surf'  filesep 'lh.volume.reg.fsaverage.mgh'],...
                    [FSDir filesep fsID filesep 'surf'  filesep 'rh.volume.reg.fsaverage.mgh'],...
                    [FSDir filesep fsID filesep 'surf'  filesep 'lh.curv.reg.fsaverage.mgh'],...
                    [FSDir filesep fsID filesep 'surf'  filesep 'rh.curv.reg.fsaverage.mgh'],...
                    [FSDir filesep fsID filesep 'surf'  filesep 'lh.sulc.reg.fsaverage.mgh'],...
                    [FSDir filesep fsID filesep 'surf'  filesep 'rh.sulc.reg.fsaverage.mgh']);
    
end
varargout{1} = OutNames;
return





