function annot2IndivSpace(IdFile)

IdfIle = '/media/Data/PROCESSING_RESULTS/HCP/5-freesurfer_processing/HCPrest_Ids.txt';
srcsubject = 'fsaverage';
annotId = 'sulcimaxprob';

setenv('SUBJECTS_DIR',SubjectsDIr);
Ids = char(textread(IdfIle,'%s'));
Ns = size(Ids,1);
% for  i = 1:Ns 
%     trgsubject = deblank(Ids(i,:));
%     
%     cad = ['mri_surf2surf --srcsubject ' srcsubject ' --trgsubject ' trgsubject ' --hemi lh --sval-annot ' SubjectsDIr filesep srcsubject filesep 'label' filesep 'lh.' annotId '.annot'...
%         ' --tval ' SubjectsDIr filesep  trgsubject filesep 'label' filesep 'lh.' annotId '.annot'];
%     system(cad);
%     
%     
%    cad = ['mri_surf2surf --srcsubject ' srcsubject ' --trgsubject ' trgsubject ' --hemi rh --sval-annot ' SubjectsDIr filesep srcsubject filesep 'label' filesep 'rh.' annotId '.annot'...
%         ' --tval ' SubjectsDIr filesep  trgsubject filesep 'label' filesep 'rh.' annotId '.annot'];
%     system(cad);
%     
%     
% end
% return;




IndSubjecsFreesDir = '/media/Data/PROCESSING_RESULTS/HCP/5-freesurfer_processing';

Slhfiles = '/media/Data/PROCESSING_RESULTS/HCP/5-freesurfer_processing/fsaverage/label/lh.sulcimaxprob.annot';
Srhfiles = '/media/Data/PROCESSING_RESULTS/HCP/5-freesurfer_processing/fsaverage/label/rh.sulcimaxprob.annot';
annotId = 'sulcimaxprob';

%% ============ Creating variables to process the average subject ====== %%
   % ----------- Setting Enviroment for Freesurfer Average ------------- %
[a,temp] = system('echo $FREESURFER_HOME');
temp = temp';temp = deblank(temp(:)');
optst.pipe.freesdir = [temp filesep 'subjects'];
% 
% optst.pipe.freesdir = FSaverageDir;
setenv('SUBJECTS_DIR',optst.pipe.freesdir);
% --------------------- Other variables --------------------------------- %
optst.pipe.subjId = 'fsaverage'; % Subject ID

%% ============ End of Setting Enviroment for Freesurfer Average ======= %%


%% ======== Creating variables to process the Individual subjects ====== %%
   % ----------- Setting Enviroment for Freesurfer Directory ------------ %
opts.pipe.freesdir = IndSubjecsFreesDir;
setenv('SUBJECTS_DIR',opts.pipe.freesdir);

%% ============ End of Setting Enviroment for Freesurfer Average ======= %%


%% ================== Subjects Processing ============================== %%
%Ids = char(textread(IdsFile,'%s'));
VolAtlasfiles = '';
FailedFiles = '';
tempd = '';
Nscales = size(Slhfiles,1);
for i = 1:Nscales
    Sfile = deblank(Slhfiles(i,:));
    [pth,nm,ext] = fileparts(Sfile);
    ind  = strfind(nm,'.');
    nmt = nm(ind+1:end);
    %% ================= Creating Color Codes File ===================== %%
    ColorFile = [pth filesep nmt '_ColorFile.txt'];
    [txt,ctab] = read_cfiles(Sfile);
    r = num2str(ctab.table(:,1));
    g = num2str(ctab.table(:,2));
    b = num2str(ctab.table(:,3));
    names = '';
    for j = 1:size(ctab.struct_names,1)
        namet = char(ctab.struct_names(j));
        inddr = isspace(namet);
        namet(inddr) = [];
        names = strvcat(names,namet);
    end
    names = char(ctab.struct_names);
    ids = num2str([0:size(ctab.table,1)-1]');
    Total = [ids repmat('   ',[size(ids,1) 1]) names repmat('             ',[size(ids,1) 1]) r repmat(' ',[size(ids,1) 1]) g repmat(' ',[size(ids,1) 1]) b repmat(' ',[size(ids,1) 1]) num2str(zeros(size(b,1),1))];
    fid = fopen(ColorFile,'wt');
    for j = 1:size(Total,1)
        fprintf(fid,'%s\n',Total(j,:));
    end
    fclose(fid);
    
    %% ============== End of Creating Color Codes File ================= %%
    
    %% ==== Creating Multiscale Atlas for template subject ============= %%
    ind  = strfind(nm,'.');
    annotId = nm(ind(1)+1:end);
    setenv('SUBJECTS_DIR',optst.pipe.freesdir);
    cadl = ['mris_ca_train -n 2 -t ' ColorFile ' lh sphere.reg ' annotId ' fsaverage ' pth filesep 'lh.' annotId '.gcs'];
    system(cadl);
    cadr = ['mris_ca_train -n 2 -t ' ColorFile ' rh sphere.reg ' annotId ' fsaverage ' pth filesep 'rh.' annotId '.gcs'];
    system(cadr);
    setenv('SUBJECTS_DIR',opts.pipe.freesdir);
    %% ==== End of Creating Multiscale Atlas for template subject ====== %%
    Nsubj = size(Ids,1);
    for j= 1:Nsubj
        Id = deblank(Ids(j,:));
         disp(['Processing =======>  Scale: ' num2str(i) ' of ' num2str(Nscales) ' . =====> Subject ID: ' Id ' . ---  ' num2str(j) ' of ' num2str(Nsubj)]);
            %% =========== Moving Surfaces to individual space ============= %%
            cadl = ['mris_ca_label -orig white -novar -t ' ColorFile ' ' Id ' lh sphere.reg ' pth filesep 'lh.' annotId '.gcs '  opts.pipe.freesdir  filesep Id filesep 'label' filesep 'lh.' annotId '.annot'];
            system(cadl);
            cadr = ['mris_ca_label -orig white -novar -t ' ColorFile ' ' Id ' rh sphere.reg ' pth filesep 'rh.' annotId '.gcs '  opts.pipe.freesdir  filesep Id filesep 'label' filesep 'rh.' annotId '.annot'];
            system(cadr);
            %% =========== End of Moving Surfaces to individual space ====== %%
    end
end
return;


