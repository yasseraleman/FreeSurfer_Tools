function Join_CH_Lobar_Results;

%OutFiles = Select_files('/media/Data/PROCESSING_RESULTS/REDES_LONGITUDINAL/5-freesurfer_processing',{'*Results_CHull_Lobar.txt'},{'Results_CHull_Lobar.txt'},'/media/COSAS/Test/Joost/LOBE_HULLPARC/Ids.txt');
names = textread('/media/Data/PROCESSING_RESULTS/REDES_LONGITUDINAL/5-freesurfer_processing/Ids.txt','%s');
names = char(names);
%OutFiles = [repmat('/media/COSAS/Test/Joost/LOBE_HULLPARC/',[size(char(names),1) 1]) char(names)    repmat('Results_CHull_Lobar.txt',[size(char(names),1) 1])];
OutFiles = '';
for i = 1:size(names,1)
    nname = ['/media/Data/PROCESSING_RESULTS/REDES_LONGITUDINAL/5-freesurfer_processing/temp/' deblank(names(i,:)) '-Results_CHull.txt' ];
    OutFiles = strvcat(OutFiles,nname);
end
Ns = size(OutFiles,1);
Cads = '';
for i = 1:Ns
    rfile = deblank(OutFiles(i,:));
    fid = fopen(rfile,'rt');
    if i == 1
        line = fgetl(fid);
        line1 = fgetl(fid);
        Cads =  strvcat(Cads,line,line1);
    else
                line = fgetl(fid);
        line1 = fgetl(fid);
        Cads =  strvcat(Cads,line1);
    end
end
Outfile = '/media/Data/PROCESSING_RESULTS/REDES_LONGITUDINAL/5-freesurfer_processing/temp/Results_CHull.txt';
fid = fopen(Outfile,'wt');
for j = 1:size(Cads,1)
    fprintf(fid,'%s\n',Cads(j,:));
end
fclose(fid);
return