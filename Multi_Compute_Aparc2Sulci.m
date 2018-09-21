function Multi_Compute_Aparc2Sulci
FreeSDir = '/media/Data/PROCESSING_RESULTS/HCP/5-freesurfer_processing';
BrainVisaDatabaseDir = '/media/MyDisk/PROCESSING_RESULTS/8-BrainVisaDataBase';
%Idfile = '/media/MyDisk/PROCESSING_RESULTS/5-freesurfer_processing/HCPrest_Ids.txt';
Idfile = '/media/Data/PROCESSING_RESULTS/HCP/5-freesurfer_processing/HCPrest_Ids.txt';
if exist(Idfile,'file')
    Ids = char(textread(Idfile,'%s'));
else
    Ids = Idfile;
end
%Ids = 'fsaverage';

%Ids =strvcat('HCP_111716-20140807-T1wMPR1','HCP_115320-20140807-T1wMPR1','HCP_118528-20140807-T1wMPR1','HCP_127630-20140807-T1wMPR1','HCP_129028-20140807-T1wMPR1','HCP_131217-20140807-T1wMPR1','HCP_133928-20140807-T1wMPR1','HCP_140925-20140807-T1wMPR1','HCP_146432-20140807-T1wMPR1');
Ns = size(Ids,1);
failed = '';
cad = '';
Ido = '';
for  i = 1:Ns
    
    Id = deblank(Ids(i,:));
    disp(['Processing Subject ' Id '. Subject ' num2str(i) ' of ' num2str(Ns)]);
    try
        delete([ FreeSDir filesep Id filesep  'surf' filesep Id 'SULCLINES_FAILED.txt']);
    end
    
    % Aparc
    lhinfl = [FreeSDir filesep Id filesep 'surf' filesep 'lh.inflated'];
    lhannot = [FreeSDir filesep Id filesep 'label' filesep 'lh.aparc.annot'];
    
    rhinfl = [FreeSDir filesep Id filesep 'surf' filesep 'rh.inflated'];
    rhannot = [FreeSDir filesep Id filesep 'label' filesep 'rh.aparc.annot'];
    
    OutAnnotFile = [FreeSDir filesep Id filesep 'label' filesep 'lh.aparc2sulci.annot'];
    %OutAnnotFile = Aparc2Sulci(lhinfl,lhannot, OutAnnotFile);
    
    OutAnnotFile = [FreeSDir filesep Id filesep 'label' filesep 'rh.aparc2sulci.annot'];
    %OutAnnotFile = Aparc2Sulci(rhinfl,rhannot, OutAnnotFile);
    
    
    % Aparc2009s to sulci
    lhannot = [FreeSDir filesep Id filesep 'label' filesep 'lh.aparc.a2009s.annot'];
    OutAnnotFile = [FreeSDir filesep Id filesep 'label' filesep 'lh.aparca2009s2sulci.annot'];
    OutAnnotFile = Aparc2009s2Sulci(lhinfl,lhannot, OutAnnotFile);
    
    rhannot = [FreeSDir filesep Id filesep 'label' filesep 'rh.aparc.a2009s.annot'];
    OutAnnotFile = [FreeSDir filesep Id filesep 'label' filesep 'rh.aparca2009s2sulci.annot'];
    OutAnnotFile = Aparc2009s2Sulci(rhinfl,rhannot, OutAnnotFile);
    
    
end
return;