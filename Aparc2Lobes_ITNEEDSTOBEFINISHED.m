
opts.pipe.freesdir = '/media/COSAS/Test/freesurfer';
opts.pipe.subjId = 'ch2';
[txt,colortable] = read_cfiles([opts.pipe.freesdir filesep opts.pipe.subjId filesep 'label' filesep 'lh.aparc.annot']);
[OutFiles, SurfF] = Exp_Surf([opts.pipe.freesdir filesep opts.pipe.subjId filesep 'surf' filesep 'lh.pial'], '0', '','', 'imp','n');
Surftl= SurfF{1};




%% ============================= Defining Lobes ======================== %%
% Fron = strvcat('superiorfrontal ','caudalmiddlefrontal ','rostralmiddlefrontal ','parsopercularis ', 'parstriangularis ','parsorbitalis ' ,'precentral ', 'paracentral ','frontalpole ');
% Par = strvcat('superiorparietal ','inferiorparietal ', 'supramarginal ', 'postcentral ', 'precuneus ');
% Temp = strvcat('superiortemporal ','middletemporal ','inferiortemporal ','bankssts ','fusiform ', 'transversetemporal ','entorhinal ','temporalpole ', 'parahippocampal ');
% Occ = strvcat('lateraloccipital ', 'lingual ','cuneus ','pericalcarine ');
% Cin = strvcat( 'rostralanteriorcingulate ','caudalanteriorcingulate ', 'posteriorcingulate ','isthmuscingulate ');
% Ins = 'insula ';

Fron = strvcat('superiorfrontal','caudalmiddlefrontal','rostralmiddlefrontal','parsopercularis', 'medialorbitofrontal','parstriangularis','parsorbitalis' ,'precentral', 'paracentral','frontalpole','rostralanteriorcingulate','caudalanteriorcingulate');
Par = strvcat('superiorparietal','inferiorparietal', 'supramarginal', 'postcentral', 'precuneus','posteriorcingulate','isthmuscingulate');
Temp = strvcat('superiortemporal','middletemporal','inferiortemporal','bankssts','fusiform', 'transversetemporal','entorhinal','temporalpole', 'parahippocampal');
Occ = strvcat('lateraloccipital', 'lingual','cuneus','pericalcarine');
Ins = 'insula';
%% =======================End of Defining Lobes ======================== %%



%% ================= Removing white spaces from Structure Names ======== %%
totnames1 = char(colortable.struct_names);
totnames = '';
for i = 1:size(totnames1,1)
    totnames = strvcat(totnames,deblank(totnames1(i,:)))
end
ind = find(sum(isspace(totnames))==size(colortable.struct_names,1));
totnames(:,ind) = [];
%% ========== End of Removing white spaces from Structure Names ======== %%

%% ================= Removing white spaces from Frontal Lobe =========== %%
ind = find(sum(isspace(Fron))==size(Fron,1));
Fron(:,ind) = [];
ind = ismember(totnames(:,1:size(Fron,2)),Fron,'rows');
FronIds = colortable.table(ind,5);
indfron = ismember(txt,FronIds);
%% ========== End of Removing white spaces from Frontal Lobe =========== %%

%% ================= Removing white spaces from Parietal Lobe ========== %%
ind = find(sum(isspace(Par))==size(Par,1));
Par(:,ind) = [];
ind = ismember(totnames(:,1:size(Par,2)),Par,'rows');
ParIds = colortable.table(ind,5);
indpar = ismember(txt,ParIds);
%% ========== End of Removing white spaces from Parietal Lobe ========== %%

%% ================= Removing white spaces from Temporal Lobe ========== %%
ind = find(sum(isspace(Temp))==size(Temp,1));
Temp(:,ind) = [];
ind = ismember(totnames(:,1:size(Temp,2)),Temp,'rows');
TempIds = colortable.table(ind,5);
indtemp = ismember(txt,TempIds);
%% ========== End of Removing white spaces from Temporal Lobe ========== %%

%% ================= Removing white spaces from Occipital Lobe ========= %%
ind = find(sum(isspace(Occ))==size(Occ,1));
Occ(:,ind) = [];
ind = ismember(totnames(:,1:size(Occ,2)),Occ,'rows');
OccIds = colortable.table(ind,5);
indocc = ismember(txt,OccIds);
%% ========== End of Removing white spaces from Occipital Lobe ========= %%

%% =============== Creating New Label Id vector ======================== %%
txt1 = txt*0;
txt1(indfron) = 1;
txt1(indpar) = 2;
txt1(indtemp) = 3;
txt1(indocc) = 4;
indcc = find(txt1==0);
%% ========== End of Creating New Label Id vector ====================== %%
%% =============== Creating Colors ===================================== %%
Surftl.Is = txt1;
Surftl.SurfData.FaceVertexCData = zeros(size(Surftl.SurfData.vertices,1),3);
color = [66 85 132]/255;
Surftl.SurfData.FaceVertexCData(indfron,:) = repmat(color,[sum(indfron) 1]);
color = [226 176 149]/255;
Surftl.SurfData.FaceVertexCData(indpar,:) = repmat(color,[sum(indpar) 1]);
color = [95 200 129]/255;
Surftl.SurfData.FaceVertexCData(indtemp,:) = repmat(color,[sum(indtemp) 1]);
color = [105 173 204]/255;
Surftl.SurfData.FaceVertexCData(indocc,:) = repmat(color,[sum(indocc) 1]);
color = [255 255 255]/255;
Surftl.SurfData.FaceVertexCData(indcc,:) = repmat(color,[length(indcc) 1]);
Surftl.SurfData.FaceColor = 'interp';

return;