function create_gyrisulciAnnot(IdFile)

IdfIle = '/media/Data/PROCESSING_RESULTS/HCP/5-freesurfer_processing/HCPrest_Ids.txt';
SubjectsDIr = '/media/Data/PROCESSING_RESULTS/HCP/5-freesurfer_processing/'
annotId = 'gyrisulci';

Ids = char(textread(IdfIle,'%s'));
Ns = size(Ids,1);
for  i = 1:Ns 
    Id = deblank(Ids(i,:));

    % Left Hemisphere
    load([SubjectsDIr filesep Id filesep 'surf' filesep 'lh.lines.crowns.depth.mat'])
    Surf = Read_Surface([SubjectsDIr filesep Id filesep 'surf' filesep 'lh.white']);
    sulgyrParcfile = [SubjectsDIr filesep  Id filesep 'label' filesep 'lh.' annotId '.annot'];
    opts.mwallind = find(WatershParcell == 4000);
    
    
    [Trip] = Vert_Neibp(double(Surf.SurfData.faces),size(Surf.SurfData.vertices,1),size(Surf.SurfData.faces,1));
    Temp = sum(Trip);
    Trip(:,Temp==0) = [];
    temp = Trip(:,3:end);
    Niter = 4;
    gyrsulparc2 = zeros(size(Surf.SurfData.vertices,1),1);
    indc = gCrowns(:,1:2);
    indc = unique(indc(:));
    for i = 1:Niter
        gyrsulparc2(indc) = 2;
        tempNeigh = temp(indc,:);
        tempNeigh = nonzeros(unique(tempNeigh(:)));
        indc = tempNeigh;
    end
    
    inds = find(gyrsulparc2 == 0);
    gyrsulparc2(inds) = 1;
    % indcorrect = find(curvlh<0);
    % gyrsulparc2(indcorrect)  = 1;
    gyrsulparc = gyrsulparc2;
    gyrsulparc(opts.mwallind) = 0;
    
    ctab.numEntries = 3;
    ctab.orig_tab = 'Gyri-Sulci Colortable';
    ctab.struct_names = {'unknown+corpuscallosum';'sulci';'gyri'};
    colors = [255 255 255;20 230 20;230 20 20 ];
    ctab.table = [colors colors(:,1)*0 colors(:,1)+colors(:,2)*2^8+colors(:,3)*2^16];
    gyrsulparctemp = gyrsulparc*0;
    gyrsulparctemp(gyrsulparc == 0) = ctab.table(1,5);
    gyrsulparctemp(gyrsulparc == 1) = ctab.table(2,5);
    gyrsulparctemp(gyrsulparc == 2) = ctab.table(3,5);
    tempAnnot = save_annotfiles(gyrsulparctemp, sulgyrParcfile,ctab);
    
    
    
     % Right Hemisphere
    load([SubjectsDIr filesep Id filesep 'surf' filesep 'rh.lines.crowns.depth.mat'])
    Surf = Read_Surface([SubjectsDIr filesep Id filesep 'surf' filesep 'rh.inflated']);
    sulgyrParcfile = [SubjectsDIr filesep  Id filesep 'label' filesep 'rh.' annotId '.annot'];
    opts.mwallind = find(WatershParcell == 4000);
    
    
    [Trip] = Vert_Neibp(double(Surf.SurfData.faces),size(Surf.SurfData.vertices,1),size(Surf.SurfData.faces,1));
    Temp = sum(Trip);
    Trip(:,Temp==0) = [];
    temp = Trip(:,3:end);
    Niter = 4;
    gyrsulparc2 = zeros(size(Surf.SurfData.vertices,1),1);
    indc = gCrowns(:,1:2);
    indc = unique(indc(:));
    for i = 1:Niter
        gyrsulparc2(indc) = 2;
        tempNeigh = temp(indc,:);
        tempNeigh = nonzeros(unique(tempNeigh(:)));
        indc = tempNeigh;
    end
    
    inds = find(gyrsulparc2 == 0);
    gyrsulparc2(inds) = 1;
    % indcorrect = find(curvlh<0);
    % gyrsulparc2(indcorrect)  = 1;
    gyrsulparc = gyrsulparc2;
    gyrsulparc(opts.mwallind) = 0;
    
    ctab.numEntries = 3;
    ctab.orig_tab = 'Gyri-Sulci Colortable';
    ctab.struct_names = {'unknown+corpuscallosum';'sulci';'gyri'};
    colors = [255 255 255;20 230 20;230 20 20 ];
    ctab.table = [colors colors(:,1)*0 colors(:,1)+colors(:,2)*2^8+colors(:,3)*2^16];
    gyrsulparctemp = gyrsulparc*0;
    gyrsulparctemp(gyrsulparc == 0) = ctab.table(1,5);
    gyrsulparctemp(gyrsulparc == 1) = ctab.table(2,5);
    gyrsulparctemp(gyrsulparc == 2) = ctab.table(3,5);
    tempAnnot = save_annotfiles(gyrsulparctemp, sulgyrParcfile,ctab);
    
    
    
end
return;
