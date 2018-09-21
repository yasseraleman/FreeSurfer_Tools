function FreeSurf_2_Circos;
%
% Syntax :
% [GMcodes] = Gray_Matter_codes(atlastype);
%
% This function extract gray matter codes from a specified atlas.
%
% Input Parameters:
%   atype        : Atlas type.
%
%
% Output Parameters:
%   GMcodes       : Gray matter codes.
%
%
% Related references:
%
%
% See also: Stats_2_Circos
%
%
%__________________________________________________
% Authors: Yasser Aleman Gomez
% LIM, HUGGM
% March 20th 2012
% Version $1.0

%% ================ Checking Input Parameters =============================%
CharFiles = strvcat(['/media/Data/PROCESSING_RESULTS/NUEVRECLU/5-freesurfer_processing/NUREC_00009__109-20070213/stats' filesep 'volum2circos.txt'],...
                    ['/media/Data/PROCESSING_RESULTS/NUEVRECLU/5-freesurfer_processing/NUREC_00009__109-20070213/stats' filesep 'thick2circos.txt'],...
                    ['/media/Data/PROCESSING_RESULTS/NUEVRECLU/5-freesurfer_processing/NUREC_00009__109-20070213/stats' filesep 'curv2circos.txt'],...
                    ['/media/Data/PROCESSING_RESULTS/NUEVRECLU/5-freesurfer_processing/NUREC_00009__109-20070213/stats' filesep 'area2circos.txt']);
ConnectFile = '/media/Data/PROCESSING_RESULTS/NUEVRECLU/7-connectome/NUREC_00009__109-20070213/tractres/probtrack/NUREC_00009__109-20070213-Connectivity_Matrix-aparc+aseg.txt';                
Outdir = '/media/Data/PROCESSING_RESULTS/NUEVRECLU/test';
atype = 'aparc+aseg';
conth = 0.05;
%% =============== End of cheking Input Parameters =======================%

switch atype
    case 'aparc+aseg'
        [GMcodes,LHMCodes,RHMCodes,Names] = Gray_Matter_codes(atype);
        % Selecting Subcortical regions
        SubC = sort([10:13 17:18 26 49:54 58]);
        SubCL = sort([10:13 17:18 26]);
        SubCR = sort([49:54 58]);
        
        % Selecting Frontal regions
        FroIds = sort([1028 1003 1027 1018 1019 1020 1012 1014 1024 1017 1032]);
        FroIdsR = FroIds+1000;
        
        % Selecting Temporal regions
        TempIds = sort([1009 1015 1030 1001 1007 1034 1006 1033 1016]);
        TempIdsR = TempIds+1000;
        
        % Selecting Parietal regions
        ParIds = sort([1029 1008 1031 1022 1025]);
        ParIdsR = ParIds+1000;
        
        % Selecting Occipital regions
        OccIds = sort([1011 1013 1005 1021]);
        OccIdsR = OccIds+1000;
        
        % Selecting Cingulate regions
        CingIds = sort([1002 1010 1023 1026]);
        CingIdsR = CingIds+1000;
        
        % Selecting Insula regions
        InsIds = [1035];
        InsIdsR = [2035];
end

Characts.Names = '';
for i = 1:size(CharFiles,1)
    Statfile = deblank(CharFiles(i,:));
    [col1, col2, col3] = textread(Statfile,'%s%s%s','delimiter',';','headerlines',1);
    Characts.Values(:,i) = str2num(char(col3));
    fid = fopen(Statfile,'rt');
    line = fgetl(fid);
    a = strread(line,'%s','delimiter',' ');
    Characts.Names = strvcat(Characts.Names,char(a{2}));
    fclose(fid);
end
Ids = str2num(char(col2));
ind = ismember(Ids,GMcodes);
indd = find(ind ==0); 
Characts.Values(indd,:) =[]; 
Characts.Ids =Ids(ind); 
Ids(indd) = [];

%% Reading Freesurfer Color LUT
[a,temp] = system('echo $FREESURFER_HOME');
txtfile =  strcat(temp,filesep,'FreeSurferColorLUT.txt');
fid = fopen(txtfile);
cont = 0;lines = '';
while 1
    cont = cont + 1;
    line = fgetl(fid);
    if ~ischar(line),   break,   end
    lines = strvcat(lines,line);
end
ind = find(sum(isspace(lines'))==size(lines,2));lines(ind,:) = [];
Codes = '';
Names = '';
RedCode = '';
GreenCode = '';
BlueCode = '';
for i = 1:size(lines,1)
    if ~strcmp(lower(deblank(lines(i,1))),'#');
        temp = strread(deblank(lines(i,:)),'%s','delimiter',' ');
        code = temp{1};
        name = temp{2};
        Codes = strvcat(Codes,code);
        Names = strvcat(Names,name);
        RedCode = strvcat(RedCode,temp{3});
        GreenCode = strvcat(GreenCode,temp{4});
        BlueCode = strvcat(BlueCode,temp{5});
    end
end
fclose(fid);
%%
fidmap = fopen([Outdir filesep 'map.txt'],'wt');
fidcbrain = fopen([Outdir filesep 'color.brain.conf'],'wt');
fidsegment = fopen([Outdir filesep 'segments.txt'],'wt');
AcronSt_Tot = '';

%% Frontal Structures selection
ind = ismember(str2num(Codes),FroIds);
FroNames = Names(ind,:);
AcronSt  = Names_2_Acron('aparc+aseg',FroNames);
AcronSt_Tot = AcronSt;
FroRedCode = RedCode(ind,:);
FroGreenCode = GreenCode(ind,:);
FroBlueCode = BlueCode(ind,:);

% Selecting Left Hemisphere
ind = ismember(Characts.Ids,FroIds);
FroVal = Characts.Values(ind,:);
indcharl = find(ind); % Selecting Ids from frontal lobe to create measure files

% Selecting Right Hemisphere
ind = ismember(Characts.Ids,FroIdsR);
indcharr = find(ind); % Selecting Ids from frontal lobe to create measure files
%
Ns = size(AcronSt,1);
fprintf(fidsegment,'%s\n',['chr - fro-l Fro 0 ' num2str(Ns*100-1) ' black']); % Saving Chromosome and Bands, Left Hemisphere
chrcad = ['chr - fro-r Fro 0 ' num2str(Ns*100-1) ' black']; % Saving Chromosome and Bands, Right Hemisphere
meascadleft = '';
meascadright = '';
for i= 1:Ns
    cad = ['Fro ' deblank(AcronSt(i,:)) ' '  deblank(FroRedCode(i,:)) ' ' deblank(FroGreenCode(i,:)) ' ' deblank(FroBlueCode(i,:))];
    % Creating lower Acronysms
    Nacro = lower(deblank(AcronSt(i,:)));ind = strfind(Nacro,'-');
    Nacro(ind) = [];
    cadc = [deblank(Nacro) ' = '  deblank(FroRedCode(i,:)) ',' deblank(FroGreenCode(i,:)) ',' deblank(FroBlueCode(i,:))];
    fprintf(fidsegment,'%s\n',['band fro-l ' deblank(AcronSt(i,:)) ' ' deblank(AcronSt(i,:)) ' ' num2str((i-1)*100) ' ' num2str((i-1)*100+99) ' ' Nacro]); % Saving Chromosome and Bands
    meascadleft = strvcat(meascadleft,['fro-l ' num2str((i-1)*100) ' ' num2str((i-1)*100+99)]); % Preparing to save Measures Files
    chrcad = strvcat(chrcad,['band fro-r ' deblank(AcronSt(i,:)) ' ' deblank(AcronSt(i,:)) ' ' num2str((i-1)*100) ' ' num2str((i-1)*100+99) ' ' Nacro]); % Saving Chromosome and Bands, Right
    meascadright = strvcat(meascadright,['fro-r ' num2str((i-1)*100) ' ' num2str((i-1)*100+99)]); % Preparing to save Measures Files
    % Creating Characteristics string
    cadn = '';
    for j = 1:size(FroVal,2)
        cadn = [cadn ' ' num2str(FroVal(i,j)) ];
    end
    cadt = [cad cadn];
    fprintf(fidmap,'%s\n',cadt);  % Saving map.txt
    fprintf(fidcbrain,'%s\n',cadc);  % Saving color.brain.conf
end
fprintf(fidmap,'\n');

%% Parietal selection
ind = ismember(str2num(Codes),ParIds);
FroNames = Names(ind,:);
AcronSt  = Names_2_Acron('aparc+aseg',FroNames);
for i = 1:size(AcronSt,1)
    AcronSt_Tot = strvcat(AcronSt_Tot,deblank(AcronSt(i,:)));
end

FroRedCode = RedCode(ind,:);
FroGreenCode = GreenCode(ind,:);
FroBlueCode = BlueCode(ind,:);

% Selecting Left Hemisphere
ind = ismember(Characts.Ids,ParIds);
FroVal = Characts.Values(ind,:);
indcharl = [indcharl;find(ind)]; % Selecting Ids from Parietal lobe to create measure files

% Selecting Right Hemisphere
ind = ismember(Characts.Ids,ParIdsR);
indcharr = [indcharr;find(ind)]; % Selecting Ids from Parietal lobe to create measure files
%

Ns = size(AcronSt,1);
fprintf(fidsegment,'%s\n',['chr - par-l Par 0 ' num2str(Ns*100-1) ' black']); % Saving Chromosome and Bands, Left Hemisphere
chrcad = strvcat(chrcad,['chr - par-r Par 0 ' num2str(Ns*100-1) ' black']); % Saving Chromosome and Bands, Right Hemisphere
for i= 1:size(AcronSt,1)
    cad = ['Par ' deblank(AcronSt(i,:)) ' '  deblank(FroRedCode(i,:)) ' ' deblank(FroGreenCode(i,:)) ' ' deblank(FroBlueCode(i,:))];
    % Creating lower Acronysms
    Nacro = lower(deblank(AcronSt(i,:)));ind = strfind(Nacro,'-');Nacro(ind) = [];
    cadc = [deblank(Nacro) ' = '  deblank(FroRedCode(i,:)) ',' deblank(FroGreenCode(i,:)) ',' deblank(FroBlueCode(i,:))];
    fprintf(fidsegment,'%s\n',['band par-l ' deblank(AcronSt(i,:)) ' ' deblank(AcronSt(i,:)) ' ' num2str((i-1)*100) ' ' num2str((i-1)*100+99) ' ' Nacro]); % Saving Chromosome and Bands
    meascadleft = strvcat(meascadleft,['par-l ' num2str((i-1)*100) ' ' num2str((i-1)*100+99)]); % Preparing to save Measures Files
    chrcad = strvcat(chrcad,['band par-r ' deblank(AcronSt(i,:)) ' ' deblank(AcronSt(i,:)) ' ' num2str((i-1)*100) ' ' num2str((i-1)*100+99) ' ' Nacro]); % Saving Chromosome and Bands Right
    meascadright = strvcat(meascadright,['par-r ' num2str((i-1)*100) ' ' num2str((i-1)*100+99)]); % Preparing to save Measures Files
    % Creating Characteristics string
    cadn = '';
    for j = 1:size(FroVal,2)
        cadn = [cadn ' ' num2str(FroVal(i,j)) ];
    end
    cadt = [cad cadn];
    fprintf(fidmap,'%s\n',cadt);  % Saving map.txt
    fprintf(fidcbrain,'%s\n',cadc);  % Saving color.brain.conf
end
fprintf(fidmap,'\n');

%% Temporal selection
ind = ismember(str2num(Codes),TempIds);
FroNames = Names(ind,:);
AcronSt  = Names_2_Acron('aparc+aseg',FroNames);
for i = 1:size(AcronSt,1)
    AcronSt_Tot = strvcat(AcronSt_Tot,deblank(AcronSt(i,:)));
end

FroRedCode = RedCode(ind,:);
FroGreenCode = GreenCode(ind,:);
FroBlueCode = BlueCode(ind,:);

% Selecting Left Hemisphere
ind = ismember(Characts.Ids,TempIds);
indcharl = [indcharl;find(ind)]; % Selecting Ids from Temporal lobe to create measure files
FroVal = Characts.Values(ind,:);

% Selecting Right Hemisphere
ind = ismember(Characts.Ids,TempIdsR);
indcharr = [indcharr;find(ind)]; % Selecting Ids from Temporal lobe to create measure files
%

Ns = size(AcronSt,1);
fprintf(fidsegment,'%s\n',['chr - tem-l Tem 0 ' num2str(Ns*100-1) ' black']); % Saving Chromosome and Bands, Left Hemisphere
chrcad = strvcat(chrcad,['chr - tem-r Tem 0 ' num2str(Ns*100-1) ' black']); % Saving Chromosome and Bands, Right Hemisphere
for i= 1:size(AcronSt,1)
    cad = ['Tem ' deblank(AcronSt(i,:)) ' '  deblank(FroRedCode(i,:)) ' ' deblank(FroGreenCode(i,:)) ' ' deblank(FroBlueCode(i,:))];
    % Creating lower Acronysms
    Nacro = lower(deblank(AcronSt(i,:)));ind = strfind(Nacro,'-');Nacro(ind) = [];
    cadc = [deblank(Nacro) ' = '  deblank(FroRedCode(i,:)) ',' deblank(FroGreenCode(i,:)) ',' deblank(FroBlueCode(i,:))];
    fprintf(fidsegment,'%s\n',['band tem-l ' deblank(AcronSt(i,:)) ' ' deblank(AcronSt(i,:)) ' ' num2str((i-1)*100) ' ' num2str((i-1)*100+99) ' ' Nacro]); % Saving Chromosome and Bands
    meascadleft = strvcat(meascadleft,['tem-l ' num2str((i-1)*100) ' ' num2str((i-1)*100+99)]); % Preparing to save Measures Files
    chrcad = strvcat(chrcad,['band tem-r ' deblank(AcronSt(i,:)) ' ' deblank(AcronSt(i,:)) ' ' num2str((i-1)*100) ' ' num2str((i-1)*100+99) ' ' Nacro]); % Saving Chromosome and Bands Right
    meascadright = strvcat(meascadright,['tem-r ' num2str((i-1)*100) ' ' num2str((i-1)*100+99)]); % Preparing to save Measures Files
    % Creating Characteristics string
    cadn = '';
    for j = 1:size(FroVal,2)
        cadn = [cadn ' ' num2str(FroVal(i,j)) ];
    end
    cadt = [cad cadn];
    fprintf(fidmap,'%s\n',cadt);  % Saving map.txt
    fprintf(fidcbrain,'%s\n',cadc);  % Saving color.brain.conf
end
fprintf(fidmap,'\n');

%% Occipital selection
ind = ismember(str2num(Codes),OccIds);
FroNames = Names(ind,:);
AcronSt  = Names_2_Acron('aparc+aseg',FroNames);
for i = 1:size(AcronSt,1)
    AcronSt_Tot = strvcat(AcronSt_Tot,deblank(AcronSt(i,:)));
end

FroRedCode = RedCode(ind,:);
FroGreenCode = GreenCode(ind,:);
FroBlueCode = BlueCode(ind,:);

% Selecting Left Hemisphere
ind = ismember(Characts.Ids,OccIds);
indcharl = [indcharl;find(ind)]; % Selecting Ids from Occipital lobe to create measure files
FroVal = Characts.Values(ind,:);

% Selecting Right Hemisphere
ind = ismember(Characts.Ids,OccIdsR);
indcharr = [indcharr;find(ind)]; % Selecting Ids from Occipital lobe to create measure files
%

Ns = size(AcronSt,1);
fprintf(fidsegment,'%s\n',['chr - occ-l Occ 0 ' num2str(Ns*100-1) ' black']); % Saving Chromosome and Bands, Left Hemisphere
chrcad = strvcat(chrcad,['chr - occ-r Occ 0 ' num2str(Ns*100-1) ' black']); % Saving Chromosome and Bands, Right Hemisphere
for i= 1:size(AcronSt,1)
    cad = ['Occ ' deblank(AcronSt(i,:)) ' '  deblank(FroRedCode(i,:)) ' ' deblank(FroGreenCode(i,:)) ' ' deblank(FroBlueCode(i,:))];
    % Creating lower Acronysms
    Nacro = lower(deblank(AcronSt(i,:)));ind = strfind(Nacro,'-');Nacro(ind) = [];
    cadc = [deblank(Nacro) ' = '  deblank(FroRedCode(i,:)) ',' deblank(FroGreenCode(i,:)) ',' deblank(FroBlueCode(i,:))];
    fprintf(fidsegment,'%s\n',['band occ-l ' deblank(AcronSt(i,:)) ' ' deblank(AcronSt(i,:)) ' ' num2str((i-1)*100) ' ' num2str((i-1)*100+99) ' ' Nacro]); % Saving Chromosome and Bands
    meascadleft = strvcat(meascadleft,['occ-l ' num2str((i-1)*100) ' ' num2str((i-1)*100+99)]); % Preparing to save Measures Files
    chrcad = strvcat(chrcad,['band occ-r ' deblank(AcronSt(i,:)) ' ' deblank(AcronSt(i,:)) ' ' num2str((i-1)*100) ' ' num2str((i-1)*100+99) ' ' Nacro]); % Saving Chromosome and Bands Right
    meascadright = strvcat(meascadright,['occ-r ' num2str((i-1)*100) ' ' num2str((i-1)*100+99)]); % Preparing to save Measures Files
    % Creating Characteristics string
    cadn = '';
    for j = 1:size(FroVal,2)
        cadn = [cadn ' ' num2str(FroVal(i,j)) ];
    end
    cadt = [cad cadn];
    fprintf(fidmap,'%s\n',cadt);  % Saving map.txt
    fprintf(fidcbrain,'%s\n',cadc);  % Saving color.brain.conf
end
fprintf(fidmap,'\n');

%% Cingulate selection
ind = ismember(str2num(Codes),CingIds);
FroNames = Names(ind,:);
AcronSt  = Names_2_Acron('aparc+aseg',FroNames);
for i = 1:size(AcronSt,1)
    AcronSt_Tot = strvcat(AcronSt_Tot,deblank(AcronSt(i,:)));
end

FroRedCode = RedCode(ind,:);
FroGreenCode = GreenCode(ind,:);
FroBlueCode = BlueCode(ind,:);

% Selecting Left Hemisphere
ind = ismember(Characts.Ids,CingIds);
indcharl = [indcharl;find(ind)]; % Selecting Ids from Cingulate lobe to create measure files
FroVal = Characts.Values(ind,:);

% Selecting Right Hemisphere
ind = ismember(Characts.Ids,CingIdsR);
indcharr = [indcharr;find(ind)]; % Selecting Ids from Cingulate lobe to create measure files
%

Ns = size(AcronSt,1);
fprintf(fidsegment,'%s\n',['chr - cin-l Cin 0 ' num2str(Ns*100-1) ' black']); % Saving Chromosome and Bands, Left Hemisphere
chrcad = strvcat(chrcad,['chr - cin-r Cin 0 ' num2str(Ns*100-1) ' black']); % Saving Chromosome and Bands, Right Hemisphere
for i= 1:size(AcronSt,1)
    cad = ['Cin ' deblank(AcronSt(i,:)) ' '  deblank(FroRedCode(i,:)) ' ' deblank(FroGreenCode(i,:)) ' ' deblank(FroBlueCode(i,:))];
    % Creating lower Acronysms
    Nacro = lower(deblank(AcronSt(i,:)));ind = strfind(Nacro,'-');Nacro(ind) = [];
    cadc = [deblank(Nacro) ' = '  deblank(FroRedCode(i,:)) ',' deblank(FroGreenCode(i,:)) ',' deblank(FroBlueCode(i,:))];
    fprintf(fidsegment,'%s\n',['band cin-l ' deblank(AcronSt(i,:)) ' ' deblank(AcronSt(i,:)) ' ' num2str((i-1)*100) ' ' num2str((i-1)*100+99) ' ' Nacro]); % Saving Chromosome and Bands
    meascadleft = strvcat(meascadleft,['cin-l ' num2str((i-1)*100) ' ' num2str((i-1)*100+99)]); % Preparing to save Measures Files
    chrcad = strvcat(chrcad,['band cin-r ' deblank(AcronSt(i,:)) ' ' deblank(AcronSt(i,:)) ' ' num2str((i-1)*100) ' ' num2str((i-1)*100+99) ' ' Nacro]); % Saving Chromosome and Bands Right
    meascadright = strvcat(meascadright,['cin-r ' num2str((i-1)*100) ' ' num2str((i-1)*100+99)]); % Preparing to save Measures Files
    % Creating Characteristics string
    cadn = '';
    for j = 1:size(FroVal,2)
        cadn = [cadn ' ' num2str(FroVal(i,j)) ];
    end
    cadt = [cad cadn];
    fprintf(fidmap,'%s\n',cadt);  % Saving map.txt
    fprintf(fidcbrain,'%s\n',cadc);  % Saving color.brain.conf
end
fprintf(fidmap,'\n');

%% Insula selection
ind = ismember(str2num(Codes),InsIds);
FroNames = Names(ind,:);
AcronSt  = Names_2_Acron('aparc+aseg',FroNames);
for i = 1:size(AcronSt,1)
    AcronSt_Tot = strvcat(AcronSt_Tot,deblank(AcronSt(i,:)));
end

FroRedCode = RedCode(ind,:);
FroGreenCode = GreenCode(ind,:);
FroBlueCode = BlueCode(ind,:);

% Selecting Left Hemisphere
ind = ismember(Characts.Ids,InsIds);
indcharl = [indcharl;find(ind)]; % Selecting Ids from Insula lobe to create measure files
FroVal = Characts.Values(ind,:);

% Selecting Right Hemisphere
ind = ismember(Characts.Ids,InsIdsR);
indcharr = [indcharr;find(ind)]; % Selecting Ids from Insula lobe to create measure files
%

Ns = size(AcronSt,1);
fprintf(fidsegment,'%s\n',['chr - ins-l Ins 0 ' num2str(Ns*100-1) ' black']); % Saving Chromosome and Bands, Left Hemisphere
chrcad = strvcat(chrcad,['chr - ins-r Ins 0 ' num2str(Ns*100-1) ' black']); % Saving Chromosome and Bands, Right Hemisphere
for i= 1:size(AcronSt,1)
    cad = ['Ins ' deblank(AcronSt(i,:)) ' '  deblank(FroRedCode(i,:)) ' ' deblank(FroGreenCode(i,:)) ' ' deblank(FroBlueCode(i,:))];
    % Creating lower Acronysms
    Nacro = lower(deblank(AcronSt(i,:)));ind = strfind(Nacro,'-');Nacro(ind) = [];
    cadc = [deblank(Nacro) ' = '  deblank(FroRedCode(i,:)) ',' deblank(FroGreenCode(i,:)) ',' deblank(FroBlueCode(i,:))];
    fprintf(fidsegment,'%s\n',['band ins-l ' deblank(AcronSt(i,:)) ' ' deblank(AcronSt(i,:)) ' ' num2str((i-1)*100) ' ' num2str((i-1)*100+99) ' ' Nacro]); % Saving Chromosome and Bands
    meascadleft = strvcat(meascadleft,['ins-l ' num2str((i-1)*100) ' ' num2str((i-1)*100+99)]); % Preparing to save Measures Files
    chrcad = strvcat(chrcad,['band ins-r ' deblank(AcronSt(i,:)) ' ' deblank(AcronSt(i,:)) ' ' num2str((i-1)*100) ' ' num2str((i-1)*100+99) ' ' Nacro]); % Saving Chromosome and Bands Right
    meascadright = strvcat(meascadright,['ins-r ' num2str((i-1)*100) ' ' num2str((i-1)*100+99)]); % Preparing to save Measures Files
    % Creating Characteristics string
    % Creating Characteristics string
    cadn = '';
    for j = 1:size(FroVal,2)
        cadn = [cadn ' ' num2str(FroVal(i,j)) ];
    end
    cadt = [cad cadn];
    fprintf(fidmap,'%s\n',cadt);  % Saving map.txt
    fprintf(fidcbrain,'%s\n',cadc);  % Saving color.brain.conf
end
fprintf(fidmap,'\n');

%% Subcortical Structures selection
ind = ismember(str2num(Codes),SubC);
FroNames = Names(ind,:);
AcronSt  = Names_2_Acron('aparc+aseg',FroNames);
for i = 1:size(AcronSt,1)
    AcronSt_Tot = strvcat(AcronSt_Tot,deblank(AcronSt(i,:)));
end

FroRedCode = RedCode(ind,:);
FroGreenCode = GreenCode(ind,:);
FroBlueCode = BlueCode(ind,:);

% Selecting Left Hemisphere
ind = ismember(Characts.Ids,SubCL);
indcharl = [indcharl;find(ind)]; % Selecting Ids from Subcortical lobe to create measure files
FroVal = Characts.Values(ind,:);

% Selecting Right Hemisphere
ind = ismember(Characts.Ids,SubCR);
indcharr = [indcharr;find(ind)]; % Selecting Ids from Subcortical lobe to create measure files
%

Ns = size(AcronSt,1);
fprintf(fidsegment,'%s\n',['chr - sbc-l SbC 0 ' num2str(Ns*100-1) ' black']); % Saving Chromosome and Bands, Left Hemisphere
chrcad = strvcat(chrcad,['chr - sbc-r SbC 0 ' num2str(Ns*100-1) ' black']); % Saving Chromosome and Bands, Right Hemisphere
for i= 1:size(AcronSt,1)
    cad = ['SbC ' deblank(AcronSt(i,:)) ' '  deblank(FroRedCode(i,:)) ' ' deblank(FroGreenCode(i,:)) ' ' deblank(FroBlueCode(i,:))];
    % Creating lower Acronysms
    Nacro = lower(deblank(AcronSt(i,:)));ind = strfind(Nacro,'-');Nacro(ind) = [];
    cadc = [deblank(Nacro) ' = '  deblank(FroRedCode(i,:)) ',' deblank(FroGreenCode(i,:)) ',' deblank(FroBlueCode(i,:))];
    fprintf(fidsegment,'%s\n',['band sbc-l ' deblank(AcronSt(i,:)) ' ' deblank(AcronSt(i,:)) ' ' num2str((i-1)*100) ' ' num2str((i-1)*100+99) ' ' Nacro]); % Saving Chromosome and Bands
    meascadleft = strvcat(meascadleft,['sbc-l ' num2str((i-1)*100) ' ' num2str((i-1)*100+99)]); % Preparing to save Measures Files
    chrcad = strvcat(chrcad,['band sbc-r ' deblank(AcronSt(i,:)) ' ' deblank(AcronSt(i,:)) ' ' num2str((i-1)*100) ' ' num2str((i-1)*100+99) ' ' Nacro]); % Saving Chromosome and Bands Right
    meascadright = strvcat(meascadright,['sbc-r ' num2str((i-1)*100) ' ' num2str((i-1)*100+99)]); % Preparing to save Measures Files
    % Creating Characteristics string
    % Creating Characteristics string
    cadn = '';
    for j = 1:size(FroVal,2)
        cadn = [cadn ' ' num2str(FroVal(i,j)) ];
    end
    cadt = [cad cadn];
    fprintf(fidmap,'%s\n',cadt);  % Saving map.txt
    fprintf(fidcbrain,'%s\n',cadc);  % Saving color.brain.conf
end
for i = 1:size(chrcad,1)
        fprintf(fidsegment,'%s\n',deblank(chrcad(i,:)));  % Saving Chromosome and Bands, Right Hemisphere
end
fclose(fidmap);
fclose(fidcbrain);
fclose(fidsegment);

%% Saving Chormosomes order
fidchr = fopen([Outdir filesep 'segment.order.conf'],'wt');
fprintf(fidchr,'%s',['chromosomes_order = fro-r,par-r,tem-r,occ-r,cin-r,ins-r,sbc-r,sbc-l,ins-l,cin-l,occ-l,tem-l,par-l,fro-l']);
fclose(fidchr);

%% Saving Measures Files
indchar = [indcharl;indcharr];
meascad = cat(1, meascadleft,meascadright);
for i = 1:size(CharFiles,1)
    fidm = fopen([Outdir filesep 'measure.' num2str(i-1) '.txt'],'wt');
    for j = 1:size(meascad,1)
        fprintf(fidm,'%s\n',[deblank(meascad(j,:)) ' ' num2str(Characts.Values(indchar(j),i))]);
    end
    fclose(fidm);
end
AcronSt_Tot = [AcronSt_Tot;AcronSt_Tot];

fidl = fopen([Outdir filesep 'structure.label.txt'],'wt');
for i = 1:size(AcronSt_Tot,1)
    fprintf(fidl,'%s\n',[deblank(meascad(i,:)) ' ' deblank(AcronSt_Tot(i,:))]);
end
fclose(fidl);
%% =============== Creating Links between different zones =============== %
Connect = Read_Connectivity_Matrix(ConnectFile,0);
Connect.Matrix = Connect.Matrix/max(Connect.Matrix(:));

% Thresholding Connectivity Matrix
ind = find(Connect.Matrix<conth);
Connect.Matrix(ind) = 0;

% Selecting the upper right corner of the matrix
ind = find(Connect.Matrix);
[X,Y] = ind2sub(size(Connect.Matrix,1),ind);
ind = find([Y-X]>0);
Y(ind) = [];
X(ind) = [];
ind = sub2ind(size(Connect.Matrix),X,Y);
Connect.Matrix(ind)=0;

% 
fidmln = fopen([Outdir filesep 'map.links.txt'],'wt');
fidln = fopen([Outdir filesep 'links.txt'],'wt');
for i = 1:size(X,1)
    typec = 1;
    inds1 = find(indchar==X(i));
    inds2 = find(indchar==Y(i));
    z1 = meascad(inds1,:);
    z2 = meascad(inds2,:);
    fprintf(fidln,'%s\n',[deblank(z1) ' ' deblank(z2) ' type=' num2str(typec) ',score=' num2str(Connect.Matrix(Y(i),X(i)))]);
    if ~isempty(strfind(z1,'-l'))
        hemi1 = 'l';
    elseif ~isempty(strfind(z1,'-r'))
        hemi1 = 'r';
    end
    if ~isempty(strfind(z2,'-l'))
        hemi2 = 'l';
    elseif ~isempty(strfind(z2,'-r'))
        hemi2 = 'r';
    end
    cad = [hemi1 ' ' deblank(AcronSt_Tot(inds1,:)) ' ' hemi2 ' ' deblank(AcronSt_Tot(inds2,:)) ' ' num2str(typec) ' ' num2str(Connect.Matrix(Y(i),X(i)))];
    fprintf(fidmln,'%s\n',cad);
end
fclose(fidln);
fclose(fidmln);
return;

%% =========== End of Creating Links between different zones ============ %
return;

function [GMcodes,LHMCodes,RHMCodes,Names] = Gray_Matter_codes(atlastype);
%
% Syntax :
% [GMcodes] = Gray_Matter_codes(atlastype);
%
% This function extract gray matter codes from a specified atlas.
%
% Input Parameters:
%   atlastype     : Atlas type.
%
%
% Output Parameters:
%   GMcodes       : Gray matter codes.
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
% March 20th 2012
% Version $1.0


switch atlastype
    case 'aparc+aseg'
        [a,temp] = system('echo $FREESURFER_HOME');
        txtfile =  strcat(temp,filesep,'FreeSurferColorLUT.txt');
        fid = fopen(txtfile);
        cont = 0;lines = '';
        while 1
            cont = cont + 1;
            line = fgetl(fid);
            if ~ischar(line),   break,   end
            lines = strvcat(lines,line);
        end
        ind = find(sum(isspace(lines'))==size(lines,2));lines(ind,:) = [];
        Codes = '';
        Names = '';
        for i = 1:size(lines,1)
            if ~strcmp(lower(deblank(lines(i,1))),'#');
                temp = strread(deblank(lines(i,:)),'%s','delimiter',' ');
                code = temp{1};
                name = temp{2};
                Codes = strvcat(Codes,code);
                Names = strvcat(Names,name);
            end
        end
        fclose(fid);
        GMcodes = [10:13 17:18 26 49:54 58 1001:1003 1005:1035 2001:2003 2005:2035];
        LHMCodes = [10:13 17:18 26 1001:1003 1005:1035];
        RHMCodes = [49:54 58 2001:2003 2005:2035];
        index = ismember(str2num(Codes),GMcodes);
        Names = Names(index,:);
    case 'a2009s+aseg'
        [a,temp] = system('echo $FREESURFER_HOME');
        txtfile =  strcat(temp,filesep,'FreeSurferColorLUT.txt');
        fid = fopen(txtfile);
        cont = 0;lines = '';
        while 1
            cont = cont + 1;
            line = fgetl(fid);
            if ~ischar(line),   break,   end
            lines = strvcat(lines,line);
        end
        ind = find(sum(isspace(lines'))==size(lines,2));lines(ind,:) = [];
        Codes = '';
        Names = '';
        for i = 1:size(lines,1)
            if ~strcmp(lower(deblank(lines(i,1))),'#');
                temp = strread(deblank(lines(i,:)),'%s','delimiter',' ');
                code = temp{1};
                name = temp{2};
                Codes = strvcat(Codes,code);
                Names = strvcat(Names,name);
            end
        end
        fclose(fid);
        GMcodes = [9:13 17:18 26 48:54 58 11101:11175 12101:12175];
        index = ismember(str2num(Codes),GMcodes);
        Names = Names(index,:);
    case 'a2005s+aseg'
        [a,temp] = system('echo $FREESURFER_HOME');
        txtfile =  strcat(temp,filesep,'FreeSurferColorLUT.txt');
        fid = fopen(txtfile);
        cont = 0;lines = '';
        while 1
            cont = cont + 1;
            line = fgetl(fid);
            if ~ischar(line),   break,   end
            lines = strvcat(lines,line);
        end
        ind = find(sum(isspace(lines'))==size(lines,2));lines(ind,:) = [];
        Codes = '';
        Names = '';
        for i = 1:size(lines,1)
            if ~strcmp(lower(deblank(lines(i,1))),'#');
                temp = strread(deblank(lines(i,:)),'%s','delimiter',' ');
                code = temp{1};
                name = temp{2};
                Codes = strvcat(Codes,code);
                Names = strvcat(Names,name);
            end
        end
        fclose(fid);
        GMcodes = [9:13 17:18 26 48:54 58 1102:1181 2102:2181];
        index = ismember(str2num(Codes),GMcodes);
        Names = Names(index,:);
    case 'yeoh7networks'
        
    case 'yeoh17networks'
        
    case 'parckmeans'
        
    case 'ibaspm116'
        txtfile = which('atlas116.cod');
        [Codes,b] = textread(txtfile,'%u%s','delimiter','=');
        %         Codes = str2num(char(a));
        Names = char(b);
        GMcodes = [1:90];
        index = ismember(Codes,GMcodes);
        Names = Names(index,:);
    case 'ibaspm71'
        txtfile = which('atlas71.cod');
        [Codes,b] = textread(txtfile,'%u%s','delimiter','=');
        %         Codes = str2num(char(a));
        Names = char(b);
        GMcodes = [1     2     4     5     6     7     9    10    11    12    14    15    16    18    19    20 ...
            23    25    26    27    32    33    36    37    38    39    41    50    52    53    54    56    60 ...
            61    62    63    64    67    69    70    72    74    75    76    80    85    88    90    97    98 ...
            99   101   102   108   110   112   114   119   125   130   132   140   145   154   159   164   165 ...
            175   196   203   251];
        index = ismember(Codes,GMcodes);
        Names = Names(index,:);
end
return

function AcronSt  = Names_2_Acron(atype,Names);
switch atype
    case 'aparc+aseg'
        AcronSt = '';
        for i = 1:size(Names,1);
            if strcmp(deblank(Names(i,:)),'Left-Thalamus')%|strcmp(deblank(Names(i,:)),'Right-Thalamus')
                AcronSt = strvcat(AcronSt,'Thalam');
            elseif strcmp(deblank(Names(i,:)),'Left-Thalamus-Proper')%|strcmp(deblank(Names(i,:)),'Right-Thalamus-Proper')
                 AcronSt = strvcat(AcronSt,'ThalProp');
            elseif strcmp(deblank(Names(i,:)),'Left-Caudate')%|strcmp(deblank(Names(i,:)),'Right-Caudate')
                 AcronSt = strvcat(AcronSt,'Caudt');
            elseif strcmp(deblank(Names(i,:)),'Left-Putamen')%|strcmp(deblank(Names(i,:)),'Right-Putamen')
                 AcronSt = strvcat(AcronSt,'Putam');
            elseif strcmp(deblank(Names(i,:)),'Left-Pallidum')%|strcmp(deblank(Names(i,:)),'Right-Pallidum')
                 AcronSt = strvcat(AcronSt,'Pallid');
            elseif strcmp(deblank(Names(i,:)),'Left-Hippocampus')%|strcmp(deblank(Names(i,:)),'Right-Hippocampus')
                 AcronSt = strvcat(AcronSt,'Hippoc');
            elseif strcmp(deblank(Names(i,:)),'Left-Amygdala')%|strcmp(deblank(Names(i,:)),'Right-Amygdala')
                 AcronSt = strvcat(AcronSt,'Amygd');
            elseif strcmp(deblank(Names(i,:)),'Left-Accumbens-area')%|strcmp(deblank(Names(i,:)),'Right-Accumbens-area')
                 AcronSt = strvcat(AcronSt,'AccumbArea');
                 
                 
            elseif strcmp(deblank(Names(i,:)),'ctx-lh-bankssts')%|strcmp(deblank(Names(i,:)),'ctx-rh-bankssts')
                 AcronSt = strvcat(AcronSt,'BanksStS');
            elseif strcmp(deblank(Names(i,:)),'ctx-lh-caudalanteriorcingulate')%|strcmp(deblank(Names(i,:)),'ctx-rh-caudalanteriorcingulate')
                 AcronSt = strvcat(AcronSt,'CaudAntCing');
            elseif strcmp(deblank(Names(i,:)),'ctx-lh-caudalmiddlefrontal')%|strcmp(deblank(Names(i,:)),'ctx-rh-caudalmiddlefrontal')
                 AcronSt = strvcat(AcronSt,'CaudMidFront');
            elseif strcmp(deblank(Names(i,:)),'ctx-lh-cuneus')%|strcmp(deblank(Names(i,:)),'ctx-rh-cuneus')
                 AcronSt = strvcat(AcronSt,'Cuneus');
            elseif strcmp(deblank(Names(i,:)),'ctx-lh-entorhinal')%|strcmp(deblank(Names(i,:)),'ctx-rh-entorhinal')
                 AcronSt = strvcat(AcronSt,'Entorh');
            elseif strcmp(deblank(Names(i,:)),'ctx-lh-fusiform')%|strcmp(deblank(Names(i,:)),'ctx-rh-fusiform')
                 AcronSt = strvcat(AcronSt,'Fusif');
            elseif strcmp(deblank(Names(i,:)),'ctx-lh-inferiorparietal')%|strcmp(deblank(Names(i,:)),'ctx-rh-inferiorparietal')
                 AcronSt = strvcat(AcronSt,'InfPariet'); 
            elseif strcmp(deblank(Names(i,:)),'ctx-lh-inferiortemporal')%|strcmp(deblank(Names(i,:)),'ctx-rh-inferiortemporal')
                 AcronSt = strvcat(AcronSt,'InfTemp'); 
            elseif strcmp(deblank(Names(i,:)),'ctx-lh-isthmuscingulate')%|strcmp(deblank(Names(i,:)),'ctx-rh-isthmuscingulate')
                 AcronSt = strvcat(AcronSt,'IsthCing'); 
            elseif strcmp(deblank(Names(i,:)),'ctx-lh-lateraloccipital')%|strcmp(deblank(Names(i,:)),'ctx-rh-lateraloccipital')
                 AcronSt = strvcat(AcronSt,'LatOccip'); 
            elseif strcmp(deblank(Names(i,:)),'ctx-lh-lateralorbitofrontal')%|strcmp(deblank(Names(i,:)),'ctx-rh-lateralorbitofrontal')
                 AcronSt = strvcat(AcronSt,'LatOrbFront'); 
            elseif strcmp(deblank(Names(i,:)),'ctx-lh-lingual')%|strcmp(deblank(Names(i,:)),'ctx-rh-lingual')
                 AcronSt = strvcat(AcronSt,'Lingual'); 
            elseif strcmp(deblank(Names(i,:)),'ctx-lh-medialorbitofrontal')%|strcmp(deblank(Names(i,:)),'ctx-rh-medialorbitofrontal')
                 AcronSt = strvcat(AcronSt,'MedOrbFront'); 
            elseif strcmp(deblank(Names(i,:)),'ctx-lh-middletemporal')%|strcmp(deblank(Names(i,:)),'ctx-rh-middletemporal')
                 AcronSt = strvcat(AcronSt,'MidTemp'); 
            elseif strcmp(deblank(Names(i,:)),'ctx-lh-parahippocampal')%|strcmp(deblank(Names(i,:)),'ctx-rh-parahippocampal')
                 AcronSt = strvcat(AcronSt,'Parahip'); 
            elseif strcmp(deblank(Names(i,:)),'ctx-lh-paracentral')%|strcmp(deblank(Names(i,:)),'ctx-rh-paracentral')
                 AcronSt = strvcat(AcronSt,'Paracent'); 
            elseif strcmp(deblank(Names(i,:)),'ctx-lh-parsopercularis')%|strcmp(deblank(Names(i,:)),'ctx-rh-parsopercularis')
                 AcronSt = strvcat(AcronSt,'ParsoPerc'); 
            elseif strcmp(deblank(Names(i,:)),'ctx-lh-parsorbitalis')%|strcmp(deblank(Names(i,:)),'ctx-rh-parsorbitalis')
                 AcronSt = strvcat(AcronSt,'ParsOrb'); 
            elseif strcmp(deblank(Names(i,:)),'ctx-lh-parstriangularis')%|strcmp(deblank(Names(i,:)),'ctx-rh-parstriangularis')
                 AcronSt = strvcat(AcronSt,'ParsTriang'); 
            elseif strcmp(deblank(Names(i,:)),'ctx-lh-pericalcarine')%|strcmp(deblank(Names(i,:)),'ctx-rh-pericalcarine')
                 AcronSt = strvcat(AcronSt,'PeriCalc'); 
            elseif strcmp(deblank(Names(i,:)),'ctx-lh-postcentral')%|strcmp(deblank(Names(i,:)),'ctx-rh-postcentral')
                 AcronSt = strvcat(AcronSt,'PostCent'); 
            elseif strcmp(deblank(Names(i,:)),'ctx-lh-posteriorcingulate')%|strcmp(deblank(Names(i,:)),'ctx-rh-posteriorcingulate')
                 AcronSt = strvcat(AcronSt,'PostCing'); 
            elseif strcmp(deblank(Names(i,:)),'ctx-lh-precentral')%|strcmp(deblank(Names(i,:)),'ctx-rh-precentral')
                 AcronSt = strvcat(AcronSt,'Precent'); 
            elseif strcmp(deblank(Names(i,:)),'ctx-lh-precuneus')%|strcmp(deblank(Names(i,:)),'ctx-rh-precuneus')
                 AcronSt = strvcat(AcronSt,'Precuneus'); 
            elseif strcmp(deblank(Names(i,:)),'ctx-lh-rostralanteriorcingulate')%|strcmp(deblank(Names(i,:)),'ctx-rh-rostralanteriorcingulate')
                 AcronSt = strvcat(AcronSt,'RostAntCing'); 
            elseif strcmp(deblank(Names(i,:)),'ctx-lh-rostralmiddlefrontal')%|strcmp(deblank(Names(i,:)),'ctx-rh-rostralmiddlefrontal')
                 AcronSt = strvcat(AcronSt,'RostMidFront'); 
            elseif strcmp(deblank(Names(i,:)),'ctx-lh-superiorfrontal')%|strcmp(deblank(Names(i,:)),'ctx-rh-superiorfrontal')
                 AcronSt = strvcat(AcronSt,'SupFront'); 
            elseif strcmp(deblank(Names(i,:)),'ctx-lh-superiorparietal')%|strcmp(deblank(Names(i,:)),'ctx-rh-superiorparietal')
                 AcronSt = strvcat(AcronSt,'SupPariet'); 
            elseif strcmp(deblank(Names(i,:)),'ctx-lh-superiortemporal')%|strcmp(deblank(Names(i,:)),'ctx-rh-superiortemporal')
                 AcronSt = strvcat(AcronSt,'SupTemp'); 
            elseif strcmp(deblank(Names(i,:)),'ctx-lh-supramarginal')%|strcmp(deblank(Names(i,:)),'ctx-rh-supramarginal')
                 AcronSt = strvcat(AcronSt,'Supramarg'); 
            elseif strcmp(deblank(Names(i,:)),'ctx-lh-frontalpole')%|strcmp(deblank(Names(i,:)),'ctx-rh-frontalpole')
                 AcronSt = strvcat(AcronSt,'FrontPole'); 
            elseif strcmp(deblank(Names(i,:)),'ctx-lh-temporalpole')%|strcmp(deblank(Names(i,:)),'ctx-rh-temporalpole')
                 AcronSt = strvcat(AcronSt,'TempPole'); 
            elseif strcmp(deblank(Names(i,:)),'ctx-lh-transversetemporal')%|strcmp(deblank(Names(i,:)),'ctx-rh-transversetemporal')
                 AcronSt = strvcat(AcronSt,'TransvTemp'); 
            elseif strcmp(deblank(Names(i,:)),'ctx-lh-insula')%|strcmp(deblank(Names(i,:)),'ctx-rh-insula')
                 AcronSt = strvcat(AcronSt,'Insula'); 
            end
        end
end
return;
