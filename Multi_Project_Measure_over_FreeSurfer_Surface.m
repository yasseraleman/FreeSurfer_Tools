function OutNames = Multi_Project_Measure_over_FreeSurfer_Surface(IdFile, FreeSDir, ConnectDir, opts);
%
% Syntax :
%  OutNames = Multi_Project_Measure_over_FreeSurfer_Surface(IdFile, FreeSDir, ConnectDir, opts);
%
% This script projects diffusion measures into freesurfer surfaces. 
% It also normalises measures maps into fsaverage template space
%
%
% Input Parameters:
%   FreeSDir          : FreeSurfer Subjects Output directory.
%   ConnectDir        : FreeSurfer Subjects Output directory.
%    IdFile           : Ids File.
%      opts           : Options
%      opts.norm      : Normalise measures into fsaverage space  
%      opts.boolcomp  : Boolean vector to select which diffusion
%                          coefficient will be computed.
%                         Examples:
%                        -- Fractional Anisotropy    [1 0 0 0 0 0 0 0 0 0] --
%                        -- Mean Diffusivity         [0 1 0 0 0 0 0 0 0 0] --
%                        -- Linear coefficient       [0 0 1 0 0 0 0 0 0 0] --
%                        -- Planar coefficient       [0 0 0 1 0 0 0 0 0 0] --
%                        -- Spherical coefficient    [0 0 0 0 1 0 0 0 0 0] --
%                        -- Volume fraction          [0 0 0 0 0 1 0 0 0 0] --
%                        -- Generalized anisotropy   [0 0 0 0 0 0 1 0 0 0] --
%                        -- Radial anisotropy        [0 0 0 0 0 0 0 1 0 0] --
%                        -- Axial Diffusivity        [0 0 0 0 0 0 0 0 1 0] --
%                        -- Radial Diffusivity       [0 0 0 0 0 0 0 0 0 1] --
%
% Output Parameters:
%     OutNames           : Output Maps
%
% Related references:
%
%
% See also: Project_Measure_over_FreeSurfer_Surface
%
%
%__________________________________________________
% Authors: Yasser Aleman Gomez and Javi Santonja
% LIM, HUGGM
% November 3rd 2014
% Version $1.0

% Idfile = '/media/Data/PROCESSING_RESULTS/5-freesurfer_processing/Connectome_Ids.txt';

if nargin < 4
    opts.norm = 1;
end
if nargin < 3
    ConnectDir = '/media/Data/PROCESSING_RESULTS/7-connectome';
    opts.norm = 1;
end
if nargin < 2
    FreeSDir = '/media/Data/PROCESSING_RESULTS/5-freesurfer_processing';
    ConnectDir = '/media/Data/PROCESSING_RESULTS/7-connectome';
    opts.norm = 1;
end
if nargin < 1
    error('Please enter a valid subject ID File');
    return
end
if exist(IdFile,'file')
    Ids = char(textread(IdFile,'%s'));
else
    Ids = IdFile;
end
setenv('SUBJECTS_DIR',FreeSDir);
Ns = size(Ids,1);
for i = 1:Ns
    subjId = deblank(Ids(i,:));
    disp([' Processing ID: ' subjId '. Subject ' num2str(i) ' of ' num2str(Ns)]);
    if exist([ConnectDir filesep subjId filesep 't1diff'],'dir');
        OutNames = Project_Measure_over_FreeSurfer_Surface(subjId, FreeSDir, ConnectDir, opts);
    end
end
return;