function [outfilename] = Unified_ctx(infilename,afile,delopt);
%
% Syntax :
% [outfilename] = Unified_ctx(infilename,delopt);
%
% This script was developed to join freesurfer ctx-GM and ctx-WM regions.
% The ctx-WM regions are generated using the command mri_aparc2aseg from
% freesurfer. This command allows the growing of gray matter regions into
% the white matter
%
% Input Parameters:
%   infilename     : Atlas filename
%   afile          : Freesurfer segmentation file
%   delopt         : Delete the original atlas file
%
% Output Parameters:
%   outfilename    : Output atlas filename.
%
% See also:
%__________________________________________________
% Authors:  Yasser Aleman Gomez
% LIM
% February 27th 2012
% Version $1.0
if nargin<3
    delopt = 0;
end
if nargin<2
    afile = 'aparc+aseg';
end
switch afile
    case 'aparc+aseg'
        V = spm_vol(infilename);
        IA = spm_read_vols(V);
        ind = find((IA(:)>=3000)&(IA(:)<3100));  %Left Hemisphere
        ind2 = find((IA(:)>=4000)&(IA(:)<4100)); %Right Hemisphere
        ind = [ind;ind2];
        It =IA;It(ind) = It(ind)-2000;
        [pth,nm,ext] = fileparts(infilename);
        if delopt == 0
            V.fname = [pth filesep nm '_gwj.nii'];
        end
        spm_write_vol(V,It);
        outfilename = V.fname;
    case 'a2005s+aseg'
        V = spm_vol(infilename);
        IA = spm_read_vols(V);
        ind = find((IA(:)>=3100)&(IA(:)<3200));  %Left Hemisphere
        ind2 = find((IA(:)>=4100)&(IA(:)<4200)); %Right Hemisphere
        ind = [ind;ind2];
        It =IA;It(ind) = It(ind)-2000;
        [pth,nm,ext] = fileparts(infilename);
        if delopt == 0
            V.fname = [pth filesep nm '_gwj.nii'];
        end
        spm_write_vol(V,It);
        outfilename = V.fname;
    case 'a2009s+aseg'
        V = spm_vol(infilename);
        IA = spm_read_vols(V);
        ind = find((IA(:)>=13100)&(IA(:)<13200));  %Left Hemisphere
        ind2 = find((IA(:)>=14100)&(IA(:)<14200)); %Right Hemisphere
        ind = [ind;ind2];
        It =IA;It(ind) = It(ind)-2000;
        [pth,nm,ext] = fileparts(infilename);
        if delopt == 0
            V.fname = [pth filesep nm '_gwj.nii'];
        end
        spm_write_vol(V,It);
        outfilename = V.fname;
    case 'yeoh7networks'
        V = spm_vol(infilename);
        IA = spm_read_vols(V);
        ind = find((IA(:)>=3000)&(IA(:)<3999));  %Left Hemisphere
        ind2 = find((IA(:)>=4000)&(IA(:)<4999)); %Right Hemisphere
        ind = [ind;ind2];
        It =IA;It(ind) = It(ind)-2000;
        [pth,nm,ext] = fileparts(infilename);
        if delopt == 0
            V.fname = [pth filesep nm '_gwj.nii'];
        end
        spm_write_vol(V,It);
        outfilename = V.fname;
    case 'yeoh17networks'
        V = spm_vol(infilename);
        IA = spm_read_vols(V);
        ind = find((IA(:)>=3000)&(IA(:)<3999));  %Left Hemisphere
        ind2 = find((IA(:)>=4000)&(IA(:)<4999)); %Right Hemisphere
        ind = [ind;ind2];
        It =IA;It(ind) = It(ind)-2000;
        [pth,nm,ext] = fileparts(infilename);
        if delopt == 0
            V.fname = [pth filesep nm '_gwj.nii'];
        end
        spm_write_vol(V,It);
        outfilename = V.fname;
    case 'parckmeans'
        V = spm_vol(infilename);
        IA = spm_read_vols(V);
        ind = find((IA(:)>=3000)&(IA(:)<3999));  %Left Hemisphere
        ind2 = find((IA(:)>=4000)&(IA(:)<4999)); %Right Hemisphere
        ind = [ind;ind2];
        It =IA;It(ind) = It(ind)-2000;
        [pth,nm,ext] = fileparts(infilename);
        if delopt == 0
            V.fname = [pth filesep nm '_gwj.nii'];
        end
        spm_write_vol(V,It);
        outfilename = V.fname;
end
return;
