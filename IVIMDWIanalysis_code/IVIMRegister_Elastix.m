%   -*- coding: utf-8 -*-
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% 
%   Image co-registration of Intravoxel Incoherent Motion (IVIM) Diffusion-weighted Images (DWI).
%   Perform two-step rigid and deformable co-registration of 3D image volumes
%   in the 4D DWI series saved in .dcm at increasing b-values using Elastix.
%   Registered series are saved as .nii and .mat files.
%
%   Not for clinical use.
%   SPDX-FileCopyrightText: 2024 University of Montreal, Montreal, CAN
%   SPDX-FileCopyrightText: 2024 Thierry L. Lefebvre
%   SPDX-FileCopyrightText: 2024 Guillaume Gilbert
%   SPDX-License-Identifier: MIT
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%


setenv("PATH", strcat(getenv("PATH"),"YOUR PATH TO ELASTIX"));

% Define path to folders containing DWI series saved as .dcm
% Path structure in this study:'YOUR PROJECT FILE PATH\PATIENT ID\DATE\IVIM\*.dcm'
PathNameIVIM = 'YOUR PROJECT FILE PATH';

files=dir(PathNameIVIM);

% Iterate over each file corresponding to a patient folder
for i = 1:length(files)
    
    nii_file=dir([PathNameIVIM files(i).name '\2*']);

    % Iterate over each file corresponding to the date of MRI scans 
    for iter = 1:length(nii_file)
        
        % Obtain directory with saved IVIM DWI data
        nii_file_sub =dir([PathNameIVIM files(i).name '\' nii_file(iter).name '\IVIM*']);

        %**************************************************************************
        % Select and read DWI series from saved DICOM files
        %**************************************************************************
        
        pathToDCM = [PathNameIVIM files(i).name '\' nii_file(iter).name '\' nii_file_sub.name '\'];
        DirContent = dir([PathNameIVIM files(i).name '\' nii_file(iter).name '\' nii_file_sub.name '\*.dcm']);       

        Counter=1;
        SliceLocation=[];
        bValue=[];
        ImageDiffusion=[];
        
        % Create directory to save co-registered DWI series in .nii and .mat
        yourFolder = [PathNameIVIM files(i).name '\' nii_file(iter).name '\' nii_file_sub.name '\' 'matFiles\'];
        disp(yourFolder)
        if ~exist(yourFolder,'dir')
            mkdir(yourFolder)
        end

        h=waitbar(0,'Reading images');

        for k=1:length(DirContent)
            waitbar(k/length(DirContent),h)
            FileName=DirContent(k).name;

            Current_info=dicominfo([pathToDCM DirContent(k).name]);
            SliceLocation(Counter)=Current_info.SliceLocation;

            if (length(Current_info.DiffusionBValue>1))
                bValue(Counter)=typecast(Current_info.DiffusionBValue,'double');
            else
                bValue(Counter)=Current_info.DiffusionBValue;
            end

            ImageDiffusion(:,:,Counter)=double(dicomread([pathToDCM DirContent(k).name]));
            Counter=Counter+1;

        end
        close(h);

        % Sort DWI series in ascending b-value order
        Diffusion=sort(unique(bValue));
        SlicePosition=sort(unique(SliceLocation));

        Nb_Slices=length(SlicePosition);
        Nb_Diffusion=length(Diffusion);

        idata=zeros(size(ImageDiffusion,1),size(ImageDiffusion,2),Nb_Slices,Nb_Diffusion);

        for k=1:Nb_Slices*Nb_Diffusion

            idata(:,:,find(SlicePosition==SliceLocation(k)),find(Diffusion==bValue(k)))=ImageDiffusion(:,:,k);

        end

        idata_Elastix=zeros(size(idata));
        
        %**************************************************************************
        % Apply 3D co-registration from Elastix model zoo (model 57)
        %**************************************************************************
        
        display(size(idata,4))
        for k=1:size(idata,4)
            display(k)
            if k==1
                idata_Elastix(:,:,:,k)=idata(:,:,:,find(Diffusion==0));
            else
                [~,out]=elastix(squeeze(idata(:,:,:,k)),squeeze(idata(:,:,:,find(Diffusion==0))),[],{'Par0057Rigid.txt','Par0057Bspline.txt'});
                idata_Elastix(:,:,:,k)=out.transformedImages{1};
            end
        end

        save([yourFolder, sprintf('%s', num2str(files(i).name)),'_',sprintf('%s',num2str(nii_file(iter).name)),'_Elastix.mat'],'idata','idata_Elastix','Diffusion')

        % Geometrical transformations before saving to NIfTI
        idata_Elastix = idata_Elastix(:,end:-1:1,end:-1:1,:);
        idata_Elastix = rot90(idata_Elastix,3);

        % Save to .nii using openly available MathWorks File Exchange "Tools for NIfTI and ANALYZE image"
        niiData = make_nii(idata_Elastix, [Current_info.PixelSpacing; Current_info.SliceThickness]',Current_info.ImagePositionPatient');
        disp([yourFolder, sprintf('%s', num2str(files(i).name)),'_',sprintf('%s',num2str(nii_file(iter).name)),'_Elastix.nii'])
        save_nii(niiData,[yourFolder, sprintf('%s', num2str(files(i).name)),'_',sprintf('%s',num2str(nii_file(iter).name)),'_Elastix.nii'])
        
        clear ImageDiffusion idata idata_Elastix bValue niiData DirContent

        disp(['Iteration ' sprintf('%s ', num2str(i)) 'sur ' num2str(length(files)) '  ' sprintf('%s',num2str(files(i).name)) 'Elastix'])

    end

end

