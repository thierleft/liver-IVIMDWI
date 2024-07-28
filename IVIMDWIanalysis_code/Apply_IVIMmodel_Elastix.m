
%   -*- coding: utf-8 -*-
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% 
%   Loop over all patients to apply voxel-wise IVIM model on DWI co-registered series within 3D segmented regions. 
%   A segmented bi-exponential model is fitted to extract IVIM parameters (f, D, and D*) 
%   and single exponential model to extract apparent diffusion coefficients (ADC).
%
%   Not for clinical use.
%   SPDX-FileCopyrightText: 2024 University of Montreal, Montreal, CAN
%   SPDX-FileCopyrightText: 2024 Thierry L. Lefebvre
%   SPDX-FileCopyrightText: 2024 Guillaume Gilbert
%   SPDX-License-Identifier: MIT
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

% Define path to folders containing co-registered DWI series saved as .mat
% Path structure in this study:'YOUR PROJECT FILE PATH\PATIENT ID\DATE\IVIM\matFiles\*.mat'
PathNameIVIM = 'YOUR PROJECT FILE PATH';
files=dir(PathNameIVIM);

% Define identifiers for both analysts having conducted manual segmentations
seg_ID = {'Analyst1','Analyst2'};

% Iterate over each file corresponding to a patient folder
for i=1:length(files)
    
    nii_file=dir([PathNameIVIM files(i).name '/2*']);

    % Iterate over each file corresponding to the date of MRI scans 
    for iter = 1:length(nii_file)
        
        % Obtain directory with co-registered DWI series
        nii_file_sub =dir([PathNameIVIM files(i).name '\' nii_file(iter).name '\IVIM*']);
        PathMat = [PathNameIVIM files(i).name '\' nii_file(iter).name '\' nii_file_sub.name '\matFiles'];

        % Iterate over each segmentation performed by different analysts
        % (here segmentations are saved as .nii and include analyst ID in the file name)
        for iterID = 1:length(seg_ID)

            dirSeg = dir([PathMat, '\*', seg_ID{iterID},'.nii']);
            file_Seg=[PathMat, '\', dirSeg.name];
            dirmat = dir([PathMat, '\*.mat']);

            % Create folder to save calculated ADC and IVIM parameters
            yourFolder = [PathNameIVIM files(i).name '\' nii_file(iter).name '\' nii_file_sub.name '\' 'Analysis\IVIM_',seg_ID{iterID},'_Elastix\'];
            if ~exist(yourFolder,'dir')
                mkdir(yourFolder)
            end

            %**************************************************************************
            % Perform voxel-wise ADC and IVIM parameters calculations within the segmentation
            %**************************************************************************
            [mean_ROI_f,std_ROI_f,mean_ROI_D,std_ROI_D,mean_ROI_P,std_ROI_P,ROI_size,mean_ROI_ADC,std_ROI_ADC]=IVIM_NII_MAP(PathMat,file_Seg,dirmat.name);               
            

            fileID=fopen([yourFolder,'IVIMparameters.txt'],'w');
            fprintf(fileID,'%s %s %s %s %s %s\r\n','mean_f  ','mean_D  ',...
                'mean_P  ','mean_ADC   ', 'ROI_size   ','name'); % Write First Line with the Type of IVIM parameters of each Column
            fprintf(fileID,'%f %f %f %f %f %f %f %s\r\n',mean_ROI_f,mean_ROI_D,mean_ROI_P,mean_ROI_ADC,ROI_size,files(i).name);
            fclose(fileID);

        end

    end
end
