
%   -*- coding: utf-8 -*-
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% 
%   Voxel-wise IVIM model on DWI co-registered series within 3D segmented regions. 
%   A segmented bi-exponential model is fitted to extract IVIM parameters (f, D, and D*) 
%   and single exponential model to extract apparent diffusion coefficients (ADC).
%
%   Not for clinical use.
%   SPDX-FileCopyrightText: 2024 University of Montreal, Montreal, CAN
%   SPDX-FileCopyrightText: 2024 Thierry L. Lefebvre
%   SPDX-FileCopyrightText: 2024 Guillaume Gilbert
%   SPDX-License-Identifier: MIT
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%


function [mean_ROI_f,std_ROI_f,mean_ROI_D,std_ROI_D,mean_ROI_P,std_ROI_P,ROI_size,mean_ROI_ADC,std_ROI_ADC]=IVIM_NII_MAP(PathDicom,PathSegmentation,matfile_Name)

% Load saved co-registered DWI series
load([PathDicom '/' matfile_Name(1:end-3) 'mat'])


% Load manual segmentation 
nii=load_untouch_nii(PathSegmentation);
SegMask=double(nii.img);

if length(size(SegMask))>3
    SegMask = squeeze(SegMask(:,:,:,1));
end

for k=1:size(SegMask,3)
    SegMask(:,:,k)=imrotate(SegMask(:,:,k),90);
end
SegMask = SegMask(:,end:-1:1,end:-1:1);

%**************************************************************************
% Perform voxel-wise ADC and IVIM parameters calculations within the segmentation
%**************************************************************************

NDiffValues=length(Diffusion);
Signal_IVIM=zeros(1,NDiffValues);

bValues=Diffusion';

idata=idata_Elastix;

idata_f=zeros(size(idata,1),size(idata,2),size(idata,3));
idata_D=zeros(size(idata,1),size(idata,2),size(idata,3));
idata_P=zeros(size(idata,1),size(idata,2),size(idata,3));
idata_MSE_Mono=zeros(size(idata,1),size(idata,2),size(idata,3));
idata_MSE_Seg=zeros(size(idata,1),size(idata,2),size(idata,3));
idata_ADC=zeros(size(idata,1),size(idata,2),size(idata,3));

imIn = zeros(size(idata,1),size(idata,2),size(idata,3));

b_threshold=150;
options=optimset('Display','off','Algorithm','levenberg-marquardt','TolFun',eps,'TolX',eps,'MaxIter',100000,'FunValCheck','off');

for k=1:size(idata,3)
    k;
    for m=1:size(idata,2)
        for n=1:size(idata,1)
            imIn(n,m,k) = idata(n,m,k,8);

            % Calculate ADC and IVIM parameters within the segmentation only
            if (SegMask(n,m,k)==1)
                try
                    Signal_IVIM=squeeze(idata(n,m,k,:))';


                    %******************************************************************
                    % Mono-exponential ADC fit 
                    %******************************************************************                    
                    Signal_ADC = Signal_IVIM./Signal_IVIM(7);
                    Signal_ADC = [Signal_ADC(7), Signal_ADC(end-1), Signal_ADC(end)];
                    bValues_ADC = [bValues(7), bValues(end-1), bValues(end)];
                    p=polyfit(bValues_ADC,log(Signal_ADC),1);
                    ADC = -p(1);

                    %******************************************************************
                    % Segmented bi-exponential IVIM (Levenberg-Marquardt)
                    %******************************************************************
                    Signal_IVIM_Mono=Signal_IVIM./Signal_IVIM(1);
                    Signal_IVIM_Mono(bValues<b_threshold)=[];
                    bValues_Mono=bValues;
                    bValues_Mono(bValues<b_threshold)=[];

                    if ~isreal(Signal_IVIM_Mono)
                        Signal_IVIM_Mono = abs(Signal_IVIM_Mono);
                    end
                    
                    % Fit 1 for IVIM parameters, f and D
                    fun_Mono= @(x,xdata)fun1(x,xdata);
                    x0_Mono=[0.1 0.001];
                    [x_Mono,~,~,~,mse_Mono] = nlinfit(bValues_Mono',log(Signal_IVIM_Mono),fun_Mono,x0_Mono,options);

                    f_Mono=1-exp(x_Mono(1));
                    D_Mono=-x_Mono(2);

                    % Fit 2 for IVIM parameter, D* (referred to as P here)
                    x0 = 0.001;
                    ub = 0.5;
                    lb = 0.00001;
                    fun_Seg= @(x,xdata)fun2(x,xdata,f_Mono,D_Mono);
                    if ~isreal(Signal_IVIM)
                        Signal_IVIM = abs(Signal_IVIM);
                    end
                    [x_Seg,mse_Seg,~,~,~] = lsqcurvefit(fun_Seg,x0,bValues',Signal_IVIM./Signal_IVIM(1),lb,ub,options);

                    % Assign calculated parameters in volume 
                    idata_f(n,m,k)=f_Mono;
                    idata_D(n,m,k)=D_Mono;
                    idata_P(n,m,k)=x_Seg;
                    idata_MSE_Mono(n,m,k)=mse_Mono;
                    idata_MSE_Seg(n,m,k)=mse_Seg;
                    idata_ADC(n,m,k) = ADC;

                catch
                    idata_f(n,m,k)=NaN;
                    idata_D(n,m,k)=NaN;
                    idata_P(n,m,k)=NaN;
                    idata_MSE_Mono(n,m,k)=NaN;
                    idata_MSE_Seg(n,m,k)=NaN;
                    idata_ADC(n,m,k) = NaN;
                    disp('--- NaN voxel returned (fit didn''t converge) ---')

                end

            end
        end

    end
end

% Apply upper and lower bounds to the calculated parametric maps
idata_f(isnan(idata_f))=0;
idata_f(isinf(idata_f))=0;
idata_f(idata_f<0)=0;

idata_ADC(isnan(idata_ADC))=0;
idata_ADC(isinf(idata_ADC))=0;
idata_ADC(idata_ADC<0)=0;

idata_D(isnan(idata_D))=0;
idata_D(isinf(idata_D))=0;
idata_D(idata_D<0)=0;

idata_P(isnan(idata_P))=0;
idata_P(isinf(idata_P))=0;
idata_P(idata_P<0)=0;

idata_MSE_Mono(isnan(idata_MSE_Mono))=0;
idata_MSE_Mono(isinf(idata_MSE_Mono))=0;
idata_MSE_Mono(idata_MSE_Mono<0)=0;

idata_MSE_Seg(isnan(idata_MSE_Seg))=0;
idata_MSE_Seg(isinf(idata_MSE_Seg))=0;
idata_MSE_Seg(idata_MSE_Seg<0)=0;

%**************************************************************************
% Define quality map and apply on manual segmentation to obtain final ROI
%**************************************************************************
MSE_prctile_threshold=0.85;

MSE_Mono=idata_MSE_Mono(:);
MSE_Mono=MSE_Mono(:);
MSE_Mono(MSE_Mono==0)=[];
prctile_MSE_Mono=quantile(MSE_Mono,MSE_prctile_threshold);

MSE_Seg=idata_MSE_Seg(:);
MSE_Seg=MSE_Seg(:);
MSE_Seg(MSE_Seg==0)=[];
prctile_MSE_Seg=quantile(MSE_Seg,MSE_prctile_threshold);

MSE_Mono_Mask=SegMask;
MSE_Mono_Mask(idata_MSE_Mono>prctile_MSE_Mono)=0;

idata_RMSE_Mono = sqrt(idata_MSE_Mono)./(max(idata_MSE_Mono(:))-min(MSE_Mono));
RMSE_Mono_Mask = SegMask;
RMSE_Mono_Mask(idata_RMSE_Mono>0.05)=0;

MSE_Seg_Mask=SegMask;
MSE_Seg_Mask(idata_MSE_Seg>prctile_MSE_Seg)=0;

idata_RMSE_Seg = sqrt(idata_MSE_Seg)./(max(idata_MSE_Seg(:))-min(MSE_Seg));
RMSE_Seg_Mask = SegMask;
RMSE_Seg_Mask(idata_RMSE_Seg>0.05)=0;

% Include only regions withing quality map for mean/SD parameters reporting
Mask=SegMask;
Mask(MSE_Mono_Mask==0)=0;
Mask(MSE_Seg_Mask==0)=0;
Mask(RMSE_Seg_Mask==0)=0;
Mask(RMSE_Mono_Mask==0)=0;

%**************************************************************************
% Calculate mean/SD of ADC and IVIM parameters within final ROI
%**************************************************************************
idata_f = abs(idata_f);
ROI_f=idata_f.*Mask;
ROI_f=ROI_f(:);
ROI_f(ROI_f==0)=[];
mean_ROI_f=mean(ROI_f);
std_ROI_f=std(ROI_f);

idata_ADC = abs(idata_ADC);
ROI_ADC=idata_ADC.*Mask;
ROI_ADC=ROI_ADC(:);
ROI_ADC(ROI_ADC==0)=[];
mean_ROI_ADC=mean(ROI_ADC);
std_ROI_ADC=std(ROI_ADC);

ROI_D=idata_D.*Mask;
ROI_D=ROI_D(:);
ROI_D(ROI_D==0)=[];
mean_ROI_D=mean(ROI_D);
std_ROI_D=std(ROI_D);

ROI_P=idata_P.*Mask;
ROI_P=ROI_P(:);
ROI_P(ROI_P==0)=[];
mean_ROI_P=mean(ROI_P);
std_ROI_P=std(ROI_P);

ROI_size=sum(Mask(:) == 1);

end

function fun_Seg=fun2(x,xdata,f_Mono,D_Mono)
    fun_Seg = (1-f_Mono)*exp(-xdata*D_Mono)+f_Mono*exp(-xdata*(x(1)));
    fun_Seg(isnan(fun_Seg))=0;
    fun_Seg(isinf(fun_Seg))=0;
end

function fun_Mono=fun1(x,xdata)
    fun_Mono = x(1)+xdata*x(2);
    fun_Mono(isnan(fun_Mono))=0;
    fun_Mono(isinf(fun_Mono))=0;
end
