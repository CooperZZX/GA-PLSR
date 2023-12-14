%% Test T02: Using A+L+CE data for training, do not reduce the effects of phase angle

% Specify the "resultpath" and "code_path" before running
% Test under Matlab R2021a
% Date: Feb. 2023
% Author: Zhenxing Zhao, China, Beijing, NSSC, Chinese Academy of Sciecnes
% (2020-2025)
% If you have any questions, please contact cooperzhaozx@gmail.com
% My website (Chinese): cooperzzx.com

clear
close all
restoredefaultpath;
code_path='*\GA-PLSR';
resultpath="*\GA-PLSR\T02\"; 
load train_dat.mat % load the data for training
addpath(code_path);
cd(resultpath);

%% model parameters
rmse_method=1; % rms_method=1-5,1:test,2:train,3:all,4:weighting,5:swithc between 1 and 2, 6: switch between case2 and case4
weight_rmse=[0.5 0.5 0.5 0.5 0.5 0.5 0.5]; % the weight of RMSE of the training set. Valid only when rmse_method=4
fitn_method=1;% fitn_method=1-5,1:rmse,2:SD of CE4 D10, 3:0.5-0.5 average,4:weighting,5:swithc between 1 and 2, 6: switch between case2 and case4
weight_fit=[0.32 0.5 0.4 0.4 0.25 0.5 0.5];% the weight of RMSE of the validation set (instead of CE4 D10 SD). Valid only when rmse_method=4
nidv = 20; % number of chorsomes
gnrt = 400; % iteration times 1500 700 400
gengap_all=[0.4 0.5 0.6 0.7 0.8 0.9];% generation gap,the proportion to generate the next generation

%% data preparation
WL=0.45:0.005:2.395;
WL_train=WL';
bands_num=length(WL_train);
CE4_D10_dat=[CE4_D10_diff;CE4_D10_SM(find(WL==ref_1band),:)]';
Spec_sele_diff=[Apo_spec_diff,Luna_spec_diff,CE3_Spec_diff(:,[2,4]),CE5_ASD_diff(:,2),CE5_insitu_diff(:,1)]; % A+L+CE
Spec_all_sele=[Apo_spec_all,Luna_spec_SM,CE3_Spec_cort(:,[2,4]),CE5_ASD_SM(:,2),CE5_Spec_inte(:,1)];
% The spectral derivative, which characterizes the shape of the spectra, is a lossy transform of the original spectra. 
% It is imperative to include a constant term that represents the absolute reflectance of the original spectra to achieve a lossless transform. 
% For this purpose, we selected the the reflectance at 1.7 ¦Ìm (r1.7), which is less affected by the varying phase angles and exhibits a relatively high signal-to-noise ratio in the spectra of CE3, CE4, and CE5.
Ref=[Spec_sele_diff;Spec_all_sele(find(WL==ref_1band),:)]; % 1st derivative spectra (the first 389 bands) + ref_1band spectral reflectance (the last 1 band)
Comp=[Apo_che_all_5;Luna_che_5;CE3_Che_5([2,4],:);CE5_Che_5;CE5_Che_5];
Samp_Names=[Apo_Name;Luna_Name;CE3_site_name([2,4],2);"CE5 ASD";"CE5 D09"];
Samp_N=["A14";"A14";"A14";"A14";"A16";"A16";"A16";"A16";"A16";"A16";"Mare";"Mare";"Mare";"Mare";"Mare";"Mare";"Mare";"Mare";"Mare";...
    "Mare";"Mare";"A14";"A16";"Mare";"A16";"Mare";"Mare";"L";"L";"L";"CE3";"CE3";"CE5";"CE5"];
%% init training
n_samp=length(Comp(:,1));
n_train=23;
n_ele=length(elements);
train_X_all=zeros(bands_num,n_train,n_ele);
test_X_all=zeros(bands_num,n_samp-n_train,n_ele);
train_Y_all=zeros(n_train,n_ele);
test_Y_all=zeros(n_samp-n_train,n_ele);
par_final_all=zeros(6,10);
chrom_final_all=zeros(6,bands_num);
train_id_all=zeros(n_train,n_ele);
test_id_all=zeros(n_samp-n_train,n_ele);
tStart=tic; % time record

% make new foders to save the training data
files_num=[length(gengap_all),nidv];
for j=1:files_num(1)
    gengap=gengap_all(j);
    for i=1:files_num(2)
        filename=strcat(resultpath,num2str(gengap),"-",num2str(i)); mkdir(filename)
    end
end
% split of training and validation set
% sort according to the oxide contents, divide all the in situ data into the training set
test_id5=[];
Ind_all=[];
split_int=3;
test_id0=2:split_int:n_samp-4;
train_id0=1:1:n_samp-4;
train_id0(test_id0)=[];
Tr_id5=zeros(23,5);
Te_id5=zeros(11,5);
CE_train_id=[34; 32; 31]; % APXS data of CE3, insitu spectrum D09 of CE5
CE_test_id=33; % CE5 ASD spectrum (high signal-to-noise ration)
Sele_map=string(zeros(n_samp,5));
for i=1:5
    [Y_sort,Ind] = sort(Comp(1:n_samp-4,i));
    Ind_all(:,i)=Ind;
    Tr_id5(:,i)=[Ind(train_id0);CE_train_id];
    Te_id5(:,i)=[Ind(test_id0);CE_test_id];
    Sele_map(Te_id5(:,i),i)="Validation";
    Sele_map(Tr_id5(:,i),i)="Training";
end
Sele_map_tab=array2table([Samp_Names,Sele_map]);
Sele_map_tab.Properties.VariableNames=["Sample",elements];
writetable(Sele_map_tab,'Train_Vali_names_T02.xlsx');
%% start training
for j=1:files_num(1)
    gengap=gengap_all(j);
    for i=1:files_num(2)
        filename=strcat(resultpath,num2str(gengap),"-",num2str(i)) % mkdir filename
        cd(filename);
        parfor (k = 1:n_ele, 5)  % training using parallel calculation
            fprintf('This is the %d th Mineral %s \n',k,elements(k));
            other=[rmse_method weight_rmse(k) fitn_method weight_fit(k)];  
            Tr_id=Tr_id5(:,k);
            Te_id=Te_id5(:,k);
            train_X=Ref(:,Tr_id);
            train_Y=Comp(Tr_id,k);
            test_X=Ref(:,Te_id);
            test_Y=Comp(Te_id,k);  
            train_id_all(:,k)=Tr_id;
            test_id_all(:,k)=Te_id;
            train_X_all(:,:,k)=train_X;
            test_X_all(:,:,k)=test_X;
            train_Y_all(:,k)=train_Y;
            test_Y_all(:,k)=test_Y;
            tic
            [par_final,chrom_final]=plsr_ga_field_WN(custo_bands_all(k,:),[WL_train,train_X],train_Y,test_X,test_Y,CE4_D10_dat,elements(k),nidv,gnrt,gengap,other);
            toc
            par_final_all(k,:)=par_final;
            chrom_final_all(k,:)=chrom_final;
            
        end
        save('result_all.mat','train_X_all','test_X_all','train_Y_all','test_Y_all','train_id_all','test_id_all',...
            'chrom_final_all', 'par_final_all', 'WL_train', 'Samp_Names', 'Samp_N', 'Comp', 'Ref','Spec_all_sele','elements');
    end
end
%%
toc(tStart)
restoredefaultpath;
