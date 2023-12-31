function [rmsep,rmsep_t,rpd,rsquare,cmpt,bands,count]=plsr_ga_mod_WN(custo_bands,train_X,train_Y,test_X,test_Y,nidv,gnrt,gengap,element,other)
% ??PLSR?GA??????????????????????GA_plsr????????????????????????RMSEP?RPD
% train_X,test_X????????????????????????????????L*??N?
% ??train_X???????test_X???????
% train_Y,test_Y??????????????????????????????????????????????????????????*??
% ??rmsep??????????????bands?GA????????????????????rmsep??????
% nidv???????????????????
% gnrt???????????????
% gengap ?????offspring?parents????
% percent ?????GA???????????????????
% count ????GA????????????????GA????????1???????2-7???
% ??????????????????????8-13????????????
% rmsep,rpd??????percent???????????????????????????rpd
% ????metal????????????
% plsr_ga_mod_4???plsr_ga_mod_3,plsr_ga_mod_2??????Cd???0??????????????????
% 2017?3?14??????plsr_ga_4?????plsr_ga_mod_6?????????????
% 2017?7?4??????plsr_ga_mod_6??????????41??64??????plsr?????plsr??????????????
% 2018?8?19???plsr_ga_mod_7()????????????plsr_ga_mod_12,??plsr_12?plsr_12_n??
% 2018?11?21???plsr_ga_mod_12????????????plsr_ga_mod_WN
% 2022.6.27 custo_bands中的0代表参加遗传变异，1代表强制遗传，-1代表强制不遗传
% 2022.12.2鉴于这个函数的功能已经在之外的函数实现，暂时抛弃这个函数
    if nargin < 4
        error('Not enough input parameter');
    end
    if nargin == 4  
        nidv=30; % number of individuals of population
    end
    if nargin >=9
        ele_str=element;
    end
    [n_ele,~]=size(ele_str);
    rmsep=[];
    bands=[];
    rpd=zeros(1,n_ele);
    rsquare=zeros(1,n_ele);
    rmsep_t=zeros(1,n_ele);
    cmpt=zeros(1,n_ele);
    nb=size(train_X,1);  % ??????????????????????
    count=zeros(2*n_ele+1,nb); % ????????????????
    count(1,:)=train_X(:,1); % 第一行是波长
    for t=1:n_ele  % n_ele可能是1或者7等，由上面的程序判断。这里是等于一
        element=ele_str(t,:);
%         fprintf('This is the %d th element %s \n',k,element)
        [temp_rmsep,temp_bands]=plsr_ga_field_WN(custo_bands,train_X,train_Y,test_X,test_Y,element,nidv,gnrt,gengap,other);%进行遗传算法训练
        % for kb=1:nb % ?????????
        %     count_temp=0;
        %     for kr=1:nidv   % ?????????
        %         if temp_bands(kr,kb)~=0
        %             count_temp=count_temp+1;
        %         end
        %     end
        %     count(t+1,kb)=count_temp; % ???????????
        %     count(t+n_ele+1,kb)=count(t+1,kb)/nidv;
        % end
        count(t+1,:) = sum(temp_bands~=0);
        count(t+n_ele+1,:) = sum(temp_bands~=0)/nidv;

        % train_X_subset=[];  % ????????GA_plsr??????
        % test_X_subset=[];  % ???????????????
        [~,ind]=min(temp_rmsep);  % ind是最小值的索引,找出rmsep最小的，然后用它来做PLSR
        % for kb=1:nb   % ?????????
        %     if temp_bands(mind,kb)~=0   % ?????????RMSEP??????
        %         temp_tri=train_X(kb,:);  % ?????????
        %         temp_tes=test_X(kb,:);
        %         train_X_subset=[train_X_subset;temp_tri];  % ?????????????????????
        %         test_X_subset=[test_X_subset;temp_tes];  
        %     end
        % end
        b_ind = temp_bands(ind,:)~=0;%b_ind是temp_bands的非0值逻辑索引
        PLSR_band=temp_bands(ind,:);
        
        train_X_subset = train_X(b_ind,:);
        test_X_subset  = test_X(b_ind,:);
        %下面将总集和作为验证集，而不是训练集
%         test_X_subset=[train_X_subset test_X_subset];
%         test_Y=[train_Y;test_Y];

        train_X_subset(:,1) = []; % ????????????????????
        [rmsep_t(t),rpd(t),rsquare(t),cmpt(t),beta,process_par]=plsr_WN_n(train_X_subset,train_Y,test_X_subset,test_Y,element,other);
        % rmsep rpd r2 主成分个数
        
%         Beta_plot(b_ind,beta,[train_X,test_X],element);%将被选波段可视化
        %被选波段，beta，数据, 元素名称
        
        rmsep=[rmsep,temp_rmsep]; % ?????????
        bands=[bands;temp_bands]; %bands中存储的是挑选的波段
        rmsep_f=process_par(4);
        Rsquare_f=process_par(3);
        cmpt_f=process_par(7);
        save(strcat(element,'_final.mat'),'PLSR_band','b_ind','beta','rmsep_f','Rsquare_f','cmpt_f')%存储用于PLSR并画图展示的波段，以及模型精度参数等
                                      
    end
end


