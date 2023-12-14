function [fitness,Rpd,Rsquare,ind_v,beta,process_par,otherout]=plsr_WN_n(train_X,train_Y,test_X,test_Y,chrome_i,CE4_D10_dat,otherin)
%Rmsep Rpd Rsquare 主成分数量 beta 过程量    
%otherin: 是辅助输入变量，如这里是1×3向量，第一个是方法，第二个是权重，第三个是画图标志位
%otherout: 是辅助输出变量，便于后期更改程序
%method：1是目标函数是test的rmsep
%        2是目标函数是train的rmsep
%        3是目标函数是all的rmsep
%        4是目标函数是train和test的rmsep的加权
%chrome_i是对应的二进制染色体
%CE4_D10_dat 390×23的矩阵，前389行是导数，最后一行是特定波段的反射率。是在主程序中幅值的全局变量，为了避免复杂的参数传递。全局变量在使用的时候要先“调用”，即通过下面语句：
otherout=0;
method_rms=otherin(1);
weight_rms=otherin(2);
method_fit=otherin(3);%fitness计算是方法1rms,2CE4D10,3平均值,4加权,5是1和2的阈值切换,6在case1和case4之间阈值切换
weight_fit=otherin(4);%fitn_method=4时有效，weight是每种元素的rms加权(vs CE4D10 stab)
CE4_D10_data=CE4_D10_dat(:,find(chrome_i>0));

    x=train_X';
    [n,p]=size(x);

    if p>15
        cmpt=15;  % 用于下面交叉验证的主成分数量
    else
        cmpt=7;
    end
    maxcmpt=min(n-1,p);
    if cmpt > maxcmpt-1  % ??????????????????????????plsregress(line 144)
        cmpt=maxcmpt-1;   
    end
    if cmpt>5
        cmpt=4;
    end
    if size(x,1)~=size(train_Y,1)
        train_Y=train_Y';
        test_Y=test_Y';
    end
    %使用n折交叉验证，n=min(训练样本数量,数据维度)
    %选用上面程序挑选的n时，程序速度较慢，但是20代迭代，每次运行的结果相对稳定
    %n=10时，程序耗时缩减了80%，但是20代迭代，每次运行的结果离散度较大，有时很好，有时较差
    [~,~,~,~,~,pctvar,plsmsecv]=plsregress(x,train_Y,cmpt,'cv',n); 
    cmpt_min=2;
    limit=cmpt_min;%解决没有合适的limit，导致limit无法识别报错，2021.2.27
    pctvar_min_persent=0.04;
    for ii=cmpt_min:cmpt
        s1=0;
        for jj=1:ii
         s1=s1+pctvar(2,jj);
        end
        s2=s1-pctvar(2,ii);
        dii=(s1-s2)/s2; %s1是第1到ii项的解释方差总和，s2是第1到ii-1项解释方差的总和
        if dii>=pctvar_min_persent   %选到某一个主成分能够解释光谱中的方差占它之前的成分解释方差总量的比例小于pctvar_min_persent
         limit=ii;
        else 
         continue;
        end
    end
    
     xx=test_X';
     [q1,r]=size(xx);
     [q2,r]=size(x);
     rmsecv=sqrt(plsmsecv);
     %indv是主成分数量，选择依据有两条：
     %1.最后一个主成分能够解释的方差占它前面总和的pctvar_min_persent以上
     %2.在这前limit+1个主成分中选择最终rmsep最小的数量的主成分
     [~,ind_v]=min(rmsecv(2,2:limit+1)); % 用于PLSR的主成分数量从最小的预测均方误差中选择
     if ind_v<3
         ind_v=2; %设置最小主成分为2
     end
     [~,~,~,~,beta,~,~,~]=plsregress(x,train_Y,ind_v);
     %验证集
     Yfit_test=[ones(q1,1),xx]*beta;
     tss=sum((test_Y-mean(test_Y)).^2);
     rss=sum((test_Y- Yfit_test).^2);
     Rsquare_test=1-rss/tss;
     Rmsep_test=sqrt(rss/q1);
     st=std(test_Y);
     Rpd_test=st/Rmsep_test;
     test_par=[Rsquare_test,Rmsep_test];
     %训练集
     Yfit_train=[ones(q2,1),x]*beta;
     tss1=sum((train_Y-mean(train_Y)).^2);
     rss1=sum((train_Y- Yfit_train).^2);
     Rsquare_train=1-rss1/tss1;
     Rmsep_train=sqrt(rss1/(q2));
     st=std(train_Y);
     Rpd_train=st/Rmsep_train;
     train_par=[Rsquare_train,Rmsep_train];
     %全部
     Yfit_all=[ones((q1+q2),1),[x;xx]]*beta;
     tss1=sum(([train_Y;test_Y]-mean([train_Y;test_Y])).^2);
     rss1=sum(([train_Y;test_Y]- Yfit_all).^2);
     Rsquare_all=1-rss1/tss1;
     Rmsep_all=sqrt(rss1/(q1+q2));
     st=std([train_Y;test_Y]);
     Rpd_all=st/Rmsep_all;
     all_par=[Rsquare_all,Rmsep_all];
     %  The stability on CE4 D10
     CE4_D10_fit=[ones(length(CE4_D10_data(:,1)),1),CE4_D10_data]*beta;
     D10_std=std(CE4_D10_fit);
     
     %对训练集和验证集效果进行加权，作为迭代依据
%      weight是训练集加权，(1-wegt)是验证集加权
     Rsquare_avg=Rsquare_train*weight_rms+Rsquare_test*(1-weight_rms);
     Rmsep_avg=Rmsep_train*weight_rms+Rmsep_test*(1-weight_rms);
     Rpd_avg=Rpd_train*weight_rms+Rpd_test*(1-weight_rms);
     
     switch method_rms
         case 1 %目标函数是test的rmsep
             Rmsep=Rmsep_test;Rpd=Rpd_test;Rsquare=Rsquare_test;
         case 2 %目标函数是train的rmsep
             Rmsep=Rmsep_train;Rpd=Rpd_train;Rsquare=Rsquare_train;
         case 3 %目标函数是all的rmsep
             Rmsep=Rmsep_all;Rpd=Rpd_all;Rsquare=Rsquare_all;
         case 4 %目标函数是train和test的rmsep的加权
             Rmsep=Rmsep_avg;Rpd=Rpd_avg;Rsquare=Rsquare_avg;
         case 5 %设置阈值，在case1和case2之间切换
             if (Rsquare_train-Rsquare_test)>0.001
                 Rmsep=Rmsep_test;Rpd=Rpd_test;Rsquare=Rsquare_test;
             else
                 Rmsep=Rmsep_train;Rpd=Rpd_train;Rsquare=Rsquare_train;
             end
         case 6 %设置阈值在case2和case4之间切换
             if (Rsquare_train-Rsquare_test)>0.001
                 Rmsep=Rmsep_avg;Rpd=Rpd_avg;Rsquare=Rsquare_avg;
             else
                 Rmsep=Rmsep_train;Rpd=Rpd_train;Rsquare=Rsquare_train;
             end
     end
    
     fit_coe=Rmsep*weight_fit+D10_std*(1-weight_fit);
     switch method_fit %1rms,2CE4D10 std,3平均值,4加权,5是1和2的阈值切换,6在case1和case4之间阈值切换
         case 1 %fitness是Rmsep
             fitness=Rmsep;
         case 2 %fitness是CE4D10 std
             fitness=D10_std;
         case 3 %fitness是平均值
             fitness=(Rmsep+D10_std)/2;
         case 4 %fitness是加权均值
             fitness=fit_coe;
         case 5 %设置阈值，在case1和case2之间切换
             if (Rmsep-D10_std)>0.001
                 fitness=Rmsep;
             else
                 fitness=D10_std;
             end
         case 6 %设置阈值在case1和case4之间切换
             if (Rmsep*weight_fit-D10_std*(1-weight_fit))>0.001
                 fitness=Rmsep;
             else
                 fitness=fit_coe;
             end
     end
     
     %汇总参数
     process_par=[train_par,test_par,all_par,ind_v,D10_std,Rmsep,fitness];
     %训练集参数2（R2,Rmsep），验证集参数2，全部数据的参数2，主成分个数, CE4D10 std，用于计算fitness的rmsep，fitness

end