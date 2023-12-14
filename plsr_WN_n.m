function [fitness,Rpd,Rsquare,ind_v,beta,process_par,otherout]=plsr_WN_n(train_X,train_Y,test_X,test_Y,chrome_i,CE4_D10_dat,otherin)
%Rmsep Rpd Rsquare ���ɷ����� beta ������    
%otherin: �Ǹ��������������������1��3��������һ���Ƿ������ڶ�����Ȩ�أ��������ǻ�ͼ��־λ
%otherout: �Ǹ���������������ں��ڸ��ĳ���
%method��1��Ŀ�꺯����test��rmsep
%        2��Ŀ�꺯����train��rmsep
%        3��Ŀ�꺯����all��rmsep
%        4��Ŀ�꺯����train��test��rmsep�ļ�Ȩ
%chrome_i�Ƕ�Ӧ�Ķ�����Ⱦɫ��
%CE4_D10_dat 390��23�ľ���ǰ389���ǵ��������һ�����ض����εķ����ʡ������������з�ֵ��ȫ�ֱ�����Ϊ�˱��⸴�ӵĲ������ݡ�ȫ�ֱ�����ʹ�õ�ʱ��Ҫ�ȡ����á�����ͨ��������䣺
otherout=0;
method_rms=otherin(1);
weight_rms=otherin(2);
method_fit=otherin(3);%fitness�����Ƿ���1rms,2CE4D10,3ƽ��ֵ,4��Ȩ,5��1��2����ֵ�л�,6��case1��case4֮����ֵ�л�
weight_fit=otherin(4);%fitn_method=4ʱ��Ч��weight��ÿ��Ԫ�ص�rms��Ȩ(vs CE4D10 stab)
CE4_D10_data=CE4_D10_dat(:,find(chrome_i>0));

    x=train_X';
    [n,p]=size(x);

    if p>15
        cmpt=15;  % �������潻����֤�����ɷ�����
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
    %ʹ��n�۽�����֤��n=min(ѵ����������,����ά��)
    %ѡ�����������ѡ��nʱ�������ٶȽ���������20��������ÿ�����еĽ������ȶ�
    %n=10ʱ�������ʱ������80%������20��������ÿ�����еĽ����ɢ�Ƚϴ���ʱ�ܺã���ʱ�ϲ�
    [~,~,~,~,~,pctvar,plsmsecv]=plsregress(x,train_Y,cmpt,'cv',n); 
    cmpt_min=2;
    limit=cmpt_min;%���û�к��ʵ�limit������limit�޷�ʶ�𱨴�2021.2.27
    pctvar_min_persent=0.04;
    for ii=cmpt_min:cmpt
        s1=0;
        for jj=1:ii
         s1=s1+pctvar(2,jj);
        end
        s2=s1-pctvar(2,ii);
        dii=(s1-s2)/s2; %s1�ǵ�1��ii��Ľ��ͷ����ܺͣ�s2�ǵ�1��ii-1����ͷ�����ܺ�
        if dii>=pctvar_min_persent   %ѡ��ĳһ�����ɷ��ܹ����͹����еķ���ռ��֮ǰ�ĳɷֽ��ͷ��������ı���С��pctvar_min_persent
         limit=ii;
        else 
         continue;
        end
    end
    
     xx=test_X';
     [q1,r]=size(xx);
     [q2,r]=size(x);
     rmsecv=sqrt(plsmsecv);
     %indv�����ɷ�������ѡ��������������
     %1.���һ�����ɷ��ܹ����͵ķ���ռ��ǰ���ܺ͵�pctvar_min_persent����
     %2.����ǰlimit+1�����ɷ���ѡ������rmsep��С�����������ɷ�
     [~,ind_v]=min(rmsecv(2,2:limit+1)); % ����PLSR�����ɷ���������С��Ԥ����������ѡ��
     if ind_v<3
         ind_v=2; %������С���ɷ�Ϊ2
     end
     [~,~,~,~,beta,~,~,~]=plsregress(x,train_Y,ind_v);
     %��֤��
     Yfit_test=[ones(q1,1),xx]*beta;
     tss=sum((test_Y-mean(test_Y)).^2);
     rss=sum((test_Y- Yfit_test).^2);
     Rsquare_test=1-rss/tss;
     Rmsep_test=sqrt(rss/q1);
     st=std(test_Y);
     Rpd_test=st/Rmsep_test;
     test_par=[Rsquare_test,Rmsep_test];
     %ѵ����
     Yfit_train=[ones(q2,1),x]*beta;
     tss1=sum((train_Y-mean(train_Y)).^2);
     rss1=sum((train_Y- Yfit_train).^2);
     Rsquare_train=1-rss1/tss1;
     Rmsep_train=sqrt(rss1/(q2));
     st=std(train_Y);
     Rpd_train=st/Rmsep_train;
     train_par=[Rsquare_train,Rmsep_train];
     %ȫ��
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
     
     %��ѵ��������֤��Ч�����м�Ȩ����Ϊ��������
%      weight��ѵ������Ȩ��(1-wegt)����֤����Ȩ
     Rsquare_avg=Rsquare_train*weight_rms+Rsquare_test*(1-weight_rms);
     Rmsep_avg=Rmsep_train*weight_rms+Rmsep_test*(1-weight_rms);
     Rpd_avg=Rpd_train*weight_rms+Rpd_test*(1-weight_rms);
     
     switch method_rms
         case 1 %Ŀ�꺯����test��rmsep
             Rmsep=Rmsep_test;Rpd=Rpd_test;Rsquare=Rsquare_test;
         case 2 %Ŀ�꺯����train��rmsep
             Rmsep=Rmsep_train;Rpd=Rpd_train;Rsquare=Rsquare_train;
         case 3 %Ŀ�꺯����all��rmsep
             Rmsep=Rmsep_all;Rpd=Rpd_all;Rsquare=Rsquare_all;
         case 4 %Ŀ�꺯����train��test��rmsep�ļ�Ȩ
             Rmsep=Rmsep_avg;Rpd=Rpd_avg;Rsquare=Rsquare_avg;
         case 5 %������ֵ����case1��case2֮���л�
             if (Rsquare_train-Rsquare_test)>0.001
                 Rmsep=Rmsep_test;Rpd=Rpd_test;Rsquare=Rsquare_test;
             else
                 Rmsep=Rmsep_train;Rpd=Rpd_train;Rsquare=Rsquare_train;
             end
         case 6 %������ֵ��case2��case4֮���л�
             if (Rsquare_train-Rsquare_test)>0.001
                 Rmsep=Rmsep_avg;Rpd=Rpd_avg;Rsquare=Rsquare_avg;
             else
                 Rmsep=Rmsep_train;Rpd=Rpd_train;Rsquare=Rsquare_train;
             end
     end
    
     fit_coe=Rmsep*weight_fit+D10_std*(1-weight_fit);
     switch method_fit %1rms,2CE4D10 std,3ƽ��ֵ,4��Ȩ,5��1��2����ֵ�л�,6��case1��case4֮����ֵ�л�
         case 1 %fitness��Rmsep
             fitness=Rmsep;
         case 2 %fitness��CE4D10 std
             fitness=D10_std;
         case 3 %fitness��ƽ��ֵ
             fitness=(Rmsep+D10_std)/2;
         case 4 %fitness�Ǽ�Ȩ��ֵ
             fitness=fit_coe;
         case 5 %������ֵ����case1��case2֮���л�
             if (Rmsep-D10_std)>0.001
                 fitness=Rmsep;
             else
                 fitness=D10_std;
             end
         case 6 %������ֵ��case1��case4֮���л�
             if (Rmsep*weight_fit-D10_std*(1-weight_fit))>0.001
                 fitness=Rmsep;
             else
                 fitness=fit_coe;
             end
     end
     
     %���ܲ���
     process_par=[train_par,test_par,all_par,ind_v,D10_std,Rmsep,fitness];
     %ѵ��������2��R2,Rmsep������֤������2��ȫ�����ݵĲ���2�����ɷָ���, CE4D10 std�����ڼ���fitness��rmsep��fitness

end