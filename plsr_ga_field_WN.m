function [process_par_final,chrom_final]=plsr_ga_field_WN(custo_bands,train_X,train_Y,test_X,test_Y,CE4_D10_dat,element,nidv,gnrt,gengap,other)
% 2022.12.2 fitness��plsr_WN_n�м��㣬process_par�洢����м����
    if nargin < 6
        error('Not enough input parameter');
    end
    if nargin == 6
        nind=30; % number of individuals of population
        maxgen=100; % maximum number of generation
        ggap=0.9; % generation gap, means how many new individuals are created
    end
    if nargin == 7 % ??????????
       nind=nidv;
       maxgen=100;
       ggap=0.9;
    end
    if nargin == 8  % ?????????????????
        nind=nidv;
        maxgen=gnrt;
        ggap=0.9;
    end
    if nargin >= 9  % ������������9�����ϣ���ô���ô���Ĳ���
        nind=nidv;
        maxgen=gnrt;
        ggap=gengap;
    end
    wavelength=train_X(:,1); % train_X�ĵ�һ���ǲ���
    train_X=train_X(:,2:end); % train_X�����������������������ǲ��θ���
    [m,~]=size(train_X);%m��X������������������
    % initialize population
    chrom_n=find(custo_bands<0);%ǿ�Ʋ��Ŵ�����
    chrom_p=find(custo_bands>0);%ǿ���Ŵ�����
    chrom_o=find(custo_bands==0);%��ͨ�������ڱ��죩
    nvar=length(chrom_o);%���б���Ļ�������
    chrom=crtbp(nind,nvar);%nind�У�Ⱦɫ���������nvar���У�������������chrome�Ƕ����Ʊ��/���� �����ʼ����Ⱥ
    gen=1;
    %process_parǰ2����ѵ����������2������֤��������2����ȫ������,���зֱ��Ӧ���̲���:Rsquare,Rmsep,��7�������ɷָ�����
    %��8��CE4D10 std����9��10�зֱ���rmsep_sele(���ڼ���fitness��rmsep)��fitness
    %[fitness, process_par, beta] beta=cell(1,20);
    [fitness, process_par, beta]=plsr_obj_field_WN(train_X,train_Y,test_X,test_Y,CE4_D10_dat,chrom,other);%fitness����Ⱥ��rmsep

    chrom_gen=false(nind,m,maxgen);
    beta_gen=cell(maxgen,nind);
    process_par_gen=zeros(nind,10,maxgen);
%     fitn_gen=zeros(maxgen,nind);

    %qual_gen_sbest��ѡ��comb_fit_bestʱ��qual,rmsep_gen_sbestͬ��
    chrom_ful=false(nind,m);
    selch_all=false(round(ggap*nind),length(chrom_p)+length(chrom_o));
    selch_all(:,chrom_p)=true;%ǿ���Ŵ�����
    while gen<=maxgen
        %��¼����
        %fitness_gen�ڼ�����fitness֮��
        process_par_gen(:,:,gen)=process_par;%ѵ��������2��R2,Rmsep������֤������2��ȫ�����ݵĲ���2�����ɷָ���
        beta_gen(gen,:)=beta;   %����������Ⱦɫ���Ӧ��Betaֵ
        chrom_ful(:,chrom_o)=chrom;
        chrom_ful(:,chrom_n)=false;
        chrom_ful(:,chrom_p)=true;
        chrom_gen(:,:,gen)=chrom_ful;%����������Ⱦɫ��ķ���

%         fitn_gen(gen,:)=fitness';    
        fitnv=ranking(fitness);
%         [~,best_ind]=min(comb_fit);
        % select individuals for breeding selch������gap*20
        selch=select('sus',chrom,fitnv,ggap);%'sus'��Stochastic Universal Sampling(���ȫ�����)������ѡ������ӽ���Ⱦɫ�� Add the chief which bears the best fitness to the selch. 2022.7.22
%         chrom_best=chrom(best_ind,:);
        % crossover/recombine individuals
        selch=recombin('xovsp',selch,0.7);%�����ӽ���ĺ����������ԭ����ͬ��'xovsp'��single-point crossover����
        %0.7�ǰ����ڳɶԸ���֮�䷢������/����ĸ��ʵı�����selch������gap*20
        
        % apply mutation�����죩
        selch=mut(selch,0.01);%�Ը����ĸ��ʶ�ÿ��Ⱦɫ����б��� selch������gap*20��Ĭ�ϸ�����0.7/bands=0.18%��mut(selch,0.1)���趨�������Ϊ10%
%         selch(1,:)=chrom_best;% 
        % evaluate offspring, call objective function
        selch_all(:,chrom_o)=selch;    
        [fitness_sel,process_par_ofs,beta_ofs]=plsr_obj_field_WN(train_X,train_Y,test_X,test_Y,CE4_D10_dat,selch_all,other);

        % reinsert offspring into population
        [chrom,~,Index]=reins(chrom,selch,1,1,fitness,fitness_sel);%Index����Ӧ������gap������Ⱦɫ���λ��
        fitness(Index)=fitness_sel;%���Ӵ���fitnessֵ���뵽ԭ����beta��
        beta(1,Index)=beta_ofs;%���Ӵ���betaֵ���뵽ԭ����beta��
        process_par(Index,:)=process_par_ofs;%���Ӵ���process_par���뵽ԭ��������

        % increment counter
        gen=gen+1;
       
    end
    save(strcat(element,'_process_par.mat'),'wavelength','process_par_gen','beta_gen','chrom_gen');
    In_process_par(process_par_gen,beta_gen,chrom_gen,wavelength,element,maxgen);%
    process_par_final_all=process_par_gen(:,:,maxgen);
    [~,I]=min(process_par_final_all(:,10));
    process_par_final=process_par_final_all(I,:);%����fitnessֵ����һ��Ⱦɫ��
    chrom_final=chrom_ful(I,:);
%     bands=chrom2bands(wavelength,chrom_ful);%����chrome�����ݣ���ȡʵ�ʶ�Ӧ�Ĳ���λ��
end