function [process_par_final,chrom_final]=plsr_ga_field_WN(custo_bands,train_X,train_Y,test_X,test_Y,CE4_D10_dat,element,nidv,gnrt,gengap,other)
% 2022.12.2 fitness在plsr_WN_n中计算，process_par存储相关中间变量
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
    if nargin >= 9  % 如果输入参数满9个以上，那么就用传入的参数
        nind=nidv;
        maxgen=gnrt;
        ggap=gengap;
    end
    wavelength=train_X(:,1); % train_X的第一列是波段
    train_X=train_X(:,2:end); % train_X的列数是样本个数，行数是波段个数
    [m,~]=size(train_X);%m是X的行数，即波段数量
    % initialize population
    chrom_n=find(custo_bands<0);%强制不遗传基因
    chrom_p=find(custo_bands>0);%强制遗传基因
    chrom_o=find(custo_bands==0);%普通基因（用于变异）
    nvar=length(chrom_o);%进行变异的基因数量
    chrom=crtbp(nind,nvar);%nind行（染色体个数），nvar是列（波段总数），chrome是二进制表达/索引 随机初始化种群
    gen=1;
    %process_par前2列是训练集，接着2列是验证集，接着2列是全部数据,两列分别对应过程参数:Rsquare,Rmsep,第7列是主成分个数，
    %第8列CE4D10 std，第9、10列分别是rmsep_sele(用于计算fitness的rmsep)，fitness
    %[fitness, process_par, beta] beta=cell(1,20);
    [fitness, process_par, beta]=plsr_obj_field_WN(train_X,train_Y,test_X,test_Y,CE4_D10_dat,chrom,other);%fitness是种群的rmsep

    chrom_gen=false(nind,m,maxgen);
    beta_gen=cell(maxgen,nind);
    process_par_gen=zeros(nind,10,maxgen);
%     fitn_gen=zeros(maxgen,nind);

    %qual_gen_sbest是选择到comb_fit_best时的qual,rmsep_gen_sbest同理
    chrom_ful=false(nind,m);
    selch_all=false(round(ggap*nind),length(chrom_p)+length(chrom_o));
    selch_all(:,chrom_p)=true;%强制遗传基因
    while gen<=maxgen
        %记录参数
        %fitness_gen在计算完fitness之后
        process_par_gen(:,:,gen)=process_par;%训练集参数2（R2,Rmsep），验证集参数2，全部数据的参数2，主成分个数
        beta_gen(gen,:)=beta;   %过程中所有染色体对应的Beta值
        chrom_ful(:,chrom_o)=chrom;
        chrom_ful(:,chrom_n)=false;
        chrom_ful(:,chrom_p)=true;
        chrom_gen(:,:,gen)=chrom_ful;%过程中所有染色体的方案

%         fitn_gen(gen,:)=fitness';    
        fitnv=ranking(fitness);
%         [~,best_ind]=min(comb_fit);
        % select individuals for breeding selch数量是gap*20
        selch=select('sus',chrom,fitnv,ggap);%'sus'是Stochastic Universal Sampling(随机全体抽样)函数，选择进行杂交的染色体 Add the chief which bears the best fitness to the selch. 2022.7.22
%         chrom_best=chrom(best_ind,:);
        % crossover/recombine individuals
        selch=recombin('xovsp',selch,0.7);%返回杂交后的后代，数量与原来相同。'xovsp'是single-point crossover函数
        %0.7是包含在成对个体之间发生重组/交叉的概率的标量。selch数量是gap*20
        
        % apply mutation（变异）
        selch=mut(selch,0.01);%以给定的概率对每个染色体进行变异 selch数量是gap*20，默认概率是0.7/bands=0.18%。mut(selch,0.1)是设定变异概率为10%
%         selch(1,:)=chrom_best;% 
        % evaluate offspring, call objective function
        selch_all(:,chrom_o)=selch;    
        [fitness_sel,process_par_ofs,beta_ofs]=plsr_obj_field_WN(train_X,train_Y,test_X,test_Y,CE4_D10_dat,selch_all,other);

        % reinsert offspring into population
        [chrom,~,Index]=reins(chrom,selch,1,1,fitness,fitness_sel);%Index是适应度最差的gap比例个染色体的位置
        fitness(Index)=fitness_sel;%将子代的fitness值插入到原来的beta中
        beta(1,Index)=beta_ofs;%将子代的beta值插入到原来的beta中
        process_par(Index,:)=process_par_ofs;%将子代的process_par插入到原来的里面

        % increment counter
        gen=gen+1;
       
    end
    save(strcat(element,'_process_par.mat'),'wavelength','process_par_gen','beta_gen','chrom_gen');
    In_process_par(process_par_gen,beta_gen,chrom_gen,wavelength,element,maxgen);%
    process_par_final_all=process_par_gen(:,:,maxgen);
    [~,I]=min(process_par_final_all(:,10));
    process_par_final=process_par_final_all(I,:);%最优fitness值的那一个染色体
    chrom_final=chrom_ful(I,:);
%     bands=chrom2bands(wavelength,chrom_ful);%根据chrome的内容，提取实际对应的波段位置
end