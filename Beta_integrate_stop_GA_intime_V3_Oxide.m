%% This program select the top 50 solutions
%% V3: Simplify the code, fix some bugs
% 2023.2.23 Optimized again
% 2023.3.19 Add statistics: 
% average iteration, 1.7 um selection times, average R2, average RMSEP, average SD of CE4 10D
% 2023.9.7 Adjust the stop iteration threshold to better avoid overfitting
% Specify the "result_path", "train_path" and "code_path" before running
% Test under Matlab R2021a
% Date: Feb. 2023
% Author: Zhenxing Zhao, China, Beijing, NSSC, Chinese Academy of Sciecnes
% (2020-2025)
% If you have any questions, please contact cooperzhaozx@gmail.com
% My website (Chinese): cooperzzx.com

close all
clear
restoredefaultpath;
code_path='*\GA-PLSR';
train_path="*\GA-PLSR\T02\"; 
result_path="*\GA-PLSR\T02_Result\"; 
addpath (code_path)
file_train=[train_path,'0.6-1'];
addpath(file_train);
load result_all
load ([code_path,'\train_data'],'ref_1band')
cd (train_path);
dirOutput=dir('*.*');
fileNames={dirOutput.name}';
date_i={dirOutput.date};
Is_dir={dirOutput.isdir}';
fileNames(1:2)=[]; % the first two of fileNames represent the parent directory and the subordinate directory
date_i(1:2)=[];
Is_dir(1:2)=[];
[~,date_num]=sort(date_i); % sort according time order. The order of the file names that MATLAB automatically fetches does not correspond to common sense
Is_dir=Is_dir(date_num);
fileNames=fileNames(date_num);
Is_dir=cell2mat(Is_dir);
fileNames=fileNames(Is_dir); 
fileNum = size(fileNames,1);
Names=elements;
sele_num=50; % the number of best solutions
thr_stopGA=[-0.040,-0.04,-0.03,-0.0130,-0.05]; % the threshold of 'fitness' to stop GA training
% T01 [-0.040,-0.022,-0.016,-0.0130,-0.030]
% T02 [-0.022,-0.022,-0.019,-0.0180,-0.017]
% T03 [-0.040,-0.022,-0.016,-0.0130,-0.030];
thr_n=20; % this indicates that after 'thr_n' iterations, if the 'fitness' change does not exceed 'thr_stop_GA', then the iteration is stopped
disp(['------thr_stopGA = ',num2str(thr_stopGA)]);

%%
par_avg=zeros(5,length(Names));
tic
for i=1:length(Names)
    % aggregate all the results in the initial training
    cd(train_path);
    cd(fileNames{1});
    dataname=join([Names(i),'process_par.mat'],'_');
    load(dataname);
    iter_n=length(beta_gen(:,1)); % iteration times
    fitn_gen=permute(process_par_gen(:,10,:),[1,3,2])';
    chrom_n=length(fitn_gen(1,:));% number of chromosomes/sulutions
    bands=length(wavelength);
    
    stop_gen_all=zeros(fileNum,1); % iteration times when stop training
    beta_all=cell(fileNum,chrom_n); % beta coefficients of PLSR
    chrom_all=zeros(chrom_n,bands,fileNum); % all the chromosomes
    fitness_all=zeros(fileNum,chrom_n); % all the 'fitness' of GA
    par_all=zeros(chrom_n,length(process_par_gen(1,:,1)),fileNum); %[training set(R2,Rmsep), validation set(R2,Rmsep), whole dataset(R2,Rmsep), Number of principal components] 
    for j=1:fileNum % traverse through all the subfolders
        cd(train_path);
        cd(fileNames{j});
        dataname=join([Names(i),'process_par.mat'],'_');
        load(dataname);
        par_best=zeros(iter_n,length(process_par_gen(1,:,1)));
        fitn_gen=permute(process_par_gen(:,10,:),[1,3,2])';
        [fit_min,bestID_gen]=min(fitn_gen,[],2);

        % strategy to stop the iteration 
        % the minimum 'fitness' is used as the criterion to determine the condition for stopping iteration. 2023.9.7
        temp=fit_min(thr_n+1:length(bestID_gen))-fit_min(1:length(bestID_gen)-thr_n);
        temp2=find(temp<thr_stopGA(i));
        temp3=1:length(temp2);
        temp4=find((temp2-temp3')==0);
        stop_gen=length(temp4)+thr_n/2;
        
       % summary parameters
        stop_gen_all(j)=stop_gen;
        beta_all(j,:)=beta_gen(stop_gen,:);
        chrom_all(:,:,j)=chrom_gen(:,:,stop_gen);
        fitness_all(j,:)=fitn_gen(stop_gen,:);
        par_all(:,:,j)=process_par_gen(:,:,stop_gen); 
    end %end of：for j=1:fileNum
    disp(['------Finish all the files! Total ',num2str(fileNum),' files!']);
    
    % plot the changes in minimum and average fitness as a function of the number of iterations
    f=figure();
    subplot(2,1,1)
    plot(stop_gen_all);
    ylabel('GA stop iterations');
    xlabel('File numbers');
    subplot(2,1,2)
    plot(min(fitness_all'));
    hold on
    plot(mean(fitness_all'));
    ylabel('Fitness');
    legend('min','mean');
    xlabel('File numbers');
    sgtitle(Names(i));
    set(gcf, 'Color', 'w','Position',[200,50,1000,700]);
    cd(result_path);
    saveas(f,strcat(Names(i),'_iter.fig'));
    exportgraphics(f,strcat(Names(i),'_iter.jpg'),'Resolution',300)
   
    % further pick solutions among all the results and drop duplicates
    u_fitn=unique(fitness_all);
    stop_gen_final=zeros(fileNum,1);
    sele_map=zeros(fileNum,chrom_n);
    sele_fit_final=zeros(sele_num,1);
    sele_chrom_final=zeros(sele_num,bands);
    sele_beta_final=cell(sele_num,1);
    sele_par_final=zeros(sele_num,length(process_par_gen(1,:,1)));
    for j=1:sele_num 
        [row_j,col_j]=find(fitness_all==u_fitn(j));
        x=row_j(1); % for duplicates, select only one
        y=col_j(1); % for duplicates, select only one
        stop_gen_final(j)=stop_gen_all(x);
        sele_map(x,y)=1; % mapping the IDs of selections among the files
        sele_fit_final(j)=u_fitn(j);
        sele_chrom_final(j,:)=chrom_all(y,:,x);
        sele_beta_final(j)=beta_all(x,y);
        sele_par_final(j,:)=par_all(y,:,x);    
    end
    
    % shown the selection map
    figure();
    sb1=subplot(1,2,1);
    imagesc(fitness_all)
    grid on
    xlabel('Chromosomes');
    ylabel('Files(gaps)')
    title('Fitness');
    colormap(sb1, 'jet');
    c1=colorbar;
    c1.Label.String='Fitness';
    sb2=subplot(1,2,2);
    imagesc(sele_map);
    colormap(sb2,[0.65 0.65 0.65;0.3 0.9 0.3]);
    grid on
    xlabel('Chromosomes');
    ylabel('Files(gaps)')
    title(join(['Top',num2str(sele_num)]));
    sgtitle(strcat(Names(i),' files'))
    c2=colorbar;
    c2.Label.String='Selection';
    set(gcf, 'Color', 'w','Position',[200,20,1000,1050]);
    saveas(gcf,strcat(Names(i),'_files.fig'));
    exportgraphics(gcf,strcat(Names(i),'_files.jpg'),'Resolution',300)
    disp('Saving Element_sele_i......');
    save(strcat(Names(i),'_sele.mat'),'sele_num','stop_gen_final','sele_map','sele_fit_final','sele_chrom_final','sele_beta_final','sele_par_final','thr_stopGA','wavelength');
    sele_par_mean=mean(sele_par_final)';
    par_avg_i=[mean(stop_gen_final);sele_par_mean([7 3 4 8])];
    par_avg(:,i)=par_avg_i;
    
    % show the training process of the best one
    cd(train_path);
    [row_j,col_j]=find(fitness_all==u_fitn(1));%the best fitness
    cd(fileNames{row_j(1)});
    dataname=join([Names(i),'process_par.mat'],'_');
    load(dataname); 
    cd(result_path);
    gen_i=stop_gen_all(row_j(1));
    chrom_loc=[char(Names(i)),'-file',fileNames{row_j(1)},'-chrom',num2str(col_j(1))];
    In_process_par(process_par_gen,beta_gen,chrom_gen,wavelength,Names(i),gen_i);
    exportgraphics(gcf,['In_process_best_chrom',chrom_loc,'.jpg'],'Resolution',300)
    savefig(gcf,['In_process_best_chrom',chrom_loc,'.fig'])
end
Par_Pro=["Iterations";"PC";"R2_Test";"RMSEP_Test";"STD"];
Par_Pro=table(Par_Pro);
par_avg_tab=array2table(par_avg);
par_avg_tab=[Par_Pro,par_avg_tab];
par_avg_tab.Properties.VariableNames=["Proteries",Names];
writetable(par_avg_tab,'par_avg.xlsx');
%% calculate the average 'beta' (weight coefficients of PLSR) of the topmost solutions
cd(train_path);
cd(fileNames{1});
load('result_all')
cd(result_path);
solu=[1 5 10 25 50]; %the top 1 5 10 25 50 solutions
ele_n=length(Names);
Beta_all_ele=zeros(length(solu),length(wavelength)+1,ele_n);
bands_all_count=zeros(length(solu),length(wavelength),ele_n);
chrom_num_id=5; %1,2,3,4,5: 1 5 10 25 50 of solutions
figure();
for i=1:ele_n
    filename=join([Names(i),'sele.mat'],'_');
    load(filename);
    indx=logical(sele_chrom_final);
    beta_all_sel=zeros(length(solu),length(wavelength)+1);
    bands_count=zeros(length(solu),length(wavelength));
    for j=1:length(solu)
        Beta_temp=zeros(solu(j),length(wavelength)+1);
        bands_temp=zeros(1,length(wavelength));
        for k=1:solu(j) % calculate the average Beta value
            Beta_temp(k,1)=sele_beta_final{k}(1);
            Beta_temp(k,[logical(0), indx(k,:)])=sele_beta_final{k}(2:end);
            bands_temp=bands_temp+sele_chrom_final(k,:);
        end
        beta_all_sel(j,:)=mean(Beta_temp,1);
        bands_count(j,:)=bands_temp;
    end
    bands_all_count(:,:,i)=bands_count;
    Beta_all_ele(:,:,i)=beta_all_sel;
    % show the selection frequency of different bands and the averaged beta coefficients
    subplot(ele_n,2,2*i-1);
    bar(wavelength,bands_count(chrom_num_id,:),'FaceColor',[0.2,0.3,0.8]);%,'color',C(1,:)
    hold on
    stem(ref_1band,bands_count(chrom_num_id,end),'Marker','s','Color',[1,0.2,0.2],...
                'MarkerEdgeColor',[1,0.2,0.2],'MarkerFaceColor',[1,0.2,0.2],'MarkerSize',6);
    ax=gca;
    axis padded
    ax.YLim=[0,50];
    ylabel('Bands count');
    set(gca,'FontName','Helvetica-Narrow','FontSize',12,'linewidth',1.2);
    title(Names(i));
    text(2.44,-5,'μm','FontSize',12,'HorizontalAlignment','left','VerticalAlignment','middle');
    RemoveSubplotWhiteArea(gca,ele_n,2,i,1,[0.97,0.99]);
    subplot(ele_n,2,2*i);
    ax=gca;
    plot(wavelength(1:end-1),beta_all_sel(chrom_num_id,2:end-1),'Color',[0.8,0.6,0.2],'LineWidth',1.4);
    hold on
    plot(ref_1band,beta_all_sel(chrom_num_id,end),'Color',[1,0.2,0.2],'Marker','s','MarkerEdgeColor',[1,0.2,0.2],'MarkerFaceColor',[1,0.2,0.2]); 
    axis padded
    ylabel('beta (wt.%)');
    ylim=ax.YLim;
    axis padded
    tx1=join(['c=',num2str(beta_all_sel(chrom_num_id,1),'%.1f'),' wt.%']);  
    text(wavelength(1),min(ylim)+range(ylim)*1.08,tx1,'FontSize',12,'Color',[0.2 0.2 0.2]);%[0 0.4 0.85]
    set(gca,'FontName','Helvetica-Narrow','FontSize',12,'linewidth',1.2);
    yline(0);
    title(Names(i));
    RemoveSubplotWhiteArea(gca,ele_n,2,i,2,[0.97,0.99]);
end
set(gcf, 'Color', 'w','Position',[200,10,1200,1000]);
saveas(gcf,'bands_count.fig');
exportgraphics(gcf,'bands_count.jpg','Resolution',300)
Beta=permute(Beta_all_ele,[3 2 1]);
Beta_all=Beta(:,:,4); % top 50
disp('------Beta_all is saving ......');
save('Beta_all.mat','Beta_all','Beta','Names','wavelength','Ref','Comp');
restoredefaultpath
toc
