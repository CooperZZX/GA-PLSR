%% This program visualize the GA-PLSR model
% V4: select according fitness
% Specify the "result_path", "train_path" and "code_path" before running
% Test under Matlab R2021a
% Date: Dec. 2022
% Author: Zhenxing Zhao, China, Beijing, NSSC, Chinese Academy of Sciecnes
% (2020-2025)
% If you have any questions, please contact cooperzhaozx@gmail.com
% My website (Chinese): cooperzzx.com

close all
clear
restoredefaultpath;
code_path='*\GA-PLSR';
result_path='*\GA-PLSR\T02_Result';
train_path='*\GA-PLSR\T02\0.6-1';
addpath(train_path);
addpath(code_path)
addpath(result_path);
cd (result_path)
load Beta_all
load result_all
names=Names;

%% plot the 'beta' (weight coefficients of PLSR)
f=figure();
for i=1:length(names)
    sp(i)=subplot(3,2,i);
end
c = turbo(240);
for i=1:length(names)
    train_X=train_X_all(:,:,i);
    train_Y=train_Y_all(:,i);
    test_X=test_X_all(:,:,i);
    test_Y=test_Y_all(:,i);
    filename=join([names(i),'sele.mat'],'_');
    load(filename);
    indx=logical(sele_chrom_final);
    [m,n]=size(sele_chrom_final);
    Beta_all=zeros(m,n);
    const=zeros(m,1);
    b_ind=indx(1,:);
    for j=1:sele_num
        const(j)=sele_beta_final{j}(1);
        Beta_all(j,indx(j,:))=sele_beta_final{j}(2:end);
        b_ind=b_ind|indx(j,:);
    end
    const_avg=mean(const);
    if sele_num==1
        beta_avg=Beta_all;
    else
        beta_avg=mean(Beta_all);
    end
    beta=beta_avg(find(beta_avg));
    figure(f);
    subplot(3,2,i);
    ax=gca;colormap(ax,'turbo');
    beta_p=zeros(1,389);
    beta_p(b_ind)=beta;
    WL1=wavelength(b_ind);beta1=beta_p(b_ind);dotsn=length(WL1);
    beta_min=min(beta1); beta_max=max(beta1); 
    offset1=20; offset2=40; 
    scale1=50; scale2=40;
    for j=1:dotsn
        if beta1(j) > 0 % Nonlinear mapping with color
            colrn=round(log1p(beta1(j))/log1p(beta_max)*scale1)+offset1;
        else
            colrn=-round(log1p(-beta1(j))/log1p(-beta_min)*scale2)-offset2;
        end
     colr=c(colrn+121,:);
     stem(WL1(j),beta1(j),'MarkerFaceColor',colr,...
     'MarkerEdgeColor',colr,'Color',colr,'MarkerSize',3);
    hold on
    end
    ylabel('Beta (wt.%)');title(names(i));
    ax.YLim=[beta_min-(beta_max-beta_min)*0.050,beta_max+(beta_max-beta_min)*0.050];
    ylim=ax.YLim;
    Ymin=min(beta);rng=range(ylim);
    if i==1 || i==3 || i==5 || i==7
        text(2.4,ylim(1)-rng*0.2,'Wavelength (μm)','FontSize',13,'color',[0.15 0.15 0.15]);
    else
        ax.YAxisLocation='left';
    end
    tx1=join(['c=',num2str(const_avg,'%.1f'),' wt.%']);  
    ylim=ax.YLim;
    text(0.5,min(ylim)+range(ylim)*1.08,tx1,'FontSize',12,'Color',[0.2 0.2 0.2]);%[0 0.4 0.85]
    ax.XLim=[wavelength(1)-0.05 wavelength(end)+0.05];
    ax.TickDir='in'; 
    ax.XMinorTick = 'on'; ax.YMinorTick = 'on';
    sgtitle('Average of 50');
    set(ax,'FontName','Helvetica-Narrow','FontSize',13,'linewidth',1.0);
    set(gcf, 'Color', 'w');
    ax.YLabel.Color='k';ax.TickLength = [0.013 0.035];
end
set(gcf,'Position',[0 0 1500 800]);
cd(result_path);
pause(3);
exportgraphics(f,'Beta_plot.jpg','Resolution',300)
savefig(f,'Beta_plot.fig')
%% plot the R-square of the model
bands_num=length(wavelength);
f1=figure();
for i=1:length(names)
    train_X=train_X_all(:,:,i);
    train_Y=train_Y_all(:,i);
    test_X=test_X_all(:,:,i);
    test_Y=test_Y_all(:,i);
    test_n=length(test_Y);
    train_n=length(train_Y);
    filename=join([names(i),'sele.mat'],'_');
    load(filename);
    indx=logical(sele_chrom_final);
    [m,n]=size(sele_chrom_final);
    Beta_all=zeros(m,n);
    const=zeros(m,1);
    b_ind=indx(1,:);
    fit_test=zeros(test_n,sele_num);
    fit_train=zeros(train_n,sele_num);
    for j=1:sele_num
        const(j)=sele_beta_final{j}(1);
        Beta_all(j,indx(j,:))=sele_beta_final{j}(2:end);
        b_ind=b_ind|indx(j,:);
        test_X_sele=test_X(indx(j,:),:)';
        train_X_sele=train_X(indx(j,:),:)';
        fit_test(:,j)=[ones(test_n,1)  test_X_sele]*sele_beta_final{j};
        fit_train(:,j)=[ones(train_n,1)  train_X_sele]*sele_beta_final{j};
    end
    fit_test(find(fit_test<0))=0; % set invalid values below 0 to 0
    fit_train(find(fit_train<0))=0; % set invalid values below 0 to 0
    Oxide_fit_test(:,i)=mean(fit_test,2);
    Oxide_fit_train(:,i)=mean(fit_train,2);
    %test
    tss=sum((test_Y-mean(test_Y)).^2);
    rss=sum((test_Y- Oxide_fit_test(:,i)).^2);
    Rsquare_test(i)=1-rss/tss;
    Rmsep_test(i)=sqrt(rss/test_n);
    st=std(test_Y);
    Rpd_test(i)=st/Rmsep_test(i);
    subplot(2,3,i);
    p1=plot(test_Y,Oxide_fit_test(:,i),'bo','LineWidth',1.2,'MarkerSize',4);
    p1.Color=[0.9 0.4 0];
    hold on
    %train
    tss=sum((train_Y-mean(train_Y)).^2);
    rss=sum((train_Y- Oxide_fit_train(:,i)).^2);
    Rsquare_train(i)=1-rss/tss;
    Rmsep_train(i)=sqrt(rss/test_n);
    st=std(test_Y);
    Rpd_train(i)=st/Rmsep_train(i);
    p2=plot(train_Y,Oxide_fit_train(:,i),'bo','LineWidth',1.2,'MarkerSize',4);
    p2.Color=[0 0.4 0.9];
    hold on
    axis tight
    axis equal
    ax=gca;x=ax.XLim;
    x_rang=x(2)-x(1);
    rang_coef=0.12;
    ax.XLim=[max(0,x(1)-x_rang*rang_coef),max(0,x(2)+x_rang*rang_coef)];
    ax.YLim=[max(0,x(1)-x_rang*rang_coef),max(0,x(2)+x_rang*rang_coef)];
    ax.XMinorTick='on';ax.YMinorTick='on';ax.XTick=ax.YTick;
    ax.TickLength = [0.03 0.04];
    x=ax.XLim;x_rang=x(2)-x(1);
    p3=plot(x,x,'--');
    p3.LineWidth=1.1;
    p3.Color=[0.5,0.5,0.5];
    ax=gca;
    cor=0.2;
    ax.XAxis.Color = [cor cor cor];ax.YAxis.Color = [cor cor cor];
    xlabel('Measured (wt.%)'); ylabel('Estimated (wt.%)');   
    txpo=[min(x)+x_rang*0.1,min(x)+x_rang*0.8];
    text(txpo(1),txpo(2),{'R^2 = ',num2str(Rsquare_test(i),'% .4f')},...
        'FontSize',11,'Color',[0.9 0.4 0],'FontName','Helvetica-Narrow','FontWeight','normal');
    txpo=[max(x)-x_rang*0.4,max(x)-x_rang*0.8];
    text(txpo(1),txpo(2),{'R^2 = ',num2str(Rsquare_train(i),'% .4f')},...
        'FontSize',11,'Color',[0 0.4 0.9],'FontName','Helvetica-Narrow','FontWeight','normal');
    title(names(i));
    hold off
    set(gca,'FontName','Helvetica-Narrow','FontSize',12,'linewidth',1)
    set(gcf, 'Color', 'w');
    axall(i)=ax;
    
end
% Calculate the R-square of all the data
Oixde_fit_test_combined=reshape(Oxide_fit_test,[],1);
test_Y_all_combined=reshape(test_Y_all,[],1);
tss=sum((test_Y_all_combined-mean(test_Y_all_combined)).^2);
rss=sum((test_Y_all_combined- Oixde_fit_test_combined).^2);
Rsquare_test_combined=1-rss/tss;
Rmsep_test_combined=sqrt(rss/length(Oixde_fit_test_combined));

legend('Validation set','Training set');
set(gcf,'Position',[50 50 900 550]);
for i=1:length(names) axall(i).XTick=axall(i).YTick; end
cd(result_path);
exportgraphics(f1,'Rsquare.jpg','Resolution',300)
savefig(f1,'Rsquare.fig')
%% display training and validation sets
f_sets=figure();
id=1:3:length(wavelength)-1;
for i=1:length(names)
    subplot(2,3,i)
    [train_X,train_Y,test_X,test_Y] = traintestsplit(Ref,Comp(:,i),3);
    plot_x=repmat(wavelength(id)',length(train_Y),1)';
    plot_y=repmat(train_Y,1,length(id))';
    plot_z=train_X((id),:);
    p1=plot3(plot_x,plot_y,plot_z,'Color',[0 0.5 0.8],'LineWidth',0.8);
    hold on
    p2=plot3(repmat(wavelength(id)',length(test_Y),1)',repmat(test_Y,1,length(id))',test_X((id),:),'Color',[0.900 0.2 0.1],'LineWidth',0.8); 
    view(3)
    ax=gca;
    ax.XLim=[0.4,2.4];ax.GridLineStyle = '-';ax.GridColor='k';
    axcolor='k';
    ax.XAxis.Color =axcolor;ax.YAxis.Color =axcolor;ax.ZAxis.Color =axcolor;
    ax.GridColor=axcolor;ax.GridAlpha=0.3;
    cont_scl=ax.YLim;
    t=title(names(i));
    t.Rotation=12;t.Position=[1.4 cont_scl(2) 0.62];
    ax.XLabel.String='(μm)';ax.XLabel.Color='k';
    ax.YLabel.String='(wt.%)';ax.YLabel.Color='k';
    ax.ZLabel.String='Ref.';ax.ZLabel.Color='k';
    ax.TickLength=[0.08 0.04];
    ax.YDir='reverse';
    ax.XMinorTick = 'on';ax.YMinorTick = 'on';ax.ZMinorTick = 'on';
    set(gca,'FontName','Helvetica-Narrow','FontSize',11.5,'FontWeight','bold','linewidth',1)
    set(gcf, 'Color', 'w');
    grid on
end
ax8=subplot(2,3,6);
legend(ax8,[p1(1),p2(1)],{'Train','Validation'}); legend('boxoff');
set(gca,'FontName','Helvetica-Narrow','FontSize',11.5,'linewidth',0.5)
axis padded
set(gcf, 'Color', 'w');
set(gcf,'Position',[40 40 1200 660]);
cd(result_path);
exportgraphics(f_sets,'Sets_disp.jpg','Resolution',300)
savefig(f_sets,'Sets_disp.fig')
%%
restoredefaultpath