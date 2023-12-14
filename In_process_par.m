%本程序用于实现对过程中的变量进行检查
%process_par_gen %训练集参数2（R2,Rmsep），验证集参数2，全部数据的参数2，主成分个数, CE4D10
%std，用于计算fitness的rmsep，fitness 20×10×1000
%beta_gen迭代过程中的beta参数，1000×20 cell
%chrom_gen迭代过程中的染色体，20×390×1000
%fitn_gen是迭代过程中的适应度，20×1000
%element对应的元素
%gen_i是标记的迭代次数
function In_process_par(process_par_gen,beta_gen,chrom_gen,WL,element,gen_i)
f=figure();
gen=length(beta_gen(:,1));
par_best=zeros(gen,length(process_par_gen(1,:,1)));
chrom_best=zeros(gen,length(chrom_gen(1,:,1)));
fitn_gen=permute(process_par_gen(:,10,:),[1,3,2]);%交换矩阵维度，变为二维矩阵
[fit_min,Ind]=min(fitn_gen,[],1);%fit_min是包含每一行的最小值的列向量
fit_mean=mean(fitn_gen,1);
fit_std=std(fitn_gen,1,1)*5;%N平均，第1维度
beta_best=cell(gen,1);
for i=1:gen
    beta_best(i)=beta_gen(i,Ind(i));
    par_best(i,:)=process_par_gen(Ind(i),:,i);
    chrom_best(i,:)=chrom_gen(Ind(i),:,i);
%     rmsep_best(i)=objv_gen(i,Ind(i));
%     std_best(i)=qual_gen(i,Ind(i));
end
subplot(4,1,1)
plot(par_best(:,9),'LineWidth',1);%rmsep_best
hold on
plot(par_best(:,8),'LineWidth',1);%std_best
plot(fit_min,'--');
xline(gen_i);
legend('rmsep best','std best','fitness min');
subplot(4,1,2)
plot(par_best(:,1));%train R2
hold on
plot(par_best(:,3));%test R2
xline(gen_i,'-','training stop');
legend('R^2 train','R^2 test','Location','east');
subplot(4,1,3)
p1=plot(fit_min);
hold on
p2=plot(fit_mean);
ax=gca;
y_lim=ax.YLim;
offset_y=mean(y_lim)-range(y_lim)*0.2-mean(fit_std);
plot(fit_std+offset_y);
plot(fit_mean,'Color',p2.Color);
plot(fit_min,'Color',p1.Color);
xline(gen_i);
yline(mean(y_lim)-range(y_lim)*0.2);
legend('fitness min','fitness mean','fitness std*5 (offseted)');
subplot(4,1,4)%绘制gen_i次迭代的beta结果
WL_id=find(chrom_best(gen_i,:));
st1=stem(WL(WL_id),beta_best{gen_i}(2:end),'filled','MarkerSize',3);
legend(st1,'beta of spectral deri');
hold on

if WL_id(end)==length(WL)
    st2=stem(WL(end),beta_best{gen_i}(end),'MarkerSize',3,'MarkerFaceColor','red','MarkerEdgeColor','red');
    legend([st1,st2],{'betas of spectral deri.','beta of spectral ref.'});
end
ax=gca;
ax.XMinorTick='on';
% ax.XLim=[0.4 2.45];
axis padded
tx1=join(['c=',num2str(beta_best{gen_i}(1),'%.2f'),' wt.%']);  
set(gcf, 'Color', 'w','Position',[200,50,800,1000]);
ylim=ax.YLim;
text(0.5,min(ylim)+range(ylim)*1.05,tx1,'Color',[0.2 0.2 0.2]);%[0 0.4 0.85] ,'FontSize',5
sgtitle(element);
saveas(f,strcat(element,'_in_pro.fig'));
exportgraphics(f,strcat(element,'_in_pro.jpg'),'Resolution',300)
end
