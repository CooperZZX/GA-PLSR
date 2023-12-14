%本程序用于呈现GA算法选择的波段和每个波段的beta值、对最终含量平均贡献值
%b_ind是波段选择id
%beta是PLSR的Beta值
%spec是包含波长信息的光谱
%element是beta值对应的元素
function Beta_plot(b_ind,beta,spec,element)
[m,n]=size(spec);
if m>n%一般样本数小于波段数（如84个样本<389个波段）
    spec=spec';%一行是一条光谱
end
wl=spec(1,:);%波长
spec(1,:)=[];
[m,n]=size(spec);
beta_p=zeros(1,n);
beta_p(b_ind)=beta(2:end);
figure();
ax_beta(1)=subplot(3,1,2);
stem(wl(b_ind),beta_p(b_ind),'MarkerFaceColor',[0 0.4 0.85],...%[0 0.4470 0.7410]
    'MarkerEdgeColor',[0 0.4 0.85],'Color',[0 0.4 0.85],'MarkerSize',3);
ylabel('Beta');

ax_beta(2)=subplot(3,1,1);
plot(wl,spec);
hold on;
spec_mean=mean(spec);
ax_beta(2).YGrid='on';
% plot(wl,spec_mean,'b','LineWidth',1.5);
ylabel('Derivative');

ax_beta(3)=subplot(3,1,3);
contrib_mean=beta_p.*spec_mean;
for i=1:m
contrib=beta_p.*spec(i,:);
stem(wl(b_ind),contrib(b_ind),'MarkerSize',3);%,'MarkerFaceColor',[0 0.4 0.85],...%[0 0.4470 0.7410]
    %'MarkerEdgeColor',[0 0.4 0.85],'Color',[0 0.4 0.85],'MarkerSize',3);
    alpha(0.5);
hold on
end
%stem(wl(b_ind),contrib_mean(b_ind),'MarkerFaceColor',[0.8500 0.3250 0.0980],...%[0 0.4470 0.7410]
%   'MarkerEdgeColor',[0.8500 0.3250 0.0980],'Color',[0.8500 0.3250 0.0980],'MarkerSize',3);
tx1=join(['cons.=',num2str(beta(1),'%.1f'),' wt.%']);
ax=gca;
Lim=ax.YLim;
text(0.42,min(Lim)+range(Lim)*0.92,tx1,'FontSize',14,'Color',[0.2 0.2 0.2]);%[0 0.4 0.85]

ylabel('Cont. (wt.%)');xlabel('Wavelength (μm)');

sgtitle(element);

for i=1:3
    subplot(3,1,i);
    ax=gca;
    ax.XLim=[0.4 2.4];
    ax.TickDir='out';
    set(ax,'FontName','Helvetica-Narrow','FontSize',14,'linewidth',1.2);%,'Fontweight','bold'
    set(gcf, 'Color', 'w');
    cor=0.2;ax.XAxis.Color = [cor cor cor];ax.YAxis.Color = [cor cor cor];
end
% saveas(gca,strcat(element,'_Beta.fig'));

end

