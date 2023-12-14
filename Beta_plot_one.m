%% 本程序的主要功能是展示某一次特定运行时的最佳适应度的染色体的波段选择和beta值
%本程序用于呈现GA算法选择的波段和每个波段的beta值、对最终含量平均贡献值
%id是所要展示的第id次训练的母代个体，id>1，因为第一次训练的母代id未保存
%chrom_gen是迭代过程中的所有染色体chrom变化
%b_ind是波段选择id
%beta是PLSR的Beta值
%spec是包含波长信息的光谱
%element是beta值对应的元素
function Beta_plot_one(id,chrom_gen,b_ind,beta,spec,element)
beta_p=zeros(1,n);
beta_p(b_ind)=beta(2:end);
figure();
ax_beta(1)=subplot(3,1,2);
stem(wl(b_ind),beta_p(b_ind),'MarkerFaceColor',[0 0.4 0.85],...%[0 0.4470 0.7410]
    'MarkerEdgeColor',[0 0.4 0.85],'Color',[0 0.4 0.85],'MarkerSize',3);
ylabel('Beta');
-
end