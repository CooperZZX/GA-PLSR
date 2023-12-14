function bands=chrom2bands(wavelengh,chrome)
% 根据chrome的内容，提取实际对应的波段位置
% cgrome的组织形式个体数*基因长度，基因长度对应波段长度
% wavelength是波长向量，与染色体的长度相同
[p,q]=size(chrome);
bands=zeros(p,q);
for i=1:p  % 循环提取每一个体的波段组合
    for j=1:q  % 遍历所有波段
        if chrome(i,j)==1
            bands(i,j)=wavelengh(j);  % 存储新提取的光谱矩阵
        else
            bands(i,j)=0;
        end
    end
end