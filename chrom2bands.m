function bands=chrom2bands(wavelengh,chrome)
% ����chrome�����ݣ���ȡʵ�ʶ�Ӧ�Ĳ���λ��
% cgrome����֯��ʽ������*���򳤶ȣ����򳤶ȶ�Ӧ���γ���
% wavelength�ǲ�����������Ⱦɫ��ĳ�����ͬ
[p,q]=size(chrome);
bands=zeros(p,q);
for i=1:p  % ѭ����ȡÿһ����Ĳ������
    for j=1:q  % �������в���
        if chrome(i,j)==1
            bands(i,j)=wavelengh(j);  % �洢����ȡ�Ĺ��׾���
        else
            bands(i,j)=0;
        end
    end
end