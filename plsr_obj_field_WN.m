function [fitn, process_par, beta]=plsr_obj_field_WN(train_X,train_Y,test_X,test_Y,CE4_D10_dat,chrome,other)
    
    [p,~]=size(chrome);
    beta=cell(1,p);
    fitn=zeros(p,1);
    process_par=zeros(p,10);
    for i=1:p  % ѭ������ΪȾɫ�����
%         temp1=[]; % ?????????????,???????????
%         temp2=[];
%         for j=1:q  % ??????
%             if chrome(i,j)==1   %chrome���ڴ洢��ѡ�Ĳ���λ��
%                 temp1=[temp1;train_X(j,:)];  %��Ϊ�µ�train_x
%                 temp2=[temp2;test_X(j,:)];   %��Ϊ�µ�test_x
%             end
%         end
        
        %��������Ĵ���֮���ٶ�������3%
        t1=find(chrome(i,:)>0);
        temp1=train_X(t1,:);
        temp2=test_X(t1,:); 
% train_X��test_X��Ҫһ��һ�����ף�train_Y��test_Y��Ҫһ��һ������
    [fitn(i),~,~,~,beta{1,i},process_par(i,:),~]=plsr_WN_n(temp1,train_Y,temp2,test_Y,chrome(i,:),CE4_D10_dat,other);%temp1��train_x������ѡ��һ����������PLSR������RMSEP��
    % plsr_WN_n���صĵ�һ��������rmesp����fitn��
    end
end