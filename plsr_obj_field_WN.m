function [fitn, process_par, beta]=plsr_obj_field_WN(train_X,train_Y,test_X,test_Y,CE4_D10_dat,chrome,other)
    
    [p,~]=size(chrome);
    beta=cell(1,p);
    fitn=zeros(p,1);
    process_par=zeros(p,10);
    for i=1:p  % 循环次数为染色体个数
%         temp1=[]; % ?????????????,???????????
%         temp2=[];
%         for j=1:q  % ??????
%             if chrome(i,j)==1   %chrome用于存储挑选的波段位置
%                 temp1=[temp1;train_X(j,:)];  %作为新的train_x
%                 temp2=[temp2;test_X(j,:)];   %作为新的test_x
%             end
%         end
        
        %改用下面的代码之后，速度提升了3%
        t1=find(chrome(i,:)>0);
        temp1=train_X(t1,:);
        temp2=test_X(t1,:); 
% train_X和test_X需要一列一条光谱，train_Y和test_Y需要一行一个样本
    [fitn(i),~,~,~,beta{1,i},process_par(i,:),~]=plsr_WN_n(temp1,train_Y,temp2,test_Y,chrome(i,:),CE4_D10_dat,other);%temp1是train_x，对挑选的一组样本进行PLSR，计算RMSEP。
    % plsr_WN_n返回的第一个参数是rmesp，即fitn。
    end
end