function [train_X,train_Y,test_X,test_Y] = traintestsplit(X,Y,interval)
%   This function splits the datasets as training set and testing set
%   Here, I first sort the chemistry value in asending order,
%   so the training and testing sets span the whole dynamic range
%   Then select the testing samples every interval 
%   Inputs:
%           X-->Spectra, L bands * N samples
%           Y-->Chemistry, N samples * 1
%           interval-->sampling interval
%   Outputs:
%           train_X-->training spectra
%           train_Y-->training chemistry value
%           test_X-->testing spectra
%           test_Y-->testing chemistry value


    N = length(Y);
    if N == size(X,2)
        %   Sort Y
        [Y_sort,Ind] = sort(Y);
        X_sort = X(:,Ind);

        test_Ind = 1:interval:N;%1/interval·ÝµÄ²âÊÔ¼¯
        test_Ind=ceil(test_Ind);
        train_Ind = 1:N;
        train_Ind(test_Ind) = [];

        train_X = X_sort(:,train_Ind);
        test_X  = X_sort(:,test_Ind);
        train_Y  = Y_sort(train_Ind);
        test_Y  = Y_sort(test_Ind);
    end

end