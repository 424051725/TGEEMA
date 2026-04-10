function [y_train,y_test,X_train,X_test,M_train,M_test] = split_data(M,X,y,n,m)
% set 24-month MRI data as test data
y_train = y;
y_train(5.*(1:n),:)=[];
y_test = y(5.*(1:n),:);
X_train = X;
X_train(5.*(1:n),:)=[];
X_test = X(5.*(1:n),:);
index=1:n*5;
index=index(mod(index, 5) ~= 0);
M_train=M(:,:,:,index);
M_test = M(:,:,:,5.*(1:n));
end