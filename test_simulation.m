function [beta_w,beta_rmse,y_rmse,y_test_rmse]=test_simulation(w,X,M,y,X_test,M_test,y_test,beta_true,beta0hat,betahat,n,n_test)

for i = 1:1:length(w)
    if i == 1
        beta_w = betahat{1,i}*w(i);
        beta0_w = beta0hat{1,i}*w(i);
    else
        beta_w = beta_w+betahat{1,i}*w(i);
        beta0_w = beta0_w+beta0hat{1,i}*w(i);
    end
end
beta_rmse = (double(beta_w)-beta_true).^2;
beta_rmse = mean(beta_rmse(:))^0.5;
y_pred=X*beta0_w + double(ttt(tensor(beta_w), M, 1:ndims(beta_w)));
y_test_pred=X_test*beta0_w + double(ttt(tensor(beta_w), M_test, 1:ndims(beta_w)));
y_rmse = ((y_pred-y)'*(y_pred-y)/n)^0.5;
y_test_rmse = ((y_test_pred-y_test)'*(y_test_pred-y_test)/n_test)^0.5;
end