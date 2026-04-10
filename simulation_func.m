%% calculate tgma weight
[w_tgma,w_equi,w_ar1,w_ind] = model_averaging(X,M,y,max_rank,fold,id,time,link);
%% calculate each model's parameter
beta0hat=cell(1,max_rank*3);
betahat=cell(1,max_rank*3);
alphahat=cell(1,max_rank*3);
devhat=zeros(1,max_rank*3);
AIC=zeros(1,max_rank*3);
BIC=zeros(1,max_rank*3);
Corr = {'equicorr','AR1','indep'};
for i=1:1:max_rank
    for j=1:1:3
        corrlation=Corr(j);
        [beta0hat{1,(i-1)*3+j},betahat{1,(i-1)*3+j},alphahat{1,(i-1)*3+j},devhat(1,(i-1)*3+j),AIC(1,(i-1)*3+j),BIC(1,(i-1)*3+j)] =tensor_gee_ma(X,M,y,i,id,time,link,corrlation,0,'scad');
    end
end

%% calculate AIC and BIC
[min_xic,where] = min(AIC);
AIC_new =AIC-min_xic ;
w_aic = zeros(1,max_rank*3);
w_aic(1,where) = 1;
[min_xic,where] = min(BIC);
BIC_new =BIC-min_xic ;
w_bic = zeros(1,max_rank*3);
w_bic(1,where) = 1;
w_saic = exp(-2*AIC_new)/sum(exp(-2*AIC_new));
w_sbic = exp(-2*BIC_new)/sum(exp(-2*BIC_new));
w_equal = ones(1,max_rank*3)/(max_rank*3);
clear min_xic;
clear where;

%% calculate each model's performance on test data
W=[w_aic;w_bic;w_saic;w_sbic;w_equal;w_tgma.';w_equi.';w_ar1.';w_ind.';eye(max_rank*3)];
for i=1:length(W)
    w_i=W(i,:);
    [beta_w{1,i},beta_rmse(i),y_rmse(i),y_test_rmse(i)]=test_simulation(w_i,X,M,y,X_test,M_test,y_test,beta,beta0hat,betahat,n,n_test);
end

%% save the result
model_name = ['AIC';'BIC';'SAIC';'SBIC';'EQMA';'TGMA';'EQUI';'AR1';'IND';cellstr(num2str((1:max_rank*3)', '%d'))]';
r = num2cell([beta_rmse', y_rmse', y_test_rmse']);
if ~exist(['2d_simulation\','sigma=',num2str(sigma),'\',num2str(n),'\',num2str(repeat)],'dir')
mkdir(['2d_simulation\','sigma=',num2str(sigma),'\',num2str(n),'\',num2str(repeat)]);
end
xlswrite(['2d_simulation\','sigma=',num2str(sigma),'\',num2str(n),'\',num2str(repeat),'\',num2str(pic_name),'w_trma.xlsx'],w_tgma);
xlswrite(['2d_simulation\','sigma=',num2str(sigma),'\',num2str(n),'\',num2str(repeat),'\',num2str(pic_name),'result.xlsx'],vertcat({'model','BETA_RMSE','Y_RMSE','Y_TEST_RMSE'},[model_name',r]));
if repeat <=1
    for i = 1:1:length(model_name)    
        figure('visible','off');  
        imagesc(double(-beta_w{1,i}));
        colormap(gray);
        axis equal;
        axis tight;
        set(gca,'FontSize',25);
        saveas(gca, ['2d_simulation\','sigma=',num2str(sigma),'\',num2str(n),'\',num2str(repeat),'\',num2str(pic_name),'_',char(model_name(i))],'png')
    end
end
