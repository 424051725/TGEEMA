function [w_hat,w_equi,w_ar1,w_ind]=model_averaging(X,M,y,max_rank,fold,id,time,link)
n=length(y);
n0 = length(y)/fold;  
theta = zeros(n,max_rank*3);
temp = zeros(n-n0,max_rank*3);
Corr={'equicorr','AR1','indep'};
for i = 1:1:max_rank
    for corr_num=1:1:3
        corr=Corr(corr_num);
        for j = 1:1:fold
            y_j = y;
            y_j((j-1)*n0+1:min(j*n0,n)) = [];
            X_j = X;
            X_j((j-1)*n0+1:min(j*n0,n),:) = [];
            M_j = double(M);
            %id_j = kron(1:n0/length(unique(time)),ones(1,length(unique(time))))';
            id_j = id(1:(n-n0));
            time_j = time;
            time_j((j-1)*n0+1:min(j*n0,n),:) = [];
            indices = repmat({':'}, 1, ndims(M_j)-1);
            M_j(indices{:},(j-1)*n0+1:min(j*n0,n)) = [];
            M_j = tensor(M_j);
            [beta0hat,betahat,~,~,~,~] =tensor_gee_ma(X_j,M_j,y_j, i,id_j,time_j,link,corr,0,'SCAD');
            M_rest = M(indices{:},(j-1)*n0+1:min(j*n0,n));
            X_rest = X((j-1)*n0+1:min(j*n0,n),:);
            theta((j-1)*n0+1:min(j*n0,n),(i-1)*3+corr_num) = X_rest*beta0hat + double(ttt(tensor(betahat), M_rest, 1:ndims(betahat)));
            %temp(:,i) = X_j*beta0hat + double(ttt(tensor(betahat), M_j, 1:ndims(betahat)));
        end
    end
end

w = ones(max_rank*3,1)/(max_rank*3);
Aeq = ones(1,max_rank*3);
beq = 1;
[w_hat,~] = fmincon(@(x) obj_func(x,theta,y),w,[],[],Aeq,beq,zeros(1,max_rank*3),ones(1,max_rank*3));

w = ones(max_rank,1)/(max_rank);
Aeq = ones(1,max_rank);
beq = 1;

theta_temp = theta(:,[1:max_rank]*3-2);
[w_equi,~] = fmincon(@(x) obj_func(x,theta_temp,y),w,[],[],Aeq,beq,zeros(1,max_rank),ones(1,max_rank));
w_equi = kron(w_equi,[1,0,0]');
theta_temp = theta(:,[1:max_rank]*3-1);
[w_ar1,~] = fmincon(@(x) obj_func(x,theta_temp,y),w,[],[],Aeq,beq,zeros(1,max_rank),ones(1,max_rank));
w_ar1 = kron(w_ar1,[0,1,0]');
theta_temp = theta(:,[1:max_rank]*3);
[w_ind,~] = fmincon(@(x) obj_func(x,theta_temp,y),w,[],[],Aeq,beq,zeros(1,max_rank),ones(1,max_rank));
w_ind = kron(w_ind,[0,0,1]');
end