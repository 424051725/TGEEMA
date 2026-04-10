beta = 1-imread(['pic_64\',num2str(pic_name),'.png']);
beta = double(beta);
[p1,p2] = size(beta);
p0 = 1;
b0 = ones(p0,1);
X = randn(n*m,p0);
M = tensor(randn(p1,p2,n*m)); 
mu = X*b0 + double(ttt(tensor(beta), M, 1:2));
if strcmp(workcorr,'equicorr') %  'equicorr', 'AR1', 'unstructured'
    S_inner=0.8.*ones(m)+0.2.*eye(m); 
elseif strcmp(workcorr,'AR1')
    S_inner = zeros(m);
    for i = 1:m
        S_inner(i,i) = 1;
        for j = 1:i-1
            S_inner(i,i-j) = 0.5^j;
            S_inner(i-j,i) = 0.5^j;
        end
    end
elseif strcmp(workcorr,'unstructured')
    S_inner = zeros(m);
    for i = 1:m
        for j = 1:m
            S_inner(i,j) = 1/(abs(i-j)+1);
        end
    end
elseif strcmp(workcorr,'indep')
    S_inner = eye(m);
else
    error('workcorr type not recogonized. equicorr, AR1, unstructured, indep accepted');
end

sigma=1;
S=sigma.*kron(eye(n),S_inner);
y=mvnrnd(mu,S)';
n_test = 100;
X_test = randn(n_test*m,p0);
M_test = tensor(randn(p1,p2,n_test*m)); 
mu_test = X_test*b0 + double(ttt(tensor(beta), M_test, 1:2));
S_test=sigma.*kron(eye(n_test),S_inner);
y_test=mvnrnd(mu_test,S_test)';

id=kron(1:n,ones(1,m))';
time=kron(ones(1,n),6:6:(6*m))';