function [beta0hat,betahat,alphahat,devhat,AIC,BIC] =tensor_gee_ma(X,M,y,r,id,time,model,workcorr,lambda,penalty)
	n=length(unique(id));
	m=length(id)/n;
	d = ndims(M)-1;             
	p = size(M);     
	p0=size(X);
	p0=p0(2);
	% pre-allocate variables
	geestats = cell(1,d+1);
	dev_final = inf;   
	
	% loop for various intial points
	for rep=1:1		
		%initial value
		beta = ktensor(arrayfun(@(j) 1-2*rand(p(j),r), 1:d, 'UniformOutput',false));    		
		% main loop
		for iter=1:100
			% update coefficients for the regular covariates
			if (iter==1)
				penidx =[false(p0,1)];
				[beta0gee,alphagee,geestats{d+1}] = gee_sparsereg(id,time,X,y,model,workcorr,lambda,'penidx',penidx,'penalty',penalty);
				dev0=inf;
			else
				eta = Xj*beta{d}(:);
				penidx =[false((p0+1),1)];
				[betatmp,alphagee,geestats{d+1}] = gee_sparsereg(id,time,[X,eta],y,model,workcorr,lambda,'penidx',penidx,'penalty',penalty);
				beta0gee = betatmp(1:end-1);
				% stopping rule
				yhat= X*beta0gee + double(ttt(tensor(beta), M, 1:d));
		        devtmp=sum((y-yhat).^2);
				diffdev = devtmp-dev0;
				dev0 = devtmp;
				if (abs(diffdev)<(1e-4)*(abs(dev0)+1))
					break;
				end
				beta = arrange(beta*betatmp(end));
				for j=1:d
					beta.U{j} = bsxfun(@times,beta.U{j},(beta.lambda').^(1/d));
				end
				beta.lambda = ones(r,1);            
			end
			% cyclic update of the array coefficients
			eta0 = X*beta0gee;
			for j=1:d
				if j==1
					cumkr = ones(1,r);
				end
				if j==d
					Xj = reshape(double(tenmat(M,[d+1,j]))*cumkr,n*m,p(j)*r);
				else
					Xj = reshape(double(tenmat(M,[d+1,j]))*khatrirao([beta.U(d:-1:j+1),cumkr]),n*m,p(j)*r);
                end

                penidx =[true(size(Xj,2),1);false];

				[betatmp,alphagee,geestats{j}] = gee_sparsereg(id,time,[Xj,eta0],y,'normal',workcorr,lambda,'penidx',penidx,'penalty',penalty);
				beta{j} = reshape(betatmp(1:end-1),p(j),r);
				eta0 = eta0*betatmp(end);
				cumkr = khatrirao(beta{j},cumkr);
			end
		end
		yhat= X*beta0gee + double(ttt(tensor(beta), M, 1:d));
		dev0=sum((y-yhat).^2);
		% record if it has a smaller deviance
		if dev0<dev_final
			beta0hat = beta0gee;
			betahat = beta;
			alphahat=alphagee;
			devhat=dev0;
			dev_final = dev0;
		end
		
	end
%output BIC
	BIC = dev_final + log(n)*(r*(p(1)+p(2)-r)+p0);
    AIC = dev_final + 2*(r*(p(1)+p(2)-r)+p0);
end