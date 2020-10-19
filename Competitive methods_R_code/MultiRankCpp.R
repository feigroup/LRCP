library(inline)
#X is the data matrix
#L is the maximum number of clusters


src = '
	using namespace Rcpp;
	int L = as<int>(L_);
	int n = as<int>(n_);
	NumericMatrix D(D_);
	NumericMatrix II(L,n);
	NumericMatrix H(L,n);
	
	for(int i=1;i<n;++i)//initialize objective function matrix
		II(0,i) = D(0,i);
	//use DP to generate the additional rows of the objective function matrix
	for(int p=1;p<L;++p){
		for(int j=p;j<n;++j){
			for(int h=p;h<=j;++h){
				if(II(p,j) < II(p-1,h-1) + D(h,j)){
					II(p,j) = II(p-1,h-1) + D(h,j);
					H(p,j) = h;
				}
			}
		}
	}
	return List::create(II,H);
'

dynamic_programming = cxxfunction(signature(L_='numeric',n_='numeric',D_='matrix'),src,plugin='Rcpp')

MultiRank = function(X,L){
	X_sort = X
	K = nrow(X)#dimension of observations
	n = ncol(X)#length of each time series
	for(i in 1:K)
		X_sort[i,] = sort(X_sort[i,])
	rank = X
	F = rank
	Cov = matrix(0,K,K)
	for(i in 1:K)#calculate ranks
		for(j in 1:n)#rank the i-th coordinate for j-th observation
			rank[i,j] = findInterval(rank[i,j],X_sort[i,])
	
	for(i in 1:K)#calculate Cov matrix
		for(j in 1:K)
			for(k in 1:n)
				Cov[i,j] = Cov[i,j] + (rank[i,k]/n-1/2)*(rank[j,k]/n-1/2)
	Cov = 4*Cov/n
	Cov = MASS::ginv(Cov)#generalized inverse
	
	DD = matrix(-Inf,n,n)
	
	for(i in 1:(n-1)){
		for(j in (i+1):n){
			Rbar = rep(0,K)
			for(k in 1:K)
				Rbar[k] = Rbar[k] + sum(rank[k,i:j])
			Rbar = Rbar/(j-i+1) - (n+1)/2#similar to average of rank of points i to j inclusive
			DD[i,j] = (t(Rbar)%*%Cov%*%Rbar)*(j-i+1)
		}
	}
	
	res = dynamic_programming(L,n,DD)
	
	II = res[[1]]
	H = res[[2]]
	
	#find number of change points
	k = 2:L
	gk = II[,n]
	gk = gk[-1]
	nk = length(k)
	Rsq = numeric(nk)
	for(p in 2:(nk-2)){
		f = as.factor(c(rep('l',p),rep('r',nk-p)))
		Rsq[p] = (summary(lm(gk ~ -1+f+k:f)))$r.sq
	}
	p.hat = which.max(Rsq) + 1
	#backtrack
	point = numeric()
	l = p.hat
	last = n
	while(l>1){
		point = c(H[l,last],point)
		last = H[l,last]
		l = l-1
	}
	point = c(1,point,n+1)
	member = numeric()
	for(i in 1:p.hat)
		member = c(member,rep(i,point[i+1]-point[i]))
	ret = NULL
	ret$clusters = member
	ret$cluster_number = p.hat
	ret$list = point
	return(ret)
}
