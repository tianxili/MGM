function adj_norm = LeeHastie(p,q,n,y,z,lambda)
X = y;
Y = z+1;
M = length(lambda);
adj_norm = struct('zz', zeros(q,q,M), 'zy', zeros(q,p,M), 'yy', zeros(p,p,M));
for m = 1:M
    L=2*ones(q,1); 
    Ltot=sum(L);
    D=[];
    for j=1:q
        Dj=zeros(n,L(j));
        for i=1:n
            Dj(i,Y(i,j))=1;
        end
        D=[D Dj];
    end
    
    theta=zeros(Ltot,p); % cts-dis params
    beta=zeros(p,p); % negative of the precision matrix
    betad=ones(p,1); % diagonal of the precision matrix
    alpha1=zeros(p,1); % cts node potential param
    alpha2=zeros(Ltot,1); % dis node potential param
    phi=zeros(Ltot,Ltot); % dis edge potential params
    Lsum=[0;cumsum(L)];
    if m==1
        x=paramToVecv5(beta,betad,theta,phi,alpha1,alpha2,L,n,p,q);
    else
        x = x_old;
    end

    lam=lambda(m);
    smoothF= @(x)lhoodTfocsv5(x,D,X,Y,L,n,p,q);
    nonsmoothH=@(varargin) tfocsProxGroupv6(lam,L,n,p,q, varargin{:} ); % only returns value of nonsmooth
    opts.alg='N83';  opts.maxIts=800; opts.printEvery=100; opts.saveHist=true;
    opts.restart=-10^4;
    opts.tol=1e-6;
    [xopt out opts]=tfocs(smoothF, {}, nonsmoothH, x,opts);
    [beta betad theta phi alpha1 alpha2]= vecToParamv5(xopt,L,n,p,q);
    adj_norm_single = adj_norm_single_fromLee(beta, phi, theta, L);
    adj_norm.zz(:,:,m) = adj_norm_single.zz;
    adj_norm.zy(:,:,m) = adj_norm_single.zy;
    adj_norm.yy(:,:,m) = adj_norm_single.yy;
    x_old = xopt;
end


