%% data reading
addpath('glmnet_matlab_MAC')

datname1 = 'Datasets/CAL500/CAL500_new.txt';

data = importdata(datname1);

p = 68;
%q = Q(i);
q = 118;


n = size(data.data,1);
z = data.data(:,1:q);
y = data.data(:,q+1:end);
ind = [1 3:5 18 20:22 35 37:39 52 54:56];
y = y(:,ind);

y = zscore(y);
size(z)
size(y)

%% data analysis

lambda =0.5.^linspace(0,5,40); 

kappa = 1.6:-0.3:0.1;

K = 100; % number of bagging sampling procedures for stability selection.

p = size(y,2);
q = size(z,2);

result = zeros(p+q,p+q,length(lambda),length(kappa));
n = size(z,1);

for k = 1:K
    k
    rand('state', k+500);
    n_sub = floor(n/2);
    ind = randsample(n,n_sub);
    z_tmp = z(ind,:);
    while(sum(sum(z_tmp,1)==n_sub)>0 || sum(sum(z_tmp,1)==0)>0) %% check if binary data are all 1 or 0
        sprintf('resample')
        ind = randsample(n,n_sub);
        z_tmp = z(ind,:);
    end
    y_tmp = y(ind,:);
    z_tmp = z(ind,:);
    for i = 1:length(kappa)
        i
	      [fitlist_post, fitlist]= sepreg_weight_flex(z_tmp,y_tmp,lambda, lambda, 'max',kappa(i));
        adj_norm = edgenorm(fitlist_post);
        adj_lambda = edgethres(adj_norm, 0);
        result(:,:,:,i) = result(:,:,:,i) + adj_lambda;
    end
end

    max_result = zeros(p+q,p+q);
for i = 1:(p+q)
	  for j = 1:(p+q)
		    max_result(i,j) = max(max(result(i,j,:,:)));
          end
end

	  final_adj = zeros(p+q,p+q);
final_adj(max_result>99) = 1;
sum(sum(final_adj(1:q,(q+1):(q+p))))

name = 'CAL500_stability_result2.txt'
dlmwrite(name,max_result);

