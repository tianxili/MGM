addpath('../glmnet_matlab') %%% This is the path for glmnet, please chanage it accordingly.

%% This is a simple example using one of the simulation setting of the paper


%%%%%%%%% simulation setup %%%%%%%%%%%
n = 200; %sample size
p = 40;
q = 10;

a = 0.1;
b=0.2;
L = 20; % number of lambda to try

lambda = 0.5.^linspace(-1,17,20);
adj = zeros(p+q, p+q);
adj(6:15, 6:15) = 1;
adj(1:5, 1:5) = 1;
adj(16:20, 16:20) = 1;
adj = adj-diag(diag(adj));
rng('default');
seed = 10;
rng(seed);
parlist = pargen_unif(adj, p, q, a, b);

partmp = parlist;
[z y] = datagen(partmp, n);                           



%%%%%%%% FIT the Weighted Lasso%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
[fitlist_post fitlist]= sepreg_weight_flex(z,y,lambda, lambda, 'max',0.1);
adj_norm = edgenorm(fitlist_post);
adj_lambda = edgethres(adj_norm, 0);
[adj_sens adj_spec total_p total_n adj_sens_zz adj_spec_zz total_p_zz total_n_zz adj_sens_yy adj_spec_yy total_p_yy total_n_yy adj_sens_zy adj_spec_zy total_p_zy total_n_zy ] = adjcompare_category(adj_lambda, adj,p,q);


 
