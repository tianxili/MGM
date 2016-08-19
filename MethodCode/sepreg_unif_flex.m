
function [fitlist_post fitlist a0_lr a0_ls beta_lr beta_ls sigma2] = sepreg_unif_flex(z,y,tune1,tune2,method,kappa)
% Input
% [z y]: data
% p: dim of y
% q: dim of z
% [tune1 tune2]: tuning parameter vector for logistic and linear
% regressions. 
% method: 'min' or 'max' postprocessing the fitlist
% Output:
% fitlist: a structure of fitted parameters path, the notations are the
% same as the paper. 
p = size(y,2);
q = size(z,2);
q1 = q*(q-1)/2;
L = length(tune1);

fitlist = struct('lambda_j', zeros(q,L), 'lambda_jk', zeros(q1,2,L), 'ita_0', zeros(p,L), 'ita_j', zeros(p,q,2,L),... 
                 'phi_0', zeros(p,p,2,L), 'phi_j', zeros(p,p,q,3,L));

%% fit models for z_j's
sprintf('fitting logistic regression')
for j = 1:q
    j
    sprintf('logistic regression %d',j);
    ytmp = (z(:,j)>0)+1;
    xtmp = design(z,y,'z',j);
    % Options for glmnet
    options = glmnetSet;
    options.maxit = 3000; %%% avoid long iterations without convergence for small lambda setting.
    options.lambda = tune1;
    options.standardize = false;
    options.penalty_factor = penalty_fac_unif_flex('z',p,q,kappa);
    if (size(xtmp,2)>size(xtmp,1))
    options.type = 'naive';
    end
    % Fit the entire solution path
    fit= glmnet(xtmp, ytmp, 'binomial', options);
    if size(fit.beta,2)<L;
        beta_tmp = [fit.beta repmat(fit.beta(:,end),1,L-size(fit.beta,2))];
        a0_tmp = [fit.a0  repmat(fit.a0(end),1,L-length(fit.a0))];
    else 
        beta_tmp = fit.beta;
        a0_tmp = fit.a0;
    end
    % Store the parameters
        fitlist = storage(fitlist, a0_tmp', beta_tmp, j);
    clear xtmp ytmp fit;
end

%% fit models for y_j's
sprintf('fitting linear regression')
for j = 1:p
    j
    sprintf('least square regression %d',j);
    ytmp = y(:,j);
    xtmp = design(z,y,'y',j);
    % Options for glmnet
    options = glmnetSet;
    options.maxit = 3000;
    options.lambda = tune2;
    options.standardize = false;
    options.penalty_factor = penalty_fac_unif_flex('y',p,q,kappa);
    if (size(xtmp,2)>size(xtmp,1))
    options.type = 'naive';
    end
     % Fit the entire solution path
    fit= glmnet(xtmp, ytmp, 'gaussian', options);
    sigma2_tmp = (1-fit.dev)*fit.nulldev;
    if size(fit.beta,2)<L;
        beta_tmp = [fit.beta repmat(fit.beta(:,end),1,L-size(fit.beta,2))];
        a0_tmp = [fit.a0; repmat(fit.a0(end),L-length(fit.a0),1)]';
        sigma2 = [sigma2_tmp; repmat(sigma2_tmp(end),L-length(sigma2_tmp),1)];
    else 
        beta_tmp = fit.beta;
        a0_tmp = fit.a0;
        sigma2 = sigma2_tmp;
    end
    % Store the parameters
    fitlist = storage(fitlist, a0_tmp, beta_tmp, j, sigma2);
    clear ytmp xtmp sigma2_tmp fit;
end
    fitlist_post = fitpost(fitlist, method);    
end



    
    
                 
