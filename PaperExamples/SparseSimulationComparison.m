cd UGM
addpath(genpath(pwd))

cd ..

addpath('TFOCS-1.3.1')

addpath('glmnet_matlab_MAC')

n = 100; % sample size
p = 90;
q = 10;
K = 100;

A = [1 0.4 4];
B = [2 0.4 4];
for ii = 1:1
    for jj = 1:1 
   a = A(ii); % absolute value of parameters
    b = B(jj); % value of off-diagonal entries of matrix
    D = [10 10 10 19]; %pth of lambda log_0.5
    L = [60 60 60 60];  % length of lambda
    kappa_seq = 0.1:0.1:1;
    rng('default');
    seed = 50;
    rng(seed);
    lambda1 = 0.5.^linspace(0,D(1),L(1));
    lambda2 = 0.5.^linspace(0,D(2),L(2));
    lambda3 = 0.5.^linspace(-4,D(3),L(3));
    lambda4 = 0.5.^linspace(-2,D(4),L(4));
    lambda_log = repmat(lambda3,[1, L(4)]);
    lambda_reg = reshape((repmat(lambda4,[L(3),1])),[1,L(3)*L(4)]);
%% Model fitting
      valid = 0;
      while(valid==0)
          adj = adjgen(n, 125, 6);
          if sum(sum(adj(1:q,1:q)))>6 && max(sum(adj))>6-1
              valid=1;
          end
      end


    parlist = pargen_unif(adj, p, q, a, b);
    model = {'both', 'main'};
Lee_X = zeros(n,p,2);
Lee_Y = zeros(n,q,2);
for i = 1:1
         partmp = parlist;
%          if (i==2)
%          partmp.ita_j = 0*partmp.ita_j;
%          partmp.phi_0 = diag(diag(partmp.phi_0));
%          else
         %if (i==2)
         partmp.phi_j = 0*partmp.phi_j;
         %end
         
         adj_par_count = edgecount(partmp);
         reg_par_count = reg_count(partmp);
         count = sum(partmp.lambda_jk~=0)+sum(sum(partmp.ita_j~=0))+...
                  sum(sum(partmp.phi_0~=0))-p+sum(sum(sum(partmp.phi_j~=0)))/2;
         idv2_sens = zeros(L(2),K);
         idv2_spec = zeros(L(2),K);
         idv3_sens = zeros(L(3),K);
         idv3_spec = zeros(L(3),K);
         idv4_sens = zeros(L(4),K);
         idv4_spec = zeros(L(4),K);
         idv5_sens = zeros(L(4),K);
         idv5_spec = zeros(L(4),K);
         adj1_sens = zeros(L(1),K);
         adj1_spec = zeros(L(1),K);
         adj2_sens = zeros(L(2),K);
         adj2_spec = zeros(L(2),K);
         adj3_sens = zeros(L(3),K);
         adj3_spec = zeros(L(3),K);
         adj4_sens = zeros(L(4),K);
         adj4_spec = zeros(L(4),K);         
         adj5_sens = zeros(L(4),K);
         adj5_spec = zeros(L(4),K);
         %% metrics for zz
         adj1_sens_zz = zeros(L(1),K);
         adj1_spec_zz = zeros(L(1),K);
         adj2_sens_zz = zeros(L(2),K);
         adj2_spec_zz = zeros(L(2),K);
         adj3_sens_zz = zeros(L(3),K);
         adj3_spec_zz = zeros(L(3),K);
         adj4_sens_zz = zeros(L(4),K);
         adj4_spec_zz = zeros(L(4),K);         
         adj5_sens_zz = zeros(L(4),K);
         adj5_spec_zz = zeros(L(4),K);
         %% metrics for zy
         adj1_sens_zy = zeros(L(1),K);
         adj1_spec_zy = zeros(L(1),K);
         adj2_sens_zy = zeros(L(2),K);
         adj2_spec_zy = zeros(L(2),K);
         adj3_sens_zy = zeros(L(3),K);
         adj3_spec_zy = zeros(L(3),K);
         adj4_sens_zy = zeros(L(4),K);
         adj4_spec_zy = zeros(L(4),K);         
         adj5_sens_zy = zeros(L(4),K);
         adj5_spec_zy = zeros(L(4),K);
         %% metrics for yy
         adj1_sens_yy = zeros(L(1),K);
         adj1_spec_yy = zeros(L(1),K);
         adj2_sens_yy = zeros(L(2),K);
         adj2_spec_yy = zeros(L(2),K);
         adj3_sens_yy = zeros(L(3),K);
         adj3_spec_yy = zeros(L(3),K);
         adj4_sens_yy = zeros(L(4),K);
         adj4_spec_yy = zeros(L(4),K);         
         adj5_sens_yy = zeros(L(4),K);
         adj5_spec_yy = zeros(L(4),K);

         for k = 1:K
              k
              [z y] = datagen(partmp, n);                           
              %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
              %%%%%%%% FIT method of regular lasso  %%%%%%%%%%%%%%%%%%

                  [fitlist_post_reg fitlist_reg]= sepreg_unif_flex(z,y,lambda1, lambda1, 'max',kappa_seq(1));
[idv1_sens(:,k) idv1_spec(:,k) idv_tp idv_tn] = fitcompare(fitlist_post_reg, parlist);
                  adj_norm = edgenorm(fitlist_post_reg);
                  adj_lambda = edgethres(adj_norm, 0);
                  
              [adj1_sens(:,k) adj1_spec(:,k) total_p total_n adj1_sens_zz(:,k) adj1_spec_zz(:,k) total_p_zz total_n_zz adj1_sens_yy(:,k) adj1_spec_yy(:,k) total_p_yy total_n_yy adj1_sens_zy(:,k) adj1_spec_zy(:,k) total_p_zy total_n_zy ] = adjcompare_category(adj_lambda, adj,p,q);

              
              %%%%%%%% FIT the Weighted Lasso%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
              
                  [fitlist_post fitlist]= sepreg_weight_flex(z,y,lambda3, lambda3, 'max',kappa_seq(1));
[idv3_sens(:,k) idv3_spec(:,k) idv_tp idv_tn] = fitcompare(fitlist_post, parlist);
                  adj_norm = edgenorm(fitlist_post);
                  adj_lambda = edgethres(adj_norm, 0);
                  
              [adj3_sens(:,k) adj3_spec(:,k) total_p total_n adj3_sens_zz(:,k) adj3_spec_zz(:,k) total_p_zz total_n_zz adj3_sens_yy(:,k) adj3_spec_yy(:,k) total_p_yy total_n_yy adj3_sens_zy(:,k) adj3_spec_zy(:,k) total_p_zy total_n_zy ] = adjcompare_category(adj_lambda, adj,p,q);
              
               %%%%%%%% FIT method of Lee & Hastie  %%%%%%%%%%%%%%%%%%
              fitlist_Lee = sepreg_Lee(z,y,lambda4, lambda4, 'max');
              [idv4_sens(:,k) idv4_spec(:,k) idv_tp idv_tn] = fitcompare(fitlist_Lee, parlist);
              adj_norm_Lee = edgenorm(fitlist_Lee);
              %adj_norm = LeeHastie(p,q,n,y,z,lambda4);
              adj_lambda_lee = edgethres(adj_norm_Lee, 0);
              %[adj4_sens(:,k) adj4_spec(:,k) adj_tp adj_tn] = adjcompare(adj_lambda_lee, adj);
              [adj4_sens(:,k) adj4_spec(:,k) total_p total_n adj4_sens_zz(:,k) adj4_spec_zz(:,k) total_p_zz total_n_zz adj4_sens_yy(:,k) adj4_spec_yy(:,k) total_p_yy total_n_yy adj4_sens_zy(:,k) adj4_spec_zy(:,k) total_p_zy total_n_zy ] = adjcompare_category(adj_lambda_lee, adj,p,q);
 %%%%%%%% FIT calibrated version of Lee & Hastie  %%%%%%%%%%%%%%%%%%
              adj_norm = LeeHastie(p,q,n,y,z,lambda4);
              adj_lambda = edgethres(adj_norm, 0);
              [adj5_sens(:,k) adj5_spec(:,k) adj_tp adj_tn] = adjcompare(adj_lambda, adj);
              [adj5_sens(:,k) adj5_spec(:,k) total_p total_n adj5_sens_zz(:,k) adj5_spec_zz(:,k) total_p_zz total_n_zz adj5_sens_yy(:,k) adj5_spec_yy(:,k) total_p_yy total_n_yy adj5_sens_zy(:,k) adj5_spec_zy(:,k) total_p_zy total_n_zy ] = adjcompare_category(adj_lambda, adj,p,q);
 
              %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
         end
              
              filename = sprintf('Simulation/Sparse_comparestudy_seeda=%.1f_b=%.1f%d_%s.mat',a,b,seed, char(model(i)));
              
              save(filename,'n','p','q','D','L','adj', 'parlist','fitlist_reg','fitlist_Lee',...
                            'adj_par_count', 'reg_par_count', 'count', 'adj1_sens', ...
                            'adj3_sens', 'adj4_sens','adj5_sens','adj1_spec', 'adj3_spec','adj4_spec','adj5_spec', ...
                             'idv3_spec',  'idv3_sens', 'idv4_spec', 'idv5_spec', 'idv4_sens', 'idv5_sens',...
                             'adj1_sens_zz','adj3_sens_zz', 'adj4_sens_zz','adj5_sens_zz',...
                         'adj1_sens_zy','adj3_sens_zy', 'adj4_sens_zy','adj5_sens_zy',....
                         'adj1_sens_yy','adj3_sens_yy', 'adj4_sens_yy','adj5_sens_yy',...
                             'adj1_spec_zz','adj3_spec_zz', 'adj4_spec_zz','adj5_spec_zz',...
                         'adj1_spec_zy','adj3_spec_zy', 'adj4_spec_zy','adj5_spec_zy',....
                         'adj1_spec_yy','adj3_spec_yy', 'adj4_spec_yy','adj5_spec_yy');
end
                  
         f = 1/5;
         g = figure;
         screen_size = get(0,'ScreenSize');
         set(g, 'Position', [0 0 1.2*screen_size(4) 0.5*screen_size(4)] );
         tt = {'Model I: main & interaction effects' 'Model II: main effects only'};
         for i = 1:1
              filename = sprintf('Simulation/Sparse_comparestudy_seeda=%.1f_b=%.1f%d_%s.mat',a,b,seed, char(model(i)));
              load(filename);
              subplot(2,2,1);
              vec_adj1_sens = adj1_sens(:);
              [vec_adj1_spec id1] = sort(adj1_spec(:));
              vec_adj3_sens = adj3_sens(:);
              [vec_adj3_spec id3] = sort(adj3_spec(:));
              vec_adj4_sens = adj4_sens(:);
              [vec_adj4_spec id4] = sort(adj4_spec(:));
              vec_adj5_sens = adj5_sens(:);
              [vec_adj5_spec id5] = sort(adj5_spec(:));
              y1 = smooth(1-vec_adj1_spec, vec_adj1_sens(id1), f,'lowess');
              y3 = smooth(1-vec_adj3_spec, vec_adj3_sens(id3), f,'lowess'); 
              y4 = smooth(1-vec_adj4_spec, vec_adj4_sens(id4), f,'lowess');
              y5 = smooth(1-vec_adj5_spec, vec_adj5_sens(id5), f,'lowess'); 
plot(1-vec_adj3_spec, y3,'g-',1-vec_adj4_spec, y4,'m--',1-vec_adj5_spec, y5, 'r-.',1-vec_adj1_spec, y1, 'b:','LineWidth', 2);
              axis([0 1 0 1]);
xlabel('Overall FPR','FontSize',20);
ylabel('Overall TPR','FontSize',20);
              % plot for zz edges
              subplot(2,2,2);
              vec_adj1_sens_zz = adj1_sens_zz(:);
              [vec_adj1_spec_zz id1] = sort(adj1_spec_zz(:));
              vec_adj3_sens_zz = adj3_sens_zz(:);
              [vec_adj3_spec_zz id3] = sort(adj3_spec_zz(:));
              vec_adj4_sens_zz = adj4_sens_zz(:);
              [vec_adj4_spec_zz id4] = sort(adj4_spec_zz(:));
              vec_adj5_sens_zz = adj5_sens_zz(:);
              [vec_adj5_spec_zz id5] = sort(adj5_spec_zz(:));
              y1 = smooth(1-vec_adj1_spec_zz, vec_adj1_sens_zz(id1), f,'lowess');
              y3 = smooth(1-vec_adj3_spec_zz, vec_adj3_sens_zz(id3), f,'lowess'); 
              y4 = smooth(1-vec_adj4_spec_zz, vec_adj4_sens_zz(id4), f,'lowess');
              y5 = smooth(1-vec_adj5_spec_zz, vec_adj5_sens_zz(id5), f,'lowess'); 
              plot(1-vec_adj3_spec_zz, y3,'g-',1-vec_adj4_spec_zz, y4,'m--',1-vec_adj5_spec_zz, y5,'r-.',1-vec_adj1_spec_zz, y1, 'b:','LineWidth', 2);
              axis([0 1 0 1]);
xlabel('ZZ-FPR','FontSize',20);
ylabel('ZZ-TPR','FontSize',20);
              % yy plot
             subplot(2,2,3);
              vec_adj1_sens_yy = adj1_sens_yy(:);
              [vec_adj1_spec_yy id1] = sort(adj1_spec_yy(:));
              vec_adj3_sens_yy = adj3_sens_yy(:);
              [vec_adj3_spec_yy id3] = sort(adj3_spec_yy(:));
              vec_adj4_sens_yy = adj4_sens_yy(:);
              [vec_adj4_spec_yy id4] = sort(adj4_spec_yy(:));
              vec_adj5_sens_yy = adj5_sens_yy(:);
              [vec_adj5_spec_yy id5] = sort(adj5_spec_yy(:));
              y1 = smooth(1-vec_adj1_spec_yy, vec_adj1_sens_yy(id1), f,'lowess');
              y3 = smooth(1-vec_adj3_spec_yy, vec_adj3_sens_yy(id3), f,'lowess'); 
              y4 = smooth(1-vec_adj4_spec_yy, vec_adj4_sens_yy(id4), f,'lowess');
              y5 = smooth(1-vec_adj5_spec_yy, vec_adj5_sens_yy(id5), f,'lowess'); 
              plot(1-vec_adj3_spec_yy, y3,'g-',1-vec_adj4_spec_yy, y4,'m--',1-vec_adj5_spec_yy, y5,'r-.',1-vec_adj1_spec_yy, y1, 'b:','LineWidth', 2);
              axis([0 1 0 1]);
xlabel('YY-FPR','FontSize',20);
ylabel('YY-TPR','FontSize',20);
              % zy plot
              subplot(2,2,4);
              vec_adj1_sens_zy = adj1_sens_zy(:);
              [vec_adj1_spec_zy id1] = sort(adj1_spec_zy(:));
              vec_adj3_sens_zy = adj3_sens_zy(:);
              [vec_adj3_spec_zy id3] = sort(adj3_spec_zy(:));
              vec_adj4_sens_zy = adj4_sens_zy(:);
              [vec_adj4_spec_zy id4] = sort(adj4_spec_zy(:));
              vec_adj5_sens_zy = adj5_sens_zy(:);
              [vec_adj5_spec_zy id5] = sort(adj5_spec_zy(:));
              y1 = smooth(1-vec_adj1_spec_zy, vec_adj1_sens_zy(id1), f,'lowess');
              y3 = smooth(1-vec_adj3_spec_zy, vec_adj3_sens_zy(id3), f,'lowess'); 
              y4 = smooth(1-vec_adj4_spec_zy, vec_adj4_sens_zy(id4), f,'lowess');
              y5 = smooth(1-vec_adj5_spec_zy, vec_adj5_sens_zy(id5), f,'lowess'); 
              plot(1-vec_adj3_spec_zy, y3,'g-',1-vec_adj4_spec_zy, y4,'m--',1-vec_adj5_spec_zy, y5,'r-.',1-vec_adj1_spec_zy, y1, 'b:','LineWidth', 2);
              axis([0 1 0 1]);
xlabel('ZY-FPR','FontSize',20);
ylabel('ZY-TPR','FontSize',20);
              % parameter plot
legend('Proposed method','Fellinghauer et al','Lee-Hastie','Regular lasso', 'Location', 'SouthEast');
         set(gcf, 'Position', [100 100 800 600])
         %mtit(sprintf('a = %.1f, \ b = %.1f', a, b), 'fontsize', 14)
         end
    end
end
         
  
