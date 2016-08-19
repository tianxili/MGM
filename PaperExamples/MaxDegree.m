
addpath('glmnet_matlab_MAC')

n = 100; % sample size
p = 90;
q = 10;
K = 100;

df = [2 6 10];
    a = 1;
    b = 2;
    D = 17; % depth of lambda log_0.5
    L = 100; % length of lambda
    lambda = 0.5.^linspace(-2,D,L);
    kappa_seq = 0.1:0.1:1;
     rng('default');
    seed = 50;
    rng(seed);
%% Model fitting 
for i = 1:3
%      name = sprintf('new_adj_nedge80_dfbd%d.txt',df(i));
%      adj = dlmread(name);
     rng('default');
    seed = 50;
    rng(seed);
          valid = 0;
     while(valid==0)
         adj = adjgen(n, 80, df(i));
         if sum(sum(adj(1:q,1:q)))>0 && max(sum(adj))>df(i)-1
             valid=1;
         end
     end  
         parlist = pargen_unif(adj, p, q, a, b);
         tmp = edgecount(parlist);
         tmp1 = reg_count(parlist);
         count = sum(parlist.lambda_jk~=0)+sum(sum(parlist.ita_j~=0))+...
             (sum(sum(parlist.phi_0~=0))-p+sum(sum(sum(parlist.phi_j~=0))))/2;
         idv_sens = zeros(L,K);
         idv_spec = zeros(L,K);
         adj_sens = zeros(L,K);
         adj_spec = zeros(L,K);
         adj_sens_zz = zeros(L,K);
         adj_spec_zz = zeros(L,K);
         adj_sens_yy = zeros(L,K);
         adj_spec_yy = zeros(L,K);
         adj_sens_zy = zeros(L,K);
         adj_spec_zy = zeros(L,K);
         %replications
        %matlabpool open;
        for k = 1:K
              k
              [z y] = datagen(parlist, n);                           
              %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
              %%%%%%%% FIT the Weighted Lasso%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
              
                  [fitlist_post fitlist]= sepreg_weight_flex(z,y,lambda, lambda, 'max',kappa_seq(1));
                  %[fitlist_post fitlist]= sepreg_weight(z,y,lambda3, 10*lambda3, 'max');
[idv_sens(:,k) idv_spec(:,k) idv_tp idv_tn] = fitcompare(fitlist_post, parlist);
                  adj_norm = edgenorm(fitlist_post);
                  adj_lambda = edgethres(adj_norm, 0);
                  %[adj3_sens(:,k) adj3_spec(:,k) adj_tp adj_tn] = adjcompare(adj_lambda, adj);
[adj_sens(:,k) adj_spec(:,k) total_p total_n adj_sens_zz(:,k) adj_spec_zz(:,k) total_p_zz total_n_zz adj_sens_yy(:,k) adj_spec_yy(:,k) total_p_yy total_n_yy adj_sens_zy(:,k) adj_spec_zy(:,k) total_p_zy total_n_zy ] = adjcompare_category(adj_lambda, adj,p,q);
              end
              name = sprintf('Simulation/MaxStudy_dgmax_dg=%d.mat', df(i));
        save(name, 'idv_sens', 'idv_spec', 'adj_sens', 'adj_spec', 'idv_tp', 'idv_tn', 'total_p',  'total_n',...
            'adj_sens_zz','adj_spec_zz','total_p_zz','total_n_zz',...
            'adj_sens_yy','adj_spec_yy','total_p_yy','total_n_yy',...
            'adj_sens_zy','adj_spec_zy','total_p_zy','total_n_zy');
              
end            
 % Plotting the results. 
 idv_tp = [80 80 81];
 idv_tn = [44920 44920 44919];
 adj_tp = 80;
 adj_tn = 4870;
         f = 1/4;
         g = figure;
         df = [2 6 10];
         screen_size = get(0,'ScreenSize');
         set(g, 'Position', [0 0 screen_size(4) screen_size(4)] );
         tt = {'Parameter' 'Parameter' 'Edge' 'Edge'};
         symb = {'r-', 'g-.', 'b--'};
         for i = 1:3
              name = sprintf('Simulation/MaxStudy_dgmax_dg=%d.mat', df(i));
              load(name);
              vec_adj_sens = adj_sens(:);
              [vec_adj_spec id1] = sort(adj_spec(:));
              vec_idv_sens = idv_sens(:);
              [vec_idv_spec id2] = sort(idv_spec(:));
              vec_adj_sens_zz = adj_sens_zz(:);
              [vec_adj_spec_zz id3] = sort(adj_spec_zz(:));
              vec_adj_sens_zy = adj_sens_zy(:);
              [vec_adj_spec_zy id4] = sort(adj_spec_zy(:));
              vec_adj_sens_yy = adj_sens_yy(:);
              [vec_adj_spec_yy id5] = sort(adj_spec_yy(:));
              y1 = smooth(1-vec_adj_spec, vec_adj_sens(id1), f,'lowess');
              y2 = smooth(1-vec_idv_spec, vec_idv_sens(id2), f,'lowess');
              y3 = smooth(1-vec_adj_spec_zz, vec_adj_sens_zz(id3), f,'lowess');
              y4 = smooth(1-vec_adj_spec_zy, vec_adj_sens_zy(id4), f,'lowess');
              y5 = smooth(1-vec_adj_spec_yy, vec_adj_sens_yy(id5), f,'lowess');
              %%%%%%%%%%%%%%%%%%%%%%%%%%%%%
              subplot(1,2,1)
              plot(1-vec_idv_spec, y2, char(symb(i)),'LineWidth', 2 );
              axis([0 0.1 0 1]);
title(char(tt(1)),'FontSize',15);
              xlabel('FPR','FontSize',15);
              ylabel('TPR','FontSize',15)
              hold on;
              subplot(1,2,2)
              plot(1-vec_adj_spec, y1, char(symb(i)),'LineWidth', 2 );
              axis([0 1 0 1]);
              title(char(tt(3)),'FontSize',15);
              xlabel('FPR','FontSize',15);
              ylabel('TPR','FontSize',15)
              hold on;              
         end
         h_legend=legend('dg\_max = 2', 'dg\_max = 6', 'dg\_max = 10', 'Location', 'SouthEast')
         set(h_legend,'FontSize',15);
         set(gcf, 'Position', [100 100 1000 300]) 
         

 
 idv_tp = [80 80 81];
 idv_tn = [44920 44920 44919];
 adj_tp = 80;
 adj_tn = 4870;
         f = 1/5;
         g = figure;
         df = [2 6 10];
         screen_size = get(0,'ScreenSize');
         set(g, 'Position', [0 0 screen_size(4) screen_size(4)] );
         tt = {'Parameter' 'Parameter (Count)' 'Edge (Rate)' 'Edge (Count)'};
         symb = {'r-', 'g-.', 'b--'};
         for i = 1:3
              name = sprintf('Simulation/MaxStudy_dgmax_dg=%d.mat', df(i));
              load(name);
              vec_adj_sens = adj_sens(:);
              [vec_adj_spec id1] = sort(adj_spec(:));
              vec_idv_sens = idv_sens(:);
              [vec_idv_spec id2] = sort(idv_spec(:));
              vec_adj_sens_zz = adj_sens_zz(:);
              [vec_adj_spec_zz id3] = sort(adj_spec_zz(:));
              vec_adj_sens_zy = adj_sens_zy(:);
              [vec_adj_spec_zy id4] = sort(adj_spec_zy(:));
              vec_adj_sens_yy = adj_sens_yy(:);
              [vec_adj_spec_yy id5] = sort(adj_spec_yy(:));
              y1 = smooth(1-vec_adj_spec, vec_adj_sens(id1), f,'lowess');
              y2 = smooth(1-vec_idv_spec, vec_idv_sens(id2), f,'lowess');
              y3 = smooth(1-vec_adj_spec_zz, vec_adj_sens_zz(id3), f,'lowess');
              y4 = smooth(1-vec_adj_spec_zy, vec_adj_sens_zy(id4), f,'lowess');
              y5 = smooth(1-vec_adj_spec_yy, vec_adj_sens_yy(id5), f,'lowess');
              %%%%%%%%%%%%%%%%%%%%%%%%%%%%%
              subplot(1,3,1)
              plot(1-vec_adj_spec_zz, y3, char(symb(i)),'LineWidth', 2 );
              axis([0 1 0 1]);
              %title(char(tt(3)));
              xlabel('ZZ-FPR','FontSize',15);
              ylabel('ZZ-TPR','FontSize',15)
              hold on;             
              subplot(1,3,2)
              plot(1-vec_adj_spec_yy, y5, char(symb(i)),'LineWidth', 2 );
              axis([0 1 0 1]);
              %title(char(tt(3)));
              xlabel('YY-FPR','FontSize',15);
              ylabel('YY-TPR','FontSize',15)
              hold on;  
              subplot(1,3,3)
              plot(1-vec_adj_spec_zy, y4, char(symb(i)),'LineWidth', 2 );
              axis([0 1 0 1]);
              %title(char(tt(3)));
              xlabel('ZY-FPR','FontSize',15);
              ylabel('ZY-TPR','FontSize',15)
              hold on;  
         end
         h_legend = legend('dg\_max = 2', 'dg\_max = 6', 'dg\_max = 10', 'Location', 'SouthEast')
         set(h_legend,'FontSize',15);
          set(gcf, 'Position', [100 100 1200 300]) 
         

        
