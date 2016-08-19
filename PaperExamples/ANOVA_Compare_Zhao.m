rng('default');
rng(1);

% case 4
%% generate beta

%% generate X
addpath('glmnet_matlab_MAC')

N = 121;
X = randn(N,55);
i = 11
for k = 2:10
	  X(:,i) = X(:,1).*X(:,k);
          i = i+1;
end
for k = 3:10
	  X(:,i) = X(:,2).*X(:,k);
          i = i+1;
end
for k = 4:10
	  X(:,i) = X(:,3).*X(:,k);
          i = i+1;
end
for k = 5:10
	  X(:,i) = X(:,4).*X(:,k);
          i = i+1;
end
for k = 6:10
	  X(:,i) = X(:,5).*X(:,k);
          i = i+1;
end
for k = 7:10
	  X(:,i) = X(:,6).*X(:,k);
          i = i+1;
end
for k = 8:10
	  X(:,i) = X(:,7).*X(:,k);
          i = i+1;
end
for k = 9:10
	  X(:,i) = X(:,8).*X(:,k);
          i = i+1;
end
for k = 10:10
	  X(:,i) = X(:,9).*X(:,k);
          i = i+1;
end
          

% case 4
beta0 = [7;2;1;1;5;0;0;4;2;0];

active_index = [1 2 3 4 11 12 13 20 21 28];

beta = zeros(55,1);

beta(active_index) = beta0;

eps = randn(N,1)*3.7;


Y = X*beta + eps;

options = glmnetSet;
options.lambda = linspace(20,0,2000);
options.standardize = false;
options.type = 'naive';
fit= glmnet(X, Y, 'gaussian', options);
glmnetPlot(fit);

w_fac = ones(10,1);
w_fac(11:55) = 3;
options.penalty_factor = w_fac;
fit2= glmnet(X, Y, 'gaussian', options);
glmnetPlot(fit2);

cd SLEP_package_4.1
root=cd;
addpath(genpath([root '/SLEP']));


opts=[];

% Starting point
opts.init=2;        % starting from a zero point

% Termination 
opts.tFlag=3;       % the relative change is less than opts.tol
opts.maxIter=5000;  % maximum number of iterations
opts.tol=1e-5;      % the tolerance parameter

% regularization
opts.rFlag=1;       % use ratio

% Normalization
opts.nFlag=0;       % without normalization


G=[1 11 12 13 14 15 16 17 18 19 2 11 20 21 22 23 24 25 26 27 3 12 20 28 29 30 31 32 33 34 4 13 21 28 35 36 37 38 39 40 5 14 22 29 35 41 42 43 44 45 6 15 23 30 36 41 46 47 48 49 7 16 24 31 37 42 46 50 51 52 8 17 25 32 38 43 47 50 53 54 9 18 26 33 39 44 48 51 53 55 10 19 27 34 40 45 49 52 54 55 11:55];
w=[ [1 10 1]', [11 20 1]', [21 30 1]', [31 40 1]',[41 50 1]', [51 60 1]', [61 70 1]', [71 80 1]', [81 90 1]', [91 100 1]' ,[ 101 101 1 ]' ,[ 102 102 1 ]' ,[ 103 103 1 ]' ,[ 104 104 1 ]' ,[ 105 105 1 ]' ,[ 106 106 1 ]' ,[ 107 107 1 ]' ,[ 108 108 1 ]' ,[ 109 109 1 ]' ,[ 110 110 1 ]' ,[ 111 111 1 ]' ,[ 112 112 1 ]' ,[ 113 113 1 ]' ,[ 114 114 1 ]' ,[ 115 115 1 ]' ,[ 116 116 1 ]' ,[ 117 117 1 ]' ,[ 118 118 1 ]' ,[ 119 119 1 ]' ,[ 120 120 1 ]' ,[ 121 121 1 ]' ,[ 122 122 1 ]' ,[ 123 123 1 ]' ,[ 124 124 1 ]' ,[ 125 125 1 ]' ,[ 126 126 1 ]' ,[ 127 127 1 ]' ,[ 128 128 1 ]' ,[ 129 129 1 ]' ,[ 130 130 1 ]' ,[ 131 131 1 ]' ,[ 132 132 1 ]' ,[ 133 133 1 ]' ,[ 134 134 1 ]' ,[ 135 135 1 ]' ,[ 136 136 1 ]' ,[ 137 137 1 ]' ,[ 138 138 1 ]' ,[ 139 139 1 ]' ,[ 140 140 1 ]' ,[ 141 141 1 ]' ,[ 142 142 1 ]' ,[ 143 143 1 ]' ,[ 144 144 1 ]' ,[ 145 145 1 ]']; 

%G=[1,5,6,7,2,5,8,9,3,6,8,10,4,7,9,10];
%w=[ [1 4 1]', [5 8 1]', [9 12 1]', [13 16 1]'];



opts.G=G;
opts.ind=w;


opts.rStartNum=100;

%----------------------- Run the code  -----------------------


glasso_lambda = linspace(1500,0,2000);
opts.maxIter2=1000;
opts.tol2=1e-8;
opts.flag2=2;
tic;
fit3_coef = zeros(55,2000);
for i=1:2000
    z = [0,glasso_lambda(i)];
    [x, funVal, ValueL]= overlapping_LeastR(X, Y, z, opts);
    fit3_coef(:,i) = x;
end
toc;

cd ..


fit1_coef = glmnetCoef(fit);
fit1_coef = fit1_coef(2:56,:);
fit2_coef = glmnetCoef(fit2);
fit2_coef = fit2_coef(2:56,:);



csvwrite('ZhaoCase4_weighted_full.csv',fit2_coef);
csvwrite('ZhaoCase4_Regular_full.csv',fit1_coef);
csvwrite('ZhaoCase4_Overlap_full.csv',fit3_coef);



%% case 5


rng('default');
rng(1);


%% generate beta

%% generate X
addpath('glmnet_matlab_MAC')

N = 121;
X = randn(N,55);
i = 11
for k = 2:10
	  X(:,i) = X(:,1).*X(:,k);
          i = i+1;
end
for k = 3:10
	  X(:,i) = X(:,2).*X(:,k);
          i = i+1;
end
for k = 4:10
	  X(:,i) = X(:,3).*X(:,k);
          i = i+1;
end
for k = 5:10
	  X(:,i) = X(:,4).*X(:,k);
          i = i+1;
end
for k = 6:10
	  X(:,i) = X(:,5).*X(:,k);
          i = i+1;
end
for k = 7:10
	  X(:,i) = X(:,6).*X(:,k);
          i = i+1;
end
for k = 8:10
	  X(:,i) = X(:,7).*X(:,k);
          i = i+1;
end
for k = 9:10
	  X(:,i) = X(:,8).*X(:,k);
          i = i+1;
end
for k = 10:10
	  X(:,i) = X(:,9).*X(:,k);
          i = i+1;
end
          

% case 5
beta0 = [7;2;1;1;7;7;7;2;2;1];

active_index = [1 2 3 4 11 12 13 20 21 28];

beta = zeros(55,1);

beta(active_index) = beta0;

eps = randn(N,1)*3.7;


Y = X*beta + eps;

options = glmnetSet;
options.lambda = linspace(20,0,2000);
options.standardize = false;
options.type = 'naive';
fit= glmnet(X, Y, 'gaussian', options);
glmnetPlot(fit);

w_fac = ones(10,1);
w_fac(11:55) = 3;
options.penalty_factor = w_fac;
fit2= glmnet(X, Y, 'gaussian', options);
glmnetPlot(fit2);

cd SLEP_package_4.1
root=cd;
addpath(genpath([root '/SLEP']));


opts=[];

% Starting point
opts.init=2;        % starting from a zero point

% Termination 
opts.tFlag=3;       % the relative change is less than opts.tol
opts.maxIter=5000;  % maximum number of iterations
opts.tol=1e-5;      % the tolerance parameter

% regularization
opts.rFlag=1;       % use ratio

% Normalization
opts.nFlag=0;       % without normalization


G=[1 11 12 13 14 15 16 17 18 19 2 11 20 21 22 23 24 25 26 27 3 12 20 28 29 30 31 32 33 34 4 13 21 28 35 36 37 38 39 40 5 14 22 29 35 41 42 43 44 45 6 15 23 30 36 41 46 47 48 49 7 16 24 31 37 42 46 50 51 52 8 17 25 32 38 43 47 50 53 54 9 18 26 33 39 44 48 51 53 55 10 19 27 34 40 45 49 52 54 55 11:55];
w=[ [1 10 1]', [11 20 1]', [21 30 1]', [31 40 1]',[41 50 1]', [51 60 1]', [61 70 1]', [71 80 1]', [81 90 1]', [91 100 1]' ,[ 101 101 1 ]' ,[ 102 102 1 ]' ,[ 103 103 1 ]' ,[ 104 104 1 ]' ,[ 105 105 1 ]' ,[ 106 106 1 ]' ,[ 107 107 1 ]' ,[ 108 108 1 ]' ,[ 109 109 1 ]' ,[ 110 110 1 ]' ,[ 111 111 1 ]' ,[ 112 112 1 ]' ,[ 113 113 1 ]' ,[ 114 114 1 ]' ,[ 115 115 1 ]' ,[ 116 116 1 ]' ,[ 117 117 1 ]' ,[ 118 118 1 ]' ,[ 119 119 1 ]' ,[ 120 120 1 ]' ,[ 121 121 1 ]' ,[ 122 122 1 ]' ,[ 123 123 1 ]' ,[ 124 124 1 ]' ,[ 125 125 1 ]' ,[ 126 126 1 ]' ,[ 127 127 1 ]' ,[ 128 128 1 ]' ,[ 129 129 1 ]' ,[ 130 130 1 ]' ,[ 131 131 1 ]' ,[ 132 132 1 ]' ,[ 133 133 1 ]' ,[ 134 134 1 ]' ,[ 135 135 1 ]' ,[ 136 136 1 ]' ,[ 137 137 1 ]' ,[ 138 138 1 ]' ,[ 139 139 1 ]' ,[ 140 140 1 ]' ,[ 141 141 1 ]' ,[ 142 142 1 ]' ,[ 143 143 1 ]' ,[ 144 144 1 ]' ,[ 145 145 1 ]']; 

%G=[1,5,6,7,2,5,8,9,3,6,8,10,4,7,9,10];
%w=[ [1 4 1]', [5 8 1]', [9 12 1]', [13 16 1]'];



opts.G=G;
opts.ind=w;


opts.rStartNum=100;

%----------------------- Run the code  -----------------------


glasso_lambda = linspace(1500,0,2000);
opts.maxIter2=1000;
opts.tol2=1e-8;
opts.flag2=2;
tic;
fit3_coef = zeros(55,2000);
for i=1:2000
    z = [0,glasso_lambda(i)];
    [x, funVal, ValueL]= overlapping_LeastR(X, Y, z, opts);
    fit3_coef(:,i) = x;
end
toc;

cd ..


fit1_coef = glmnetCoef(fit);
fit1_coef = fit1_coef(2:56,:);
fit2_coef = glmnetCoef(fit2);
fit2_coef = fit2_coef(2:56,:);



csvwrite('ZhaoCase5_weighted_full.csv',fit2_coef);
csvwrite('ZhaoCase5_Regular_full.csv',fit1_coef);
csvwrite('ZhaoCase5_Overlap_full.csv',fit3_coef);





%% case 3





rng('default');
rng(1);


%% generate beta

%% generate X
addpath('glmnet_matlab_MAC')

N = 121;
X = randn(N,55);
i = 11
for k = 2:10
	  X(:,i) = X(:,1).*X(:,k);
          i = i+1;
end
for k = 3:10
	  X(:,i) = X(:,2).*X(:,k);
          i = i+1;
end
for k = 4:10
	  X(:,i) = X(:,3).*X(:,k);
          i = i+1;
end
for k = 5:10
	  X(:,i) = X(:,4).*X(:,k);
          i = i+1;
end
for k = 6:10
	  X(:,i) = X(:,5).*X(:,k);
          i = i+1;
end
for k = 7:10
	  X(:,i) = X(:,6).*X(:,k);
          i = i+1;
end
for k = 8:10
	  X(:,i) = X(:,7).*X(:,k);
          i = i+1;
end
for k = 9:10
	  X(:,i) = X(:,8).*X(:,k);
          i = i+1;
end
for k = 10:10
	  X(:,i) = X(:,9).*X(:,k);
          i = i+1;
end
          

% case 3
beta0 = [7;2;1;1;1;0;0;0.5;0.4;0.1];

active_index = [1 2 3 4 11 12 13 20 21 28];

beta = zeros(55,1);

beta(active_index) = beta0;

eps = randn(N,1)*3.7;


Y = X*beta + eps;

options = glmnetSet;
options.lambda = linspace(20,0,2000);
options.standardize = false;
options.type = 'naive';
fit= glmnet(X, Y, 'gaussian', options);
glmnetPlot(fit);

w_fac = ones(10,1);
w_fac(11:55) = 3;
options.penalty_factor = w_fac;
fit2= glmnet(X, Y, 'gaussian', options);
glmnetPlot(fit2);

cd SLEP_package_4.1
root=cd;
addpath(genpath([root '/SLEP']));


opts=[];

% Starting point
opts.init=2;        % starting from a zero point

% Termination 
opts.tFlag=3;       % the relative change is less than opts.tol
opts.maxIter=5000;  % maximum number of iterations
opts.tol=1e-5;      % the tolerance parameter

% regularization
opts.rFlag=1;       % use ratio

% Normalization
opts.nFlag=0;       % without normalization


G=[1 11 12 13 14 15 16 17 18 19 2 11 20 21 22 23 24 25 26 27 3 12 20 28 29 30 31 32 33 34 4 13 21 28 35 36 37 38 39 40 5 14 22 29 35 41 42 43 44 45 6 15 23 30 36 41 46 47 48 49 7 16 24 31 37 42 46 50 51 52 8 17 25 32 38 43 47 50 53 54 9 18 26 33 39 44 48 51 53 55 10 19 27 34 40 45 49 52 54 55 11:55];
w=[ [1 10 1]', [11 20 1]', [21 30 1]', [31 40 1]',[41 50 1]', [51 60 1]', [61 70 1]', [71 80 1]', [81 90 1]', [91 100 1]' ,[ 101 101 1 ]' ,[ 102 102 1 ]' ,[ 103 103 1 ]' ,[ 104 104 1 ]' ,[ 105 105 1 ]' ,[ 106 106 1 ]' ,[ 107 107 1 ]' ,[ 108 108 1 ]' ,[ 109 109 1 ]' ,[ 110 110 1 ]' ,[ 111 111 1 ]' ,[ 112 112 1 ]' ,[ 113 113 1 ]' ,[ 114 114 1 ]' ,[ 115 115 1 ]' ,[ 116 116 1 ]' ,[ 117 117 1 ]' ,[ 118 118 1 ]' ,[ 119 119 1 ]' ,[ 120 120 1 ]' ,[ 121 121 1 ]' ,[ 122 122 1 ]' ,[ 123 123 1 ]' ,[ 124 124 1 ]' ,[ 125 125 1 ]' ,[ 126 126 1 ]' ,[ 127 127 1 ]' ,[ 128 128 1 ]' ,[ 129 129 1 ]' ,[ 130 130 1 ]' ,[ 131 131 1 ]' ,[ 132 132 1 ]' ,[ 133 133 1 ]' ,[ 134 134 1 ]' ,[ 135 135 1 ]' ,[ 136 136 1 ]' ,[ 137 137 1 ]' ,[ 138 138 1 ]' ,[ 139 139 1 ]' ,[ 140 140 1 ]' ,[ 141 141 1 ]' ,[ 142 142 1 ]' ,[ 143 143 1 ]' ,[ 144 144 1 ]' ,[ 145 145 1 ]']; 

%G=[1,5,6,7,2,5,8,9,3,6,8,10,4,7,9,10];
%w=[ [1 4 1]', [5 8 1]', [9 12 1]', [13 16 1]'];



opts.G=G;
opts.ind=w;


opts.rStartNum=100;

%----------------------- Run the code  -----------------------


glasso_lambda = linspace(1500,0,2000);
opts.maxIter2=1000;
opts.tol2=1e-8;
opts.flag2=2;
tic;
fit3_coef = zeros(55,2000);
for i=1:2000
    z = [0,glasso_lambda(i)];
    [x, funVal, ValueL]= overlapping_LeastR(X, Y, z, opts);
    fit3_coef(:,i) = x;
end
toc;

cd ..


fit1_coef = glmnetCoef(fit);
fit1_coef = fit1_coef(2:56,:);
fit2_coef = glmnetCoef(fit2);
fit2_coef = fit2_coef(2:56,:);


csvwrite('ZhaoCase3_weighted_full.csv',fit2_coef);
csvwrite('ZhaoCase3_Regular_full.csv',fit1_coef);
csvwrite('ZhaoCase3_Overlap_full.csv',fit3_coef);





% case 2




rng('default');
rng(1);


%% generate beta

%% generate X
addpath('glmnet_matlab_MAC')

N = 121;
X = randn(N,55);
i = 11
for k = 2:10
	  X(:,i) = X(:,1).*X(:,k);
          i = i+1;
end
for k = 3:10
	  X(:,i) = X(:,2).*X(:,k);
          i = i+1;
end
for k = 4:10
	  X(:,i) = X(:,3).*X(:,k);
          i = i+1;
end
for k = 5:10
	  X(:,i) = X(:,4).*X(:,k);
          i = i+1;
end
for k = 6:10
	  X(:,i) = X(:,5).*X(:,k);
          i = i+1;
end
for k = 7:10
	  X(:,i) = X(:,6).*X(:,k);
          i = i+1;
end
for k = 8:10
	  X(:,i) = X(:,7).*X(:,k);
          i = i+1;
end
for k = 9:10
	  X(:,i) = X(:,8).*X(:,k);
          i = i+1;
end
for k = 10:10
	  X(:,i) = X(:,9).*X(:,k);
          i = i+1;
end
          

% case 2
beta0 = [7;2;1;1;0.5;0;0;0.1;0.1;0];

active_index = [1 2 3 4 11 12 13 20 21 28];

beta = zeros(55,1);

beta(active_index) = beta0;

eps = randn(N,1)*3.7;


Y = X*beta + eps;

options = glmnetSet;
options.lambda = linspace(20,0,2000);
options.standardize = false;
options.type = 'naive';
fit= glmnet(X, Y, 'gaussian', options);
glmnetPlot(fit);

w_fac = ones(10,1);
w_fac(11:55) = 3;
options.penalty_factor = w_fac;
fit2= glmnet(X, Y, 'gaussian', options);
glmnetPlot(fit2);

cd SLEP_package_4.1
root=cd;
addpath(genpath([root '/SLEP']));


opts=[];

% Starting point
opts.init=2;        % starting from a zero point

% Termination 
opts.tFlag=3;       % the relative change is less than opts.tol
opts.maxIter=5000;  % maximum number of iterations
opts.tol=1e-5;      % the tolerance parameter

% regularization
opts.rFlag=1;       % use ratio

% Normalization
opts.nFlag=0;       % without normalization


G=[1 11 12 13 14 15 16 17 18 19 2 11 20 21 22 23 24 25 26 27 3 12 20 28 29 30 31 32 33 34 4 13 21 28 35 36 37 38 39 40 5 14 22 29 35 41 42 43 44 45 6 15 23 30 36 41 46 47 48 49 7 16 24 31 37 42 46 50 51 52 8 17 25 32 38 43 47 50 53 54 9 18 26 33 39 44 48 51 53 55 10 19 27 34 40 45 49 52 54 55 11:55];
w=[ [1 10 1]', [11 20 1]', [21 30 1]', [31 40 1]',[41 50 1]', [51 60 1]', [61 70 1]', [71 80 1]', [81 90 1]', [91 100 1]' ,[ 101 101 1 ]' ,[ 102 102 1 ]' ,[ 103 103 1 ]' ,[ 104 104 1 ]' ,[ 105 105 1 ]' ,[ 106 106 1 ]' ,[ 107 107 1 ]' ,[ 108 108 1 ]' ,[ 109 109 1 ]' ,[ 110 110 1 ]' ,[ 111 111 1 ]' ,[ 112 112 1 ]' ,[ 113 113 1 ]' ,[ 114 114 1 ]' ,[ 115 115 1 ]' ,[ 116 116 1 ]' ,[ 117 117 1 ]' ,[ 118 118 1 ]' ,[ 119 119 1 ]' ,[ 120 120 1 ]' ,[ 121 121 1 ]' ,[ 122 122 1 ]' ,[ 123 123 1 ]' ,[ 124 124 1 ]' ,[ 125 125 1 ]' ,[ 126 126 1 ]' ,[ 127 127 1 ]' ,[ 128 128 1 ]' ,[ 129 129 1 ]' ,[ 130 130 1 ]' ,[ 131 131 1 ]' ,[ 132 132 1 ]' ,[ 133 133 1 ]' ,[ 134 134 1 ]' ,[ 135 135 1 ]' ,[ 136 136 1 ]' ,[ 137 137 1 ]' ,[ 138 138 1 ]' ,[ 139 139 1 ]' ,[ 140 140 1 ]' ,[ 141 141 1 ]' ,[ 142 142 1 ]' ,[ 143 143 1 ]' ,[ 144 144 1 ]' ,[ 145 145 1 ]']; 

%G=[1,5,6,7,2,5,8,9,3,6,8,10,4,7,9,10];
%w=[ [1 4 1]', [5 8 1]', [9 12 1]', [13 16 1]'];



opts.G=G;
opts.ind=w;


opts.rStartNum=100;

%----------------------- Run the code  -----------------------


glasso_lambda = linspace(1500,0,2000);
opts.maxIter2=1000;
opts.tol2=1e-8;
opts.flag2=2;
tic;
fit3_coef = zeros(55,2000);
for i=1:2000
    z = [0,glasso_lambda(i)];
    [x, funVal, ValueL]= overlapping_LeastR(X, Y, z, opts);
    fit3_coef(:,i) = x;
end
toc;

cd ..


fit1_coef = glmnetCoef(fit);
fit1_coef = fit1_coef(2:56,:);
fit2_coef = glmnetCoef(fit2);
fit2_coef = fit2_coef(2:56,:);


csvwrite('ZhaoCase2_weighted_full.csv',fit2_coef);
csvwrite('ZhaoCase2_Regular_full.csv',fit1_coef);
csvwrite('ZhaoCase2_Overlap_full.csv',fit3_coef);





% case 1




rng('default');
rng(1);


%% generate beta

%% generate X
addpath('glmnet_matlab_MAC')

N = 121;
X = randn(N,55);
i = 11
for k = 2:10
	  X(:,i) = X(:,1).*X(:,k);
          i = i+1;
end
for k = 3:10
	  X(:,i) = X(:,2).*X(:,k);
          i = i+1;
end
for k = 4:10
	  X(:,i) = X(:,3).*X(:,k);
          i = i+1;
end
for k = 5:10
	  X(:,i) = X(:,4).*X(:,k);
          i = i+1;
end
for k = 6:10
	  X(:,i) = X(:,5).*X(:,k);
          i = i+1;
end
for k = 7:10
	  X(:,i) = X(:,6).*X(:,k);
          i = i+1;
end
for k = 8:10
	  X(:,i) = X(:,7).*X(:,k);
          i = i+1;
end
for k = 9:10
	  X(:,i) = X(:,8).*X(:,k);
          i = i+1;
end
for k = 10:10
	  X(:,i) = X(:,9).*X(:,k);
          i = i+1;
end
          

% case 1
beta0 = [7;2;1;1;0;0;0;0;0;0];

active_index = [1 2 3 4 11 12 13 20 21 28];

beta = zeros(55,1);

beta(active_index) = beta0;

eps = randn(N,1)*3.7;


Y = X*beta + eps;

options = glmnetSet;
options.lambda = linspace(20,0,2000);
options.standardize = false;
options.type = 'naive';
fit= glmnet(X, Y, 'gaussian', options);
glmnetPlot(fit);

w_fac = ones(10,1);
w_fac(11:55) = 3;
options.penalty_factor = w_fac;
fit2= glmnet(X, Y, 'gaussian', options);
glmnetPlot(fit2);

cd SLEP_package_4.1
root=cd;
addpath(genpath([root '/SLEP']));


opts=[];

% Starting point
opts.init=2;        % starting from a zero point

% Termination 
opts.tFlag=3;       % the relative change is less than opts.tol
opts.maxIter=5000;  % maximum number of iterations
opts.tol=1e-5;      % the tolerance parameter

% regularization
opts.rFlag=1;       % use ratio

% Normalization
opts.nFlag=0;       % without normalization


G=[1 11 12 13 14 15 16 17 18 19 2 11 20 21 22 23 24 25 26 27 3 12 20 28 29 30 31 32 33 34 4 13 21 28 35 36 37 38 39 40 5 14 22 29 35 41 42 43 44 45 6 15 23 30 36 41 46 47 48 49 7 16 24 31 37 42 46 50 51 52 8 17 25 32 38 43 47 50 53 54 9 18 26 33 39 44 48 51 53 55 10 19 27 34 40 45 49 52 54 55 11:55];
w=[ [1 10 1]', [11 20 1]', [21 30 1]', [31 40 1]',[41 50 1]', [51 60 1]', [61 70 1]', [71 80 1]', [81 90 1]', [91 100 1]' ,[ 101 101 1 ]' ,[ 102 102 1 ]' ,[ 103 103 1 ]' ,[ 104 104 1 ]' ,[ 105 105 1 ]' ,[ 106 106 1 ]' ,[ 107 107 1 ]' ,[ 108 108 1 ]' ,[ 109 109 1 ]' ,[ 110 110 1 ]' ,[ 111 111 1 ]' ,[ 112 112 1 ]' ,[ 113 113 1 ]' ,[ 114 114 1 ]' ,[ 115 115 1 ]' ,[ 116 116 1 ]' ,[ 117 117 1 ]' ,[ 118 118 1 ]' ,[ 119 119 1 ]' ,[ 120 120 1 ]' ,[ 121 121 1 ]' ,[ 122 122 1 ]' ,[ 123 123 1 ]' ,[ 124 124 1 ]' ,[ 125 125 1 ]' ,[ 126 126 1 ]' ,[ 127 127 1 ]' ,[ 128 128 1 ]' ,[ 129 129 1 ]' ,[ 130 130 1 ]' ,[ 131 131 1 ]' ,[ 132 132 1 ]' ,[ 133 133 1 ]' ,[ 134 134 1 ]' ,[ 135 135 1 ]' ,[ 136 136 1 ]' ,[ 137 137 1 ]' ,[ 138 138 1 ]' ,[ 139 139 1 ]' ,[ 140 140 1 ]' ,[ 141 141 1 ]' ,[ 142 142 1 ]' ,[ 143 143 1 ]' ,[ 144 144 1 ]' ,[ 145 145 1 ]']; 

%G=[1,5,6,7,2,5,8,9,3,6,8,10,4,7,9,10];
%w=[ [1 4 1]', [5 8 1]', [9 12 1]', [13 16 1]'];



opts.G=G;
opts.ind=w;


opts.rStartNum=100;

%----------------------- Run the code  -----------------------


glasso_lambda = linspace(1500,0,2000);
opts.maxIter2=1000;
opts.tol2=1e-8;
opts.flag2=2;
tic;
fit3_coef = zeros(55,2000);
for i=1:2000
    z = [0,glasso_lambda(i)];
    [x, funVal, ValueL]= overlapping_LeastR(X, Y, z, opts);
    fit3_coef(:,i) = x;
end
toc;

cd ..


fit1_coef = glmnetCoef(fit);
fit1_coef = fit1_coef(2:56,:);
fit2_coef = glmnetCoef(fit2);
fit2_coef = fit2_coef(2:56,:);


csvwrite('ZhaoCase1_weighted_full.csv',fit2_coef);
csvwrite('ZhaoCase1_Regular_full.csv',fit1_coef);
csvwrite('ZhaoCase1_Overlap_full.csv',fit3_coef);
