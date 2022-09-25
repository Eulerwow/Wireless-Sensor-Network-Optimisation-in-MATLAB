%% script to call Matlab built-in Optimiser and Penalty method with BFGS (3D), 
% with .
clear all;
clc;
numSen = 2; %number of sensors
X = (-10+20*rand(numSen,2)); %list of sensors
numInitials = 10; %number of initial points
S0 = [-10+20*rand(numInitials,2),100.*rand(numInitials,1)]; %list of initial points
tolerance1 = 10^(-2); %tolerance for BFGS
tolerance2 = 10^(-5); % tolerance for Golden section search
tolerance3 = 10^(-4); % tolerance for Penalty method 
T = 10; %step size
%%

%% Matlab Optimiser
for i = 1:numInitials
    MLtStart = tic; %start recording time
    s0 = S0(i,1:2);
    sf = @(s)f(s,X);
    %options = optimoptions(@fminunc, 'OptimalityTolerance', 1e-5);
    [MLmins(i,:),MLminf(i)]=fminunc(sf,s0);%,options);

    MLtElapsed(i) = toc(MLtStart); %end recording time
end
MLminf=round(MLminf,2);
MLmatrix = [MLminf',round(MLmins,2), MLtElapsed']; %combine the results data into a single matrix
MLmatrix = sortrows(MLmatrix);
MLSmins = unique(MLmatrix(:,2:3),'rows','stable'); %extract unique minimisers

%calculate f value appearance times, average iterations and search times
MLrep = zeros(1,size(MLSmins,1));
MLtime = zeros(1,size(MLSmins,1));
for n = 1:size(MLSmins,1)
    MLfs(n) = feval(sf,MLSmins(n,:));
    for m = 1:numInitials
        if MLmatrix(m,2:3) == MLSmins(n,:)
            MLrep(n) = MLrep(n) + 1;
            MLtime(n) = MLtime(n)+ MLmatrix(m,4);
        end
    end
end
aveMLtime = MLtime./MLrep;

%print performance analysis in a table
varNames = {'P_value','Minimiser','No_of_times','Ave_search_time'};
MLTb = table(MLfs',MLSmins,MLrep',aveMLtime','VariableNames',varNames);
disp('Matlab Optimiser');
MLTb



%% l2 penalty
for i = 1:numInitials
    L2tStart = tic;  % internal timer
    pk=1; %penalty parameter
    s0 = S0(i,:);
    Pval_old = feval(@P,pk,s0,X); % P function value at starting point s0
    [BFGSs(pk,:),Pval(pk),~] = BFGS('P', 'gradP', pk, s0, X, eye(3), tolerance1, tolerance2, T); %use BFGS to find next point
    Pval_new = Pval(pk);
    while ( abs(Pval_old-Pval_new) >= tolerance3 ) 
        pk=pk+1; %increase penalty parameter when P value not converge enough
        Pval_old = Pval_new;
        [BFGSs(pk,:),Pval(pk),~] = BFGS('P', 'gradP', pk, BFGSs(pk-1,:), X, eye(3), tolerance1, tolerance2, T) ;
        Pval_new = Pval(pk); 
    end
    L2tElapsed(i) = toc(L2tStart);
    L2pk(i) = pk;%count pk
    L2smin(i,:) = BFGSs(pk,:);
    L2Pval(i) = Pval_new; %
      
end

L2matrix = [L2Pval',round(L2smin,2),L2pk',L2tElapsed']; %combine the results data into a single 'numInstances x 6' matrix
L2smins = unique(L2matrix(:,2:4),'rows'); %extract unique minimum P values

%calculate P value appearance times, average iterations and search times
L2rep = zeros(1,size(L2smins,1));
L2iter = zeros(1,size(L2smins,1));
L2time = zeros(1,size(L2smins,1));

for n = 1:size(L2smins,1)
    L2Ps(n) = feval(@P,L2pk(n),L2smins(n,:),X);
    L2timeDetail = [];
    for m = 1:numInitials
        if L2matrix(m,2:4) == L2smins(n,:)
            L2rep(n) = L2rep(n) + 1;
            L2iter(n) = L2iter(n)+ L2matrix(m,5);
            L2time(n) = L2time(n)+ L2matrix(m,6);
            L2timeDetail = [L2timeDetail; L2matrix(m,6)];       
        end
    end
    stdL2time(n) = std( L2timeDetail);
end
aveL2iter = L2iter./L2rep;
aveL2time = L2time./L2rep;
%print performance analysis in a table
varNames = {'P_value','Minimiser','No_of_times','Ave_iterations','Ave_search_time','Search_time_STD'};
L2Tb = table(L2Ps',L2smins,L2rep',aveL2iter',aveL2time',stdL2time','VariableNames',varNames);
disp('l2 Penalty Method');
L2Tb