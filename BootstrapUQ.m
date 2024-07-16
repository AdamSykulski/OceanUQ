%% GENERATING SOME TOY DATA FROM A STUDENT T-DISTRIBUTION WITH MEAN 0.5 %%%
close all
clear all
N = 100; % sample size
dof = 5; % degrees of freedom
rng(4); X = 0.5+trnd(dof,[N,1]); % sample from a t-distribution with mean 0.5
figure; histogram(X,20,'Normalization','pdf') % let's look at a plot
hold on; plot(min(X):0.01:max(X),tpdf((min(X):0.01:max(X))-0.5,dof),'LineWidth',2)
hold on; plot(min(X):0.01:max(X),normpdf((min(X):0.01:max(X))-0.5),'LineWidth',2)
legend('data','student-t','normal')
mean_of_sample = mean(X) % the mean, but what is the uncertainty of our mean estimate?



%% IDEAL SCENARIO: WE KNOW THE TRUE STANDARD DEVIATION OF THE DATA
standard_error_theory = sqrt(dof/(dof-2))/sqrt(N) % the standard error of the mean












%% 2nd IDEAL SCENARIO: REPEAT EXPERIMENTS
MM = 1000; % number of repeats
Y = 0.5+trnd(dof,[N,MM]); % Monte Carlo repeats
m = mean(Y); % lots of Monte Carlo means
standard_error_MC=std(m) % the Monte Carlo standard error









%% OK but we only have X, let's try the bootstrap!
BB = 1000; mb = zeros(1,BB);
for jj = 1:BB
    mb(jj) = mean(datasample(X,N));
end
standard_error_bootstrap=std(mb) % the bootstrap standard error








%% and how about confidence intervals?
figure; histogram(mb,20,'Normalization','pdf') % let's look at a plot
hold on; plot(min(mb):0.01:max(mb),normpdf(min(mb):0.01:max(mb),mean(mb),std(mb)),'LineWidth',2)
legend('boostrap estimates','normal')
%% order the data and take percentiles!
bootstrap_confidence_interval=[prctile(mb,2.5), prctile(mb,97.5)]








%% OK but did we get lucky? let's try running it on all the samples
for ii = 1:MM % Monte Carlo
    for jj = 1:BB % bootstrap
        mb(jj) = mean(datasample(Y(:,ii),N)); % the mean
    end
    bste(ii) = std(mb);
end
figure; histogram(bste,100); hold on; line([standard_error_theory standard_error_theory],[0 100],'linewidth',2); xlabel("standard error estimates"); legend("bootstrap","truth") % a figure comparing all the boostrap estimates
fontsize(gca,20,"pixels")
legend('bootstrap','truth')




%% OK what if the data is correlated?
N = 100;
corr = 0.3; % correlation of neighbouring data points
A = toeplitz([1 corr zeros(1,N-2)]); % Cholesky method for generating correlated data
AA = chol(A);
for ii = 1:MM % Monte Carlo
    X = 0.5+trnd(dof,[N,1]);
    X = AA'*(X-0.5) + 0.5;
    m(ii) = mean(X); % the mean
    for jj = 1:BB % Bootstrap
        mb(jj) = mean(datasample(X,N));
    end
    bste(ii) = std(mb);
end
figure; histogram(bste,100); hold on; line([std(m) std(m)],[0 100],'linewidth',2); xlabel("standard error estimates"); legend("bootstrap","truth")
fontsize(gca,20,"pixels")
% time series bootstrap
for ii = 1:MM % Monte Carlo
X = 0.5+trnd(dof,[N,1]);
X = AA'*(X-0.5) + 0.5;
m(ii) = mean(X); % the mean
NB = 5; % block length
for jj = 1:BB % The bootstrap
    bb=datasample(1:N/NB,N/NB);
    for kk = 1:length(bb)
        XB((kk-1)*NB+1:kk*NB) = X((bb(kk)-1)*NB+1:bb(kk)*NB);
    end
    mb(jj) = mean(XB);
end
bste(ii) = std(mb);
end
figure; histogram(bste,100); hold on; line([std(m) std(m)],[0 100],'linewidth',2); xlabel("standard error estimates"); legend("block bootstrap","truth")
fontsize(gca,20,"pixels")