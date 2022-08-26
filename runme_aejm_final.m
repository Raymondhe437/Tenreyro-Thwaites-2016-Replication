%% File to produce results in published version of Tenreyro & Thwaites
% GDT 10 December 2015

clear all
close all

%% construct the data and shocks
construct_data

% construct the transition function and the shocks using different state variables
generate_shocks

%% Create the baseline set of data and parameters

% Load data
clear all
close all
load data.mat R LHS names shocksandstates startS endS samplelength
runstate=1;                 % run state-dependent regressions as well as linear ones? 1=yes
if runstate==1              % State variable for state-dependent regression
    M=shocksandstates(1).M(startS:endS,1);     % baseline state variable
end

% Parameters
includetrend=1;             % 1 includes a trend
useNLshocks=1;              % 1 uses non-linearly identified shocks
H=20;                       % Max horizon of IRF
blocklength=H;              % Size of blocks for bootstrap
CI=.90;                     % Confidence intervals in charts (as decimal fraction)
Xlags=1;                    % Number of lags of control variables
Rlags=1;                    % Number of lags of policy variables
B=10000;                    % number of bootstrap iterations
scaling=1;                  % whether scaled parameters are calculated (set to zero to save time if you don't need them)

% weights for smoothing the results
weight=zeros(H+1,H+1);
weight(1,1)=1; % End-points are not MAs
weight(H+1,H+1)=1;
for h=2:H
    weight(h,h-1:h+1)=[1/3 1/3 1/3]; % Centred three-period moving average
end
cumweight=tril(ones(H+1,H+1));
trend=linspace(1,samplelength,samplelength)';

save baseline.mat           % baseline workspace

%% Run the baseline simulation
clear all
close all
load baseline.mat
LHSindices=[1 2 7 10 13 22 25 28 15 34 35 38];
LHS=LHS(:,LHSindices);
[results]=stlpm(R,shocksandstates(1).E,shocksandstates(1).M_Z_t,LHS,trend,Xlags,Rlags,runstate,includetrend,useNLshocks,H,CI,startS,endS,scaling,weight,cumweight)
stlpm_inference
block_bootstrap
save baseline_results.mat

%% Linear identification
clear all
close all
load baseline.mat
LHSindices=[1 2];
LHS=LHS(:,LHSindices);
useNLshocks=0;
[results]=stlpm(R,shocksandstates(1).E,shocksandstates(1).M_Z_t,LHS,trend,Xlags,Rlags,runstate,includetrend,useNLshocks,H,CI,startS,endS,scaling,weight,cumweight)
stlpm_inference
block_bootstrap
save linearid_results.mat

%% VAR identification
clear all
close all
load baseline.mat
LHSindices=[1 2];
LHS=LHS(:,LHSindices);
[results]=stlpm(R,shocksandstates(7).E,shocksandstates(7).M_Z_t,LHS,trend,Xlags,Rlags,runstate,includetrend,useNLshocks,H,CI,startS,endS,scaling,weight,cumweight)
stlpm_inference
block_bootstrap
save var_results.mat

%% NBER recessions
clear all
close all
load baseline.mat
LHSindices=[1 2];
LHS=LHS(:,LHSindices);
[results]=stlpm(R,shocksandstates(2).E,shocksandstates(2).M_Z_t,LHS,trend,Xlags,Rlags,runstate,includetrend,useNLshocks,H,CI,startS,endS,scaling,weight,cumweight)
stlpm_inference
block_bootstrap
save nber_results.mat

%% HP-filtered output gap
clear all
close all
load baseline.mat
LHSindices=[1 2];
LHS=LHS(:,LHSindices);
[results]=stlpm(R,shocksandstates(3).E,shocksandstates(3).M_Z_t,LHS,trend,Xlags,Rlags,runstate,includetrend,useNLshocks,H,CI,startS,endS,scaling,weight,cumweight)
stlpm_inference
block_bootstrap
save hpfilterog_results.mat

%% Centred state variable
clear all
close all
load baseline.mat
LHSindices=[1 2];
LHS=LHS(:,LHSindices);
[results]=stlpm(R,shocksandstates(4).E,shocksandstates(4).M_Z_t,LHS,trend,Xlags,Rlags,runstate,includetrend,useNLshocks,H,CI,startS,endS,scaling,weight,cumweight)
stlpm_inference
block_bootstrap
save centred_results.mat

%% c=50
clear all
close all
load baseline.mat
LHSindices=[1 2];
LHS=LHS(:,LHSindices);
[results]=stlpm(R,shocksandstates(4).E,shocksandstates(4).M_Z_t,LHS,trend,Xlags,Rlags,runstate,includetrend,useNLshocks,H,CI,startS,endS,scaling,weight,cumweight)
stlpm_inference
block_bootstrap
save c50_results.mat

%% theta=1
clear all
close all
load baseline.mat
LHSindices=[1 2];
LHS=LHS(:,LHSindices);
[results]=stlpm(R,shocksandstates(5).E,shocksandstates(5).M_Z_t,LHS,trend,Xlags,Rlags,runstate,includetrend,useNLshocks,H,CI,startS,endS,scaling,weight,cumweight)
stlpm_inference
block_bootstrap
save theta1_results.mat

%% theta=10
clear all
close all
load baseline.mat
LHSindices=[1 2];
LHS=LHS(:,LHSindices);
[results]=stlpm(R,shocksandstates(6).E,shocksandstates(6).M_Z_t,LHS,trend,Xlags,Rlags,runstate,includetrend,useNLshocks,H,CI,startS,endS,scaling,weight,cumweight)
stlpm_inference
block_bootstrap
save theta10_results.mat

%% No trend
clear all
close all
load baseline.mat
LHSindices=[1 2];
LHS=LHS(:,LHSindices);
includetrend=0;             % 1 includes a trend
[results]=stlpm(R,shocksandstates(1).E,shocksandstates(1).M_Z_t,LHS,trend,Xlags,Rlags,runstate,includetrend,useNLshocks,H,CI,startS,endS,scaling,weight,cumweight)
stlpm_inference
block_bootstrap
save notrend_results.mat

%% Calculate information criteria for lag length
clear all
close all
load baseline.mat
LHSindices=[1 2];
LHS=LHS(:,LHSindices);
Xlagsvec=[1 2 3 4];
Rlagsvec=[1 2 3 4];

% Loop over the lag length parameters
for a1=1:length(Xlagsvec)
    for a2=1:length(Rlagsvec)
        Xlags=Xlagsvec(a1);                    % Number of lags of control variables
        Rlags=Rlagsvec(a2);                    % Number of lags of policy variables

        [results]=stlpm(R,shocksandstates(1).E,shocksandstates(1).M_Z_t,LHS,trend,Xlags,Rlags,runstate,includetrend,useNLshocks,H,CI,startS,endS,scaling,weight,cumweight)
        
        for a3=1:size(results,2)-1
            AIC(a1,a2,a3)=sum(log(var(results(a3).NL.res))+2*(size(results(a3).NL.XNL,2))/size(results(a3).NL.XNL,1));
            SBC(a1,a2,a3)=sum((size(results(a3).NL.XNL,2))*log(size(results(a3).NL.XNL,1))+size(results(a3).NL.XNL,1)*log(var(results(a3).NL.res)));
        end
        a3=a3+1;
            AIC(a1,a2,a3)=sum(log(var(results(a3).NL.res))+2*(size(results(a3).NL.XRNL,2))/size(results(a3).NL.XRNL,1));
            SBC(a1,a2,a3)=sum((size(results(a3).NL.XRNL,2))*log(size(results(a3).NL.XRNL,1))+size(results(a3).NL.XRNL,1)*log(var(results(a3).NL.res)));
        
    end
end

AICsum=sum(AIC(:,:,1:end-1),3);
[M,I]=min(AICsum(:));
[AIC_row,AIC_col]=ind2sub(size(AICsum),I);
display(strcat('AIC criterion indicates ',num2str(Xlagsvec(AIC_row)),' lags of X and ',num2str(Rlagsvec(AIC_col)),' lags of R'))
% test=exp((min(min(AIC(:,:,1)))-AIC(:,:,1))/(2*21));

SBCsum=sum(SBC(:,:,1:end-1),3);
[M,I]=min(SBCsum(:));
[SBC_row,SBC_col]=ind2sub(size(SBCsum),I);
display(strcat('SBC criterion indicates ',num2str(Xlagsvec(SBC_row)),' lags of X and ',num2str(Rlagsvec(SBC_col)),' lags of R'))

save laglength_tests.mat

%% Run alternative lag lengths
clear all
close all
load baseline.mat
LHSindices=[1 2];
LHS=LHS(:,LHSindices);
Xlags=2;                    % Number of lags of control variables
[results]=stlpm(R,shocksandstates(1).E,shocksandstates(1).M_Z_t,LHS,trend,Xlags,Rlags,runstate,includetrend,useNLshocks,H,CI,startS,endS,scaling,weight,cumweight)
stlpm_inference
block_bootstrap
save xlags2_results.mat

%% Positive v negative
clear all
close all
load baseline.mat
useNLshocks=0;              % 1 uses non-linearly identified shocks
LHSindices=[1 2];
LHS=LHS(:,LHSindices);
[results]=stlpm2(R,shocksandstates(1).E,shocksandstates(1).M_Z_t,LHS,trend,Xlags,Rlags,includetrend,useNLshocks,H,CI,startS,endS,scaling,weight,cumweight)
stlpm_inference2
block_bootstrap2
save posvneg_results.mat

%% Large v small
clear all
close all
load baseline.mat
useNLshocks=0;              % 1 uses non-linearly identified shocks
LHSindices=[1 2];
LHS=LHS(:,LHSindices);
[results]=stlpm3(R,shocksandstates(1).E,shocksandstates(1).M_Z_t,LHS,trend,Xlags,Rlags,includetrend,useNLshocks,H,CI,startS,endS,scaling,weight,cumweight)
stlpm_inference3
block_bootstrap3
save largevsmall_results.mat

%% Tables and charts only
clear all
close all
load baseline_results.mat
standardcharts_fig2(results,[1 2 size(results,2)],'Baseline results, headline variables','Baseline results, headline variables',[names(LHSindices(1:2)) 'Fed Funds rate'],{'Log points' 'Log points' 'Percentage points'},H,CI,CIparam,{'smoothed','cum','smoothed'})
standardcharts_fig3(results,3:5,'Baseline results, expenditure volumes','Baseline results, expenditure volumes',names(LHSindices(3:5)),{'Log points' 'Log points' 'Log points'},H,CI,CIparam,{'smoothed','smoothed','smoothed'})
standardcharts5(results,9:12,'Baseline results, fiscal and credit variables','Baseline results, fiscal and credit variables',names(LHSindices(9:12)),{'Log points' 'Log points' 'Log points' 'Log points'},H,CI,CIparam,{'smoothed','smoothed','smoothed','smoothed'})
tables_aejm_randr
filename='table_baseline_aejm.csv';
csvwrite(filename,table,3,4);
% Unsmoothed baseline IRFs
standardcharts(results,[1 2 size(results,2)],'Baseline results, unsmoothed IRFs, headline variables','Baseline results, unsmoothed IRFs, headline variables',[names(LHSindices(1:2)) 'Fed Funds rate'],{'Log points' 'Log points' 'Percentage points'},H,CI,CIparam,{'scaled','scaled','scaled'})

% Linear identification
clear all
close all
load linearid_results.mat
standardcharts_fig8(results,[1 2 size(results,2)],'Linear identification, headline variables','Linear identification, headline variables',[names(LHSindices(1:2)) 'Fed Funds rate'],{'Log points' 'Log points' 'Percentage points'},H,CI,CIparam,{'smoothed','cum','smoothed'})
tables_aejm_randr
filename='table_linearID_aejm.csv';
csvwrite(filename,table,3,4);

% VAR identification
clear all
close all
load var_results.mat
standardcharts_fig9(results,[1 2 size(results,2)],'VAR identification, headline variables','VAR identification, headline variables',[names(LHSindices(1:2)) 'Fed Funds rate'],{'Log points' 'Log points' 'Percentage points'},H,CI,CIparam,{'smoothed','cum','smoothed'})
tables_aejm_randr
filename='table_VARID_aejm.csv';
csvwrite(filename,table,3,4);

% No trend
clear all
close all
load notrend_results.mat
standardcharts_fig10(results,[1 2 size(results,2)],'No trend, headline variables','No trend, headline variables',[names(LHSindices(1:2)) 'Fed Funds rate'],{'Log points' 'Log points' 'Percentage points'},H,CI,CIparam,{'smoothed','cum','smoothed'})
tables_aejm_randr
filename='table_notrend_aejm.csv';
csvwrite(filename,table,3,4);

% Run alternative lag lengths
clear all
close all
load xlags2_results.mat
standardcharts_fig11(results,[1 2 size(results,2)],'2 lags of dependent variable, headline variables','2 lags of dependent variable, headline variables',[names(LHSindices(1:2)) 'Fed Funds rate'],{'Log points' 'Log points' 'Percentage points'},H,CI,CIparam,{'smoothed','cum','smoothed'})
tables_aejm_randr
filename='table_2lags_aejm.csv';
csvwrite(filename,table,3,4);

% NBER recessions
clear all
close all
load nber_results.mat
standardcharts4(results,[1 2 size(results,2)],'NBER recessions, headline variables','NBER recessions, headline variables',[names(LHSindices(1:2)) 'Fed Funds rate'],{'Log points' 'Log points' 'Percentage points'},H,CI,CIparam,{'smoothed','cum','smoothed'})
tables_aejm_randr
filename='table_NBER_aejm.csv';
csvwrite(filename,table,3,4);

% HP-filtered output gap
clear all
close all
load hpfilterog_results.mat
standardcharts4(results,[1 2 size(results,2)],'HP filter, headline variables','HP filter, headline variables',[names(LHSindices(1:2)) 'Fed Funds rate'],{'Log points' 'Log points' 'Percentage points'},H,CI,CIparam,{'smoothed','cum','smoothed'})
tables_aejm_randr
filename='table_HPfilter_aejm.csv';
csvwrite(filename,table,3,4);

% Centred state variable
clear all
close all
load centred_results.mat
standardcharts(results,[1 2 size(results,2)],'Centered state variable, headline variables','Centered state variable, headline variables',[names(LHSindices(1:2)) 'Fed Funds rate'],{'Log points' 'Log points' 'Percentage points'},H,CI,CIparam,{'smoothed','cum','smoothed'})
tables_aejm_randr
filename='table_centered_aejm.csv';
csvwrite(filename,table,3,4);

% c=50
clear all
close all
load c50_results.mat
standardcharts(results,[1 2 size(results,2)],'Recession at 50th percentile, headline variables','Recession at 50th percentile, headline variables',[names(LHSindices(1:2)) 'Fed Funds rate'],{'Log points' 'Log points' 'Percentage points'},H,CI,CIparam,{'smoothed','cum','smoothed'})
tables_aejm_randr
filename='table_c50_aejm.csv';
csvwrite(filename,table,3,4);

% theta=1
clear all
close all
load theta1_results.mat
standardcharts(results,[1 2 size(results,2)],'Theta=1, headline variables','Theta=1, headline variables',[names(LHSindices(1:2)) 'Fed Funds rate'],{'Log points' 'Log points' 'Percentage points'},H,CI,CIparam,{'smoothed','cum','smoothed'})
tables_aejm_randr
filename='table_theta1_aejm.csv';
csvwrite(filename,table,3,4);

% theta=10
clear all
close all
load theta10_results.mat
standardcharts(results,[1 2 size(results,2)],'Theta=10, headline variables','Theta=10, headline variables',[names(LHSindices(1:2)) 'Fed Funds rate'],{'Log points' 'Log points' 'Percentage points'},H,CI,CIparam,{'smoothed','cum','smoothed'})
tables_aejm_randr
filename='table_theta10_aejm.csv';
csvwrite(filename,table,3,4);

% Positive v negative
clear all
close all
load posvneg_results.mat
[test]=load('baseline_results.mat','results')
results(1).lin=test.results(1).lin
results(2).lin=test.results(1).lin
results(end).lin=test.results(end).lin
standardcharts2(results,[1 2 size(results,2)],'Positive and negative shocks, headline variables','Positive and negative shocks, headline variables',[names(LHSindices(1:2)) 'Fed Funds rate'],{'Log points' 'Log points' 'Percentage points'},H,CI,CIparam,{'smoothed','cum','smoothed'})
tables_aejm_randr
filename='table_posvneg_aejm.csv';
csvwrite(filename,table,3,4);

% Large v small
clear all
close all
load largevsmall_results.mat
[test]=load('baseline_results.mat','results')
results(1).lin=test.results(1).lin
results(2).lin=test.results(1).lin
results(end).lin=test.results(end).lin
standardcharts3(results,[1 2 size(results,2)],'Large and small shocks, headline variables','Large and small shocks, headline variables',[names(LHSindices(1:2)) 'Fed Funds rate'],{'Log points' 'Log points' 'Percentage points'},H,CI,CIparam,{'smoothed','cum','scaled'})
tables_aejm_randr_3cols
filename='table_bigvsmall_aejm.csv';
csvwrite(filename,table,3,4);


%% Chart 7
clear all
close all
load('baseline_results.mat')
figure('Name','Chart 7: regime-specific pdf and cdfs of shocks')
subplot(2,1,1)
[boomdens,xi]=ksdensity(shocksandstates(1).E(:,1),'weights',shocksandstates(1).M_Z_t,'function','pdf');
[recdens,yi]=ksdensity(shocksandstates(1).E(:,1),'weights',1-shocksandstates(1).M_Z_t,'function','pdf');
[avgdens,zi]=ksdensity(shocksandstates(1).E(:,1),'function','pdf');
plot(xi,boomdens,'--g',yi,recdens,'r-.',zi,avgdens,'-k', 'LineWidth', 1.5);
handle=legend('expansion','recession','average');
legend(handle,'boxoff');
box on
set(gcf,'Color', [1 1 1])
title('PDF', 'FontName', 'Arial Black')

subplot(2,1,2)
[boomdens,xi]=ksdensity(shocksandstates(1).E(:,1),'weights',shocksandstates(1).M_Z_t,'function','cdf');
[recdens,yi]=ksdensity(shocksandstates(1).E(:,1),'weights',1-shocksandstates(1).M_Z_t,'function','cdf');
[avgdens,zi]=ksdensity(shocksandstates(1).E(:,1),'function','cdf');
plot(xi,boomdens,'--g',yi,recdens,'r-.',zi,avgdens,'-k','LineWidth', 1.5);
hold on
plot(xi,0.5*ones(size(xi)),'-k');
handle=legend('expansion','recession','average');
legend(handle,'boxoff');
box on
title('CDF','FontName', 'Arial Black')
hold off
saveas(gcf,'chart 6 regime-specific pdfs and cdfs','pdf');

%% Chart 1
% Import the original Romer shocks
[~, ~, raw] = xlsread('RomerandRomerDataAppendix.xlsx','DATA BY MONTH','T16:T127');
Romer_original= reshape([raw{:}],size(raw));
dates=linspace(1969,2002.75,136);
threeshocks=[[[Romer_original; NaN(length(shocksandstates(1).E(:,1))-length(Romer_original),1)] shocksandstates(1).E(:,1) shocksandstates(1).E(:,2)]]
figure('name','Chart 1: monetary policy shocks and state variable')
[AX,H1,H2]=plotyy(dates,threeshocks,dates,shocksandstates(1).M_Z_t)
set(AX(2),'YLim',[-0.5 1.5])
handle=legend('Romer shocks','Extended linear shocks','Nonlinearly identified shocks','Probability of an expansion');
legend(handle,'boxoff');
box on
set(gcf,'Color', [1 1 1])
hold off
saveas(gcf,'chart 1 mp shocks and state variable merged','pdf');
