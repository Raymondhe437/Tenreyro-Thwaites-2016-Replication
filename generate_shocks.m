%% This script generates Romer&Romer shock (both linearly and non-linearly) 
%% Thus, there is nothing to do with the authors' datasets. This script only uses Data >> "romer" data.
%% Generate the first day in the subsequent quarter then subtract one day from it
%% Let's see if this shows up on Ray's github!!!
%% hmmm need to view and commit changes to push again.
%% Otherwise, the push button is not clickable.

%hsk% Data construction need to be done beforehand
settings.sample0 = startS:1:endS; 
%hsk% create field "sample0" under the structure "settings" and store
%values from startS to endS.

%hsk% recall that our dates are initially expressed as "1947.0", "1947.25", "1947.5","1947.75"  
datestring=datenum(floor(dates(startS:endS)),3+12*(dates(startS:endS)-floor(dates(startS:endS))),1);
%hsk% floor(number) returns the integer part of the number.
%3+12*(dates(startS:endS)-floor(dates(startS:endS))) generates the subsequent month
%1 indicates the first dat.
%datestring is an array of dates in general number format.


for i=1:length(datestring)
    datestring(i)=addtodate(addtodate(datestring(i),1,'month'),-1,'day');
end
%hsk% addtodate(addtodate(datestring(i),1,'month'),-1,'day') adds one month
%(approximately 31 days) to each cell and then subtracts 1 day.

test=datestr(datestring);
%hsk% Transform (MATLAB serial) numbers into the usual date format  (DD-MM-YYYY). 
%hsk% The date ranges from 31-Mar-1969 to 31-Dec-2002 (Quarterly)
%hsk% The date length is equivalent to Data > "Romer"
%hsk% test variable is used in the main script.
%% Import the data, extracting spreadsheet dates in MATLAB serial date number format (datenum)
[~, ~, raw, dateNums] = xlsread('C:\Users\HUSANG KIM\Desktop\Papers\3. MPU&MPP\6. Local_Projection_MATLAB\data.xlsx','romer','A2:U352','',@convertSpreadsheetDates);
%hsk% The only difference b/w raw and dateNums is how the date column (from Data > Romer) is formatted
%hsk% The ~ represents the discarded outputs from function xlsread.

% Replace date strings by MATLAB serial date numbers (datenum)
R1 = ~cellfun(@isequalwithequalnans,dateNums,raw) & cellfun('isclass',raw,'char'); 
%hsk% Find spreadsheet dates by finding cells that are NOT equal between dateNums and raw. The only difference between dateNumbs and raw was the format of dates. So R1 is a logical array that consists of 1, for the dates and 0 for other values.
raw(R1) = dateNums(R1);
%hsk% Replace the dates in raw with dates in dateNums.

%format longG
% Create output variable
RomerandRomer20041 = reshape([raw{:}],size(raw));
%hsk% I am not sure what format reshape() leads to. But RomerandRomer20041
%is a matrix that contains data from sheet "romer".

% Clear temporary variables
clearvars raw dateNums;

%% Create the D matrix.
% now create a matrix to get M into meeting space
Qdate=repmat(datestring,[1 size(RomerandRomer20041,1)]);
%hsk% size(RomerandRomer20041,1) returns the number of 'rows' of RomerandRomer20041.
%hsk% Creates 351 repeat copies of datestring.
%hsk% repmat(datestring,1,351) creates 351 copies stacked horizontally (351 columns).
%hsk% repmat(datestring,351) creates 351 copies stacked vertically.

Qdatelag=repmat([0; datestring(1:end-1)],[1 size(RomerandRomer20041,1)]);
%hsk% create repeat copies of "lagged" dates (lag = 1). The initial empty date is set to 0

meetdate=repmat(RomerandRomer20041(:,1)',[length(datestring) 1]);
%hsk% RomerandRomer20041 dates: 719226 ~ 733387 (351 obs)
%hsk% datestring dates: 719253 ~ 731581 (136 obs)
%hsk% meetdate: 136 rows of 719226 ~ 733387 (351 obs)
D=(Qdate>=meetdate & Qdatelag<meetdate);    % D is Q x m = 136 x 351
%hsk% 136:31-Mar-1969 to 31-Dec-2002 (Quarterly) - This is the sample period the paper uses.
%hsk% 351: 4-Mar-1969 to 11-Dec-2007 (FOMC meetings monthly or bimonthly)
%hsk% D is an indicator matrix that tells whether a FOMC meeting is in between two quarters.
%% regime variable - output growth
i=1;
shocksandstates(i).Zlag=0; % how many lags for the state variable - AG use 1.  Fully lagged would be maband+1/8
shocksandstates(i).maband=7; % how many periods over which moving average is taken for the state variable - AG use 7
shocksandstates(i).percentile=20; % percentile of state variable cutoff
shocksandstates(i).theta=3; % intensity of regime switching

shocksandstates(i).statename='GDP growth';

pma = filter((1/shocksandstates(i).maband)*ones(1,shocksandstates(i).maband),1,[zeros(1,1); diff(LHS(:,1),1)],[]); %Trailing moving average
%hsk% I wonder what the last [] is doing?


%filter function example%
%x = [1,1,1,1,1,1,1,1,1,1];
%windowSize = 3; 
%b = (1/windowSize)*ones(1,windowSize);
%a = 1;
%y = filter(b,a,x,[]);
%x;
%y;

Z=pma(startS:endS); %hsk% extracts 136 observations.
Z0=Z-prctile(Z,shocksandstates(i).percentile); % Set zero at the percentile we specified
%hsk% Find the lower 20% value of state variable Z

Z0=Z0/std(Z0); % Give the variable unit standard deviation


shocksandstates(i).M_Z_t=exp(shocksandstates(i).theta*Z0(:,1))./(1+exp(shocksandstates(i).theta*Z0(:,1))); 
% Logistic function of the state variable to use in regression

%hsk% the code below shows the logistic function
%plot(sort(Z0),sort(shocksandstates(i).M_Z_t))


shocksandstates(i).M=D'*shocksandstates(i).M_Z_t;
% M is the state of the world for every FOMC meeting
% M_Z_t is the state of the world for every quarter.
x = [ones(351,1) RomerandRomer20041(:,3:end-1)];
%hsk% end-1 because end is column u which contains residual. (see sheet "romer", column U.)

y = RomerandRomer20041(:,2); %hsk% y variable is changes in federal funds rate
p=x*inv(x'*x)*x'; %hsk% projection matrix
m=eye(size(p))-p;
r=m*y; % these are the linearly identified shocks %hsk% Obtaining residual after regressing y onto p
r_Q=D*r;
%multiplication by D reduces the length from 351 to 136. Something to do with the sample period%


% Identify the shocks using a nonlinear analogue of the R&R equation
Mrep=repmat(shocksandstates(i).M,[1 size(x,2)]);
xnl=[Mrep.*x (1-Mrep).*x]; %This is just [F(z_t)X_t, (1-F(z_{t}))X_{t}] in eq (3)%
pnl=xnl*inv(xnl'*xnl)*xnl';
mnl=eye(size(pnl))-pnl;
rnl=mnl*y; % these are the nonlinearly identified shocks
rnl_Q=D*rnl;

%hsk% 19 explanatory variables * 2 states = 38 beta coefficients.
%beta = inv(xnl'*xnl)*xnl'*y;

shocksandstates(i).E=[r_Q rnl_Q]; %Y(settings.startS:settings.endS,shockindex+settings.endoravg);

%% regime variable - nber dates 
i=2;
shocksandstates(i).Zlag=0; % how many lags for the state variable - AG use 1.  Fully lagged would be maband+1/8
shocksandstates(i).maband=7; % how many periods over which moving average is taken for the state variable - AG use 7
shocksandstates(i).percentile=20; % percentile of state variable cutoff
shocksandstates(i).theta=3; % intensity of regime switching

shocksandstates(i).statename='NBER dates'
pma=filter((1/shocksandstates(i).maband)*ones(1,shocksandstates(i).maband),1,alt_states(:,3),[]); %Trailing moving average
Z=pma(startS:endS);

Z=1-Z;
Z0=Z-prctile(Z,shocksandstates(i).percentile); % Set zero at the percentile we specified
Z0=Z0/std(Z0); % Give the variable unit standard deviation
shocksandstates(i).M_Z_t=exp(shocksandstates(i).theta*Z0(:,1))./(1+exp(shocksandstates(i).theta*Z0(:,1))); % Logistic function of the state variable to use in regression

shocksandstates(i).M=D'*shocksandstates(i).M_Z_t;
x=[ones(351,1) RomerandRomer20041(:,3:end-1)];
y=RomerandRomer20041(:,2);
p=x*inv(x'*x)*x';
m=eye(size(p))-p;
r=m*y; % these are the linearly identified shocks
r_Q=D*r;

Mrep=repmat(shocksandstates(i).M,[1 size(x,2)]);
xnl=[Mrep.*x (1-Mrep).*x];
pnl=xnl*inv(xnl'*xnl)*xnl';
mnl=eye(size(pnl))-pnl;
rnl=mnl*y; % these are the nonlinearly identified shocks
rnl_Q=D*rnl;

shocksandstates(i).E=[r_Q rnl_Q]; %Y(settings.startS:settings.endS,shockindex+settings.endoravg);

%% regime variable - HP filtered output gap
i=3;
shocksandstates(i).Zlag=0; % how many lags for the state variable - AG use 1.  Fully lagged would be maband+1/8
shocksandstates(i).maband=7; % how many periods over which moving average is taken for the state variable - AG use 7
shocksandstates(i).percentile=20; % percentile of state variable cutoff
shocksandstates(i).theta=3; % intensity of regime switching

shocksandstates(i).statename='HP filter'
pma=filter((1/shocksandstates(i).maband)*ones(1,shocksandstates(i).maband),1,alt_states(:,1),[]); %Trailing moving average
Z=pma(startS:endS);
Z0=Z-prctile(Z,shocksandstates(i).percentile); % Set zero at the percentile we specified
Z0=Z0/std(Z0); % Give the variable unit standard deviation
shocksandstates(i).M_Z_t=exp(shocksandstates(i).theta*Z0(:,1))./(1+exp(shocksandstates(i).theta*Z0(:,1))); % Logistic function of the state variable to use in regression

shocksandstates(i).M=D'*shocksandstates(i).M_Z_t;
x=[ones(351,1) RomerandRomer20041(:,3:end-1)];
y=RomerandRomer20041(:,2);
p=x*inv(x'*x)*x';
m=eye(size(p))-p;
r=m*y; % these are the linearly identified shocks
r_Q=D*r;

Mrep=repmat(shocksandstates(i).M,[1 size(x,2)]);
xnl=[Mrep.*x (1-Mrep).*x];
pnl=xnl*inv(xnl'*xnl)*xnl';
mnl=eye(size(pnl))-pnl;
rnl=mnl*y; % these are the nonlinearly identified shocks
rnl_Q=D*rnl;

% shocks
% Create a (possibly lagged) matrix of the Romer shocks
shocksandstates(i).E=[r_Q rnl_Q]; %Y(settings.startS:settings.endS,shockindex+settings.endoravg);

%% regime variable - phase shift - centred
i=4; % 
shocksandstates(i).Zlag=0; % how many lags for the state variable - AG use 1.  Fully lagged would be maband+1/8
shocksandstates(i).maband=7; % how many periods over which moving average is taken for the state variable - AG use 7
shocksandstates(i).percentile=20; % percentile of state variable cutoff
shocksandstates(i).theta=3; % intensity of regime switching

shocksandstates(i).statename='GDP growth, centred'
pma=filter((1/shocksandstates(i).maband)*ones(1,shocksandstates(i).maband),1,[zeros(1,1); diff(LHS(:,1),1)],[]); %Trailing moving average
Z=pma(startS+4:endS+4);
Z0=Z-prctile(Z,shocksandstates(i).percentile); % Set zero at the percentile we specified
Z0=Z0/std(Z0); % Give the variable unit standard deviation
shocksandstates(i).M_Z_t=exp(shocksandstates(i).theta*Z0(:,1))./(1+exp(shocksandstates(i).theta*Z0(:,1))); % Logistic function of the state variable to use in regression

shocksandstates(i).M=D'*shocksandstates(i).M_Z_t;
x=[ones(351,1) RomerandRomer20041(:,3:end-1)];
y=RomerandRomer20041(:,2);
p=x*inv(x'*x)*x';
m=eye(size(p))-p;
r=m*y; % these are the linearly identified shocks
r_Q=D*r;

Mrep=repmat(shocksandstates(i).M,[1 size(x,2)]);
xnl=[Mrep.*x (1-Mrep).*x];
pnl=xnl*inv(xnl'*xnl)*xnl';
mnl=eye(size(pnl))-pnl;
rnl=mnl*y; % these are the nonlinearly identified shocks
rnl_Q=D*rnl;

% shocks
% Create a (possibly lagged) matrix of the Romer shocks
shocksandstates(i).E=[r_Q rnl_Q]; %Y(settings.startS:settings.endS,shockindex+settings.endoravg);

%% regime variable - theta=1
i=5; % 
shocksandstates(i).Zlag=0; % how many lags for the state variable - AG use 1.  Fully lagged would be maband+1/8
shocksandstates(i).maband=7; % how many periods over which moving average is taken for the state variable - AG use 7
shocksandstates(i).percentile=20; % percentile of state variable cutoff
shocksandstates(i).theta=1; % intensity of regime switching

shocksandstates(i).statename='GDP growth, theta=1'
pma=filter((1/shocksandstates(i).maband)*ones(1,shocksandstates(i).maband),1,[zeros(1,1); diff(LHS(:,1),1)],[]); %Trailing moving average
Z=pma(startS:endS);
Z0=Z-prctile(Z,shocksandstates(i).percentile); % Set zero at the percentile we specified
Z0=Z0/std(Z0); % Give the variable unit standard deviation
shocksandstates(i).M_Z_t=exp(shocksandstates(i).theta*Z0(:,1))./(1+exp(shocksandstates(i).theta*Z0(:,1))); % Logistic function of the state variable to use in regression

shocksandstates(i).M=D'*shocksandstates(i).M_Z_t;
x=[ones(351,1) RomerandRomer20041(:,3:end-1)];
y=RomerandRomer20041(:,2);
p=x*inv(x'*x)*x';
m=eye(size(p))-p;
r=m*y; % these are the linearly identified shocks
r_Q=D*r;

Mrep=repmat(shocksandstates(i).M,[1 size(x,2)]);
xnl=[Mrep.*x (1-Mrep).*x];
pnl=xnl*inv(xnl'*xnl)*xnl';
mnl=eye(size(pnl))-pnl;
rnl=mnl*y; % these are the nonlinearly identified shocks
rnl_Q=D*rnl;

% shocks
% Create a (possibly lagged) matrix of the Romer shocks
shocksandstates(i).E=[r_Q rnl_Q]; %Y(settings.startS:settings.endS,shockindex+settings.endoravg);

%% regime variable - theta=10
i=6; % theta=10
shocksandstates(i).Zlag=0; % how many lags for the state variable - AG use 1.  Fully lagged would be maband+1/8
shocksandstates(i).maband=7; % how many periods over which moving average is taken for the state variable - AG use 7
shocksandstates(i).percentile=20; % percentile of state variable cutoff
shocksandstates(i).theta=10; % intensity of regime switching

shocksandstates(i).statename='GDP growth, theta=10'
pma=filter((1/shocksandstates(i).maband)*ones(1,shocksandstates(i).maband),1,[zeros(1,1); diff(LHS(:,1),1)],[]); %Trailing moving average

Z=pma(startS:endS);
Z0=Z-prctile(Z,shocksandstates(i).percentile); % Set zero at the percentile we specified
Z0=Z0/std(Z0); % Give the variable unit standard deviation
shocksandstates(i).M_Z_t=exp(shocksandstates(i).theta*Z0(:,1))./(1+exp(shocksandstates(i).theta*Z0(:,1))); % Logistic function of the state variable to use in regression

shocksandstates(i).M=D'*shocksandstates(i).M_Z_t;
x=[ones(351,1) RomerandRomer20041(:,3:end-1)];
y=RomerandRomer20041(:,2);
p=x*inv(x'*x)*x';
m=eye(size(p))-p;
r=m*y; % these are the linearly identified shocks
r_Q=D*r;

Mrep=repmat(shocksandstates(i).M,[1 size(x,2)]);
xnl=[Mrep.*x (1-Mrep).*x];
pnl=xnl*inv(xnl'*xnl)*xnl';
mnl=eye(size(pnl))-pnl;
rnl=mnl*y; % these are the nonlinearly identified shocks
rnl_Q=D*rnl;

% shocks
% Create a (possibly lagged) matrix of the Romer shocks
shocksandstates(i).E=[r_Q rnl_Q]; %Y(settings.startS:settings.endS,shockindex+settings.endoravg);

%% regime variable - proportion of sample in recession is 50%
i=7; % 
shocksandstates(i).Zlag=0; % how many lags for the state variable - AG use 1.  Fully lagged would be maband+1/8
shocksandstates(i).maband=7; % how many periods over which moving average is taken for the state variable - AG use 7
shocksandstates(i).percentile=50; % percentile of state variable cutoff
shocksandstates(i).theta=3; % intensity of regime switching

shocksandstates(i).statename='GDP growth, recession percentile 50';
pma=filter((1/shocksandstates(i).maband)*ones(1,shocksandstates(i).maband),1,[zeros(1,1); diff(LHS(:,1),1)],[]); %Trailing moving average
Z=pma(startS:endS);
Z0=Z-prctile(Z,shocksandstates(i).percentile); % Set zero at the percentile we specified
Z0=Z0/std(Z0); % Give the variable unit standard deviation
shocksandstates(i).M_Z_t=exp(shocksandstates(i).theta*Z0(:,1))./(1+exp(shocksandstates(i).theta*Z0(:,1))); % Logistic function of the state variable to use in regression

shocksandstates(i).M=D'*shocksandstates(i).M_Z_t;
x=[ones(351,1) RomerandRomer20041(:,3:end-1)];
y=RomerandRomer20041(:,2);
p=x*inv(x'*x)*x';
m=eye(size(p))-p;
r=m*y; % these are the linearly identified shocks
r_Q=D*r;

Mrep=repmat(shocksandstates(i).M,[1 size(x,2)]);
xnl=[Mrep.*x (1-Mrep).*x];
pnl=xnl*inv(xnl'*xnl)*xnl';
mnl=eye(size(pnl))-pnl;
rnl=mnl*y; % these are the nonlinearly identified shocks
rnl_Q=D*rnl;

% shocks
% Create a (possibly lagged) matrix of the Romer shocks
shocksandstates(i).E=[r_Q rnl_Q]; %Y(settings.startS:settings.endS,shockindex+settings.endoravg);

%% regime variable - var shocks
i=8; % 
shocksandstates(i).statename='GDP growth, VAR shocks';
y=[LHS(:,1:2) R];
x=[ones(length(y(startS:endS,:)),1) y(startS-1:endS-1,:) y(startS-2:endS-2,:) y(startS-3:endS-3,:) y(startS-4:endS-4,:)];
p=x*inv(x'*x)*x';
m=eye(size(p))-p;
e=m*y(startS:endS,:); % these are the linearly identified shocks
A=chol(cov(e),'Lower');
u=((A^-1)*e')';

Mrep=repmat(shocksandstates(1).M_Z_t,[1 size(x,2)]);
xnl=[Mrep.*x (1-Mrep).*x];
pnl=xnl*inv(xnl'*xnl)*xnl';
mnl=eye(size(pnl))-pnl;
enl=mnl*y(startS:endS,:); % these are the linearly identified shocks

for j=1:size(enl,1)
    sigma_t(:,:,j)=enl(j,:)'*enl(j,:);
end

A_b=zeros(size(sigma_t(:,:,1)));
A_r=A_b;

rhs=[shocksandstates(1).M_Z_t 1-shocksandstates(1).M_Z_t];
b=inv(rhs'*rhs)*rhs';
for k=1:size(enl,2)
    for j=1:size(enl,2)
        beta=b*squeeze(sigma_t(k,j,:));
        A_b(k,j)=beta(1);
        A_r(k,j)=beta(2);
    end
end

% now do Choleski orthogonalisation
B_b=chol(A_b,'Lower');
B_r=chol(A_r,'Lower');
for t=1:length(shocksandstates(1).M_Z_t)
A_t(:,:,t)=(shocksandstates(1).M_Z_t(t)*B_b+(1-shocksandstates(1).M_Z_t(t))*B_r);
u_nl(:,t)=inv(A_t(:,:,t))*(enl(t,:)');
end

% shocks
% Create a (possibly lagged) matrix of the Romer shocks
shocksandstates(i).E=[u(:,3) u_nl(3,:)']; %Y(settings.startS:settings.endS,shockindex+settings.endoravg);
%% save the result
save ('data.mat','shocksandstates','-append')