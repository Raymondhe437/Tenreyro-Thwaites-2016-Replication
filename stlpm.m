function [results]=stlpm(R,E,M,LHS,trend,Xlags,Rlags,runstate,includetrend,useNLshocks,H,CI,startS,endS,scaling,weight,cumweight)

%hsk% Spits out a structure "results" with 13 sub-structures inside.
% 13 substructures = 12 explanatory variables + 1 policy rates (ffr)
% footnote 4 says "In the baseline specification, X_t contains one lag each of the dependent variable and federal funds rate."
% But when ffr itself is dependent variable, not X_t just contains one lag of ffr. - This is done first. 

%E = shocksandstates(1).E;
%M = shocksandstates(1).M_Z_t;



% estimates an STLPM as a function

% set up regression matrices
N=size(LHS,2);              % Number of endogenous variables = 12
XR=[];
if Rlags>0
    XR=R(startS-1:endS-1); 
    if Rlags>1
        for l=2:Rlags
            XR=[XR R(startS-l:endS-l)]; %   lagged interest rates                                                                
        end
    end
end

%hsk% R is interest rates. XR is a subset of interest rates for the sample period.
%hsk% E is two series of shocks. Linear and Non-linear,

XR=[ones(length(XR),1) XR E(:,1)];   % construct common RHS projection matrix using linear shock.
if includetrend==1
    XR=[trend XR];
end
%hsk% XR consists of (trend, intercept, interest rates, linear shocks)

XRproj=inv(XR'*XR)*XR'; %hsk% not exactly projection matrix. But used to calculate coefficient later
KR=size(XR,2); %hsk% Number of variables in XRproj

%hsk% Setting up non-linear projection matrix.
if runstate==1;

    if useNLshocks==0
        XRNL=[repmat(M,[1 size(XR(:,1+includetrend:end),2)]).*XR(:,1+includetrend:end) (1-repmat(M,[1 size(XR(:,1+includetrend:end),2)])).*XR(:,1+includetrend:end)];
    else %hsk% defualt is NLshocks == 1, M is probability of expansion for each FOMC meeting. This is just [F(z_t)X_t, (1-F(z_{t}))X_{t}] but X_{t} excludes time trend. 
        XRNL=[repmat(M,[1 size(XR(:,1+includetrend:end),2)]).*[XR(:,1+includetrend:end-1) E(:,2)] (1-repmat(M,[1 size(XR(:,1+includetrend:end),2)])).*[XR(:,1+includetrend:end-1) E(:,2)]];
    end
    if includetrend==1
        XRNL=[trend XRNL];
    end
    
    XRNLproj=inv(XRNL'*XRNL)*XRNL';  
    KRNL=size(XRNL,2);
    
end

%% Main regression loop

% loop over H, store coefficients and residuals
% policy variable regressions 
%hsk% dependent variable is ffr
%daf% So one set of regressions is run for every set of lags h. However,
%daf% The other set of regressions is only performed when runstate is 1.

betaR=[];
resR=[];
betafull=[];
betafullNL=[];
betaRNL=zeros(H+1,2);
resRNL=[];
for h=0:H
    betaRh=XRproj*R(startS+h:endS+h); %teporary storage of beta_{t+h}%
    results(N+1).lin.betafull=betaRh;
    betaR=[betaR; betaRh(end)]; %store response to the shock%
    betafull(:,h+1)=betaRh; %need to add 1 because there is no column with index 0%
    resRh=R(startS+h:endS+h)-XR*betaRh; %extract residual - maybe used for inference.%
    resR(:,h+1)=resRh;
    results(N+1).lin.variance(:,:,h+1)=var(resR(:,h+1))*inv(XR'*XR);
    
    if runstate==1
        betaRhNL=XRNLproj*R(startS+h:endS+h);
        betafullNL(:,h+1)=betaRhNL;
        betaRNL(h+1,1)=betaRhNL(1+(end-includetrend)/2);
        betaRNL(h+1,2)=betaRhNL(end);
        resRhNL=R(startS+h:endS+h)-XRNL*betaRhNL;
        resRNL(:,h+1)=resRhNL;
        results(N+1).NL.variance(:,:,h+1)=var(resRNL(:,h+1))*inv(XRNL'*XRNL);
    end

%daf% This is getting the fitted estimates and residuals across both states

end
results(N+1).lin.res=resR;
results(N+1).lin.beta=betaR;
results(N+1).lin.betafull=betafull;
results(N+1).lin.XR=XR;
if runstate==1
    results(N+1).NL.res=resRNL;
    results(N+1).NL.beta=betaRNL;
    results(N+1).NL.betafullNL=betafullNL;
    results(N+1).NL.XRNL=XRNL;
end

%daf% This is only run if scaling is chosen. Then the coefficients are weighted by
%daf% first period impact of Fed Funds, which seems to be called "weight".

if scaling==1
    results(N+1).lin.scaled.coef=betaR/betaR(1); % Scale by first period impact on Fed Funds
    results(N+1).lin.smoothed.coef=weight*betaR/betaR(1); % Scale by first period impact on Fed Funds
    results(N+1).lin.cum.coef=cumweight*betaR/betaR(1);

    if runstate==1
        results(N+1).NL.scaled.coef(:,1)=betaRNL(:,1)/betaRNL(1,1);
        results(N+1).NL.scaled.coef(:,2)=betaRNL(:,2)/betaRNL(1,2);
        results(N+1).NL.scaled.coef_diff=results(N+1).NL.scaled.coef(:,1)-results(N+1).NL.scaled.coef(:,2);

        results(N+1).NL.smoothed.coef(:,1)=weight*results(N+1).NL.scaled.coef(:,1);
        results(N+1).NL.smoothed.coef(:,2)=weight*results(N+1).NL.scaled.coef(:,2);
        results(N+1).NL.smoothed.coef_diff=results(N+1).NL.smoothed.coef(:,1)-results(N+1).NL.smoothed.coef(:,2);

        results(N+1).NL.cum.coef(:,1)=cumweight*results(N+1).NL.scaled.coef(:,1);
        results(N+1).NL.cum.coef(:,2)=cumweight*results(N+1).NL.scaled.coef(:,2);        
        results(N+1).NL.cum.coef_diff=results(N+1).NL.cum.coef(:,1)-results(N+1).NL.cum.coef(:,2);       
    end
end

% response variable regressions
%hsk% 12 other dependent variables. The only difference from above is there is one more explanatory variable, ffr.
for n=1:N; % loop over LHS variables
    %n = 1; %hsk% GDP volume
    results(N+1).lin.beta=[];
    res=[];
    X=[];
    beta=[];
    betafull=[];
    betafullNL=[];
    
    if Xlags>0 %hsk% default is 1.
        X=LHS(startS-1:endS-1,n); %hsk% One lag of each dependenet variable (see footnote 4)
        if Xlags>1
            for l=2:Xlags
                X=[X LHS(startS-l:endS-l,n)];
            end
        end
    end
%    X2=X;
    if includetrend==0
        X=[X XR]; 
    else %hsk% default is includetrend == 1
        X=[trend X XR(:,2:end)]; %hsk% X consists of (trend, lagged GDP volume (if n =1), interest rates, linear shocks)
    end
    results(n).lin.X=X;
    Xproj=inv(X'*X)*X';
    K=size(X,2);
    linshockindex=K;

    results(n).NL.beta=[];
    resNL=[];
    XNL=[];
    betaNL=[];
    betaNL=zeros(H+1,2);
    %daf% if a trend is included, then it is made the initial column of XNL, with
    %daf% apparently some repositioning of the other (repeated) columns to drop
    %daf% their first columns with "2:end", which drops the 1st col?
    if runstate==1;
%        XNL=[M.*X2 (1-M).*X2 XRNL];
        if includetrend==0
            XNL=[repmat(M,[1 size(X,2)]).*X (1-repmat(M,[1 size(X,2)])).*X];
        else
            XNL=[trend repmat(M,[1 size(X(:,2:end),2)]).*X(:,2:end) (1-repmat(M,[1 size(X(:,2:end),2)])).*X(:,2:end)];
        end
        results(n).NL.XNL=XNL;
        XNLproj=inv(XNL'*XNL)*XNL';
    end 
    KNL=size(XNL,2);
    recshockindex=KNL;
    boomshockindex=KNL-(KNL-includetrend)/2; %hsk% Index that indicates the shock column?
    
    % loop over H, store coefficients and residuals
    for h=0:H
        betah=Xproj*LHS(startS+h:endS+h,n); %hsk% (n=1) Regress GDP Volume on (trend, lagged GDP volume (if n =1), interest rates, linear shocks).
        betafull(:,h+1)=betah; %hsk% need to +1 because h = 0 cannot be used as a column index.
        beta=[beta; betah(end)]; %hsk% store response to shock
        resh=LHS(startS+h:endS+h,n)-X*betah;
        res(:,h+1)=resh;
%        results(n).lin.variance(:,:,h+1)=var(resR(:,h+1))*inv(XR'*XR);

        if runstate==1;
            betahNL=XNLproj*LHS(startS+h:endS+h,n);
            betafullNL(:,h+1)=betahNL;
            betaNL(h+1,1)=betahNL(KNL-(KNL-includetrend)/2);
            betaNL(h+1,2)=betahNL(KNL);
            reshNL=LHS(startS+h:endS+h,n)-XNL*betahNL;
            resNL(:,h+1)=reshNL;
%            results(n).NL.variance(:,:,h+1)=var(resRNL(:,h+1))*inv(XRNL'*XRNL);

        end
    end
    
    results(n).lin.beta=beta;
    results(n).lin.res=res;
    results(n).lin.betafull=betafull;
    
    if runstate==1;
        results(n).NL.betaNL=betaNL;
        results(n).NL.res=resNL;            
        results(n).NL.betafullNL=betafullNL;
    end
    
    if scaling==1
    
    % Construct smoothed and cumulative scaled IRFs
    % matrix of weights for the moving average process
    % Linear model MAs and CIs
    %daf% It looks like this is a preparation of the creation of CIs (coef difference)
        results(n).lin.scaled.coef=results(n).lin.beta/betaR(1); % Scale by first period impact on Fed Funds
        results(n).lin.smoothed.coef=weight*results(n).lin.beta/betaR(1); % Scale by first period impact on Fed Funds
        results(n).lin.cum.coef=cumweight*results(n).lin.beta/betaR(1);

        if runstate==1;
            results(n).NL.scaled.coef(:,1)=results(n).NL.betaNL(:,1)/betaRNL(1);
            results(n).NL.scaled.coef(:,2)=results(n).NL.betaNL(:,2)/betaRNL(2);
            results(n).NL.scaled.coef_diff=results(n).NL.scaled.coef(:,1)-results(n).NL.scaled.coef(:,2);
            
            %daf% smoothed is just scaled times weight
            
            results(n).NL.smoothed.coef(:,1)=weight*results(n).NL.scaled.coef(:,1);
            results(n).NL.smoothed.coef(:,2)=weight*results(n).NL.scaled.coef(:,2);
            results(n).NL.smoothed.coef_diff=results(n).NL.smoothed.coef(:,1)-results(n).NL.smoothed.coef(:,2);
            
            %daf% cum is just scaled times cumulative weight

            results(n).NL.cum.coef(:,1)=cumweight*results(n).NL.scaled.coef(:,1);
            results(n).NL.cum.coef(:,2)=cumweight*results(n).NL.scaled.coef(:,2);        
            results(n).NL.cum.coef_diff=results(n).NL.cum.coef(:,1)-results(n).NL.cum.coef(:,2);       
        end
    end
    
end
%display('Estimation complete')

end

