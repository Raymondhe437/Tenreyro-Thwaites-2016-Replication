function [results]=stlpm2(R,E,M,LHS,trend,Xlags,Rlags,includetrend,useNLshocks,H,CI,startS,endS,scaling,weight,cumweight)
% estimates an STLPM as a function

% set up regression matrices
N=size(LHS,2);              % Number of endogenous variables
XR=[];
if Rlags>0
    XR=R(startS-1:endS-1);
    if Rlags>1
        for l=2:Rlags
            XR=[XR R(startS-l:endS-l)]; %   lagged interest rates                                                                
        end
    end
end

XR=[ones(length(XR),1) XR E(:,1) E(:,1).^3];   % construct common RHS projection matrix
if includetrend==1
    XR=[trend XR];
end
XRproj=inv(XR'*XR)*XR';
KR=size(XR,2);


%% Main regression loop
% loop over H, store coefficients and residuals

% policy variable regressions
betaR=[];
resR=[];
betafull=[];
betafullNL=[];
betaRNL=zeros(H+1,2);
resRNL=[];
for h=0:H

    betaRh=XRproj*R(startS+h:endS+h);
    results(N+1).lin.betafull=betaRh;
    betaR=[betaR; betaRh(end-1:end)];
    betafull(h+1,:)=betaRh;
    resRh=R(startS+h:endS+h)-XR*betaRh;
    resR(:,h+1)=resRh;
    
end
results(N+1).NL.res=resR;
results(N+1).NL.beta=betaR;
results(N+1).NL.betafull=betafull;
results(N+1).NL.XR=XR;

if scaling==1
    results(N+1).NL.scaled.coef(:,1)=betafull(:,end-1)/betafull(1,end-1); % Scale by first period impact on Fed Funds
    results(N+1).NL.scaled.coef(:,2)=betafull(:,end); % Scale by first period impact on Fed Funds
    results(N+1).NL.scaled.coef_diff=results(N+1).NL.scaled.coef(:,1)-results(N+1).NL.scaled.coef(:,2);

    results(N+1).NL.smoothed.coef(:,1)=weight*betafull(:,end-1)/betafull(1,end-1); % Scale by first period impact on Fed Funds
    results(N+1).NL.smoothed.coef(:,2)=weight*betafull(:,end); % Scale by first period impact on Fed Funds
    results(N+1).NL.smoothed.coef_diff=results(N+1).NL.smoothed.coef(:,1)-results(N+1).NL.smoothed.coef(:,2);
    
    results(N+1).NL.cum.coef(:,1)=cumweight*betafull(:,end-1)/betafull(1,end-1);
    results(N+1).NL.cum.coef(:,2)=cumweight*betafull(:,end);
    results(N+1).NL.cum.coef_diff=results(N+1).NL.cum.coef(:,1)-results(N+1).NL.cum.coef(:,2);       
end

% response variable regressions
for n=1:N; % loop over LHS variables
    
    results(N+1).NL.beta=[];
    res=[];
    X=[];
    beta=[];
    betafull=[];
    betafullNL=[];
    
    if Xlags>0
        X=LHS(startS-1:endS-1,n);
        if Xlags>1
            for l=2:Xlags
                X=[X LHS(startS-l:endS-l,n)];
            end
        end
    end
%    X2=X;
    if includetrend==0
        X=[X XR]; 
    else
        X=[trend X XR(:,2:end)];
    end
    results(n).NL.X=X;
    Xproj=inv(X'*X)*X';
    K=size(X,2);
    linshockindex=K;

    results(n).NL.beta=[];
    resNL=[];
    XNL=[];
    betaNL=[];
    betaNL=zeros(H+1,2);
    KNL=size(XNL,2);
    recshockindex=KNL;
    boomshockindex=KNL-(KNL-includetrend)/2;
    
    % loop over H, store coefficients and residuals
    for h=0:H
        betah=Xproj*LHS(startS+h:endS+h,n);
        betafull(h+1,:)=betah;
        beta=[beta; betah(end)];
        resh=LHS(startS+h:endS+h,n)-X*betah;
        res(:,h+1)=resh;

    end
    
    results(n).NL.beta=beta;
    results(n).NL.res=res;
    results(n).NL.betafull=betafull;
        
if scaling==1
    results(n).NL.scaled.coef(:,1)=results(n).NL.betafull(:,end-1)/results(N+1).NL.betafull(1,end-1); % Scale by first period impact on Fed Funds
    results(n).NL.scaled.coef(:,2)=results(n).NL.betafull(:,end); % Scale by first period impact on Fed Funds
    results(n).NL.scaled.coef_diff=results(n).NL.scaled.coef(:,1)-results(n).NL.scaled.coef(:,2);
    
    results(n).NL.smoothed.coef(:,1)=weight*results(n).NL.betafull(:,end-1)/results(N+1).NL.betafull(1,end-1); % Scale by first period impact on Fed Funds
    results(n).NL.smoothed.coef(:,2)=weight*results(n).NL.betafull(:,end); % Scale by first period impact on Fed Funds
    results(n).NL.smoothed.coef_diff=results(n).NL.smoothed.coef(:,1)-results(n).NL.smoothed.coef(:,2);

    results(n).NL.cum.coef(:,1)=cumweight*results(n).NL.betafull(:,end-1)/results(N+1).NL.betafull(1,end-1);
    results(n).NL.cum.coef(:,2)=cumweight*results(n).NL.betafull(:,end);
    results(n).NL.cum.coef_diff=results(n).NL.cum.coef(:,1)-results(n).NL.cum.coef(:,2);       
    
end

end
%display('Estimation complete')

end

