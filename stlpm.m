function [results]=stlpm(R,E,M,LHS,trend,Xlags,Rlags,runstate,includetrend,useNLshocks,H,CI,startS,endS,scaling,weight,cumweight)
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

XR=[ones(length(XR),1) XR E(:,1)];   % construct common RHS projection matrix
if includetrend==1
    XR=[trend XR];
end
XRproj=inv(XR'*XR)*XR';
KR=size(XR,2);

if runstate==1;

    if useNLshocks==0
        XRNL=[repmat(M,[1 size(XR(:,1+includetrend:end),2)]).*XR(:,1+includetrend:end) (1-repmat(M,[1 size(XR(:,1+includetrend:end),2)])).*XR(:,1+includetrend:end)];
    else
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
betaR=[];
resR=[];
betafull=[];
betafullNL=[];
betaRNL=zeros(H+1,2);
resRNL=[];
for h=0:H

    betaRh=XRproj*R(startS+h:endS+h);
    results(N+1).lin.betafull=betaRh;
    betaR=[betaR; betaRh(end)];
    betafull(:,h+1)=betaRh;
    resRh=R(startS+h:endS+h)-XR*betaRh;
    resR(:,h+1)=resRh;
%    results(N+1).lin.variance(:,:,h+1)=var(resR(:,h+1))*inv(XR'*XR);
    
    if runstate==1
        betaRhNL=XRNLproj*R(startS+h:endS+h);
        betafullNL(:,h+1)=betaRhNL;
        betaRNL(h+1,1)=betaRhNL(1+(end-includetrend)/2);
        betaRNL(h+1,2)=betaRhNL(end);
        resRhNL=R(startS+h:endS+h)-XRNL*betaRhNL;
        resRNL(:,h+1)=resRhNL;
%        results(N+1).NL.variance(:,:,h+1)=var(resRNL(:,h+1))*inv(XRNL'*XRNL);
    end
    
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
for n=1:N; % loop over LHS variables
    
    results(N+1).lin.beta=[];
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
    results(n).lin.X=X;
    Xproj=inv(X'*X)*X';
    K=size(X,2);
    linshockindex=K;

    results(n).NL.beta=[];
    resNL=[];
    XNL=[];
    betaNL=[];
    betaNL=zeros(H+1,2);
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
    boomshockindex=KNL-(KNL-includetrend)/2;
    
    % loop over H, store coefficients and residuals
    for h=0:H
        betah=Xproj*LHS(startS+h:endS+h,n);
        betafull(:,h+1)=betah;
        beta=[beta; betah(end)];
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
        results(n).lin.scaled.coef=results(n).lin.beta/betaR(1); % Scale by first period impact on Fed Funds
        results(n).lin.smoothed.coef=weight*results(n).lin.beta/betaR(1); % Scale by first period impact on Fed Funds
        results(n).lin.cum.coef=cumweight*results(n).lin.beta/betaR(1);

        if runstate==1;
            results(n).NL.scaled.coef(:,1)=results(n).NL.betaNL(:,1)/betaRNL(1);
            results(n).NL.scaled.coef(:,2)=results(n).NL.betaNL(:,2)/betaRNL(2);
            results(n).NL.scaled.coef_diff=results(n).NL.scaled.coef(:,1)-results(n).NL.scaled.coef(:,2);

            results(n).NL.smoothed.coef(:,1)=weight*results(n).NL.scaled.coef(:,1);
            results(n).NL.smoothed.coef(:,2)=weight*results(n).NL.scaled.coef(:,2);
            results(n).NL.smoothed.coef_diff=results(n).NL.smoothed.coef(:,1)-results(n).NL.smoothed.coef(:,2);

            results(n).NL.cum.coef(:,1)=cumweight*results(n).NL.scaled.coef(:,1);
            results(n).NL.cum.coef(:,2)=cumweight*results(n).NL.scaled.coef(:,2);        
            results(n).NL.cum.coef_diff=results(n).NL.cum.coef(:,1)-results(n).NL.cum.coef(:,2);       
        end
    end
    
end
%display('Estimation complete')

end

