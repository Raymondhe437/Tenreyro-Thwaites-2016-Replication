KX=size(results(1).lin.X,2);
KR=size(results(end).lin.XR,2);
K=KX+KR;
if runstate==1;
    KNL=size(results(1).NL.XNL,2);
    KRNL=size(results(end).NL.XRNL,2);
    KXRNL=KNL+KRNL;
end
recshockindex=KNL;
boomshockindex=KNL-(KNL-includetrend)/2;
CIparam=norminv(1-(1-CI)/2,0,1);

for n=1:size(results,2)-1; % loop over LHS variables

for h=0:H
    % heteroscedastic but uncorrelated case
    results(n).lin.variance(:,:,h+1)=var(results(n).lin.res(:,h+1))*inv(results(n).lin.X'*results(n).lin.X);
    if runstate==1;
        results(n).NL.variance(:,:,h+1)=var(results(n).NL.res(:,h+1))*inv(results(n).NL.XNL'*results(n).NL.XNL);
    end
end

    results(n).lin.variance_beta=squeeze(results(n).lin.variance(end,end,:));
    results(n).lin.varDK=NeweyWestDK([results(n).lin.res; repmat(results(end).lin.res(:,1),[1 H+1])],blkdiag(results(n).lin.X,results(end).lin.XR),H+1,0);
    results(n).lin.varDK_beta=results(n).lin.varDK([KX+K*[0:H] K*(H+1)],[KX+K*[0:H] K*(H+1)]);

    if runstate==1;
        KXRNL=KNL+KRNL;
        results(n).NL.variance_beta=squeeze(results(n).NL.variance([(KNL-(KNL-includetrend)/2) KNL],[(KNL-(KNL-includetrend)/2) KNL],:));
        results(n).NL.varDK=NeweyWestDK([results(n).NL.res; repmat(results(end).NL.res(:,1),[1 H+1])],blkdiag(results(n).NL.XNL,results(end).NL.XRNL),H+1,0);
        indices=[(boomshockindex+KXRNL*[0:H]) (recshockindex+KXRNL*[0:H]) (recshockindex+KXRNL*H+(KR-(KRNL-includetrend)/2)) (recshockindex+KXRNL*H+KRNL)];
        results(n).NL.varDK_beta=results(n).NL.varDK(indices,indices);
    end

end

for h=0:H
    results(end).lin.variance(:,:,h+1)=var(results(end).lin.res(:,h+1))*inv(results(end).lin.XR'*results(end).lin.XR);
end

results(end).lin.variance_beta=squeeze(results(end).lin.variance(end,end,:));
results(end).lin.varDK=NeweyWestDK(results(end).lin.res,results(end).lin.XR,H+1,0);
results(end).lin.varDK_beta=results(end).lin.varDK(KR*[1:H+1],KR*[1:H+1]);
if runstate==1;
    results(end).NL.variance=var(results(end).NL.res(:,h+1))*inv(results(end).NL.XRNL'*results(end).NL.XRNL);
    betaindices=[KRNL-(KRNL-includetrend)/2 KRNL];
    results(end).NL.variance_beta=squeeze(results(end).NL.variance(betaindices,betaindices,:));
    results(end).NL.varDK=NeweyWestDK([results(end).NL.res; repmat(results(end).NL.res(:,1),[1 H+1])],results(end).NL.XRNL,H+1,0);
    indicesR=[KRNL*[1:H+1]-(KRNL-includetrend)/2 KRNL*[1:H+1]];
    results(end).NL.varDK_beta=results(end).NL.varDK(indicesR,indicesR);
end


% delta method 
for n=1:size(results,2)-1; % loop over LHS variables

results(n).lin.scaled.delta=[eye(H+1)/results(end).lin.betafull(end,1) -results(n).lin.beta/(results(end).lin.betafull(end,1)^2)];
results(n).lin.scaled.var=results(n).lin.scaled.delta*results(n).lin.varDK_beta*results(n).lin.scaled.delta';
results(n).lin.scaled.SE=sqrt(diag(results(n).lin.scaled.var));
results(n).lin.scaled.tstats=results(n).lin.scaled.coef./results(n).lin.scaled.SE; % T stats defined as standard
results(n).lin.scaled.CIhigh =  results(n).lin.scaled.coef + CIparam*results(n).lin.scaled.SE; % 90% confidence intervals
results(n).lin.scaled.CIlow  =  results(n).lin.scaled.coef - CIparam*results(n).lin.scaled.SE; % 90% confidence intervals

results(n).lin.smoothed.delta=[weight/results(end).lin.betafull(end,1) -weight*results(n).lin.beta/(results(end).lin.betafull(end,1)^2)];
results(n).lin.smoothed.var=results(n).lin.smoothed.delta*results(n).lin.varDK_beta*results(n).lin.smoothed.delta';
results(n).lin.smoothed.SE=sqrt(diag(results(n).lin.smoothed.var));
results(n).lin.smoothed.tstats=results(n).lin.smoothed.coef./results(n).lin.smoothed.SE; % T stats defined as standard
results(n).lin.smoothed.CIhigh =  results(n).lin.smoothed.coef + CIparam*results(n).lin.smoothed.SE; % 90% confidence intervals
results(n).lin.smoothed.CIlow =  results(n).lin.smoothed.coef - CIparam*results(n).lin.smoothed.SE; % 90% confidence intervals

results(n).lin.cum.delta=[cumweight/results(end).lin.betafull(end,1) -cumweight*results(n).lin.beta/(results(end).lin.betafull(end,1)^2)];
results(n).lin.cum.var=results(n).lin.cum.delta*results(n).lin.varDK_beta*results(n).lin.cum.delta';
results(n).lin.cum.SE=sqrt(diag(results(n).lin.cum.var));
results(n).lin.cum.tstats=results(n).lin.cum.coef./results(n).lin.cum.SE; % T stats defined as standard
results(n).lin.cum.CIhigh =  results(n).lin.cum.coef + CIparam*results(n).lin.cum.SE; % 90% confidence intervals
results(n).lin.cum.CIlow  =  results(n).lin.cum.coef - CIparam*results(n).lin.cum.SE; % 90% confidence intervals

if runstate==1; 
    % scaled impulse responses
    results(n).NL.scaled.delta_boom=[eye(H+1)/results(end).NL.beta(1,1) zeros(H+1) -results(n).NL.betaNL(:,1)/(results(end).NL.beta(1,1)^2) zeros(H+1,1)];
    results(n).NL.scaled.delta_recession=[zeros(H+1) eye(H+1)/results(end).NL.beta(1,2) zeros(H+1,1) -results(n).NL.betaNL(:,2)/(results(end).NL.beta(1,2)^2)];  
    results(n).NL.scaled.delta_diff=[eye(H+1)/results(end).NL.beta(1,1) -eye(H+1)/results(end).NL.beta(1,2) -results(n).NL.betaNL(:,1)/(results(end).NL.beta(1,1)^2) results(n).NL.betaNL(:,2)/(results(end).NL.beta(1,2)^2)];

    results(n).NL.scaled.var_boom=results(n).NL.scaled.delta_boom*results(n).NL.varDK_beta*results(n).NL.scaled.delta_boom';
    results(n).NL.scaled.var_recession=results(n).NL.scaled.delta_recession*results(n).NL.varDK_beta*results(n).NL.scaled.delta_recession';
    results(n).NL.scaled.var_diff=results(n).NL.scaled.delta_diff*results(n).NL.varDK_beta*results(n).NL.scaled.delta_diff';

    results(n).NL.scaled.SE_boom=sqrt(diag(results(n).NL.scaled.var_boom));
    results(n).NL.scaled.SE_recession=sqrt(diag(results(n).NL.scaled.var_recession));
    results(n).NL.scaled.SE_diff=sqrt(diag(results(n).NL.scaled.var_diff));

    results(n).NL.scaled.tstats_boom=results(n).NL.scaled.coef(:,1)./results(n).NL.scaled.SE_boom; % T stats defined as standard
    results(n).NL.scaled.tstats_recession=results(n).NL.scaled.coef(:,2)./results(n).NL.scaled.SE_recession; % T stats defined as standard
    results(n).NL.scaled.tstats_diff=results(n).NL.scaled.coef_diff./sqrt(diag(results(n).NL.scaled.var_diff));

    results(n).NL.scaled.CIhigh_boom =  results(n).NL.scaled.coef(:,1) + CIparam*results(n).NL.scaled.SE_boom; % 90% confidence intervals
    results(n).NL.scaled.CIhigh_recession =  results(n).NL.scaled.coef(:,2) + CIparam*results(n).NL.scaled.SE_recession; % 90% confidence intervals
    results(n).NL.scaled.CIlow_boom =  results(n).NL.scaled.coef(:,1) - CIparam*results(n).NL.scaled.SE_boom; % 90% confidence intervals
    results(n).NL.scaled.CIlow_recession =  results(n).NL.scaled.coef(:,2) - CIparam*results(n).NL.scaled.SE_recession; % 90% confidence intervals

    % smoothed impulse responses
    results(n).NL.smoothed.delta_boom=[weight/results(end).NL.beta(1,1) zeros(H+1) -weight*results(n).NL.betaNL(:,1)/(results(end).NL.beta(1,1)^2) zeros(H+1,1)];
    results(n).NL.smoothed.delta_recession=[zeros(H+1) weight/results(end).NL.beta(1,2) zeros(H+1,1) -weight*results(n).NL.betaNL(:,2)/(results(end).NL.beta(1,2)^2)];  
    results(n).NL.smoothed.delta_diff=[weight/results(end).NL.beta(1,1) -weight/results(end).NL.beta(1,2) -weight*results(n).NL.betaNL(:,1)/(results(end).NL.beta(1,1)^2) weight*results(n).NL.betaNL(:,2)/(results(end).NL.beta(1,2)^2)];

    results(n).NL.smoothed.var_boom=results(n).NL.smoothed.delta_boom*results(n).NL.varDK_beta*results(n).NL.smoothed.delta_boom';
    results(n).NL.smoothed.var_recession=results(n).NL.smoothed.delta_recession*results(n).NL.varDK_beta*results(n).NL.smoothed.delta_recession';
    results(n).NL.smoothed.var_diff=results(n).NL.smoothed.delta_diff*results(n).NL.varDK_beta*results(n).NL.smoothed.delta_diff';

    results(n).NL.smoothed.SE_boom=sqrt(diag(results(n).NL.smoothed.var_boom));
    results(n).NL.smoothed.SE_recession=sqrt(diag(results(n).NL.smoothed.var_recession));
    results(n).NL.smoothed.SE_diff=sqrt(diag(results(n).NL.smoothed.var_diff));

    results(n).NL.smoothed.tstats_boom=results(n).NL.smoothed.coef(:,1)./results(n).NL.smoothed.SE_boom; % T stats defined as standard
    results(n).NL.smoothed.tstats_recession=results(n).NL.smoothed.coef(:,2)./results(n).NL.smoothed.SE_recession; % T stats defined as standard
    results(n).NL.smoothed.tstats_diff=results(n).NL.smoothed.coef_diff./sqrt(diag(results(n).NL.smoothed.var_diff));

    results(n).NL.smoothed.CIhigh_boom =  results(n).NL.smoothed.coef(:,1) + CIparam*results(n).NL.smoothed.SE_boom; % 90% confidence intervals
    results(n).NL.smoothed.CIhigh_recession =  results(n).NL.smoothed.coef(:,2) + CIparam*results(n).NL.smoothed.SE_recession; % 90% confidence intervals
    results(n).NL.smoothed.CIlow_boom =  results(n).NL.smoothed.coef(:,1) - CIparam*results(n).NL.smoothed.SE_boom; % 90% confidence intervals
    results(n).NL.smoothed.CIlow_recession =  results(n).NL.smoothed.coef(:,2) - CIparam*results(n).NL.smoothed.SE_recession; % 90% confidence intervals

    % cumulated impulse responses
    results(n).NL.cum.delta_boom=[cumweight/results(end).NL.beta(1,1) zeros(H+1) -cumweight*results(n).NL.betaNL(:,1)/(results(end).NL.beta(1,1)^2) zeros(H+1,1)];
    results(n).NL.cum.delta_recession=[zeros(H+1) cumweight/results(end).NL.beta(1,2) zeros(H+1,1) -cumweight*results(n).NL.betaNL(:,2)/(results(end).NL.beta(1,2)^2)];  
    results(n).NL.cum.delta_diff=[cumweight/results(end).NL.beta(1,1) -cumweight/results(end).NL.beta(1,2) -cumweight*results(n).NL.betaNL(:,1)/(results(end).NL.beta(1,1)^2) cumweight*results(n).NL.betaNL(:,2)/(results(end).NL.beta(1,2)^2)];

    results(n).NL.cum.var_boom=results(n).NL.cum.delta_boom*results(n).NL.varDK_beta*results(n).NL.cum.delta_boom';
    results(n).NL.cum.var_recession=results(n).NL.cum.delta_recession*results(n).NL.varDK_beta*results(n).NL.cum.delta_recession';
    results(n).NL.cum.var_diff=results(n).NL.cum.delta_diff*results(n).NL.varDK_beta*results(n).NL.cum.delta_diff';

    results(n).NL.cum.SE_boom=sqrt(diag(results(n).NL.cum.var_boom));
    results(n).NL.cum.SE_recession=sqrt(diag(results(n).NL.cum.var_recession));
    results(n).NL.cum.SE_diff=sqrt(diag(results(n).NL.cum.var_diff));

    results(n).NL.cum.tstats_boom=results(n).NL.cum.coef(:,1)./results(n).NL.cum.SE_boom; % T stats defined as standard
    results(n).NL.cum.tstats_recession=results(n).NL.cum.coef(:,2)./results(n).NL.cum.SE_recession; % T stats defined as standard
    results(n).NL.cum.tstats_diff=results(n).NL.cum.coef_diff./sqrt(diag(results(n).NL.cum.var_diff));

    results(n).NL.cum.CIhigh_boom =  results(n).NL.cum.coef(:,1) + CIparam*results(n).NL.cum.SE_boom; % 90% confidence intervals
    results(n).NL.cum.CIhigh_recession =  results(n).NL.cum.coef(:,2) + CIparam*results(n).NL.cum.SE_recession; % 90% confidence intervals
    results(n).NL.cum.CIlow_boom =  results(n).NL.cum.coef(:,1) - CIparam*results(n).NL.cum.SE_boom; % 90% confidence intervals
    results(n).NL.cum.CIlow_recession =  results(n).NL.cum.coef(:,2) - CIparam*results(n).NL.cum.SE_recession; % 90% confidence intervals
end

end

% policy variable
results(end).lin.scaled.delta=eye(H+1)/results(end).lin.betafull(end,1);
results(end).lin.scaled.delta(:,1)=results(end).lin.scaled.delta(:,1)-results(end).lin.betafull(end,:)'/(results(end).lin.betafull(end,1)^2);
results(end).lin.scaled.var=results(end).lin.scaled.delta*results(end).lin.varDK_beta*results(end).lin.scaled.delta';
results(end).lin.scaled.SE=sqrt(diag(results(end).lin.scaled.var));
results(end).lin.scaled.tstats=results(end).lin.scaled.coef./results(end).lin.scaled.SE; % T stats defined as standard
results(end).lin.scaled.CIhigh =  results(end).lin.scaled.coef + CIparam*results(end).lin.scaled.SE; % 90% confidence intervals
results(end).lin.scaled.CIlow  =  results(end).lin.scaled.coef - CIparam*results(end).lin.scaled.SE; % 90% confidence intervals

results(end).lin.smoothed.delta=weight/results(end).lin.betafull(end,1);
results(end).lin.smoothed.delta(:,1)=results(end).lin.smoothed.delta(:,1)-weight*results(end).lin.betafull(end,:)'/(results(end).lin.betafull(end,1)^2);
results(end).lin.smoothed.var=results(end).lin.smoothed.delta*results(end).lin.varDK_beta*results(end).lin.smoothed.delta';
results(end).lin.smoothed.SE=sqrt(diag(results(end).lin.smoothed.var));
results(end).lin.smoothed.tstats=results(end).lin.smoothed.coef./results(end).lin.smoothed.SE; % T stats defined as standard
results(end).lin.smoothed.CIhigh =  results(end).lin.smoothed.coef + CIparam*results(end).lin.smoothed.SE; % 90% confidence intervals
results(end).lin.smoothed.CIlow =  results(end).lin.smoothed.coef - CIparam*results(end).lin.smoothed.SE; % 90% confidence intervals

results(end).lin.cum.delta=cumweight/results(end).lin.betafull(end,1);
results(end).lin.cum.delta(:,1)=results(end).lin.cum.delta(:,1)-cumweight*results(end).lin.betafull(end,:)'/(results(end).lin.betafull(end,1)^2);
results(end).lin.cum.var=results(end).lin.cum.delta*results(end).lin.varDK_beta*results(end).lin.cum.delta';
results(end).lin.cum.SE=sqrt(diag(results(end).lin.cum.var));
results(end).lin.cum.tstats=results(end).lin.cum.coef./results(end).lin.cum.SE; % T stats defined as standard
results(end).lin.cum.CIhigh =  results(end).lin.cum.coef + CIparam*results(end).lin.cum.SE; % 90% confidence intervals
results(end).lin.cum.CIlow  =  results(end).lin.cum.coef - CIparam*results(end).lin.cum.SE; % 90% confidence intervals

if runstate==1; 
    % scaled impulse responses
    results(end).NL.scaled.delta_boom=[eye(H+1)/results(end).NL.beta(1,1) zeros(H+1)];
    results(end).NL.scaled.delta_boom(:,1)=results(end).NL.scaled.delta_boom(:,1)-results(end).NL.beta(:,1)/(results(end).NL.beta(1,1)^2);
    results(end).NL.scaled.delta_recession=[zeros(H+1) eye(H+1)/results(end).NL.beta(1,2)];
    results(end).NL.scaled.delta_recession(:,H+2)=results(end).NL.scaled.delta_recession(:,H+2)-results(end).NL.beta(:,2)/(results(end).NL.beta(1,2)^2);
    results(end).NL.scaled.delta_diff=results(end).NL.scaled.delta_boom-results(end).NL.scaled.delta_recession;

    results(end).NL.scaled.var_boom=results(end).NL.scaled.delta_boom*results(end).NL.varDK_beta*results(end).NL.scaled.delta_boom';
    results(end).NL.scaled.var_recession=results(end).NL.scaled.delta_recession*results(end).NL.varDK_beta*results(end).NL.scaled.delta_recession';
    results(end).NL.scaled.var_diff=results(end).NL.scaled.delta_diff*results(end).NL.varDK_beta*results(end).NL.scaled.delta_diff';

    results(end).NL.scaled.SE_boom=sqrt(diag(results(end).NL.scaled.var_boom));
    results(end).NL.scaled.SE_recession=sqrt(diag(results(end).NL.scaled.var_recession));
    results(end).NL.scaled.SE_diff=sqrt(diag(results(end).NL.scaled.var_diff));

    results(end).NL.scaled.tstats_boom=results(end).NL.scaled.coef(:,1)./results(end).NL.scaled.SE_boom; % T stats defined as standard
    results(end).NL.scaled.tstats_recession=results(end).NL.scaled.coef(:,2)./results(end).NL.scaled.SE_recession; % T stats defined as standard
    results(end).NL.scaled.tstats_diff=results(end).NL.scaled.coef_diff./sqrt(diag(results(end).NL.scaled.var_diff));

    results(end).NL.scaled.CIhigh_boom =  results(end).NL.scaled.coef(:,1) + CIparam*results(end).NL.scaled.SE_boom; % 90% confidence intervals
    results(end).NL.scaled.CIhigh_recession =  results(end).NL.scaled.coef(:,2) + CIparam*results(end).NL.scaled.SE_recession; % 90% confidence intervals
    results(end).NL.scaled.CIlow_boom =  results(end).NL.scaled.coef(:,1) - CIparam*results(end).NL.scaled.SE_boom; % 90% confidence intervals
    results(end).NL.scaled.CIlow_recession =  results(end).NL.scaled.coef(:,2) - CIparam*results(end).NL.scaled.SE_recession; % 90% confidence intervals

    % smoothed impulse responses
    results(end).NL.smoothed.delta_boom=[weight/results(end).NL.beta(1,1) zeros(H+1)];
    results(end).NL.smoothed.delta_boom(:,1)=results(end).NL.smoothed.delta_boom(:,1)-weight*results(end).NL.beta(:,1)/(results(end).NL.beta(1,1)^2);
    results(end).NL.smoothed.delta_recession=[zeros(H+1) weight/results(end).NL.beta(1,2)];
    results(end).NL.smoothed.delta_recession(:,H+2)=results(end).NL.smoothed.delta_recession(:,H+2)-weight*results(end).NL.beta(:,2)/(results(end).NL.beta(1,2)^2);
    results(end).NL.smoothed.delta_diff=results(end).NL.smoothed.delta_boom-results(end).NL.smoothed.delta_recession;

    results(end).NL.smoothed.var_boom=results(end).NL.smoothed.delta_boom*results(end).NL.varDK_beta*results(end).NL.smoothed.delta_boom';
    results(end).NL.smoothed.var_recession=results(end).NL.smoothed.delta_recession*results(end).NL.varDK_beta*results(end).NL.smoothed.delta_recession';
    results(end).NL.smoothed.var_diff=results(end).NL.smoothed.delta_diff*results(end).NL.varDK_beta*results(end).NL.smoothed.delta_diff';

    results(end).NL.smoothed.SE_boom=sqrt(diag(results(end).NL.smoothed.var_boom));
    results(end).NL.smoothed.SE_recession=sqrt(diag(results(end).NL.smoothed.var_recession));
    results(end).NL.smoothed.SE_diff=sqrt(diag(results(end).NL.smoothed.var_diff));

    results(end).NL.smoothed.tstats_boom=results(end).NL.smoothed.coef(:,1)./results(end).NL.smoothed.SE_boom; % T stats defined as standard
    results(end).NL.smoothed.tstats_recession=results(end).NL.smoothed.coef(:,2)./results(end).NL.smoothed.SE_recession; % T stats defined as standard
    results(end).NL.smoothed.tstats_diff=results(end).NL.smoothed.coef_diff./sqrt(diag(results(end).NL.smoothed.var_diff));

    results(end).NL.smoothed.CIhigh_boom =  results(end).NL.smoothed.coef(:,1) + CIparam*results(end).NL.smoothed.SE_boom; % 90% confidence intervals
    results(end).NL.smoothed.CIhigh_recession =  results(end).NL.smoothed.coef(:,2) + CIparam*results(end).NL.smoothed.SE_recession; % 90% confidence intervals
    results(end).NL.smoothed.CIlow_boom =  results(end).NL.smoothed.coef(:,1) - CIparam*results(end).NL.smoothed.SE_boom; % 90% confidence intervals
    results(end).NL.smoothed.CIlow_recession =  results(end).NL.smoothed.coef(:,2) - CIparam*results(end).NL.smoothed.SE_recession; % 90% confidence intervals

    % cumulated impulse responses
    results(end).NL.cum.delta_boom=[cumweight/results(end).NL.beta(1,1) zeros(H+1)];
    results(end).NL.cum.delta_boom(:,1)=results(end).NL.cum.delta_boom(:,1)-cumweight*results(end).NL.beta(:,1)/(results(end).NL.beta(1,1)^2);
    results(end).NL.cum.delta_recession=[zeros(H+1) cumweight/results(end).NL.beta(1,2)];
    results(end).NL.cum.delta_recession(:,H+2)=results(end).NL.cum.delta_recession(:,H+2)-cumweight*results(end).NL.beta(:,2)/(results(end).NL.beta(1,2)^2);
    results(end).NL.cum.delta_diff=results(end).NL.cum.delta_boom-results(end).NL.cum.delta_recession;

    results(end).NL.cum.var_boom=results(end).NL.cum.delta_boom*results(end).NL.varDK_beta*results(end).NL.cum.delta_boom';
    results(end).NL.cum.var_recession=results(end).NL.cum.delta_recession*results(end).NL.varDK_beta*results(end).NL.cum.delta_recession';
    results(end).NL.cum.var_diff=results(end).NL.cum.delta_diff*results(end).NL.varDK_beta*results(end).NL.cum.delta_diff';

    results(end).NL.cum.SE_boom=sqrt(diag(results(end).NL.cum.var_boom));
    results(end).NL.cum.SE_recession=sqrt(diag(results(end).NL.cum.var_recession));
    results(end).NL.cum.SE_diff=sqrt(diag(results(end).NL.cum.var_diff));

    results(end).NL.cum.tstats_boom=results(end).NL.cum.coef(:,1)./results(end).NL.cum.SE_boom; % T stats defined as standard
    results(end).NL.cum.tstats_recession=results(end).NL.cum.coef(:,2)./results(end).NL.cum.SE_recession; % T stats defined as standard
    results(end).NL.cum.tstats_diff=results(end).NL.cum.coef_diff./sqrt(diag(results(end).NL.cum.var_diff));

    results(end).NL.cum.CIhigh_boom =  results(end).NL.cum.coef(:,1) + CIparam*results(end).NL.cum.SE_boom; % 90% confidence intervals
    results(end).NL.cum.CIhigh_recession =  results(end).NL.cum.coef(:,2) + CIparam*results(end).NL.cum.SE_recession; % 90% confidence intervals
    results(end).NL.cum.CIlow_boom =  results(end).NL.cum.coef(:,1) - CIparam*results(end).NL.cum.SE_boom; % 90% confidence intervals
    results(end).NL.cum.CIlow_recession =  results(end).NL.cum.coef(:,2) - CIparam*results(end).NL.cum.SE_recession; % 90% confidence intervals
end