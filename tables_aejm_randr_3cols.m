
%Create the table for horizons 4,8,12 and 16 for GDP variable considering
%three different s.e.
table=zeros(8,3);
table(1,1)=results(1).NL.cum.coef(5,2);
table(1,2)=min(normcdf(results(1).NL.cum.tstats_neg(5)),1-normcdf(results(1).NL.cum.tstats_neg(5)));
table(1,3)=1-results(1).NL.cum.betadiffshare(5);
table(2,1)=results(1).NL.cum.coef(9,2);
table(2,2)=min(normcdf(results(1).NL.cum.tstats_neg(9)),1-normcdf(results(1).NL.cum.tstats_neg(9)));
table(2,3)=1-results(1).NL.cum.betadiffshare(9);
table(3,1)=results(1).NL.cum.coef(13,2);
table(3,2)=min(normcdf(results(1).NL.cum.tstats_neg(13)),1-normcdf(results(1).NL.cum.tstats_neg(13)));
table(3,3)=1-results(1).NL.cum.betadiffshare(13);
table(4,1)=results(1).NL.cum.coef(17,2);
table(4,2)=min(normcdf(results(1).NL.cum.tstats_neg(17)),1-normcdf(results(1).NL.cum.tstats_neg(17)));
table(4,3)=1-results(1).NL.cum.betadiffshare(17);

%Continue the table for horizons 4,8,12 and 16 for Inflation variable,
%considering three different s.e.
table(5,1)=results(2).NL.cum.coef(5,2);
table(5,2)=min(normcdf(results(2).NL.cum.tstats_neg(5)),1-normcdf(results(2).NL.cum.tstats_neg(5)));
table(5,3)=1-results(2).NL.cum.betadiffshare(5);
table(6,1)=results(2).NL.cum.coef(9,2);
table(6,2)=min(normcdf(results(2).NL.cum.tstats_neg(9)),1-normcdf(results(2).NL.cum.tstats_neg(9)));
table(6,3)=1-results(2).NL.cum.betadiffshare(9);
table(7,1)=results(2).NL.cum.coef(13,2);
table(7,2)=min(normcdf(results(2).NL.cum.tstats_neg(13)),1-normcdf(results(2).NL.cum.tstats_neg(13)));
table(7,3)=1-results(2).NL.cum.betadiffshare(13);
table(8,1)=results(2).NL.cum.coef(17,2);
table(8,2)=min(normcdf(results(2).NL.cum.tstats_neg(17)),1-normcdf(results(2).NL.cum.tstats_neg(17)));
table(8,3)=1-results(2).NL.cum.betadiffshare(17);
