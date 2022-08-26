function [] = standardcharts3(results,varindices,chartname,filename,rownames,colnames,H,CI,CIparam,scaling)
h=linspace(1,H+1,H+1);
zero=zeros(1,H+1);
tstathigh=CIparam*ones(1,H+1);
tstatlow=-CIparam*ones(1,H+1);

figure('Name',chartname,'NumberTitle','off','Units','normalized','Position',[0,0,.5,.5])
rows=length(varindices);
columns=2;

% First row
index=varindices(1);
subplot(rows,columns,1)
patch([h fliplr(h)], [[results(index).NL.(scaling{1}).CIhigh_neg(1:21)]',fliplr([results(index).NL.(scaling{1}).CIlow_neg(1:21)]')],'c');
hold on
plot([h fliplr(h)], [[results(index).NL.(scaling{1}).coef(1:21,2)]' zero]);
set(gca,'YTick',[-0.01 -0.005 0 0.005 0.01])
set(gca,'XTick', [0 10 20])
xlim([0 20])
box on
hold off
ylabel(colnames(1),'FontName','Arial Black')
set(gcf,'Color', [1 1 1])
xlim([0 20])

subplot(rows,columns,2)
patch([h fliplr(h)], [tstathigh,fliplr(tstatlow)],'c');
hold on
plot(1:H+1,results(index).NL.(scaling{1}).tstats_neg(1:21),'black-',1:H+1,results(index).NL.(scaling{1}).betadifft(1:21),'red-');%,1:H,results(index).pos.tstatsequalDK(2:21),'green-');
set(gca,'XTick', [0 10 20])
xlim([0 20])
box on
hold off
ylabel('t statistic','FontName','Arial Black')

% Second row
index=varindices(2);
subplot(rows,columns,columns+1)
patch([h fliplr(h)], [[results(index).NL.(scaling{2}).CIhigh_neg(1:21)]',fliplr([results(index).NL.(scaling{2}).CIlow_neg(1:21)]')],'cyan');
hold on
plot([h fliplr(h)], [[results(index).NL.(scaling{2}).coef(1:21,2)]' zero]);
xlim([0 20])
set(gca,'XTick', [0 10 20])
box on
hold off
ylabel(colnames(2),'FontName','Arial Black');
xlim([0 20])
set(gca,'XTick', [0 10 20])
box on

subplot(rows,columns,columns+2)
patch([h fliplr(h)], [tstathigh,fliplr(tstatlow)],'cyan');
hold on
plot(1:H+1,results(index).NL.(scaling{2}).tstats_neg(1:21),'black-',1:H+1,results(index).NL.(scaling{2}).betadifft(1:21),'red-');%,1:H,results(index).pos.tstatsequalcumDK(2:21),'green-')
xlim([0 20])
set(gca,'XTick', [0 10 20])
box on
hold off
ylabel('t statistic','FontName','Arial Black')

% Bottom row
index=varindices(end);
subplot(rows,columns,(length(varindices)-1)*columns+1)
ylabel(colnames(end),'FontName','Arial Black');
xlabel('Cubed policy rate','FontName','Arial Black')
patch([h fliplr(h)], [[results(index).NL.(scaling{end}).CIhigh_neg(1:21)]',fliplr([results(index).NL.(scaling{end}).CIlow_neg(1:21)]')],'cyan');
hold on
plot([h fliplr(h)], [[results(index).NL.(scaling{end}).coef(1:21,2)]' zero]);
xlim([0 20])
set(gca,'XTick', [0 10 20])
box on
hold off

subplot(rows,columns,(length(varindices)-1)*columns+2)
patch([h fliplr(h)], [tstathigh,fliplr(tstatlow)],'cyan');
xlabel({'Cubed policy rate = 0'},'FontName','Arial Black')
hold on
plot(1:H+1,results(index).NL.(scaling{end}).tstats_neg(1:21),'black-',1:H+1,results(index).NL.(scaling{end}).betadifft(1:21),'red-');%,1:H,results(index).pos.tstatsequalDK(2:21),'red-');
xlim([0 20])
set(gca,'XTick', [0 10 20])
box on
hold off
ylabel('t statistic','FontName','Arial Black')

saveas(gcf,filename,'pdf');

end

