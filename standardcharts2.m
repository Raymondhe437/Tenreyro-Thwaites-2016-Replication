function [] = standardcharts(results,varindices,chartname,filename,rownames,colnames,H,CI,CIparam,scaling)
h=linspace(1,H+1,H+1);
zero=zeros(1,H+1);
tstathigh=CIparam*ones(1,H+1);
tstatlow=-CIparam*ones(1,H+1);

figure('Name',chartname,'NumberTitle','off','Units','normalized','Position',[0,0,.5,.5])
rows=length(varindices);
columns=5;

% First row
index=varindices(1);
subplot(rows,columns,1)
plot(1:H+1,results(index).lin.(scaling{1}).coef(1:21),'-b',1:H+1,results(index).NL.(scaling{1}).coef(1:21,1),'--g',1:H+1,results(index).NL.(scaling{1}).coef(1:21,2),'r-.',1:H+1,zeros(1,H+1),'-black','LineWidth', 1.2);
ylabel(colnames(1),'FontName','Arial Black')
set(gcf,'Color', [1 1 1])
xlim([0 20])

subplot(rows,columns,2)
patch([h fliplr(h)], [[results(index).lin.(scaling{1}).CIhigh(1:21)]',fliplr([results(index).lin.(scaling{1}).CIlow(1:21)]')],'c');
hold on
plot([h fliplr(h)], [[results(index).lin.(scaling{1}).coef(1:21)]' zero]);
xlim([0 20])
set(gca,'XTick', [0 10 20])
box on
hold off

subplot(rows,columns,3)
patch([h fliplr(h)], [[results(index).NL.(scaling{1}).CIhigh_pos(1:21)]',fliplr([results(index).NL.(scaling{1}).CIlow_pos(1:21)]')],'c');
hold on
plot([h fliplr(h)], [[results(index).NL.(scaling{1}).coef(1:21,1)]' zero]);
title(rownames(1),'FontName','Arial Black')
xlim([0 20])
set(gca,'XTick', [0 10 20])
box on
hold off

subplot(rows,columns,4)
patch([h fliplr(h)], [[results(index).NL.(scaling{1}).CIhigh_neg(1:21)]',fliplr([results(index).NL.(scaling{1}).CIlow_neg(1:21)]')],'c');
hold on
plot([h fliplr(h)], [[results(index).NL.(scaling{1}).coef(1:21,2)]' zero]);
set(gca,'YTick',[-0.01 -0.005 0 0.005 0.01])
set(gca,'XTick', [0 10 20])
xlim([0 20])
box on
hold off

subplot(rows,columns,5)
patch([h fliplr(h)], [tstathigh,fliplr(tstatlow)],'c');
hold on
plot(1:H+1,results(index).NL.(scaling{1}).tstats_diff(1:21),'black-',1:H+1,results(index).NL.(scaling{1}).betadifft(1:21),'red-');%,1:H,results(index).pos.tstatsequalDK(2:21),'green-');
set(gca,'XTick', [0 10 20])
xlim([0 20])
box on
hold off
ylabel('t statistic','FontName','Arial Black')

% Second row
index=varindices(2);
subplot(rows,columns,columns+1)
plot(1:H+1,results(index).lin.(scaling{2}).coef(1:21),'b-',1:H+1,results(index).NL.(scaling{2}).coef(1:21,1),'--g',1:H+1,results(index).NL.(scaling{2}).coef(1:21,2),'r-.',1:H+1,zeros(1,H+1),'black-','LineWidth', 1.2);
ylabel(colnames(2),'FontName','Arial Black');
xlim([0 20])
set(gca,'XTick', [0 10 20])
box on

subplot(rows,columns,columns+2)
patch([h fliplr(h)], [[results(index).lin.(scaling{2}).CIhigh(1:21)]',fliplr([results(index).lin.(scaling{2}).CIlow(1:21)]')],'cyan');
hold on
plot([h fliplr(h)], [[results(index).lin.(scaling{2}).coef(1:21)]' zero]);
set(gca,'YTick',[-0.03 0 0.02])
set(gca,'XTick', [0 10 20])
xlim([0 20])
box on
hold off

subplot(rows,columns,columns+3)
patch([h fliplr(h)], [[results(index).NL.(scaling{2}).CIhigh_pos(1:21)]',fliplr([results(index).NL.(scaling{2}).CIlow_pos(1:21)]')],'cyan');
hold on
plot([h fliplr(h)], [[results(index).NL.(scaling{2}).coef(1:21,1)]' zero]);
title(rownames(2),'FontName','Arial Black')
xlim([0 20])
set(gca,'XTick', [0 10 20])
box on
hold off

subplot(rows,columns,columns+4)
patch([h fliplr(h)], [[results(index).NL.(scaling{2}).CIhigh_neg(1:21)]',fliplr([results(index).NL.(scaling{2}).CIlow_neg(1:21)]')],'cyan');
hold on
plot([h fliplr(h)], [[results(index).NL.(scaling{2}).coef(1:21,2)]' zero]);
xlim([0 20])
set(gca,'XTick', [0 10 20])
box on
hold off

subplot(rows,columns,columns+5)
patch([h fliplr(h)], [tstathigh,fliplr(tstatlow)],'cyan');
hold on
plot(1:H+1,results(index).NL.(scaling{2}).tstats_diff(1:21),'black-',1:H+1,results(index).NL.(scaling{2}).betadifft(1:21),'red-');%,1:H,results(index).pos.tstatsequalcumDK(2:21),'green-')
xlim([0 20])
set(gca,'XTick', [0 10 20])
box on
hold off
ylabel('t statistic','FontName','Arial Black')

if rows==4
    % Third row
    index=varindices(3);
    subplot(rows,columns,2*columns+1)
    plot(1:H+1,results(index).lin.(scaling{3}).coef(1:21),'b-',1:H+1,results(index).NL.(scaling{3}).coef(1:21,1),'--g',1:H+1,results(index).NL.(scaling{3}).coef(1:21,2),'r-.',1:H+1,zeros(1,H+1),'black-','LineWidth', 1.2);
    ylabel(colnames(3),'FontName','Arial Black');
    xlim([0 20])
    set(gca,'XTick', [0 10 20])
    box on

    subplot(rows,columns,2*columns+2)
    patch([h fliplr(h)], [[results(index).lin.(scaling{3}).CIhigh(1:21)]',fliplr([results(index).lin.(scaling{3}).CIlow(1:21)]')],'cyan');
    hold on
    plot([h fliplr(h)], [[results(index).lin.(scaling{3}).coef(1:21)]' zero]);
    set(gca,'YTick',[-0.03 0 0.02])
    set(gca,'XTick', [0 10 20])
    xlim([0 20])
    box on
    hold off

    subplot(rows,columns,2*columns+3)
    patch([h fliplr(h)], [[results(index).NL.(scaling{3}).CIhigh_pos(1:21)]',fliplr([results(index).NL.(scaling{3}).CIlow_pos(1:21)]')],'cyan');
    hold on
    plot([h fliplr(h)], [[results(index).NL.(scaling{3}).coef(1:21,1)]' zero]);
    title(rownames(3),'FontName','Arial Black')
    xlim([0 20])
    set(gca,'XTick', [0 10 20])
    box on
    hold off

    subplot(rows,columns,2*columns+4)
    patch([h fliplr(h)], [[results(index).NL.(scaling{3}).CIhigh_neg(1:21)]',fliplr([results(index).NL.(scaling{3}).CIlow_neg(1:21)]')],'cyan');
    hold on
    plot([h fliplr(h)], [[results(index).NL.(scaling{3}).coef(1:21,2)]' zero]);
    xlim([0 20])
    set(gca,'XTick', [0 10 20])
    box on
    hold off

    subplot(rows,columns,2*columns+5)
    patch([h fliplr(h)], [tstathigh,fliplr(tstatlow)],'cyan');
    hold on
    plot(1:H+1,results(index).NL.(scaling{3}).tstats_diff(1:21),'black-',1:H+1,results(index).NL.(scaling{3}).betadifft(1:21),'red-');%,1:H,results(index).pos.tstatsequalcumDK(2:21),'green-')
    xlim([0 20])
    set(gca,'XTick', [0 10 20])
    box on
    hold off
    ylabel('t statistic','FontName','Arial Black')
    
end

% Bottom row
index=varindices(end);
subplot(rows,columns,(length(varindices)-1)*columns+1)
plot(1:H+1,results(index).lin.(scaling{end}).coef(1:21),'b-',1:H+1,results(index).NL.(scaling{end}).coef(1:21,1),'g--',1:H+1,results(index).NL.(scaling{end}).coef(1:21,2),'r-.',1:H,zeros(1,H),'black-','LineWidth', 1.2);
ylabel(colnames(end),'FontName','Arial Black');
xlabel('Three models','FontName','Arial Black')
xlim([0 20])
set(gca,'XTick', [0 10 20])
box on

subplot(rows,columns,(length(varindices)-1)*columns+2)
patch([h fliplr(h)], [[results(index).lin.(scaling{end}).CIhigh(1:21)]',fliplr([results(index).lin.(scaling{end}).CIlow(1:21)]')],'cyan');
xlabel('Linear model','FontName','Arial Black')
hold on
plot([h fliplr(h)], [[results(index).lin.(scaling{end}).coef(1:21)]' zero]);
xlim([0 20])
set(gca,'XTick', [0 10 20])
box on
hold off

subplot(rows,columns,(length(varindices)-1)*columns+3)
patch([h fliplr(h)], [[results(index).NL.(scaling{end}).CIhigh_pos(1:21)]',fliplr([results(index).NL.(scaling{end}).CIlow_pos(1:21)]')],'cyan');
xlabel('Tightening','FontName','Arial Black')
title(rownames(end),'FontName','Arial Black')
hold on
plot([h fliplr(h)], [[results(index).NL.(scaling{end}).coef(1:21,1)]' zero]);
xlim([0 20])
set(gca,'XTick', [0 10 20])
box on
hold off

subplot(rows,columns,(length(varindices)-1)*columns+4)
patch([h fliplr(h)], [[results(index).NL.(scaling{end}).CIhigh_neg(1:21)]',fliplr([results(index).NL.(scaling{end}).CIlow_neg(1:21)]')],'cyan');
xlabel('Loosening','FontName','Arial Black')
hold on
plot([h fliplr(h)], [[results(index).NL.(scaling{end}).coef(1:21,2)]' zero]);
xlim([0 20])
set(gca,'XTick', [0 10 20])
box on
hold off

subplot(rows,columns,(length(varindices)-1)*columns+5)
patch([h fliplr(h)], [tstathigh,fliplr(tstatlow)],'cyan');
xlabel({'Tightening =';'Loosening'},'FontName','Arial Black')
hold on
plot(1:H+1,results(index).NL.(scaling{end}).tstats_diff(1:21),'black-',1:H+1,results(index).NL.(scaling{end}).betadifft(1:21),'red-');%,1:H,results(index).pos.tstatsequalDK(2:21),'red-');
xlim([0 20])
set(gca,'XTick', [0 10 20])
box on
hold off
ylabel('t statistic','FontName','Arial Black')

saveas(gcf,filename,'pdf');

end

