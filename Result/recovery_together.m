
fh = figure();
strtitle = ['The recovery rates of propane in 1 year', char(13,10)'];
h = title(strtitle);
set(h, 'fontsize', 12, 'FontWeight', 'bold', 'HorizontalAlignment', 'center');
h = xlabel('Time(s)');
set(h, 'fontsize', 12, 'FontWeight', 'bold', 'HorizontalAlignment', 'center');
h = ylabel('Recovery Rate(%)');
set(h, 'fontsize', 12, 'FontWeight', 'bold', 'VerticalAlignment', 'middle');
hold on;

nt = 1*365*24;
timeEnd = 1*365*24*3600;
ts = (0:nt)*timeEnd/nt;

mfh = zeros(nt+1, 1);
temp = load('case1/soln_1PhFlw_moleHistory.txt');
k = 1;
for i = 2 : nt+1
    mfh(i) = temp(k);
    k = k + 1;
end    

plot(ts, mfh*100, 'r', 'linewidth', 2);
hold on;

temp = load('case2/soln_1PhFlw_moleHistory.txt');
k = 1;
for i = 2 : nt+1
    mfh(i) = temp(k);
    k = k + 1;
end    

plot(ts, mfh*100, 'b', 'linewidth', 2);
hold on;

temp = load('case3/soln_1PhFlw_moleHistory.txt');
k = 1;
for i = 2 : nt+1
    mfh(i) = temp(k);
    k = k + 1;
end    

plot(ts, mfh*100, 'k', 'linewidth', 2);
hold on;

temp = load('case4/soln_1PhFlw_moleHistory.txt');
k = 1;
for i = 2 : nt+1
    mfh(i) = temp(k);
    k = k + 1;
end    

plot(ts, mfh*100, 'y', 'linewidth', 2);
hold on;

temp = load('case5/soln_1PhFlw_moleHistory.txt');
k = 1;
for i = 2 : nt+1
    mfh(i) = temp(k);
    k = k + 1;
end    

plot(ts, mfh*100, 'g', 'linewidth', 2);
hold on;

temp = load('case6/soln_1PhFlw_moleHistory.txt');
k = 1;
for i = 2 : nt+1
    mfh(i) = temp(k);
    k = k + 1;
end    

plot(ts, mfh*100, 'c', 'linewidth', 2);
hold on;

legend('1','2','3','4','5','6','location','southeast');
hold off;

saveas(fh,'recovery_toghther');

