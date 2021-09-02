
dat = readtable('ocount_1_type.txt');

x = table2array(dat(:,5));

y = table2array(dat(:,4));            

plot(x, y);                           

typestring = table2array(dat(1, 3));

th = title(typestring);

set(th,'Interpreter','none','FontSize',14,'FontWeight','bold');

