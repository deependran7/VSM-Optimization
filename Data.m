Rf = 0.01;
Cf=1e-4;
Lf=0.1;
Rc = 0.01;
Lc = 0.00125;

wbase = 2*pi*50;
Vbase =400;
Sbase = 1000;
Vbase = Vbase/sqrt(3);
Ibase = Sbase/Vbase;
Vbase = Vbase*sqrt(2);
Ibase = Ibase*sqrt(2);
Zbase = Vbase/Ibase;
Rbase = Vbase/Ibase;
Lbase = Zbase/(wbase);
RLoad = 0.8*Zbase ;           %115.47
XLoad = 0.6*Zbase ;           %153.96
LLoad = XLoad/wbase;
Vdc = sqrt(2)*100;

FileData = load('data_values.mat');
csvwrite('data_optimal.csv', FileData.Data);