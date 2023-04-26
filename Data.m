Rf = 0.01;
Cf=1e-4;
Lf=0.1;
Rc = 0.01;
Lc = 0.00125;

wbase = 2*pi*50;
Vdc = sqrt(2)*100;
Vnom = 400;
Pbase = 1000;
Vbase = Vnom*sqrt(2/3);
Ibase = 2*(Pbase/Vbase); %vdid +vqiq = 2P

%FileData = load('data_values.mat');
%csvwrite('data_optimal.csv', FileData.Data);
