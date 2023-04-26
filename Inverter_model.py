import numpy as np
pi = np.pi
Vnom = 400
Pbase = 1000
Vbase = Vnom*np.sqrt(2/3)
Ibase = 2*(Pbase/Vbase) #vdid +vqiq = 2P
Zbase = Vbase/Ibase



#Inverter Variables(in pu)
w = 2*pi*50
Lf = 0.1
Rf = 0.01
Cf= 1*10**(-4)
Rc = 0.01
Lc = 0.00125
Vdc = np.sqrt(2)*100
#Load Variables

PLoad = 800;
QLoad = 600
Inom = Pbase/(np.sqrt(3)*Vnom)
RLoad = PLoad/(3*Inom**2)
XLoad = QLoad/(3*Inom**2)
T_inv = 0.0001
K_inv = 0.5


#Machine Variables(in pu)
xd =1.56
xdd=0.296
xq =2
xqd=0.4
Td0=3.7;
Tq0=0.6;
J=40;
D=0.2;
Rsg=0.00;
Ka=24;
Ta=0.05;
Tsv=0.05;
Tch=0.5;
Rd=0.005;
Pref=1.0;
Vref=1.0;

def Inverter(t,x):
    Ifd,Ifq,Vcd,Vcq,Iod,Ioq,Vid,Viq,Eqd,Edd,delta,d_w,Efd,Psv,Pm = x
    Vod = Iod*RLoad - Ioq*XLoad
    Voq = Iod*XLoad + Ioq*RLoad
    
    md = (Edd - Ioq*xqd/Ibase - Iod*Rsg/Ibase)
    mq = (Eqd - Ioq*Rsg/Ibase + Iod*xdd/Ibase)

    Vid_dot = (1/T_inv)*(-Vid + md*Vdc*K_inv)
    Viq_dot = (1/T_inv)*(-Viq + mq*Vdc*K_inv)
    Ifd_dot = (-Rf/Lf)*Ifd + w*Ifq + (Vid - Vcd)/Lf
    Ifq_dot = -w*Ifd + (-Rf/Lf)*Ifq + (Viq - Vcq)/Lf
    Vcd_dot = w*Vcq + (Ifd - Iod)/Cf
    Vcq_dot = -w*Vcd + (Ifq - Ioq)/Cf
    Iod_dot = (-(Rc)/(Lc))*Iod + w*Ioq + (Vcd -Vod)/Lc
    Ioq_dot = -w*Iod + (-(Rc)/(Lc))*Ioq + (Vcq -Voq)/Lc

    ##Machine Equations
    Eqd_dot = (-Eqd - (xd-xdd)*Iod/Ibase + Efd)/Td0 
    Edd_dot = (-Edd + (xq-xqd)*Ioq/Ibase)/Tq0 
    del_dot =  d_w
    d_w_dot = (Pm - (Eqd*Ioq/Ibase + Edd*Iod/Ibase - (xdd-xqd)*(Iod)*(Ioq)/(Ibase*Ibase)) + D*(d_w))/J 
    Efd_dot = (-Efd +  Ka*(Vref - np.sqrt((Vod/Vbase)**2 + (Voq/Vbase)**2)))/Ta
    Psv_dot = (-Psv + Pref - (d_w)/Rd)/Tsv 
    Pm_dot =  (-Pm + Psv)/Tch 
    
    return np.array([Ifd_dot,Ifq_dot,Vcd_dot,Vcq_dot,Iod_dot,Ioq_dot,Vid_dot,Viq_dot,Eqd_dot,Edd_dot,del_dot,d_w_dot,Efd_dot,Psv_dot,Pm_dot])
