import numpy as np
pi = np.pi
Vbase = 400/np.sqrt(3)
Sbase = 1000
Ibase = Sbase/Vbase
Vbase = np.sqrt(2)*Vbase
Ibase = np.sqrt(2)*Ibase
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
RLoad = 0.8*Zbase            #115.47
XLoad = 0.6*Zbase            #153.96
LLoad = XLoad/(w)

T_inv = 0.0001
K_inv = 1.2


Pref=1.0;
Vref=1.4;
X_optimal = [1.47527927e+00, 4.03544283e-01, 6.35738918e-02, 2.94491065e+00,
 5.82967049e-03, 1.73255221e-02, 5.71779349e+00, 9.79004084e-01,
 9.69122020e-01, 1.70774306e+00, 1.03728095e-04, 1.43135922e-04,
 1.22065566e-04, 5.79434404e-03, 141.42]

Td0,Tq0,xd,xq,xdd,xqd,J,D,Rsg,Ka,Ta,Tsv,Tch,Rd,Vdc = X_optimal


def Inverter_optimal(t,x):
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
    Iod_dot = (-(Rc+RLoad)/(Lc+LLoad))*Iod + w*Ioq + (Vcd -Vod)/Lc
    Ioq_dot = -w*Iod + (-(Rc+RLoad)/(Lc+LLoad))*Ioq + (Vcq -Voq)/Lc

    ##Machine Equations
    Eqd_dot = (-Eqd - (xd-xdd)*Iod/Ibase + Efd)/Td0 
    Edd_dot = (-Edd + (xq-xqd)*Ioq/Ibase)/Tq0 
    del_dot =  d_w
    d_w_dot = (Pm - (Eqd*Ioq/Ibase + Edd*Iod/Ibase - (xdd-xqd)*(Iod)*(Ioq)/(Ibase*Ibase)) + D*(d_w))/J 
    Efd_dot = (-Efd +  Ka*(Vref - np.sqrt((Vod/Vbase)**2 + (Voq/Vbase)**2)))/Ta
    Psv_dot = (-Psv + Pref - (d_w)/Rd)/Tsv 
    Pm_dot =  (-Pm + Psv)/Tch 
    
    return np.array([Ifd_dot,Ifq_dot,Vcd_dot,Vcq_dot,Iod_dot,Ioq_dot,Vid_dot,Viq_dot,Eqd_dot,Edd_dot,del_dot,d_w_dot,Efd_dot,Psv_dot,Pm_dot])
