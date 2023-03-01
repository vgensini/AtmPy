
__all__ = ['stp','scp','ship']

#sig tornado parameter
def stp(cape,s06,lcl,cin,srh):
    term1 = cape/1500.
    term2 = (2000.-lcl)/1000.
    term3 = srh/150.
    term4 = s06/20.
    term5 = (200+cin)/150.
    term2[lcl<1000.] = 1.0
    term2[lcl>2000.] = 0.0
    term5[cin>-50.] = 1.0
    term5[cin<-200.] = 0.0
    term4[s06>30.] = 1.5
    term4[s06<12.5] = 0.0
    stp = term1*term2*term3*term4*term5
    stp[stp<0.] = 0.0
    return stp

#supercell composite parameter
def scp(cape,srh,bwd,cin):
    term1 = cape/1000.
    term2 = srh/100.
    term3 = bwd/20.
    term3[bwd<10.] = 0.
    term3[bwd>20.] = 1.
    term4 = -40/cin
    term4[cin>=-40.] = 1.
    scp = term1*term2*term3*term4
    scp[scp<0.] = 0.0
    return scp

#sig hail parameter
def ship(mucape,mixr2m,lr57,t500,bs06,frzlvl):
    mumixr = mixr2m.copy()
    t500[t500>-5.] = -5.5
    bs06[bs06<7.] = 7.
    bs06[bs06>27.] = 27.
    mumixr[mumixr<11.] = 11.
    mumixr[mumixr>13.6] = 13.6
    ship = (mucape*mumixr*lr57*(-1.*t500)*bs06)/42000000.
    ship[mucape<1300.] = ship[mucape<1300.] * (mucape[mucape<1300.]/1300.)
    ship[lr57<5.8] = ship[lr57<5.8] * (lr57[lr57<5.8]/5.8)
    ship[frzlvl<2400.] = ship[frzlvl<2400.] * (frzlvl[frzlvl<2400.]/2400.)
    return ship