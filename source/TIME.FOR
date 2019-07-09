* Source file TIME.FOR |||||||||||||||||||||||||||||||||||||||||||||||||

      subroutine TmCont(dt,dtMaxW,dtOpt,dMul,dMul2,dtMin,Iter,tPrint,
     !                  tAtm,t,tMax,dtMaxC,ItMin,ItMax,lMinStep,dtInit)
      logical lMinStep
      double precision t,tPrint,tAtm,tMax,tFix

      if(lMinStep) then
        dtMax=amin1(dtMaxW,dtMaxC,dtInit,dtOpt)
        dtOpt=dtMax
        lMinStep=.false.
      else
        dtMax=amin1(dtMaxW,dtMaxC)
      end if
      tFix=dmin1(tPrint,tAtm,tMax)
      if(Iter.le.ItMin.and.(tFix-t).ge.dMul*dtOpt) 
     !  dtOpt=amin1(dtMax,dMul*dtOpt)
      if(Iter.ge.ItMax)
     !  dtOpt=amax1(dtMin,dMul2*dtOpt)
      dt=amin1(dtOpt,sngl(tFix-t))
      iStep=1
      if(dt.gt.0.) iStep=anint(sngl(tFix-t)/dt)
      if(iStep.ge.1.and.iStep.le.10) 
     !  dt=amin1(sngl(tFix-t)/iStep,dtMax)
      if(iStep.eq.1) then
        dt=sngl(tFix-t)
        if(dt-dtMax.gt.dtMin) dt=dt/2.
      end if
      if(dt.le.0.0) dt=dtMin/3.

      return
      end

************************************************************************

      double precision function RTime(iMonth,iDay,iHours,iMins,iSecs,
     !                                i100th)
      
      integer*2 iMonth,iDay,iHours,iMins,iSecs,i100th

      if(iMonth.eq.1.or.iMonth.eq.3.or.iMonth.eq.5.or.iMonth.eq.7.or.
     !   iMonth.eq.8.or.iMonth.eq.10.or.iMonth.eq.12) then
        NoDay=31
      else if(iMonth.eq.4.or.iMonth.eq.6.or.iMonth.eq.9.or.iMonth.eq.11) 
     !                                                then
        NoDay=30
      else if(iMonth.eq.2) then
        NoDay=28
      end if
      nMonth=NoDay*24.*60.*60.
      RTime=nMonth+iDay*24.*60.*60.+iHours*60.*60.+iMins*60.+iSecs+
     !      i100th/100.

      return
      end

************************************************************************

      subroutine SetBC(tMax,tAtm,rTop,rRoot,rBot,hCritA,hBot,hTop,GWL0L,
     !                 TopInF,BotInF,cT,cBot,NS,tTop,tBot,Ampl,lTemp,
     !                 lChem,KodTop,lVarBC,ierr,lMinStep,lMeteo,Prec,
     !                 rSoil,lLAI,rExtinct,lCentrif,CosAlf,xConv,
     !                 tConv,iModel,hTopN,iRootIn,xRoot,WLayer,lLinear,
     !                 lActRSU,SPot)

      logical TopInF,BotInF,lTemp,lChem,lMinStep,lVarBC,lMeteo,lCentrif,
     !        WLayer,lLAI,lVarGr,lLinear(NS),lActRSU
      dimension cBot(NS),cT(NS),cB(11)
      character*3 TheEnd
      double precision tAtm,tMax

      KodTOld=KodTop
      read(31,101) TheEnd
      if(TheEnd.eq.'end') then
        tMax=tAtm
        return
      else
        backspace 31
        if(iRootIn.ne.0) then
          if(.not.lChem.and..not.lTemp) then
            read(31,*,err=901) tAtm,Prec,rSoil,rR,hCA,rB,hB,hT
          else if(lTemp.and..not.lChem) then
            read(31,*,err=901) tAtm,Prec,rSoil,rR,hCA,rB,hB,hT,tTop,
     !                         tBot,Ampl
          else
            read(31,*,err=901) tAtm,Prec,rSoil,rR,hCA,rB,hB,hT,tTop,
     !                         tBot,Ampl,(cT(jj),cB(jj),jj=1,NS)
          end if
        else
          if(.not.lChem.and..not.lTemp) then
            read(31,*,err=901) tAtm,Prec,rSoil,rR,hCA,rB,hB,hT,xRoot
          else if(lTemp.and..not.lChem) then
            read(31,*,err=901) tAtm,Prec,rSoil,rR,hCA,rB,hB,hT,tTop,
     !                         tBot,Ampl,xRoot
          else
            read(31,*,err=901) tAtm,Prec,rSoil,rR,hCA,rB,hB,hT,tTop,
     !                         tBot,Ampl,(cT(jj),cB(jj),jj=1,NS),xRoot
          end if
        end if
        if(lLAI) then
          rLAI=rR
          rPET=rSoil
          SC=1.
          rR=0.
          if(rLAI.gt.0.) 
     !      rR=rPET*amax1(0.,1.-exp(-amax1(rExtinct,0.1)*rLAI))*SC
          rSoil=rPET-rR
        end if
      end if

*     Top of the profile
      if(TopInF) then
        lVarGr=.false. ! Variable gravity field - for Scott Jones
        if(lVarGr) then
          if(CosAlf.ne.Prec) lMinStep=.true.
          CosAlf=Prec
          Prec=0.
        end if
        rTopOld=rTop
        hCritA=-abs(hCA)
        if(lVarBC) then
          rTop=Prec
          if(abs(rTopOld-rTop).gt.abs(rTop)*0.2.and.rTop.lt.0.)
     !      lMinStep=.true.
          KodTop=int(rSoil)
          rSoil=0.
          if(KodTop.eq.-1.and.KodTOld.eq.+1.and.Prec.gt.0.
     !      .and.hNewT.gt.0) hNewT=-0.01*xConv
        else
          if(.not.lMeteo) rTop=abs(rSoil)-abs(Prec)
          if(lMeteo)      rTop=          -abs(Prec)
          if(abs(rTopOld-rTop).gt.abs(rTop)*0.2.and.rTop.lt.0.)
     !      lMinStep=.true.
          if(rTop.gt.0..and.rTopOld.lt.0..and..not.WLayer) then
            xLimit=0.0
            if(iModel.eq.3) xLimit=-0.03*xConv
            if(KodTop.eq.4.or.hTopN.gt.xLimit) then 
              if(iModel.ne.3) xLimit=-0.01*xConv
              hTopN=xLimit
              KodTop=-4
            end if
          end if
        end if
        if(KodTop.eq.3.or.lVarBC) then
          hTop=hT
          if(abs(hTop-hT).gt.abs(hTop)*0.2) lMinStep=.true.
        end if
        rRoot=abs(rR)
      end if

*     Bottom of the profile
      if(BotInF) then
        if(lCentrif) then
          g=9.80665*xConv/tConv/tConv 
          CosAlf=hB*hB/g
          hB=0.
          lMinStep=.true.
        end if
        if(abs(rBot-rB).gt.abs(rBot)*0.2) lMinStep=.true.
        rBot=rB
        if(abs(hBot-hB-GWL0L).gt.abs(hBot)*0.2) lMinStep=.true.
        hBot=hB+GWL0L
      end if

      if(lChem) then
        do 11 jj=1,NS
          if(.not.lLinear(JJ).and.cT(jj).gt.0.) lMinStep=.true.
          cBot(jj)=cB(jj)
          if(lActRSU.and.NS.eq.1) then
            cBot(jj)=0.
            SPot=cB(jj)
          end if
11      continue
      end if
      return

*     Error when reading from an input file 
901   ierr=1
      return

101   format(a3)
      end

************************************************************************

      subroutine SetChemBC(Prec,rSoil,NS,cTop,cT,WLayer,
     !                     hNewT,KodTop,kTopCh)

      logical WLayer
      dimension cTop(NS),cT(NS)

      do 11 jj=1,NS
        if(WLayer.and.hNewT.gt.0.) then
          cTop(jj)=cTop(jj)    ! this is handled in the main program   
        else
          cTop(jj)=cT(jj)
          if(abs(KodTop).eq.4.and.kTopCh.le.0) then
            if(Prec-rSoil.gt.0.) then
              cTop(jj)=cT(jj)*Prec/(Prec-rSoil)
            else if(rSoil.gt.0.) then
              cTop(jj)=0.
            end if
          end if
        end if
11    continue


      return
      end

************************************************************************

      subroutine DailyVar(tConv,t,rRoot,rRootD)

*     Temperature, max at 1. p.m.
*     Radiation, max at noon, hourly values between 0-6 a.m. and 18-24 p.m.
*     represent 1% of daily value, sinusoid in between

      double precision t

      PI=3.141592654
      tPeriod=1.                 ! one day  !24.*60.*60.*tConv
      tDay=sngl(t)/tConv/86400   ! time in day units

c      if(tPeriod.gt.0.) tTopA=tTop+Ampl*sin(2.*PI*sngl(t)/tPeriod-7.*PI/12.)
      tRemainder=amod(tDay,tPeriod)
      if(tRemainder.le.0.264.or.tRemainder.ge.0.736) then
        rRoot=0.24*rRootD
      else
        rRoot=2.75*rRootD*sin(2.*PI*tDay/tPeriod-6.*PI/12.)
      end if
      
      return
      end

************************************************************************

      subroutine SinPrec(t,t1,t2,rPrec,rPrecD)

*     Cosinusoidal distribution of precipitation

      double precision t,t1,t2

      PI=3.141592654
      dt=sngl(t2-t1)
      if(rPrecD.gt.0.) then
        rPrec=rPrecD*(1.+1.*cos(2.*PI*sngl(t-t1)/dt-PI))
      else
        rPrec=0.
      end if

      return
      end

************************************************************************

      subroutine Snow(Prec,dt,Temp,SnowMF,SnowLayer,rEvap,xConv,
     !                lMinStep,cTop,cT,NS)

      logical lMinStep
      dimension cTop(NS),cT(NS)

      PrecOld=Prec
      rEvapOld=rEvap
      Q=1.
      if(SnowLayer.lt.0.001*xConv) then
        if(Temp.lt.-2.0) then
          Q=1.
        else if(Temp.lt.2.0) then
          Q=1.-((Temp+2.)/4.)
        else
          Q=0.
        end if
      end if

      rTop=Prec*(1.-Q)
      SnowF=Prec*Q

      if(Temp.gt.0..and.SnowLayer.gt.0.) then
        SnowMelt=Temp*SnowMF*dt
      else
        SnowMelt=0.
      end if

      SnowLayer=SnowLayer+SnowF*dt-SnowMelt
      if(SnowLayer.lt.0.) then
        SnowMelt=SnowMelt+SnowLayer
        SnowLayer=0.
      else if(SnowLayer.gt.0..and.rEvap.gt.0.) then
        SnowLayerO=SnowLayer
        if(rEvap*dt.lt.SnowLayer) then
          SnowLayer=SnowLayer-rEvap*dt
          rEvap=0.
        else
          rEvap=(rEvap*dt-SnowLayer)/max(dt,1e-8)
          SnowLayer=0.
        end if
      end if

      Prec=rTop+SnowMelt/max(dt,1e-8)
      if(abs(PrecOld-Prec).gt.abs(Prec)*0.2.and.Prec.gt.0.)
     !      lMinStep=.true.

      if(SnowLayer.gt.0.001*xConv) then
        do 11 jj=1,NS
          if((SnowLayer+dt*(PrecOld-rEvapOld)).gt.0.) 
     !       cTop(jj)=(SnowLayer*cTop(jj)+dt*PrecOld*cT(jj))/
     !                (SnowLayer+dt*(PrecOld-rEvapOld))
11      continue
      end if

      return
      end

************************************************************************

      subroutine Meteo(iKod,lMetDaily,lDayVar,t,dt,tInit,tMax,tAtm2,
     !                 tAtmN,tAtm2O,dtMax,rLat,rAlt,ShWRadA,ShWRadB,
     !                 rLWRadA,rLWRadB,rLWRadA1,rLWRadB1,WindHeight,
     !                 TempHeight,iCrop,iLAI,rRoot,xConv,tConv,rGrowth,
     !                 nGrowth,iInterc,rInterc,aInterc,ExcesInt,lEnBal,
     !                 rExtinct,lPrint,lHargr,iRadiation,iSunSh,iRelHum,
     !                 iMetHour,CloudF_Ac,CloudF_Bc,Prec,Precc,rSoil,
     !                 EvapP,TransP,Rns,Rnl,RadTerm,AeroTerm,Rst,ETcomb,
     !                 Rad,RadN,RadO,Wind,WindN,WindO,Albedo,AlbedoN,
     !                 xLAI,xLAIN,xRoot,xRootN,CropHeight,CropHeightN,
     !                 Ampl,tTop,TMaxAN,TMinAN,TMax1,TMaxN,TMaxO,TMin1,
     !                 TMinN,TMinO,TempA,TMaxA,TMaxAO,TMinA,TMinAO,
     !                 SunHours,SunHoursN,SunHoursO,RHMean,RHMeanN,
     !                 RHMeanO,RHMax,RHMaxN,RHMaxO,RHMin,RHMinN,RHMinO,
     !                 RH_A,EaMean,EaMeanN,rTop,ierr)

      implicit real(A-H,L-Z)
      logical lMetDaily,lDayVar,lEnBal,lPrint,lHargr
      dimension rGrowth(1000,5)
      integer nGrowth,iCrop,iRelHum
      double precision t,tInit,tMax,tAtm2,tAtm2O,tAtmN

      ierr=0
      if(iKod.eq.1) then
        if(lMetDaily) then
          dtMax=min(3600.*tConv,dtMax) ! Maximum time step is less than 1 hour
          lDayVar=.false.
          call DailyMet(1,t,dt,tInit,tMax,tAtm2,tAtmN,TMax1,TMaxN,TMaxO,
     !                  TMin1,TMinN,TMinO,TempA,RHMax,RHMaxN,RHMaxO,
     !                  RHMin,RHMinN,RHMinO,RH_A,iRelHum,EaMean,EaMeanN,
     !                  Rad,RadN,Wind,WindN,SunHours,SunHoursN,
     !                  CropHeight,CropHeightN,Albedo,AlbedoN,xLAI,
     !                  xLAIN,xRoot,xRootN,iCrop,xConv,ierr)
          if(ierr.ne.0) goto 901
          call SetDayMet(rLat,rAlt,ShWRadA,ShWRadB,rLWRadA,rLWRadB,
     !                   rLWRadA1,rLWRadB1,WindHeight,TempHeight,iCrop,
     !                   iLAI,CropHeight,Albedo,xLAI,xRoot,tAtm2,rRoot,
     !                   xConv,tConv,iInterc,aInterc,nGrowth,rGrowth,
     !                   ExcesInt,rExtinct,Prec,rSoil,iRadiation,Wind,
     !                   Rad,SunHours,iSunSh,CloudF_Ac,CloudF_Bc,TempA,
     !                   RH_A,t,lEnBal,Rst,ETcomb,EvapP,TransP,Rns,Rnl,
     !                   RadTerm,AeroTerm,rInterc,Precc,ierr)
          if(ierr.eq.3) goto 903
        else 
          call SetMeteo(rLat,rAlt,ShWRadA,ShWRadB,rLWRadA,rLWRadB,
     !                  rLWRadA1,rLWRadB1,WindHeight,TempHeight,iCrop,
     !                  iLAI,CropHeight,Albedo,xLAI,xRoot,tAtm2,rRoot,
     !                  xConv,tConv,iInterc,aInterc,nGrowth,rGrowth,
     !                  ExcesInt,rExtinct,lEnBal,Prec,rSoil,iRadiation,
     !                  TMaxAN,TMinAN,WindN,RHMeanN,RadN,SunHoursN,
     !                  lPrint,iSunSh,iRelHum,tTop,Ampl,CloudF_Ac,
     !                  CloudF_Bc,lHargr,tInit,1,iMetHour,ierr)
          if(ierr.ne.0) goto (901,903) ierr
          call MeteoInt(1,tInit,tAtm2O,tAtm2,Rad,RadO,RadN,TMaxA,TMaxAO,
     !                  TMaxAN,TMinA,TMinAO,TMinAN,Wind,WindO,WindN,
     !                  RHMean,RHMeanO,RHMeanN,SunHours,SunHoursO,
     !                  SunHoursN,lEnBal)
          if(ierr.ne.0) goto (902,903) ierr
        end if
      else if(iKod.eq.2) then
        if(lMetDaily) then
          call DailyMet(2,t,dt,tInit,tMax,tAtm2,tAtmN,TMax1,TMaxN,TMaxO,
     !                  TMin1,TMinN,TMinO,TempA,RHMax,RHMaxN,RHMaxO,
     !                  RHMin,RHMinN,RHMinO,RH_A,iRelHum,EaMean,EaMeanN,
     !                  Rad,RadN,Wind,WindN,SunHours,SunHoursN,
     !                  CropHeight,CropHeightN,Albedo,AlbedoN,xLAI,
     !                  xLAIN,xRoot,xRootN,iCrop,xConv,ierr)		 
          if(ierr.ne.0) goto 901
        else
          call MeteoInt(2,t,tAtm2O,tAtm2,Rad,RadO,RadN,TMaxA,TMaxAO,
     !                  TMaxAN,TMinA,TMinAO,TMinAN,Wind,WindO,WindN,
     !                  RHMean,RHMeanO,RHMeanN,SunHours,SunHoursO,
     !                  SunHoursN,lEnBal)
          call SetMeteo(rLat,rAlt,ShWRadA,ShWRadB,rLWRadA,rLWRadB,
     !                  rLWRadA1,rLWRadB1,WindHeight,TempHeight,iCrop,
     !                  iLAI,CropHeight,Albedo,xLAI,xRoot,tAtm2,rRoot,
     !                  xConv,tConv,iInterc,aInterc,nGrowth,rGrowth,
     !                  ExcesInt,rExtinct,lEnBal,Prec,rSoil,iRadiation,
     !                  TMaxAN,TMinAN,WindN,RHMeanN,RadN,SunHoursN,
     !                  lPrint,iSunSh,iRelHum,tTop,Ampl,CloudF_Ac,
     !                  CloudF_Bc,lHargr,tInit,2,iMetHour,ierr)
          if(ierr.ne.0) goto (901,903) ierr
        end if
      else if(iKod.eq.3) then
        if(lMetDaily) then
          call DailyMet(4,t,dt,tInit,tMax,tAtm2,tAtmN,TMax1,TMaxN,TMaxO,
     !                  TMin1,TMinN,TMInO,TempA,RHMax,RHMaxN,RHMaxO,
     !                  RHMin,RHMinN,RHMinO,RH_A,iRelHum,EaMean,EaMeanN,
     !                  Rad,RadN,Wind,WindN,SunHours,SunHoursN,
     !                  CropHeight,CropHeightN,Albedo,AlbedoN,xLAI,
     !                  xLAIN,xRoot,xRootN,iCrop,xConv,ierr)			     	 
          if(ierr.ne.0) goto 901
          call SetDayMet(rLat,rAlt,ShWRadA,ShWRadB,rLWRadA,rLWRadB,
     !                   rLWRadA1,rLWRadB1,WindHeight,TempHeight,iCrop,
     !                   iLAI,CropHeight,Albedo,xLAI,xRoot,tAtm2,rRoot,
     !                   xConv,tConv,iInterc,aInterc,nGrowth,rGrowth,
     !                   ExcesInt,rExtinct,Prec,rSoil,iRadiation,Wind,
     !                   Rad,SunHours,iSunSh,CloudF_Ac,CloudF_Bc,TempA,
     !                   RH_A,t,lEnBal,Rst,ETcomb,EvapP,TransP,Rns,Rnl,
     !                   RadTerm,AeroTerm,rInterc,Precc,ierr)
          if(ierr.eq.3) goto 903
          rTop=abs(rSoil)-abs(Prec)
        end if
        if(lEnBal.and..not.lMetDaily)
     !    call MeteoInt(3,t,tAtm2O,tAtm2,Rad,RadO,RadN,TMaxA,TMaxAO,
     !                  TMaxAN,TMinA,TMinAO,TMinAN,Wind,WindO,WindN,
     !                  RHMean,RHMeanO,RHMeanN,SunHours,SunHoursO,
     !                  SunHoursN,lEnBal)
      end if
      return

901   ierr=1 !932
      return
902   ierr=2 !913
      return
903   ierr=3 !933
      return

      return
      end

************************************************************************

*     Subroutine reading meteorological input data and calculating potential
*     evapotranspiration using either Penman-Montheith or Hargreaves equations.

      subroutine SetMeteo(Latitude,Altitude,ShortWaveRadA,ShortWaveRadB,
     !                    LongWaveRadA,LongWaveRadB,LongWaveRadA1,
     !                    LongWaveRadB1,WindHeight,TempHeight,iCrop,
     !                    iLAI,CropHeight,Albedo,LAI,xRoot,tAtm,rRoot,
     !                    xConv,tConv,iInterc,aInterc,nGrowth,rGrowth,
     !                    ExcesInt,rExtinct,lEnBal,Prec,rSoil,
     !                    iRadiation,TMax,TMin,Wind_ms,RHMean,Rad,
     !                    SunHours,lPrint,iSunSh,iRelHum,TAver,Ampl,
     !                    CloudF_Ac,CloudF_Bc,lHargr,tInit,iInit,
     !                    iMetHour,ierr)

*     DayNo   - Day number
*     Rad     - Net radiation flux at the surface [MJ/m2/d]
*     TMax    - Maximum temperature [C]
*     TMin    - Minimum temperature [C]
*     RHMean  - Relative humidity [%]
*     Wind_kmd- Average daily wind speed [km/d]
*     Albedo  - Albedo [-]
*     xRoot   - Rooting depth [L]
*     iRadiation - = 0: potential radiation, = 1: solar radiation, = 2: net radiation
*     iInterc - interception
*               =1: uses LAI and Soil Cover (SC)
*               =2; uses only SC

      implicit real(A-H,L-Z)
      logical lEnBal,lPrint,lHargr
      dimension rGrowth(1000,5)
      integer nGrowth
      double precision tInit,tAtm,tAtmOld
      character*3 TheEnd

      PI=3.141592654
      raa=0.
      rc=60.
      SCF=0.

*     Conversion to mm/d from L/T
      rConv=0.001*xConv
      TTConv=24.*60.*60.*tConv

      tAtmOld=tAtm
      read(33,101) TheEnd
      if(TheEnd.eq.'end') then
c        tMax=tAtm
        return
      else
        backspace 33
        if(iCrop.eq.3) then
          read(33,*,err=901) tAtm,Rad,TMax,TMin,RHMean,Wind_kmd,
     !                       SunHours,CropHeight,Albedo,LAI,xRoot
          CropHeight=CropHeight*100./xConv      ! conversion to cm
        else
          read(33,*,err=901) tAtm,Rad,TMax,TMin,RHMean,Wind_kmd,SunHours
        end if
      end if
      DayNo=sngl(tAtm)/TTConv                   ! Conversion to [d]

      if(iInit.eq.1) then                     ! Check time interval of meteo data  
        if((tAtm-tInit).le.0.9999*TTConv) then	 ! short interval 
          iMetHour=1
        else
          iMetHour=0							               ! daily data
        end if
      end if

      if(iCrop.eq.2) then
        i=1000
        j=5
        call Table(nGrowth,rGrowth,i,j,tAtm,CropHeight,Albedo,LAI,	
     !             xRoot)
      end if

      Wind_ms=Wind_kmd/86.4                   ! Conversion, [m/s]
      if(lEnBal) return
      call CropRes(rc,iLAI,LAI,rExtinct,iCrop,CropHeight,SCF)

      TAver=(TMax+TMin)/2.
      if(.not.lHargr) then
        Ampl=(TMax-TMin)/2.
        Ea_TMax=0.6108*exp((17.27*TMax)/(TMax+237.3))           ! Equation 10
        Ea_TMin=0.6108*exp((17.27*TMin)/(TMin+237.3))           ! Equation 10
        if(iRelHum.eq.0) then
*         When average relative humidity is the input
          EaDew=RHMean/(50./Ea_TMin+50./Ea_TMax)                ! Equation 14
        else if(iRelHum.eq.1) then
*         When average daily vapor pressure is the input
          EaDew=RHMean
          RHMean=EaDew*(50./Ea_TMin+50./Ea_TMax)                ! Equation 13
        end if
      end if

*     Net shortwave radiation
      DayNo=mod(DayNo,365.)
      if(iRadiation.ne.2.or.lHargr) then
        call RadGlobal(Ra,Latitude,DayNo,Omega,xx,yy,SC)
        if(lHargr) goto 11
        call Cloudiness(CloudF,iSunSh,SunHours,Omega,LongWaveRadB,
     !                  LongWaveRadA,n_N,Cover,Tt)

        if(iRadiation.eq.0.or.iSunSh.eq.3) then
          if(iRadiation.eq.0) Rad=Ra*(ShortWaveRadA+ShortWaveRadB*n_N) ! Equation 52
          if(iSunSh.eq.3) then
            Tt=Rad/Ra
            Cover=max(0.1,min(1.,2.330-3.330*Tt))
            n_N=1.-Cover
          end if
          if(iRadiation.eq.1.and.iSunSh.eq.3)	then
            Rad_cs=Ra*(ShortWaveRadA+ShortWaveRadB*1.0)         ! solar radiation of clear-sky		 
            if(iMetHour.eq.1) then	  ! Short term interval, calculate daily variation of Rad_cs
              Sum=0.
              do 10 k=1,24
                sine1=xx+yy*cos(2.*PI/24.*(k-12.))
                Sum=Sum+max(sine1,0.)/24.        
10            continue
              HourNo=24.*mod(DayNo,1.)
              sine=xx+yy*cos(2.*PI/24.*(HourNo-12.))     ! Hourly variations eq.12.8 (Basic)
              Rad_csH=max(sine*Rad_cs/Sum,0.)
              if     (Rad_csH.le.0.0001) then
                CloudF=CloudF_Ac*0.6+CloudF_Bc                  ! average Rad/Rad_csH=0.6 in night time 	   
              else if(Rad.ge.Rad_csH) then 				 
                CloudF=CloudF_Ac*1.+CloudF_Bc                   ! Equation 57	   
              else
                CloudF=max(CloudF_Ac*Rad/Rad_csH+CloudF_Bc,0.01) ! Equation 57	        
              end if
            else					                               ! Daily interval
              CloudF=CloudF_Ac*Rad/Rad_cs+CloudF_Bc             ! Equation 57	   
            end if
          end if
        end if
        Rns=(1.-Albedo)*Rad                                     ! Equation 51

*       Calculate net longwave radiation
        call RadLongNet(Rnl,TMax,TMin,LongWaveRadA1,LongWaveRadB1,
     !                  EaDew,CloudF)

*       Net radiation
        Rn=Rns-Rnl                                              ! Equation 50
      else
        Rn=Rad
      end if

*     Calculate aerodynamic and radiation terms of the Penman-Montheith equation
      call Aero(AeroTerm,RadTerm,Rn,CropHeight,WindHeight,TempHeight,
     !          Wind_ms,Altitude,TAver,TMax,TMin,rc,raa,Ea_TMax,Ea_TMin,
     !          EaDew,ierr)
      if(ierr.eq.3) return

*     Evapotranspiration
      ETcomb=amax1(0.,RadTerm+AeroTerm)                         ! Equation 69, 31
      row=1000.   !(ms) water density [kg/m3]
      row=(1.-7.37e-6*(TAver-4.)**2+3.79e-8*(TAver-4.)**3)*1000.  
      ETcomb=ETcomb/row*1000.      ! conversion from [kg/m2/d] to [mm/d]
	
11    if(lHargr) then	 ! Hargreaves Formula
        ETcomb=0.0023*0.408*Ra*(TAver+17.8)*sqrt(TMax-Tmin)	![mm/d], 0.408 - conversin from MJ/m2/d to mm/d
      end if

*     Potential Evaporation and Transpirations [mm/d]
      EvapP=ETcomb*(1.-SCF)
      TransP=ETcomb*SCF

*     Calculate interception
      Prec=Prec/rConv*TTConv ! to mm/d
      Precc=Prec
      call Intercept(rInterc,iInterc,LAI,aInterc,SCF,Prec,ExcesInt,
     !               TransP)

*     conversion from mm/d to L/T
      rRoot=(amax1(TransP-rInterc,0.))*rConv/TTConv
      rSoil=EvapP*rConv/TTConv
      Prec=Prec*rConv/TTConv

      if(DayNo.eq.0) DayNo=365
      if(lPrint) write(43,120) DayNo,ETcomb,EvapP,TransP,Rns,Rnl,
     !                         RadTerm,AeroTerm,Precc,rInterc,ExcesInt
      return

*     Error when reading from an input file 
901   ierr=1
      return

101   format(a3)
120   format(f8.2,10f10.3)
      end

************************************************************************

*     Calculate aerodynamic and radiation terms of the Penman-Montheith equation

      subroutine Aero(AeroTerm,RadTerm,Rn,CropHeight,WindHeight,
     !                TempHeight,Wind_ms,Altitude,TAver,TMax,TMin,rc,
     !                raa,Ea_TMax,Ea_TMin,EaDew,ierr)

      real lambda

      dl0=0.667*CropHeight
      if(dl0.ge.WindHeight.or.dl0.ge.TempHeight) goto 901
      AerDynRes=alog((WindHeight-dl0)/(0.123*CropHeight))*
     !          alog((TempHeight-dl0)/(0.0123*CropHeight))/
     !          0.41**2                                         ! Equation 36,37,38,39
      if(Wind_ms.gt.0.) raa=AerDynRes/Wind_ms                   ! Equation 36
      AeroTCff=0.622*3.486*86400./AerDynRes/1.01                ! Equation 47b
      PAtm=101.3*((293.-0.0065*Altitude)/293.)**5.253           ! Equation 6
      lambda=2.501-0.002361*TAver                               ! Equation 1
      gamma=0.0016286*PAtm/lambda                               ! Equation 4
      gamma1=gamma
      if(raa.gt.0.) gamma1=gamma*(1+rc/raa)                     ! Equation 41
      Dlt=2049.*Ea_TMax/(TMax+237.3)**2+
     !    2049.*Ea_TMin/(TMin+237.3)**2                         ! Equation 3
      dl_dl=Dlt/(Dlt+gamma1)                                    ! Equation 49a
      gm_dl=gamma/(Dlt+gamma1)                                  ! Equation 47a
      EaMean=(Ea_TMax+Ea_TMin)/2.                               ! Equation 11
      AeroTerm=gm_dl*AeroTCff/(TAver+273.)*Wind_ms*(EaMean-EaDew)! Equation 45, 47   

*     Radiation term
      RadTerm=0.
      if(lambda.gt.0.) RadTerm=dl_dl*Rn/lambda                  ! Equation 49

      return

901   ierr=3
      return

      end
      
************************************************************************

*     Calculate global radiation

      subroutine RadGlobal(Ra,Latitude,DayNo,Omega,xx,yy,SC)

      real Lat1,Lat2,Latitude

      PI=3.141592654
      SC=118.08                                             ! Solar Constant MJ/m2/d

      Lat1=(Latitude-aint(Latitude))*(5./3.)+aint(Latitude)
      Lat2=Lat1*PI/180.
      SolDeclin=sin(2.*PI/365.*DayNo-1.39)*0.4093           ! Equation 22
      Omega=acos(-tan(SolDeclin)*tan(Lat2))                 ! Equation 20
      xx=sin(SolDeclin)*sin(Lat2)                           ! Equation 18a
      yy=cos(SolDeclin)*cos(Lat2)                           ! Equation 18a
      dr=(1.+0.033*cos((2.*Pi/365.)*DayNo))                 ! Equation 21
      Ra=SC/PI*dr*(Omega*xx+sin(Omega)*yy)                  ! Equation 18

      return
      end

************************************************************************

*     Calculate cloudiness factor

*     CloudF	- Cloudiness factor [-]
*     Cover	  - Cloud cover fraction 

      subroutine Cloudiness(CloudF,iSunSh,SunHours,Omega,LongWaveRadB,
     !                      LongWaveRadA,n_N,Cover,Tt)

      implicit real(A-H,L-Z)
      PI=3.141592654

      if(iSunSh.eq.0) then                ! sunshine hours
        NN=24./PI*Omega                                     ! Equation 25
        n_N=min(SunHours/NN,1.)                             ! Equation 52a
        CloudF=LongWaveRadB+LongWaveRadA*n_N                ! Equation 59
        Cover=1.-n_N
      else if(iSunSh.eq.1) then           ! cloudiness
        CloudF=SunHours
        n_N=(CloudF-LongWaveRadB)/LongWaveRadA
        Cover=1.-n_N
      else if(iSunSh.eq.2) then           ! transmission coeff
        Tt=SunHours
        Cover=max(0.0001,min(1.,2.330-3.330*Tt))
        n_N=1.-Cover
        CloudF=LongWaveRadB+LongWaveRadA*n_N                ! Equation 59
      end if

      return
      end

************************************************************************

*     Calculate crop canopy resistance, rc

*     SunHours- Bright sunshine hours per day [hr]
*               alternatively, it can be Tt instead of SunHours
*     iLAI    - LAI is calculated from crop height for grass (1), for alfalfa (2), and soil cover (3)
*     CropHeight - Crop Height [cm]
*     LAI     - Leaf Area Index [-]

      subroutine CropRes(rc,iLAI,LAI,rExtinct,iCrop,CropHeight,SCF)

      real LAI
      logical lCrop

      if(iCrop.eq.0.or.CropHeight.le.0) then
        CropHeight=0.1 ! cm
      	lCrop=.false.
      else
      	lCrop=.true.
      end if

      rc=0.
      if(lCrop) then
        if(iLAI.eq.1) then                    ! only clipped grass
          LAI=0.24*CropHeight                                   ! Equation 35
        else if(iLAI.eq.2) then               ! alfalfa
          LAI=1.5*alog(CropHeight)+5.5                          ! Equation 34
        else if(iLAI.eq.3) then               ! SCF
          if(LAI.lt.1.) then
            LAI=-alog(1.-LAI)/amax1(rExtinct,0.1)
          else 
            LAI=10.
          end if
        end if
        if(LAI.gt.0) rc=200./LAI                                ! Equation 32
        if(LAI.gt.0) SCF=amax1(0.,1.-exp(-amax1(rExtinct,0.1)*LAI))
      else
        rc=0.
      end if

      return
      end

************************************************************************

*     Calculate net longwave radiation

      subroutine RadLongNet(Rnl,TMax,TMin,LongWaveRadA1,LongWaveRadB1,
     !                      Ea,CloudF)

      real LongWaveRadA1,LongWaveRadB1

      sigma=0.00000000245*((TMax+273.16)**4+(TMin+273.16)**4)   ! Equation 63a
      Emissivity=LongWaveRadA1+LongWaveRadB1*sqrt(Ea)           ! Equation 60
      Rnl=sigma*CloudF*Emissivity                               ! Equation 63

      return
      end

************************************************************************

*     Calculate interception [mm]

      subroutine Intercept(rInterc,iInterc,LAI,aInterc,SCF,Prec,
     !                     ExcesInt,TransP)

      real LAI

      rInterc=0.
      if(iInterc.gt.0) then
        if(iInterc.eq.1.and.LAI.gt.0.) then
          rInterc=amin1(aInterc*LAI*(1.-1./(1+SCF*Prec/aInterc/LAI)),
     !                  Prec) ! Newly intercepted
	    if(rInterc+ExcesInt.gt.aInterc*LAI) ! can not be more than max interception
     !      rInterc=aInterc*LAI-ExcesInt
        end if
        Prec=amax1(Prec-rInterc,0.)
        rInterc=ExcesInt+rInterc ! Old interception + new interception
        ExcesInt=0.
        if(TransP-rInterc.lt.0.) ExcesInt=-TransP+rInterc
      end if

      return
      end

************************************************************************

*     Interpolate meteo variables in time
*     i = 1 (initialize), i = 2 (move new into old), i = 3 (interpolate)

      subroutine MeteoInt(i,t,tAtm2O,tAtm2,Rad,RadO,RadN,TMaxA,TMaxAO,
     !                    TMaxAN,TMinA,TMinAO,TMinAN,Wind,WindO,WindN,
     !                    RHMean,RHMeanO,RHMeanN,SunHours,SunHoursO,
     !                    SunHoursN,lEnBal)

      double precision t,tAtm2O,tAtm2
      logical lEnBal ! for lEnBal, there is no need to interpolate Rad and SunHours

      if(i.eq.1) then
        tAtm2O=t
        RadO=RadN
        Rad =RadN
        TMaxAO=TMaxAN
        TMaxA =TMaxAN
        TMinAO=TMinAN
        TMinA =TMinAN
        WindO =WindN
        Wind  =WindN
        RHMeanO=RHMeanN
        RHMean =RHMeanN
        SunHoursO=SunHoursN
        SunHours =SunHoursN
      else if(i.eq.2) then
        tAtm2O=tAtm2
        RadO=RadN
        TMaxAO=TMaxAN
        TMinAO=TMinAN
        WindO =WindN
        RHMeanO=RHMeanN
        SunHoursO=SunHoursN
      else if(i.eq.3) then
        Rad  =RadO  +(RadN  -RadO)  *(sngl(t)-tAtm2O)/(tAtm2-tAtm2O)
        if(lEnBal) Rad=RadN
        TMaxA=TMaxAO+(TMaxAN-TMaxAO)*(sngl(t)-tAtm2O)/(tAtm2-tAtm2O)
        TMinA=TMinAO+(TMinAN-TMinAO)*(sngl(t)-tAtm2O)/(tAtm2-tAtm2O)
        Wind =WindO +(WindN -WindO) *(sngl(t)-tAtm2O)/(tAtm2-tAtm2O)
        RHMean=RHMeanO+(RHMeanN-RHMeanO)*(sngl(t)-tAtm2O)/(tAtm2-tAtm2O)
        SunHours=SunHoursO+(SunHoursN-SunHoursO)*(sngl(t)-tAtm2O)/
     !                                           (tAtm2  -tAtm2O)
        if(lEnBal) SunHours=SunHoursN
      end if

      return
      end

************************************************************************

*     Evaluate surface energy balance

      subroutine Evapor(t,TempS,TMaxA,TMinA,Rad,hTop,TempHeight,
     !                  WindHeight,Wind,RHMean,HeatFlux,rTop,Prec,
     !                  tPeriod,Latitude,Albedo,SunHours,ThetaT,xConv,
     !                  tConv,iRadiation,Rns,Rnl,Rn,SensFlux,Evap,Lat,
     !                  Const,iSunSh,r_H,lMetDaily,Rst,TempA,RH_A,
     !                  ShortWaveRadA,ShortWaveRadB,LongWaveRadA,
     !                  LongWaveRadB,iMetHour)

      implicit real(A-H,L-Z)
      real Lat,Latitude,n_N
      double precision t
      logical lMetDaily

*     Rn      - Net radiation [MJ/m2/d]
*     Ra      - Potential radiation [MJ/m2/d]
*     Rad     - Solar radiation flux at the surface [MJ/m2/d]
*     Rns     - Net short wave radiation [MJ/m2/d]
*     Rst     - Daily variated net short wave radiation [MJ/m2/d]
*     Rnl     - Net long wave radiation [MJ/m2/d]
*     Rnlu    - Outgoing long wave radiation [MJ/m2/d]
*     Rnld    - Incoming long wave radiation [MJ/m2/d]
*     SC      - Solar Constant [MJ/m2/d] (=118.08)
*     g       - gravitational acceleration [m/s2] (9.81)
*     xMol    - molecular weight of water [kg/mol] (0.018015)
*     R       - universal gas constant [J/mol/K] (8.314)
*     row     - density of soil water [kg/m3]
*     Lat     - mass latent heat of vaporization of water [J/kg, m2/s2]
*     xLatent - volumetric latent heat of vaporization of water [J/m3,kg/m/s2]
*     TempS   - soil temperature
*     TempA   - atmosphere temperature
*     rovs    - saturated vapor  density [kg/m3]
*     rov     - vapor density [kg/m3]
*     Hr      - relative humidity [-]
*     uu      - friction velocity
*     Wind    - wind speed [m/s]
*     SensFlux - sensible flux [W/m2]
*     HeatFlux - heat flux fo soil [W/m2]
*     Evap    - evaporation flux [kg/s/m2]
*     sigma   - Stephan-Boltzmann constant [5.6697e-8 J/s/m2/K4], [4.899e-9 MJ/d/m2/K4]
*     Tt      - Transmission coefficient (either given or calculated from potential and 
*               measured solar radiation)
*     r_v		  - Aerodynamic resistance to vapor flow [s/m]
*     r_h		  - Aerodynamic resistance to heat flow [s/m] (=r_v)
*     r_s		  - Soil surface resistance [s/m]
*     CloudF	- Cloudiness factor [-]
*     Cover	  - Cloud cover fraction 

      Pi=3.141592654
      g=9.81
      xMol=0.018015
      R=8.314
      Ca=1200.0

*     Hourly variations of atmospheric temperature
      if(.not.lMetDaily) then
        TMean=(TMaxA+TMinA)/2.
        TempA=TMean
        if(tPeriod.gt.0.) TempA=TMean+(TMaxA-TMinA)/2.*
     !                          sin(2.*PI*sngl(t)/tPeriod-7.*PI/12.)
      end if
      TKelvS=TempS+273.15
      TKelvA=TempA+273.15
      rovsA=0.001*exp(31.3716-6014.79/TKelvA-0.00792495*TkelvA)/TKelvA
      if(lMetDaily) RHMean=RH_A
      roa=RHMean/100.*rovsA

*     Net shortwave radiation (Cloud cover fraction to calculate Rnl)
      if(iRadiation.ne.2) then
        DayNo=mod(sngl(t),365.)
        call RadGlobal(Ra,Latitude,DayNo,Omega,xx,yy,SC)
        call Cloudiness(CloudF,iSunSh,SunHours,Omega,LongWaveRadB,
     !                  LongWaveRadA,n_N,Cover,Tt)
        if(iRadiation.eq.0.or.iSunSh.eq.3) then
          if(iRadiation.eq.0) Rad=Ra*(ShortWaveRadA+ShortWaveRadB*n_N) ! Equation 52
          if(iSunSh.eq.3) then         
            if(iMetHour.eq.1) then    ! short term interval radiation data
              HourNo=24.*mod(DayNo,1.)
              sine=xx+yy*cos(2.*PI/24.*(HourNo-12.))
              RaH=max(SC*sine,0.)     ! Hourly variations of extraterrestorial radiation, eq.12.8 (Basic)
              if(RaH.le.0.0001) then  ! night time
                Cover=0.6             ! average value calculated from Rad/Rad_csH=0.6 	   
              else
                Tt=min(Rad/RaH,1.)	        
                Cover=max(0.1,min(1.,2.330-3.330*Tt))
              end if	      
            else                      ! daily interval radiation data 
              Tt=Rad/Ra 		
              Cover=max(0.1,min(1.,2.330-3.330*Tt))
            end if
          end if
        end if
	  if(iRadiation.eq.0.or.iMethour.eq.0)then
          if(iSunSh.ne.2)Tt=Rad/Ra
          HourNo=24.*mod(DayNo,1.)
          sine=xx+yy*cos(2.*Pi/24.*(HourNo-12.))                ! Hourly variations eq.12.8 (Basic)
          RadH=max(SC*Tt*sine,0.)                               ! Equation 12.7
        else
          RadH=Rad
        end if

*       Surface Albedo van Bavel and Hillel (1976)
        Albedo=0.25
        if(ThetaT.gt.0.25) Albedo=0.1
        if(ThetaT.gt.0.1.and.ThetaT.le.0.25) Albedo=0.35-ThetaT

        if(lMetDaily) RadH=Rst	   
        Rns=(1.-Albedo)*RadH                                    ! Equation 51

*       Net longwave radiation
        call RadLongNet1(Rnl,TempA,RHMean,Cover,ThetaT,TkelvA,TkelvS)

        Rn=Rns+Rnl                                              ! Equation 50
      else
        Rn=Rad
      end if

*     Calculate aerodynamic resistance to vapor flow
      call AeroRes(r_v,TempHeight,WindHeight,Wind,TKelvS,TKelvA,Ca,g)

*     Soil surface resistance (van de Griend and Owe, 1994)
     	r_s=10.0
	if(ThetaT.lt.0.15) r_s=10.0*exp(35.63*(0.15-ThetaT))

*     Aerodynamic resistances to heat transfer
      r_h=r_v
      SensFlux=Ca*(TempS-TempA)/r_h

      h=hTop/xConv              ! conversions to m
      Hr=exp(h*xMol*g/R/TkelvS)
c      if(hTop.lt.0.999*hCritA.and.TempS.gt.TempS1) Hr=0.0001
      rovsS=0.001*exp(31.3716-6014.79/TKelvS-0.00792495*TkelvS)/TKelvS
      rov=rovsS*Hr
      Evap=max(0.,(rov-roa)/(r_v+r_s))
      row=(1.-7.37e-6*(TempS-4.)**2+3.79e-8*(TempS-4.)**3)*1000.
      Lat=xLatent(TempA)/row

*     Conversion from [MJ/m2/d] to [J/m2/s,W/m2]
      rConv=1000000./86400.
      HeatFlux=Rn*rConv-SensFlux-Lat*Evap   ! [W/m2]=[MJ/m2/d]*conv-[W/m2]-[J/kg][kg/m2/s]

*     Conversion to HYDRUS units
      Evap1=Evap/row*xConv/tConv            ! [L/T]=[kg/m2/s][m3/kg]*conv
      Const=Lat*row/xConv*tConv/tConv/tConv/tConv  ! [J/kg][kg/m3]*conv=[J/m3]*conv=[kg/m/s2]*conv
      rTop=-Prec+Evap1

*     Unit conversion from W/m2 [kg/s3] to [kg/T3]
      HeatFlux=-HeatFlux/tConv/tConv/tConv      ! Negative is in the profile

      return
      end

************************************************************************

*     Net longwave radiation

      subroutine RadLongNet1(Rnl,TempA,RHMean,Cover,ThetaT,TkelvA,
     !                       TkelvS)

      sigma=4.899e-09
      Es=0.6108*exp((17.27*TempA)/(TempA+237.3))              ! Equation 10
      Ea=RHMean/100.*Es
      Epsi=1.24*(Ea/TkelvA)**(1./7.)                          ! Brutsaert (1975)      
      EpsiA=max(0.,min((1.-0.84*Cover )*Epsi+0.84*Cover,1.))
      EpsiS=min(0.9+0.18*ThetaT,1.)
      Rlu=EpsiS*sigma*TKelvS**4
      Rld=EpsiA*sigma*TKelvA**4
      Rnl=Rld-Rlu

      return
      end

************************************************************************

*     Aerodynamic resistance to vapor flow (Camillo and Gurney, 1986)

      subroutine AeroRes(r_v,TempHeight,WindHeight,Wind,TKelvS,TKelvA,
     !                   Ca,g)

      real Mo

*     rK      - von Karman constant (=0.41)
*     psim    - atmospheric stability factor for momentum
*     psih    - atmospheric stability factor for heat
*     zm		  - roughness parameter for momentum [m]
*     zh		  - roughness parameter for heat transport [m]
*     dl		  - displacement level for heat transport [m] (0)
*     Mo		  - Monin-Obukhov scaling length
*     zeta	  - Unitless height

      rK=0.41
      Pi=3.141592654

      zm=0.001                  ! roughness parameter for momentum [m]
      zh=zm                     ! roughness parameter for heat transport [m]
      dl=0.0                    ! displacement level for heat transport [m]
      THeight=TempHeight/100.   ! conversions to m
      WHeight=WindHeight/100.

      if(Wind.gt.0.) then
        if(abs(TKelvS-TKelvA).lt.0.01) then  ! The atmosphere is neutral
          r_v=1./Wind/rK/rK*(alog((THeight-dl)/zh))*
     !                      (alog((WHeight-dl)/zm))
          goto 12  
        end if
        psim=0.0	! The atmosphere is not neutral
        psih=0.0
        do 11 i=1,6
          uu=Wind*rK/(alog((WHeight-dl+zm)/zm)+psim)
          r_v=1./uu/rK*(alog((THeight-dl+zh)/zh)+psih)
          if(i.eq.1) then
            rvMin= 0.1*r_v
            rvMax=10.0*r_v
          else
            if(r_v.gt.rvMax) then
              r_v=rvMax
              goto 12
            else if(r_v.lt.rvMin) then
              r_v=rvMin
              goto 12
            end if
          end if
          Mo=-Ca*TKelvA*uu*uu*uu/rK/g/(Ca*(TKelvS-TKelvA)/r_v)  ! Monin-Obukhov scaling length
          zeta=(THeight-dl)/Mo	    ! unitless height
          if(zeta.lt.0.) then		    ! unstable
            xx=(max(0.,1.-16*zeta))**0.25
          	psih=-2.*alog((1+xx*xx)/2.)		
            psim=-2.*alog((1+xx)/2.)-alog((1+xx*xx)/2.)+
     !           2.*atan(xx)-pi/2.
          else if(zeta.gt.0) then	  ! stable
            if(zeta.lt.1) then		
              psih=5.*zeta
              psim=psih
            else if(zeta.ge.1) then
              psih=5.
              psim=psih
            end if
          end if
11      continue
      else       ! no wind
        Diff0=2.12e-5
        DiffT=Diff0*(TKelvA/273.15)**2
        r_v=THeight/DiffT
      end if
12    continue

      return
      end

************************************************************************

*     Adjusts evaporation and heat flux based on the difference between 
*     potential and actual evaporation

      subroutine UpdateEnergy(t,vTop,rTop,HeatFl,TempS,Rns,Rnl,Rn,Evap,
     !                        Lat,SensFlux,xConv,tConv,Const,TLevel,
     !                        nPrStep,r_h,dz,rLamb,iTemp,lPrint,
     !                        lMetDaily,TempA,RH_A,Rst)

*     DeltaE    - Difference of Latent heat term [kg/s3]
*     DeltaT    - Approximate temperature change at surface [C]
*     DeltaSens - Difference of Sensible heat term [kg/s3]
*     rLamb	    - Thermal conductivity of soil surface [kg/m/s3/K]
*     dz        - Distance between 1st and 2nd lattice [L]

      integer TLevel
      real Lat
      logical lPrint,lMetDaily
      double precision t 

      DeltaE=0.
      if(rTop.gt.0..and.rTop.gt.vTop) then
        row=(1.-7.37e-6*(TempS-4.)**2+3.79e-8*(TempS-4.)**3)*1000.
        DeltaE=(rTop-amax1(0.,vTop))*Const                   ! [kg/T3]
        DeltaE=DeltaE*tConv*tConv*tConv                      ! to [kg/s3]
        if(DeltaE.gt.0..and.abs(vTop).gt.0.) then
          Ca=1200.0     ! [J/m3/K]
          Cw=4180000.0
          rLamb=rLamb/xConv*tConv*tConv*tConv                ! to [kg/m/s3/K]
          DeltaT=DeltaE/(Ca/r_h+rLamb/(dz/xConv)-Cw*(vTop*tConv/xConv))   ! [C]
          DeltaSens=Ca*DeltaT/r_h                            ! [kg/s3]
          SensFlux=SensFlux+DeltaSens                        ! [kg/s3]
          HeatFl=HeatFl-(DeltaE-DeltaSens)/tConv/tConv/tConv ! [kg/T3]
          Evap=Evap-(rTop-amax1(0.,vTop))/xConv*tConv*row
        end if
      end if

*     Conversion from [J/m2/s,W/m2] to [MJ/m2/d]
      rConv=1000000./86400.
      Balance=Rn-SensFlux/rConv-Lat*Evap/rConv+
     !        HeatFl/rConv*tConv*tConv*tConv
      if(lPrint.and.abs(float((TLevel+nPrStep-1)/nPrStep)-
     ! (TLevel+nPrStep-1)/float(nPrStep)).lt.0.0001.and.iTemp.eq.1) then
        if(.not.lMetDaily) then
          write(43,110) t,Rns,Rnl,Rn,SensFlux/rConv,Lat*Evap/rConv,
     !                  HeatFl/rConv*tConv*tConv*tConv,Balance
        else
          write(43,120) t,Rns,Rnl,Rn,SensFlux/rConv,Lat*Evap/rConv,
     !                HeatFl/rConv*tConv*tConv*tConv,Balance,TempA,RH_A,
     !                Rst
        end if
      end if
      return

110   format(12e13.5)
120   format(12e13.5)						
      end

************************************************************************

*     Read Meteo.in and generate daily variations of air temperature
*     and relative humidity (April 2008 Masaru Sakai)  

      subroutine DailyMet(iKod,t,dt,tInit,tEnd,tAtm,tAtmN,TMax,TMaxN,
     !                    TMaxO,TMin,TMinN,TMinO,TempA,RHMax,RHMaxN,
     !                    RHMaxO,RHMin,RHMinN,RHMinO,RH_A,iRelHum,
     !                    EaMean,EaMeanN,Rad,RadN,Wind_ms,Wind_msN,
     !                    SunHours,SunHoursN,CropHeight,CropHeightN,
     !                    Albedo,AlbedoN,LAI,LAIN,xRoot,xRootN,iCrop,
     !                    xConv,ierr)
      implicit real(A-H,L-Z)
      double precision t,tInit,tEnd,tAtm,tAtmN
      integer iDaily,iCrop,iRelHum

*     iDaily    - =1: initial, =2: update, =3: tFinal-1 (last day), =4: 0:00<t<24:00
*     tInit     - Initial time
*     tEnd	    - Final time
*     tAtm,N    - Current time, next DOY
*     TMax,N(O)	- Max temperature at current, next, and previous DOY
*     TMin,N(O)	- Min temperature at current, next, and previous DOY
*     TempA		  - Air temperature at Current Time Step
*     RHMax,N(O)- Max relative humidity at current, next, previous DOY
*     RHMin,N(O)- Min relative humidity at current, next, previous DOY
*     RH_A		  - Relative humidity at the current time step
*     EaMean,N  - Daily averaged actual vapor pressure
*     Rad,N     - Radiation data at the current and next DOY
*     Wind_ms,N	- Wind speed (m/s) at the current and next DOY
*     SunHours,N- Sunshine hours at the current and next DOY
*     CropHeight,N - Crop height at the current and next DOY
*     Albedo,N  - Albedo at the current and next DOY
*     LAI,N		  - LAI at the current and next DOY
*     xRoot,N	  - Root depth  at the current and next DOY

*     Read meteorological data from Meteo.In
*     iDaily=1: Initialy read 1st and 2nd day's data
*     iDaily=2: In the middle of calculations, update data and then read the next day's data
*     iDaily=3: At the tFinal-1 (last day) or 0:00<t<24:00 (iDaily=4), do not read data
*     At the end of calculations, skip this subroutine

      iDaily=iKod
      if(dabs(tEnd-t).le.0.001*dt) return
      if(dabs((tEnd-1.)-t).le.0.001*dt) iDaily=3
      if(iDaily.eq.1) then
        if(iCrop.eq.3) then
          read(33,*,err=901) tAtm,Rad,TMax,TMin,RHMean,Wind_kmd,
     !                       SunHours,CropHeight,Albedo,LAI,xRoot
          CropHeight=CropHeight*100./xConv  ! conversion to cm
        else
          read(33,*,err=901) tAtm,Rad,TMax,TMin,RHMean,Wind_kmd,SunHours
        end if
        call Humidity(TMax,TMin,RHMean,iRelHum,RHMax,RHMin,EaMean)  ! Estimate RHMax,RHmin
        Wind_ms=Wind_kmd/86.4               ! Conversion to m/s
      end if

*     Update Data
      if(iDaily.eq.2.or.iDaily.eq.3) then
        tAtm=tAtmN
        TMaxO=TMax
        TMax=TMaxN
        TMinO=Tmin
        Tmin=TminN
        RHMaxO=RHMax
        RHMax=RHMaxN
        RHMinO=RHMin
        RHMin=RHMinN
        EaMean=EaMeanN
        Rad=RadN
        Wind_ms=Wind_msN
        SunHours=SunHoursN
        if(iCrop.eq.3)then
          CropHeight=CropHeightN
          Albedo=AlbedoN
          LAI=LAIN
          xRoot=xRootN
        end if
      end if

*     Read Next DOY's data
      if(iDaily.le.2) then
        if(iCrop.eq.3) then
          read(33,*,err=901) tAtmN,RadN,TMaxN,TMinN,RHMeanN,Wind_kmdN,
     !                       SunHoursN,CropHeightN,AlbedoN,LAIN,xRootN
          CropHeightN=CropHeightN*100./xConv  ! conversion to cm
        else
          read(33,*,err=901) tAtmN,RadN,TMaxN,TMinN,RHMeanN,Wind_kmdN,
     !                       SunHoursN
        end if
        call Humidity(TMaxN,TMinN,RHMeanN,iRelHum,RHMaxN,RHMinN,
     !                EaMeanN)                ! Estimate RHMax,RHmin
        Wind_msN=Wind_kmdN/86.4               ! Conversion to m/s
      end if

      if(iDaily.eq.2) return                  ! do not calculate temp and RH
      	
*     Calculate temperature, TempA, at the current time step
      call DailyTemp(t,dt,tInit,tEnd,tAtm,TMax,TMaxN,TMaxO,TMin,TMinN,
     !               TMinO,TempA)    
*     Calculate relative humidity, RH_A, at the current time step
      call DailyTemp(t,dt,tInit,tEnd,tAtm,RHMin,RHMinN,RHMinO,RHMax,
     !               RHMaxN,RHMaxO,RH_A) 
      return

*     Error  
901   ierr=1

      return
      end

************************************************************************

*     Calculate temperature (or relative humidity) from daily max and min values

      subroutine DailyTemp(t,dt,tInit,tEnd,tAtm,TMax,TMaxN,TMaxO,TMin,
     !                   	 TMinN,TMinO,TempA)

      implicit real(A-H,L-Z)
      double precision t,tInit,tEnd,tAtm

*     maxTime	- The time of maximum temperature 13:00
*     minTime	- The time of minimum temperature 1:00

      PI=3.141592654
      maxTime=13./24.
      minTime=1./24.

*     Calculate maximum temperature, TMaxA, for current time
      if(t.lt.tInit+maxTime) then		                        ! Initial - 13:00 of 1st day
        TMaxA=TMax
      else if(t.ge.tEnd-1.+maxTime)	then                    ! 13:00 of last day - final time
        TMaxA=TMax
      else if((t.ge.tAtm-1.+maxTime).and.((t.lt.tAtm)
     !               .or.(abs(t-tAtm).le.dt*0.001))) then   ! 13:00-24:00
        TMaxA=(TMaxN-TMax)*(t-(tAtm-1.+maxTime))+TMax
      else if(((t.gt.tAtm-1.).or.(abs(t-(tAtm-1.))).le.dt*0.001)
     !               .and.(t.lt.tAtm-1.+maxTime)) then      ! 0:00-13:00
        TMaxA=(TMax-TMaxO)*(t-(tAtm-2.+maxTime))+TMaxO
      end if

*     Calculate minimum temperature, TMinA, for current time
      if(t.lt.tInit+minTime) then		                        ! Initial - 1:00 of 1st day
        TMinA=TMin
      else if(t.ge.tEnd-1.+minTime)	then                    ! 1:00 of last day - final time
        TMinA=TMin
      else if((t.ge.tAtm-1.+minTime).and.((t.lt.tAtm)
     !                .or.(abs(t-tAtm).le.dt*0.001))) then  ! 1:00-24:00
        TMinA=(TMinN-TMin)*(t-(tAtm-1.+minTime))+TMin
      else if(((t.gt.tAtm-1.).or.(abs(t-(tAtm-1.))).le.dt*0.001)
     !              	.and.(t.lt.tAtm-1.+minTime)) then     ! 0:00-1:00
        TMinA=(TMin-TMinO)*(t-(tAtm-2.+minTime))+TMinO
      end if

*     Calculate current temperature
      TempA=(TMaxA+TminA)/2.+(TMaxA-TminA)/2.*Cos(2.*PI*(t-maxTime))

      return
      end

************************************************************************

*     Calculate max and min relative humidity from daily max and min 
*     temperatures and average RH 

      subroutine Humidity(TMax,TMin,RHMean,iRelHum,RHMax,RHMin,EaMean)

      implicit real(A-H,L-Z)
      integer iRelHum

*     iRelHum	- =0: Input RH, =1: Input Vapor pressure for RHMean
*     EaMean	- Average Saturation vapor pressure
*     Es_TMax	- Maximum Saturation vapor pressure
*     Es_TMin	- Minimum Saturation vapor pressure
*     RHMax   - Estimated Max RH
*     RHMin   - Estimated Min RH

      Es_TMax=0.6108*exp((17.27*TMax)/(TMax+237.3))             ! Equation 10
      Es_TMin=0.6108*exp((17.27*TMin)/(TMin+237.3))             ! Equation 10

      if(iRelHum.eq.0) then
        RHMin=2.*RHMean*Es_TMin/(Es_TMax+Es_TMin)	
        RHMax=2.*RHMean*Es_TMax/(Es_TMax+Es_TMin)
        EaMean=RHMin/100.*Es_TMax	
      else
        EaMean=RHMean
        RHMin=EaMean/Es_TMax*100.
        RHMax=EaMean/Es_TMin*100.
      end if    

      if(RHMax.ge.100.) then
        if(iRelHum.eq.1) RHMean=EaMean*(50./Es_TMin+50./Es_TMax)   ! Equation 13
        RHMin=RHMean-(100.-RHMean)
        RHMax=100.        
      end if

      return
      end

************************************************************************

      subroutine SetDayMet(Latitude,Altitude,ShortWaveRadA,
     !                     ShortWaveRadB,LongWaveRadA,LongWaveRadB,
     !                     LongWaveRadA1,LongWaveRadB1,WindHeight,
     !                     TempHeight,iCrop,iLAI,CropHeight,Albedo,LAI,
     !                     xRoot,tAtm,rRoot,xConv,tConv,iInterc,aInterc,
     !                     nGrowth,rGrowth,ExcesInt,rExtinct,Prec,rSoil,
     !                     iRadiation,Wind_ms,Rad,SunHours,iSunSh,
     !                     CloudF_Ac,CloudF_Bc,TempA,RH_A,t,lEnBal,Rst,
     !                     ETcomb,EvapP,TransP,Rns,Rnl,RadTerm,AeroTerm,
     !                     rInterc,Precc,ierr)

*     tAtm	- Current DOY
*     TempA - Air tempeature [C]
*     RH_A  - Air relative humidity [%]
*     Ea    - Daily averaged vapor pressure [kPa]
*     Es		- Saturation vapor pressure at current time step [kPa]
*     Rad   - Incoming solar radiation [MJ/m2/d]
*     Rst   - Solar radiation at current time step [MJ/m2/d]

      implicit real(A-H,L-Z)
      logical lEnBal
      dimension rGrowth(1000,5)
      integer nGrowth
      double precision t,tAtm

      PI=3.141592654
      raa=0.
      rc=60.
      SCF=0.

*     Conversion to mm/d from L/T
      rConv=0.001*xConv
      TTConv=24.*60.*60.*tConv

      if(iCrop.eq.2) then
        i=1000
        j=5
        call Table(nGrowth,rGrowth,i,j,t,CropHeight,Albedo,LAI,	
     !             xRoot)
      end if

      call CropRes(rc,iLAI,LAI,rExtinct,iCrop,CropHeight,SCF)

*     Calculate vapor pressure
      Es=0.6108*exp((17.27*TempA)/(TempA+237.3))		            ! Equation 10
      Ea=Es*RH_A/100.

*     Daily Variated Net shortwave radiation, Rst
      DayNo=mod(sngl(tAtm),365.)
      if(iRadiation.ne.2) then
        call RadGlobal(Ra,Latitude,DayNo,Omega,xx,yy,SC)

*       Calculate daily average incoming shortwave radiation
        if(iRadiation.eq.0) then
          call Cloudiness(CloudF,iSunSh,SunHours,Omega,LongWaveRadB,
     !                    LongWaveRadA,n_N,Cover,Tt)
          Rad=Ra*(ShortWaveRadA+ShortWaveRadB*n_N)            ! Equation 52
        end if

*       Average of Rst of every hour should be daily averaged Rs value 
        Sum=0.
        do 10 k=1,24
          sine1=xx+yy*cos(2.*PI/24.*(k-12.))
          Sum=Sum+max(sine1,0.)/24.        
10      continue
        t2=(t-DayNo)*24.
        sine=xx+yy*cos(2.*PI/24.*(t2-12.))
        Rst=max(sine*Rad/Sum,0.)
        Rns=(1.-Albedo)*Rst                                     ! Equation 51

*       Calculate Cloudiness Factor for net long wave radiation       
        call Cloudiness(CloudF,iSunSh,SunHours,Omega,LongWaveRadB,
     !                  LongWaveRadA,n_N,Cover,Tt)
        if(iRadiation.eq.1.and.iSunSh.eq.3) then	      ! measured solar radiation
          Rad_cs=Ra*(ShortWaveRadA+ShortWaveRadB*1.0)   ! solar radiation of clear-sky		 
          CloudF=CloudF_Ac*Rad/Rad_cs+CloudF_Bc                 ! Equation 57
        end if
*       Net longwave radiation
        call RadLongNet(Rnl,TempA,TempA,LongWaveRadA1,LongWaveRadB1,
     !                  Ea,CloudF)

*       Net radiation
        Rn=Rns-Rnl                                              ! Equation 50
      else	  ! Use Measured Net Radiation data
        Rn=Rad
      end if

      if(lEnBal) return

*     Calculate aerodynamic and radiation terms of the Penman-Montheith equation
      call Aero(AeroTerm,RadTerm,Rn,CropHeight,WindHeight,TempHeight,
     !          Wind_ms,Altitude,TempA,TempA,TempA,rc,raa,Es,Es,Ea,ierr)
      if(ierr.eq.3) return

*     Evapotranspiration
      ETcomb=amax1(0.,RadTerm+AeroTerm)                         ! Equation 69, 31
      row=1000.   !(ms) water density [kg/m3]
      row=(1.-7.37e-6*(TempA-4.)**2+3.79e-8*(TempA-4.)**3)*1000.  
      ETcomb=ETcomb/row*1000.      !(ms) conversion [kg/m2/d] to [mm/d]
	
*     Potential Evaporation and Transpirations [mm/d]
      EvapP=ETcomb*(1.-SCF)
      TransP=ETcomb*SCF

*     Calculate interception
      Prec=Prec/rConv*TTConv ! to mm/d
      Precc=Prec
      call Intercept(rInterc,iInterc,LAI,aInterc,SCF,Prec,ExcesInt,
     !               TransP)

*     conversion from mm/d to L/T
      rRoot=(amax1(TransP-rInterc,0.))*rConv/TTConv
      rSoil=EvapP*rConv/TTConv
      Prec=Prec*rConv/TTConv

      return

!120   format(f8.2,14f10.3)

      end

************************************************************************

*	    Output for Penman-Monteith with Daily variated meteorological information

      subroutine DayMeteoOut(t,ETcomb,EvapP,TransP,Rns,Rnl,RadTerm,
     !                       AeroTerm,Precc,rInterc,ExcesInt,TempA,
     !                       RH_A,Rst,lPrint)

      logical lPrint
      double precision t

      if(lPrint) write(43,120) t,ETcomb,EvapP,TransP,Rns,Rnl,RadTerm,
     !                         AeroTerm,Precc,rInterc,ExcesInt,TempA,
     !                         RH_A,Rst
      return
120   format(f8.2,14f10.3)

      end

************************************************************************

      subroutine Table(n1,rTable,n,m,t,a,b,c,d)

      dimension rTable(n,m),y(4)
      double precision t

      if(t.le.rTable(1,1)) then
        do 11 i=1,m-1
          y(i)=rTable(1,1+i)
11      continue
      else if(t.ge.rTable(n1,1)) then
        do 12 i=1,m-1
          y(i)=rTable(n1,1+i)
12      continue
      else
        do 14 j=2,n1
          if(t.gt.rTable(j-1,1).and.t.le.rTable(j,1)) then
            dt=t-rTable(j-1,1)
            do 13 i=1,m-1
              y(i)=rTable(j-1,1+i)+(rTable(j,1+i)-rTable(j-1,1+i))*dt/
     !                             (rTable(j,1)  -rTable(j-1,1))
13          continue
          end if
14      continue
      end if
      a=y(1)
      b=y(2)
      c=y(3)
      d=y(4)

      return
      end

* ||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||
