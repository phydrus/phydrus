*||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||*
*                                                                      *
*     HYDRUS   - Numerical model of one-dimensional variably saturated *
*                water flow, heat transport, and transport of solutes  *
*                involved in sequential first-order decay reactions    *
*                version 7.0                                           *
*                                                                      *
*                Version coupled with user-friendly interface for      *
*                MS Windows environment                                *
*                                                                      *
*     Designed by J.Simunek, M. Sejna, and M. Th. van Genuchten (1996) *
*                                                                      *
*                                        Last modified: October, 2009  *
*                                                                      *
*                                                                      *
*     COPYRIGHT (c) 2005-9, Jiri Simunek, PC-Progress                  *
*                                                                      *
*||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||*

      program HYDRUS

*#######################################################################
*     Procedura pro osetreni nasledujicich signalu: 
*         - signal CTRL+C            
*         - signal CTRK+BREAK   
*         - abnormal termination 
*         - floating point error 
*-----------------------------------------------------------------------
*      use MSFLIB
*
*      interface
*        function h_sig (signum)
*          !MS$ATTRIBUTES C :: h_sig
*          integer(4) h_sig
*          integer(2) signum
*        end function
*      end interface
*
*     interface
*        function hand_fpe (sigid, except)
*           !MS$ATTRIBUTES C :: hand_fpe
*          integer(4) hand_fpe
*          integer(2)  sigid, except
*        end function
*      end interface
*#######################################################################

      parameter (NumNPD=1001,
     !           NMatD =20,
     !           NTabD =100,
     !           NObsD =100,
     !           NSD   =11,
     !           NUnitD=7,
     !           NPD   =30000)

      integer PLevel,Alevel,TLevel,err
      logical SinkF,WLayer,qGWLF,TopInF,ShortO,lWat,lChem,lTemp,ConvgF,
     !        FreeD,SeepF,BotInF,AtmBC,lRoot,lUpW,lWTDep,lSTDep,lEquil,
     !        lLinear(NSD),lArtD,lMoSink,lSolRed,lSolAdd,lScreen,lTable,
     !        lMinStep,qDrain,lMsSink,lTort,lInitW,lVarBC,lPrint,lFiltr,
     !        lMobIm(NMatD),lMeteo,lBact,lVapor,lDayVar,lEnter,lEnBal,
     !        lVaporOut,lExtrap,lPrintD,lLAI,lDensity,lSnow,lCentrif,
     !        lEqInit,lSinPrec,lDualNEq,lFlux,lHargr,lMetDaily,lEnd,
     !        lMassIni,lStopConv,lActRSU,lOmegaW,lFluxOut
      double precision P,R,S,Q,RTime1,t,tInit,tOld,RTime,tMax,tAtm,
     !                 tPrint1,tAtm1,tAtm2,tAtmOld,tAtmN,tAtm2O,
     !                 TPrint,tPrintInt
      character cFileName*260,cDataPath*260
      integer*2 status
      integer*2 i2,iYear,iMonth,iDay,iHours,iMins,iSecs,i100th

      dimension x(NumNPD),hNew(NumNPD),hOld(NumNPD),hTemp(NumNPD),
     !  MatNum(NumNPD),ParD(11,NMatD),TPrint(NPD),Sink(NumNPD),
     !  POptm(NMatD),Beta(NumNPD),LayNum(NumNPD),CumQ(12),SubVol(10),
     !  ConTab(NTabD,NMatD),CapTab(NTabD,NMatD),TheTab(NTabD,NMatD),
     !  Con(NumNPD),Cap(NumNPD),ConSat(NMatD),ths(NMatD),P(NumNPD),
     !  R(NumNPD),S(NumNPD),hSat(NMatD),ParW(11,NMatD),Area(10),
     !  ChPar(NSD*16+4,NMatD),Conc(NSD,NumNPD),vNew(NumNPD),
     !  vOld(NumNPD),Disp(NumNPD),Retard(NumNPD),g0(NumNPD),g1(NumNPD),
     !  Q(NumNPD),wc(NumNPD),ThNew(NumNPD),ThOld(NumNPD),WatIn(NumNPD),
     !  SolIn(NumNPD),Node(NObsD),iUnit(NUnitD),cvCh0(NSD),cvCh1(NSD),
     !  cTop(NSD),cBot(NSD),cRoot(NSD),CumCh(10,NSD),cVOlI(NSD),
     !  cCumA(NSD),cCumT(NSD),cvTop(NSD),cvBot(NSD),cNew(NumNPD),
     !  cTemp(NumNPD),cPrevO(NumNPD),TDep(NSD*16+4),TPar(10,NMatD),
     !  TempO(NumNPD),TempN(NumNPD),SorbN(NumNPD),Sorb(NSD,NumNPD),
     !  q0(NumNPD),q1(NumNPD),aOsm(NSD),thr(NMatD),sSink(NumNPD),
     !  Ah(NumNPD),AK(NumNPD),ATh(NumNPD),AhW(NMatD),AThW(NMatD),
     !  AKW(NMatD),Kappa(NumNPD),AThS(NumNPD),ThRR(NumNPD),cvChR(NSD),
     !  ConR(NumNPD),ConO(NumNPD),AKS(NumNPD),KappaO(NumNPD),
     !  cRootMax(NSD),cvChIm(NSD),hTab(NTabD,NMatD),NTab(NMatD),
     !  rGrowth(1000,5),ConLT(NumNPD),ConVT(NumNPD),ConVh(NumNPD),
     !  ThEq(NumNPD),ThVOld(NumNPD),ThVNew(NumNPD),vVOld(NumNPD),
     !  vVNew(NumNPD),Sorb2(NSD,NumNPD),SorbN2(NumNPD),
     !  ThNewIm(NumNPD),ThOldIm(NumNPD),SinkIm(NumNPD),STrans(NumNPD),
     !  DMoist(NMatD,NSD,13,6),WDep(2+NMatD,NSD*9),cT(NSD)

      data lStopConv,iNonConv /.true. , 0 /

*-----------------------------------------------------------------------
*      i4ret = SIGNALQQ(SIG$ABORT, h_sig)
*      i4ret = SIGNALQQ(SIG$BREAK, h_sig)
*      i4ret = SIGNALQQ(SIG$INT,   h_sig)
*      i4ret = SIGNALQQ(SIG$FPE,   hand_fpe)
*-----------------------------------------------------------------------

*      iCount = NARGS()
      iCount = IARGC()
      if(iCount.gt.1.or.iCount.eq.0) then
        i2=1
        call GETARG(i2, cDataPath)
      else
        cFileName = 'LEVEL_01.DIR'
        open(10, file=cFileName, status='old', err=901)
        read(10, 101, err=904) cDataPath
        close (10)
      end if
      iLengthPath = Len_Trim(cDataPath)
      if(iLengthPath.gt.260-13) goto 930
      cFileName = cDataPath(1:iLengthPath)//'/Options.in'
      open(35,file=cFileName, status='old',err=998)
998   continue

*     Initialization
      call Init(CosAlf,NTab,ItCum,TLevel,ALevel,PLevel,hRoot,vRoot,
     !          IterW,IterC,dtMaxC,wCumT,wCumA,err,lVarBC,NSD,cRoot,
     !          cCumT,cCumA,NumNPD,Sink,wc,CumQ,lMeteo,lBact,lVapor,
     !          lEnBal,lDayVar,lEnter,lFiltr,TauW,nPrStep,nTabMod,
     !          iDualPor,dtMaxT,lExtrap,lPrintD,lLAI,rExtinct,lDensity,
     !          ExcesInt,lMinstep,lPrint,lSnow,SnowMF,SnowLayer,
     !          cTemp,iSunSh,iRelHum,xRoot,lCentrif,Radius,GWL0L,
     !          tPrintInt,iTort,iEnhanc,hSeep,lEqInit,lSinPrec,
     !          iMoistDep,OmegaC,WTransf,lDualNEq,lFlux,iCrop,iRootIn,
     !          lMassIni,lMetDaily,Sorb2,lActRSU,OmegaS,SPot,lOmegaW,
     !          OmegaW,lEnd,lVaporOut,lFluxOut)
      data iUnit /50,70,71,75,76,77,78/

      ! Force lEnter to be False for the Python implementation
      lEnter = .false.

      cFileName = cDataPath(1:iLengthPath)//'/SELECTOR.IN'
      open(30,file=cFileName, status='old',err=901)
      cFileName = cDataPath(1:iLengthPath)//'/PROFILE.DAT'
      open(32,file=cFileName, status='old',err=901)
      cFileName = cDataPath(1:iLengthPath)//'/I_CHECK.OUT'
      open(50,file=cFileName, status='unknown',err=902)
      if(lPrint) then
        cFileName = cDataPath(1:iLengthPath)//'/RUN_INF.OUT'
        open(70,file=cFileName, status='unknown',err=902)
        cFileName = cDataPath(1:iLengthPath)//'/T_LEVEL.OUT'
        open(71,file=cFileName, status='unknown',err=902)
        cFileName = cDataPath(1:iLengthPath)//'/NOD_INF.OUT'
        open(75,file=cFileName, status='unknown',err=902)
        cFileName = cDataPath(1:iLengthPath)//'/BALANCE.OUT'
        open(76,file=cFileName, status='unknown',err=902)
        cFileName = cDataPath(1:iLengthPath)//'/OBS_NODE.OUT'
        open(77,file=cFileName, status='unknown',err=902)
        cFileName = cDataPath(1:iLengthPath)//'/PROFILE.OUT'
        open(78,file=cFileName, status='unknown',err=902)
      end if

*     Read input data --------------------------------------------------

      call BasInf(CosAlf,MaxIt,TolTh,TolH,TopInF,BotInF,ShortO,lWat,
     !            lChem,SinkF,WLayer,qGWLF,FreeD,SeepF,AtmBC,KodTop,
     !            KodBot,rTop,rRoot,rBot,hCritS,hCritA,GWL0L,Aqh,Bqh,
     !            kTOld,kBOld,NUnitD,iUnit,NMat,NMatD,NLay,lRoot,lTemp,
     !            lWTDep,lEquil,lScreen,qDrain,zBotDr,BaseGW,rSpacing,
     !            iPosDr,rKhTop,rKhBot,rKvTop,rKvBot,Entres,WetPer,
     !            zInTF,GeoFac,lInitW,lVarBC,xConv,tConv,lMeteo,lVapor,
     !            iVer,lPrint,lCentrif,lSnow,hSeep,lFlux,lActRSU,err)
      if(err.ne.0) goto (905,906,903,916) err

      if(lPrint) then
        if(lFluxOut) then
          cFileName = cDataPath(1:iLengthPath)//'/T_Level1.OUT'
          open(44,file=cFileName, status='unknown',err=902)
        end if
        if(lVaporOut.and.lVapor) then
          cFileName = cDataPath(1:iLengthPath)//'/NOD_INF_V.OUT'
          open(45,file=cFileName, status='unknown',err=902)
        end if
      end if
      if(TopInF.or.BotInF.or.AtmBC) then
        cFileName = cDataPath(1:iLengthPath)//'/ATMOSPH.IN'
        open(31,file=cFileName, status='old',err=901)
        if(lPrint) then
          cFileName = cDataPath(1:iLengthPath)//'/A_LEVEL.OUT'
          open(72,file=cFileName, status='unknown',err=902)
        end if
        if(lMeteo) then
          cFileName = cDataPath(1:iLengthPath)//'/Meteo.in'
          open(33,file=cFileName, status='old',err=902)
          if(lPrint) then
            cFileName = cDataPath(1:iLengthPath)//'/Meteo.out'
            open(43,file=cFileName, status='unknown',err=902)
          end if
        end if
      end if

      call NodInf(NumNPD,NumNP,NObsD,NObs,hTop,hBot,x,hNew,hOld,
     !            MatNum,hTemp,LayNum,Beta,Ah,AK,ATh,Conc,Sorb,
     !            TempN,TempO,Node,NSD,NS,xSurf,lChem,lTemp,lEquil,
     !            lScreen,lBact,Sorb2,err,lPrint,lDualNEq,lFlux)
      if(err.ne.0) goto (912,903,914,915,929) err
      call MatIn (NMat,ParD,ParW,hTab(1,1),hTab(NTab(1),1),lScreen,err,
     !            NumNP,Ah,iHyst,AhW,AThW,AKW,MatNum,hNew,Kappa,AThS,
     !            ThRR,ConR,AKS,KappaO,iModel,xConv,lTable,IKappa,
     !            nTabMod,iDualPor)
      if(err.ne.0) goto (906,903) err
      if(iHyst.gt.0.and.iHyst.ne.3)
     !call HysterIn(NumNP,NMat,hOld,MatNum,ParD,ParW,ThNew,ThOld,Kappa,
     !              AThS,ThRR,ConO,ConR,AKS,KappaO,Ah,AK,iHyst,iModel,
     !              cDataPath)      
      if(iModel.eq.nTabMod) then
        cFileName = cDataPath(1:iLengthPath)//'/Mater.in'
        open(36,file=cFileName, status='old',err=901)
      end if
      if(lInitW.and.iModel.lt.nTabMod) then
        call InitW(NumNP,NMat,Matnum,Kappa,hNew,hOld,hTemp,ParD,ParW,
     !             iModel,hTop,hBot,iDualPor,ThNewIm,err)
        if(err.ne.0) goto 924
      end if
      call GenMat(NTab,NTabD,NMat,thr,ths,hSat,ParD,hTab,ConTab,CapTab,
     !            ConSat,TheTab,iModel,lScreen,nTabMod,ConSMax,xConv,
     !            tConv,err)
      if(err.eq.1) goto 903
      if(iHyst.ne.3) then
       call SetMat(NumNP,NTab,NTabD,NMat,hTab,ConTab,CapTab,hNew,MatNum,
     !            ParD,Con,Cap,ConSat,Ah,AK,ATh,hSat,hTemp,TheTab,ThOld,
     !            thr,ths,lWTDep,TempN,IterW,ConO,Kappa,AThS,ThRR,ConR,
     !            AKS,AhW,AThW,AKW,iModel,lTable,lVapor,ThVOld,ConLT,
     !            ConVT,ConVh,xConv,tConv,nTabMod,hCritA,lDensity,Conc,
     !            NSD,iEnhanc)
      else ! Bob Lenhard hysteresis
        call Hyst(NumNP,NMatD,ParD,ParW,MatNum,Kappa,hNew,hOld,ThOld,
     !            Con,Cap,IKappa,1)
      end if
      call InitDualPor(NumNP,NMat,Matnum,ParD,ThOld,iDualPor,ThNewIm,
     !                 ThOldIm,SinkIm,hNew,STrans,lInitW)
      call TmIn(tInit,tMax,tAtm,tOld,dt,dtMax,dMul,dMul2,dtMin,TPrint,
     !          t,dtOpt,TopInF,BotInF,lScreen,ItMin,ItMax,MaxAL,hCritS,
     !          NPD,AtmBC,iVer,lPrintD,nPrStep,tPrintInt,lEnter,lDayVar,
     !          lSinPrec,lLAI,rExtinct,err)
      if(err.ne.0) goto (907,928) err
      tPrint1=tMax
      if(lPrintD) TPrint1=tInit+tPrintInt
      if(lMeteo) then
        call MeteoIn(rLat,rAlt,ShWRadA,ShWRadB,rLWRadA,rLWRadB,rLWRadA1,
     !               rLWRadB1,WindHeight,TempHeight,iCrop,iLAI,
     !               CropHeight,Albedo,xLAI,xRoot,iInterc,aInterc,
     !               nGrowth,rGrowth,rExtinct,iRadiation,lEnBal,lPrint,
     !               iSunSh,iRelHum,CloudF_Ac,CloudF_Bc,lHargr,
     !               lMetDaily,xConv,err)
        if(err.ne.0) goto (932) err
        if(lEnBal) iCrop=0
      end if
      tAtm2=tMax
      dtInit=dt

      if(lRoot) then
        call RootIn(tRMin,tRHarv,xRMin,xRMax,RGR,lScreen,iver,iRootIn,
     !              nGrowth,rGrowth,tRPeriod,err)
        if(err.ne.0) goto (908,903) err
	end if
      if(lRoot.or.iCrop.gt.1)
     !  call SetRG(NumNP,x,Beta,t,tRMin,tRHarv,xRMin,xRMax,RGR,xRoot,
     !             lRoot,iRootIn,nGrowth,rGrowth,tRPeriod)
      if(lTemp) then
        call TempIn(NMat,TPar,Ampl,tPeriod,kTopT,tTop,kBotT,tBot,TopInF,
     !              BotInF,iCampbell,iVer,SnowMF,lScreen,ierr)
        if(err.ne.0) goto (910,903) err
        if(lEnBal) kTopT=-1
      end if

      if(lChem) then
        call ChemIn(lUpW,lSTDep,NMat,NS,NSD,MaxItC,ChPar,TDep,kTopCh,
     !              cTop,kBotCh,cBot,epsi,tPulse,CumCh,cTolA,cTolR,
     !              lLinear,lEquil,lArtD,PeCr,lScreen,dSurf,cAtm,lTort,
     !              lMobIm,lBact,lFiltr,iMoistDep,WDep,NMatD,iModel,
     !              ParD,iVer,lDualNEq,lMassIni,lEqInit,iTort,err)
        if(err.ne.0) goto (911,903,931) err
        if(iDualPor.gt.0) lEquil=.false.
        if(iMoistDep.eq.2)
     !    call MoistDepIn(cDataPath,cFileName,NMat,NMatD,NS,NSD,DMoist,
     !                    iMoistDep)
        if(lMassIni)
     !    call MassInit(NumNP,NSD,NS,NMat,MatNum,TDep,TempN,ChPar,Conc,
     !                  ThOld,ThOldIm,ths,lLinear,lBact)
        if(.not.lEquil.and.(lEqInit.or.lMassIni))
     !    call NonEqInit(NumNP,NSD,NS,NMat,MatNum,TDep,TempO,ChPar,
     !                   Conc,Sorb,lLinear,lMobIm,iDualPor,lBact,Sorb2,
     !                   ThOld)
        if(lPrint) then 
          call OpenSoluteFiles(NS,cDataPath,iLengthPath,cFileName,err)
          if(err.eq.1) goto 902
        end if
      end if
      if(TopInF.or.BotInF.or.AtmBC) then
        call SetBC(tMax,tAtm1,rTop,rRoot,rBot,hCritA,hBot,hTop,GWL0L,
     !             TopInF,BotInF,cT,cBot,NS,tTop,tBot,Ampl,lTemp,lChem,
     !             KodTop,lVarBC,err,lMinStep,lMeteo,Prec,rSoil,lLAI,
     !             rExtinct,lCentrif,CosAlf,xConv,tConv,iModel,
     !             hNew(NumNP),iRootIn,xRoot,WLayer,lLinear,lActRSU,
     !             SPot)
        if(err.eq.1) goto 913
        if(lMeteo) then
          call Meteo(1,lMetDaily,lDayVar,t,dt,tInit,tMax,tAtm2,
     !               tAtmN,tAtm2O,dtMax,rLat,rAlt,ShWRadA,ShWRadB,
     !               rLWRadA,rLWRadB,rLWRadA1,rLWRadB1,WindHeight,
     !               TempHeight,iCrop,iLAI,rRoot,xConv,tConv,rGrowth,
     !               nGrowth,iInterc,rInterc,aInterc,ExcesInt,lEnBal,
     !               rExtinct,lPrint,lHargr,iRadiation,iSunSh,iRelHum,
     !               iMetHour,CloudF_Ac,CloudF_Bc,Prec,Precc,rSoil,
     !               EvapP,TransP,Rns,Rnl,RadTerm,AeroTerm,Rst,ETcomb,
     !               Rad,RadN,RadO,Wind,WindN,WindO,Albedo,AlbedoN,
     !               xLAI,xLAIN,xRoot,xRootN,CropHeight,CropHeightN,
     !               Ampl,tTop,TMaxAN,TMinAN,TMax1,TMaxN,TMaxO,TMin1,
     !               TMinN,TMinO,TempA,TMaxA,TMaxAO,TMinA,TMinAO,
     !               SunHours,SunHoursN,SunHoursO,RHMean,RHMeanN,
     !               RHMeanO,RHMax,RHMaxN,RHMaxO,RHMin,RHMinN,RHMinO,
     !               RH_A,EaMean,EaMeanN,rTop,err)
          if(err.ne.0) goto (932,913,933) err
        end if
        if(lChem.and.SnowLayer.le.0)
     !    call SetChemBC(Prec,rSoil,NS,cTop,cT,WLayer,hNew(NumNP),
     !                   KodTop,kTopCh)
        if(lSnow) then
          dtSnow=sngl(tAtm1-tInit) 
          call Snow(Prec,dtSnow,tTop,SnowMF,SnowLayer,rSoil,xConv,
     !              lMinStep,cTop,cT,NS)
        end if
        tAtm=dmin1(tAtm1,tAtm2)
        if(lDayVar) then
          rRootD=rRoot
          rSoilD=rSoil
          call DailyVar(tConv,t,rRoot,rRootD)
          call DailyVar(tConv,t,rSoil,rSoilD)
        end if
        if(lSinPrec) then
          PrecD=Prec
          tAtmOld=sngl(tInit)
          call SinPrec(t,tAtmOld,tAtm1,Prec,PrecD)
        end if
        if(KodTop.eq.-4) rTop=abs(rSoil)-abs(Prec)
      end if
      if(SinkF) then
        call SinkIn(NMat,lChem,lMoSink,lSolRed,lSolAdd,P0,POptm,P2H,
     !              P2L,P3,r2H,r2L,aOsm,c50,P3c,NS,lMsSink,cRootMax,
     !              iVer,OmegaC,lActRSU,OmegaS,SPot1,rKM,cMin,lOmegaW,
     !              lScreen,err)
        if(err.ne.0) goto (909,903) err
        if(.not.TopInF.and..not.BotInF.and..not.AtmBC) SPot=SPot1
        call SetSnk(NumNP,NMat,MatNum,x,hRoot,vRoot,Sink,rRoot,hNew,
     !              lMoSink,lSolRed,lSolAdd,P0,POptm,P2H,P2L,P3,r2H,r2L,
     !              aOsm,c50,P3c,Beta,lChem,NS,NSD,Conc,cRoot,lMsSink,
     !              ThOld,ParD,dt,OmegaC,iModel,Con,lOmegaW,OmegaW,rBot)
      end if
      close(30)
      close(32)
      close(50)

      if(lPrint) then
        call Profil (NumNP,NMat,x,MatNum,xSurf,Beta,Ah,AK,ATh,thr,ths,
     !               ConSat,hSat,lScreen,err)
        if(err.eq.1) goto 923
        call NodOut (NumNP,NMat,hNew,ThOld,Con,x,xSurf,CosAlf,tInit,
     !               MatNum,Cap,AK,Sink,ConSat,NS,NSD,Conc,TempO,Sorb,
     !               Kappa,lBact,Sorb2,lVapor,lWTDep,ConLT,ConVT,ConVh,
     !               ThOld(NumNP),dt,iDualPor,ThNewIm,SinkIm,STrans,
     !               lDensity,lCentrif,Radius,lVaporOut,lDualNEq,err)    
        if(err.eq.1) goto 920
      end if
      call SubReg (NumNP,NMat,NLay,hNew,ThOld,ThOld,x,MatNum,LayNum,
     !             t-dt,dt,CosAlf,Con,lChem,Conc,ChPar,0,ths,wCumT,
     !             wCumA,cCumT,cCumA,wVolI,cVolI,WatIn,SolIn,lWat,lTemp,
     !             TempN,TPar,TDep,NS,NSD,Sorb,lLinear,lEquil,lMobIm,
     !             err,SubVol,Area,lPrint,lBact,Sorb2,lVapor,ThVOld,
     !             ThVNew,lWTDep,ConLT,ConVh,ConVT,iDualPor,ThNewIm,
     !             ThOldIm,lDensity,lCentrif,Radius,lDualNEq,cPrevO)
      if(err.eq.1) goto 921
      if(lChem.or.lTemp.or.lFlux)
     !  call Veloc(NumNP,hOld,Con,x,CosAlf,vOld,ThOld,ThOld,Sink,dt,
     !             lVapor,lWTDep,ConLT,ConVT,ConVh,TempO,vVOld,ThVOld,
     !             ThVOld,lDensity,Conc,NSD,lCentrif,Radius)
      do 11 i=1,NumNP
        vNew(i) =vOld(i)
        ThNew(i)=ThOld(i)
        if(lVapor) ThVNew(i)=ThVOld(i)
        if(lVapor) vVNew(i) =vVOld(i)
11    continue
      if(lBact.and..not.lWat) 
     !  call Exclusion(NumNP,NMat,NSD,ParD,ChPar,ThNew,vNew,ThOld,vOld)

      if(lScreen) write(*,*) 'beginning of numerical solution'
      ! using keyword arguments

      Rtime1=RTime(iMonth,iDay,iHours,iMins,iSecs,i100th)

*     Time stepping calculations ---------------------------------------
12    continue

*     Loop between water flow and heat and vapor transport
      iTemp=0
13    continue

      if(lEnBal) then
        EpsiT=0.5
        TempS=EpsiT*TempN(NumNP)+(1.-EpsiT)*TempO(NumNP)
        call Evapor(t,TempS,TMaxA,TMinA,Rad,hNew(NumNP),TempHeight,
     !              WindHeight,Wind,RHMean,HeatFl,rTop,Prec,tPeriod,
     !              rLat,Albedo,SunHours,ThNew(NumNP),xConv,tConv,
     !              iRadiation,Rns,Rnl,Rn,SensFlux,Evap,xLat,Const,
     !              iSunSh,r_H,lMetDaily,Rst,TempA,RH_A,ShWRadA,ShWRadB,
     !              rLWRadA,rLWRadB,iMetHour)
      end if

*     Solve water movement ---------------------------------------------
      if(lWat) then
        ItCumO=ItCum
        call WatFlow(NumNP,NTab,NTabD,NMat,hTab,ConTab,CapTab,hNew,hOld,
     !             MatNum,ParD,ParW,Con,Cap,ConSat,Ah,AK,ATh,hSat,hTemp,
     !               KodTop,KodBot,rTop,rBot,CosAlf,t,dt,x,Sink,P,R,S,
     !               FreeD,SeepF,qGWLF,Aqh,Bqh,GWL0L,hTop,hBot,hCritA,
     !               hCritS,WLayer,IterW,ItCum,TopInf,kTOld,kBOld,TolTh,
     !               TolH,MaxIt,dtMin,tOld,dtOpt,ConvgF,TheTab,ThNew,
     !               ThOld,thr,ths,lWTDep,TempN,Kappa,KappaO,AThS,ThRR,
     !               ConO,ConR,AKS,AhW,AThW,AKW,iHyst,iModel,qDrain,
     !               zBotDr,BaseGW,rSpacing,iPosDr,rKhTop,rKhBot,rKvTop,
     !               rKvBot,Entres,WetPer,zInTF,GeoFac,lTable,lVapor,
     !               xConv,tConv,ConLT,ConVT,ConVh,TauW,ThEq,ThVNew,
     !               ThVOld,nTabMod,iDualPor,ThNewIm,ThOldIm,SinkIm,
     !               vTop,TempO,iTemp,WTransf,lDensity,Conc,NSD,iEnhanc,
     !               lCentrif,Radius,hSeep)
        if(.not.ConvgF) then
          iNonConv=iNonConv+1
        else
          iNonConv=0
        end if
        if(lStopConv.and..not.ConvgF.and.iNonConv.ge.10) then
          call CloseOutput(RTime1,NS,TopInF,BotInF,lChem,lScreen,lMeteo,
     !                     lPrint)
          if(lEnter) then
            write(*,*) 'Press Enter to continue'
            read(*,*)
          end if
          stop
        end if
        if(lMeteo.and.lMetDaily.and..not.lEnBal) ! Output for daily variated meteo information
     !    call DayMeteoOut(t,ETcomb,EvapP,TransP,Rns,Rnl,RadTerm,
     !                     AeroTerm,Precc,rInterc,ExcesInt,TempA,RH_A,
     !                     Rst,lPrint)
        if(lEnBal) then
          if(ItCum.gt.ItCumO+MaxIt.and..not.lMetDaily)
     !    call MeteoInt(3,t,tAtm2O,tAtm2,Rad,RadO,RadN,TMaxA,TMaxAO,
     !                  TMaxAN,TMinA,TMinAO,TMinAN,Wind,WindO,WindN,
     !                  RHMean,RHMeanO,RHMeanN,SunHours,SunHoursO,
     !                  SunHoursN,lEnBal)
          TempS=EpsiT*TempN(NumNP)+(1.-EpsiT)*TempO(NumNP)
          M=MatNum(NumNP)
          rLamb=amax1(0.,TPar(4,M)+TPar(5,M)*ThNew(NumNP)+
     !                TPar(6,M)*sqrt(ThNew(NumNP)))
          dz=x(NumNP)-x(NumNP-1)
          call UpdateEnergy(t,vTop,rTop,HeatFl,TempS,Rns,Rnl,Rn,Evap,
     !                      xLat,SensFlux,xConv,tConv,Const,TLevel,
     !                      nPrStep,r_h,dz,rLamb,iTemp,lPrint,lMetDaily,
     !                      TempA,RH_A,Rst)
        end if
      else
        iterW=1
        ItCum=ItCum+1
      end if

*     To calculate the velocities --------------------------------------
      if(lWat.and.(lTemp.or.lChem.or.lFlux))
     !  call Veloc(NumNP,hNew,Con,x,CosAlf,vNew,ThNew,ThOld,Sink,dt,
     !             lVapor,lWTDep,ConLT,ConVT,ConVh,TempN,vVNew,ThVNew,
     !             ThVOld,lDensity,Conc,NSD,lCentrif,Radius)

*     Root zone calculations
      if(lRoot.or.iCrop.gt.1)
     !  call SetRG(NumNP,x,Beta,t,tRMin,tRHarv,xRMin,xRMax,RGR,xRoot,
     !             lRoot,iRootIn,nGrowth,rGrowth,tRPeriod)
      if(SinkF)
     !  call SetSnk(NumNP,NMat,MatNum,x,hRoot,vRoot,Sink,rRoot,hNew,
     !              lMoSink,lSolRed,lSolAdd,P0,POptm,P2H,P2L,P3,r2H,
     !              r2L,aOsm,c50,P3c,Beta,lChem,NS,NSD,Conc,cRoot,
     !              lMsSink,ThNew,ParD,dt,OmegaC,iModel,Con,lOmegaW,
     !              OmegaW,rBot)

*     Calculation of heat transport ------------------------------------
      if(lTemp) then
        call Temper(NumNP,NMat,x,dt,t,MatNum,TempO,TempN,TPar,Ampl,P,R,
     !              S,Q,vOld,vNew,ThOld,ThNew,Retard,Disp,Sink,tPeriod,
     !              kTopT,tTop,kBotT,tBot,lVapor,ThVOld,ThVNew,vVOld,
     !              vVNew,g0,lEnBal,HeatFl,xConv,tConv,dtMaxT,iCampbell,
     !              iTemp)
        if((lVapor.or.lEnBal).and.iTemp.lt.1) then
          iTemp=iTemp+1
          goto 13
        end if
      end if

*     Calculations of the solute transport -----------------------------
      if(lChem) then
        iKod=0
        if(lFlux) iKod=2
        call Solute(NumNP,NMat,NS,NSD,x,dt,t,tPulse,ChPar,MatNum,ThOld,
     !              ThNew,vOld,vNew,Disp,epsi,kTopCh,cTop,kBotCh,cBot,
     !              Conc,P,R,S,Q,g0,g1,Retard,cvTop,cvBot,cvCh0,cvCh1,
     !              lUpW,wc,Peclet,Courant,dtMaxC,TempO,TempN,cNew,
     !              cPrevO,cTemp,TDep,ths,cTolA,cTolR,IterC,MaxItC,
     !              hTemp,Sorb,SorbN,lLinear,lEquil,lArtD,PeCr,q0,q1,
     !              dSurf,cAtm,lTort,Sink,cRootMax,sSink,cvChR,lMobIm,
     !              cvChIm,TLevel,lBact,Sorb2,SorbN2,dtMin,dtOpt,lWat,
     !              lFiltr,iDualPor,ThOldIm,ThNewIm,SinkIm,STrans,
     !              iTort,xConv,tConv,lVapor,rBot,err,iMoistDep,NMatD,
     !              DMoist,WDep,iKod,Beta,lDualNEq,AtmBC,SinkF,lActRSU,
     !              OmegaS,OmegaW,SPot,rKM,cMin,lDensity)
        if(err.ne.0) goto 927
      end if

*     Output ------------------------------------------------------------
*     T-level information 
      if(abs(t-tMax).le.0.5*dtMin.or.t.gt.tMax) lEnd=.true.
      jPrint=0
      if((abs(TPrint(PLevel)-t).lt.0.001*dt.or.
     !    (lPrintD.and.(abs(TPrint1-t).lt.0.001*dt)).or. 
     !    (.not.ShortO.and.abs(float((TLevel+nPrStep-1)/nPrStep)-
     !    (TLevel+nPrStep-1)/float(nPrStep)).lt.0.0001))) jPrint=1

      call TLInf(NumNP,Con,x,CosAlf,t,dt,IterW,IterC,TLevel,rTop,rRoot,
     !           vRoot,hNew,hRoot,CumQ,ItCum,KodTop,KodBot,ConvgF,lWat,
     !           lChem,cRoot,NS,NSD,Conc,cvTop,cvBot,cvCh0,cvCh1,Peclet,
     !           Courant,wCumT,wCumA,cCumT,cCumA,CumCh,ThNew,ThOld,Sink,
     !           lScreen,err,cvChR,cvChIm,lPrint,lVapor,lWTDep,ConLT,
     !           ConVh,ConVT,TempN,rSoil,Prec,nPrStep,ThVOld,ThVNew,
     !           xConv,iDualPor,SinkIm,WTransf,lDensity,SnowLayer,
     !           lCentrif,Radius,ThNewIm,WLayer,hCritS,lEnd,lFluxOut,
     !           jPrint,dummy1,dummy2,lFlux,NObs,Node,vNew,cNew,CTop)
      if(err.ne.0) goto (919,918,926) err

      if(NObs.gt.0.and.lPrint.and.jPrint.eq.1) then
        call ObsNod(t,NumNP,NObs,NS,NSD,Node,Conc,hNew,ThNew,TempN,
     !              lChem,ThNewIm,vNew,vVNew,lFlux,err)
        if(err.eq.1) goto 922
      end if
      if(lPrintD.and.dabs(TPrint1-t).lt.0.001*dt)
     !                                         TPrint1=TPrint1+tPrintInt

*     P-level information ----------------------------------------------
      if(abs(TPrint(PLevel)-t).lt.0.001*dt) then
        if(lPrint) then
          call NodOut(NumNP,NMat,hNew,ThNew,Con,x,xSurf,CosAlf,
     !                TPrint(PLevel),MatNum,Cap,AK,Sink,ConSat,NS,NSD,
     !                Conc,TempN,Sorb,Kappa,lBact,Sorb2,lVapor,lWTDep,
     !                ConLT,ConVT,ConVh,ThOld(NumNP),dt,iDualPor,
     !                ThNewIm,SinkIm,STrans,lDensity,lCentrif,Radius,
     !                lVaporOut,lDualNEq,err)
          if(err.eq.1) goto 920
        end if
        call SubReg(NumNP,NMat,NLay,hNew,ThNew,ThOld,x,MatNum,LayNum,t,
     !              dt,CosAlf,Con,lChem,Conc,ChPar,PLevel,ths,wCumT,
     !              wCumA,cCumT,cCumA,wVolI,cVolI,WatIn,SolIn,lWat,
     !              lTemp,TempN,TPar,TDep,NS,NSD,Sorb,lLinear,lEquil,
     !              lMobIm,err,SubVol,Area,lPrint,lBact,Sorb2,lVapor,
     !              ThVOld,ThVNew,lWTDep,ConLT,ConVh,ConVT,iDualPor,
     !              ThNewIm,ThOldIm,lDensity,lCentrif,Radius,lDualNEq,
     !              cPrevO)
        if(err.eq.1) goto 921
        PLevel=PLevel+1
      end if

*     A-level information ----------------------------------------------
      if(dabs(t-tAtm).le.0.001*dt.and.(TopInF.or.BotInF.or.AtmBC)) then
        if(lPrint) then
          call ALInf(t,CumQ,hNew(NumNP),hRoot,hNew(1),ALevel,err)
          if(err.ne.0) goto (925) err
        end if
        if(dabs(t-tAtm1).le.0.001*dt) then
          tAtmOld=tAtm1
          call SetBC(tMax,tAtm1,rTop,rR,rBot,hCritA,hBot,hTop,GWL0L,
     !               TopInF,BotInF,cT,cBot,NS,tTop,tBot,Ampl,lTemp,
     !               lChem,KodTop,lVarBC,err,lMinStep,lMeteo,Prec,rS,
     !               lLAI,rExtinct,lCentrif,CosAlf,xConv,tConv,
     !               iModel,hNew(NumNP),iRootIn,xRoot,WLayer,lLinear,
     !               lActRSU,SPot)
          if(err.eq.1) goto 913
          if(.not.lMeteo) then
            rRoot=rR
            rSoil=rS
          end if
        end if
        if(lMeteo.and.dabs(t-tAtm2).le.0.001*dt) then
          call Meteo(2,lMetDaily,lDayVar,t,dt,tInit,tMax,tAtm2,
     !               tAtmN,tAtm2O,dtMax,rLat,rAlt,ShWRadA,ShWRadB,
     !               rLWRadA,rLWRadB,rLWRadA1,rLWRadB1,WindHeight,
     !               TempHeight,iCrop,iLAI,rRoot,xConv,tConv,rGrowth,
     !               nGrowth,iInterc,rInterc,aInterc,ExcesInt,lEnBal,
     !               rExtinct,lPrint,lHargr,iRadiation,iSunSh,iRelHum,
     !               iMetHour,CloudF_Ac,CloudF_Bc,Prec,Precc,rSoil,
     !               EvapP,TransP,Rns,Rnl,RadTerm,AeroTerm,Rst,ETcomb,
     !               Rad,RadN,RadO,Wind,WindN,WindO,Albedo,AlbedoN,
     !               xLAI,xLAIN,xRoot,xRootN,CropHeight,CropHeightN,
     !               Ampl,tTop,TMaxAN,TMinAN,TMax1,TMaxN,TMaxO,TMin1,
     !               TMinN,TMinO,TempA,TMaxA,TMaxAO,TMinA,TMinAO,
     !               SunHours,SunHoursN,SunHoursO,RHMean,RHMeanN,
     !               RHMeanO,RHMax,RHMaxN,RHMaxO,RHMin,RHMinN,RHMinO,
     !               RH_A,EaMean,EaMeanN,rTop,err)
          if(err.ne.0) goto (932,913,933) err
        end if
        tAtm=dmin1(tAtm1,tAtm2)
        if(lChem.and.SnowLayer.le.0)
     !    call SetChemBC(Prec,rSoil,NS,cTop,cT,WLayer,hNew(NumNP),
     !                   KodTop,kTopCh)
        if(lSnow) then
          dtSnow=sngl(tAtm1-tAtmOld)
          call Snow(Prec,dtSnow,tTop,SnowMF,SnowLayer,rSoil,xConv,
     !              lMinStep,cTop,cT,NS)
        end if
        if(lDayVar) then
          rRootD=rRoot
          rSoilD=rSoil
        end if
        if(lSinPrec) PrecD=Prec
        if(.not.lVarBC.and.KodTop.eq.-4) rTop=abs(rSoil)-abs(Prec)
        if(.not.lVarBC) rTop=abs(rSoil)-abs(Prec)
        ALevel=ALevel+1
      end if

      if(WLayer.and.hNew(NumNP).gt.0.) then ! mass balance in the surface layer
        hT=hNew(NumNP)
        do 15 jj=1,NS
          if((hT+dt*(Prec-rSoil)).gt.0.) 
     !      cTop(jj)=(hT*cTop(jj)+dt*Prec*cT(jj))/(hT+dt*(Prec-rSoil))
15      continue
      end if

*     Time governing ---------------------------------------------------
      if(abs(t-tMax).le.0.5*dtMin.or.t.gt.tMax) then
        call CloseOutput(RTime1,NS,TopInF,BotInF,lChem,lScreen,lMeteo,
     !                   lPrint)
        if(lEnter) then
          write(*,*) 'Press Enter to continue'
          read(*,*)
        end if
        stop
      else
        tOld=t
        dtOld=dt
        kTOld=KodTop
        kBOld=KodBot
        if(.not.lWat) IterW=1
        Iter=max0(IterW,IterC)
        dtMaxA=min(dtMaxC,dtMaxT)
        if(lSinPrec.and.PrecD.gt.0) 
     !    dtMaxA=min(dtMaxA,sngl(tAtm1-tAtmOld)/20.)
        call TmCont(dt,dtMax,dtOpt,dMul,dMul2,dtMin,Iter,
     !              min(TPrint(PLevel),TPrint1),tAtm,t,tMax,dtMaxA,
     !              ItMin,ItMax,lMinStep,dtInit)
        t=t+dt
        if(lDayVar) then
          call DailyVar(tConv,t,rRoot,rRootD)
          call DailyVar(tConv,t,rSoil,rSoilD)
          rTop=abs(rSoil)-abs(Prec)
        end if
        if(lSinPrec) then
          call SinPrec(t,tAtmOld,tAtm1,Prec,PrecD)
          rTop=abs(rSoil)-abs(Prec)
        end if
        if(lMeteo) then
          call Meteo(3,lMetDaily,lDayVar,t,dt,tInit,tMax,tAtm2,
     !               tAtmN,tAtm2O,dtMax,rLat,rAlt,ShWRadA,ShWRadB,
     !               rLWRadA,rLWRadB,rLWRadA1,rLWRadB1,WindHeight,
     !               TempHeight,iCrop,iLAI,rRoot,xConv,tConv,rGrowth,
     !               nGrowth,iInterc,rInterc,aInterc,ExcesInt,lEnBal,
     !               rExtinct,lPrint,lHargr,iRadiation,iSunSh,iRelHum,
     !               iMetHour,CloudF_Ac,CloudF_Bc,Prec,Precc,rSoil,
     !               EvapP,TransP,Rns,Rnl,RadTerm,AeroTerm,Rst,ETcomb,
     !               Rad,RadN,RadO,Wind,WindN,WindO,Albedo,AlbedoN,
     !               xLAI,xLAIN,xRoot,xRootN,CropHeight,CropHeightN,
     !               Ampl,tTop,TMaxAN,TMinAN,TMax1,TMaxN,TMaxO,TMin1,
     !               TMinN,TMinO,TempA,TMaxA,TMaxAO,TMinA,TMinAO,
     !               SunHours,SunHoursN,SunHoursO,RHMean,RHMeanN,
     !               RHMeanO,RHMax,RHMaxN,RHMaxO,RHMin,RHMinN,RHMinO,
     !               RH_A,EaMean,EaMeanN,rTop,err)
          if(err.ne.0) goto (932,913,933) err
        end if
        TLevel=TLevel+1
        if(TLevel.gt.1000000) TLevel=2
      end if

*     New updated values
      call Update(NumNP,lWat,lChem,lTemp,lVapor,iDualPor,lExtrap,dt,
     !            dtOld,hTemp,hNew,hOld,ThOld,ThNew,vOld,vNew,ThVOld,
     !            ThVNew,vVOld,vVNew,ThOldIm,ThNewIm,TempO,TempN,rTop,
     !            xConv,ConSMax,KodTop,KodBot)

      goto 12

* --- End of time loop -------------------------------------------------

*     Error messages
901   ierr=1
      goto 1000
902   ierr=2
      goto 1000
903   ierr=3
      goto 1000
904   ierr=4
      goto 1000
905   ierr=5
      goto 1000
906   ierr=6
      goto 1000
907   ierr=7
      goto 1000
908   ierr=8
      goto 1000
909   ierr=9
      goto 1000
910   ierr=10
      goto 1000
911   ierr=11
      goto 1000
912   ierr=12
      goto 1000
913   ierr=13
      goto 1000
914   ierr=14
      goto 1000
915   ierr=15
      goto 1000
916   ierr=16
      goto 1000
917   ierr=17
      goto 1000
918   ierr=18
      goto 1000
919   ierr=19
      goto 1000
920   ierr=20
      goto 1000
921   ierr=21
      goto 1000
922   ierr=22
      goto 1000
923   ierr=23
      goto 1000
924   ierr=24
      goto 1000
925   ierr=25
      goto 1000
926   ierr=26
      goto 1000
927   ierr=27
      goto 1000
928   ierr=28
      goto 1000
929   ierr=29
      goto 1000
930   ierr=30
      goto 1000
931   ierr=31
      goto 1000
932   ierr=32
      goto 1000
933   ierr=33
      goto 1000

1000  call ErrorOut(ierr,cFileName,cDataPath,iLengthPath,lScreen)
      if(lEnter) then
        write(*,*) 'Press Enter to continue'
        read(*,*)
      end if
      stop

101   format(a)
      end

************************************************************************

      subroutine ErrorOut(ierr,cFileName,cDataPath,iLengthPath,lScreen)

      character*260 cErr(33),cFileName,cDataPath,cFileNameErr
      logical lScreen

      cErr( 1)='Open file error in file :'
      cErr( 2)='File already exists or hard disk is full ! Open file err
     !or in output file : '
      cErr( 3)='Error when writing to an output file !'
      cErr( 4)='Error when reading from an input file Level_01.dir data
     !pathway !'
      cErr( 5)='Error when reading from an input file Selector.in Basic
     !Informations !'
      cErr( 6)='Error when reading from an input file Selector.in Water
     !Flow Informations !'
      cErr( 7)='Error when reading from an input file Selector.in Time I
     !nformations !'
      cErr( 8)='Error when reading from an input file Selector.in Root G
     !rowth Informations !'
      cErr( 9)='Error when reading from an input file Selector.in Sink I
     !nformations !'
      cErr(10)='Error when reading from an input file Selector.in Heat T
     !ransport Informations !'
      cErr(11)='Error when reading from an input file Selector.in Solute
     ! Transport Informations !'
      cErr(12)='Error when reading from an input file Profile.dat !'
      cErr(13)='Error when reading from an input file Atmosph.in !'
      cErr(14)='Dimension in NumNPD is exceeded !'
      cErr(15)='Dimension in NObsD is exceeded !'
      cErr(16)='Dimension in NMatD or NLay is exceeded !'
      cErr(17)='Error when writing into an output file I_CHECK.OUT !'
      cErr(18)='Error when writing into an output file RUN_INF.OUT !'
      cErr(19)='Error when writing into an output file T_LEVEL.OUT !'
      cErr(20)='Error when writing into an output file NOD_INF.OUT !'
      cErr(21)='Error when writing into an output file BALANCE.OUT !'
      cErr(22)='Error when writing into an output file OBS_NODE.OUT !'
      cErr(23)='Error when writing into an output file PROFILE.OUT !'
      cErr(24)='Initial water content condition is lower than Qr !'
      cErr(25)='Error when writing into an output file A_LEVEL.OUT !'
      cErr(26)='Error when writing into an output file SOLUTE.OUT !'
      cErr(27)='Does not converge in the solute transport module !'
      cErr(28)='Number of Print-Times is exceeded !'
      cErr(29)='Dimension in NSD is exceeded !'
      cErr(30)='The path to the project is too long !!!'
      cErr(31)='Bulk density can not be equal to zero !'
      cErr(32)='Error when reading from an input file Meteo.in !'
      cErr(33)='Crop Height must be smaller than the height of wind and 
     !temperature measurements!'

      cFileNameErr = cDataPath(1:iLengthPath)//'/Error.msg'
      open(99,file=cFileNameErr,status='unknown',err=901)
      if(ierr.le.2) then
        if(lScreen) write( *,*) cErr(ierr),cFileName
        write(99,*) cErr(ierr),cFileName
      else
        if(lScreen) write( *,*) cErr(ierr)
        write(99,*) cErr(ierr)
      end if
      close(99)
      return

901   write(*,*) 'Folder with input data of the specified project does n
     !ot exist or pathway is too long or corrupted'
      write(*,*) cFileName
      return
      end

************************************************************************

      subroutine CloseOutput(RTime1,NS,TopInF,BotInF,lChem,lScreen,
     !                       lMeteo,lPrint)

*      use MSFLIB
      logical TopInF,BotInF,lChem,lScreen,lMeteo,lPrint
      integer*2 iYear,iMonth,iDay,iHours,iMins,iSecs,i100th
      double precision RTime1,RTime2,RTime
*---- START gettim and getdat gnu alternative function ----
      integer dtv(8)
      character(8)  :: date
      character(10) :: time
      character(5)  :: zone
      integer hundred
      
      call date_and_time(date,time,zone,dtv)
      hundred=100
      iYear=dtv(1)
      iMonth=dtv(2)
      iDay=dtv(3)
      iHours=dtv(5)
      iMins=dtv(6)
      iSecs=dtv(7)
      i100th=mod(dtv(8),hundred)
*      call getdat(iYear,iMonth,iDay)
*      call gettim(iHours,iMins,iSecs,i100th)
*---- END gettim and getdat gnu alternative function ----

      Rtime2=RTime(iMonth,iDay,iHours,iMins,iSecs,i100th)
      if(lScreen) write( *,*) 'Real time [sec]',Rtime2-RTime1
      if(lPrint) then
        write(70,'(''end'')')
        write(76,*)
        write(76,*) 'Calculation time [sec]',Rtime2-RTime1
        write(71,'(''end'')')
        if(TopInF.or.BotInF) then
          write(72,'(''end'')')
          if(lMeteo) write(43,'(''end'')')
        end if
        write(77,'(''end'')')
        if(lChem) then
          do 13 jj=1,NS
            write(80+jj,'(''end'')')
13        continue
        end if
      end if
      return
      end

* ||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||

*-----------------------------------------------------------------------
*     Signal handler routines 
*-----------------------------------------------------------------------
*      integer(4) function h_sig (signum)
*        use MSFLIB
*        !MS$ATTRIBUTES C :: h_sig
*        integer(2) signum             
*        select case (signum)
*          case (SIG$ABORT)
*            write(*,*) 'Abnormal termination !'
*          case (SIG$INT) 
*            write(*,*) 'Program terminated !'
*          case (SIG$FPE)
*            write(*,*) 'Floating point error !'
*          case (SIG$BREAK) 
*            write(*,*) 'Program terminated !'
*          case default
*            write(*,*) 'Unknown signal handler !'
*        end select
*        h_sig = 1
*        call CloseFiles
*        stop
*      end
*
!      function hand_fpe (signum, excnum)
!        !MS$ATTRIBUTES C :: hand_fpe
!        use MSFLIB
!        integer(2)  signum, excnum
!        select case(excnum)
!          case(FPE$INVALID )
!            write(*,*) 'Floating Point Error - Invalid number !'
!          case( FPE$DENORMAL )
!            write(*,*) 'Floating Point Error - Denormalized number !'
!          case( FPE$ZERODIVIDE )
!            write(*,*) 'Floating Point Error - Zero divide !'
!          case( FPE$OVERFLOW )
!            write(*,*) 'Floating Point Error - Overflow !'
!          case( FPE$UNDERFLOW )
!            write(*,*) 'Floating Point Error - Underflow !'
!          case( FPE$INEXACT )
!            write(*,*) 'Floating Point Error - Inexact precision !'
!          case default
!            write(*,*) 'Floating Point Error - Non-IEEE type !'
!        end select
!        hand_fpe = 1
!        call CloseFiles
!      end
*
************************************************************************

      subroutine CloseFiles
      logical lOpen

      inquire(unit=43,opened=lOpen)
      if(lOpen) then
        write(43,'(''end'')')
        close(43)
      end if
      inquire(unit=44,opened=lOpen)
      if(lOpen) then
        write(44,'(''end'')')
        close(44)
      end if
      inquire(unit=70,opened=lOpen)
      if(lOpen) then
        write(70,'(''end'')')
        close(70)
      end if
      inquire(unit=71,opened=lOpen)
       if(lOpen) then
        write(71,'(''end'')')
        close(71)
      end if
      inquire(unit=72,opened=lOpen)
      if(lOpen) then
        write(72,'(''end'')')
        close(72)
      end if
      inquire(unit=77,opened=lOpen)
      if(lOpen) then
        write(77,'(''end'')')
        close(77)
      end if
      inquire(unit=75,opened=lOpen)
      if(lOpen) close(75)
      inquire(unit=76,opened=lOpen)
      if(lOpen) close(76)
      inquire(unit=78,opened=lOpen)
      if(lOpen) close(78)
      inquire(unit=31,opened=lOpen)
      if(lOpen) close(31)
      inquire(unit=81,opened=lOpen)
      if(lOpen) then
        write(81,'(''end'')')
        close(81)
      end if
      inquire(unit=82,opened=lOpen)
      if(lOpen) then
        write(82,'(''end'')')
        close(82)
      end if
      inquire(unit=83,opened=lOpen)
      if(lOpen) then
        write(83,'(''end'')')
        close(83)
      end if
      inquire(unit=41,opened=lOpen)
      if(lOpen) close(41)
      inquire(unit=42,opened=lOpen)
      if(lOpen) close(42)

      return
      end

************************************************************************
