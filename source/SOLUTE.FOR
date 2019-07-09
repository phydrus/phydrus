* Source file SOLUTE.FOR |||||||||||||||||||||||||||||||||||||||||||||||

*     To assemble and solve the solute transport equation
*     Mass-lumping finite elements

      subroutine Solute(N,NMat,NS,NSD,x,dt,t,tPulse,ChPar,MatNum,thO,
     !                  thN,vO,vN,Disp,epsi,kTopCh,cTop,kBotCh,cBot,
     !                  Conc,B,D,E,F,g0,g1,Retard,cvTop,cvBot,cvCh0,
     !                  cvCh1,lUpW,wc,Peclet,Courant,dtMaxC,TempO,TempN,
     !                  cNew,cPrevO,cTemp,TDep,thSat,cTolA,cTolR,IterC,
     !                  MaxItC,vCorr,Sorb,SorbN,lLinear,lEquil,lArtD,
     !                  PeCr,q0,q1,dSurf,cAtm,lTort,Sink,cRootMax,sSink,
     !                  cvChR,lMobIm,cvChIm,TLevel,lBact,Sorb2,SorbN2,
     !                  dtMin,dtOpt,lWat,lFiltr,iDualPor,ThOIm,ThNIm,
     !                  SinkIm,STrans,iTort,xConv,tConv,lVapor,rBot,
     !                  ierr,iMoistDep,NMatD,DMoist,WDep,iConcType,Beta,
     !                  lDualNEq,AtmBC,SinkF,lActRSU,OmegaS,OmegaW,SPot,
     !                  rKM,cMin,lDensity)

      logical lUpW,lConv,lLinear(NSD),lEquil,lArtD,lTort,lMobIm(NMat),
     !        lBact,lWat,lFiltr,lVapor,SinkF,lNEquil,lDualNEq,AtmBC,
     !        lActRSU,lDensity
      double precision B,D,E,F,t
      integer TLevel
      dimension x(N),ChPar(NSD*16+4,NMat),MatNum(N),thO(N),thN(N),vO(N),
     !          vN(N),Disp(N),cTop(NS),cBot(NS),Conc(NSD,N),B(N),D(N),
     !          E(N),F(N),g0(N),g1(N),Retard(N),cvTop(NS),cvBot(NS),
     !          cvCh0(NS),cvCh1(NS),wc(N),TempO(N),TempN(N),cNew(N),
     !          cPrevO(N),cTemp(N),TDep(NSD*16+4),thSat(NMat),sSink(N),
     !          vCorr(N),Sorb(NSD,N),SorbN(N),q0(N),q1(N),Sink(N),
     !          cvChR(NS),cRootMax(NS),cvChIm(NS),Sorb2(NSD,N),SorbN2(N)
     !         ,ThOIm(N),ThNIm(N),SinkIm(N),STrans(N),Beta(N),
     !          DMoist(NMatD,NSD,13,6),WDep(2+NMatD,NSD*9)

      alf=1.-epsi
      IterC=1.
      NLevel=2
      Peclet=0.
      Courant=0.
      dtMaxC=1.e+30
      rMin=1.e-30
      Tr=293.15
      R=8.314
*     Sequential first order decay goes into equilibrium phase (lNEquil=.false.) or nonequilbrium phase (lNEquil=.true.) 
      lNEquil=.false.

10    continue

*     Loop on species in the chain

      do 18 jS=1,NS
        Iter=0
        jjj=(jS-1)*16
        cvTop(jS)=0.
        cvBot(jS)=0.
        cvCh0(jS)=0.
        cvCh1(jS)=0.
        cvChR(jS)=0.
        cvChIm(jS)=0.
        if(t-tPulse.gt.dtMin.and..not.AtmBC) then
          cTop(jS)=0.
          cBot(jS)=0.
        end if
        if(kBotCh.lt.0) then
          if(vO(1).ge.0.) cvBot(jS)=alf*cBot(jS)  *vO(1)
          if(vO(1).lt.0.) cvBot(jS)=alf*Conc(jS,1)*vO(1)
          if(lVapor.and.rBot.eq.0.) cvBot(jS)=0.
        else if(kBotCh.eq.0) then
          cvBot(jS)=alf*Conc(jS,1)*vO(1)
        end if
        if(kTopCh.lt.0..and.TLevel.ne.1) then 
          if(vO(N).lt.0.) cvTop(jS)=alf*cTop(jS)*vO(N)
        end if
        if(kTopCh.eq.-2) then
          M=MatNum(N)
          Tr=293.15
          R=8.314
          TT=(TempO(N)+273.15-Tr)/R/(TempO(N)+273.15)/Tr
          Dg   =ChPar(jjj+ 6,M)*exp(TDep(jjj+ 6)*TT)
          Henry=ChPar(jjj+10,M)*exp(TDep(jjj+10)*TT)
          dSurfT=dSurf*exp(TDep(jjj+9)*TT)
          cvTop(jS)=cvTop(jS)+alf*Dg/dSurfT*Henry*Conc(jS,N)-
     !              Dg/dSurfT*cAtm
        end if
        if(.not.lLinear(jS)) then
          do 11 i=1,N
            cNew(i)=Conc(jS,i)
            if(.not.lEquil)       SorbN(i) =Sorb(jS,i)
            if(lBact.or.lDualNEq) SorbN2(i)=Sorb2(jS,i)
11        continue
        end if

*       Root Solute Uptake
        if(SinkF) 
     !    call SetSSnk(jS,NS,N,t,x,Beta,Sink,sSink,NSD,Conc,OmegaW,
     !                 cRootMax(jS),lActRSU,OmegaS,SPot,rKM,cMin)

*       Iterative loop for a nonlinear adsorption isotherm
12      Iter=Iter+1
        if(.not.lLinear(jS)) then
          do 13 i=1,N
            cTemp(i)=cNew(i)
13        continue
        end if

*       To construct the matrix equation
        do 15 Level=1,NLevel

*         Calculate the dispersion coefficients, retardation factors, source/
*         decay coefficients, Peclet and Courant numbers, upstream weighting
*         factors
          call Coeff(jS,Level,NLevel,N,NMat,NSD,x,Disp,vO,vN,thO,thN,
     !               thSat,ChPar,MatNum,TempN,TempO,TDep,g0,g1,Retard,
     !               Conc,cNew,cPrevO,dt,Pecl,Cour,dtMxC,lLinear,lEquil,
     !               lUpW,lArtD,Iter,wc,vCorr,Sorb,SorbN,epsi,PeCr,q0,
     !               q1,lTort,sSink,lMobIm,lBact,Sorb2,SorbN2,lFiltr,
     !               iDualPor,ThOIm,ThNIm,SinkIm,iTort,xConv,tConv,
     !               iMoistDep,NMatD,DMoist,WDep,lNEquil,lDualNEq)
          Peclet=amax1(Peclet,Pecl)
          Courant=amax1(Courant,Cour)
          dtMaxC=amin1(dtMaxC,dtMxC)

*         Set up the matrix equation
          call MatSet(jS,N,NS,NSD,Level,epsi,alf,dt,kBotCh,kTopCh,cBot,
     !                cTop,x,thO,thN,vO,vN,Conc,Disp,Retard,wc,g0,g1,B,
     !                D,E,F,E1,D1,F1,BN,DN,FN,NMat,ChPar,TempO,TempN,
     !                TDep,dSurfT,cAtm,MatNum,lMobIm,iDualPor,lVapor,
     !                rBot,lBact)
          do 14 i=1,N
            if(Level.eq.1) vO(i)=vO(i)+vCorr(i)
            if(Level.eq.2) vN(i)=vN(i)+vCorr(i)
14        continue

*         Calculate mass-transfer fluxes at the beginning of the time interval
          if(Level.eq.1.and.Iter.eq.1)
     !      call MassTran(jS,NS,NSD,N,MatNum,TempO,lMobIm,lEquil,ChPar,
     !                    TDep,Sorb,Conc,NMat,x,cvCh0,cvCh1,cvChR,
     !                    cvChIm,alf,q0,q1,sSink,lBact,ThO,Sorb2,lFiltr,
     !                    vO,iDualPor,SinkIm,xConv,tConv,lDualNEq,
     !                    STrans,lLinear)
15      continue

*       Solve matrix equation
        call BanSol(N,B,D,E,F)

*       Test for convergence for nonlinear problem
        lConv=.true.
        do 16 i=1,N
          if((NS.gt.1.and.Iter.eq.1).or.lDensity) cPrevO(i)=Conc(jS,i)
          if(lLinear(jS)) then
            Conc(jS,i)=amax1(sngl(F(i)),0.)
            if(Conc(jS,i).lt.1.e-30.and.Conc(jS,i).gt.0.) Conc(jS,i)=0.
          else
            cNew(i)=sngl(F(i))
            if(cNew(i).lt.1.0e-30) cNew(i)=0.
            if(abs(cNew(i)-cTemp(i)).gt.cTolA+cTolR*Conc(jS,i))
     !        lConv=.false.
          end if
16      continue
        if(.not.lLinear(jS)) then
          if(.not.lConv) then
            if(iter.lt.MaxItC) then
              goto 12
            else if(dt.gt.dtMin.and..not.lWat) then
c              ierr=1
              dtOld=dt
              dt=amax1(dt/3.,dtMin)
              dtOpt=dt
              t=t-dtOld+dt
              goto 10
            else
              ierr=1
            end if
          end if
          do 17 i=1,N
            Conc(jS,i)=cNew(i)
            if(.not.lEquil)       Sorb(jS,i) =SorbN(i)
            if(lBact.or.lDualNEq) Sorb2(jS,i)=SorbN2(i)
17        continue
        end if

*       Calculate sorbed concentration for linear noneq. adsorption or 
*       concentration in the imobile water.
        if(.not.lEquil.and.lLinear(jS)) 
     !    call SorbConc(jS,NSD,N,MatNum,TempN,lMobIm,ChPar,TDep,Sorb,
     !                  Conc,dt,NMat,lBact,thN,Sorb2,lFiltr,vN,iDualPor,
     !                  ThNIm,ThOIm,SinkIm,STrans,iMoistDep,NMatD,
     !                  DMoist,WDep,xConv,tConv,lDualNEq)

*       Calculate mass-transfer fluxes at the end of the time interval
        call MassTran(jS,NS,NSD,N,MatNum,TempN,lMobIm,lEquil,ChPar,TDep,
     !                Sorb,Conc,NMat,x,cvCh0,cvCh1,cvChR,cvChIm,epsi,
     !                q0,q1,sSink,lBact,ThN,Sorb2,lFiltr,vN,iDualPor,
     !                SinkIm,xConv,tConv,lDualNEq,STrans,lLinear)

*       Set up mass fluxes
        if(kTopCh.lt.0) then
          if(TLevel.ne.1) then
            if(vN(N).lt.0.) cvTop(jS)=cvTop(jS)+epsi*vN(N)*cTop(jS)
          else
            if(vN(N).lt.0.) cvTop(jS)=cvTop(jS)+     vN(N)*cTop(jS)
          end if
        else 
          cvTop(jS)=FN-BN*Conc(jS,N-1)-DN*Conc(jS,N)
        end if
        if(kTopCh.eq.-2) then
          M=MatNum(N)
          TT=(TempN(N)+273.15-Tr)/R/(TempN(N)+273.15)/Tr
          Dg   =ChPar(jjj+ 6,M)*exp(TDep(jjj+ 6)*TT)
          Henry=ChPar(jjj+10,M)*exp(TDep(jjj+10)*TT)
          cvTop(jS)=cvTop(jS)+epsi*Dg/dSurfT*Henry*Conc(jS,N)-
     !              Dg/dSurfT*cAtm
        end if
        if(kBotCh.lt.0) then
          if(vN(1).ge.0.) cvBot(jS)=cvBot(jS)+epsi*cBot(jS)  *vN(1)
          if(vN(1).lt.0.) cvBot(jS)=cvBot(jS)+epsi*Conc(jS,1)*vN(1)
          if(lVapor.and.rBot.eq.0.) cvBot(jS)=0.
        else if(kBotCh.eq.0) then 
          cvBot(jS)=cvBot(jS)+epsi*Conc(jS,1)*vN(1)
        else
          cvBot(jS)=D1*Conc(jS,1)+E1*Conc(jS,2)-F1
        end if
        IterC=max0(IterC,Iter)
        if(abs(cvTop(jS)).lt.rMin) cvTop(jS)=0.
        if(abs(cvBot(jS)).lt.rMin) cvBot(jS)=0.
18    continue

*     Calculate flux concentrations
      if(iConcType.eq.2)
     !  call FluxConc(N,NMat,NSD,x,vN,thN,thSat,ChPar,MatNum,TempN,TDep,
     !                Conc,cNew,lTort,lMobIm,iDualPor,ThNIm,1)
      return
      end

************************************************************************    

*     Calculate the dispersion coefficients, retardation factors, source/
*     decay coefficients, Peclet and Courant numbers, upstream weighting
*     factors

      subroutine Coeff(jS,Level,NLevel,NumNP,NMat,NSD,x,Disp,vO,vN,thO,
     !                 thN,thSat,ChPar,MatNum,TempN,TempO,TDep,g0,g1,
     !                 Retard,Conc,cNew,cPrevO,dt,Peclet,Courant,dtMaxC,
     !                 lLinear,lEquil,lUpW,lArtD,Iter,wc,vCorr,Sorb,
     !                 SorbN,epsi,PeCr,q0,q1,lTort,sSink,lMobIm,lBact,
     !                 Sorb2,SorbN2,lFiltr,iDualPor,ThOIm,ThNIm,SinkIm,
     !                 iTort,xConv,tConv,iMoistDep,NMatD,DMoist,WDep,
     !                 lNEquil,lDualNEq)

      logical lUpW,lLinear(NSD),lEquil,lArtD,lTort,lMobIm(NMat),lBact,
     !        lFiltr,lNEquil,
     !        lDualNEq
      dimension x(NumNP),Disp(NumNP),vO(NumNP),vN(NumNP),thO(NumNP),
     !          thN(NumNP),thSat(NMat),ChPar(NSD*16+4,NMat),g0(NumNP),
     !          g1(NumNP),MatNum(NumNP),Conc(NSD,NumNP),TempO(NumNP),
     !          TempN(NumNP),TDep(NSD*16+4),cNew(NumNP),cPrevO(NumNP),
     !          Retard(NumNP),wc(NumNP),Sorb(NSD,NumNP),SorbN(NumNP),
     !          vCorr(NumNP),q0(NumNP),q1(NumNP),
     !          sSink(NumNP),Sorb2(NSD,NumNP),SorbN2(NumNP),
     !          ThOIm(NumNP),ThNIm(NumNP),SinkIm(NumNP),
     !          DMoist(NMatD,NSD,13,6),WDep(2+NMatD,NSD*9)

*     Inicialization
      jjj=(jS-1)*16
      if(jS.gt.1) jj1=jjj-16
      Peclet=0.
      Courant=0.
      CourMax=1.
      dtMaxC=1.e+30
      Tr=293.15
      R=8.314

      do 11 i=NumNP,1,-1
        j=i+1
        k=i-1
        M=MatNum(i)
        if(Level.eq.NLevel) then
          ThW=ThN(i)
          ThWO=ThO(i)
          ThG=amax1(0.,thSat(M)-ThW)
          if(lMobIm(M).and.iDualPor.eq.0.or.lBact) then 
            ThImob=ChPar(4,M)
            ThImobO=ThImob
            ThW=amax1(ThW-ThImob,0.001)
          end if
          if(iDualPor.gt.0) then
            ThImob=ThNIm(i)
            ThImobO=ThOIm(i)
          end if
          v=vN(i)
          if(i.ne.NumNP) then
            vj=vN(j)
            Thj=ThN(j)
            if(lMobIm(M).and.iDualPor.eq.0.or.lBact) 
     !                                       Thj=amax1(Thj-ThImob,0.001)
          end if
          TT=(TempN(i)+273.15-Tr)/R/(TempN(i)+273.15)/Tr
          if(jS.gt.1) cPrev=Conc(jS-1,i)
        else
          ThW=ThO(i)
          ThG=amax1(0.,thSat(M)-ThW)
          if(lMobIm(M).and.iDualPor.eq.0.or.lBact) then
            ThImob=ChPar(4,M)
            ThW=amax1(ThW-ThImob,0.001)
          end if
          if(iDualPor.gt.0) then
            ThImob=ThOIm(i)
          end if
          v=vO(i)
          if(i.ne.NumNP) then
            vj=vO(j)
            Thj=ThO(j)
            if(lMobIm(M).and.iDualPor.eq.0.or.lBact) 
     !                                       Thj=amax1(Thj-ThImob,0.001)
          end if
          TT=(TempO(i)+273.15-Tr)/R/(TempO(i)+273.15)/Tr
          if(jS.gt.1) cPrev=cPrevO(i)
        end if

*       Temperature dependence
        f1=1.
        ro   =ChPar(1,     M)*exp(TDep(1)     *TT)
        Frac =ChPar(3,     M)*exp(TDep(3)     *TT)
        Dw   =ChPar(jjj+ 5,M)*exp(TDep(jjj+ 5)*TT)
        Dg   =ChPar(jjj+ 6,M)*exp(TDep(jjj+ 6)*TT)
        xKs  =ChPar(jjj+ 7,M)*exp(TDep(jjj+ 7)*TT)
        xNu  =ChPar(jjj+ 8,M)*exp(TDep(jjj+ 8)*TT)
        fExp =ChPar(jjj+ 9,M)
        Henry=ChPar(jjj+10,M)*exp(TDep(jjj+10)*TT)
        if(iMoistDep.gt.0)  
     !    f1=rMD(NMatD,NSD,M,jS,1,DMoist,1,WDep,ThW,iMoistDep)
        GamL  =ChPar(jjj+11,M)*exp(TDep(jjj+11)*TT)*f1
        if(iMoistDep.gt.0) 
     !    f1=rMD(NMatD,NSD,M,jS,10,DMoist,1,WDep,ThImob,iMoistDep)
        GamLi =ChPar(jjj+11,M)*exp(TDep(jjj+11)*TT)*f1    ! reaction in the immobile phase
        if(iMoistDep.gt.0) 
     !    f1=rMD(NMatD,NSD,M,jS,2,DMoist,2,WDep,ThW,iMoistDep)
        GamS  =ChPar(jjj+12,M)*exp(TDep(jjj+12)*TT)*f1
        if(iMoistDep.gt.0) 
     !    f1=rMD(NMatD,NSD,M,jS,11,DMoist,2,WDep,ThImob,iMoistDep)
        GamSi =ChPar(jjj+12,M)*exp(TDep(jjj+12)*TT)*f1
        if(iMoistDep.gt.0) 
     !    f1=rMD(NMatD,NSD,M,jS,3,DMoist,3,WDep,ThW,iMoistDep)
        GamG  =ChPar(jjj+13,M)*exp(TDep(jjj+13)*TT)*f1
        GamGi=GamG
        if(iMoistDep.gt.0) 
     !    f1=rMD(NMatD,NSD,M,jS,4,DMoist,4,WDep,ThW,iMoistDep)
        GamL1 =ChPar(jjj+14,M)*exp(TDep(jjj+14)*TT)*f1
        if(iMoistDep.gt.0) 
     !    f1=rMD(NMatD,NSD,M,jS,12,DMoist,4,WDep,ThImob,iMoistDep)
        GamL1i=ChPar(jjj+14,M)*exp(TDep(jjj+14)*TT)*f1
        if(iMoistDep.gt.0) 
     !    f1=rMD(NMatD,NSD,M,jS,5,DMoist,5,WDep,ThW,iMoistDep)
        GamS1 =ChPar(jjj+15,M)*exp(TDep(jjj+15)*TT)*f1
        if(iMoistDep.gt.0) 
     !    f1=rMD(NMatD,NSD,M,jS,13,DMoist,5,WDep,ThImob,iMoistDep)
        GamS1i=ChPar(jjj+15,M)*exp(TDep(jjj+15)*TT)*f1
        if(iMoistDep.gt.0) 
     !    f1=rMD(NMatD,NSD,M,jS,6,DMoist,6,WDep,ThW,iMoistDep)
        GamG1 =ChPar(jjj+16,M)*exp(TDep(jjj+16)*TT)*f1
        if(iMoistDep.gt.0) 
     !    f1=rMD(NMatD,NSD,M,jS,7,DMoist,7,WDep,ThW,iMoistDep)
        xMuL  =ChPar(jjj+17,M)*exp(TDep(jjj+17)*TT)*f1
        if(iMoistDep.gt.0) 
     !    f1=rMD(NMatD,NSD,M,jS,8,DMoist,8,WDep,ThW,iMoistDep)
        xMuS  =ChPar(jjj+18,M)*exp(TDep(jjj+18)*TT)*f1
        if(iMoistDep.gt.0) 
     !    f1=rMD(NMatD,NSD,M,jS,9,DMoist,9,WDep,ThW,iMoistDep)
        xMuG  =ChPar(jjj+19,M)*exp(TDep(jjj+19)*TT)*f1
        Omega =ChPar(jjj+20,M)*exp(TDep(jjj+20)*TT)
        f_em=1.
        if(lDualNEq) f_em  =ChPar(jjj+13,M)*exp(TDep(jjj+13)*TT)
        if(lDualNEq) OmegaS=ChPar(jjj+16,M)*exp(TDep(jjj+16)*TT)
        if(lBact) then
          Dg=0. 
          GamG =0.
          GamGi=0.
          GamL1=0.
          GamL1i=0.
          GamS1=0.
          GamS1i=0.
          GamG1=0.
          GamG1i=0.
          xMuL =0.
          xMuS =0.
          xMuG =0.
          Omega=0.
          SMax2 =ChPar(jjj+15,M)*exp(TDep(jjj+15)*TT)
          rKa2  =ChPar(jjj+16,M)*exp(TDep(jjj+16)*TT)
          rKd2  =ChPar(jjj+17,M)*exp(TDep(jjj+17)*TT)
          SMax1 =ChPar(jjj+18,M)*exp(TDep(jjj+18)*TT)
          rKa1  =ChPar(jjj+19,M)*exp(TDep(jjj+19)*TT)
          rKd1  =ChPar(jjj+20,M)*exp(TDep(jjj+20)*TT)
          iPsi1=0
          iPsi2=0
          if(.not.lFiltr) iPsi2=int(ChPar(jjj+13,M))
          if(.not.lFiltr) iPsi1=int(ChPar(jjj+14,M))
          if(iPsi1.eq.0.and.SMax1.gt.0.) iPsi1=1
          if(iPsi2.eq.0.and.SMax2.gt.0.) iPsi2=1
          if(iPsi1.ge.3.or.iPsi2.ge.3) 
     !      Dc=ChPar(jjj+6,M)*exp(TDep(jjj+6)*TT)
          if(iPsi1.eq.5.or.iPsi2.eq.5) aa=ChPar(jjj+15,M)
          if(Level.eq.NLevel) then
            ss1=SorbN(i)
            ss2=SorbN2(i)
          else
            ss1=Sorb(jS,i)
            ss2=Sorb2(jS,i)
          end if
          psi1=1.
          psi2=1.
          if(iPsi1.gt.0) call Blocking(iPsi1,SMax1,psi1,x(i),ss1,dc,aa)
          if(iPsi2.gt.0) call Blocking(iPsi2,SMax2,psi2,x(i),ss2,dc,aa)

*         recalculate ka1 and ka2 based on filtration theory
          if(lFiltr) then
            GamG =0.
            GamL1=0.
            Dc=ChPar(jjj+13,M)*exp(TDep(jjj+13)*TT)
            Dp=ChPar(jjj+14,M)*exp(TDep(jjj+14)*TT)
            Alfa1=rKa1
            Alfa2=rKa2
            call Deposit(rKa1,rKa2,Dc,Dp,Alfa1,Alfa2,ThW,v,TempN(i),
     !                   xConv,tConv)
          end if
        end if

        if(jS.gt.1) then
          xKsP  =ChPar(jj1+ 7,M)*exp(TDep(jj1+ 7)*TT)
          xNuP  =ChPar(jj1+ 8,M)*exp(TDep(jj1+ 8)*TT)
          fExpP =ChPar(jj1+ 9,M) !*exp(TDep(jj1+ 9)*TT)
          HenryP=ChPar(jj1+10,M)*exp(TDep(jj1+10)*TT)
          f1=1.
          if(iMoistDep.gt.0) 
     !      f1=rMD(NMatD,NSD,M,jS-1,4,DMoist,4,WDep,ThW,iMoistDep)
          GamL1P =ChPar(jj1+14,M)*exp(TDep(jj1+14)*TT)*f1
          if(iMoistDep.gt.0) 
     !      f1=rMD(NMatD,NSD,M,jS-1,12,DMoist,4,WDep,ThImob,iMoistDep)
          GamL1Pi=ChPar(jj1+14,M)*exp(TDep(jj1+14)*TT)*f1
          if(iMoistDep.gt.0) 
     !      f1=rMD(NMatD,NSD,M,jS-1,5,DMoist,5,WDep,ThW,iMoistDep)
          GamS1P =ChPar(jj1+15,M)*exp(TDep(jj1+15)*TT)*f1
          if(iMoistDep.gt.0) 
     !      f1=rMD(NMatD,NSD,M,jS-1,13,DMoist,5,WDep,ThImob,iMoistDep)
          GamS1Pi=ChPar(jj1+15,M)*exp(TDep(jj1+15)*TT)*f1
          if(iMoistDep.gt.0) 
     !      f1=rMD(NMatD,NSD,M,jS-1,6,DMoist,6,WDep,ThW,iMoistDep)
          GamG1P=ChPar(jj1+16,M)*exp(TDep(jj1+16)*TT)*f1
          if(lBact) then
            GamL1P=0.
            GamL1Pi=0.
            GamS1P=0.
            GamS1Pi=0.
            GamG1P=0.
          end if
        end if
        if(Level.eq.NLevel) then
          TTO=(TempO(i)+273.15-Tr)/R/(TempO(i)+273.15)/Tr
          xKsO  =ChPar(jjj+ 7,M)*exp(TDep(jjj+ 7)*TTO)
          xNuO  =ChPar(jjj+ 8,M)*exp(TDep(jjj+ 8)*TTO)
          fExpO =ChPar(jjj+ 9,M)
          HenryO=ChPar(jjj+10,M)*exp(TDep(jjj+10)*TTO)
          f1=1.
          if(iMoistDep.gt.0) 
     !      f1=rMD(NMatD,NSD,M,jS,1,DMoist,1,WDep,ThO(i),iMoistDep)
          GamLO =ChPar(jjj+11,M)*exp(TDep(jjj+11)*TTO)*f1
          if(iMoistDep.gt.0) 
     !      f1=rMD(NMatD,NSD,M,jS,10,DMoist,1,WDep,ThImobO,iMoistDep)
          GamLOi=ChPar(jjj+11,M)*exp(TDep(jjj+11)*TTO)*f1
          if(iMoistDep.gt.0) 
     !      f1=rMD(NMatD,NSD,M,jS,2,DMoist,2,WDep,ThO(i),iMoistDep)
          GamSO =ChPar(jjj+12,M)*exp(TDep(jjj+12)*TTO)*f1
          if(iMoistDep.gt.0) 
     !      f1=rMD(NMatD,NSD,M,jS,11,DMoist,2,WDep,ThImobO,iMoistDep)
          GamSOi=ChPar(jjj+12,M)*exp(TDep(jjj+12)*TTO)*f1
          if(iMoistDep.gt.0) 
     !      f1=rMD(NMatD,NSD,M,jS,4,DMoist,4,WDep,ThO(i),iMoistDep)
          GamL1O =ChPar(jjj+14,M)*exp(TDep(jjj+14)*TTO)*f1
          if(iMoistDep.gt.0) 
     !      f1=rMD(NMatD,NSD,M,jS,12,DMoist,4,WDep,ThImobO,iMoistDep)
          GamL1Oi=ChPar(jjj+14,M)*exp(TDep(jjj+14)*TTO)*f1
          if(iMoistDep.gt.0) 
     !      f1=rMD(NMatD,NSD,M,jS,5,DMoist,5,WDep,ThO(i),iMoistDep)
          GamS1O=ChPar(jjj+15,M)*exp(TDep(jjj+15)*TTO)*f1
          if(iMoistDep.gt.0) 
     !      f1=rMD(NMatD,NSD,M,jS,13,DMoist,5,WDep,ThImobO,iMoistDep)
          GamS1Oi=ChPar(jjj+15,M)*exp(TDep(jjj+15)*TTO)*f1
          if(iMoistDep.gt.0) 
     !      f1=rMD(NMatD,NSD,M,jS,7,DMoist,7,WDep,ThO(i),iMoistDep)
          xMuLO =ChPar(jjj+17,M)*exp(TDep(jjj+17)*TTO)*f1
          if(iMoistDep.gt.0) 
     !      f1=rMD(NMatD,NSD,M,jS,8,DMoist,8,WDep,ThO(i),iMoistDep)
          xMuSO =ChPar(jjj+18,M)*exp(TDep(jjj+18)*TTO)*f1
          OmegaO=ChPar(jjj+20,M)*exp(TDep(jjj+20)*TTO)
          dKs   =(xKs  -  xKsO)/dt
          dNu   =(xNu  -  xNuO)/dt
          ddExp =(fExp - fExpO)/dt
          dHenry=(Henry-HenryO)/dt
          if(i.ne.1)     TTi=(TempN(k)+273.15-Tr)/R/(TempN(k)+273.15)/Tr
          if(i.ne.NumNP) TTj=(TempN(j)+273.15-Tr)/R/(TempN(j)+273.15)/Tr
          if(lBact) then
            GamS1O=0.
            GamS1Oi=0.
            xMuLO =0.
            xMuSO =0.
            xMuGO =0.
            OmegaO=0.
            SMax2O =ChPar(jjj+15,M)*exp(TDep(jjj+15)*TTO)
            rKa2O  =ChPar(jjj+16,M)*exp(TDep(jjj+16)*TTO)
            rKd2O  =ChPar(jjj+17,M)*exp(TDep(jjj+17)*TTO)
            SMax1O =ChPar(jjj+18,M)*exp(TDep(jjj+18)*TTO)
            rKa1O  =ChPar(jjj+19,M)*exp(TDep(jjj+19)*TTO)
            rKd1O  =ChPar(jjj+20,M)*exp(TDep(jjj+20)*TTO)
            iPsi1=0
            iPsi2=0
            if(.not.lFiltr) iPsi2=int(ChPar(jjj+13,M))
            if(.not.lFiltr) iPsi1=int(ChPar(jjj+14,M))
            if(iPsi1.eq.0.and.SMax1O.gt.0.) iPsi1=1
            if(iPsi2.eq.0.and.SMax2O.gt.0.) iPsi2=1
            if(iPsi1.ge.3.or.iPsi2.ge.3) 
     !        Dc=ChPar(jjj+6,M)*exp(TDep(jjj+6)*TT)
            if(iPsi1.eq.5.or.iPsi2.eq.5) aa=ChPar(jjj+15,M)
            psi1O=1.
            psi2O=1.
            if(iPsi1.gt.0) 
     !        call Blocking(iPsi1,SMax1O,psi1O,x(i),Sorb(jS,i),dc,aa)
            if(iPsi2.gt.0) 
     !        call Blocking(iPsi2,SMax2O,psi2O,x(i),Sorb2(jS,i),dc,aa)
            if(lFiltr) then
              GamL1O=0.
              Dc=ChPar(jjj+13,M)*exp(TDep(jjj+13)*TTO)
              Dp=ChPar(jjj+14,M)*exp(TDep(jjj+14)*TTO)
              Alfa1=rKa1O
              Alfa2=rKa2O
              call Deposit(rKa1O,rKa2O,Dc,Dp,Alfa1,Alfa2,ThWO,vO(i),
     !                     TempO(i),xConv,tConv)
            end if
          end if
        else
          TTN=(TempN(i)+273.15-Tr)/R/(TempN(i)+273.15)/Tr
          xKsN  =ChPar(jjj+ 7,M)*exp(TDep(jjj+ 7)*TTN)
          xNuN  =ChPar(jjj+ 8,M)*exp(TDep(jjj+ 8)*TTN)
          fExpN =ChPar(jjj+ 9,M) !*exp(TDep(jjj+ 9)*TTN)
          HenryN=ChPar(jjj+10,M)*exp(TDep(jjj+10)*TTN)
          dKs   =(xKsN  -  xKs)/dt
          dNu   =(xNuN  -  xNu)/dt
          ddExp =(fExpN - fExp)/dt
          dHenry=(HenryN-Henry)/dt
          if(i.ne.1)     TTi=(TempO(k)+273.15-Tr)/R/(TempO(k)+273.15)/Tr
          if(i.ne.NumNP) TTj=(TempO(j)+273.15-Tr)/R/(TempO(j)+273.15)/Tr
        end if
        if(i.ne.1) Henryi=ChPar(jjj+10,MatNum(k))*exp(TDep(jjj+10)*TTi)
        if(i.ne.NumNP)
     !             Henryj=ChPar(jjj+10,MatNum(j))*exp(TDep(jjj+10)*TTj)

        dSConc=1.
        dConc=1.
        SConcP=1.
        SConc=1.
        SConcO=1.
        dRetard=0.

        SConcS=1.
        SConcOS=1.
        dSConcS=1.
        dConcS=1.
        SConcPS=1.
        dRetardS=0.

*       Effects of nonlinear adsorption
        if(.not.lLinear(jS)) then
          cc=Conc(jS,i)
          cMid=(Conc(jS,i)+cNew(i))/2.
          if(Level.eq.NLevel) cc=cNew(i)
          if(cc.gt.0.) then
            dSConc=fExp*cc**(fExp-1.)/(1.+xNu*cc**fExp)**2
            SConc =     cc**(fExp-1.)/(1.+xNu*cc**fExp)
          end if
          if(cMid.gt.0.) then
            dConc=fExp*cMid**(fExp-1.)/(1.+xNu*cMid**fExp)**2
            dRetard=cMid**fExp/(1.+xNu*cMid**fExp)*dKs-
     !           xKs*cMid**(2.*fExp)/(1.+xNu*cMid**fExp)**2*dNu+
     !           xKs*alog(cMid)*cMid**fExp/(1.+xNu*cMid**fExp)**2*ddExp
          end if
          if(Level.eq.NLevel.and..not.lEquil.and.Conc(jS,i).gt.0.)
     !      SConcO=Conc(jS,i)**(fExpO-1.)/(1.+xNuO*Conc(jS,i)**fExpO)
          if(lMobIm(M).or.iDualPor.gt.0) then     ! mobile-immobile model
            ss=Sorb(jS,i)
            sMid=(Sorb(jS,i)+SorbN(i))/2.
            if(Level.eq.NLevel) ss=SorbN(i)
            if(ss.gt.0.) then
              dSConcS=fExp*ss**(fExp-1.)/(1.+xNu*ss**fExp)**2
              SConcS =     ss**(fExp-1.)/(1.+xNu*ss**fExp)
            end if
            if(sMid.gt.0.) then
              dConcS=fExp*sMid**(fExp-1.)/(1.+xNu*sMid**fExp)**2
              dRetardS=sMid**fExp/(1.+xNu*sMid**fExp)*dKs-
     !            xKs*sMid**(2.*fExp)/(1.+xNu*sMid**fExp)**2*dNu+
     !            xKs*alog(sMid)*sMid**fExp/(1.+xNu*sMid**fExp)**2*ddExp
            end if
            if(Level.eq.NLevel.and..not.lEquil.and.Sorb(jS,i).gt.0.)
     !        SConcOS=Sorb(jS,i)**(fExpO-1.)/(1.+xNuO*Sorb(jS,i)**fExpO)
          end if
        else
          if(Conc(jS,i).gt.0.) dRetard=Conc(jS,i)*dKs
          if(lMobIm(M).or.iDualPor.gt.0) then     ! mobile-immobile model
            if(Sorb(jS,i).gt.0) dRetardS=Sorb(jS,i)*dKs
          end if
        end if
        if(jS.gt.1) then
          if(.not.lLinear(jS-1)) then
            if(cPrev.gt.0.)
     !        SConcP=cPrev**(fExpP-1.)/(1.+xNuP*cPrev**fExpP)
            if(Sorb(jS-1,i).gt.0.)
     !        SConcPS=Sorb(jS-1,i)**(fExpP-1.)/
     !                                     (1.+xNuP*Sorb(jS-1,i)**fExpP)
          end if
        end if

*       Calculate the retardation factors
        Retard(i)=(ro*Frac*f_em*xKs*dConc+ThG*Henry)/ThW+1.

*       Calculate the dispersion coefficients
        call Disper(i,M,NMat,NSD,NumNP,dt,lTort,lArtD,lUpW,lMobIm,
     !              iDualPor,Level,NLevel,ChPar,Retard,Disp,thSat,
     !              ThImob,ThW,ThG,v,Dw,Dg,Henry,ro,Frac,xKs,fExp,xNu,
     !              cMid,dSConc,PeCr,TauG,iTort,lBact)

*       Calculate the adsorbed concentration on kinetic sites or
*       the concentration in an imobile zone, before solving matrix equation
        if(.not.lEquil)
     !    call NEquil(i,jS,NSD,NumNP,NMat,M,Conc,Sorb,Sorb2,cNew,SorbN,
     !                SorbN2,SSorb,SSorb2,lMobIm,lLinear,lBact,Level,
     !                NLevel,dt,epsi,ro,xKs,xKsO,cc,SConc,SConcO,SConcS,
     !                SConcOS,dSConcS,xMuL,xMuLO,xMuS,xMuSO,dRetardS,
     !                GamLi,GamL1i,GamLOi,GamL1Oi,GamSi,GamS1i,GamSOi,
     !                GamS1Oi,Omega,OmegaO,rKa1,rKa1O,rKa2,rKa2O,rKd1,
     !                rKd1O,rKd2,rKd2O,ThW,ThWO,psi1,psi1O,psi2,psi2O,
     !                DMobI,Frac,ThImob,ThImobO,iDualPor,SinkIm,FlMacro,
     !                lNEquil,GamL1Pi,GamS1Pi,xKsP,SConcPS,lDualNEq,
     !                f_em,OmegaS)

*       Calculate zero-order coefficient g0
        g0(i)=xMuL*ThW+Frac*f_em*ro*xMuS+ThG*xMuG-sSink(i)
        q0(i)=xMuL*ThW+          ro*xMuS+ThG*xMuG
        if(.not.lEquil) then
          if((lMobIm(M).or.iDualPor.gt.0).and..not.lBact) then
            g0(i)=g0(i)+Omega*SSorb
            if(iDualPor.gt.0.and.SinkIm(i).le.0) g0(i)=g0(i)-FlMacro
            if(lDualNEq) g0(i)=g0(i)+OmegaS*ro*SSorb2
          else if(.not.lBact) then
            g0(i)=g0(i)+Omega*ro*SSorb
          else if(lBact) then
            g0(i)=g0(i)+rKd1*ro*SSorb+rKd2*ro*SSorb2
          end if
        end if
        if(jS.gt.1) then
          cG=cPrev*(GamL1P*ThW+ro*Frac*f_em*xKsP*GamS1P*SConcP+
     !              ThG*HenryP*GamG1P)
          cG1=cG
          if(.not.lEquil) then
            if((lMobIm(M).or.iDualPor.gt.0).and..not.lBact) then
              aa=Sorb(jS-1,i)*(ThImob*GamL1Pi+
     !                               (1.-Frac)*ro*GamS1Pi*xKsP*SConcPS)
              if(.not.lNEquil) cG=cG+aa
	        cG1=cG1+aa    
              if(lDualNEq) then
                aa=GamS1Pi*ro*Sorb2(jS-1,i)
                if(.not.lNEquil) cG=cG+aa
	          cG1=cG1+aa
              end if
            else if(.not.lBact) then
              aa=GamS1Pi*ro*Sorb(jS-1,i)
              if(.not.lNEquil) cG=cG+aa
	        cG1=cG1+aa
            else if(lBact) then
              write(*,*) 'Attachment/dettachment model is implemented 
     !only for one solute'
              write(*,*)'Press Enter to continue'
              read(*,*)
              stop
            end if
          end if
          g0(i)=g0(i)+cG
          q0(i)=q0(i)+cG1
        end if
        if(cMid.gt.0.) g0(i)=g0(i)-ro*Frac*f_em*dRetard
     
*       Calculate first-order coefficient g1
        g1(i)=-(GamL+GamL1)*ThW-(GamS+GamS1)*ro*Frac*f_em*xKs*SConc-
     !         (GamG+GamG1)*ThG*Henry
c        if(Level.eq.NLevel) g1(i)=g1(i)-ThG*dHenry-Henry*(ThWO-ThW)/dt
        if(.not.lEquil) then
          if((lMobIm(M).or.iDualPor.gt.0).and..not.lBact) then ! mobile-immobile model
            g1(i)=g1(i)-Omega
            if(iDualPor.gt.0.and.SinkIm(i).gt.0) g1(i)=g1(i)-SinkIm(i)
            if(Level.eq.NLevel.and.lLinear(jS))
     !         g1(i)=g1(i)+Omega*dt*Omega/DMobI
            if(lDualNEq) then
              g1(i)=g1(i)-OmegaS*ro*Frac*(1.-f_em)*SConc*xKs
              if(Level.eq.NLevel.and.lLinear(jS)) g1(i)=g1(i)+OmegaS*ro*
     !                  (dt*OmegaS*Frac*(1.-f_em)*xKs/
     !                  (2.+dt*(OmegaS+GamSi+GamS1i)))
            end if
          else if(.not.lBact) then                             ! two-site sorption model
            g1(i)=g1(i)-Omega*ro*(1.-Frac)*SConc*xKs
            if(Level.eq.NLevel.and.lLinear(jS)) g1(i)=g1(i)+Omega*ro*
     !             (dt*Omega*(1.-Frac)*xKs/(2.+dt*(Omega+GamSi+GamS1i)))
          else if(lBact) then                                  ! filtration model
            g1(i)=g1(i)-ThW*(rKa1*psi1+rKa2*psi2)
            if(Level.eq.NLevel.and.lLinear(jS)) g1(i)=g1(i)+dt*ThW*
     !               (rKd1*rKa1/(2.+dt*(rKd1+GamSi))+
     !                rKd2*rKa2/(2.+dt*(rKd2+GamSi)))
          end if
        end if
        q1(i)=(-(GamL+GamL1)*ThW-(GamS+GamS1)*ro*Frac*f_em*xKs*SConc-
     !          (GamG+GamG1)*ThG*Henry)*Conc(jS,i)
        if(.not.lEquil) then
          if((lMobIm(M).or.iDualPor.gt.0).and..not.lBact) then
            q1(i)=q1(i)-Sorb(jS,i)*(ThImob*(GamLi+GamL1i)+
     !                  (1.-Frac)*ro*xKs*SConcS*(GamSi+GamS1i))
            if(lDualNEq) q1(i)=q1(i)-(GamSi+GamS1i)*ro*Sorb2(jS,i)
          else if(.not.lBact) then
            q1(i)=q1(i)-(GamSi+GamS1i)*ro*Sorb(jS,i)
          else if(lBact) then
            q1(i)=q1(i)-ro*(GamSi+GamS1i)*(Sorb(jS,i)+Sorb2(jS,i))
          end if
        end if

*       Velocity corrections
        if(i.eq.1) then
          dx=x(2)-x(1)
          derK=(Henryj-Henry)/dx
        else if(i.eq.NumNP) then
          dx=x(NumNP)-x(NumNP-1)
          derK=(Henry-Henryi)/dx
        else
          dx=(x(j)-x(k))/2.
          derK=(Henryj-Henryi)/dx
        end if
        vCorr(i)=ThG*Dg*TauG*derK
        if(Level.eq.1)      vO(i)=vO(i)-vCorr(i)
        if(Level.eq.NLevel) vN(i)=vN(i)-vCorr(i)

*       Calculate the maximum local Peclet and Courant numbers
        call PeCour(i,j,NumNP,Level,NLevel,lUpW,lArtD,dt,x,v,wc,ThW,vj,
     !              Thj,Disp,Retard,Peclet,Courant,CourMax,PeCr,dtMaxC,
     !              Iter,epsi)
11    continue
      return
      end

************************************************************************

      subroutine MatSet(jS,N,NS,NSD,Level,epsi,alf,dt,kBotCh,kTopCh,
     !                  cBot,cTop,x,thO,thN,vO,vN,Conc,Disp,Retard,wc,
     !                  g0,g1,B,D,E,F,E1,D1,F1,BN,DN,FN,NMat,ChPar,
     !                  TempO,TempN,TDep,dSurf,cAtm,MatNum,lMobIm,
     !                  iDualPor,lVapor,rBot,lBact)

      logical lMobIm(NMat),lVapor,lBact
      double precision B,D,E,F
      dimension cBot(NS),cTop(NS),x(N),thO(N),thN(N),vO(N),vN(N),
     !          Conc(NSD,N),Disp(N),Retard(N),wc(N),g0(N),g1(N),B(N),
     !          D(N),E(N),F(N),ChPar(NSD*16+4,NMat),TempO(N),TempN(N),
     !          TDep(NSD*16+4),MatNum(N)

      do 10 i=1,N
        M=MatNum(i)
        if(lMobIm(M).and.iDualPor.eq.0.or.lBact) then
          ThImob=ChPar(4,M)
          if(ThImob.gt.thO(i)) write(*,*) "Warning !!! ThImob > Theta"
          thN(i)=max(thN(i)-ThImob,0.001)
          thO(i)=max(thO(i)-ThImob,0.001)
        end if
10    continue

*     Lower boundary condition
      b1=x(2)-x(1)
      if(Level.eq.1) then
        F1=           Conc(jS,1)*
     !     (b1/2./dt*thO(1)*Retard(1)+
     !      alf*(-(thO(1)*Disp(1)+thO(2)*Disp(2))/b1/2.-
     !           ((2.+3.*wc(1))*vO(1)+vO(2))/6.+
     !           b1/12.*(3.*g1(1)+g1(2))))+
     !                Conc(jS,2)*
     !     alf*((thO(1)*Disp(1)+thO(2)*Disp(2))/b1/2.-
     !          (vO(1)+(2.-3.*wc(1))*vO(2))/6.+b1/12.*(g1(1)+g1(2)))+
     !                alf*b1/6.*(2.*g0(1)+g0(2))

*       3. type  BC
        if(kBotCh.eq.-1) F(1)=F1+alf*cBot(jS)*vO(1)
      else
        E1=epsi*(-(thN(1)*Disp(1)+thN(2)*Disp(2))/b1/2.+
     !           (vN(1)+(2.-3.*wc(1))*vN(2))/6.-b1/12.*(g1(1)+g1(2)))
        D1=b1/2./dt*thN(1)*Retard(1)+
     !     epsi*((thN(1)*Disp(1)+thN(2)*Disp(2))/b1/2.+
     !           ((2.+3.*wc(1))*vN(1)+vN(2))/6.-b1/12.*(3.*g1(1)+g1(2)))
        F2=epsi*b1/6.*(2.*g0(1)+g0(2))
        F1=F1+F2

*       1.type BC
        if(kBotCh.eq.1) then
          D(1)=1.
          E(1)=0.
          F(1)=cBot(jS)
        end if

*       3. type  BC
        if(kBotCh.eq.-1) then
          if(vN(1).gt.0..or.(lVapor.and.rBot.eq.0.)) then
            E(1)=E1
            D(1)=D1
            F(1)=F(1)+F2+epsi*cBot(jS)*vN(1)
          else
            D(1)=-1.
            E(1)=1.
            F(1)=0.
          end if
        end if

*       Free drainage
        if(kBotCh.eq.0) then
          D(1)=-1.
          E(1)=1.
          F(1)=0.
        end if
      end if

      do 11 i=2,N-1
        a1=b1
        b1=x(i+1)-x(i)
        dx=(x(i+1)-x(i-1))/2.
        if(Level.eq.1) then
          F(i)=       Conc(jS,i-1)*
     !       alf*((thO(i-1)*Disp(i-1)+thO(i)*Disp(i))/a1/2.+
     !            ((2.+3.*wc(i-1))*vO(i-1)+vO(i))/6.+
     !            a1/12.*(g1(i-1)+g1(i)))+
     !                Conc(jS,i)*
     !      (dx/dt*thO(i)*Retard(i)+
     !      alf*(-(thO(i-1)*Disp(i-1)+thO(i)*Disp(i))/a1/2.-
     !           (thO(i+1)*Disp(i+1)+thO(i)*Disp(i))/b1/2.-
     !           (vO(i+1)+3.*(wc(i-1)+wc(i))*vO(i)-vO(i-1))/6.+
     !           (a1*(g1(i-1)+3.*g1(i))+b1*(3.*g1(i)+g1(i+1)))/12.))+
     !                Conc(jS,i+1)*
     !      alf*((thO(i+1)*Disp(i+1)+thO(i)*Disp(i))/b1/2.-
     !           (vO(i)+(2.-3.*wc(i))*vO(i+1))/6.+
     !           b1/12.*(g1(i)+g1(i+1)))+
     !              alf*(a1*(g0(i-1)+2.*g0(i))+b1*(2.*g0(i)+g0(i+1)))/6.
        else
          B(i)=epsi*(-(thN(i-1)*Disp(i-1)+thN(i)*Disp(i))/a1/2.-
     !               ((2.+3.*wc(i-1))*vN(i-1)+vN(i))/6.-
     !               a1/12.*(g1(i-1)+g1(i)))
          D(i)=dx/dt*thN(i)*Retard(i)+
     !         epsi*((thN(i-1)*Disp(i-1)+thN(i)*Disp(i))/a1/2.+
     !               (thN(i+1)*Disp(i+1)+thN(i)*Disp(i))/b1/2.+
     !               (vN(i+1)+3.*(wc(i-1)+wc(i))*vN(i)-vN(i-1))/6.-
     !               (a1*(g1(i-1)+3.*g1(i))+b1*(3.*g1(i)+g1(i+1)))/12.)
          E(i)=epsi*(-(thN(i+1)*Disp(i+1)+thN(i)*Disp(i))/b1/2.+
     !               (vN(i)+(2.-3.*wc(i))*vN(i+1))/6.-
     !               b1/12.*(g1(i)+g1(i+1)))
          F(i)=F(i)+epsi*(a1*(g0(i-1)+2.*g0(i))+
     !                    b1*(2.*g0(i)+g0(i+1)))/6.
        end if
11    continue

*     Upper boundary condition
      if(Level.eq.1) then
        FN=           Conc(jS,N-1)*
     !    alf*((thO(N-1)*Disp(N-1)+thO(N)*Disp(N))/b1/2.+
     !       ((2.+3.*wc(N-1))*vO(N-1)+vO(N))/6.+b1/12.*(g1(N-1)+g1(N)))+
     !                Conc(jS,N)*
     !     (b1/2./dt*thO(N)*Retard(N)+
     !      alf*(-(thO(N-1)*Disp(N-1)+thO(N)*Disp(N))/b1/2.+
     !           (vO(N-1)+(2.-3.*wc(N-1))*vO(N))/6.+
     !           b1/12.*(g1(N-1)+3*g1(N))))+
     !                alf*b1/6.*(g0(N-1)+2.*g0(N))

*       3. type BC
        if(kTopCh.le.0) then
          F(N)=FN
          if(vO(N).lt.0.) F(N)=F(N)-alf*vO(N)*cTop(jS)
          if(kTopCh.eq.-2) then
            M=MatNum(N)
            Tr=293.15
            R=8.314
            jjj=(jS-1)*16
            TT=(TempO(N)+273.15-Tr)/R/(TempO(N)+273.15)/Tr
            Dg=ChPar(jjj+6,M)*exp(TDep(jjj+6)*TT)
            Henry=ChPar(jjj+10,M)*exp(TDep(jjj+10)*TT)
            F(N)=F(N)-alf*Dg/dSurf*Henry*Conc(jS,N)+Dg/dSurf*cAtm
          end if
        end if
      else
        BN=epsi*(-(thN(N-1)*Disp(N-1)+thN(N)*Disp(N))/b1/2.-
     !        ((2.+3.*wc(N-1))*vN(N-1)+vN(N))/6.-b1/12.*(g1(N-1)+g1(N)))
        DN=b1/2./dt*thN(N)*Retard(N)+
     !     epsi*((thN(N-1)*Disp(N-1)+thN(N)*Disp(N))/b1/2.-
     !           (vN(N-1)+(2.-3.*wc(N-1))*vN(N))/6.-
     !           b1/12.*(g1(N-1)+3.*g1(N)))
        FE=epsi*b1/6.*(g0(N-1)+2.*g0(N))
        FN=FN+FE

*       1. type BC
        if(kTopCh.gt.0) then
          B(N)=0.
          D(N)=1.
          F(N)=cTop(jS)

*       3. type BC
        else
          B(N)=BN
          D(N)=DN
          F(N)=F(N)+FE
          if(vN(N).lt.0.) F(N)=F(N)-epsi*vN(N)*cTop(jS)
          if(kTopCh.eq.-2) then
            M=MatNum(N)
            Tr=293.15
            R=8.314
            jjj=(jS-1)*16
            TT=(TempN(N)+273.15-Tr)/R/(TempN(N)+273.15)/Tr
            Dg=ChPar(jjj+6,M)*exp(TDep(jjj+6)*TT)
            Henry=ChPar(jjj+10,M)*exp(TDep(jjj+10)*TT)
            D(N)=D(N)+epsi*Dg/dSurf*Henry
          end if
        end if
      end if

      do 12 i=1,N
        M=MatNum(i)
        if(lMobIm(M).and.iDualPor.eq.0.or.lBact) then
          ThImob=ChPar(4,M)
          thN(i)=thN(i)+ThImob
          thO(i)=thO(i)+ThImob
        end if
12    continue

      return
      end

*************************************************************************

*     Solve matrix equation

      subroutine BanSol(N,A,B,C,F)
      double precision A,B,C,F
      dimension A(N),B(N),C(N),F(N)

      do 11 i=2,N
        B(i)=B(i)-A(i)*C(i-1)/B(i-1)
        F(i)=F(i)-A(i)*F(i-1)/B(i-1)
11    continue
      F(N)=F(N)/B(N)
      do 12 i=2,N
        j=N-i+1
        F(j)=(F(j)-C(j)*F(j+1))/B(j)
12    continue
      return
      end

*************************************************************************

*     Calculate the maximum local Peclet and Courant numbers

      subroutine PeCour(i,j,NumNP,Level,NLevel,lUpW,lArtD,dt,x,v,wc,ThW,
     !                  vj,Thj,Disp,Retard,Peclet,Courant,CourMax,PeCr,
     !                  dtMaxC,Iter,epsi)

      logical lUpW,lArtD
      dimension x(NumNP),Disp(NumNP),Retard(NumNP),wc(NumNP)

      TanH(z)=(exp(z)-exp(-z))/(exp(z)+exp(-z))

      if(i.ne.NumNP) then
        dx=x(j)-x(i)
        vv=0.
        if(ThW.gt.1.e-6.and.Thj.gt.1.e-6) vv=(abs(v)/ThW+abs(vj)/Thj)/2.
        vv1=0.
        if(ThW.gt.1.e-6.and.Thj.gt.1.e-6) vv1=(v/ThW+vj/Thj)/2.
        DD=(Disp(i)+Disp(j))/2.
        if(Level.eq.NLevel) then
          Pec=99999.
          dtMax=1.e+30
c          vMax=amax1(abs(v)/ThW,abs(vj)/Thj)
          vMax=(abs(v)+abs(vj))/(ThW+Thj)
          RMin=amin1(Retard(i),Retard(j))
          if(DD.gt.0.) Pec=abs(vv)*dx/DD
          Cour=vMax*dt/dx/RMin
          Peclet=amax1(Peclet,Pec)
          Courant=amax1(Courant,Cour)
          Cour1=CourMax
          if(.not.lUpW.and..not.lArtD) then
            if(Pec.ne.99999.) Cour1=amin1(1.,PeCr/amax1(0.5,Pec))
          end if
          if(epsi.lt.1..and.vMax.gt.1.e-20) dtMax=Cour1*dx*RMin/vMax
*         the von Neumann time step limit
c          RThE=(ThW+thj)/2.*RMin
c          if(abs(DD).gt.1.e-20)dtMax=amin1(dtMax,10.*RThE*dx*dx/2./DD)
          dtMaxC=amin1(dtMaxC,dtMax)

*       Calculate upstream weighting factors
        else if(lUpW.and.Iter.eq.1) then
          Pe2=11.
          if(DD.gt.0.) Pe2=dx*vv1/DD/2.
          if(abs(vv).lt.1.e-30) then
            wc(i)=0.
          else if(abs(Pe2).gt.10.) then
            if(vv1.gt.0.) wc(i)=1.
            if(vv1.lt.0.) wc(i)=-1
          else
            wc(i)=1./TanH(Pe2)-1./Pe2
            wc(i)=amin1( 1.,wc(i))
            wc(i)=amax1(-1.,wc(i))
          end if
        end if
      end if

      return
      end

*************************************************************************

*     Calculate the dispersion coefficients

      subroutine Disper(i,M,NMat,NSD,NumNP,dt,lTort,lArtD,lUpW,lMobIm,
     !                  iDualPor,Level,NLevel,ChPar,Retard,Disp,thSat,
     !                  ThImob,ThW,ThG,v,Dw,Dg,Henry,ro,Frac,xKs,fExp,
     !                  xNu,cMid,dSConc,PeCr,TauG,iTort,lBact)

      logical lTort,lMobIm(NMat),lArtD,lUpW,lBact
      dimension thSat(NMat),ChPar(NSD*16+4,NMat),Disp(NumNP),
     !          Retard(NumNP)

      if(lTort) then
        ThS=thSat(M)
        if(lMobIm(M).and.iDualPor.eq.0.or.lBact) 
     !                                  ThS=max(thSat(M)-ThImob,0.001)
        if(              iDualPor.gt.0) ThS=    thSat(M)+ThImob
        if(iTort.eq.0) then
          TauW=ThW**(7./3.)/ThS**2
          TauG=ThG**(7./3.)/ThS**2
        else
          TauW=0.66*(ThW/ThS)**(8./3.)
          TauG=ThG**1.5/ThS
        end if
      else
        TauW=1.
        TauG=1.
      end if
      Disp(i)=ChPar(2,M)*abs(v)/ThW+Dw*TauW+ThG/ThW*Dg*Henry*TauG
      if(.not.lArtD.and..not.lUpW) then
        fi=0.
        if(cMid.gt.0.)
     !    fi=6.*ThW*ro*xKs*cMid**(fExp-1.)*
     !         (fExp/(1.+xNu*cMid**fExp)**2-1./(1.+xNu*cMid**fExp))
        DPom=amax1(dt/(6.*ThW*(ThW+ro*Frac*xKs*dSConc+ThG*Henry)+fi),0.)
        if(Level.ne.NLevel) then
          Disp(i)=Disp(i)+v*v*DPom
        else
          Disp(i)=amax1(Disp(i)-v*v*DPom,Disp(i)/2.)
        end if
      end if
      if(lArtD) then
        DD=0.
        if(PeCr.ne.0.and.abs(v).gt.1.e-15) DD=v*v*dt/thW/thW/
     !       Retard(i)/PeCr
        if(DD.gt.Disp(i)) Disp(i)=DD
      end if

      return
      end

*************************************************************************

*     Calculate the adsorbed concentration on kinetic sites or
*     the concentration in an imobile zone, before solving matrix equation

      subroutine NEquil(i,jS,NSD,NumNP,NMat,M,Conc,Sorb,Sorb2,cNew,
     !                  SorbN,SorbN2,SSorb,SSorb2,lMobIm,lLinear,lBact,
     !                  Level,NLevel,dt,epsi,ro,xKs,xKsO,cc,SConc,
     !                  SConcO,SConcS,SConcOS,dSConcS,xMuL,xMuLO,xMuS,
     !                  xMuSO,dRetardS,GamL,GamL1,GamLO,GamL1O,GamS,
     !                  GamS1,GamSO,GamS1O,Omega,OmegaO,rKa1,rKa1O,rKa2,
     !                  rKa2O,rKd1,rKd1O,rKd2,rKd2O,ThW,ThWO,psi1,psi1O,
     !                  psi2,psi2O,DMobI,Frac,ThImob,ThImobO,iDualPor,
     !                  SinkIm,FlMacro,lNEquil,GamL1Pi,GamS1Pi,xKsP,
     !                  SConcPS,lDualNEq,f_em,OmegaS)

      logical lMobIm(NMat),lLinear(NSD),lBact,lNEquil,lDualNEq
      dimension Conc(NSD,NumNP),Sorb(NSD,NumNP),Sorb2(NSD,NumNP),
     !          SinkIm(NumNP),cNew(NumNP),SorbN(NumNP),SorbN2(NumNP)


      FlMacro=0.
      if(iDualPor.gt.0) then            ! mobile-immobile model
        if(SinkIm(i).gt.0) then
          FlMacro=SinkIm(i)*Conc(jS,i)
        else
          FlMacro=SinkIm(i)*Sorb(jS,i)
        end if
      end if
      SSorb =Sorb (jS,i)
      SSorb2=Sorb2(jS,i)
      if(Level.eq.NLevel) then

*       mobile-immobile model
        if((lMobIm(M).or.iDualPor.gt.0).and..not.lBact) then
          AMobI=(ThImob+ThImobO)/2.+(1.-Frac)*ro*xKs*dSConcS
          dTheta=ThImob-ThImobO
          EMobI =ThImob*xMuL+(1.-Frac)*ro*xMuS-(1.-Frac)*ro*dRetardS
          EMobIO=ThImobO*xMuLO+(1.-Frac)*ro*xMuSO-(1.-Frac)*ro*dRetardS
          EMobI =EMobI +FlMacro
          EMobIO=EMobIO+FlMacro
          if(lNEquil) then
            cS=Sorb(jS-1,i)*(ThImob*GamL1Pi+
     !                               (1.-Frac)*ro*GamS1Pi*xKsP*SConcPS)
            EMobI =EMobI +cS
            EMobIO=EMobIO+cS
          end if
          BMobI =ThImob *(GamL +GamL1 )+
     !                       (1.-Frac)*ro*(GamS +GamS1 )*xKs *SConcS
          BMobIO=ThImobO*(GamLO+GamL1O)+
     !                       (1.-Frac)*ro*(GamSO+GamS1O)*xKsO*SConcOS
          if(lLinear(jS)) then
            DMobI=  2.*AMobI+dt*(Omega +BMobI )+dTheta
            GMobI0=(2.*AMobI-dt*(OmegaO+BMobIO)-dTheta)/DMobI
            Sorb(jS,i)=Sorb(jS,i)*GMobI0+
     !                      dt*(OmegaO*Conc(jS,i)+EMobI+EMobIO)/DMobI
            SSorb=Sorb(jS,i)
          else
            SorbN(i)=Sorb(jS,i)+dt/AMobI*  (   epsi *(Omega *
     !              (cNew(i)   -SorbN(i))  -BMobI *SorbN(i)  +EMobI)+
     !                                     (1.-epsi)*(OmegaO*
     !              (Conc(jS,i)-Sorb(jS,i))-BMobIO*Sorb(jS,i)+EMobIO))
            SSorb=SorbN(i)
          end if
          if(lDualNEq) then
            cS=0.
            if(lNEquil) cS=GamS1Pi*Sorb2(jS-1,i)*dt*2.
            if(lLinear(jS)) then
              Sorb2(jS,i)=((2.-(OmegaS+GamSO+GamS1O)*dt)*Sorb2(jS,i)+
     !                  dt*Frac*(1.-f_em)*OmegaS*xKsO*Conc(jS,i)+
     !                  dt*(1.-f_em)*(xMuSO+xMuS)+cS)/
     !                  (2.+dt*(OmegaS+GamS+GamS1))
              SSorb2=Sorb2(jS,i)
            else
              SorbN2(i)=Sorb2(jS,i)+dt*
     !    (epsi* (OmegaS*(Frac*(1.-f_em)*SConc *xKs *cc     -SorbN2(i))-
     !                         (GamS+GamS1)*SorbN2(i)+(1.-f_em)*xMuS)+
     ! (1.-epsi)*(OmegaS*(Frac*(1.-f_em)*SConcO*xKsO*Conc(jS,i)-SSorb2)-
     !                         (GamSO+GamS1O)*SSorb2+(1.-f_em)*xMuSO))
              SSorb2=SorbN2(i)
            end if
          end if

*       two-site sorption model
        else if(.not.lBact) then
          cS=0.
          if(lNEquil) cS=GamS1Pi*Sorb(jS-1,i)*dt*2.
          if(lLinear(jS)) then
            Sorb(jS,i)=((2.-(OmegaO+GamSO+GamS1O)*dt)*Sorb(jS,i)+
     !                  dt*(1.-Frac)*OmegaO*xKsO*Conc(jS,i)+
     !                  dt*(1.-Frac)*(xMuSO+xMuS)+cS)/
     !                  (2.+dt*(Omega+GamS+GamS1))
            SSorb=Sorb(jS,i)
          else
            SorbN(i)=Sorb(jS,i)+dt*
     !          (epsi* (Omega* ((1.-Frac)*SConc *xKs *cc     -SorbN(i))-
     !                         (GamS+GamS1)*SorbN(i)+(1.-Frac)*xMuS)+
     !       (1.-epsi)*(OmegaO*((1.-Frac)*SConcO*xKsO*Conc(jS,i)-SSorb)-
     !                         (GamSO+GamS1O)*SSorb+(1.-Frac)*xMuSO))
            SSorb=SorbN(i)
          end if

*       filtration model
        else if(lBact) then
          if(lLinear(jS)) then
            Sorb(jS,i)=((2.-dt*(rKd1O+GamSO+GamS1O))*Sorb(jS,i)+
     !                      dt*rKa1O*ThW*Conc(jS,i)/ro)/
     !                      (2.+dt*(rKd1+GamS+GamS1))
            SSorb=Sorb(jS,i)
            Sorb2(jS,i)=((2.-dt*(rKd2O+GamSO+GamS1O))*Sorb2(jS,i)+
     !                      dt*rKa2O*ThW*Conc(jS,i)/ro)/
     !                      (2.+dt*(rKd2+GamS+GamS1))
            SSorb2=Sorb2(jS,i)
          else
            SorbN(i)=Sorb(jS,i)+dt*
     !               (epsi*    (rKa1*ThW/ro*psi1*cc-
     !                         (rKd1+GamS+GamS1)*SorbN(i))+
     !               (1.-epsi)*(rKa1O*ThWO/ro*psi1O*Conc(jS,i)-
     !                         (rKd1O+GamSO+GamS1O)*Sorb(jS,i)))
            SSorb=SorbN(i)
            SorbN2(i)=Sorb2(jS,i)+dt*
     !               (epsi*    (rKa2*ThW/ro*psi2*cc-
     !                         (rKd2+GamS+GamS1)*SorbN2(i))+
     !               (1.-epsi)*(rKa2O*ThWO/ro*psi2O*Conc(jS,i)-
     !                         (rKd2O+GamSO+GamS1O)*Sorb2(jS,i)))
            SSorb2=SorbN2(i)
          end if
        end if
      end if

      return
      end

*************************************************************************

*     Calculate sorbed concentration for linear noneq. adsorption or 
*     concentration in the imobile water. 
*     At the end of the time step. After solving matrix equation

      subroutine SorbConc(jS,NSD,N,MatNum,Temp,lMobIm,ChPar,TDep,Sorb,
     !                    Conc,dt,NMat,lBact,ThW,Sorb2,lFiltr,Veloc,
     !                    iDualPor,ThIm,ThOIm,SinkIm,STrans,iMoistDep,
     !                    NMatD,DMoist,WDep,xConv,tConv,lDualNEq)

      logical lMobIm(NMat),lBact,lFiltr,lDualNEq
      dimension MatNum(N),Temp(N),ChPar(NSD*16+4,NMat),TDep(NSD*16+4),
     !          Conc(NSD,N),Sorb(NSD,N),ThW(N),Sorb2(NSD,N),Veloc(N),
     !          ThIm(N),ThOIm(N),SinkIm(N),STrans(N),
     !          DMoist(NMatD,NSD,13,6),WDep(2+NMatD,NSD*9)

      Tr=293.15
      R=8.314
      jjj=(jS-1)*16
      do 11 i=1,N
        M=MatNum(i)
        TT=(Temp(i)+273.15-Tr)/R/(Temp(i)+273.15)/Tr
        Frac  =ChPar(3,     M)*exp(TDep(3)     *TT)
        xKs   =ChPar(jjj+ 7,M)*exp(TDep(jjj+ 7)*TT)
        f1=1.
        if(iMoistDep.gt.0) 
     !    f1=rMD(NMatD,NSD,M,jS,11,DMoist,2,WDep,ThW(i),iMoistDep)
        GamS  =ChPar(jjj+12,M)*exp(TDep(jjj+12)*TT)*f1
        if(iMoistDep.gt.0) 
     !    f1=rMD(NMatD,NSD,M,jS,13,DMoist,5,WDep,ThW(i),iMoistDep)
        GamS1 =ChPar(jjj+15,M)*exp(TDep(jjj+15)*TT)*f1
        Omega =ChPar(jjj+20,M)*exp(TDep(jjj+20)*TT)

        if((lMobIm(M).or.iDualPor.gt.0).and..not.lBact) then ! mobile-immobile model
          ro    =ChPar(1,     M)*exp(TDep(1)     *TT)
          if(iDualPor.eq.0) then
            ThImob=ChPar(4,M)*exp(TDep(4)*TT)
            ThImobO=ThImob
          else if(iDualPor.gt.0) then
            ThImob =ThIm(i)
            ThImobO=ThOIm(i)
          end if
          if(iMoistDep.gt.0) 
     !      f1=rMD(NMatD,NSD,M,jS,10,DMoist,1,WDep,ThImob,iMoistDep)
          GamL  =ChPar(jjj+11,M)*exp(TDep(jjj+11)*TT)*f1
          if(iMoistDep.gt.0) 
     !      f1=rMD(NMatD,NSD,M,jS,12,DMoist,4,WDep,ThImob,iMoistDep)
          GamL1 =ChPar(jjj+14,M)*exp(TDep(jjj+14)*TT)*f1
          if(iMoistDep.gt.0) 
     !      f1=rMD(NMatD,NSD,M,jS,11,DMoist,2,WDep,ThImob,iMoistDep)
          GamS  =ChPar(jjj+12,M)*exp(TDep(jjj+12)*TT)*f1
          if(iMoistDep.gt.0) 
     !      f1=rMD(NMatD,NSD,M,jS,13,DMoist,5,WDep,ThImob,iMoistDep)
          GamS1 =ChPar(jjj+15,M)*exp(TDep(jjj+15)*TT)*f1
          dTheta=ThImob-ThImobO
          AMobI=(ThImob+ThImobO)/2.+(1.-Frac)*ro*xKs
          BMobI=ThImob*(GamL+GamL1)+(GamS+GamS1)*ro*(1.-Frac)*xKs
          DMobI=2.*AMobI+dt*(Omega+BMobI)+dTheta
          Sorb(jS,i)=Sorb(jS,i)+dt*Omega*Conc(jS,i)/DMobI
          FlMacro=0.
          if(iDualPor.gt.0) then
            if(SinkIm(i).gt.0) then
              FlMacro=SinkIm(i)*Conc(jS,i)
            else
              FlMacro=SinkIm(i)*Sorb(jS,i)
            end if
          end if
          if(jS.eq.1) STrans(i)=Omega*(Conc(jS,i)-Sorb(jS,i))+FlMacro
          if(lDualNEq) then
            f_em  =ChPar(jjj+13,M)*exp(TDep(jjj+13)*TT)
            OmegaS=ChPar(jjj+16,M)*exp(TDep(jjj+16)*TT)
            Sorb2(jS,i)=
     !              Sorb2(jS,i)+dt*OmegaS*Frac*(1.-f_em)*xKs*Conc(jS,i)/
     !                         (2.+dt*(OmegaS+GamS+GamS1))
          end if

        else if(.not.lBact) then                 ! two-site sorption model
          Sorb(jS,i)=Sorb(jS,i)+dt*Omega*(1.-Frac)*xKs*Conc(jS,i)/
     !                       (2.+dt*(Omega+GamS+GamS1))

        else if(lBact) then                      ! filtration model
          ro    =ChPar(1,     M)*exp(TDep(1)     *TT)
          ThImob=ChPar(4,M)
          Theta=ThW(i)-ThImob
          GamS  =ChPar(jjj+12,M)*exp(TDep(jjj+12)*TT)
          GamS1 =0.
          rKa1  =ChPar(jjj+19,M)*exp(TDep(jjj+19)*TT)
          rKd1  =ChPar(jjj+20,M)*exp(TDep(jjj+20)*TT)
          rKa2  =ChPar(jjj+16,M)*exp(TDep(jjj+16)*TT)
          rKd2  =ChPar(jjj+17,M)*exp(TDep(jjj+17)*TT)
          if(lFiltr) then
            Dc=ChPar(jjj+13,M)*exp(TDep(jjj+13)*TT)
            Dp=ChPar(jjj+14,M)*exp(TDep(jjj+14)*TT)
            Alfa1=rKa1
            Alfa2=rKa2
            call Deposit(rKa1,rKa2,Dc,Dp,Alfa1,Alfa2,Theta,Veloc(i),
     !                   Temp(i),xConv,tConv)
          end if
          Sorb(jS,i) =Sorb(jS,i) +dt*rKa1*Theta*Conc(jS,i)/ro/
     !                          (2.+dt*(rKd1+GamS+GamS1))
          Sorb2(jS,i)=Sorb2(jS,i)+dt*rKa2*Theta*Conc(jS,i)/ro/
     !                          (2.+dt*(rKd2+GamS+GamS1))
        end if
11    continue

      return
      end

*************************************************************************

*     Calculate mass-transfer fluxes at the end of the time interval

      subroutine MassTran(jS,NS,NSD,N,MatNum,Temp,lMobIm,lEquil,ChPar,
     !                    TDep,Sorb,Conc,NMat,x,cvCh0,cvCh1,cvChR,
     !                    cvChIm,epsi,q0,q1,sSink,lBact,theta,Sorb2,
     !                    lFiltr,Veloc,iDualPor,SinkIm,xConv,tConv,
     !                    lDualNEq,STrans,lLinear)

      logical lEquil,lMobIm(NMat),lBact,lFiltr,lDualNEq,lLinear(NS)
      dimension x(N),ChPar(NSD*16+4,NMat),MatNum(N),Conc(NSD,N),
     !          cvCh0(NS),cvCh1(NS),Temp(N),TDep(NSD*16+4),sSink(N),
     !          Sorb(NSD,N),q0(N),q1(N),cvChR(NS),cvChIm(NS),theta(N),
     !          Sorb2(NSD,N),Veloc(N),SinkIm(N),STrans(N)

      Tr=293.15
      R=8.314
      jjj=(jS-1)*16
      do 11 i=1,N
        Mi=MatNum(i)
        TTi=(Temp(i)+273.15-Tr)/R/(Temp(i)+273.15)/Tr
        if(i.eq.N) goto 10
        j=i+1
        dx=x(j)-x(i)
        cvCh0(jS)=cvCh0(jS)+epsi*dx*(q0(i)+q0(j))/2.
        cvCh1(jS)=cvCh1(jS)+epsi*dx*(q1(i)+q1(j))/2.
        cvChR(jS)=cvChR(jS)+epsi*dx*(sSink(i)+sSink(j))/2.
        if(.not.lEquil) then
          Mj=MatNum(j)
          TTj=(Temp(j)+273.15-Tr)/R/(Temp(j)+273.15)/Tr
          Omegai=ChPar(jjj+20,Mi)*exp(TDep(jjj+20)*TTi)
          Omegaj=ChPar(jjj+20,Mj)*exp(TDep(jjj+20)*TTj)

*         mobile-immobile model
          if((lMobIm(Mi).or.iDualPor.gt.0).and..not.lBact) then
            FlMacroi=0.
            FlMacroj=0.
            if(iDualPor.gt.0) then
              if(SinkIm(i).gt.0) then
                FlMacroi=SinkIm(i)*Conc(jS,i)
              else
                FlMacroi=SinkIm(i)*Sorb(jS,i)
              end if
              if(SinkIm(j).gt.0) then
                FlMacroj=SinkIm(j)*Conc(jS,j)
              else
                FlMacroj=SinkIm(j)*Sorb(jS,j)
              end if
            end if
            cvChIm(jS)=cvChIm(jS)+epsi*dx/2.*
     !                  (Omegai*(Conc(jS,i)-Sorb(jS,i))+
     !                   Omegaj*(Conc(jS,j)-Sorb(jS,j))+
     !                   FlMacroi+FlMacroj)
            if(jS.eq.1.and..not.lLinear(jS)) 
     !          STrans(i)=epsi*(Omegai*(Conc(jS,i)-Sorb(jS,i))+
     !                          FlMacroi)
            if(lDualNEq) then
              roi   =ChPar(1,     Mi)*exp(TDep(1)     *TTi)
              Fraci =ChPar(3,     Mi)*exp(TDep(3)     *TTi)
              xKsi  =ChPar(jjj+ 7,Mi)*exp(TDep(jjj+ 7)*TTi)
              xNui  =ChPar(jjj+ 8,Mi)*exp(TDep(jjj+ 8)*TTi)
              fExpi =ChPar(jjj+ 9,Mi) !*exp(TDep(jjj+ 9)*TTi)
              roj   =ChPar(1,     Mj)*exp(TDep(1)     *TTj)
              Fracj =ChPar(3,     Mj)*exp(TDep(3)     *TTj)
              xKsj  =ChPar(jjj+ 7,Mj)*exp(TDep(jjj+ 7)*TTj)
              xNuj  =ChPar(jjj+ 8,Mj)*exp(TDep(jjj+ 8)*TTj)
              fExpj =ChPar(jjj+ 9,Mj) !*exp(TDep(jjj+ 9)*TTj)
              f_emi =ChPar(jjj+13,Mi)*exp(TDep(jjj+13)*TTi)
              f_emj =ChPar(jjj+13,Mj)*exp(TDep(jjj+13)*TTj)
              OmegaSi=ChPar(jjj+16,Mi)*exp(TDep(jjj+16)*TTi)
              OmegaSj=ChPar(jjj+16,Mj)*exp(TDep(jjj+16)*TTj)
              cci   =Conc(jS,i)
              ccj   =Conc(jS,j)
              SorbEi=0.
              SorbEj=0.
              if(cci.gt.0.) SorbEi=
     !             Fraci*(1.-f_emi)*xKsi*cci**fExpi/(1.+xNui*cci**fExpi)
              if(ccj.gt.0.) SorbEj=
     !             Fracj*(1.-f_emj)*xKsj*ccj**fExpj/(1.+xNuj*ccj**fExpj)
              cvChIm(jS)=cvChIm(jS)+epsi*dx/2.*
     !                      (roi*OmegaSi*(SorbEi-Sorb2(jS,i))+
     !                       roj*OmegaSj*(SorbEj-Sorb2(jS,j)))
              if(jS.eq.1.and..not.lLinear(jS)) 
     !         STrans(i)=STrans(i)+epsi*roi*OmegaSi*(SorbEi-Sorb2(jS,i))
            end if

*         two-site sorption model
      	  else if(.not.lBact) then
            roi   =ChPar(1,     Mi)*exp(TDep(1)     *TTi)
            Fraci =ChPar(3,     Mi)*exp(TDep(3)     *TTi)
            xKsi  =ChPar(jjj+ 7,Mi)*exp(TDep(jjj+ 7)*TTi)
            xNui  =ChPar(jjj+ 8,Mi)*exp(TDep(jjj+ 8)*TTi)
            fExpi =ChPar(jjj+ 9,Mi) !*exp(TDep(jjj+ 9)*TTi)
            roj   =ChPar(1,     Mj)*exp(TDep(1)     *TTj)
            Fracj =ChPar(3,     Mj)*exp(TDep(3)     *TTj)
            xKsj  =ChPar(jjj+ 7,Mj)*exp(TDep(jjj+ 7)*TTj)
            xNuj  =ChPar(jjj+ 8,Mj)*exp(TDep(jjj+ 8)*TTj)
            fExpj =ChPar(jjj+ 9,Mj) !*exp(TDep(jjj+ 9)*TTj)
            cci   =Conc(jS,i)
            ccj   =Conc(jS,j)
            SorbEi=0.
            SorbEj=0.
            if(cci.gt.0.)
     !          SorbEi=(1.-Fraci)*xKsi*cci**fExpi/(1.+xNui*cci**fExpi)
            if(ccj.gt.0.)
     !          SorbEj=(1.-Fracj)*xKsj*ccj**fExpj/(1.+xNuj*ccj**fExpj)
            cvChIm(jS)=cvChIm(jS)+epsi*dx/2.*
     !                      (roi*Omegai*(SorbEi-Sorb(jS,i))+
     !                       roj*Omegaj*(SorbEj-Sorb(jS,j)))
            if(jS.eq.1) 
     !        STrans(i)=epsi*roi*Omegai*(SorbEi-Sorb(jS,i))

*         filtration model
      	  else if(lBact) then
            roi   =ChPar(1,     Mi)*exp(TDep(1)     *TTi)
            roj   =ChPar(1,     Mj)*exp(TDep(1)     *TTj)
            SMax1i=ChPar(jjj+18,Mi)*exp(TDep(jjj+18)*TTi)
            SMax1j=ChPar(jjj+18,Mj)*exp(TDep(jjj+18)*TTj)
            rKa1i =ChPar(jjj+19,Mi)*exp(TDep(jjj+19)*TTi)
            rKa1j =ChPar(jjj+19,Mj)*exp(TDep(jjj+19)*TTj)
            rKd1i =ChPar(jjj+20,Mi)*exp(TDep(jjj+20)*TTi)
            rKd1j =ChPar(jjj+20,Mj)*exp(TDep(jjj+20)*TTj)
            SMax2i=ChPar(jjj+15,Mi)*exp(TDep(jjj+15)*TTi)
            SMax2j=ChPar(jjj+15,Mj)*exp(TDep(jjj+15)*TTj)
            rKa2i =ChPar(jjj+16,Mi)*exp(TDep(jjj+16)*TTi)
            rKa2j =ChPar(jjj+16,Mj)*exp(TDep(jjj+16)*TTj)
            rKd2i =ChPar(jjj+17,Mi)*exp(TDep(jjj+17)*TTi)
            rKd2j =ChPar(jjj+17,Mj)*exp(TDep(jjj+17)*TTj)
            ThImobi=ChPar(4,Mi)
            ThImobj=ChPar(4,Mj)
            Thetai=theta(i)-ThImobi
            Thetaj=theta(j)-ThImobj
            iPsi1=0
            iPsi2=0
            if(.not.lFiltr) iPsi2=int(ChPar(jjj+13,Mi))
            if(.not.lFiltr) iPsi1=int(ChPar(jjj+14,Mj))
            if(iPsi1.eq.0.and.SMax1i.gt.0.) iPsi1=1
            if(iPsi2.eq.0.and.SMax2i.gt.0.) iPsi2=1
            if(iPsi1.ge.3.or.iPsi2.ge.3) 
     !        Dc=ChPar(jjj+6,Mi)*exp(TDep(jjj+6)*TTi)
            if(iPsi1.eq.5.or.iPsi2.eq.5) aa=ChPar(jjj+15,Mi)
            psi1i=1.
            psi1j=1.
            psi2i=1.
            psi2j=1.
            if(iPsi1.gt.0) then
              call Blocking(iPsi1,SMax1i,psi1i,x(i),Sorb(jS,i),Dc,aa)
              call Blocking(iPsi1,SMax1j,psi1j,x(j),Sorb(jS,j),Dc,aa)
            end if
            if(iPsi2.gt.0) then
              call Blocking(iPsi2,SMax2i,psi2i,x(i),Sorb2(jS,i),Dc,aa)
              call Blocking(iPsi2,SMax2j,psi2j,x(j),Sorb2(jS,j),Dc,aa)
            end if
            if(lFiltr) then
              Dc=ChPar(jjj+13,Mi)*exp(TDep(jjj+13)*TTi)
              Dp=ChPar(jjj+14,Mi)*exp(TDep(jjj+14)*TTi)
              Alfa1=rKa1i
              Alfa2=rKa2i
              call Deposit(rKa1,rKa2,Dc,Dp,Alfa1,Alfa2,Thetai,
     !                     Veloc(i),Temp(i),xConv,tConv)
              rKa1i=rKa1
              rKa1j=rKa1
              rKa2i=rKa2
              rKa2j=rKa2
            end if
            cvChIm(jS)=cvChIm(jS)+epsi*dx/2.*
     !            (Conc(jS,i)*Thetai*(psi1i*rKa1i+psi2i*rKa2i)+
     !             Conc(jS,j)*Thetaj*(psi1j*rKa1j+psi2j*rKa2j)-
     !             roi*(Sorb(jS,i)*rKd1i+Sorb2(jS,i)*rKd2i)-
     !             roj*(Sorb(jS,j)*rKd1j+Sorb2(jS,j)*rKd2j))
            if(jS.eq.1) STrans(i)=epsi*
     !                     (Conc(jS,i)*Thetai*(psi1i*rKa1i+psi2i*rKa2i)-
     !                      roi*(Sorb(jS,i)*rKd1i+Sorb2(jS,i)*rKd2i))
       	  end if
        end if
10      continue
11    continue

      return
      end

*************************************************************************

*     Calculate flux concentration for the first solute

      subroutine FluxConc(NumNP,NMat,NSD,x,v,theta,thSat,ChPar,MatNum,
     !                    Temp,TDep,Conc,ConcF,lTort,lMobIm,iDualPor,
     !                    ThIm,jS)

      logical lTort,lMobIm(NMat)
      dimension x(NumNP),v(NumNP),theta(NumNP),thSat(NMat),ConcF(NumNP),
     !          ChPar(NSD*16+4,NMat),Conc(NSD,NumNP),MatNum(NumNP),
     !          Temp(NumNP),TDep(NSD*16+4),ThIm(NumNP)

      jjj=(jS-1)*16
      Tr=293.15
      R=8.314

      do 11 i=1,NumNP
        M=MatNum(i)
        ThW=Theta(i)
        ThG=amax1(0.,thSat(M)-ThW)
        if(lMobIm(M)) then
          if(iDualPor.eq.0) then
            ThImob=ChPar(4,M)
            ThW=max(ThW-ThImob,0.001)
          else if(iDualPor.gt.0) then
            ThImob=ThIm(i)
          end if
        end if
        TT=(Temp(i)+273.15-Tr)/R/(Temp(i)+273.15)/Tr
        Dw   =ChPar(jjj+ 5,M)*exp(TDep(jjj+ 5)*TT)
        Dg   =ChPar(jjj+ 6,M)*exp(TDep(jjj+ 6)*TT)
        Henry=ChPar(jjj+10,M)*exp(TDep(jjj+10)*TT)
        if(lTort) then
          ThS=thSat(M)
          if(lMobIm(M).and.iDualPor.eq.0) ThS=max(thSat(M)-ThImob,0.001)
          if(              iDualPor.gt.0) ThS=    thSat(M)+ThImob
          TauW=ThW**(7./3.)/ThS**2
          TauG=ThG**(7./3.)/ThS**2
        else
          TauW=1.
          TauG=1.
        end if
        qW=v(i)
        Disp=ChPar(2,M)*abs(qW)/ThW+Dw*TauW+ThG/ThW*Dg*Henry*TauG
        cGrad=0.
        if(i.eq.1) then
          cGrad=(Conc(jS,i+1)-Conc(jS,i))/(x(i+1)-x(i))
        else if(i.eq.NumNP) then
          cGrad=(Conc(jS,i)-Conc(jS,i-1))/(x(i)-x(i-1))
        else
          cGrad=(Conc(jS,i+1)-Conc(jS,i-1))/(x(i+1)-x(i-1))
        end if
        ConcF(i)=Conc(jS,i)
        if(qW.ne.0.) ConcF(i)=Conc(jS,i)-Disp*ThW/qW*cGrad
11    continue

      return
      end

*************************************************************************

*     Calculate blocking coefficient for the attachment process

      subroutine Blocking(iPsi,SMax,psi,x,ss,Dc,SMax2)

      real Minf

      psi=1.
      if     (iPsi.eq.1) then
        if(SMax.gt.0.) psi=1.-ss/SMax
      else if(iPsi.eq.2) then
        if(SMax.gt.0.) psi=max(ss**SMax,psi)
      else if(iPsi.eq.3) then
        Binf=1./SMax
        Sinf=.546
        Minf=Dc
        const=Sinf*Binf*ss  
        if(ss.le.(0.8*SMax))
     !    psi=1.-(4.*const)+(3.08*const**2.)+(1.4069*const**3.)
        if(ss.gt.(0.8*SMax))
     !    psi=((1.-Binf*ss)**3.)/(2.*(Minf**2.)*(Binf**3.))
      else if(iPsi.eq.4) then
        if(SMax.gt.0..and.Dc.gt.0.) psi=((abs(x)+Dc)/Dc)**(-SMax)
      else if(iPsi.eq.5) then
        if(SMax.gt.0..and.Dc.gt.0.) psi=((abs(x)+Dc)/Dc)**(-SMax)
        if(SMax2.gt.0.) psi=psi*(1.-ss/SMax2)
      end if

      return
      end

*************************************************************************

*     Calculate the deposition coefficient for the bacteria transport,
*     All calculations within this subroutines are in meters and seconds
*     Conversions are needed

      subroutine Deposit(Ka1,Ka2,Dc1,Dp1,Alfa1,Alfa2,Theta,q,Temp,xConv,
     !                   tConv)

      real Ka1,Ka2,mu,N_Pe,N_Lo,N_R,N_G

*     Ka       - deposition coefficient (output) [1/T]
*     Dc       - diameter of the sand grains (m)
*     Dp       - diameter of the bacteria (0.95 microm) (m)
*     Alfa     - sticking efficiency (-)
*     Theta    - porosity (-)
*     q        - Darcy flux [L/T]
*     Temp     - Temperature in Celcius

      Dc=Dc1/xConv
      Dp=Dp1/xConv
      if(Dp.le.0.and.Dc.le.0.) then
        write(*,*) 'Both Dp and Dc are equal to zero !!!'
        write(*,*) 'Press Enter to continue'
        read(*,*)
        stop
      end if
      PI=3.1415                ! Ludolf's number
      mu=0.00093               ! fluid viscosity (Pa s)
      Bk=1.38048e-23           ! Boltzman constatnt (J/K)
      H=1.e-20                 ! Hamaker constant (J)
      g=9.81                   ! gravitational acceleration (m/s2)
      rop=1080.                ! bacterial density (kg/m3)
      rof=998.                 ! fluid density (kg/m3)
      Veloc=abs(q/xConv*tConv) ! absolute value of Darcy flux (converted to m/s)
      PVeloc=Veloc/Theta       ! pore velocity (converted to m/s)
      Dc=Dc1/xConv             ! conversion to m
      Dp=Dp1/xConv             ! conversion to m
      
      if(Veloc.gt.0.) then
        gamma=(1.-Theta)**(1./3.)
        As=2.*(1.-gamma**5)/(2.-3.*gamma+3.*gamma**5-2.*gamma**6) ! Correct.factor
        N_Pe=3.*PI*mu*Dp*Dc*Veloc/(Bk*(Temp+273.15)) ! Peclet number
        e_diff=4.*As**(1./3.)*N_Pe**(-2./3.)         ! removal by diffusion

        N_Lo=4.*H/(9.*PI*mu*Dp**2*Veloc)             ! London number
        N_R=Dp/Dc                                    ! Interception number
        e_inter=As*N_Lo**(1./8.)*N_R**(15./8.)       ! removal interception

        N_G=g*(rop-rof)*Dp**2/(18.*mu*Veloc)         ! Gravitation number
        e_grav=0.00338*As*N_G**1.2*N_R**(-0.4)       ! removal by gravitational
*                                                      sedimentation
      else
        e_diff =0.
        e_inter=0.
        e_grav =0.
      end if
      eta=e_diff+e_inter+e_grav                      ! single-collector efficiency

*     Original Filtration Theory
      Ka1=3.*(1.-Theta)/2./dc*eta*Alfa1*PVeloc
      Ka1=Ka1/tConv
      Ka2=3.*(1.-Theta)/2./dc*eta*Alfa2*PVeloc
      Ka2=Ka2/tConv

      return
      end

*************************************************************************

*     Nonequilibrium phase is initially in equilibrium with liquid phase

      subroutine NonEqInit(NumNP,NSD,NS,NMat,MatNum,TDep,Temp,ChPar,
     !                     Conc,Sorb,lLinear,lMobIm,iDualPor,lBact,
     !                     Sorb2,Theta)

      logical lLinear(NSD),lMobIm(NMat),lBact
      dimension ChPar(NSD*16+4,NMat),MatNum(NumNP),TDep(NSD*16+4),
     !          Conc(NSD,NumNP),Sorb(NSD,NumNP),Sorb2(NSD,NumNP),
     !          Temp(NumNP),Theta(NumNP)

      Tr=293.15
      R=8.314
      do 12 jS=1,NS
        jjj=(jS-1)*16
        do 11 i=1,NumNP
          M=MatNum(i)
          TT=(Temp(i)+273.15-Tr)/R/(Temp(i)+273.15)/Tr
          if(lMobIm(M).or.iDualPor.gt.0) then
            Sorb(jS,i)=Conc(jS,i)
          else
            ro   =ChPar(1,     M)*exp(TDep(1)     *TT)
            Frac =ChPar(3,     M)*exp(TDep(3)     *TT)
            xKs  =ChPar(jjj+ 7,M)*exp(TDep(jjj+ 7)*TT)
            xNu  =ChPar(jjj+ 8,M)*exp(TDep(jjj+ 8)*TT)
            fExp =ChPar(jjj+ 9,M) !*exp(TDep(jjj+ 9)*TT)
            SConc=1.
            cc=Conc(jS,i)
            if(.not.lLinear(jS).and.cc.gt.0.) 
     !        SConc=cc**(fExp-1.)/(1.+xNu*cc**fExp)
            Sorb(jS,i)=(1.-Frac)*SConc*xKs*cc
            if(lBact) then
              Frac=0.
              rKa1=ChPar(jjj+19,M)*exp(TDep(jjj+19)*TT)
              rKd1=ChPar(jjj+20,M)*exp(TDep(jjj+20)*TT)
              xKs1=Theta(i)*rKa1/ro/rKd1 
              Sorb(jS,i)=xKs1*Conc(jS,i)
              rKa2=ChPar(jjj+16,M)*exp(TDep(jjj+16)*TT)
              rKd2=ChPar(jjj+17,M)*exp(TDep(jjj+17)*TT)
              xKs2=Theta(i)*rKa1/ro/rKd1 
              Sorb2(jS,i)=xKs2*Conc(jS,i)
            end if
          end if
11      continue
12    continue

      return
      end

************************************************************************

*     Subroutine calculating accesible water content for colloids and colloid velocity
*     Only for steady-state water flow and homogeneous soil profile

      subroutine Exclusion(NumNP,NMat,NSD,Par,ChPar,ThNew,vNew,ThOld,
     !                      vOld)

      double precision thr,ths,swr,alpha,vgn1,vgn2,vgm1,vgm2,Pc_c,
     !                 sw_c_eff,sw_c,sw,sw_eff,kr_c,krw,r_c
      dimension Par(11,NMat),ChPar(NSD*16+4,NMat),ThNew(NumNP),
     !          vNew(NumNP),ThOld(NumNP),vOld(NumNP)

      M=1
      th_c=ChPar(4,M)   ! Water Content from which colloids are excluded
      if(abs(th_c).lt.1.e-20) return
      ChPar(4,M)=0.
      thr=Par(1,M)
      ths=Par(2,M)
      swr=thr/ths
      sw_c=th_c/ths
      alpha=Par(3,M)
      vgn1=Par(4,M)
      xL=Par(6,M)
      vgm1=1.0-1.0/vgn1
      vgn2=vgn1+1
      vgm2=1.0-2.0/vgn2
      if(sw_c.lt.0.or.sw_c.gt.1.) then
        write(*,*) 'Problem with size exclusion!'
        write(*,*) 'Press Enter to continue'
        read(*,*)
      end if

*     Accessible Water Content to Colloid
      ThC=ThNew(1)-ths*sw_c

      sw=ThNew(1)/ths
      sw_eff=(sw-swr)/(1.0-swr)
      sw_c_eff=(sw_c-swr)/(1.0-swr)

*     Colloid Permeability according to Burdine Model
      if(sw_eff.gt.sw_c_eff) then
        kr_c=(sw_eff**2.0)*(((1.0-sw_c_eff**(1.0/vgm2))**vgm2)-
     !                      ((1.0-sw_eff  **(1.0/vgm2))**vgm2)) 
      else
        kr_c=0.0
      end if

*     Mualem Water Relative Permeability
      if(sw_eff.gt.0.0) then
      	krw=(sw_eff**xL)*(1-(1-sw_eff**(1/vgm1))**vgm1)**2.0
        VelC=vNew(1)*kr_c/krw
      else
        krw=0.0
        VelC=vNew(1)
      end if

*     Convert sw_c to r_c
      Pc_c=(1.0/alpha)*(sw_c_eff**(-1./vgm1)-1)**(1.0/vgn1) 
      r_c=(2.0*72.0)/(Pc_c*981.0)
      PorVelSolute=vNew(1)/ThNew(1)
      PorVelColloid=VelC/ThC
c     write(*,*) "r_c microns", 10000.*r_c

      do 11 i=1,NumNP
        ThNew(i)=ThC
        ThOld(i)=ThC
        vNew(i)=VelC
        vOld(i)=VelC
11    continue
      return
      end

************************************************************************

*     Reads parameters for the function expressing reaction rate dependence 
*     on the water content

      subroutine MoistDepIn(cDataPath,cFileName,NMat,NMatD,NS,NSD,
     !                      DMoist,iMoistDep)

      character cFileName*260,cDataPath*260
      dimension DMoist(NMatD,NSD,13,6)

      iLengthPath = Len_Trim(cDataPath)
      cFileName = cDataPath(1:iLengthPath)//'\MoistDep.in'
      open(15,file=cFileName, status='unknown',err=901)

      read(15,*,err=902)
      do 13 M=1,NMat
        read(15,*,err=902)
        do 12 jS=1,NS
          read(15,*,err=902)
          do 11 jReact=1,13
            read(15,*,err=902) (DMoist(M,jS,jReact,i),i=1,6)
11        continue
12      continue
13    continue
      close(15)
      return

*     Error opening an input file 
901   if(iMoistDep.eq.2) iMoistDep=0
      return
*     Error reading from an input file 
902   if(iMoistDep.eq.2) iMoistDep=0
      write(*,*) 'Error reading from an input file MoistDep.in !!!!'
      close(15)
      return
      end

************************************************************************

      real function rMD(NMatD,NSD,M,jS,jReact,DMoist,iReact,WDep,Theta,
     !                  iMoistDep)

*     Function expressing reaction rate dependence on the water content
*     ReacMin0  - relative minimum rate of reaction at low water contents
*     Theta0    - water content at which reaction rate start increasing
*     Theta1    - water content at which reaction rate stops increasing
*     Theta2    - water content at which reaction rate start decreasing  
*     Theta3    - water content at which reaction rate stops decreasing 
*     ReacMin1  - relative minimum rate of reaction at high water contents
*     If theta2=theta3=thetaS -> Anaerobic process
*     If theta0=theta1=0      -> Aerobic process
*     If theta2=0 -> no reduction

      dimension DMoist(NMatD,NSD,13,6),WDep(2+NMatD,NSD*9)

      rMD=1.
      if(iMoistDep.eq.2) then
        if(jReact.eq.0) return
        ReacMin0=DMoist(M,jS,jReact,1)
        Theta0  =DMoist(M,jS,jReact,2)
        Theta1  =DMoist(M,jS,jReact,3)
        Theta2  =DMoist(M,jS,jReact,4)
        Theta3  =DMoist(M,jS,jReact,5)
        ReacMin1=DMoist(M,jS,jReact,6)
        if(abs(Theta2).lt.0.001) return
        if     (Theta.le.Theta0) then
          rMD=ReacMin0
        else if(Theta.le.Theta1) then
          rMD=ReacMin0+(Theta-Theta0)/(Theta1-Theta0)*(1.-ReacMin0)
        else if(Theta.le.Theta2) then
          rMD=1.
        else if(Theta.le.Theta3) then
          rMD=ReacMin1+(Theta-Theta3)/(Theta2-Theta3)*(1.-ReacMin1)
        else
          rMD=ReacMin1
        end if
      else if(iMoistDep.eq.1) then ! Walker's formula
        jjj=(jS-1)*9
        if(WDep(2+M,jjj+iReact).gt.Theta.and.WDep(2+M,jjj+iReact).gt.0.) 
     !    rMD=(Theta/WDep(2+M,jjj+iReact))**WDep(1,jjj+iReact)
      end if

      return
      end

************************************************************************

      real function rMD1(NMatD,NSD,M,jS,jReact,DMoist,Theta)

*     Function expressing reaction rate dependence on the water content
*     Type  =1 or -1: Rate increases or decreases with water content, respectively)
*     Theta0 - water content at which reaction rate start increasing or decreasing  
*     Theta1 - water content at which reaction rate stops increasing or decreasing 
*     ReacMin	- relative minimum rate of reaction

      dimension DMoist(NMatD,NSD,9,4)

      rType  =DMoist(M,jS,jReact,1)
      Theta0 =DMoist(M,jS,jReact,2)
      Theta1 =DMoist(M,jS,jReact,3)
      ReacMin=DMoist(M,jS,jReact,4)
      rMD1=1.
      if(abs(rType).lt.0.1) return
      if(rType.gt.0) then                 ! increasing rate
        if(Theta.ge.Theta1) then
          rMD1=1.
        else if(Theta.le.Theta0) then
          rMD1=ReacMin
        else
          if(abs(Theta1-Theta0).gt.0.)
     !      rMD1=ReacMin+(Theta-Theta0)/(Theta1-Theta0)*(1.-ReacMin)
        end if
      else                                ! decreasing rate
        if(Theta.ge.Theta1) then
          rMD1=ReacMin
        else if(Theta.le.Theta0) then
          rMD1=1.
        else
          if(abs(Theta1-Theta0).gt.0.)
     !      rMD1=ReacMin+(Theta-Theta1)/(Theta0-Theta1)*(1.-ReacMin)
        end if
      end if

      return
      end

*************************************************************************

*     Distribute mass into different phases

      subroutine MassInit(NumNP,NSD,NS,NMat,MatNum,TDep,Temp,ChPar,
     !                    Conc,Theta,ThetaIm,ThSat,lLinear,lBact)

      logical lLinear(NS),lBact
      dimension ChPar(NSD*16+4,NMat),MatNum(NumNP),TDep(NSD*16+4),
     !          Conc(NSD,NumNP),Temp(NumNP),Theta(NumNP),ThSat(NMat),
     !          ThetaIm(NumNP),Par(10)

      Tr=293.15
      R=8.314
      do 12 jS=1,NS
        jjj=(jS-1)*16
        do 11 i=1,NumNP
          if(Conc(jS,i).gt.0.) then
            M=MatNum(i)
            TT=(Temp(i)+273.15-Tr)/R/(Temp(i)+273.15)/Tr
            Par(1)=ChPar(1,     M)*exp(TDep(1)     *TT) !ro
            Par(2)=ChPar(3,     M)*exp(TDep(3)     *TT) !frac
            Par(3)=ChPar(jjj+ 7,M)*exp(TDep(jjj+ 7)*TT) !xKs
            Par(4)=ChPar(jjj+ 8,M)*exp(TDep(jjj+ 8)*TT) !xNu
            Par(5)=ChPar(jjj+ 9,M) !*exp(TDep(jjj+ 9)*TT) !fExp
            Par(6)=ChPar(jjj+10,M)*exp(TDep(jjj+ 7)*TT) !xKH
            Par(7)=Theta(i)+ThetaIm(i)
            Par(8)=ThSat(M)-Par(7)
            if(lBact) then
              rKa1=ChPar(jjj+19,M)*exp(TDep(jjj+19)*TT)
              rKd1=ChPar(jjj+20,M)*exp(TDep(jjj+20)*TT)
              xKs1=Theta(i)*rKa1/Par(1)/rKd1 
              rKa2=ChPar(jjj+16,M)*exp(TDep(jjj+16)*TT)
              rKd2=ChPar(jjj+17,M)*exp(TDep(jjj+17)*TT)
              xKs2=Theta(i)*rKa1/Par(1)/rKd1
              Par(3)=xKs1+xKs2 
            end if
            if(lLinear(jS)) then
              Conc(jS,i)=Conc(jS,i)/(Par(7)+Par(1)*Par(3)+Par(8)*Par(6))
            else
              Conc(jS,i)=cInit(Conc(jS,i),Par,NPar)
            end if
          end if
11      continue
12    continue

      return
      end


************************************************************************

*     Evaluate Liquid concentration from the total solute mass

      real function cInit(xMass,Par,NPar)

      dimension Par(NPar)

      x1=1.e-3
      x2=1.e+3
      call ZBRAK1(X1,X2,XB1,XB2,xMass,Par,NPar)
      cInit=ZBRENT1(XB1,XB2,xMass,Par,NPar)

      return
      end

************************************************************************

      function SolMass(Conc,xMass,Par,NPar)

*     Calculate total solute mass for concentration Conc

      dimension Par(NPar)

      ro=Par(1)
      frac=Par(2)
      xKs=Par(3)
      xNu=Par(4)
      fExp=Par(5)
      xKH=Par(6)
      Theta=Par(7)
      ThetaA=Par(8)

      yMass=Theta*Conc+ro*xKs*Conc**fExp/(1.+xNu*Conc**fExp)+
     !      ThetaA*xKH*Conc
      SolMass=xMass-yMass

      return
      end

************************************************************************

*     Bracketing of the root, Numerical recepies (345)

      subroutine ZBRAK1(X1,X2,XB1,XB2,xMass,Par,NPar)

      dimension Par(NPar)

      NBB=1
      NB=1000

      dlh=(alog10(X2)-alog10(X1))/(NB-1)
      FP=SolMass(X1,xMass,Par,NPar)
      do 11 i=1,NB
        dx2=alog10(X1)+(i)*dlh
        X2=10**dx2
        FC=SolMass(X2,xMass,Par,NPar)
        if(FC*FP.lt.0.) then
          NBB=NBB+1
          XB1=X1
          XB2=X2
          return
        end if
        FP=FC
        X1=X2
        if(NBB.eq.NB) return
11    continue

      return
      end

************************************************************************

*     Brent method of finding root that lies between x1 and x2, 
*     Numerical recepies (354)

      function ZBRENT1(X1,X2,xMass,Par,NPar)

      parameter (ITMAX=100,EPS=3.E-8,TOL=1.e-6)
      dimension Par(NPar)

      A=X1
      B=X2
      FA=SolMass(A,xMass,Par,NPar)
      FB=SolMass(B,xMass,Par,NPar)
      IF(FB*FA.GT.0.) PAUSE 'Root must be bracketed for ZBRENT1.'
      FC=FB
      DO 11 ITER=1,ITMAX
        IF(FB*FC.GT.0.) THEN
          C=A
          FC=FA
          D=B-A
          E=D
        ENDIF
        IF(ABS(FC).LT.ABS(FB)) THEN
          A=B
          B=C
          C=A
          FA=FB
          FB=FC
          FC=FA
        ENDIF
        TOL1=2.*EPS*ABS(B)+0.5*TOL
        XM=.5*(C-B)
        IF(ABS(XM).LE.TOL1 .OR. FB.EQ.0.)THEN
          ZBRENT1=B
          RETURN
        ENDIF
        IF(ABS(E).GE.TOL1 .AND. ABS(FA).GT.ABS(FB)) THEN
          S=FB/FA
          IF(A.EQ.C) THEN
            P=2.*XM*S
            Q=1.-S
          ELSE
            Q=FA/FC
            R=FB/FC
            P=S*(2.*XM*Q*(Q-R)-(B-A)*(R-1.))
            Q=(Q-1.)*(R-1.)*(S-1.)
          ENDIF
          IF(P.GT.0.) Q=-Q
          P=ABS(P)
          IF(2.*P .LT. MIN(3.*XM*Q-ABS(TOL1*Q),ABS(E*Q))) THEN
            E=D
            D=P/Q
          ELSE
            D=XM
            E=D
          ENDIF
        ELSE
          D=XM
          E=D
        ENDIF
        A=B
        FA=FB
        IF(ABS(D) .GT. TOL1) THEN
          B=B+D
        ELSE
          B=B+SIGN(TOL1,XM)
        ENDIF
        FB=SolMass(B,xMass,Par,NPar)
11    CONTINUE
      PAUSE 'ZBRENT1 exceeding maximum iterations.'
      ZBRENT1=B

      RETURN
      END

* ||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||