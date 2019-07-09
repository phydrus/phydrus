* Source file OUTPUT.FOR |||||||||||||||||||||||||||||||||||||||||||||||

      subroutine TLInf(N,Con,x,CosAlf,t,dt,IterW,IterC,TLevel,rTop,
     !                 rRoot,vRoot,hNew,hRoot,CumQ,ItCum,KodTop,KodBot,
     !                 ConvgF,lWat,lChem,cRoot,NS,NSD,Conc,cvTop,cvBot,
     !                 cvCh0,cvCh1,Peclet,Courant,wCumT,wCumA,cCumT,
     !                 cCumA,CumCh,ThNew,ThOld,Sink,lScreen,ierr,cvChR,
     !                 cvChIm,lPrint,lVapor,lWTDep,ConLT,ConVh,ConVT,
     !                 Temp,rSoil,Prec,nPrStep,ThVOld,ThVNew,xConv,
     !                 iDualPor,SinkIm,WTransf,lDensity,SnowLayer,
     !                 lCentrif,Radius,ThNewIm,WLayer,hCritS,lEnd,
     !                 lFluxOut,jPrint,vTop,vBot,lFlux,NObs,Node,vNew,
     !                 cNew,cTop)
      integer TLevel
      double precision t
      logical ConvgF,lWat,lChem,lScreen,lPrint,lVapor,lWTDep,
     !        lFluxOut,lDensity,lCentrif,WLayer,lEnd,lFlux
      dimension CumQ(12),cRoot(NS),Conc(NSD,N),cvTop(NS),cvBot(NS),
     !          cvCh0(NS),cvCh1(NS),cCumT(NS),cCumA(NS),CumCh(10,NS),
     !          cvChR(NS),cvChIm(NS),Con(N),ThOld(N),ThNew(N),hNew(N),
     !          Sink(N),x(N),ConLT(N),ConVh(N),ConVT(N),Temp(N),
     !          ThVNew(N),ThVOld(N),SinkIm(N),ThNewIm(N),Node(NObs),
     !          vNew(N),cNew(N),cTop(NS),cRunOff(15),cGWL(15)

      fRE=1.
      Grav=CosAlf
      M=N-1
      if(lDensity) fRE=fRo(1,Conc(1,N))
      if(lCentrif) Grav=CosAlf*(Radius+abs((x(N)+x(M))/2.))
      dxN=x(N)-x(M)
      vT=-(Con(N)+Con(M))/2.*((hNew(N)-hNew(M))/dxN+Grav*fRE)-
     !    (ThNew(N)-ThOld(N))*fRE*dxN/2./dt-Sink(N)*dxN/2.
      if(lDensity) fRE=fRo(1,Conc(1,1))
      dx1=x(2)-x(1)
      if(lCentrif) Grav=CosAlf*(Radius+abs((x(2)+x(1))/2.))
      vB=-(Con(1)+Con(2))/2.*((hNew(2)-hNew(1))/dx1+Grav*fRE)+
     !  (ThNew(1)-ThOld(1))*fRE*dx1/2./dt+Sink(1)*dx1/2.
      if(iDualPor.gt.0) vT=vT-SinkIm(N)*dxN/2.
      if(lWTDep) then
        vT=vT-(ConLT(N)+ConLT(M))/2.*(Temp(N)-Temp(M))/dxN
        vB=vB-(ConLT(1)+ConLT(2))/2.*(Temp(2)-Temp(1))/dx1
      end if
      vTopW=vT
      if(lVapor) then
        vT=vT-(ConVh(N)+ConVh(M))/2.*(hNew(N)-hNew(M))/dxN
        vB=vB-(ConVh(1)+ConVh(2))/2.*(hNew(2)-hNew(1))/dx1
        vT=vT-(ConVT(N)+ConVT(M))/2.*(Temp(N)-Temp(M))/dxN-
     !        (ThVNew(N)-ThVOld(N))*dxN/2./dt
        vB=vB-(ConVT(1)+ConVT(2))/2.*(Temp(2)-Temp(1))/dx1+
     !        (ThVNew(1)-ThVOld(1))*dx1/2./dt
      end if
      vTopV=vT-vTopW
      vTop=vT
      vBot=vB
      vRunOff=0.

      if((.not.WLayer.or.(WLayer.and.hNew(N).ge.hCritS)).and.rTop.lt.0.)
     !                     vRunOff=abs(rTop-vTop)
      if(vRunOff.lt.1.e-5) vRunOff=0.
      rInfil=0.
      rEvap=0.
      if(vTop.lt.0..and.(Prec.gt.0.or.(WLayer.and.hNew(N).gt.0.))) 
     !                             rInfil=-vTop+rSoil
      if(vTop.ge.0..and.Prec.gt.0) rInfil=Prec
      if(vTop.gt.0.)               rEvap=vTop+Prec
      if(vTop.le.0..and.rSoil.gt.0.and.Prec.gt.0) rEvap=rSoil
      if(vTop.lt.0..and.WLayer.and.hNew(N).gt.0.) rEvap=rSoil
	if(vTop.lt.0..and.KodTop.gt.0) rInfil=-vTop
      CumQ(1)=CumQ(1)+rTop *dt
      CumQ(2)=CumQ(2)+rRoot*dt
      CumQ(3)=CumQ(3)+vTop *dt
      CumQ(4)=CumQ(4)+vRoot*dt
      CumQ(5)=CumQ(5)+vBot *dt
      CumQ(6)=CumQ(6)+vRunOff*dt
      CumQ(7)=CumQ(7)+rInfil*dt
      CumQ(8)=CumQ(8)+rEvap*dt
      if(lFluxOut) then
        CumQ(9) =CumQ(9) +Prec *dt
        CumQ(10)=CumQ(10)+rSoil*dt
      end if
      CumQ(11)=CumQ(11)+WTransf*dt
      wCumT=wCumT+(vBot-vTop-vRoot)*dt
      wCumA=wCumA+(abs(vBot)+abs(vTop)+abs(vRoot))*dt
      if(lChem) then
        do 11 jS=1,NS
          CumCh(1,jS)=CumCh(1,jS)-cvTop(jS)*dt
          CumCh(2,jS)=CumCh(2,jS)+cvBot(jS)*dt
          CumCh(3,jS)=CumCh(3,jS)+cvCh0(jS)*dt
          CumCh(4,jS)=CumCh(4,jS)+cvCh1(jS)*dt
          CumCh(5,jS)=CumCh(5,jS)+cvChR(jS)*dt
          CumCh(6,jS)=CumCh(6,jS)+cvChIm(jS)*dt
          cCumT(jS)=cCumT(jS)+(cvTop(jS)-cvBot(jS)-cvCh0(jS)-cvCh1(jS)+
     !                         cvChR(jS))*dt
          cCumA(jS)=cCumA(jS)+(abs(cvBot(jS))+abs(cvTop(jS))+
     !                         abs(cvCh0(jS))+abs(cvCh1(jS))+
     !                         abs(cvChR(jS)))*dt
          if(lFlux) then
            if(jS.eq.1) then ! using flux concentration (available)
              if(NObs.ge.1) 
     !          CumCh(7,jS)=CumCh(7,jS)+vNew(Node(1))*cNew(Node(1))*dt
              if(NObs.ge.2) 
     !          CumCh(8,jS)=CumCh(8,jS)+vNew(Node(2))*cNew(Node(2))*dt
              if(NObs.ge.3) 
     !          CumCh(9,jS)=CumCh(9,jS)+vNew(Node(3))*cNew(Node(3))*dt
            else             ! using resident concentration (flux conc unavailable)
              if(NObs.ge.1) 
     !         CumCh(7,jS)=CumCh(7,jS)+vNew(Node(1))*Conc(jS,Node(1))*dt
              if(NObs.ge.2) 
     !         CumCh(8,jS)=CumCh(8,jS)+vNew(Node(2))*Conc(jS,Node(2))*dt
              if(NObs.ge.3) 
     !         CumCh(9,jS)=CumCh(9,jS)+vNew(Node(3))*Conc(jS,Node(3))*dt
            end if
          end if
          cRunOff(jS)=vRunOff*cTop(jS)
          CumCh(10,jS)=CumCh(10,jS)+cRunOff(jS)*dt

*         Average GWL concentration
          cGWL(jS)=0.
          dGWL=0
          iBreak=0
          do 14 i=1,N-1
            j=i+1
            dx=x(j)-x(i)
            if(hNew(j).gt.0..and.iBreak.eq.0) then
              cGWL(jS)=cGWL(jS)+(Conc(jS,i)+Conc(jS,j))/2.*dx
              dGWL=dGWL+dx          
            else
              iBreak=1
            end if
14        continue
          if(dGWL.gt.0.) cGWL(jS)=cGWL(jS)/dGWL

11      continue

      end if

      Volume=0.
      do 10 i=N-1,1,-1
        j=i+1
        dx=x(j)-x(i)
        VNewi=dx*(ThNew(i)+ThNew(j))/2.
        if(lDensity) VNewi=dx*(ThNew(i)*fRo(1,(Conc(1,i)))+
     !                         ThNew(j)*fRo(1,(Conc(1,j))))/2.
        Volume=Volume+VNewi
        if(iDualPor.gt.0.) then
          VNewImi=dx*(ThNewIm(i)+ThNewIm(j))/2.
c          if(lDensity) VNewImi=dx*(ThNewIm(i)*fRo(1,(Sorb(1,i)))+
c     !                             ThNewIm(j)*fRo(1,(Sorb(1,j))))/2.
          Volume=Volume+VNewImi
        end if
10    continue

      if(lScreen) then
        if(jPrint.eq.1.and.abs(float((TLevel+20*nPrStep-1)/20/nPrStep)-
     !            (TLevel+20*nPrStep-1)/float(20*nPrStep)).lt.0.0001)
     !  write(*,110)
        if(jPrint.eq.1) then
          if(t.lt.9999999.) then
            if(xConv.eq.1.) then
              write(*,122) t,IterW,ItCum,vTop,CumQ(3),CumQ(4),CumQ(5),
     !                     hNew(N),hRoot,hNew(1)
            else
              write(*,120) t,IterW,ItCum,vTop,CumQ(3),CumQ(4),CumQ(5),
     !                     hNew(N),hRoot,hNew(1)
            end if
          else
            if(xConv.eq.1.) then
              write(*,123) t,IterW,ItCum,vTop,CumQ(3),CumQ(4),CumQ(5),
     !                     hNew(N),hRoot,hNew(1)
            else
              write(*,121) t,IterW,ItCum,vTop,CumQ(3),CumQ(4),CumQ(5),
     !                     hNew(N),hRoot,hNew(1)
            end if
          end if
        end if
      end if
      if(TLevel.eq.1.and.lPrint) then
        write(71,130,err=901)
        if(lFluxOut) write(44,111)
        if(lChem) then
          write(70,150,err=902)
          do 12 jS=1,NS
            write(80+jS,160,err=903)
12        continue
        else
          write(70,140,err=902)
        end if
      end if
      if(lPrint.and.jPrint.eq.1) then
        if(lWat.or.TLevel.eq.1.or.lEnd) then
          if(t.lt.9999999.) then
c            write(71,170,err=901) t,rTop,vTopW,vTop,vTopV,vBot,(CumQ(i),
            write(71,170,err=901) t,rTop,rRoot,vTop,vRoot,vBot,(CumQ(i),
     !               i=1,5),hNew(N),hRoot,hNew(1),vRunOff,CumQ(6),
     !               Volume,CumQ(7),CumQ(8),TLevel,CumQ(11),SnowLayer
	    else
            write(71,171,err=901) t,rTop,rRoot,vTop,vRoot,vBot,(CumQ(i),
     !               i=1,5),hNew(N),hRoot,hNew(1),vRunOff,CumQ(6),
     !               Volume,CumQ(7),CumQ(8),TLevel,CumQ(11),SnowLayer
	    end if
	  end if
        if(lChem) then
          write(70,180,err=902) TLevel,t,dt,IterW,IterC,ItCum,KodTop,
     !                          KodBot,ConvgF,Peclet,Courant
          do 13 jS=1,NS
            mObs=NObs
            if(.not.lFlux) mObs=0
            if(t.lt.99999999.) then
              write(80+jS,190,err=903) t,-cvTop(jS),cvBot(jS),
     !                  CumCh(1,jS),CumCh(2,jS),CumCh(3,jS),CumCh(4,jS),
     !                  Conc(jS,N),cRoot(jS),Conc(jS,1),cvChR(jS),
     !                  CumCh(5,jS),CumCh(6,jS),TLevel,cGWL(jS),
     !                  cRunOff(jS),CumCh(10,jS),
     !                  (vNew(Node(j))*Conc(jS,Node(j)),
     !                  CumCh(6+j,jS),j=1,min(mObs,3))
            else
              write(80+jS,191,err=903) t,-cvTop(jS),cvBot(jS),
     !                  CumCh(1,jS),CumCh(2,jS),CumCh(3,jS),CumCh(4,jS),
     !                  Conc(jS,N),cRoot(jS),Conc(jS,1),cvChR(jS),
     !                  CumCh(5,jS),CumCh(6,jS),TLevel,cGWL(jS),
     !                  cRunOff(jS),CumCh(10,jS),
     !                  (vNew(Node(j))*Conc(jS,Node(j)),
     !                  CumCh(6+j,jS),j=1,min(mObs,3))
            end if
13        continue
        else
          write(70,200,err=902) TLevel,t,dt,IterW,ItCum,KodTop,KodBot,
     !                          ConvgF
        end if
      end if
      if(lPrint.and.lFluxOut.and.jPrint.eq.1) then
        write(44,222) t,Prec,rSoil,rTop,vTop,rInfil,rEvap,rRoot,vRoot,
     !                vBot,vRunOff,CumQ(9),CumQ(10),CumQ(1),CumQ(3),
     !                CumQ(7),CumQ(8),CumQ(2),CumQ(4),CumQ(5),CumQ(6),
     !                Volume,SnowLayer,vTopW,vTopV
222     format(f14.4,24e13.5)
      end if
      return

*     Error when writing into an output file 
901   ierr=1
      return
902   ierr=2
      return
903   ierr=3
      return

110   format(/
     !'         Time ItW   ItCum  vTop    SvTop    SvRoot   SvBot   ',
     !' hTop hRoot hBot'/)
111   format(/'          time     PrecipP      EvaporP     FluxTopP     
     !FluxTopA     InfiltrA      EvaporA      TranspP      TranspA      
     !FluxBot       RunOff     Sum(PrecP)   Sum(EvapP)  Sum(FlTopP)  Sum
     !(FlTopA)   Sum(InfA)   Sum(EvapA)  Sum(TransP)  Sum(TransA)  Sum(F
     !luxBot)  Sum(RunOff)    Storage   SnowLayer      vTopW        vTop
     !V'/)
120   format(f13.4,i3,i7,4e9.2,f8.1,2f6.0)
121   format(e14.7,i3,i7,4e9.2,f7.0,2f6.0)
122   format(f14.4,i3,i7,4e9.2,f7.2,2f5.2)
123   format(e14.7,i3,i7,4e9.2,f7.2,2f5.2)
c120   format('+',f12.3,2i3,i6,4e9.2,3f6.0)  ! writing at one line
130   format(/
     !'       Time          rTop        rRoot        vTop         vRoot 
     !       vBot       sum(rTop)   sum(rRoot)    sum(vTop)   sum(vRoot)
     !    sum(vBot)      hTop         hRoot        hBot        RunOff   
     ! sum(RunOff)     Volume     sum(Infil)    sum(Evap) TLevel Cum(WTr
     !ans)  SnowLayer'/
     !'        [T]         [L/T]        [L/T]        [L/T]        [L/T] 
     !       [L/T]         [L]          [L]          [L]         [L]    
     !       [L]         [L]           [L]         [L]          [L/T]   
     !      [L]          [L]          [L]          [L]'/)
140   format(//'    TLevel      Time          dt      Iter    ItCum  K',
     !'odT  KodB  Convergency'/)
150   format(//'    TLevel      Time           dt     ItrW ItrC    ItC',
     !'um  KodT  KodB Converg  Peclet   Courant'/)
160   format(' All solute fluxes and cumulative solute fluxes are positi
     !ve into the region'//
     !'       Time         cvTop        cvBot      Sum(cvTop)   Sum(cvBo
     !t)     cvCh0        cvCh1         cTop        cRoot         cBot  
     !      cvRoot    Sum(cvRoot)  Sum(cvNEql) TLevel      cGWL        c
     !RunOff   Sum(cRunOff)    (cv(i),    Sum(cv(i)), i=1,NObs)'/
     !'        [T]        [M/L2/T]     [M/L2/T]      [M/L2]       [M/L2]
     !       [M/L2]      [M/L2]        [M/L3]      [M/L3]        [M/L3] 
     !     [M/L2/T]      [M/L2]       [M/L2]              [M/L3]        
     ![M/L2]      [M/L3]      [M/L2/T]      [M/L2]')
170   format(f13.4,11e13.5,2e13.5,5e13.5,i7,e13.5,f11.3)
171   format(e14.8,11e13.5,2e13.5,5e13.5,i7,e13.5,f11.3)
180   format(i9,e15.7,e13.5,2i5,i9,2i6,l6,2f10.3)
190   format(f14.4,12e13.5,i8,e13.5,8e13.5)
191   format(e15.8,12e13.5,i8,e13.5,8e13.5)
200   format(i9,e15.7,e13.5,i5,i9,2i6,l6,2f10.3)
      end

************************************************************************

      subroutine ALInf(t,CumQ,hNewN,hRoot,hNew1,ALevel,ierr)

      integer ALevel
      double precision t
      dimension CumQ(12)

      if(ALevel.eq.1) write(72,110,err=901)
      if(t.lt.999999.) then
        write(72,120,err=901) t,(CumQ(i),i=1,5),hNewN,hRoot,hNew1,ALevel
      else
        write(72,121,err=901) t,(CumQ(i),i=1,5),hNewN,hRoot,hNew1,ALevel
      end if
      return

*     Error when writing into an output file 
901   ierr=1
      return

110   format(//
     !'   Time         sum(rTop)     sum(rRoot)    sum(vTop)     sum(vRo
     !ot)     sum(vBot)    hTop       hRoot      hBot      A-level'/
     !'    [T]           [L]           [L]           [L]           [L]  
     !          [L]        [L]         [L]       [L] '/)
120   format(f12.5,5e14.6,3f11.3,i8)
121   format(e14.7,5e14.6,3f11.3,i8)
      end

************************************************************************

      subroutine SubReg(N,NMat,NLay,hNew,ThN,ThO,x,MatNum,LayNum,t,dt,
     !                  CosAlf,Con,lChem,Conc,ChPar,PLevel,ths,wCumT,
     !                  wCumA,cCumT,cCumA,wVolI,cVolI,WatIn,SolIn,lWat,
     !                  lTemp,Temp,TPar,TDep,NS,NSD,Sorb,lLinear,lEquil,
     !                  lMobIm,ierr,SubVol,Area,lPrint,lBact,Sorb2,
     !                  lVapor,ThVOld,ThVNew,lWTDep,ConLT,ConVh,ConVT,
     !                  iDualPor,ThNewIm,ThOldIm,lDensity,lCentrif,
     !                  Radius,lDualNEq,cPrev)

      logical lWat,lChem,lTemp,lLinear(NS),lEquil,lPrint,lMobIm(NMat),
     !        lVapor,lWTDep,lBact,lCentrif,lDensity,lDualNEq
      integer PLevel
      double precision t
      dimension hNew(N),ThN(N),ThO(N),x(N),MatNum(N),LayNum(N),
     !          Conc(NSD,N),ChPar(NSD*16+4,NMat),ths(NMat),cCumA(NS),
     !          cCumT(NS),cVolI(NS),WatIn(N),SolIn(N),Temp(N),
     !          TPar(10,NMat),TDep(NSD*16+4),Sorb(NSD,N),Con(N),
     !          ConLT(N),ConVh(N),ConVT(N),Sorb2(NSD,N),ThNewIm(N),
     !          ThOldIm(N),ThVOld(N),ThVNew(N),hMean(10),
     !          cMean(11,10),TMean(10),SubVol(10),SubCha(10),ConVol(11),
     !          ConSub(11,10),cTot(11),SubT(10),Area(10),ConVolIm(11),
     !          ConSubIm(11,10),cMeanIm(11,10),cTotIm(11),ConVolIm2(11),
     !          ConSubIm2(11,10),cPrev(N)

      fRE=1.
      Grav=CosAlf
      ATot=0.
      Tr=293.15
      R=8.314
      if(lWat.or.PLevel.eq.0) then
        Volume=0.
        VolumeIm=0.
        Change=0.
        hTot=0.
        DeltW=0.
      end if
      if(lTemp) then
        TTot=0.
        TVol=0.
      end if
      if(lChem) then
        if(NS.gt.11) then
          write(*,*) 'Dimensions in the Subreg subroutine need to be inc
     !reased !'
          stop
        end if
        do 11 jS=1,NS
          cTot(jS)=0.
          ConVol(jS)=0.
          if(.not.lEquil) ConVolIm(jS)=0.
          if(.not.lEquil) cTotIm(jS)=0.
          if(lBact.or.lDualNEq) ConVolIm2(jS)=0.
11      continue
        DeltC=0.
      end if
      do 13 Lay=1,NLay
        Area(Lay)=0.
        if(lWat.or.PLevel.eq.0) then
          SubVol(Lay)=0.
          SubCha(Lay)=0.
          hMean(Lay)=0.
        end if
        if(lTemp) then
          SubT(Lay)=0.
          TMean(Lay)=0.
        end if
        if(lChem) then
          do 12 jS=1,NS
            ConSub(jS,Lay)=0.
            cMean(jS,Lay)=0.
            if(.not.lEquil)       ConSubIm(jS,Lay)=0.
            if(lBact.or.lDualNEq) ConSubIm2(jS,Lay)=0.
            if(.not.lEquil)       cMeanIm(jS,Lay)=0.
12        continue
        end if
13    continue

      do 15 i=N-1,1,-1
        j=i+1
        cEl=0.
        Mi=MatNum(i)
        Mj=MatNum(j)
        Lay=LayNum(i)
        dx=x(j)-x(i)
        Area(Lay)=Area(Lay)+dx
        ATot=ATot+dx
        TT=(Temp(i)+Temp(j))/2.+273.15
        if(lWat.or.PLevel.eq.0) then
          hE=(hNew(i)+hNew(j))/2.
          VNewi=dx*(ThN(i)+ThN(j))/2.
          VOldi=dx*(ThO(i)+ThO(j))/2.
          if(lDensity) then
            VNewi=dx*(ThN(i)*fRo(1,(Conc(1,i)))+
     !                ThN(j)*fRo(1,(Conc(1,j))))/2.
            VOldi=dx*(ThO(i)*fRo(1,(cPrev(i)))+
     !                ThO(j)*fRo(1,(cPrev(j))))/2.
          end if
          if(lVapor) then
            VNewi=VNewi+dx*(ThVNew(i)+ThVNew(j))/2.
            VOldi=VOldi+dx*(ThVOld(i)+ThVOld(j))/2.
          end if
          Volume=Volume+VNewi
          Change=Change+(VNewi-VOldi)/dt
          SubCha(Lay)=SubCha(Lay)+(VNewi-VOldi)/dt
          SubVol(Lay)=SubVol(Lay)+VNewi
          hTot=hTot+hE*dx
          hMean(Lay)=hMean(Lay)+hE*dx
          if(iDualPor.gt.0) then
            VNewImi=dx*(ThNewIm(i)+ThNewIm(j))/2.
            VOldImi=dx*(ThOldIm(i)+ThOldIm(j))/2.
            if(lDensity) then
              VNewImi=dx*(ThNewIm(i)*fRo(1,(Sorb(1,i)))+
     !                    ThNewIm(j)*fRo(1,(Sorb(1,j))))/2.
              VOldImi=dx*(ThOldIm(i)*fRo(1,(Sorb(1,i)))+
     !                    ThOldIm(j)*fRo(1,(Sorb(1,j))))/2.
            end if
            VolumeIm=VolumeIm+dx*(ThNewIm(i)+ThNewIm(j))/2.
            Change=Change+(VNewImi-VOldImi)/dt
            SubCha(Lay)=SubCha(Lay)+(VNewImi-VOldImi)/dt
            SubVol(Lay)=SubVol(Lay)+VNewImi
          end if
        end if
        if(lTemp) then
          TE=(Temp(i)+Temp(j))/2.
          TNewE=dx*((Temp(i)+273.15)*(TPar(1,Mi)*TPar(7,Mi)+
     !                    TPar(2,Mi)*TPar(8,Mi)+TPar(9,Mi)*ThN(i))+
     !              (Temp(j)+273.15)*(TPar(1,Mj)*TPar(7,Mj)+
     !                    TPar(2,Mj)*TPar(8,Mj)+TPar(9,Mj)*ThN(j)))/2.
          TVol=TVol+TNewE
          SubT(Lay)=SubT(Lay)+TNewE
          TTot=TTot+TE*dx
          TMean(Lay)=TMean(Lay)+TE*dx
        end if
        if(lChem) then
          do 14 jS=1,NS
            jjj=(jS-1)*16
            cE=(Conc(jS,i)+Conc(jS,j))/2.
            TTi=(Temp(i)+273.15-Tr)/R/TT/Tr
            xKsi  =ChPar(jjj+ 7,Mi)*exp(TDep(jjj+ 7)*TTi)
            xNui  =ChPar(jjj+ 8,Mi)*exp(TDep(jjj+ 8)*TTi)
            fExpi =ChPar(jjj+ 9,Mi)*exp(TDep(jjj+ 9)*TTi)
            Henryi=ChPar(jjj+10,Mi)*exp(TDep(jjj+10)*TTi)
            TTj=(Temp(j)+273.15-Tr)/R/TT/Tr
            xKsj  =ChPar(jjj+ 7,Mj)*exp(TDep(jjj+ 7)*TTj)
            xNuj  =ChPar(jjj+ 8,Mj)*exp(TDep(jjj+ 8)*TTj)
            fExpj =ChPar(jjj+ 9,Mj)*exp(TDep(jjj+ 9)*TTj)
            Henryj=ChPar(jjj+10,Mj)*exp(TDep(jjj+10)*TTj)
            C1=1.
            C2=1.
            if(.not.lLinear(jS)) then
              if(Conc(jS,i).gt.0.) C1=Conc(jS,i)**(fExpi-1.)/
     !                       (1.+xNui*Conc(jS,i)**fExpi)
              if(Conc(jS,j).gt.0.) C2=Conc(jS,j)**(fExpj-1.)/
     !                       (1.+xNuj*Conc(jS,j)**fExpj)
            end if
            ThWi=ThN(i)
            ThWj=ThN(j)
            ThImobi=ChPar(4,Mi)
            ThImobj=ChPar(4,Mj)
            ThGi=amax1(0.,ths(Mi)-ThWi)
            ThGj=amax1(0.,ths(Mj)-ThWj)
            if(iDualPor.gt.0) then
              ThImobi=ThNewIm(i)
              ThImobj=ThNewIm(j)
c              ThGi=amax1(0.,ths(Mi)-ThWi+thSIm(Mi)-ThImobi)
c              ThGj=amax1(0.,ths(Mj)-ThWj+thSIm(Mj)-ThImobj)
            end if
            if(lMobIm(Mi).and.iDualPor.eq.0.or.lBact)
     !                                      ThWi=max(ThWi-ThImobi,0.001)
            if(lMobIm(Mj).and.iDualPor.eq.0.or.lBact)
     !                                      ThWj=max(ThWj-ThImobj,0.001)
            f_em=1.
            if(lDualNEq) f_em=ChPar(jjj+13,Mi)
            cNewi=dx/2.*(Conc(jS,i)*
     !          (thWi+f_em*ChPar(3,Mi)*ChPar(1,Mi)*xKsi*C1+thGi*Henryi)+
     !                   Conc(jS,j)*
     !          (thWj+f_em*ChPar(3,Mj)*ChPar(1,Mj)*xKsj*C2+thGj*Henryj))
            ConVol(jS)=ConVol(jS)+cNewi
            ConSub(jS,Lay)=ConSub(jS,Lay)+cNewi
            if(.not.lEquil) then
              if(lMobIm(Mi).or.iDualPor.gt.0) then ! mobile-immobile model
                S1=1.
                S2=1.
                if(.not.lLinear(jS)) then
                  if(Sorb(jS,i).gt.0.) S1=Sorb(jS,i)**(fExpi-1.)/
     !                                 (1.+xNui*Sorb(jS,i)**fExpi)
                  if(Sorb(jS,j).gt.0.) S2=Sorb(jS,j)**(fExpj-1.)/
     !                                 (1.+xNuj*Sorb(jS,j)**fExpj)
	              end if
                cNewiIm=dx/2.*(Sorb(jS,i)*
     !               (ThImobi+(1.-ChPar(3,Mi))*ChPar(1,Mi)*xKsi*S1)+
     !                         Sorb(jS,j)*
     !               (ThImobj+(1.-ChPar(3,Mj))*ChPar(1,Mj)*xKsj*S2))
                if(lDualNEq) then
                  cNewiIm2=dx/2.*(ChPar(1,Mi)*Sorb2(jS,i)+
     !                            ChPar(1,Mj)*Sorb2(jS,j))
                  ConVolIm2(jS)    =ConVolIm2(jS)    +cNewiIm2
                  ConSubIm2(jS,Lay)=ConSubIm2(jS,Lay)+cNewiIm2
                end if
              else                                 ! two-site sorption model
                cNewiIm=dx/2.*(ChPar(1,Mi)*Sorb(jS,i)+
     !                         ChPar(1,Mj)*Sorb(jS,j))
                if(lBact) then
                  cNewiIm2=dx/2.*(ChPar(1,Mi)*Sorb2(jS,i)+
     !                            ChPar(1,Mj)*Sorb2(jS,j))
                  ConVolIm2(jS)    =ConVolIm2(jS)    +cNewiIm2
                  ConSubIm2(jS,Lay)=ConSubIm2(jS,Lay)+cNewiIm2
                end if
              end if
              ConVolIm(jS)    =ConVolIm(jS)    +cNewiIm
              ConSubIm(jS,Lay)=ConSubIm(jS,Lay)+cNewiIm
              cEIm=(Sorb(jS,i)+Sorb(jS,j))/2.
              cMeanIm(jS,Lay)=cMeanIm(jS,Lay)+cEIm*dx
              cTotIm(jS)=cTotIm(jS)+cEIm*dx
            end if
            cTot(jS)=cTot(jS)+cE*dx
            cMean(jS,Lay)=cMean(jS,Lay)+cE*dx
            if(jS.eq.1) cEl=cNewi
14        continue
        end if
        if(PLevel.eq.0) then
          if(lWat)  WatIn(i)=vNewi
          if(lWat.and.iDualPor.gt.0) WatIn(i)=vNewi+vNewImi
          if(lChem) SolIn(i)=cEl
        else
          if(lWat) then
            if(iDualPor.le.0) DeltW=DeltW+abs(WatIn(i)-vNewi)
            if(iDualPor.gt.0) DeltW=DeltW+abs(WatIn(i)-vNewi-vNewImi)
          end if
          if(lChem) DeltC=DeltC+abs(SolIn(i)-cEl)
        end if
15    continue
      do 17 Lay=1,NLay
        if(Area(Lay).gt.0.)then
          if(lWat.or.PLevel.eq.0) hMean(Lay)=hMean(Lay)/Area(Lay)
          if(lTemp)               TMean(Lay)=TMean(Lay)/Area(Lay)
          do 16 jS=1,NS
            if(lChem) then
                              cMean  (jS,Lay)=cMean  (jS,Lay)/Area(Lay)
              if(.not.lEquil) cMeanIm(jS,Lay)=cMeanIm(jS,Lay)/Area(Lay)
            end if
16        continue
        end if
17    continue
      if(ATot.gt.0.) then
        if(lWat.or.PLevel.eq.0) hTot=hTot/ATot
        if(lTemp)               TTot=TTot/ATot
      end if
      do 18 jS=1,NS
        if(lChem.and.ATot.gt.0.) then
          cTot(jS)=cTot(jS)/ATot
          if(.not.lEquil) cTotIm(jS)=cTotIm(jS)/ATot
        end if
18    continue
      if(lDensity) fRE=fRo(1,Conc(1,1))
      if(lCentrif) Grav=CosAlf*(Radius+abs((x(2)+x(1))/2.))
      dx1=x(2)-x(1)
      v1=-(Con(1)+Con(2))/2.*((hNew(2)-hNew(1))/dx1+Grav*fRE)
      if(lDensity) fRE=fRo(1,Conc(1,N))
      if(lCentrif) Grav=CosAlf*(Radius+abs((x(N)+x(N-1))/2.))
      dxN=x(N)-x(N-1)
      vN=-(Con(N)+Con(N-1))/2*((hNew(N)-hNew(N-1))/dxN+Grav*fRE)
      if(lWTDep) then
        v1=v1-(ConLT(1)+ConLT(2)  )/2.*(Temp(2)-Temp(1)  )/dx1
        vN=vN-(ConLT(N)+ConLT(N-1))/2.*(Temp(N)-Temp(N-1))/dxN
      end if
      if(lVapor) then
        v1=v1-(ConVh(1)+ConVh(2)  )/2.*(hNew(2)-hNew(1)  )/dx1
        vN=vN-(ConVh(N)+ConVh(N-1))/2.*(hNew(N)-hNew(N-1))/dxN
        v1=v1-(ConVT(1)+ConVT(2)  )/2.*(Temp(2)-Temp(1)  )/dx1
        vN=vN-(ConVT(N)+ConVT(N-1))/2.*(Temp(N)-Temp(N-1))/dxN
      end if

      if(lPrint) then
        if(t.lt.99999999.) then
          write(76,110,err=901) t
        else
          write(76,111,err=901) t
        end if
        write(76,120,err=901) (i,i=1,NLay)
        write(76,130,err=901)
        write(76,140,err=901)   ATot,  (Area(i),i=1,NLay)
        if(lWat.or.PLevel.eq.0) then
          write(76,150,err=901) Volume,(SubVol(i),i=1,NLay)
          if(iDualPor.gt.0) write(76,151,err=901) VolumeIm
          write(76,160,err=901) Change,(SubCha(i),i=1,NLay)
          write(76,170,err=901) hTot,  ( hMean(i),i=1,NLay)
        end if
        if(lTemp) then
          write(76,180,err=901) TVol,  (  SubT(i),i=1,NLay)
          write(76,190,err=901) TTot,  ( TMean(i),i=1,NLay)
        end if
        if(lChem) then
          do 19 jS=1,NS
            write(76,200,err=901) jS,ConVol(jS),(ConSub(jS,i),i=1,NLay)
            write(76,210,err=901) jS,cTot(jS),  ( cMean(jS,i),i=1,NLay)
            if(.not.lEquil) then
              if(lMobIm(1).or.iDualPor.gt.0) then
                write(76,201,err=901) 
     !                       jS,ConVolIm(jS),(ConSubIm(jS,i),i=1,NLay)
                write(76,211,err=901) 
     !                       jS,cTotIm(jS)  ,(cMeanIm(jS,i) ,i=1,NLay)
              else
                write(76,202,err=901) 
     !                       jS,ConVolIm(jS),(ConSubIm(jS,i),i=1,NLay)
                write(76,212,err=901) 
     !                       jS,cTotIm(jS)  ,(cMeanIm(jS,i) ,i=1,NLay)
              end if
            end if
            if(lBact.or.lDualNEq) write(76,203,err=901) 
     !                       jS,ConVolIm2(jS),(ConSubIm2(jS,i),i=1,NLay)
19        continue
        end if
        if(lWat.or.PLevel.eq.0) write(76,220,err=901) vN,v1
      end if

*     Mass balance calculation
      if(PLevel.eq.0) then
        wVolI=Volume
        if(iDualPor.gt.0) wVolI=Volume+VolumeIm
        if(lChem) then
          do 20 jS=1,NS
            cVolI(jS)=ConVol(jS)
            if(.not.lEquil)     cVolI(jS)=cVolI(jS)+ConVolIm(jS)
            if(lBact.or.lDualNEq) 
     !                          cVolI(jS)=cVolI(jS)+ConVolIm2(jS)
20        continue
        end if
      else
        if(lWat) then
          wBalT=Volume-wVolI-wCumT
          if(iDualPor.gt.0) wBalT=Volume-wVolI-wCumT+VolumeIm
          if(lPrint) write(76,230,err=901) wBalT
          ww=amax1(DeltW,wCumA)
          if(ww.gt.1.e-25) then
            wBalR=abs(wBalT)/ww*100.
            if(lPrint) write(76,240,err=901) wBalR
          end if
        end if
        if(lChem) then
          do 21 jS=1,NS
            cBalT=ConVol(jS)-cVolI(jS)+cCumT(jS)
            if(.not.lEquil)       cBalT=cBalT+ConVolIm (jS)
            if(lBact.or.lDualNEq) cBalT=cBalT+ConVolIm2(jS)
            if(lPrint) write(76,250,err=901) jS,cBalT
            cc=amax1(DeltC,cCumA(jS))
            if(cc.gt.1.e-25) then
              cBalR=abs(cBalT)/cc*100.
              if(lPrint) write(76,260,err=901) jS,cBalR
            end if
21        continue
        end if
      end if
      if(lPrint) write(76,130,err=901)
      return

*     Error when writing into an output file 
901   ierr=1
      return

110   format(
     !    /'----------------------------------------------------------'/
     !        ' Time       [T]',f14.4/
     !     '----------------------------------------------------------')
111   format(
     !    /'----------------------------------------------------------'/
     !        ' Time       [T]',e15.8/
     !     '----------------------------------------------------------')
120   format( ' Sub-region num.               ',9(I7,6x))
130   format( 
     !     '----------------------------------------------------------')
140   format( ' Area     [L]      ',e13.5,9e13.5)
150   format( ' W-volume [L]      ',e13.5,9e13.5)
151   format( ' W-volumeI[L]      ',e13.5,9e13.5)
160   format( ' In-flow  [L/T]    ',e13.5,9e13.5)
170   format( ' h Mean   [L]      ',e13.5,9e13.5)
180   format( ' HeatVol  [M/T2]   ',e13.5,10e13.5)
190   format( ' tMean    [K]      ',f13.3,10f13.3)
200   format( ' ConcVol  [M/L2] ',i1,1x,e13.5,10e13.5)
201   format( ' ConcVolIm[M/L2] ',i1,1x,e13.5,10e13.5)
202   format( ' SorbVolIm[M/L2] ',i1,1x,e13.5,10e13.5)
203   format( ' SorbVolIm2[M/L2]',i1,1x,e13.5,10e13.5)
210   format( ' cMean    [M/L3] ',i1,1x,e13.5,10e13.5)
211   format( ' cMeanIm  [M/L3] ',i1,1x,e13.5,10e13.5)
212   format( ' sMeanIm  [-]    ',i1,1x,e13.5,10e13.5)
220   format( ' Top Flux [L/T]    ',e13.5/
     !        ' Bot Flux [L/T]    ',e13.5)
230   format( ' WatBalT  [L]      ',e13.5)
240   format( ' WatBalR  [%]      ',f13.3)
250   format( ' CncBalT  [M]    ',i1,1x,e13.5)
260   format( ' CncBalR  [%]    ',i1,1x,f13.3)
      end

***********************************************************************

      subroutine NodOut(N,NMat,hNew,thN,Con,x,xSurf,CosAlf,TPrint,
     !                  MatNum,Cap,Bxz,Sink,ConS,NS,NSD,Conc,Temp,Sorb,
     !                  Kappa,lBact,Sorb2,lVapor,lWTDep,ConLT,ConVT,
     !                  ConVh,ThOldT,dt,iDualPor,ThNewIm,SinkIm,STrans,
     !                  lDensity,lCentrif,Radius,lVaporOut,lDualNEq,
     !                  ierr)

      dimension hNew(N),thN(N),Con(N),x(N),MatNum(N),Cap(N),Bxz(N),
     !          Sink(N),ConS(NMat),Conc(NSD,N),Temp(N),Sorb(NSD,N),
     !          Kappa(N),Sorb2(NSD,N),ConLT(N),ConVT(N),ConVh(N),
     !          ThNewIm(N),SinkIm(N),STrans(N)
      logical lBact,lVapor,lWTDep,lDensity,lCentrif,
     !        lDualNEq,lVaporOut
      double precision TPrint

      fRE=1.
      Grav=CosAlf
      xS=xSurf+Radius
      if(TPrint.lt.99999999.) then
        write(75,110,err=901) TPrint
        if(lVaporOut.and.lVapor) write(45,110,err=901) TPrint
      else
        write(75,111,err=901) TPrint
        if(lVaporOut.and.lVapor) write(45,111,err=901) TPrint
      end if
      if(NS.eq.0) then
        if(iDualPor.eq.0) then
          write(75,112,err=901)
        else
          write(75,113,err=901)
        end if
      else
        if(iDualPor.eq.0) then
          write(75,114,err=901)
        else
          write(75,115,err=901)
        end if
      end if
      if(lVaporOut.and.lVapor) write(45,116,err=901)
      do 11 i=N,1,-1
        Mat=MatNum(i)
        if(i.eq.1) then
          if(lDensity) fRE=fRo(1,Conc(1,1))
          if(lCentrif) Grav=CosAlf*(Radius+abs((x(2)+x(1))/2.))
          dx=x(2)-x(1)
          vi=-(Con(1)+Con(2))/2.*((hNew(2)-hNew(1))/dx+Grav*fRE)
          if(lWTDep)
     !    vi =vi-(ConLT(1)+ConLT(2))/2.*(Temp(2)-Temp(1))/dx
          vVi=0.
          if(lVapor) then
            vVii=-(ConVh(1)+ConVh(2))/2.*(hNew(2)-hNew(1))/dx
            vVTi=-(ConVT(1)+ConVT(2))/2.*(Temp(2)-Temp(1))/dx
            vVi=vVii+vVTi
          end if
        else if(i.eq.N) then
          ConSN=ConS(Mat)*Bxz(N)
          N1=N-1
          dx=(x(N)-x(N-1))
          if(lDensity) fRE=fRo(1,Conc(1,N))
          if(lCentrif) Grav=CosAlf*(Radius+abs((x(N)+x(N1))/2.))
          vi=-(Con(N)+Con(N1))/2.*((hNew(N)-hNew(N1))/dx+Grav*fRE)-
     !        (ThN(N)-ThOldT)*fRE*dx/2./dt-Sink(N)*dx/2.
          if(lWTDep) vi=vi-(ConLT(N)+ConLT(N1))/2.*(Temp(N)-Temp(N1))/dx
          vVi=0.
          if(lVapor) then
            vVii=-(ConVh(N)+ConVh(N1))/2.*(hNew(N)-hNew(N1))/dx
            vVTi=-(ConVT(N)+ConVT(N1))/2.*(Temp(N)-Temp(N1))/dx
            vVi=vVTi+vVii
          end if
        else
          dxA=x(i+1)-x(i)
          dxB=x(i)-x(i-1)
          if(lDensity) fRE=(fRo(1,Conc(1,i))+fRo(1,Conc(1,i+1)))/2.
          if(lCentrif) Grav=CosAlf*(Radius+abs((x(i+1)+x(i))/2.))
          vA=-(Con(i)+Con(i+1))/2.*((hNew(i+1)-hNew(i))/dxA+Grav*fRE)
          if(lDensity) fRE=(fRo(1,Conc(1,i))+fRo(1,Conc(1,i-1)))/2.
          if(lCentrif) Grav=CosAlf*(Radius+abs((x(i)+x(i-1))/2.))
          vB=-(Con(i)+Con(i-1))/2.*((hNew(i)-hNew(i-1))/dxB+Grav*fRE)
          vi= (vA*dxA+vB*dxB)/(dxA+dxB)
          if(lWTDep) then
            vAT =-(ConLT(i)+ConLT(i+1))/2.*(Temp(i+1)-Temp(i))/dxA
            vBT =-(ConLT(i)+ConLT(i-1))/2.*(Temp(i)-Temp(i-1))/dxB
            vi  =vi+(vAT*dxA+vBT*dxB)/(dxA+dxB)
          end if
          vVi=0.
          if(lVapor) then
            vVA =-(ConVh(i)+ConVh(i+1))/2.*(hNew(i+1)-hNew(i))/dxA
            vVB =-(ConVh(i)+ConVh(i-1))/2.*(hNew(i)-hNew(i-1))/dxB
            vVii= (vVA*dxA+vVB*dxB)/(dxA+dxB)
            vVAT=-(ConVT(i)+ConVT(i+1))/2.*(Temp(i+1)-Temp(i))/dxA
            vVBT=-(ConVT(i)+ConVT(i-1))/2.*(Temp(i)-Temp(i-1))/dxB
            vVTi= (vVAT*dxA+vVBT*dxB)/(dxA+dxB)
            vVi =vVii+vVTi
          end if
        end if
        if(lVaporOut.and.lVapor) write(45,140,err=901) N-i+1,x(i)-xS,
     !         Con(i),ConLT(i),ConVh(i),ConVT(i),vi,vVi,vi+vVi,vVii,vVTi
        if(hNew(i).gt.-9.9e+05) then
          if(.not.lBact) then
            if(iDualPor.eq.0) then
              if(.not.lDualNEq) then
                write(75,120,err=901) N-i+1,x(i)-xS,hNew(i),thN(i),
     !                Con(i),Cap(i),vi+vVi,Sink(i),Kappa(i),vi/ConSN,
     !                Temp(i),(Conc(jS,i),jS=1,NS),(Sorb(jS,i),jS=1,NS)
              else
                write(75,120,err=901) N-i+1,x(i)-xS,hNew(i),thN(i),
     !                Con(i),Cap(i),vi+vVi,Sink(i),Kappa(i),vi/ConSN,
     !                Temp(i),(Conc(jS,i),jS=1,NS),(Sorb(jS,i),jS=1,NS),
     !                (Sorb2(jS,i),jS=1,NS)
              end if
            else
              if(.not.lDualNEq) then
                write(75,120,err=901) N-i+1,x(i)-xS,hNew(i),thN(i),
     !              SinkIm(i),ThNewIm(i),vi+vVi,STrans(i),Kappa(i),
     !              vi/ConSN,Temp(i),(Conc(jS,i),jS=1,NS),
     !              (Sorb(jS,i),jS=1,NS)
              else
                write(75,120,err=901) N-i+1,x(i)-xS,hNew(i),thN(i),
     !              SinkIm(i),ThNewIm(i),vi+vVi,STrans(i),Kappa(i),
     !              vi/ConSN,Temp(i),(Conc(jS,i),jS=1,NS),
     !              (Sorb(jS,i),jS=1,NS),(Sorb2(jS,i),jS=1,NS)
              end if
            end if
          else
            write(75,120,err=901) N-i+1,x(i)-xS,hNew(i),thN(i),
     !            Con(i),Cap(i),vi+vVi,Sink(i),Kappa(i),vi/ConSN,
     !            Temp(i),(Conc(jS,i),jS=1,NS),(Sorb(jS,i),jS=1,NS),
     !            (Sorb2(jS,i),jS=1,NS)
          end if
        else
          if(.not.lBact) then
            if(iDualPor.eq.0) then
              if(.not.lDualNEq) then
                write(75,130,err=901) N-i+1,x(i)-xS,hNew(i),thN(i),
     !                Con(i),Cap(i),vi+vVi,Sink(i),Kappa(i),vi/ConSN,
     !                Temp(i),(Conc(jS,i),jS=1,NS),(Sorb(jS,i),jS=1,NS)
              else
                write(75,130,err=901) N-i+1,x(i)-xS,hNew(i),thN(i),
     !                Con(i),Cap(i),vi+vVi,Sink(i),Kappa(i),vi/ConSN,
     !                Temp(i),(Conc(jS,i),jS=1,NS),(Sorb(jS,i),jS=1,NS),
     !                (Sorb2(jS,i),jS=1,NS)
              end if
            else
              if(.not.lDualNEq) then
                write(75,130,err=901) N-i+1,x(i)-xS,hNew(i),thN(i),
     !                SinkIm(i),ThNewIm(i),vi+vVi,STrans(i),Kappa(i),
     !                vi/ConSN,Temp(i),(Conc(jS,i),jS=1,NS),
     !                (Sorb(jS,i),jS=1,NS)
              else
                write(75,130,err=901) N-i+1,x(i)-xS,hNew(i),thN(i),
     !              SinkIm(i),ThNewIm(i),vi+vVi,STrans(i),Kappa(i),
     !              vi/ConSN,Temp(i),(Conc(jS,i),jS=1,NS),
     !              (Sorb(jS,i),jS=1,NS),(Sorb2(jS,i),jS=1,NS)
              end if
            end if
          else
            write(75,130,err=901) N-i+1,x(i)-xS,hNew(i),thN(i),
     !            Con(i),Cap(i),vi+vVi,Sink(i),Kappa(i),vi/ConSN,
     !            Temp(i),(Conc(jS,i),jS=1,NS),(Sorb(jS,i),jS=1,NS),
     !            (Sorb2(jS,i),jS=1,NS)
          end if
        end if
11    continue
      write(75,'(''end'')',err=901)
      if(lVaporOut.and.lVapor) write(45,'(''end'')',err=901)
      return

*     Error when writing into an output file 
901   ierr=1
      return

110   format(//' Time:',f14.4//)
111   format(//' Time:',e15.8//)
112   format(
     !' Node      Depth      Head Moisture       K          C         ',
     !'Flux        Sink         Kappa   v/KsTop   Temp'/
     !'           [L]        [L]    [-]        [L/T]      [1/L]      [',
     !'L/T]        [1/T]         [-]      [-]      [C]'/)
113   format(
     !' Node      Depth      Head Moisture    WTrans    Im.Moist.     ',
     !'Flux       STrans        Kappa   v/KsTop   Temp'/
     !'           [L]        [L]    [-]        [1/T]       [-]       [',
     !'L/T]     [M/L*3/T]        [-]      [-]      [C]'/)
114   format(
     !' Node      Depth      Head Moisture       K          C         ',
     !'Flux        Sink         Kappa   v/KsTop   Temp   Conc(1..NS) Sor
     !b(1...NS)'/
     !'           [L]        [L]    [-]        [L/T]      [1/L]      [',
     !'L/T]        [1/T]         [-]      [-]      [C]      [M/L*3]'/)
115   format(
     !' Node      Depth      Head Moisture    WTrans    Im.Moist.     ',
     !'Flux       STrans        Kappa   v/KsTop   Temp   Conc(1..NS) Sor
     !b(1...NS)'/
     !'           [L]        [L]    [-]        [1/T]       [-]       [',
     !'L/T]     [M/L*3/T]        [-]      [-]      [C]      [M/L*3]'/)
116   format(
     !' Node    Depth        Con        ConLT        ConVh        ConVT 
     !      vLiquid       vVapor       vTotal       vVapIso     vVapTerm
     !'/
     !'          [L]        [L/T]     [L2/K/T]       [L/T]      [L2/K/T]
     !       [L/T]         [L/T]        [L/T]        [L/T]        [L/T]'
     !/)
120   format(i4,1x,f10.4,1x,f11.3,1x,f6.4,1x,4e12.4,i8,1x,e11.3,f8.2,
     !       30e12.4)
130   format(i4,1x,f10.4,1x,e11.4,1x,f6.4,1x,4e12.4,i8,2x,f10.3,f8.2,
     !       30e12.4)
140   format(i4,1x,f10.4,1x,9e13.5)
      end

***********************************************************************

      subroutine ObsNod(t,N,NObs,NS,NSD,Node,Conc,hNew,ThNew,TempN,
     !                  lChem,ThNewIm,vNew,vVNew,lFlux,ierr)

      dimension Node(NObs),Conc(NSD,N),ThNew(N),TempN(N),hNew(N),
     !          ThNewIm(N),ECa(10),vNew(N),vVNew(N),Th(1000)
      double precision t
      logical lChem,lEC,lFlux

      do 10 i=1,NObs
        Th(i)=ThNew(Node(i))+ThNewIm(Node(i))
10    continue
      if(.not.lChem) then
        if(t.lt.99999999.) then
          if(lFlux) then
            write(77,102,err=901) t,(hNew(Node(i)),Th(i),
     !                 vNew(Node(i))+vVNew(Node(i)),i=1,NObs)
          else
            write(77,100,err=901) t,(hNew(Node(i)),Th(i),
     !                               TempN(Node(i)),i=1,NObs)
          end if
        else
          if(lFlux) then
            write(77,103,err=901) t,(hNew(Node(i)),Th(i),
     !                 vNew(Node(i))+vVNew(Node(i)),i=1,NObs)
          else
           write(77,101,err=901) t,(hNew(Node(i)),Th(i),
     !                              TempN(Node(i)),i=1,NObs)
          end if
        end if
      else
        lEC=.false.
        if(lEC) then
          do 11 i=1,NObs
            ECw=(Conc(1,Node(i))/0.008465)**(1./1.073)
            Thw=ThNew(Node(i))
            ECa(i)=1.45*ECw*Thw*Thw+0.102*Thw
11        continue
        end if
        if(t.lt.99999999.) then
        if(lFlux) then
          if(NS.eq.1) write(77,112,err=901) t,(hNew(Node(i)),Th(i),
     !                  vNew(Node(i)),(Conc(j,Node(i)),j=1,NS),i=1,NObs)
          if(NS.eq.2) write(77,122,err=901) t,(hNew(Node(i)),Th(i),
     !                  vNew(Node(i)),(Conc(j,Node(i)),j=1,NS),i=1,NObs)
          if(NS.eq.3) write(77,132,err=901) t,(hNew(Node(i)),Th(i),
     !                  vNew(Node(i)),(Conc(j,Node(i)),j=1,NS),i=1,NObs)
          if(NS.eq.4) write(77,140,err=901) t,(hNew(Node(i)),Th(i),
     !                  vNew(Node(i)),(Conc(j,Node(i)),j=1,NS),i=1,NObs)
          if(NS.eq.5) write(77,150,err=901) t,(hNew(Node(i)),Th(i),
     !                  vNew(Node(i)),(Conc(j,Node(i)),j=1,NS),i=1,NObs)
          if(NS.eq.6) write(77,160,err=901) t,(hNew(Node(i)),Th(i),
     !                  vNew(Node(i)),(Conc(j,Node(i)),j=1,NS),i=1,NObs)
          if(NS.eq.7) write(77,170,err=901) t,(hNew(Node(i)),Th(i),
     !                  vNew(Node(i)),(Conc(j,Node(i)),j=1,NS),i=1,NObs)
          if(NS.eq.8) write(77,180,err=901) t,(hNew(Node(i)),Th(i),
     !                  vNew(Node(i)),(Conc(j,Node(i)),j=1,NS),i=1,NObs)
          if(NS.eq.9) write(77,190,err=901) t,(hNew(Node(i)),Th(i),
     !                  vNew(Node(i)),(Conc(j,Node(i)),j=1,NS),i=1,NObs)
          if(NS.eq.10)write(77,200,err=901) t,(hNew(Node(i)),Th(i),
     !                  vNew(Node(i)),(Conc(j,Node(i)),j=1,NS),i=1,NObs)
          if(NS.eq.11)write(77,210,err=901) t,(hNew(Node(i)),Th(i),
     !                  vNew(Node(i)),(Conc(j,Node(i)),j=1,NS),i=1,NObs)
          if(NS.eq.12)write(77,220,err=901) t,(hNew(Node(i)),Th(i),
     !                  vNew(Node(i)),(Conc(j,Node(i)),j=1,NS),i=1,NObs)
	  else
          if(NS.eq.1) write(77,110,err=901) t,(hNew(Node(i)),Th(i),
     !                 TempN(Node(i)),(Conc(j,Node(i)),j=1,NS),i=1,NObs)
          if(NS.eq.2) write(77,120,err=901) t,(hNew(Node(i)),Th(i),
     !                 TempN(Node(i)),(Conc(j,Node(i)),j=1,NS),i=1,NObs)
          if(NS.eq.3) write(77,130,err=901) t,(hNew(Node(i)),Th(i),
     !                 TempN(Node(i)),(Conc(j,Node(i)),j=1,NS),i=1,NObs)
          if(NS.eq.4) write(77,140,err=901) t,(hNew(Node(i)),Th(i),
     !                 TempN(Node(i)),(Conc(j,Node(i)),j=1,NS),i=1,NObs)
          if(NS.eq.5)write(77,150,err=901) t,(hNew(Node(i)),Th(i),
     !                 TempN(Node(i)),(Conc(j,Node(i)),j=1,NS),i=1,NObs)
          if(NS.eq.6) write(77,160,err=901) t,(hNew(Node(i)),Th(i),
     !                 TempN(Node(i)),(Conc(j,Node(i)),j=1,NS),i=1,NObs)
          if(NS.eq.7) write(77,170,err=901) t,(hNew(Node(i)),Th(i),
     !                 TempN(Node(i)),(Conc(j,Node(i)),j=1,NS),i=1,NObs)
          if(NS.eq.8) write(77,180,err=901) t,(hNew(Node(i)),Th(i),
     !                 TempN(Node(i)),(Conc(j,Node(i)),j=1,NS),i=1,NObs)
          if(NS.eq.9) write(77,190,err=901) t,(hNew(Node(i)),Th(i),
     !                 TempN(Node(i)),(Conc(j,Node(i)),j=1,NS),i=1,NObs)
          if(NS.eq.10)write(77,200,err=901) t,(hNew(Node(i)),Th(i),
     !                 TempN(Node(i)),(Conc(j,Node(i)),j=1,NS),i=1,NObs)
          if(NS.eq.11)write(77,210,err=901) t,(hNew(Node(i)),Th(i),
     !                 TempN(Node(i)),(Conc(j,Node(i)),j=1,NS),i=1,NObs)
          if(NS.eq.12)write(77,220,err=901) t,(hNew(Node(i)),Th(i),
     !                 TempN(Node(i)),(Conc(j,Node(i)),j=1,NS),i=1,NObs)
	  end if
        else
        if(lFlux) then
          if(NS.eq.1) write(77,111,err=901) t,(hNew(Node(i)),Th(i),
     !                  vNew(Node(i)),(Conc(j,Node(i)),j=1,NS),i=1,NObs)
          if(NS.eq.2) write(77,121,err=901) t,(hNew(Node(i)),Th(i),
     !                  vNew(Node(i)),(Conc(j,Node(i)),j=1,NS),i=1,NObs)
          if(NS.eq.3) write(77,131,err=901) t,(hNew(Node(i)),Th(i),
     !                  vNew(Node(i)),(Conc(j,Node(i)),j=1,NS),i=1,NObs)
          if(NS.eq.4) write(77,141,err=901) t,(hNew(Node(i)),Th(i),
     !                  vNew(Node(i)),(Conc(j,Node(i)),j=1,NS),i=1,NObs)
          if(NS.eq.5) write(77,151,err=901) t,(hNew(Node(i)),Th(i),
     !                  vNew(Node(i)),(Conc(j,Node(i)),j=1,NS),i=1,NObs)
          if(NS.eq.6) write(77,161,err=901) t,(hNew(Node(i)),Th(i),
     !                  vNew(Node(i)),(Conc(j,Node(i)),j=1,NS),i=1,NObs)
          if(NS.eq.7) write(77,171,err=901) t,(hNew(Node(i)),Th(i),
     !                  vNew(Node(i)),(Conc(j,Node(i)),j=1,NS),i=1,NObs)
          if(NS.eq.8) write(77,181,err=901) t,(hNew(Node(i)),Th(i),
     !                  vNew(Node(i)),(Conc(j,Node(i)),j=1,NS),i=1,NObs)
          if(NS.eq.9) write(77,191,err=901) t,(hNew(Node(i)),Th(i),
     !                  vNew(Node(i)),(Conc(j,Node(i)),j=1,NS),i=1,NObs)
          if(NS.eq.10)write(77,201,err=901) t,(hNew(Node(i)),Th(i),
     !                  vNew(Node(i)),(Conc(j,Node(i)),j=1,NS),i=1,NObs)
          if(NS.eq.11)write(77,211,err=901) t,(hNew(Node(i)),Th(i),
     !                  vNew(Node(i)),(Conc(j,Node(i)),j=1,NS),i=1,NObs)
          if(NS.eq.12)write(77,221,err=901) t,(hNew(Node(i)),Th(i),
     !                  vNew(Node(i)),(Conc(j,Node(i)),j=1,NS),i=1,NObs)
	  else
          if(NS.eq.1) write(77,111,err=901) t,(hNew(Node(i)),Th(i),
     !                 TempN(Node(i)),(Conc(j,Node(i)),j=1,NS),i=1,NObs)
          if(NS.eq.2) write(77,121,err=901) t,(hNew(Node(i)),Th(i),
     !                 TempN(Node(i)),(Conc(j,Node(i)),j=1,NS),i=1,NObs)
          if(NS.eq.3) write(77,131,err=901) t,(hNew(Node(i)),Th(i),
     !                 TempN(Node(i)),(Conc(j,Node(i)),j=1,NS),i=1,NObs)
          if(NS.eq.4) write(77,141,err=901) t,(hNew(Node(i)),Th(i),
     !                 TempN(Node(i)),(Conc(j,Node(i)),j=1,NS),i=1,NObs)
          if(NS.eq.5) write(77,151,err=901) t,(hNew(Node(i)),Th(i),
     !                 TempN(Node(i)),(Conc(j,Node(i)),j=1,NS),i=1,NObs)
          if(NS.eq.6) write(77,161,err=901) t,(hNew(Node(i)),Th(i),
     !                 TempN(Node(i)),(Conc(j,Node(i)),j=1,NS),i=1,NObs)
          if(NS.eq.7) write(77,171,err=901) t,(hNew(Node(i)),Th(i),
     !                 TempN(Node(i)),(Conc(j,Node(i)),j=1,NS),i=1,NObs)
          if(NS.eq.8) write(77,181,err=901) t,(hNew(Node(i)),Th(i),
     !                 TempN(Node(i)),(Conc(j,Node(i)),j=1,NS),i=1,NObs)
          if(NS.eq.9) write(77,191,err=901) t,(hNew(Node(i)),Th(i),
     !                 TempN(Node(i)),(Conc(j,Node(i)),j=1,NS),i=1,NObs)
          if(NS.eq.10)write(77,201,err=901) t,(hNew(Node(i)),Th(i),
     !                 TempN(Node(i)),(Conc(j,Node(i)),j=1,NS),i=1,NObs)
          if(NS.eq.11)write(77,211,err=901) t,(hNew(Node(i)),Th(i),
     !                 TempN(Node(i)),(Conc(j,Node(i)),j=1,NS),i=1,NObs)
          if(NS.eq.12)write(77,221,err=901) t,(hNew(Node(i)),Th(i),
     !                 TempN(Node(i)),(Conc(j,Node(i)),j=1,NS),i=1,NObs)
	  end if
        end if
      end if

100   format(2x,f14.4,100(f12.2,f8.4,f9.3,       2x))
101   format( x,e15.8,100(f12.2,f8.4,f9.3,       2x))
102   format(2x,f14.4,100(f12.2,f8.4,e11.3,      2x))
103   format( x,e15.8,100(f12.2,f8.4,e11.3,      2x))
110   format(2x,f14.4,100(f12.2,f8.4,f9.3, e12.4,2x))
111   format( x,e15.8,100(f12.2,f8.4,f9.3, e12.4,2x))
112   format(2x,f14.4,100(f12.2,f8.4,     2e12.4,2x))
120   format(2x,f14.4,100(f12.2,f8.4,f9.3,2e12.4,2x))
121   format( x,e15.8,100(f12.2,f8.4,f9.3,2e12.4,2x))
122   format(2x,f14.4,100(f12.2,f8.4,     3e12.4,2x))
130   format(2x,f14.4,100(f12.2,f8.4,f9.3,3e12.4,2x))
131   format( x,e15.8,100(f12.2,f8.4,f9.3,3e12.4,2x))
132   format(2x,f14.4,100(f12.2,f8.4,     4e12.4,2x))
140   format(2x,f14.4,100(f12.2,f8.4,f9.3,4e12.4,2x))
141   format( x,e15.8,100(f12.2,f8.4,f9.3,4e12.4,2x))
150   format(2x,f14.4,100(f12.2,f8.4,f9.3,5e12.4,2x))
151   format( x,e15.8,100(f12.2,f8.4,f9.3,5e12.4,2x))
160   format(2x,f14.4,100(f12.2,f8.4,f9.3,6e12.4,2x))
161   format( x,e15.8,100(f12.2,f8.4,f9.3,6e12.4,2x))
170   format(2x,f14.4,100(f12.2,f8.4,f9.3,7e12.4,2x))
171   format( x,e15.8,100(f12.2,f8.4,f9.3,7e12.4,2x))
180   format(2x,f14.4,100(f12.2,f8.4,f9.3,8e12.4,2x))
181   format( x,e15.8,100(f12.2,f8.4,f9.3,8e12.4,2x))
190   format(2x,f14.4,100(f12.2,f8.4,f9.3,9e12.4,2x))
191   format( x,e15.8,100(f12.2,f8.4,f9.3,9e12.4,2x))
200   format(2x,f14.4,100(f12.2,f8.4,f9.3,10e12.4,2x))
201   format( x,e15.8,100(f12.2,f8.4,f9.3,10e12.4,2x))
210   format(2x,f14.4,100(f12.2,f8.4,f9.3,11e12.4,2x))
211   format( x,e15.8,100(f12.2,f8.4,f9.3,11e12.4,2x))
220   format(2x,f14.4,100(f12.2,f8.4,f9.3,12e12.4,2x))
221   format( x,e15.8,100(f12.2,f8.4,f9.3,12e12.4,2x))
      return

*     Error when writing into an output file 
901   ierr=1
      return

      end

* ||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||