* Source file HYSTER.FOR ||||||||||||||||||||||||||||||||||||||||||||||

*     Written/revised by RJ Lenhard, Nov 2004
*     Program to evaluate hysteretic water-wet K-S-P routines

      subroutine Hyst(NumNP,NMatD,ParD,ParW,MatNum,Kappa,hNew,hOld,
     !                Theta,Con,Cap,IKappa,iKod)

      implicit real*8(A-H,O-Z)
      implicit integer*4 (I-N)
      logical lPrint
      real ParD(11,NMatD),ParW(11,NMatD),hNew(NumNP),hOld(NumNP),
     !     Theta(NumNP),Con(NumNP),Cap(NumNP)
      dimension MatNum(NumNP),Kappa(NumNP)
      common /properties/ ALPHAD,ALPHAI,XN,XM,XXM,SM,SARWI,IPATH,IHYST,
     !                    XNW,XMW,XXMW
      common /glob/ RHSW(1001,7),RASW(1001,7),SARW(1001),RSWAW(1001),
     !              PHSW(1001),MPSW(1001),JJH(1001),IPSW(1001) 
      common /local/ SW,SRAW,SAT,RSW,RAW,PERMW,MPSW1,JJJH,IPSW1,ILSW,
     !               ESW,ESAT,ASW,DWW

*     SW  - actual water content
*     ESW - effective water content
*     ASW - apparent water content
*     PERMW - relative conductivity
*     ALPHAD = 'alpha for main drainage'
*     ALPHAI = 'alpha for main imbibition'
*     XN = 'VG n parameter' for drying
*     XNW= 'VG n parameter' for wetting
*     XM = 1.d+0 - 1.d+0/XN
*     XXM= 1.d+0/XM
*     SM = 'Sr'

      lPrint=.false.
*     iPath - number of reversal points the code remembers (<=7)
*     =1: Fluid entrapment only, non hysteretic, i.e., alphaw=alphad
      iPath=7
*     iHyst - Location of the initial condition
*     =1: Main drainage curve
*     =2: Main inhibition curve
*     =3: Primary drainage curve (with maximum entrapped air)
      iHyst=1
      if(iKod.eq.1) then
        if(IKappa.eq.-1) iHyst=1
        if(IKappa.eq. 1) iHyst=2
        call HysterIni(NumNP)
      end if

      do 11 n=1,NumNP
        m=MatNum(n)
        AlphaD=ParD(3,m)
        AlphaI=ParW(3,m)
        xn=ParD(4,m)
        xm=1.d+0-1.d+0/xn
        xxm=1.d+0/xm
*       different n values for the drying curve
        xnw=ParW(4,m)
        xmw=1.d+0-1.d+0/xnw
        xxmw=1.d+0/xmw
*       Residual water saturation
        sm=ParD(1,m)/ParD(2,m)
c        sm=ParD(1,m)/(ParD(2,m)-ParD(1,m))
*       Max amount of air that get traped on main inhibition branch (saturation)
c        SARWI=(ParD(2,m)-ParW(2,m))/ParD(2,m)
        SARWI=(ParD(2,m)-ParW(2,m))/(ParD(2,m)-ParD(1,m))

*       Initialize arrays for hysteretic saturations
        hW=hNew(n)
        hA=0.
        PHSW(n)=-hOld(n)

*       subprogram for hysteretic routines
        call Path(n,hW,hA)
 
*       after convergence
        if(iKod.eq.1.or.iKod.eq.3) call UpdateHyst(n,hW,hA)
        Theta(n)=SW*ParD(2,m)
c        Theta(n)=SW*(ParD(2,m)-ParD(1,m))+ParD(1,m)
        Con(n)=ParD(5,m)*sngl(PermW)
        Cap(n)=sngl(DWW)*ParD(2,m)
        if(JJH(n).eq.0) Kappa(n)=-1
        if(JJH(n).eq.1) Kappa(n)= 1
      
        if(lPrint.and.iKod.eq.3) then
          write(*,1000)
          write(*,1002)
          write(*,1003) RASW(1,1),RASW(1,2),RASW(1,3),RASW(1,4),
     !                  RASW(1,5),RASW(1,6),RASW(1,7)
          write(*,1004) RHSW(1,1),RHSW(1,2),RHSW(1,3),RHSW(1,4),
     !                  RHSW(1,5),RHSW(1,6),RHSW(1,7)
          write(*,1005) PHSW(1),SARW(1),SRAW,SARWI,RSWAW(1),RSW,JJJH,
     !                  JJH(1),IPSW(1),IPSW1,MPSW(1),MPSW1,RAW
          write(*,1006) ASW,ESW,SW,ILSW,ESAT,SAT,PERMW
          write(*,1007)
        end if
11    continue      
      return

1000  format(//)
1001  format(/,10X,'CURRENT PRESSURE = ',F8.3,7X,'PREVIOUS HAW = ',F8.3)
1002  format(3X,'CONVERGED RESULTS')
1003  format(/5X,'RASW-1 = ',F5.3,5X,'RASW-2 = ',F5.3,5X,'RASW-3 = ',
     !F5.3,5X,'RASW-4 = ',F5.3,/5X,'RASW-5 = ',F5.3,5X,'RASW-6 = ',F5.3,
     !5X,'RASW-7 = ',F5.3)
1004  format(/5X,'RHSW-1 = ',F5.1,5X,'RHSW-2 = ',F5.1,5X,'RHSW-3 = ',
     !F5.1,5X,'RHSW-4 = ',F5.1,/5X,'RHSW-5 = ',F5.1,5X,'RHSW-6 = ',F5.1,
     !5X,'RHSW-7 = ',F5.1)
1005  format(//5X,'PHSW = ',F7.1,5X,'SARW = ',F5.3,5X,'SRAW = ',F5.3,5X,
     !'SARWI = ',F5.3,/5X,'RSWAW = ',F5.3,4X,'RSW = ',F5.3,6X,'JJJH = ',
     !I1,9X,'JJH = ',I1,/5X,'IPSW = ',I1,9X,'IPSW1 = ',I1,8X,'MPSW = ',
     !I1,9X,'MPSW1 = ',I1,/5X,'RAW = ',F5.3)
1006  format(//5X,'ASW = ',F5.3,6X,'ESW = ',F5.3,6X,'SW = ',F5.3,5X,
     !'ILSW = ',I1,/5X'ESAT = ',F5.3,5X,'SAT = ',F5.3,6X,'PERMW = ',
     !E12.3)  
1007  format(3X,'**********************************************')
 
      end
      
************************************************************************

*     Initialize arrays for hysteretic saturations.
*     Written/revised by RJ Lenhard, Nov. 2004

      subroutine HysterIni(N)
      implicit real*8 (A-H,O-Z)
      implicit integer*4 (I-N)
      common /properties/ ALPHAD,ALPHAI,XN,XM,XXM,SM,SARWI,IPATH,IHYST,
     !                    XNW,XMW,XXMW 
      common /glob/ RHSW(1001,7),RASW(1001,7),SARW(1001),RSWAW(1001),
     !              PHSW(1001),MPSW(1001),JJH(1001),IPSW(1001) 
      common /local/ SW,SRAW,SAT,RSW,RAW,PERMW,MPSW1,JJJH,IPSW1,ILSW,
     !               ESW,ESAT,ASW,DWW

*     Initializing reversal points, depending on number of saturation
*     paths being considered (odd versus even number)

c      IRES = restart option
c      IF( IRES .GE. 1) GO TO 210
      ZERO=0.D+0
      ONE =1.D+0
      DO 100, I = 1,IPATH
        IF( IPATH-(I*2) .EQ. 0) GO TO 110
  100 CONTINUE
      GO TO 160
  110 CONTINUE
      DO 130, IJ = 4,IPATH,2
        DO 120, JI = 1,N
          RASW(JI,IJ) = ONE
  120   CONTINUE
  130 CONTINUE
      DO 150, IJ = 3,IPATH-1,2
        DO 140, JI = 1,N
          RASW(JI,IJ) = ZERO
  140   CONTINUE
  150 CONTINUE
      GO TO 210
  160 CONTINUE
      DO 180, IJ = 4,IPATH-1,2
        DO 170, JI = 1,N
          RASW(JI,IJ) = ONE
  170   CONTINUE
  180 CONTINUE
      DO 200, IJ = 3,IPATH,2
        DO 190, JI = 1,N
          RASW(JI,IJ) = ZERO
  190   CONTINUE
  200 CONTINUE
  210 CONTINUE
  
      DO 220, I = 1,N
        RASW(I,1) = ONE
        RHSW(I,1) = ZERO
        MPSW(I) = 3
        SARW(I) = ZERO
        RSWAW(I) = ONE
  220 CONTINUE

*     Setting historical values for given initial conditions

*     For nonhysteresis and fluid entrapment options
      IF( IPATH .EQ. 1 ) THEN
        IF( SARWI .EQ. ZERO ) GO TO 260
        IF( IHYST .EQ. 2) RSWAW(I) = ZERO
        GO TO 260
      ENDIF

*     For hysteresis options         
      IF( IHYST .EQ. 1 ) THEN

*---    Starting from the main drainage branch

        DO 230, I = 1,N
          RHSW(I,2) = ZERO
          RASW(I,2) = ONE
          JJH(I) = 0
          PHSW(I) = ZERO 
          RSWAW(I) = ONE   
          IPSW(I) = 1
  230   CONTINUE
      ELSE IF ( IHYST .EQ. 2 ) THEN

*       Starting from the main imbibition branch

        DO 240, I = 1,N
          RHSW(I,2) = 1.D+3
          RASW(I,2) = ZERO
          PHSW(I) = 1.D+3
          JJH(I) = 1
          RSWAW(I) = ZERO
          IPSW(I) = 2
  240   CONTINUE
      ELSE IF ( IHYST .EQ. 3 ) THEN

*       Starting from drainage scanning curve when all depths where
*       previously apparent saturated

        DO 250, I = 1,N
          RHSW(I,2) = 1.D+5
          RHSW(I,3) = ZERO
          RASW(I,2) = ZERO
          RASW(I,3) = ONE
          PHSW(I) = ZERO
          JJH(I) = 0
          RSWAW(I) = ZERO
          IPSW(I) = 3
  250   CONTINUE
      ENDIF
  260 CONTINUE    

*    End of HYINI group

      RETURN
      END
      
************************************************************************

*     Beginning subprogram for hysteretic routines
*     Written/revised by RJ Lenhard, Nov 2004

      subroutine PATH(N,HW,HA)

      implicit real*8 (A-H,O-Z)
      implicit integer*4 (I-N)
      common /properties/ ALPHAD,ALPHAI,XN,XM,XXM,SM,SARWI,IPATH,IHYST,
     !                    XNW,XMW,XXMW 
      common /glob/ RHSW(1001,7),RASW(1001,7),SARW(1001),RSWAW(1001),
     !              PHSW(1001),MPSW(1001),JJH(1001),IPSW(1001) 
      common /local/ SW,SRAW,SAT,RSW,RAW,PERMW,MPSW1,JJJH,IPSW1,ILSW,
     !               ESW,ESAT,ASW,DWW

      ZERO=0.D+0
      ONE =1.D+0

*     Assigning property values

      IF( SARWI .EQ. ZERO ) THEN
        RAW = ZERO
        SRAW = ZERO
      ELSE  
        RAW = (ONE/SARWI) - ONE
      ENDIF  
      HAW = MAX(ZERO,(HA-HW))

*     Initializing common block variables

      ESAT = ZERO
      ASW = ZERO

*     Setting global variables to local variables

      JJJH = JJH(N)
      IPSW1 = IPSW(N)
      MPSW1 = MPSW(N)
      RSW = RSWAW(N)
      IF( RSW .EQ. 0. ) SARW(N) = SARWI
      SRAW = SARW(N)

*     Two phase K-S-P relations

      CALL HAWPATH(N,HAW)
      RETURN

*     End of PATH group

      END

************************************************************************

*     Determining saturation path history & calling k-S-P routines
*     Written/revised by RJ Lenhard, Nov 2004

      subroutine HAWPATH( N,HAW )
      implicit real*8 (A-H,O-Z)
      implicit integer*4 (I-N)
      common /properties/ ALPHAD,ALPHAI,XN,XM,XXM,SM,SARWI,IPATH,IHYST,
     !                    XNW,XMW,XXMW 
      common /glob/ RHSW(1001,7),RASW(1001,7),SARW(1001),RSWAW(1001),
     !              PHSW(1001),MPSW(1001),JJH(1001),IPSW(1001) 
      common /local/ SW,SRAW,SAT,RSW,RAW,PERMW,MPSW1,JJJH,IPSW1,ILSW,
     !               ESW,ESAT,ASW,DWW

      ZERO=0.D+0
      ONE =1.D+0

*     For fluid entrapment only option

      IF( IPATH .EQ. 1 ) THEN
        CALL DRAIN( N,HAW )
        ILSW=0
        RETURN
      ENDIF

*     For hysteresis option
*     Main air-water drainage branch
    
  100 CONTINUE
      IF( HAW .GE. RHSW(N,2) ) THEN
        CALL DRAIN( N,HAW )
        ILSW=0
        JJJH=0
        IPSW1=1
        MPSW1=3
        RETURN
      ENDIF

*     Wetting scanning paths

  110 CONTINUE
      IF(((((HAW .LT. PHSW(N)) .AND. (JJJH.EQ.0)) .OR. 
     & ((HAW .le. PHSW(N)) .AND. (JJJH.EQ.1))) .AND. (MPSW1.EQ.3)) .OR. 
     & (MPSW1 .EQ. 1) .OR. ((JJJH .EQ. 1) .AND. (HAW .EQ. ZERO))) THEN
     
*       On max. sat. path, but passed reveral point and became a drying path
        IF( (IPSW1 .EQ. IPATH) .AND. (HAW .GT. RHSW(N,IPATH)) ) THEN
          MPSW1=0
          IPSW1=IPATH-1
          GOTO 100
        ENDIF

*       Close of path 2 at HAW=0.
        if(HAW.eq.0.d0) then
          if(IPSW1.lt.iPath) MPSW1=3
          IPSW1=2
          JJJH=1
          call WET( N, HAW)
          return
        end if

*       Switching from drying to wetting paths      
        IF( (JJJH.EQ.0) .AND. (MPSW1.EQ.3) ) THEN
          IPSW1=IPSW1+1
     
*         Switched to max. sat. path, which is a wetting path     
          IF( IPSW1 .GE. IPATH ) THEN
            MPSW1=1
            IPSW1=IPATH
          ENDIF
        ENDIF
        
*       Evaluating the current wetting path number and whether
*       saturation paths have closed  
        IPSW2 = IPSW1
        DO 120, L = 0,IPSW2-1,2
          IF( (HAW .LT. RHSW(N,IPSW1-L)) .AND. 
     &     (HAW .GT. RHSW(N,IPSW1-1-L)) ) THEN
            IPSW1=IPSW1-L
            IF( IPSW1 .LT. IPATH ) MPSW1=3  
            CALL WET( N,HAW )
            ILSW=L          
            JJJH=1
            RETURN
          ENDIF
  120   CONTINUE
        WRITE( *,'(A)') 'ERROR: Passed Through Wetting Scanning Path Loo
     &p Without Entering subroutine WET2P'
        GOTO 200

*     Drying scanning paths

      ELSEIF(((((HAW .GT. PHSW(N)) .AND. (JJJH .EQ. 1)) .OR.
     &  ((HAW .ge. PHSW(N)) .AND. (JJJH .EQ. 0))) .AND. (MPSW1 .EQ. 3))       
     &  .OR. (MPSW1 .EQ. 0)) THEN

*       K-S-P Relations corresponding to the main drainage branch

        IF( IPSW1 .EQ. 1 ) THEN
          CALL DRAIN( N,HAW )
          ILSW=0
          JJJH=0
          IPSW1=1
          MPSW1=3
          RETURN
        ENDIF
        
*       On max. sat. path, but passed reveral point and became 
*       a wetting path        
        IF( (IPSW1 .EQ. IPATH) .AND. (HAW .LT. RHSW(N,IPATH)) ) THEN
          MPSW1=1
          IPSW1=IPATH-1
          GOTO 110
        ENDIF
        
*       Switching from wetting to drying paths        
        IF( (JJJH .EQ. 1) .AND. (MPSW1 .EQ. 3) ) THEN
          IPSW1=IPSW1+1
          IF( IPSW1 .GE. IPATH ) THEN
            MPSW1=0
            IPSW1=IPATH
          ENDIF
        ENDIF
        
*       Evaluating the current drying path number and whether
*       saturation paths have closed     
        IPSW2 = IPSW1
        DO 140, L = 0,IPSW2-2,2
          IF( (HAW .GT. RHSW(N,IPSW1-L)) .AND.
     &      (HAW .LT. RHSW(N,IPSW1-1-L)) )THEN
            IPSW1=IPSW1-L
            IF( IPSW1 .LT. IPATH ) MPSW1=3
            CALL DRY( N,HAW )
            ILSW=L
            JJJH=0
            RETURN
          ENDIF
  140   CONTINUE
        WRITE( *,'(A)') 'ERROR: Passed Through Three-Phase Liquid Wettin
     &g Loop Without Entering subroutine DRY2P'
        GOTO 200
      ELSE
        WRITE( *,'(A)') 'ERROR: Passed Through subroutine HAWPATH Withou
     &t Entering a Two-Phase Liquid Saturation Path'
      ENDIF
  200 CONTINUE

      RETURN 
      END     

************************************************************************

*     Main drainage k-S-P relations 
*     Options: nonhysteretic, fluid entrapment only, hysteretic
*     Written/revised by RJ Lenhard, Nov 2004

      subroutine DRAIN( N,HAW )
      implicit real*8 (A-H,O-Z)
      implicit integer*4 (I-N)
      common /properties/ ALPHAD,ALPHAI,XN,XM,XXM,SM,SARWI,IPATH,IHYST,
     !                    XNW,XMW,XXMW 
      common /glob/ RHSW(1001,7),RASW(1001,7),SARW(1001),RSWAW(1001),
     !              PHSW(1001),MPSW(1001),JJH(1001),IPSW(1001) 
      common /local/ SW,SRAW,SAT,RSW,RAW,PERMW,MPSW1,JJJH,IPSW1,ILSW,
     !               ESW,ESAT,ASW,DWW

*     Apparent water saturation for all options

      ZERO=0.D+0
      ONE =1.D+0
      ASW = MIN(((1.D+0+(ALPHAD*HAW)**XN)**(-XM)),ONE)
      IF( ASW .LT. RSW ) THEN
        RSW = ASW
        IF( SARWI .NE. 0. ) THEN
          SRAW = (ONE-ASW)/(ONE+RAW*(ONE-ASW))
        ELSE
          SRAW = 0.D+0
        ENDIF    
      ENDIF
      IF( IPATH .EQ. 1 ) THEN
        ESAT = SRAW*(ASW-RSW)/(ONE-RSW)
      ELSE
        ESAT= 0.D+0
      ENDIF

*     Effective and actual saturations for all options

      ESW = ASW-ESAT
      SW = ESW*(ONE-SM)+SM
      SAT = ESAT*(ONE-SM)
      IF( IPATH .NE. 1 ) GOTO 100

*     Derivative terms for nonhysteretic/fluid entrapment options

      X1 = ASW**XXM
      DWWD = (ONE-SM)*ALPHAD*(XN-ONE)*X1*(ONE-X1)**XM
      IF( ESAT .EQ. 0. ) THEN
        DWW = DWWD
      ELSE
        CSATW = SRAW/(ONE-RSW)
        X7 = MAX(ZERO,ONE-CSATW)
        DWW = X7*DWWD
      ENDIF

*     Permeability terms for nonhysteretic/fluid entrapment options

      IF( ESAT .EQ. 0. ) THEN
        CSATW = 0.D+0
        P3 = 0.D+0
      ELSE
        CSATW = SRAW/(ONE-RSW)
        P3 = CSATW*(ONE-RSW**XXM)**XM
      ENDIF  
      P0 = ONE-(ONE-CSATW)*(ONE-ASW**XXM)**XM
      PERMW = (ESW**.5)*(P0-P3)**2.
      GOTO 110

*     Derivative terms for hysteresis option

  100 CONTINUE
      X1 = ASW**XXM
      DWW = (ONE-SM)*ALPHAD*(XN-ONE)*X1*(ONE-X1)**XM

*     Permeability terms for hysteresis option

      P0 = ONE-(ONE-ESW**XXM)**XM
      PERMW = (ESW**.5)*(P0**2.)
  110 CONTINUE

*     Relative permeability error message

      IF( (PERMW.LT.0.D+0) .OR. (PERMW.GT.ONE) ) THEN
        WRITE(*,'(A,I6)') 'ERROR: Water Relative Permeability in Subro-
     &utine DRAIN @ Node ',N
        STOP
      ENDIF

*     End of DRAIN group

      RETURN
      END

************************************************************************

*     Hysteretic k-S-P relations of scanning drying path
*     Written/revised by RJ Lenhard, Nov 2004

      subroutine DRY( N,HAW )

      implicit real*8 (A-H,O-Z)
      implicit integer*4 (I-N)
      common /properties/ ALPHAD,ALPHAI,XN,XM,XXM,SM,SARWI,IPATH,IHYST,
     !                    XNW,XMW,XXMW 
      common /glob/ RHSW(1001,7),RASW(1001,7),SARW(1001),RSWAW(1001),
     !              PHSW(1001),MPSW(1001),JJH(1001),IPSW(1001) 
      common /local/ SW,SRAW,SAT,RSW,RAW,PERMW,MPSW1,JJJH,IPSW1,ILSW,
     !               ESW,ESAT,ASW,DWW

*     Saturations of drying scanning curves for hysteresis option

      ZERO=0.D+0
      ONE =1.D+0
      RIDSWD=(ONE+(ALPHAD*RHSW(N,IPSW1))**XN)**(-XM)
      RDISWD=(ONE+(ALPHAD*RHSW(N,IPSW1-1))**XN)**(-XM)
      SWD=(ONE+(ALPHAD*HAW)**XN)**(-XM)
      ASW=(((SWD-RDISWD)*(RASW(N,IPSW1)-RASW(N,IPSW1-1)))/
     &  (RIDSWD-RDISWD))+RASW(N,IPSW1-1)
      ASW=MIN(ASW,ONE)
      ESAT=SRAW*((ASW-RSW)/(ONE-RSW))
      ESW=ASW-ESAT
      SW=ESW*(ONE-SM)+SM
      SAT=ESAT*(ONE-SM)

*     Permeability terms for hysteresis option

      P0=ONE-((ONE-ASW**XXM)**XM)
      CSATW=SRAW/(ONE-RSW)
      P3=((ONE-RSW**XXM)**XM)-((ONE-ASW**XXM)**XM)
      P3=CSATW*P3
      PERMW=(ESW**.5)*((P0-P3)**2.)

*     Derivative terms for hysteresis option

      X1=SWD**XXM
      X7=MAX(ZERO,ONE-CSATW)
      X8=(RASW(N,IPSW1)-RASW(N,IPSW1-1))/(RIDSWD-RDISWD)
      DWW=(ONE-SM)*ALPHAD*(XN-ONE)*X1*X7*X8*(ONE-X1)**XM

*     Relative permeability error message

      IF( (PERMW.LT.0.D+0) .OR. (PERMW.GT.ONE) ) THEN
        WRITE(*,'(A,I6)') 'ERROR: Water Relative Permeability in Subrout-
     &ine DRY @ Node ',N
        STOP
      ENDIF

*     End of DRY group

      RETURN
      END

************************************************************************

*     Hysteretic k-S-P relations of scanning wetting path
*     Written/revised by RJ Lenhard, Nov 2004

      subroutine WET( N,HAW )
      implicit real*8 (A-H,O-Z)
      implicit integer*4 (I-N)
      common /properties/ ALPHAD,ALPHAI,XN,XM,XXM,SM,SARWI,IPATH,IHYST,
     !                    XNW,XMW,XXMW 
      common /glob/ RHSW(1001,7),RASW(1001,7),SARW(1001),RSWAW(1001),
     !              PHSW(1001),MPSW(1001),JJH(1001),IPSW(1001) 
      common /local/ SW,SRAW,SAT,RSW,RAW,PERMW,MPSW1,JJJH,IPSW1,ILSW,
     !               ESW,ESAT,ASW,DWW

*     Saturations of wetting scanning paths for hysteresis option

      ZERO=0.D+0
      ONE =1.D+0
      RIDSWI = (ONE+(ALPHAI*RHSW(N,IPSW1-1))**XNw)**(-XMw)
      RDISWI = (ONE+(ALPHAI*RHSW(N,IPSW1))**XNw)**(-XMw)
      SWI = (ONE+(ALPHAI*HAW)**XNw)**(-XMw)
      ASW = (((SWI-RIDSWI)*(RASW(N,IPSW1)-RASW(N,IPSW1-1)))/
     &  (RDISWI-RIDSWI))+RASW(N,IPSW1-1)
      ASW = MIN(ASW,ONE)
      ESAT = SRAW*((ASW-RSW)/(ONE-RSW))
      ESW = ASW-ESAT
      SW = ESW*(ONE-SM)+SM
      SAT = ESAT*(ONE-SM)

*     Permeability terms for hysteresis option

      P0 = ONE-((ONE-ASW**XXM)**XM)
      CSATW = SRAW/(ONE-RSW)
      P3 = ((ONE-RSW**XXM)**XM)-((ONE-ASW**XXM)**XM)
      P3 = CSATW*P3
      PERMW = (ESW**.5)*((P0-P3)**2.)

*     Derivative terms for hysteresis option

      X1 = SWI**XXMw
      X7 = MAX(ZERO,ONE-CSATW)
      X8 = (RASW(N,IPSW1)-RASW(N,IPSW1-1))/(RDISWI-RIDSWI)
      DWW = (ONE-SM)*ALPHAI*(XNw-ONE)*X1*X8*X7*(ONE-X1)**XMw

*     Relative permeability error message

      IF( PERMW.LT.ZERO .OR. PERMW.GT.ONE ) THEN
        WRITE(*,'(A,I6)') 'ERROR: Water Relative Permeability in Subrout-
     &ine WET @ Node ',N
        STOP
      ENDIF

*     End of WET group

      RETURN
      END

************************************************************************

*     Sets hysteresis variables following convergence
*     Written/revised by RJ Lenhard, Nov 2004

      subroutine UPDATEHYST(N,HW,HA)

      implicit real*8 (A-H,O-Z)
      implicit integer*4 (I-N)
      common /properties/ ALPHAD,ALPHAI,XN,XM,XXM,SM,SARWI,IPATH,IHYST,
     !                    XNW,XMW,XXMW 
      common /glob/ RHSW(1001,7),RASW(1001,7),SARW(1001),RSWAW(1001),
     !              PHSW(1001),MPSW(1001),JJH(1001),IPSW(1001) 
      common /local/ SW,SRAW,SAT,RSW,RAW,PERMW,MPSW1,JJJH,IPSW1,ILSW,
     !               ESW,ESAT,ASW,DWW

      ZERO=0.D+0
      ONE =1.D+0

*       HW = 'converged water pressure for the node'
*       HA = 'converged air pressure for the node'        
        HAW = MAX( ZERO,(HA-HW) )
        IF( ASW .LT. RSWAW(N) ) THEN
          RSWAW(N) = ASW
          IF( SARWI .NE. ZERO ) SARW(N) = (ONE-ESW)/
     &     (ONE+RAW*(ONE-ESW))
        ENDIF
        IF( IPATH.EQ.1 ) GOTO 200

*       Setting reversal points when HAW = 0 on path 2  
        IF( (HAW .EQ. ZERO) .AND. (IPSW1 .NE. 1) ) THEN
          IPSW1 = 2
          JJJH = 1
          DO 102, M = IPSW1+3,IPATH,2
            RASW(N,M) = ZERO
  102     CONTINUE
          DO 103, M = IPSW1+2,IPATH,2
            RASW(N,M) = ONE
  103     CONTINUE
          RHSW(N,IPSW1+1) = HAW
          RASW(N,IPSW1+1) = ASW
          PHSW(N) = HAW
        ENDIF
        
*       A drying scanning path has closed - resetting reversal points 
        IF( (ILSW .GT. 0) .AND. (JJJH .EQ. 0) ) THEN
          DO 120, M = IPSW1+3,IPATH,2
            RASW(N,M) = ONE
  120     CONTINUE
          DO 121, M = IPSW1+2,IPATH,2
            RASW(N,M) = ZERO
  121     CONTINUE
          RHSW(N,IPSW1+1) = HAW
          RASW(N,IPSW1+1) = ASW
          PHSW(N) = HAW
        
*       A wetting scanning path has closed - resetting reversal points 
        ELSEIF( ILSW.GT.0 .AND. JJJH.EQ.1 ) THEN
          DO 122, M = IPSW1+3,IPATH,2
            RASW(N,M) = ZERO
  122     CONTINUE
          DO 123, M = IPSW1+2,IPATH,2
            RASW(N,M) = ONE
  123     CONTINUE
          RHSW(N,IPSW1+1) = HAW
          RASW(N,IPSW1+1) = ASW
          PHSW(N) = HAW
        ENDIF

*       Drying scanning path closed with main drainage - resetting reversals  
        IF( (HAW .GE. RHSW(N,2)) .AND. 
     &    (RASW(N,3) .GT. ZERO) ) THEN
          IPSW1 = 1
          JJJH = 0
          DO 124, M = IPSW1+2,IPATH,2
            RASW(N,M) = ZERO
  124     CONTINUE
          DO 125, M = IPSW1+3,IPATH,2
            RASW(N,M) = ONE
  125     CONTINUE
          RHSW(N,IPSW1+1) = HAW
          RASW(N,IPSW1+1) = ASW
          PHSW(N) = HAW
        ENDIF

*       Setting gobal variables from converged local variables  
        JJH(N) = JJJH
        IPSW(N) = IPSW1
        MPSW(N) = MPSW1
        IF( IPSW1 .EQ. IPATH ) GOTO 200

*       No scanning paths closed, only continued  
        IF( (JJJH.EQ.1) .AND. (ASW .GT. RASW(N,IPSW1+1)) ) THEN
          RHSW(N,IPSW1+1) = HAW
          RASW(N,IPSW1+1) = ASW
          PHSW(N) = HAW
        ENDIF
        IF( (JJJH.EQ.0) .AND. (ASW .LT. RASW(N,IPSW1+1)) ) THEN
          RHSW(N,IPSW1+1) = HAW
          RASW(N,IPSW1+1) = ASW
          PHSW(N) = HAW
        ENDIF
        IF( IPSW1.EQ.1 ) THEN
          PHSW(N) = RHSW(N,2)
        ENDIF
  200 CONTINUE

      RETURN
      END

* ||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||