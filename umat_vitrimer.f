      SUBROUTINE UMAT(STRESS,STATEV,DDSDDE,SSE,SPD,SCD,
     1 RPL,DDSDDT,DRPLDE,DRPLDT,
     2 STRAN,DSTRAN,TIME,DTIME,TEMP,DTEMP,PREDEF,DPRED,CMNAME,
     3 NDI,NSHR,NTENS,NSTATV,PROPS,NPROPS,COORDS,DROT,PNEWDT,
     4 CELENT,DFGRD0,DFGRD1,NOEL,NPT,LAYER,KSPT,KSTEP,KINC)
C
      INCLUDE 'ABA_PARAM.INC'
C
      CHARACTER*80 CMNAME
      DIMENSION STRESS(NTENS),STATEV(NSTATV),
     1 DDSDDE(NTENS,NTENS),DDSDDT(NTENS),DRPLDE(NTENS),
     2 STRAN(NTENS),DSTRAN(NTENS),TIME(2),PREDEF(1),DPRED(1),
     3 PROPS(NPROPS),COORDS(3),DROT(3,3),DFGRD0(3,3),DFGRD1(3,3)
C
      REAL*8 I1(NDI,NDI), Fv_old(NDI,NDI), invFv_old(NDI,NDI),Fv3_old(NDI,NDI),invFv3_old(NDI,NDI),
     +    Fv_new(NDI,NDI), invFv_new(NDI,NDI), Fv3_new(NDI,NDI), invFv3_new(NDI,NDI),
     +    Fe_old(NDI,NDI), Fe_int(NDI,NDI), Fe_new(NDI,NDI), Fe3_old(NDI,NDI), Fe3_int(NDI,NDI), Fe3_new(NDI,NDI),
     +    Fve_old(NDI,NDI), Fve_int(NDI,NDI), Fve_new(NDI,NDI),
     +    Fb_old(NDI,NDI), invFb_old(NDI,NDI), Fcure(NDI,NDI),Fthm(NDI,NDI),FBER(NDI,NDI),
     +    Fb_new(NDI,NDI), invFb_new(NDI,NDI), Fcure_old(NDI,NDI),Fthm_old(NDI,NDI),FBER_old(NDI,NDI),
     +    Fcure_tau(NDI,NDI), Fthm_tau(NDI,NDI), FBER_tau(NDI,NDI),
     +    Falpha_old(NDI,NDI),Falpha_tau(NDI,NDI), Falpha(NDI,NDI),
     +    invFalpha(NDI,NDI), F_eff(NDI,NDI), F_new(NDI,NDI), 
     +    Ee(NDI,NDI), Ee_d(NDI,NDI), Eve(NDI,NDI), Eve_d(NDI,NDI), 
     +    Te(NDI,NDI), Te_d, Tve(NDI,NDI), Tve_d(NDI,NDI), Tve3(NDI,NDI),Tve3_d(NDI,NDI),
     +    T_old(NDI,NDI), T_int(NDI,NDI), T_int_d(NDI,NDI),
     +    T_new(NDI,NDI), Tp(NDI,NDI), Tp_d(NDI,NDI),
     +    Dv(NDI,NDI), Dvp(NDI,NDI), Dv_exp(NDI,NDI),Dvp3(NDI,NDI), Dv3_exp(NDI,NDI),
     +    Db(NDI,NDI), Dbp(NDI,NDI), Db_exp(NDI,NDI),
     +    Qve(NDI,NDI), QveT(NDI,NDI), Qe(NDI,NDI), Qe3(NDI,NDI),Ve3(NDI,NDI),Ee3(NDI,NDI),Ee3_d(NDI,NDI),
     +    Q(NDI,NDI), QT(NDI,NDI), SIGP(NDI), Vve(NDI,NDI), Ve(NDI,NDI),SIGP3(NDI),Qve3(NDI,NDI),QveT3(NDI,NDI),
     +    D1(NTENS,NTENS), D2(NTENS,NTENS),D3(NTENS,NTENS)
C
      REAL*8 gdotv, gdotb, sve, sb, kappa1, kappa2,kappa3, mu1,mu2,mu3,gdotv3,sve3,
     +    tau_T0, tau_B0, tau_T, tau_B, Tf_T0, Tf_B0, Tf_T, Tf_B, Tf3_T,Tf3_T0
     +    Tg, Tv, C1T, C2T, alpha_NT, Ea, beta_B, alpha_NB,trTve3,
     +    phidot, alpha_C, phi0, E1, E2, phi, phi_old, dphi, deps_M,deps_M0,
     +    alpha_T, alpha_B, TEMP0, deps_C, deps_T, deps_B,deps_Tg,deps_Tg0,deps_T0,deps_Tv,deps_Tv0,
     +    trTve, trT_int, trEe,trEe3, trEve,Tcure,theta_tau, alpha_NT0, AF, temp_ref,E3, tau3_T0,tau3_T,tau2_T,tau2_T0

C
      REAL*8 ZERO,ONE,TWO,THREE,HALF,THIRD,R,nu
      PARAMETER (ZERO = 0.D0, ONE=1.D0, TWO=2.D0, THREE = 3.D0, 
     +    THIRD=1.D0/3.D0, TEN = 10.D0, R=8.3145, nu=0.49D0)
      !----------------------------------------------------------------------
      ! ELASTIC PARAMETERS
      !----------------------------------------------------------------------  
c     kappa = PROPS(1)       ! nominal bulk parameter
      !E1 = PROPS(2)         ! nominal storage modulus of elastic leg
      !E2 = PROPS(3)         ! nominal storage modulus of viscoelastic leg1	  	  
	  theta_tau = TEMP + DTEMP	  

	  CALL COMPUTEMU_E1(TEMP,E1)
      CALL COMPUTEMU_E2(TEMP,E2)
      CALL COMPUTEMU_E3(TEMP,E3) 	  

	  CALL COMPUTEMU_TAU2(TEMP,tau2_T)
	  CALL COMPUTEMU_TAU3(TEMP,tau3_T)
      !----------------------------------------------------------------------
      ! RELAXATION PARAMETERS
      !----------------------------------------------------------------------
      tau_T0 = PROPS(4)      ! thermal relaxation timescale
      tau_B0 = PROPS(5)      ! BER relaxation timescale 
      !----------------------------------------------------------------------
      ! THERMAL EXPANSION PARAMETERS
      !----------------------------------------------------------------------
      C1T = PROPS(6)         ! WLF parameter 1
      C2T = PROPS(7)         ! WLF parameter 2
      alpha_NT0 = PROPS(8)    ! thermal expansion coefficient
      !----------------------------------------------------------------------
      ! BER EXPANSION PARAMETERS
      !----------------------------------------------------------------------
      Ea = PROPS(9)          ! BER activation energy
      beta_B = PROPS(10)     ! BER expansion fitting parameter 1
      alpha_NB = PROPS(11)   ! BER expansion fitting parameter 2
      !----------------------------------------------------------------------
      ! TRANSITION TEMPS
      !----------------------------------------------------------------------
      Tg = PROPS(12)         ! glass transition temp
      Tv = PROPS(13)         ! vitrimeric transition temperature
      !----------------------------------------------------------------------
      ! CURE PARAMETERS
      !----------------------------------------------------------------------
      phidot = PROPS(14)     ! cure rate
      alpha_C = PROPS(15)    ! cure expansion coefficient
	  AF = PROPS(16)         ! The parameter for shift factor when it below the reference temp 
	  temp_ref = PROPS(17)   ! Reference temperature for shift factor-fitting paramter
      phi0 = 0
	  Tcure = 163.1D0
	  
	  E3 = PROPS(18)
	  tau3_T0 = PROPS(19)
C

C
      !--------------------------------------------------------------------------
      ! INITIALIZE PARAMETERS/CONSTANTS
      !--------------------------------------------------------------------------
      CALL ONEM(I1)
C      DO I=1,NDI
C          DO J=1,NDI
C              I1(I,J) = ZERO
C          ENDDO
C          I1(I,I) = ONE
C      ENDDO
C
      !--------------------------------------------------------------------------
      ! DEFINE TANGENT MATRIX - LINEAR ELASTIC FOR SMALL STRAIN
      !--------------------------------------------------------------------------
	  mu1 = E1/2/(1+nu)
	  mu2 = E2/2/(1+nu)
	  mu3 = E3/2/(1+nu)
      kappa1 = E1/(THREE*(ONE-TWO*nu))
      kappa2 = E2/(THREE*(ONE-TWO*nu))
	  kappa3 = E3/(THREE*(ONE-TWO*nu))

      CALL ASSEMBLEDDSDDE(E1,nu,D1)
      CALL ASSEMBLEDDSDDE(E2,nu,D2)
	  CALL ASSEMBLEDDSDDE(E3,nu,D3)
	  
      DDSDDE = (D1+D2+D3)
      !----------------------------------------------------------------------
      ! DEFINE STATE VARIABLES
      !----------------------------------------------------------------------
      IF (TIME(2)==ZERO) THEN
          Fv_old = I1
		  Fv3_old = I1
          Fb_old = I1
          Falpha_old = I1
          CI = ZERO
          CR = ZERO
          phi_old = ZERO
          Tf_B0 = Tcure
          Tf_T0 = Tcure
		  Fcure_old = I1
		  Fthm_old = I1
		  FBER_old = I1
		  alpha_B0 = ONE
		  alpha_T0 = ONE
		  deps_T0 = ZERO
		  deps_Tg0 = ZERO
		  deps_Tv0 = ZERO
		  alpha_NT0 = PROPS(8)
      ELSE
          Fv_old(1,1)=STATEV(1)
          Fv_old(2,2)=STATEV(2)
          Fv_old(3,3)=STATEV(3)
          Fv_old(1,2)=STATEV(4)
          Fv_old(2,3)=STATEV(5)
          Fv_old(3,1)=STATEV(6)
          Fv_old(2,1)=STATEV(7)
          Fv_old(3,2)=STATEV(8)
          Fv_old(1,3)=STATEV(9)

          Fb_old(1,1)=STATEV(13)
          Fb_old(2,2)=STATEV(14)
          Fb_old(3,3)=STATEV(15)
          Fb_old(1,2)=STATEV(16)
          Fb_old(2,3)=STATEV(17)
          Fb_old(3,1)=STATEV(18)
          Fb_old(2,1)=STATEV(19)
          Fb_old(3,2)=STATEV(20)
          Fb_old(1,3)=STATEV(21)
C        
          CI=STATEV(22)
          CR=STATEV(23)
C
          phi_old=STATEV(24)
C
          Falpha_old(1,1)=STATEV(26)
          Falpha_old(2,2)=STATEV(27)
          Falpha_old(3,3)=STATEV(28)
          Falpha_old(1,2)=STATEV(29)
          Falpha_old(2,3)=STATEV(30)
          Falpha_old(3,1)=STATEV(31)
          Falpha_old(2,1)=STATEV(32)
          Falpha_old(3,2)=STATEV(33)
          Falpha_old(1,3)=STATEV(34)
C
          Tf_B0 = STATEV(35)
          Tf_T0 = STATEV(36)
C
          Fcure_old(1,1) = STATEV(95)
          Fcure_old(2,2) = STATEV(96)
          Fcure_old(3,3) = STATEV(97)
          Fcure_old(1,2) = STATEV(98)  
          Fcure_old(2,3) = STATEV(99) 
          Fcure_old(3,1) = STATEV(100) 
          Fcure_old(2,1) = STATEV(101)  
          Fcure_old(3,2) = STATEV(102) 
          Fcure_old(1,3) = STATEV(103) 
C
          FBER_old(1,1) = STATEV(104) 
          FBER_old(2,2) = STATEV(105) 
          FBER_old(3,3) = STATEV(106)  
          FBER_old(1,2) = STATEV(107) 
          FBER_old(2,3) = STATEV(108) 
          FBER_old(3,1) = STATEV(109) 
          FBER_old(2,1) = STATEV(110)  
          FBER_old(3,2) = STATEV(111)  
          FBER_old(1,3) = STATEV(112) 
C
          Fthm_old(1,1) = STATEV(113)  
          Fthm_old(2,2) = STATEV(114)  
          Fthm_old(3,3) = STATEV(115)  
          Fthm_old(1,2) = STATEV(116)  
          Fthm_old(2,3) = STATEV(117)  
          Fthm_old(3,1) = STATEV(118)  
          Fthm_old(2,1) = STATEV(119)  
          Fthm_old(3,2) = STATEV(120)  
          Fthm_old(1,3) = STATEV(121) 
          alpha_B0 = STATEV(122)
		  alpha_T0 = STATEV(123)
	      alpha_NT0 = STATEV(130)		  
		  deps_T0 = STATEV(131)		  
C		  
          Fv3_old(1,1)=STATEV(134)
          Fv3_old(2,2)=STATEV(135)
          Fv3_old(3,3)=STATEV(136)
          Fv3_old(1,2)=STATEV(137)
          Fv3_old(2,3)=STATEV(138)
          Fv3_old(3,1)=STATEV(139)
          Fv3_old(2,1)=STATEV(140)
          Fv3_old(3,2)=STATEV(141)
          Fv3_old(1,3)=STATEV(142)		  
      ENDIF
C
      CALL COMPUTEALPHANT(TEMP,Tg,Tv,alpha_NT)! remember to Turn off BER for 0% Catalyst      
      ! Denote new step for Element 1, IP 1
      IF (NPT==ONE .AND. NOEL==ONE) THEN
      PRINT*, '%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%'
      PRINT*,'NEW STEP: Time = ',TIME(2)
      PRINT*,'KINC = ', KINC,', NPT = ',NPT
      PRINT*, '  '
      ENDIF
C
      !----------------------------------------------------------------------
      ! DEFINE CURE CONTRACTION STRAIN (DEGREE OF CURE GOES TOWARDS MAX 
      ! DEGREE, phi = 1 - INDUCES CONTRACTION)
      !----------------------------------------------------------------------
      phi = (DTIME*phidot+phi_old)/(1+DTIME*phidot)
	  deps_C = -alpha_C*phi   
C      
      !----------------------------------------------------------------------
      ! STEP 1: compute shift factors
      !----------------------------------------------------------------------
      TEMP0 = TEMP-DTEMP
  
      IF (TEMP>=Tv) THEN
          alpha_B = 10**(Ea/R*(1/TEMP-1/Tv))
	  ELSE 
	      alpha_B = one
      ENDIF      
C      
      !----------------------------------------------------------------------
      ! STEP 2: compute retardation times
      !----------------------------------------------------------------------	  
      tau_B = tau_B0*alpha_B      
C         
      !--------------------------------------------------------------------------
      ! STEP 3: COMPUTE THERMAL AND BER STRAINS AS FUNCTION OF dt AND dT
	  !   update with current value 
      !--------------------------------------------------------------------------

	  deps_T = deps_T0 + alpha_NT*(theta_tau-TEMP)	  
	  deps_B = beta_B*(EXP(alpha_NB*(Tf_B-Tv))-1)
C
      !----------------------------------------------------------------------
      ! DEFINE THE INCREMENTAL, AND TOTAL EFFECTIVE DEFORMATIONS DUE TO CURE,
      ! BERS, AND DTEMP - GRADIENTS ARE FROM TIME t->t+dt
      !----------------------------------------------------------------------
      Fcure_tau = EXP(deps_C)*I1
	  Fcure = Fcure_tau

      FBER_tau = EXP(deps_B)*I1
	  FBER = FBER_tau   
	  
      Fthm_tau = EXP(deps_T)*I1
	  Fthm = Fthm_tau
	  
	  Falpha_tau = MATMUL(Fcure_tau,Fthm_tau)
	  Falpha = Falpha_tau
C
      !----------------------------------------------------------------------
      ! USE CURE/THERMAL DEFORMATION AND TRIAL DEFORMATION TO COMPUTE 
      ! THE EFFECTIVE INTERNAL DEFORMATIONS (TOTAL, VISCOELASTIC AND ELASTIC)
      !----------------------------------------------------------------------
      ! Clean up DFGRD1 - eliminate numerical rounding error
C      DO I=1,NDI
C          DO J=1,NDI
C              IF (DFGRD1(I,J)<1E-15) THEN
C                  DFGRD1(I,J) = ZERO
C              ENDIF
C          ENDDO
C      ENDDO 
C
      CALL M3INV(Falpha,invFalpha)
      F_eff = MATMUL(DFGRD1,invFalpha)
      CALL M3INV(Fb_old,invFb_old)
      Fve_old = MATMUL(F_eff,invFb_old)
      CALL M3INV(Fv_old,invFv_old)
      Fe_old = MATMUL(Fve_old,invFv_old)
      CALL M3INV(Fv3_old,invFv3_old)      ! calculate the initial F for leg 3
      Fe3_old = MATMUL(Fve_old,invFv3_old)	  
C
      !----------------------------------------------------------------------
      ! COMPUTE PRE-RELAXATION TRIAL STRESS AND STRAIN
      !----------------------------------------------------------------------
      ! In elastic leg (single spring)
      ! Compute strain in reference frame (must intermediately compute V in priniciple orientation)
      CALL LEFTDECOMP(Fve_old,Qve,Vve,Eve)      
C      
      ! Compute stress in elastic leg  1
      trEve = Eve(1,1) + Eve(2,2) + Eve(3,3)
      Eve_d = Eve - THIRD*trEve*I1
      Te = TWO*mu1*Eve_d + kappa1*trEve*I1
C         
      ! In viscoelastic leg (spring + damper)
      ! Compute strain in reference frame (must intermediately compute V in priniciple orientation)
      CALL LEFTDECOMP(Fe_old,Qe,Ve,Ee)      
C
      ! Compute stress in viscoelastic leg 2
      trEe = Ee(1,1) + Ee(2,2) + Ee(3,3)
      Ee_d = Ee - THIRD*trEe*I1
      Tve = TWO*mu2*Ee_d + kappa2*trEe*I1  
C	  
      CALL LEFTDECOMP(Fe3_old,Qe3,Ve3,Ee3) 
      ! Compute stress in viscoelastic leg 3
      trEe3 = Ee3(1,1) + Ee3(2,2) + Ee3(3,3)
      Ee3_d = Ee3 - THIRD*trEe3*I1
      Tve3 = TWO*mu3*Ee3_d + kappa3*trEe3*I1  	  
	  
C
      ! Compute total stress
  
      T_old = Te + Tve + Tve3      
C
      !----------------------------------------------------------------------
      ! UPDATE VISCOUS DEFORMATION, Fv, USING LINEAR VISCOELASTICITY  leg 2
      !----------------------------------------------------------------------
      ! Put stress into principle orientation to compute flow
      CALL SPECTRAL(Tve,SIGP,Qve)
      CALL MTRANS(Qve,QveT)
      CALL ZEROM(Tve)
      Tve(1,1) = SIGP(1)
      Tve(2,2) = SIGP(2)
      Tve(3,3) = SIGP(3)
C      
      trTve = Tve(1,1) + Tve(2,2) + Tve(3,3)
      Tve_d = Tve - THIRD*trTve*I1      
C
      ! Compute flow rate
      sve = (ONE/SQRT(TWO))*SQRT(SUM(MATMUL(Tve_d,Tve_d))) 
      gdotv = sve/(mu2*ABS(tau2_T))     
C
      ! Compute rate of deformation tensor
      CALL ZEROM(Dvp)
      IF (sve>ZERO) THEN
          Dvp = gdotv*Tve_d/(SQRT(TWO)*sve) 
      ELSE
          Dvp = ZERO*I1
      ENDIF  
C 
      ! Define exponential factor
      CALL ZEROM(Dv_exp)
      Dv_exp(1,1) = EXP(Dvp(1,1)*DTIME)
      Dv_exp(2,2) = EXP(Dvp(2,2)*DTIME)
      Dv_exp(3,3) = EXP(Dvp(3,3)*DTIME)   
C

      Dv_exp = MATMUL(MATMUL(Qve,Dv_exp),QveT)
C
      ! Update viscous deformation
      Fv_new = MATMUL(Dv_exp,Fv_old)
C
      ! Update other deformation gradients given Fv_new
      CALL M3INV(Fb_old,invFb_old)
      Fve_int = MATMUL(F_eff,invFb_old)
      CALL M3INV(Fv_new,invFv_new)
      Fe_int = MATMUL(Fve_int,invFv_new)
C
C      IF (NPT==ONE .AND. NOEL==ONE) THEN
C          PRINT*, 'Fv_old = ', Fv_old
C          PRINT*, 'Fv_new = ', Fv_new
C          PRINT*, 'Fb_old = ', Fb_old
C          PRINT*, 'invFb_old = ', invFb_old
C          PRINT*, 'Fve_int = ', Fve_int
C          PRINT*, 'F_eff = ', F_eff
C          PRINT*, 'Fv_new = ', Fv_new
C          PRINT*, 'invFv_new = ', invFv_new
C          PRINT*, 'Fve_int = ', Fve_int
C          PRINT*, 'Fe_int = ', Fe_int
C      ENDIF
C
      !----------------------------------------------------------------------
      ! UPDATE VISCOUS DEFORMATION, Fv3, USING LINEAR VISCOELASTICITY  leg 3
      !----------------------------------------------------------------------
      ! Put stress into principle orientation to compute flow
      CALL SPECTRAL(Tve3,SIGP3,Qve3)
      CALL MTRANS(Qve3,QveT3)
      CALL ZEROM(Tve3)
      Tve3(1,1) = SIGP3(1)
      Tve3(2,2) = SIGP3(2)
      Tve3(3,3) = SIGP3(3)
C      
      trTve3 = Tve3(1,1) + Tve3(2,2) + Tve3(3,3)
      Tve3_d = Tve3 - THIRD*trTve3*I1      
C
      ! Compute flow rate
      sve3 = (ONE/SQRT(TWO))*SQRT(SUM(MATMUL(Tve3_d,Tve3_d))) 
      gdotv3 = sve3/(mu3*ABS(tau3_T))     
C
      ! Compute rate of deformation tensor
      CALL ZEROM(Dvp3)
      IF (sve3>ZERO) THEN
          Dvp3 = gdotv3*Tve3_d/(SQRT(TWO)*sve3)   ! sqrt(two) or just two
		  !Dvp = ZERO*I1   ! turn off the plastic flow
      ELSE
          Dvp3 = ZERO*I1
      ENDIF  
C 
      ! Define exponential factor
      CALL ZEROM(Dv3_exp)
      Dv3_exp(1,1) = EXP(Dvp3(1,1)*DTIME)
      Dv3_exp(2,2) = EXP(Dvp3(2,2)*DTIME)
      Dv3_exp(3,3) = EXP(Dvp3(3,3)*DTIME)   
C
      ! Rotate back into reference orientation
      Dv3_exp = MATMUL(MATMUL(Qve3,Dv3_exp),QveT3)
C
      ! Update viscous deformation
      Fv3_new = MATMUL(Dv3_exp,Fv3_old)
C
      ! Update other deformation gradients given Fv_new
      CALL M3INV(Fb_old,invFb_old)
      Fve_int = MATMUL(F_eff,invFb_old)
      CALL M3INV(Fv3_new,invFv3_new)
      Fe3_int = MATMUL(Fve_int,invFv3_new)
C
C      IF (NPT==ONE .AND. NOEL==ONE) THEN
C          PRINT*, 'Fv_old = ', Fv_old
C          PRINT*, 'Fv_new = ', Fv_new
C          PRINT*, 'Fb_old = ', Fb_old
C          PRINT*, 'invFb_old = ', invFb_old
C          PRINT*, 'Fve_int = ', Fve_int
C          PRINT*, 'F_eff = ', F_eff
C          PRINT*, 'Fv_new = ', Fv_new
C          PRINT*, 'invFv_new = ', invFv_new
C          PRINT*, 'Fve_int = ', Fve_int
C          PRINT*, 'Fe_int = ', Fe_int
C      ENDIF
C
      !---------------------------------------------------------------------
      ! FIND INTTERMEDIATE STRESS AND STRAIN GIVEN VISCOUS DISSIPATION
      !----------------------------------------------------------------------
      ! In elastic leg (single spring)
      ! Compute strain in reference frame (must intermediately compute V in priniciple orientation)
      CALL LEFTDECOMP(Fve_int,Qve,Vve,Eve)      
      ! Compute stress in elastic leg
      trEve = Eve(1,1) + Eve(2,2) + Eve(3,3)
      Eve_d = Eve - THIRD*trEve*I1
      Te = TWO*mu1*Eve_d + kappa1*trEve*I1      
C    
      ! In viscoelastic leg (spring + damper)
      ! Compute strain in reference frame (must intermediately compute V in priniciple orientation)
      CALL LEFTDECOMP(Fe_int,Qe,Ve,Ee)      
      ! Compute stress in viscoelastic leg
      trEe = Ee(1,1) + Ee(2,2) + Ee(3,3)
      Ee_d = Ee - THIRD*trEe*I1
      Tve = TWO*mu2*Ee_d + kappa2*trEe*I1    
C
      CALL LEFTDECOMP(Fe3_int,Qe3,Ve3,Ee3)      
      ! Compute stress in viscoelastic leg
      trEe3 = Ee3(1,1) + Ee3(2,2) + Ee3(3,3)
      Ee3_d = Ee3 - THIRD*trEe3*I1
      Tve3 = TWO*mu3*Ee3_d + kappa3*trEe3*I1 
C	  
      T_int = Te + Tve + Tve3     
C      
      IF (NPT==ONE .AND. NOEL==ONE) THEN
          PRINT*, 'Te = ', Te
          PRINT*, 'Tve = ', Tve
          PRINT*, 'T_int = ', T_int
          PRINT*, 'Ee_d = ', Ee_d
          PRINT*, 'Ee = ', Ee
          PRINT*, 'Eve_d = ', Ee_d
          PRINT*, 'Eve = ', Ee
      ENDIF
C 
      !----------------------------------------------------------------------
      ! NOW USE INTERMEDIATE STRESS TO UPDATE THE EFFECTIVE BER RELAXATION
      !----------------------------------------------------------------------
      ! Put stress into principle orientation to compute flow
      CALL SPECTRAL(T_int,SIGP,Q)
      CALL MTRANS(Q,QT)
      CALL ZEROM(Tp)
      Tp(1,1) = SIGP(1)
      Tp(2,2) = SIGP(2)
      Tp(3,3) = SIGP(3)
C     
      trT_int = Tp(1,1) + Tp(2,2) + Tp(3,3)
      Tp_d = Tp - THIRD*trT_int*I1      
C
      ! Compute flow rate
      sb = (1/SQRT(TWO))*SQRT(SUM(MATMUL(Tp_d,Tp_d)))
      gdotb = sb/((mu1+mu2+mu3)*tau_B)    
C
      ! Compute rate of deformation tensor for BER
	  ! Turn on this for BER case
      IF (sb>0 .AND. TEMP>500) THEN
          Dbp = gdotb*Tp_d/(SQRT(TWO)*sb) !Delta F viscous - correction made to denominator for units
      ELSE
          Dbp = ZERO*I1
      ENDIF
C 
      ! Define exponential factor
      CALL ZEROM(Db_exp)
      Db_exp(1,1) = EXP(Dbp(1,1)*DTIME)
      Db_exp(2,2) = EXP(Dbp(2,2)*DTIME)
      Db_exp(3,3) = EXP(Dbp(3,3)*DTIME)   
C
      ! Rotate back into reference orientation
      Db_exp = MATMUL(MATMUL(Q,Db_exp),QT)
C
      ! Update BER deformation
      Fb_new = MATMUL(Db_exp,Fb_old)
C      
      ! Update other deformation gradients given Fv_new
      CALL M3INV(Fb_new,invFb_new)
      Fve_new = MATMUL(F_eff,invFb_new)
      CALL M3INV(Fv_new,invFv_new)
      Fe_new = MATMUL(Fve_new,invFv_new)
      CALL M3INV(Fv3_new,invFv3_new)
      Fe3_new = MATMUL(Fve_new,invFv3_new)	  
C
      !----------------------------------------------------------------------
      ! COMPUTE NEW STRESS & DEF. FOR BOTH LEGS GIVEN VISCOUS + BER DISSIPATION
      !----------------------------------------------------------------------
      ! In elastic leg (single spring)
      ! Compute strain in reference frame (must intermediately compute V in priniciple orientation)
      CALL LEFTDECOMP(Fve_new,Qve,Vve,Eve)      
      ! Compute stress in elastic leg
      trEve = Eve(1,1) + Eve(2,2) + Eve(3,3)
      Eve_d = Eve - THIRD*trEve*I1
      Te = TWO*mu1*Eve_d + kappa1*trEve*I1      
C    
      ! In viscoelastic leg (spring + damper)
      ! Compute strain in reference frame (must intermediately compute V in priniciple orientation)
      CALL LEFTDECOMP(Fe_new,Qe,Ve,Ee)      
      ! Compute stress in viscoelastic leg
      trEe = Ee(1,1) + Ee(2,2) + Ee(3,3)
      Ee_d = Ee - THIRD*trEe*I1
      Tve = TWO*mu2*Ee_d + kappa2*trEe*I1    
C
      ! Output F_new for state var
      F_new = MATMUL(MATMUL(Fb_new,Fve_new),Falpha)	  
C
	  ! Leg 3 	
      CALL LEFTDECOMP(Fe3_new,Qe3,Ve3,Ee3)      
      ! Compute stress in viscoelastic leg3
      trEe3 = Ee3(1,1) + Ee3(2,2) + Ee3(3,3)
      Ee3_d = Ee3 - THIRD*trEe3*I1
      Tve3 = TWO*mu3*Ee3_d + kappa3*trEe3*I1    
C
      !----------------------------------------------------------------------
      ! COMPUTE TOTAL STRESS
      !----------------------------------------------------------------------
      T_new = Te + Tve + Tve3
C
      !----------------------------------------------------------------------
      ! PRINT DEBUGGING OUTPUTS
      !----------------------------------------------------------------------
      IF (NPT==ONE .AND. NOEL==ONE) THEN
c      PRINT*, 'D_v = ', Dv
c      PRINT*, 'depsC = ',deps_C
c      PRINT*, 'STRAN11 = ',STRAN(1)
c      PRINT*, 'STRAN22 = ',STRAN(2)
c      PRINT*, 'STRAN33 = ',STRAN(3)
c      PRINT*, 'STRAN12 = ',STRAN(4)
c      PRINT*, 'STRAN23 = ',STRAN(5)
c      PRINT*, 'STRAN31 = ',STRAN(6)
c      PRINT*, 'DSTRAN11 = ',DSTRAN(1)
c      PRINT*, 'DSTRAN22 = ',DSTRAN(2)
c      PRINT*, 'DSTRAN33 = ',DSTRAN(3)
c      PRINT*, 'DSTRAN12 = ',DSTRAN(4)
c      PRINT*, 'DSTRAN23 = ',DSTRAN(5)
c      PRINT*, 'DSTRAN31 = ',DSTRAN(6)
c      PRINT*, 'F0_11 = ',DFGRD0(1,1)
c      PRINT*, 'F0_22 = ',DFGRD0(2,2)
c      PRINT*, 'F0_33 = ',DFGRD0(3,3)
c      PRINT*, 'F0_12 = ',DFGRD0(1,2)
c      PRINT*, 'F0_23 = ',DFGRD0(2,3)
c      PRINT*, 'F0_31 = ',DFGRD0(3,1)
c      PRINT*, 'DFGRD1_11 = ',DFGRD1(1,1)
c      PRINT*, 'DFGRD1_22 = ',DFGRD1(2,2)
c      PRINT*, 'DFGRD1_33 = ',DFGRD1(3,3)
c      PRINT*, 'DFGRD1_12 = ',DFGRD1(1,2)
c      PRINT*, 'DFGRD1_23 = ',DFGRD1(2,3)
c      PRINT*, 'DFGRD1_31 = ',DFGRD1(3,1)
c      PRINT*, 'Tv = ',Tv
c      PRINT*, 'Tg = ',Tg
c      PRINT*, 'TfT = ',Tf_T
c      PRINT*, 'TfB = ',Tf_B
c      PRINT*, 'phi = ',phi,', phi0 = ',phi_old
C
c      PRINT*, 'Falph = ', Falpha(1,1)
c      PRINT*, 'D_v = ', Dv(1,1)
c      PRINT*, 'D_b = ', Db(1,1)
c      PRINT*, 'Falph = ', Falpha(1,1)
C     
c      PRINT*, 'T11 = ', T_new(1,1)
c      PRINT*, 'T22 = ', T_new(2,2)
c      PRINT*, 'T33 = ', T_new(3,3)
c      PRINT*, 'T12 = ', T_new(1,2)
c      PRINT*, 'T23 = ', T_new(2,3)
c      PRINT*, 'T31 = ', T_new(3,1)
c      PRINT*, 'T12 = ', T_new(1,2)
c      PRINT*, 'T23 = ', T_new(2,3)
c      PRINT*, 'T32 = ', T_new(3,1)
C      
c      PRINT*, 'TL1_11 = ', Te(1,1)
c      PRINT*, 'TL1_22 = ', Te(2,2)
c      PRINT*, 'TL1_33 = ', Te(3,3)
c      PRINT*, 'TL1_12 = ', Te(1,2)
c      PRINT*, 'TL1_23 = ', Te(2,3)
c      PRINT*, 'TL1_31 = ', Te(3,1)
C
c      PRINT*, 'TL2_11 = ', Tve(1,1)
c      PRINT*, 'TL2_22 = ', Tve(2,2)
c      PRINT*, 'TL2_33 = ', Tve(3,3)
c      PRINT*, 'TL2_12 = ', Tve(1,2)
c      PRINT*, 'TL2_23 = ', Tve(2,3)
c      PRINT*, 'TL2_31 = ', Tve(3,1)
C      
c      PRINT*, 'F11 = ', F_eff(1,1)
c      PRINT*, 'Fve = ', Fve_new(1,1)
c      PRINT*, 'Fe11 = ', Fe_new(1,1)
c      PRINT*, 'Ehe11 = ', Ehe(1,1)
c      PRINT*, 'Eve11 = ', Eve(1,1)
c      PRINT*, 'Fe33 = ', Fe_new(3,3)
c      PRINT*, 'Ehe33 = ', Ehe(3,3)
c      PRINT*, 'Eve33 = ', Eve(3,3)
c      PRINT*, 'Fv = ', Fv_new(1,1)  
      PRINT*, '%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%'
      ENDIF
      
      !----------------------------------------------------------------------
      ! REDEFINE STRESS
      !----------------------------------------------------------------------
      STRESS(1) = T_new(1,1)
      STRESS(2) = T_new(2,2)
      STRESS(3) = T_new(3,3)
      STRESS(4) = T_new(1,2)
      STRESS(5) = T_new(1,3)
      STRESS(6) = T_new(2,3)      
C
      !----------------------------------------------------------------------
      ! REDEFINE STATE VARIABLES
      !----------------------------------------------------------------------
      STATEV(1) = Fv_new(1,1)
      STATEV(2) = Fv_new(2,2)
      STATEV(3) = Fv_new(3,3)
      STATEV(4) = Fv_new(1,2)
      STATEV(5) = Fv_new(2,3)
      STATEV(6) = Fv_new(3,1)
      STATEV(7) = Fv_new(2,1)
      STATEV(8) = Fv_new(3,2)
      STATEV(9) = Fv_new(1,3)
C
      STATEV(10) = F_new(1,1)
      STATEV(11) = F_new(2,2)
      STATEV(12) = F_new(3,3)
C
      STATEV(13) = Fb_new(1,1)
      STATEV(14) = Fb_new(2,2)
      STATEV(15) = Fb_new(3,3)
      STATEV(16) = Fb_new(1,2)
      STATEV(17) = Fb_new(2,3)
      STATEV(18) = Fb_new(3,1)
      STATEV(19) = Fb_new(2,1)
      STATEV(20) = Fb_new(3,2)
      STATEV(21) = Fb_new(1,3)
C
      STATEV(22) = C_I
      STATEV(23) = C_R
C
      STATEV(24) = phi
C
      STATEV(26) = Falpha(1,1)
      STATEV(27) = Falpha(2,2)
      STATEV(28) = Falpha(3,3)
      STATEV(29) = Falpha(1,2)
      STATEV(30) = Falpha(2,3)
      STATEV(31) = Falpha(3,1)
      STATEV(32) = Falpha(2,1)
      STATEV(33) = Falpha(3,2)
      STATEV(34) = Falpha(1,3)
C
      STATEV(35) = Tf_B
      STATEV(36) = Tf_T
C
      STATEV(37) = Te(1,1)
      STATEV(38) = Te(2,2)
      STATEV(39) = Te(3,3)
      STATEV(40) = Te(1,2)
      STATEV(41) = Te(2,3)
      STATEV(42) = Te(3,1)
C
      STATEV(43) = Tve(1,1)
      STATEV(44) = Tve(2,2)
      STATEV(45) = Tve(3,3)
      STATEV(46) = Tve(1,2)
      STATEV(47) = Tve(2,3)
      STATEV(48) = Tve(3,1)
C
      STATEV(49) = sve
      STATEV(50) = sb
C
      STATEV(51) = tau_T
      STATEV(52) = tau_B
C
      STATEV(53) = Fve_new(1,1)
      STATEV(54) = Fve_new(2,2)
      STATEV(55) = Fve_new(3,3)
      STATEV(56) = Fve_new(1,2)
      STATEV(57) = Fve_new(2,3)
      STATEV(58) = Fve_new(3,1)
      STATEV(59) = Fve_new(2,1)
      STATEV(60) = Fve_new(3,2)
      STATEV(61) = Fve_new(1,3)
C
      STATEV(62) = Fe_new(1,1)
      STATEV(63) = Fe_new(2,2)
      STATEV(64) = Fe_new(3,3)
      STATEV(65) = Fe_new(1,2)
      STATEV(66) = Fe_new(2,3)
      STATEV(67) = Fe_new(3,1)
      STATEV(68) = Fe_new(2,1)
      STATEV(69) = Fe_new(3,2)
      STATEV(70) = Fe_new(1,3)
C
      STATEV(71) = Ee(1,1)
      STATEV(72) = Ee(2,2)
      STATEV(73) = Ee(3,3)
      STATEV(74) = Ee(1,2)
      STATEV(75) = Ee(2,3)
      STATEV(76) = Ee(3,1)
      STATEV(77) = Ee(2,1)
      STATEV(78) = Ee(3,2)
      STATEV(79) = Ee(1,3)
C
      STATEV(80) = Eve(1,1)
      STATEV(81) = Eve(2,2)
      STATEV(82) = Eve(3,3)
      STATEV(83) = Eve(1,2)
      STATEV(84) = Eve(2,3)
      STATEV(85) = Eve(3,1)
      STATEV(86) = Eve(2,1)
      STATEV(87) = Eve(3,2)
      STATEV(88) = Eve(1,3)
C
      STATEV(89) = Te(1,1) + Tve(1,1) + Tve3(1,1)
      STATEV(90) = Te(2,2) + Tve(2,2) + Tve3(2,2)
      STATEV(91) = Te(3,3) + Tve(3,3) + Tve3(3,3)
      STATEV(92) = Te(1,2) + Tve(1,2) + Tve3(1,2)
      STATEV(93) = Te(2,3) + Tve(2,3) + Tve3(2,3)
      STATEV(94) = Te(3,1) + Tve(3,1) + Tve3(3,1)
C
      STATEV(95) = Fcure(1,1)
      STATEV(96) = Fcure(2,2)
      STATEV(97) = Fcure(3,3)
      STATEV(98) = Fcure(1,2)
      STATEV(99) = Fcure(2,3)
      STATEV(100) = Fcure(3,1)
      STATEV(101) = Fcure(2,1)
      STATEV(102) = Fcure(3,2)
      STATEV(103) = Fcure(1,3)
C
      STATEV(104) = FBER(1,1)
      STATEV(105) = FBER(2,2)
      STATEV(106) = FBER(3,3)
      STATEV(107) = FBER(1,2)
      STATEV(108) = FBER(2,3)
      STATEV(109) = FBER(3,1)
      STATEV(110) = FBER(2,1)
      STATEV(111) = FBER(3,2)
      STATEV(112) = FBER(1,3)
C
      STATEV(113) = Fthm(1,1)
      STATEV(114) = Fthm(2,2)
      STATEV(115) = Fthm(3,3)
      STATEV(116) = Fthm(1,2)
      STATEV(117) = Fthm(2,3)
      STATEV(118) = Fthm(3,1)
      STATEV(119) = Fthm(2,1)
      STATEV(120) = Fthm(3,2)
      STATEV(121) = Fthm(1,3)
	  STATEV(122) = alpha_B
	  STATEV(123) = alpha_T
C
      STATEV(124) = mu1
      STATEV(125) = mu2
      STATEV(126) = E1
      STATEV(127) = E2
      STATEV(128) = E3
      STATEV(129) = kappa2
	  STATEV(130) = alpha_NT
	  STATEV(131) = deps_T
	  STATEV(132) = gdotv
	  !STATEV(133) = deps_Tv
      STATEV(134) = Fv3_new(1,1)
      STATEV(135) = Fv3_new(2,2)
      STATEV(136) = Fv3_new(3,3)
      STATEV(137) = Fv3_new(1,2)
      STATEV(138) = Fv3_new(2,3)
      STATEV(139) = Fv3_new(3,1)
      STATEV(140) = Fv3_new(2,1)
      STATEV(141) = Fv3_new(3,2)
      STATEV(142) = Fv3_new(1,3)
      STATEV(143) = Fe3_new(1,1)
      STATEV(144) = Fe3_new(2,2)
      STATEV(145) = Fe3_new(3,3)
      STATEV(146) = Fe3_new(1,2)
      STATEV(147) = Fe3_new(2,3)
      STATEV(148) = Fe3_new(3,1)
      STATEV(149) = Fe3_new(2,1)
      STATEV(150) = Fe3_new(3,2)
      STATEV(151) = Fe3_new(1,3)	
C
      STATEV(152) = Tve3(1,1)
      STATEV(153) = Tve3(2,2)
      STATEV(154) = Tve3(3,3)
      STATEV(155) = Tve3(1,2)
      STATEV(156) = Tve3(2,3)
      STATEV(157) = Tve3(3,1)	  
      RETURN
      END
      
C**********************************************************************
      SUBROUTINE M3INV(A,AINV)

C 	  THIS SUBROUTINE CALCULATES THE THE INVERSE OF A 3 BY 3 MATRIX
C	      [A] AND PLACES THE RESULT IN [AINV]. 
C 	  IF DET(A) IS ZERO, THE CALCULATION
C 	  IS TERMINATED AND A DIAGNOSTIC STATEMENT IS PRINTED.
C**********************************************************************

      REAL*8  A(3,3), AINV(3,3), DET, ACOFAC(3,3), AADJ(3,3)

C	    A(3,3)	        -- THE MATRIX WHOSE INVERSE IS DESIRED.
C	    DET		-- THE COMPUTED DETERMINANT OF [A].
C	    ACOFAC(3,3)	-- THE MATRIX OF COFACTORS OF A(I,J).
C	    		   THE SIGNED MINOR (-1)**(I+J)*M_IJ
C	    		   IS CALLED THE COFACTOR OF A(I,J).
C	    AADJ(3,3)	-- THE ADJOINT OF [A]. IT IS THE MATRIX
C	    		   OBTAINED BY REPLACING EACH ELEMENT OF
C	    		   [A] BY ITS COFACTOR, AND THEN TAKING
C	    		   TRANSPOSE OF THE RESULTING MATRIX.
C	    AINV(3,3)	-- RETURNED AS INVERSE OF [A].
C	    		   [AINV] = [AADJ]/DET.
C----------------------------------------------------------------------
      CALL MDET(A,DET)
C      IF ( DET .EQ. 0.D0 ) THEN
C        WRITE(80,100)
C        WRITE(*,100)
C        STOP
C      ENDIF
      CALL MCOFAC(A,ACOFAC)
      CALL MTRANS(ACOFAC,AADJ)
      DO I = 1,3
          DO J = 1,3
              AINV(I,J) = AADJ(I,J)/DET
          ENDDO
      ENDDO
!1      CONTINUE
!2    CONTINUE

 !    100  FORMAT(5X,'--ERROR IN M3INV--- THE MATRIX IS SINGULAR',/,
 !    +       10X,'PROGRAM TERMINATED')

      RETURN
      END
C
C**********************************************************************
      SUBROUTINE MCOFAC(A,ACOFAC)
 
C 	  THIS SUBROUTINE CALCULATES THE COFACTOR OF A 3 BY 3 MATRIX [A],
C 	  AND PLACES THE RESULT IN [ACOFAC]. 
C**********************************************************************
      REAL*8  A(3,3), ACOFAC(3,3)

      ACOFAC(1,1) = A(2,2)*A(3,3) - A(3,2)*A(2,3)
      ACOFAC(1,2) = -(A(2,1)*A(3,3) - A(3,1)*A(2,3))
      ACOFAC(1,3) = A(2,1)*A(3,2) - A(3,1)*A(2,2)
      ACOFAC(2,1) = -(A(1,2)*A(3,3) - A(3,2)*A(1,3))
      ACOFAC(2,2) = A(1,1)*A(3,3) - A(3,1)*A(1,3)
      ACOFAC(2,3) = -(A(1,1)*A(3,2) - A(3,1)*A(1,2))
      ACOFAC(3,1) = A(1,2)*A(2,3)  - A(2,2)*A(1,3)
      ACOFAC(3,2) = -(A(1,1)*A(2,3) - A(2,1)*A(1,3))
      ACOFAC(3,3) = A(1,1)*A(2,2) - A(2,1)*A(1,2)

      RETURN
      END
C
C**********************************************************************
      SUBROUTINE MTRANS(A,ATRANS)
 
C	    THIS SUBROUTINE CALCULATES THE TRANSPOSE OF A 3 BY 3 
C	    MATRIX [A], AND PLACES THE RESULT IN ATRANS. 
C**********************************************************************
      REAL*8 A(3,3),ATRANS(3,3)

      DO 2 I=1,3
        DO 1 J=1,3
          ATRANS(J,I) = A(I,J)
 1      CONTINUE
 2    CONTINUE

      RETURN
      END
C
C**********************************************************************
      SUBROUTINE MDET(A,DET)
 
C 	  THIS SUBROUTINE CALCULATES THE DETERMINANT
C 	  OF A 3 BY 3 MATRIX [A].
C**********************************************************************
      REAL*8  A(3,3), DET

      DET = A(1,1)*A(2,2)*A(3,3) 
     +      + A(1,2)*A(2,3)*A(3,1)
     +      + A(1,3)*A(2,1)*A(3,2)
     +      - A(3,1)*A(2,2)*A(1,3)
     +      - A(3,2)*A(2,3)*A(1,1)
     +      - A(3,3)*A(2,1)*A(1,2)

      RETURN
      END
      
C**********************************************************************
      SUBROUTINE ASSEMBLEDDSDDE(E,nu,D)
 
C	    THIS SUBROUTINE CALCULATES THE TRANSPOSE OF A 3 BY 3 
C	    MATRIX [A], AND PLACES THE RESULT IN ATRANS. 
C**********************************************************************

      REAL*8 E, nu, D(6,6), PF
      
      PARAMETER (ZERO = 0.D0)

      PF = E/((1+nu)*(1-2*nu))
      DO I=1,6
          DO J=1,6
              D(I,J) = ZERO
          ENDDO
      ENDDO
      
      D(1,1) = 1-nu
      D(2,2) = 1-nu
      D(3,3) = 1-nu
      D(1,2) = nu
      D(2,3) = nu
      D(3,1) = nu
      D(2,1) = nu
      D(3,2) = nu
      D(1,3) = nu
      D(4,4) = (1-2*nu)/2
      D(5,5) = (1-2*nu)/2
      D(6,6) = (1-2*nu)/2
      D = PF*D

c      PRINT*, '0 = ',ZERO
c      PRINT*, 'nu = ',nu
c      PRINT*, 'E = ',E
c      PRINT*, 'PF = ',PF
c      PRINT*, 'D = ',D

      RETURN
      END
      
C**********************************************************************
	SUBROUTINE LEFTDECOMP(F,R,V,E)

C	THIS SUBROUTINE PERFORMS THE LEFT POLAR DECOMPOSITION
C	[F] = [V][R] OF THE DEFORMATION GRADIENT [F] INTO
C	A ROTATION [R] AND THE LEFT STRETCH TENSOR [V].
C	THE EIGENVALUES AND EIGENVECTORS OF [V] AND
C	THE LOGARITHMIC STRAIN [E] = LN [V]
C	ARE ALSO RETURNED.
C**********************************************************************

	IMPLICIT REAL*8 (A-H,O-Z)
	REAL*8 F(3,3),FTRANS(3,3), B(3,3), OMEGA(3),
     +	          VEIGVAL(3),EIGVEC(3,3), EIGVECT(3,3), 
     +            V(3,3),E(3,3),VINV(3,3),R(3,3),TEMPM(3,3)
     
	COMMON/ERRORINFO/UMERROR     

C	F(3,3)	-- THE DEFORMATION GRADIENT MATRIX WHOSE
C		   POLAR DECOMPOSITION IS DESIRED.
C	DETF	-- THE DETRMINANT OF [F]; DETF > 0.
C	FTRANS(3,3)	-- THE TRANSPOSE OF [F].
C	R(3,3)	-- THE ROTATION MATRIX; [R]^T [R] = [I];
C		   OUTPUT.
C	V(3,3)	-- THE LEFT STRETCH TENSOR; SYMMETRIC
C		   AND POSITIVE DEFINITE; OUTPUT.
C	VINV(3,3)	-- THE INVERSE OF [V].
C	B(3,3)	-- THE LEFT CAUCHY-GREEN TENSOR = [V][V];
C		   SYMMETRIC AND POSITIVE DEFINITE.
C	OMEGA(3)-- THE SQUARES OF THE PRINCIPAL STRETCHES.
C 	VEIGVAL(3)	-- THE PRINCIPAL STRETCHES; OUTPUT.
C	EIGVEC(3,3)	-- MATRIX OF EIGENVECTORS OF [V];OUTPUT.
C	EIGVECT(3,3)    -- TRANSPOSE OF THE ABOVE.
C	E(3,3)	-- THE LOGARITHMIC STRAIN TENSOR, [E]=LN[V];
C		   OUTPUT.
C**********************************************************************

C	STORE THE IDENTITY MATRIX IN  [R], [V], AND [VINV]

	CALL ONEM(R)
	CALL ONEM(V)
	CALL ONEM(VINV)
     
C	STORE THE ZERO MATRIX IN [E]

	CALL ZEROM(E)

C      	CHECK IF THE DETERMINANT OF [F] IS GREATER THAN ZERO.
C	IF NOT, THEN PRINT DIAGNOSTIC AND STOP.

        CALL MDET(F,DETF)
c        IF (DETF .LE. 0.D0) THEN
c          WRITE(*,100)
c          UMERROR=5.
c          RETURN
c        ENDIF
         
C      	CALCULATE THE RIGHT CAUCHY GREEN STRAIN TENSOR [B]

        CALL  MTRANS(F,FTRANS)
        CALL  MPROD(F,FTRANS,B)
 
C	CALCULATE THE EIGENVALUES AND EIGENVECTORS OF [B]

	CALL SPECTRAL(B,OMEGA,EIGVEC)
          
C	CALCULATE THE PRINCIPAL VALUES OF [V] AND [E]

	VEIGVAL(1) = SQRT(OMEGA(1))
	VEIGVAL(2) = SQRT(OMEGA(2))
	VEIGVAL(3) = SQRT(OMEGA(3))
          
	V(1,1) = VEIGVAL(1)
	V(2,2) = VEIGVAL(2)
	V(3,3) = VEIGVAL(3)
          
	E(1,1) = LOG( VEIGVAL(1) )
	E(2,2) = LOG( VEIGVAL(2) )
	E(3,3) = LOG( VEIGVAL(3) )

C	TRANSFORM [V] AND [E] BACK TO REFERENCE FRAME

	CALL MTRANS(EIGVEC,EIGVECT)
      TEMPM = MATMUL(EIGVEC,V)
      V = MATMUL(TEMPM,EIGVECT)
      TEMPM = MATMUL(EIGVEC,E)
      E = MATMUL(TEMPM,EIGVECT)
C      
c	CALL MPROD(EIGVEC,V,TEMPM)
c	CALL MPROD(TEMPM,EIGVECT,V)
c	CALL MPROD(EIGVEC,E,TEMPM)
c	CALL MPROD(TEMPM,EIGVECT,E)

C	CALCULATE [VINV]

	CALL M3INV(V,VINV)
          
C	CALCULATE [R]

      CALL MPROD(VINV,F,R)

!100       FORMAT(5X,'--ERROR IN KINEMATICS-- THE DETERMINANT OF [F]',
!     +         ' IS NOT GREATER THAN 0')

	RETURN
      END

C**********************************************************************
	SUBROUTINE RIGHTDECOMP(F,R,U,E)

C	THIS SUBROUTINE PERFORMS THE RIGHT POLAR DECOMPOSITION
C	[F] = [R][U] OF THE DEFORMATION GRADIENT [F] INTO
C	A ROTATION [R] AND THE RIGHT  STRETCH TENSOR [U].
C	THE EIGENVALUES AND EIGENVECTORS OF [U] AND
C	THE LOGARITHMIC STRAIN [E] = LN [U]
C	ARE ALSO RETURNED.
C**********************************************************************

	IMPLICIT REAL*8 (A-H,O-Z)
	DIMENSION F(3,3),FTRANS(3,3), C(3,3), OMEGA(3),
     +	          UEIGVAL(3),EIGVEC(3,3), EIGVECT(3,3), 
     +            U(3,3),E(3,3),UINV(3,3),R(3,3),TEMPM(3,3)
     
	COMMON/ERRORINFO/UMERROR     

C	F(3,3)	-- THE DEFORMATION GRADIENT MATRIX WHOSE
C		   POLAR DECOMPOSITION IS DESIRED.
C	DETF	-- THE DETRMINANT OF [F]; DETF > 0.
C	FTRANS(3,3)	-- THE TRANSPOSE OF [F].
C	R(3,3)	-- THE ROTATION MATRIX; [R]^T [R] = [I];
C		   OUTPUT.
C	U(3,3)	-- THE RIGHT STRETCH TENSOR; SYMMETRIC
C		   AND POSITIVE DEFINITE; OUTPUT.
C	UINV(3,3)	-- THE INVERSE OF [U].
C	C(3,3)	-- THE RIGHT CAUCHY-GREEN TENSOR = [U][U];
C		   SYMMETRIC AND POSITIVE DEFINITE.
C	OMEGA(3)-- THE SQUARES OF THE PRINCIPAL STRETCHES.
C 	UEIGVAL(3)	-- THE PRINCIPAL STRETCHES; OUTPUT.
C	EIGVEC(3,3)	-- MATRIX OF EIGENVECTORS OF [U];OUTPUT.
C	EIGVECT(3,3)    -- TRANSPOSE OF THE ABOVE.
C	E(3,3)	-- THE LOGARITHMIC STRAIN TENSOR, [E]=LN[U];
C		   OUTPUT.
C**********************************************************************

C	STORE THE IDENTITY MATRIX IN  [R], [U], AND [UINV]

	CALL ONEM(R)
	CALL ONEM(U)
	CALL ONEM(UINV)

C	STORE THE ZERO MATRIX IN [E]

	CALL ZEROM(E)

C      	CHECK IF THE DETERMINANT OF [F] IS GREATER THAN ZERO.
C	IF NOT, THEN PRINT DIAGNOSTIC AND STOP.

        CALL MDET(F,DETF)
c        IF (DETF .LE. 0.D0) THEN
c          WRITE(*,100)
c          UMERROR=5.
c          RETURN
c        ENDIF

C      	CALCULATE THE RIGHT CAUCHY GREEN STRAIN TENSOR [C]

        CALL  MTRANS(F,FTRANS)
        CALL  MPROD(FTRANS,F,C)
 
C	CALCULATE THE EIGENVALUES AND EIGENVECTORS OF  [C]

	CALL SPECTRAL(C,OMEGA,EIGVEC)

C	CALCULATE THE PRINCIPAL VALUES OF [U] AND [E]

	UEIGVAL(1) = SQRT(OMEGA(1))
	UEIGVAL(2) = SQRT(OMEGA(2))
	UEIGVAL(3) = SQRT(OMEGA(3))

	U(1,1) = UEIGVAL(1)
	U(2,2) = UEIGVAL(2)
	U(3,3) = UEIGVAL(3)

	E(1,1) = LOG( UEIGVAL(1) )
	E(2,2) = LOG( UEIGVAL(2) )
	E(3,3) = LOG( UEIGVAL(3) )

C	CALCULATE THE COMPLETE TENSORS [U] AND [E]

	CALL MTRANS(EIGVEC,EIGVECT)
	CALL MPROD(EIGVEC,U,TEMPM)
	CALL MPROD(TEMPM,EIGVECT,U)
	CALL MPROD(EIGVEC,E,TEMPM)
	CALL MPROD(TEMPM,EIGVECT,E)

C	CALCULATE [UINV]

	CALL M3INV(U,UINV)

C	CALCULATE [R]

      CALL MPROD(F,UINV,R)
100     FORMAT(5X,'--ERROR IN KINEMATICS-- THE DETERMINANT OF [F]',
     +         ' IS NOT GREATER THAN 0')

	RETURN
      END

C**********************************************************************
C	THE FOLLOWING SUBROUTINES CALCULATE THE SPECTRAL
C	DECOMPOSITION OF A SYMMETRIC THREE BY THREE MATRIX
C**********************************************************************
	SUBROUTINE SPECTRAL(A,D,V)
C
C	THIS SUBROUTINE CALCULATES THE EIGENVALUES AND EIGENVECTORS OF
C	A SYMMETRIC 3 BY 3 MATRIX [A]. 
C
C	THE OUTPUT CONSISTS OF A VECTOR D CONTAINING THE THREE
C	EIGENVALUES IN ASCENDING ORDER, AND
C	A MATRIX [V] WHOSE COLUMNS CONTAIN THE CORRESPONDING
C	EIGENVECTORS.
C**********************************************************************

	IMPLICIT REAL*8 (A-H,O-Z)
	PARAMETER(NP=3)
	DIMENSION D(NP),V(NP,NP)
	DIMENSION A(3,3),E(NP,NP)

	DO 2 I = 1,3
          DO 1 J= 1,3
            E(I,J) = A(I,J)
1	  CONTINUE
2	CONTINUE

	CALL JACOBI(E,3,NP,D,V,NROT)
	CALL EIGSRT(D,V,3,NP)

	RETURN
      END

C**********************************************************************
	SUBROUTINE JACOBI(A,N,NP,D,V,NROT)

C	COMPUTES ALL EIGENVALUES AND EIGENVECTORS OF A REAL SYMMETRIC
C	MATRIX [A], WHICH IS OF SIZE N BY N, STORED IN A PHYSICAL 
C	NP BY BP ARRAY. ON OUTPUT, ELEMENTS OF [A] ABOVE THE DIAGONAL 
C	ARE DESTROYED, BUT THE DIAGONAL AND SUB-DIAGONAL ARE UNCHANGED
C	AND GIVE FULL INFORMATION ABOUT THE ORIGINAL SYMMETRIC MATRIX.
C	VECTOR D RETURNS THE EIGENVALUES OF [A] IN ITS FIRST N ELEMENTS.
C	[V] IS A MATRIX WITH THE SAME LOGICAL AND PHYSICAL DIMENSIONS AS
C	[A] WHOSE COLUMNS CONTAIN, ON OUTPUT, THE NORMALIZED
C	EIGENVECTORSOF [A]. NROT RETURNS THE NUMBER OF JACOBI ROTATIONS
C	WHICH WERE REQUIRED.

C	THIS SUBROUTINE IS TAKEN FROM "NUMERICAL RECIPES", PAGE 346.
C**********************************************************************

	IMPLICIT REAL*8 (A-H,O-Z)
	PARAMETER (NMAX =100)
	DIMENSION A(NP,NP),D(NP),V(NP,NP),B(NMAX),Z(NMAX)

C	INITIALIZE [V] TO THE IDENTITY MATRIX

	DO 12 IP = 1,N	
	  DO 11 IQ = 1,N
	    V(IP,IQ) = 0.D0
11        CONTINUE
          V(IP,IP) = 1.D0
12	CONTINUE

C	INITIALIZE [B] AND [D] TO THE DIAGONAL OF [A], AND Z TO ZERO.
C	THE VECTOR Z WILL ACCUMULATE TERMS OF THE FORM T*A_PQ AS
C	IN EQUATION (11.1.14)

	DO 13 IP = 1,N
	  B(IP) = A(IP,IP)
	  D(IP) = B(IP)
	  Z(IP) = 0.D0
13	CONTINUE
C
	NROT = 0
	DO 24 I = 1,50

C	SUM OFF-DIAGONAL ELEMENTS

          SM = 0.D0
          DO 15 IP = 1, N-1
            DO 14 IQ = IP + 1, N
	      SM = SM + DABS ( A(IP,IQ ))
14          CONTINUE
15        CONTINUE

C	IF SUM = 0., THEN RETURN. THIS IS THE NORMAL RETURN
C	WHICH RELIES ON QUADRATIC CONVERGENCE TO MACHINE 
C	UNDERFLOW.

          IF ( SM .EQ. 0.D0) RETURN
C
C	  IF ( SM .LT. 1.0D-15) RETURN

C	IN THE FIRST THREE SWEEPS CARRY OUT THE PQ ROTATION ONLY IF
C	|A_PQ| > TRESH, WHERE TRESH IS SOME THRESHOLD VALUE, 
C	SEE EQUATION (11.1.25). THEREAFTER TRESH = 0.

          IF ( I .LT. 4) THEN
            TRESH = 0.2D0*SM/N**2
          ELSE
            TRESH = 0.D0
          ENDIF
C
          DO 22 IP = 1, N-1
            DO 21 IQ = IP+1,N
              G = 100.D0*DABS(A(IP,IQ))

C	AFTER FOUR SWEEPS, SKIP THE ROTATION IF THE
C	OFF-DIAGONAL ELEMENT IS SMALL.

	      IF ((I .GT. 4) .AND. (DABS(D(IP))+G .EQ. DABS(D(IP)))
     +            .AND. ( DABS(D(IQ))+G .EQ. DABS(D(IQ)))) THEN
                A(IP,IQ) = 0.D0
              ELSE IF ( DABS(A(IP,IQ)) .GT. TRESH) THEN
                H = D(IQ) - D(IP)
                IF (DABS(H)+G .EQ. DABS(H)) THEN

C	T = 1./(2.*THETA), EQUATION(11.1.10)

	          T =A(IP,IQ)/H
	        ELSE
	          THETA = 0.5D0*H/A(IP,IQ)
	          T =1.D0/(DABS(THETA)+SQRT(1.D0+THETA**2))
	          IF (THETA .LT. 0.D0) T = -T
	        ENDIF
	        C = 1.D0/SQRT(1.D0 + T**2)
	        S = T*C
	        TAU = S/(1.D0 + C)
	        H = T*A(IP,IQ)
	        Z(IP) = Z(IP) - H
	        Z(IQ) = Z(IQ) + H
	        D(IP) = D(IP) - H
	        D(IQ) = D(IQ) + H
	        A(IP,IQ) = 0.D0

C	CASE OF ROTATIONS 1 <= J < P
				
	        DO 16 J = 1, IP-1
	          G = A(J,IP)
	          H = A(J,IQ)
	          A(J,IP) = G - S*(H + G*TAU)
	          A(J,IQ) = H + S*(G - H*TAU)
16	        CONTINUE

C	CASE OF ROTATIONS P < J < Q

	        DO 17 J = IP+1, IQ-1
	          G = A(IP,J)
	          H = A(J,IQ)
	          A(IP,J) = G - S*(H + G*TAU)
	          A(J,IQ) = H + S*(G - H*TAU)
17	        CONTINUE

C	CASE OF ROTATIONS Q < J <= N

	        DO 18 J = IQ+1, N
                  G = A(IP,J)
	          H = A(IQ,J)
	          A(IP,J) = G - S*(H + G*TAU)
	          A(IQ,J) = H + S*(G - H*TAU)
18	        CONTINUE
	        DO 19 J = 1,N
	          G = V(J,IP)
	          H = V(J,IQ)
	          V(J,IP) = G - S*(H + G*TAU)
	          V(J,IQ) = H + S*(G - H*TAU)
19	        CONTINUE
	        NROT = NROT + 1
              ENDIF
21	    CONTINUE
22	  CONTINUE

C	UPDATE D WITH THE SUM OF T*A_PQ, AND REINITIALIZE Z

	  DO 23 IP = 1, N
	    B(IP) = B(IP) + Z(IP)
	    D(IP) = B(IP)
	    Z(IP) = 0.D0
23	  CONTINUE
24	CONTINUE

C	IF THE ALGORITHM HAS REACHED THIS STAGE, THEN
C	THERE ARE TOO MANY SWEEPS, PRINT A DIAGNOSTIC
C	AND STOP.

	WRITE (*,'(/1X,A/)') '50 ITERS IN JACOBI SHOULD NEVER HAPPEN'

	RETURN
	END

C**********************************************************************
	SUBROUTINE EIGSRT(D,V,N,NP)

C	GIVEN THE EIGENVALUES [D] AND EIGENVECTORS [V] AS OUTPUT FROM
C	JACOBI, THIS ROUTINE SORTS THE EIGENVALUES INTO ASCENDING ORDER, 
C	AND REARRANGES THE COLUMNS OF [V] ACCORDINGLY.

C	THIS SUBROUTINE IS TAKEN FROM "NUMERICAL RECIPES", P. 348.
C**********************************************************************

	IMPLICIT REAL*8 (A-H,O-Z)
	DIMENSION D(NP),V(NP,NP)

	DO 13 I = 1,N-1
	  K = I
	  P = D(I)
	  DO 11 J = I+1,N
	    IF (D(J) .GE. P) THEN
	      K = J
	      P = D(J)
	    END IF
11	  CONTINUE
	  IF (K .NE. I) THEN
	    D(K) = D(I)
	    D(I) = P
	    DO 12 J = 1,N
	      P = V(J,I)
	      V(J,I) = V(J,K)
	      V(J,K) = P
12	    CONTINUE
  	  ENDIF
13	CONTINUE

	RETURN
      END
      
C**********************************************************************
	SUBROUTINE MPROD(A,B,C)
 
C 	THIS SUBROUTINE MULTIPLIES TWO 3 BY 3 MATRICES [A] AND [B],
C 	AND PLACE THEIR PRODUCT IN MATRIX [C]. 
C**********************************************************************

	REAL*8 A(3,3),B(3,3),C(3,3)

	DO 2 I = 1, 3
	  DO 2 J = 1, 3
	    C(I,J) = 0.D0
	    DO 1 K = 1, 3
	      C(I,J) = C(I,J) + A(I,K) * B(K,J)                       
1	    CONTINUE
2	CONTINUE
C
	RETURN
      END
      
C**********************************************************************
      	SUBROUTINE ZEROM(A)
C
C	THIS SUBROUTINE SETS ALL ENTRIES OF A 3 BY 3 MATRIX TO 0.D0.
C**********************************************************************

        REAL*8 A(3,3)

	DO 1 I=1,3
	  DO 1 J=1,3
	    A(I,J) = 0.D0
1	CONTINUE
C	
	RETURN
      END

C**********************************************************************
	SUBROUTINE ONEM(A)

C	THIS SUBROUTINE STORES THE IDENTITY MATRIX IN THE 
C	3 BY 3 MATRIX [A]
C**********************************************************************

        REAL*8 A(3,3)
        DATA ZERO/0.D0/
        DATA ONE/1.D0/

	DO 1 I=1,3
	  DO 1 J=1,3
	    IF (I .EQ. J) THEN
              A(I,J) = 1.0
            ELSE
              A(I,J) = 0.0
            ENDIF
1       CONTINUE

	RETURN
	END
C**********************************************************************
      SUBROUTINE COMPUTEALPHANT(theta_tau,Tg,Tv,alpha_NT)
      
C	THIS SUBROUTINE COMPUTES THE COEFFICIENT OF THERMAL EXPANSION (CTE)
C     AS A FUNCTION OF TEMPERTURE (PIECEWISE) 
C**********************************************************************

      REAL*8 TEMP, alpha_NT, Tg, Tv,theta_tau
      ! CTE for 0% cataly, please update for 5% catalyst
      IF (theta_tau<Tg) THEN	  
          alpha_NT = 0.0001761D0
      ELSEIF (theta_tau>=Tg) THEN
          alpha_NT = 0.0002486D0     
      ENDIF

	RETURN
      END
      
C**********************************************************************
      SUBROUTINE COMPUTEMU_E1(theta_tau,E1)

C	THIS SUBROUTINE COMPUTES MU1 (ELASTIC LEG) AS A FUNCTION OF TEMPERATURE
C     BASED ON THE EXPERIMENTAL DATA
C**********************************************************************
      PARAMETER(MIDTEMP = 40.D0, T1 = 39.9999D0, T2 = 40.9639D0,
     +    Y1 = 1568.88D0, Y2 = 1311D0)

      REAL*8 TEMP, mu1, A1, A2, A3, A4, K1, K2, K3, K4, T, theta_tau, E1
      
      A1 = 1.635E11
      A2 = -0.9486
      A3 = 25.13
	  A4 = -0.009538
      K1 = -0.01
      K2 = 16.1
      K3 = 0.06719
      K4 = 10.44	  
	  IF (theta_tau < 50D0) THEN
	      E1 = A1*EXP(A2*theta_tau) + A3*EXP(A4*theta_tau)		  
	  ELSEIF (theta_tau >= 50D0 .AND. theta_tau <= 70D0) THEN
          E1 = K1*theta_tau+K2 	  
	  ELSEIF (theta_tau>70.D0) THEN
          E1 = K3*theta_tau+K4 	  
	  ENDIF       
	RETURN
      END
      
C**********************************************************************
      SUBROUTINE COMPUTEMU_E2(theta_tau,E2)

C	THIS SUBROUTINE COMPUTES MU2 (VISCOELASTIC LEG) AS A FUNCTION OF TEMPERATURE
C     BASED ON THE EXPERIMENTAL DATA
C**********************************************************************
      PARAMETER(MIDTEMP = 47D0)

      REAL*8 TEMP, mu2, TERM1, TERM2, TERM3, TERM4, TERM5, W, A0, 
     +    A1, A2, A3, A4, A5, A6, B1, B2, B3, B4, B5, B6, 
     +    K1, K2, K3, K4, T, Y1, Y2, T1, T2,theta_tau, E2
      

      A1 = 6.315E9
      A2 = -4.705
      A3 = -12.18
      K1 = -0.3143
      K2 = 342.5

	  IF (theta_tau < 35D0) THEN
          E2 = K1*theta_tau+K2	
	  ELSEIF (theta_tau >= 35D0 .AND. theta_tau <= 70D0) THEN
          E2 = A1*theta_tau**A2+A3  	  
	  ELSEIF (theta_tau>70.D0) THEN
          E2 = 1  	  
	  ENDIF 
        
	RETURN
      END
C**********************************************************************
      SUBROUTINE COMPUTEMU_E3(theta_tau,E3)

C	THIS SUBROUTINE COMPUTES MU2 (VISCOELASTIC LEG) AS A FUNCTION OF TEMPERATURE
C     BASED ON THE EXPERIMENTAL DATA
C**********************************************************************
      PARAMETER(MIDTEMP = 47D0)

      REAL*8 TEMP, mu1, A1, B1, C1, A2, B2, C2, A3, B3, C3, A4, B4, C4,
     +    TERM1, TERM2, TERM3, TERM4,theta_tau, E3     
     
      A1 = -1.885
      A2 = 127.3
      A3 = -1603
	  B1 = -2.582
	  B2 = 181.7	  
	  IF (theta_tau < 50D0) THEN
          E3 = A1*theta_tau**2+A2*theta_tau+A3
	  ELSEIF (theta_tau >= 50D0 .AND. theta_tau <= 70D0) THEN
		  E3 = B1*theta_tau+B2		  		  
	  ELSEIF (theta_tau>70.D0) THEN
          E3 = 1.D0  	  
	  ENDIF 

	RETURN
      END	  

C**********************************************************************
      SUBROUTINE COMPUTEMU_TAU(theta_tau,tau_T)

C	THIS SUBROUTINE COMPUTES MU1 (ELASTIC LEG) AS A FUNCTION OF TEMPERATURE
C     BASED ON THE EXPERIMENTAL DATA
C**********************************************************************
      PARAMETER(MIDTEMP = 40.D0, T1 = 39.9999D0, T2 = 40.9639D0,
     +    Y1 = 1568.88D0, Y2 = 1311D0)

      REAL*8 TEMP, mu1, A1, A2, A3, A4, K1, K2, K3, K4, T, C1, tau_T, theta_tau
          
      A1 = 2.529E10
      A2 = -0.9406D0
      A3 = 20D0
      A4 = 5.713E-12
	  tau_T = A1*EXP(A2*theta_tau) + A3*EXP(A4*theta_tau)+20	  
	RETURN
      END
C**********************************************************************
      SUBROUTINE COMPUTEMU_TAU2(theta_tau,tau2_T)

C	THIS SUBROUTINE COMPUTES MU1 (ELASTIC LEG) AS A FUNCTION OF TEMPERATURE
C     BASED ON THE EXPERIMENTAL DATA
C**********************************************************************
      PARAMETER(MIDTEMP = 40.D0, T1 = 39.9999D0, T2 = 40.9639D0,
     +    Y1 = 1568.88D0, Y2 = 1311D0)

      REAL*8 TEMP, mu1, A1, B1, A3, A4, K1, K2, K3, K4, T, C1, tau2_T, theta_tau
      
      A1 = 1.395E8
      A2 = -2.9
      A3 = -1639	  
	  IF (theta_tau < 50D0) THEN
          tau2_T = A1*theta_tau**A2+A3			  
	  ELSE
          tau2_T = 10D0  	  
	  ENDIF 	  
	RETURN
      END	  


C**********************************************************************
      SUBROUTINE COMPUTEMU_TAU3(theta_tau,tau3_T)

C	THIS SUBROUTINE COMPUTES MU1 (ELASTIC LEG) AS A FUNCTION OF TEMPERATURE
C     BASED ON THE EXPERIMENTAL DATA
C**********************************************************************
      PARAMETER(MIDTEMP = 40.D0, T1 = 39.9999D0, T2 = 40.9639D0,
     +    Y1 = 1568.88D0, Y2 = 1311D0)

      REAL*8 TEMP, mu1, A1,A2, B1, A3, A4, K1, K2, K3, K4, T, C1, tau3_T, theta_tau
      
      A1 = 2.184E12
      A2 = -6.594D0
      A3 = -3.664	  
	  IF (theta_tau < 50D0) THEN
          tau3_T = A1*theta_tau**A2+A3
	  ELSE
          tau3_T = 10D0  	  
	  ENDIF 	  
	RETURN
      END	  