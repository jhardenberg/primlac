! version 22
! Fixed issues which prevented compilation when ECOSYSTEM is used 

! Remember to fix immbound bottom condition also for scalars!!!
! -------------------------------------------------------------
!
!  lac.F -- Primitive equations model for shallow lakes
!  by Jost von Hardenberg - ISAC-CNR - 3/2009 - 2/2010
!     Maurice Rossi       - UPMC-CNRS    
!
!  Uses the Wu and Tsanis 1995 parametrization of vertical eddy viscosity
!
!  Arakawa C staggered grid                       
!  Vertical sigma coordinates
!  Second order interpolation/differencing 
!  Leapfrog + Adams-Moulton P-C timestepping                             
!
!  This code contains some very long lines, 
!  please use the 132 columns option of your compiler! 

!   gfortran -ffixed-line-length-none lac16.F 
!   ifort lac16.F -132
!
!
!  When the horizontal resolution NY is set to 1               
!  the preprocessor automatically removes not relevant lines 
!
!  BCs: no-stress+impermeability on the sides,
!       no-slip at the bottom and wind stress at the top 
!
! -----------------------------------------------------------------

#define TIMESPLIT
#define IMPLICITDIFF 

!#define DEBUGCONSTSCAL

! This option sets slc at half-points between the zlws
!#define DEBUGTANPAR

! USE IMMERSED BOUNDARIES
!#define IMMBOUND

! newlago8 : version having immersed boundaries and not allowing yet for lateral no-slip BCs 
! newlago10 : adding implicit vertical diffusion
! newlago11  : more detailed parabolic nabla2 at the bottom (removed  in newlago16)
! newlago12 : diffusion is now completely implicit :)
! newlago13  : cleaned up version
! newlago14 :  Optional split time-stepping (External(Barotropic)+Internal(Baroclinic) modes)
! newlago15 : working version with both timesplitting and implicit diffusion
! newlago16 : last version with INSIMPSON ETC. Fixed issues with IMMBOUND. Removed parabolic deriv at bottom.
! newlago17 : Fixed issues with computation of taubu+IMMBOUND; Removed INTSIMPSON etc
! newlago18 : Fixed yet another immbound bug!
! newlago19 : All scalars using Patankard....
! newlago20 : fixed LMIXING for the Immersed Boundary cqs 

! Notes: dutop used also in barotropic part if INTSQUARE is not defined 

! %%%%%%%%%%%%DEFINES AND MACROS %%%%%%%%%%%%%%%%%%%%%

! %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
!  MODEL DESCRIPTION 
! %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

! Uncomment #defines to activate the flags
! For parameters Check params.dat and paramseco.dat

! Radiation (acts on temperature field and/or light for ecosystem) ***
!#define RADIATION

! Ecosystem
#define ECOSYSTEM

! Second Phyto species
!#define PHYTO2

!  Wind 
! 
! WIND  is constant if VARIABLEWIND is undefined (using params.dat for the constant values)
! WIND  varies  if VARIABLEWIND  is defined (using params.dat for amplitude and then check the subroutine updatewind)

!#define VARIABLEWIND
! Example: regular morning wind
!#define WINDCYCLE(ntime,ddt)  sin(ntime*ddt/86400*pi+pi/8)**16 
! This is a wind with maximum at 1d1h and lasting 1d14h
!#define WINDCYCLE(ntime,ddt) (0.1+exp(-(ntime*ddt/86400-(1.+1./24))**2/(1.+14./24.)**2))  

! Temperature 
c-#define TEMPERATURE

! Initialize with a fixed temperature profile (tanh((xz-0.666)*10)+10)
c-#define TEMPTANH

! Diurnal cycle for light fluxes (the alternative is CONSTANT if not defined)
!#define DIURNALCYCLE
! Function used to modulate light (w/ DIURNALCYCLE) when RADIATION is not used (only for ecosystem)
!#define LIGHTCYCLE(ntime,ddt) (1+tanh(6*sin(2*pi*ntime*ddt/86400.d0-pi/2)))/2  

! Define this to fix the velocity profile to the Tsanis et al. solution 
! (used for tests with a 1D column model)
!#define FIXED_VEL_PROF

!#define NUTCYCLE(ntime,ddt) (0.+exp(-(ntime*ddt/86400-5)**2/(1.)**2))  

! Add small island in the center
!#define MODBATHY
!#define ISLANDRADIUS 2
!#define ISLANDBANKS 5

! -----------------------------------------------
! *** Vertical Diffusion ***
! Choose one. Priorities: LMIXING>KVPAR>KVERT
! -----------------------------------------------
! Use mixing length model
! This case uses parameter: kvmin
! KARMAN is the von Karman constant used for LMIXING
c-#define LMIXING 
#define KARMAN 0.4
c-#define KARMAN 0.

! If density effect :Parametrization based on gradient Richardson number (Pacanowski et al)
! used only for LMIXING case with TEMPERATURE
#define PACANOWSKI
#define ALPHAPAC 5
#define NPAC -2

! Parabolic vertical diffusion (Tsanis) 
! this case uses parameters:   lambda,  zsh,  zbh . see Chao et al.
!#define KVPAR

#define ZSH 2.2e-4
#define ZBH 2e-3 
#define LAMBDA 0.4

! Constant vertical diffusivity 
#define KVERT kvmin

!Convective adjustment. Relevand only when TEMPERATURE is active.
! Defines how often (in steps) to perform it
#define CONVADJ 600

! ---------------------------------------------------

! %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
! SPATIAL RESOLUTION OF THE GRID
! %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
!NNT   rescales NX, NY, NZ   
!NNT   also rescales the timestep defined in params.dat

#define NFASTSTEPS 100

#define NNT 1
#define NX  256*NNT
#define NZ  32*NNT
!Use  Define NY = 1 to have a 2D code (vertical slice) If NOT 3D CODE
!#define NY  40*NNT
#define NY 1
 
!Vertical irregular Tanh grid parameter. If you want a regular grid UNdefine the following (0.2 default)
#define TANPAR 0.4

! If to save scalar statistics
#define CHECKSCAL 

!SEDIMENT TRANSPORT:Flux Corrected Transport scheme to guarantee positive definite scalar fields
!#define FCT_SCAL


! ---------------------------------------------------
 
!%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
! TIME and position  active with variablewind and Radiation
!%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
 
! HOUR, First of time step 
! Used both with RADIATION or simply for light modulation for ECOSYSTEM
! or simply for the WIND modulation
#define HOURSTART 0

! DAY   (Meaningful ONLY WHEN Radiation is active) here July 1st  
#define STARTDAY 181

! LAKE POSITION (Lon, Lat, timezone)    (Meaningful ONLY WHEN Radiation is active)
! Viverone
#define LONGITUDE 8.03d0 
#define LATITUDE 45.42d0
#define LTM 1

! Air pressure at lake position (here 2000 m)
!#define PAIR 760
#define PAIR 1000

! %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
! PHYSICS CONSTANTS  :  MOST  model parameters are   read from params.dat
! %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

! %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
!PHYSICS CONSTANTS  :  wind Density  
! %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
! Gravity
#define GG 9.81

!Temperature, Density 
!Values at 20 C 
#define RHOAIR 1.2041
#define RHOWAT 998.2071
!Equation of state for freshwater between 0 and 20 °C
#define RHOEOS(T)  (999.84298+1e-3*(65.4891*T-8.56272*T**2+0.059385*T**3))

!SPECIFIC HEAT  
#define CAIR 1005
#define CWAT 4193

! Wind drag 1e-3 from Chao et al 2008   2.6e-3 from ???
#define Cd 1.0e-3 
 
!%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
! PHYSICS CONSTANTS  :Radiation Constants 
!%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
! 0.218 J = 1 muMol light
#define MUMOL2JOULE 0.218

!Solar constant w/m^2
#define PHI0 1366

! Sky neatness (no pollution) and cloud cover
#define AASKY 0.77
#define BBSKY 1.22
#define AACONVEC 0.00117
#define BBCONVEC 0.00115
#define CLOUDCOV 0.1
#define ALBEDO 0.06
#define EPSWAT 0.97
#define STEFAN 5.67e-8
#define CEVAP 1.3e-3

!%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
! Derived quantities
!%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

#define WINDRAD windangle*6.28318530717959/360.

#define TAUS   (Cd*RHOAIR*UWIND*UWIND)
! Used for Tsanis
#define USTARS  sqrt(TAUS/RHOWAT)

!%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%% 
!  ------------  No need to touch anything beyond this point  ---------- 
!%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%% 
#define double real*8

#define NNX (NX+2)
#define NNZ (NZ+1) 
#if (NY!=1) 
#define c2d          
! GENERATES A3D CODE by activating all lines with c2d in front 
#define NYMIN 0
#define NYMAX NY+1
#define NNY (NY+2)
#else
#define NYMIN 1 
#define NYMAX 1
#define NNY 1
#endif

#define RANGE 0:(NX+1),NYMIN:NYMAX,0:NZ
#define RANGE2D 0:(NX+1),NYMIN:NYMAX
! ------------------------------
!Definition of points
#define HGRID 0
#define UGRID 1
#define VGRID 2
#define OGRID 3

!%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
! SCALAR FIELDS
!%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
!Total Number of scalars fields and what they are
#ifdef ECOSYSTEM
#ifndef PHYTO2
#define NSCAL 6
#else
#define NSCAL 7
#endif
#else
#define NSCAL 2
#endif

!Field identifiers
#define SCAL_NEUTRAL  1
#define SCAL_SEDIMENT 2
#define SCAL_DETRITUS 3
#define SCAL_NUTRIENT 4
#define SCAL_ZOO      5
#define SCAL_PHYTO1   6
#define SCAL_PHYTO2   7

#define SED(i,j,k) c(i,j,k,SCAL_SEDIMENT)

#ifdef ECOSYSTEM
#define DET(i,j,k)  c(i,j,k,SCAL_DETRITUS)
#define NUT(i,j,k)  c(i,j,k,SCAL_NUTRIENT)
#define ZOO(i,j,k)  c(i,j,k,SCAL_ZOO)
#define PHY1(i,j,k) c(i,j,k,SCAL_PHYTO1)
#define PHY2(i,j,k) c(i,j,k,SCAL_PHYTO2)

#endif

!%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
!Various combinations of fields used calling functions
!%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

#define FIELDSUVH u(1,1,1,now),v(1,1,1,now),h(1,1,now)
#define FIELDSUV u(1,1,1,now),v(1,1,1,now)

#ifdef TEMPERATURE

#define TENDFIELDS2 utold(0,NYMIN,0,nt2),vtold(0,NYMIN,0,nt2),htold(0,NYMIN,nt2),ctold(0,NYMIN,0,1,nt2),ttold(0,NYMIN,0,nt2)
#define TENDFIELDS1 utold(0,NYMIN,0,nt1),vtold(0,NYMIN,0,nt1),htold(0,NYMIN,nt1),ctold(0,NYMIN,0,1,nt1),ttold(0,NYMIN,0,nt1)
#define TENDFIELDS ut,vt,ht,ct,tt

#define FIELDLIST u,v,h,c,t
#define FIELDS u(1,1,1,now),v(1,1,1,now),h(1,1,now),c(1,1,1,1,now),t(1,1,1,now)
#define OLDFLD u(1,1,1,old),v(1,1,1,old),h(1,1,old),c(1,1,1,1,old),t(1,1,1,old)
#define EXTFIELDS u(1,1,1,now),v(1,1,1,now),w(1,1,1),h(1,1,now),c(1,1,1,1,now),t(1,1,1,now)

#else

#define TENDFIELDS2 utold(0,NYMIN,0,nt2),vtold(0,NYMIN,0,nt2),htold(0,NYMIN,nt2),ctold(0,NYMIN,0,1,nt2)
#define TENDFIELDS1 utold(0,NYMIN,0,nt1),vtold(0,NYMIN,0,nt1),htold(0,NYMIN,nt1),ctold(0,NYMIN,0,1,nt1)
#define TENDFIELDS ut,vt,ht,ct

#define FIELDLIST u,v,h,c
#define FIELDS u(1,1,1,now),v(1,1,1,now),h(1,1,now),c(1,1,1,1,now)
#define OLDFLD u(1,1,1,old),v(1,1,1,old),h(1,1,old),c(1,1,1,1,old)
#define EXTFIELDS u(1,1,1,now),v(1,1,1,now),w(1,1,1),h(1,1,now),c(1,1,1,1,now)

#endif

#ifdef TIMESPLIT 

#define TENDFIELDS2FAST umtold(0,NYMIN,nt2fast),vmtold(0,NYMIN,nt2fast),htold(0,NYMIN,nt2fast)
#define TENDFIELDS1FAST umtold(0,NYMIN,nt1fast),vmtold(0,NYMIN,nt1fast),htold(0,NYMIN,nt1fast)
#define TENDFIELDSFAST umt,vmt,ht

#define FIELDSFAST um(1,1,nowfast),vm(1,1,nowfast),hm(1,1,nowfast)
#define OLDFLDFAST um(1,1,oldfast),vm(1,1,oldfast),hm(1,1,oldfast)

#endif

!%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
! MACROS TO FIX HORIZONTAL Boundary Conditions
!%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

#define FIXBC2D(F,GTYPE) where(llm(:,:,GTYPE)==0); F(:,:)=0.; end where 
!#define FIXBC4D(F,GTYPE) do L=1,NZ; do K=1,NZ; where(llm(:,:,GTYPE)==0); F(:,:,K,L)=0.; end where; enddo; enddo ! NOT USED

!#define FIXBC3D(F,GTYPE) do K=1,NZ; where(llm(:,:,GTYPE)==0); F(:,:,K)=0.; end where; enddo
#define FIXBC3D(F,GTYPE) where(llm3d(:,:,:,GTYPE)==0); F(:,:,:)=0.; end where

!%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
!  Second order Horizontal  differencing operators on staggered grid 
!  cffx = 1/dx and cffy = 1/dy must be defined.  
!  All these operators are defined at (I+1/2,J,K), (I,J+1/2,K)   for points H(I,J,K)

!  E.g.: DX2 is the (horizontal)differencing operator in the x direction for 2D fields
!        DX3 is the same operator for 3D  fields   
!   	 DX4 is the same for 4D fields (the scalars)
!%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

#define DX2(F,I,J)     (CFFX*(F(I+1,J)-F(I,J)))
#define DY2(F,I,J)     (CFFY*(F(I,J+1)-F(I,J)))

#define DX3(F,I,J,K)   (CFFX*(F(I+1,J,K)-F(I,J,K)))
#define DY3(F,I,J,K)   (CFFY*(F(I,J+1,K)-F(I,J,K)))

#define DX4(F,I,J,K,L) (CFFX*(F(I+1,J,K,L)-F(I,J,K,L)))
#define DY4(F,I,J,K,L) (CFFY*(F(I,J+1,K,L)-F(I,J,K,L)))
 

!%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
!   IX2 IX3 IX4 are Horizontal  interpolation operators
!   All these operators are defined at (I+1/2,J,K), (I,J+1/2,K)  for point H(I,J,K)
!%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

#define IX2(F,I,J)     (0.5*(F(I+1,J)+F(I,J)))
#define IX3(F,I,J,K)   (0.5*(F(I+1,J,K)+F(I,J,K)))
#define IX4(F,I,J,K,L) (0.5*(F(I+1,J,K,L)+F(I,J,K,L)))

#define IY2(F,I,J)     (0.5*(F(I,J+1)+F(I,J)))
#define IY3(F,I,J,K)   (0.5*(F(I,J+1,K)+F(I,J,K)))
#define IY4(F,I,J,K,L) (0.5*(F(I,J+1,K,L)+F(I,J,K,L)))


!%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
! Block  horizontal interpolation and  differencing macros
! NOFF is an offset 
!%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

#define BLOCKDX2(F,FD,NOFF)     FD(-(NOFF):(NX-(NOFF)),:)=(CFFX*(F(1:(NX+1),:)-F(0:NX,:)))
#define BLOCKDY2(F,FD,NOFF)     FD(:,-(NOFF):(NY-(NOFF)))=(CFFY*(F(:,1:(NY+1))-F(:,0:NY)))
#define BLOCKIX2(F,FD,NOFF)     FD(-(NOFF):(NX-(NOFF)),:)=(0.5*(F(1:(NX+1),:)+F(0:NX,:)))
#define BLOCKIY2(F,FD,NOFF)     FD(:,-(NOFF):(NY-(NOFF)))=(0.5*(F(:,1:(NY+1))+F(:,0:NY)))


#define BLOCKDX3(F,FD,NOFF)     FD(-(NOFF):(NX-(NOFF)),:,:)=(CFFX*(F(1:(NX+1),:,:)-F(0:NX,:,:)))
#define BLOCKDY3(F,FD,NOFF)     FD(:,-(NOFF):(NY-(NOFF)),:)=(CFFY*(F(:,1:(NY+1),:)-F(:,0:NY,:)))
#define BLOCKIX3(F,FD,NOFF)     FD(-(NOFF):(NX-(NOFF)),:,:)=(0.5*(F(1:(NX+1),:,:)+F(0:NX,:,:)))
#define BLOCKIY3(F,FD,NOFF)     FD(:,-(NOFF):(NY-(NOFF)),:)=(0.5*(F(:,1:(NY+1),:)+F(:,0:NY,:)))


#define BLOCKDX4(F,FD,NOFF,L)   FD(-(NOFF):(NX-(NOFF)),:,:)=(CFFX*(F(1:(NX+1),:,:,L)-F(0:NX,:,:,L)))
#define BLOCKDY4(F,FD,NOFF,L)   FD(:,-(NOFF):(NY-(NOFF)),:)=(CFFY*(F(:,1:(NY+1),:,L)-F(:,0:NY,:,L)))
#define BLOCKIX4(F,FD,NOFF,L)   FD(-(NOFF):(NX-(NOFF)),:,:)=(0.5*(F(1:(NX+1),:,:,L)+F(0:NX,:,:,L)))
#define BLOCKIY4(F,FD,NOFF,L)   FD(:,-(NOFF):(NY-(NOFF)),:)=(0.5*(F(:,1:(NY+1),:,L)+F(:,0:NY,:,L)))


!%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
!Vertical Second order  differencing operators on staggered grid in   sigma coordinates.
!cffsigma   must be defined in  sigma 

!All these operators are defined at  (I,J,K+1/2)   for points C(I,J,K)
! Hence on plane of points w
!        DZC3 is the differencing vertical operator in  z direction using the points located on   plane of  point C
!        IZC3 is the  interpolation vertical operator in  z direction starting at the plane of  point C

!All these operators are defined at  (I,J,K+1/2)   for points W(I,J,K)
! Hence on plane of points C

!        DZW3 is the differencing vertical operator in  z direction starting at the plane of points W
!        IZW3 is the interpolation vertical operator in z direction starting at the plane of points W

!   	 DZC4 IZC4 are the same for 4D fields  
!   	 DZW4 IZW4 are the same for 4D fields  
!%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
        
#ifndef TANPAR
! REGULAR GRID

#define DZC3(F,I,J,K)   (cffsigma*(F(I,J,K+1)-F(I,J,K)))
#define IZC3(F,I,J,K)   (0.5*(F(I,J,K+1)+F(I,J,K)))
#define DZW3(F,I,J,K)   (cffsigma*(F(I,J,K+1)-F(I,J,K)))
#define IZW3(F,I,J,K)   (0.5*(F(I,J,K+1)+F(I,J,K)))
#define DZC4(F,I,J,K,L) (cffsigma*(F(I,J,K+1,L)-F(I,J,K,L)))
#define IZC4(F,I,J,K,L) (0.5*(F(I,J,K+1,L)+F(I,J,K,L)))
#define DZW4(F,I,J,K,L) (cffsigma*(F(I,J,K+1,L)-F(I,J,K,L)))
#define IZW4(F,I,J,K,L) (0.5*(F(I,J,K+1,L)+F(I,J,K,L)))

#else
! IRREGULAR GRID

#define DZC3(F,I,J,K)   (cfDZC(K,2)*F(I,J,K+1)+cfDZC(K,1)*F(I,J,K))
#define IZC3(F,I,J,K)   (cfIZC(K,2)*F(I,J,K+1)+cfIZC(K,1)*F(I,J,K))

#define DZW3(F,I,J,K)   (cfDZW(K,2)*F(I,J,K+1)+cfDZW(K,1)*F(I,J,K))
#define IZW3(F,I,J,K)   (cfIZW(K,2)*F(I,J,K+1)+cfIZW(K,1)*F(I,J,K))

#define DZC4(F,I,J,K,L) (cfDZC(K,2)*F(I,J,K+1,L)+cfDZC(K,1)*F(I,J,K,L))
#define IZC4(F,I,J,K,L) (cfIZC(K,2)*F(I,J,K+1,L)+cfIZC(K,1)*F(I,J,K,L))

#define DZW4(F,I,J,K,L) (cfDZW(K,2)*F(I,J,K+1,L)+cfDZW(K,1)*F(I,J,K,L))
#define IZW4(F,I,J,K,L) (cfIZW(K,2)*F(I,J,K+1,L)+cfIZW(K,1)*F(I,J,K,L))

#endif


!%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
! Block VERTICAL interpolation and    differencing macros
! NOFF is an offset 
!%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
 
#ifndef TANPAR

#define BLOCKDZC3(F,FD,NOFF)      do K=0,NZ-1; FD(:,:,(K-(NOFF)))=(cffsigma*(F(:,:,K+1)-F(:,:,K))); enddo
#define BLOCKIZC3(F,FD,NOFF)      do K=0,NZ-1; FD(:,:,(K-(NOFF)))=(0.5*(F(:,:,K+1)+F(:,:,K))); enddo

#define BLOCKDZW3(F,FD,NOFF)      do K=0,NZ-1; FD(:,:,(K-(NOFF)))=(cffsigma*(F(:,:,K+1)-F(:,:,K))); enddo
#define BLOCKIZW3(F,FD,NOFF)      do K=0,NZ-1; FD(:,:,(K-(NOFF)))=(0.5*(F(:,:,K+1)+F(:,:,K))); enddo


#define BLOCKDZC4(F,FD,NOFF,L)    do K=0,NZ-1; FD(:,:,(K-(NOFF)))=(cffsigma*(F(:,:,K+1,L)-F(:,:,K,L))); enddo
#define BLOCKIZC4(F,FD,NOFF,L)    do K=0,NZ-1; FD(:,:,(K-(NOFF)))=(0.5*(F(:,:,K+1,L)+F(:,:,K,L))); enddo

#define BLOCKDZW4(F,FD,NOFF,L)    do K=0,NZ-1; FD(:,:,(K-(NOFF)))=(cffsigma*(F(:,:,K+1,L)-F(:,:,K,L))); enddo
#define BLOCKIZW4(F,FD,NOFF,L)    do K=0,NZ-1; FD(:,:,(K-(NOFF)))=(0.5*(F(:,:,K+1,L)+F(:,:,K,L))); enddo

#else

#define BLOCKDZC3(F,FD,NOFF)     do K=0,NZ-1; FD(:,:,(K-(NOFF)))=(cfDZC(K,2)*F(:,:,K+1)+cfDZC(K,1)*F(:,:,K)); enddo
#define BLOCKIZC3(F,FD,NOFF)     do K=0,NZ-1; FD(:,:,(K-(NOFF)))=(cfIZC(K,2)*F(:,:,K+1)+cfIZC(K,1)*F(:,:,K)); enddo

#define BLOCKDZW3(F,FD,NOFF)     do K=0,NZ-1; FD(:,:,(K-(NOFF)))=(cfDZW(K,2)*F(:,:,K+1)+cfDZW(K,1)*F(:,:,K)); enddo
#define BLOCKIZW3(F,FD,NOFF)     do K=0,NZ-1; FD(:,:,(K-(NOFF)))=(cfIZW(K,2)*F(:,:,K+1)+cfIZW(K,1)*F(:,:,K)); enddo


#define BLOCKDZC4(F,FD,NOFF,L)   do K=0,NZ-1; FD(:,:,(K-(NOFF)))=(cfDZC(K,2)*F(:,:,K+1,L)+cfDZC(K,1)*F(:,:,K,L)); enddo
#define BLOCKIZC4(F,FD,NOFF,L)   do K=0,NZ-1; FD(:,:,(K-(NOFF)))=(cfIZC(K,2)*F(:,:,K+1,L)+cfIZC(K,1)*F(:,:,K,L)); enddo

#define BLOCKDZW4(F,FD,NOFF,L)   do K=0,NZ-1; FD(:,:,(K-(NOFF)))=(cfDZW(K,2)*F(:,:,K+1,L)+cfDZW(K,1)*F(:,:,K,L)); enddo
#define BLOCKIZW4(F,FD,NOFF,L)   do K=0,NZ-1; FD(:,:,(K-(NOFF)))=(cfIZW(K,2)*F(:,:,K+1,L)+cfIZW(K,1)*F(:,:,K,L)); enddo

#endif
!%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%% END OF DEFINES AND MACROS  %%%%%%%%%%% 



 
!%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%% 
!  -------------------Main --------------------------
!%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
      program lac

      implicit none
#include "commonlac.h"

      integer fout,i,j,k
      integer maxtime,mintime,mstar
      integer now,old,ntemp,nt1,nt2,nt1fast,nt2fast,nowfast,oldfast
	  
! Main output file
      parameter(fout=20)
! ------------------------------
! FIELDS 
      double u(NNX,NNY,NNZ,2)
      double v(NNX,NNY,NNZ,2)
      double w(NNX,NNY,NNZ)
      double h(NNX,NNY,2)
      double c(NNX,NNY,NNZ,NSCAL,2)
#ifdef TEMPERATURE
      double t(NNX,NNY,NNZ,2)
#endif

#ifdef TIMESPLIT

! ----------------------------
! Stuff for time-splitting
	  double aw(NFASTSTEPS*2)
	  double hm(NNX,NNY,2)
	  double hmean(NNX,NNY),umean(NNX,NNY),vmean(NNX,NNY)
	  double um(NNX,NNY,2)
	  double vm(NNX,NNY,2)
	  double um3d(NNX,NNY)
	  double vm3d(NNX,NNY)
	  double bm3du(NNX,NNY),bm3dv(NNX,NNY),bm3dx(NNX,NNY),bm3dy(NNX,NNY)
	  double rslowu(RANGE2D)
	  double rslowv(RANGE2D)
	  double umt(RANGE2D)
      double vmt(RANGE2D)
      double umtmean(RANGE2D)
      double vmtmean(RANGE2D)
      double umtold(RANGE2D,2)
      double vmtold(RANGE2D,2)
#endif 

! These are used only for IMPLICITDIFF but we keep them so that computerslow has the same structure
#if defined(IMPLICITDIFF) || defined(TIMESPLIT) 
	  double hhi(RANGE2D)
	  double hhui(RANGE2D)
	  double hhvi(RANGE2D)
	  double taub(RANGE2D)
	  double taubu(RANGE2D)
	  double taubv(RANGE2D)
	  double tautu(RANGE2D)
	  double tautv(RANGE2D)
#endif

! ------------------------------
!FIELD TIME DERIVATIVES
      double ut(RANGE)
      double vt(RANGE)
      double ht(RANGE2D)
      double ct(RANGE,NSCAL)
#ifdef TEMPERATURE    
      double tt(RANGE)
#endif
! ------------------------------
! OLD TIME DERIVATIVES
      double utold(RANGE,2)
      double vtold(RANGE,2)
      double htold(RANGE2D,2)
      double ctold(RANGE,NSCAL,2)
#ifdef TEMPERATURE    
      double ttold(RANGE,2)
#endif

!%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
! Initial ordering of the fields (values 1 or 2 are not important)
      now=1
      old=2  
! Initial ordering of tendencies (time derivative fields) (values 1 or 2 are not important)
      nt1=1
      nt2=2

! ------------------------------
! FIELD TIME DERIVATIVES INITIALIZATION
      utold(:,:,:,:)=0.
      vtold(:,:,:,:)=0.
      htold(:,:,:)=0.
      ctold(:,:,:,:,:)=0.
#ifdef TEMPERATURE    
      ttold(:,:,:,:)=0.
#endif

!%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
!READ PHYSICAL PARAMETERS 
      call readparameters()
	
! ------------------------------
! COMPUTE THE BASIC ELEMENTS FOR THE GRID,  Interpolation and derivatives in Z
      call initparameters()
      call printparameters()
! ------------------------------
!READ ECOLOGICAL PARAMETERS 
      call readparameterseco()
      call printparameterseco()
! ------------------------------
      print *,outfile
      open(unit=fout, file=outfile, form='unformatted' )

!%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
! Reads a bathymetry file and   generates  the mask  for h points
! OR generates  a bathymetry file  and read an external mask for h points
! OR generates  a bathymetry file  and makes simple 'rectangular swimming pool'
! --------------------------------------------------------------

! freadbathy logical true or false in params.dat 
      if(freadbathy) then
!Reads bathymetry from file bathyfile and Computes  mask for  h (or c or w points)  from bathymetry 
        call readbathymetry()
      else
!generates a bathymetry for  h points  and create a file mybathy  
        call initbathymetry()
! freadmask logical true or false in params.dat 
           if(freadmask) then
	   call readllmask()
           else
!Makes simple 'rectangular swimming pool' land-lake         
           call initmask(llm)
           endif
      endif   

      izoff(:,:,:)=0 ! JJNXX	
#ifdef IMMBOUND
      if(freadimmbound) then
!Reads immersed boundary from file immbound 
        call readimmbound()
!      else
!	zimm(:,:)=0.
      endif   
#endif
 
! ----------------------------------
! Compute mask for u and v points (still to be used? XXX)
! ----------------------------------
      call makevecmask(llm)

! ----------------------------------
! Compute 3D mask for u and v points
! ----------------------------------
      call makemask3d(llm,llm3d,izoff,zlc)

! %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
! Initialize fields with initial conditions 
! -------------------------------------------
! Initial time in the day
      ntime=(HOURSTART*3600/ddt)
      print *,'Starting at time ',HOURSTART,' hours = ',ntime,' steps'
      feuler=.true. ! first step should be Euler   
   
! frestart logical true or false in params.dat 

      if(frestart) then
         call restart(FIELDS)
      else
         call initfields(FIELDLIST)
      endif

      mintime=ntime
      maxtime=ntime+numtime

     
#ifdef VARIABLEWIND
      uwind=uwindamp*WINDCYCLE(ntime,ddt)
#endif

!%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
! COMPUTE Horizontal diffusion parameters (uses llm for sponge layer  so it must be here) 
! ----------------------------------
      call inithordiffusion()

!%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
! COMPUTE  VERTICAL diffusion parameters  ONLY ACTIVE when LMIXING IS NOT ACTIVE
!FOR PARABOLIC DIFFUSIOn TSANIS this coeff is K/H where K is the diffusion at point w in the vertical
!FOR CONSTANT  DIFFUSION  this coeff is K where K is the diffusion at point w in the vertical
!For LMIXING IT IS NOT USED
! ----------------------------------
  
#ifndef LMIXING    
      call initvertdiffusion()
#endif

! Fix u,v,h,c on ghost points      
      call boundaries(FIELDS)      
      call boundaries(OLDFLD)

! ----------------------      
      call saveoutput(fout,EXTFIELDS)                                                                                                                                                                                                                                                               
      call checkCFL(EXTFIELDS)

#ifdef TIMESPLIT
! %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
! Initialize time-splitting fields
! %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

	  rslowu(:,:)=0.
	  rslowv(:,:)=0.
	  vmt(:,:)=0.
	  umt(:,:)=0

! Initial ordering of the fields (values 1 or 2 are not important)
      nowfast=1
      oldfast=2  
! Initial ordering of tendencies (time derivative fields) (values 1 or 2 are not important)
      nt1fast=1
      nt2fast=2
! ------------------------------
! FIELD TIME DERIVATIVES INITIALIZATION
	  umtold(:,:,:)=0.
	  vmtold(:,:,:)=0.
	  htold(:,:,:)=0
      umtmean(:,:)=0.d0	  	
c2d      vmtmean(:,:)=0.d0
      feulerfast=.true.  ! we start with an euler step 

! Initialize weights for external mode integration	
	  call makeweightsbaro(NFASTSTEPS,aw,mstar)
	
  	  hm(:,:,1)=h(:,:,1)
	  hm(:,:,2)=h(:,:,2)
		
#ifdef TEMPERATURE
	  call averagevert(u(1,1,1,now),v(1,1,1,now),t(1,1,1,now),um3d,vm3d,bm3du,bm3dv,bm3dx,bm3dy)
#else
	  call averagevert(u(1,1,1,now),v(1,1,1,now),um3d,vm3d)
#endif		
   	
	  hmean(:,:)=0.d0
	  umean(:,:)=0.d0
c2d      vmean(:,:)=0.d0
	  umtmean(:,:)=0.d0	  	
c2d      vmtmean(:,:)=0.d0	  	

	  htold(:,:,:)=0.
	  umtold(:,:,:)=0.
c2d  	  vmtold(:,:,:)=0.

	  feulerfast=.true.
#endif	! TIMESPLIT	

!%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
!  MAIN TIME LOOP
!%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
   
      do while(ntime.le.maxtime) 

! ------------------------------------------------------------
#ifdef VARIABLEWIND
      uwind=uwindamp*WINDCYCLE(ntime,ddt)
      call initvertdiffusion()
#endif
! ------------------------------------------------------------

!utold(nt2) is changed and will contain the derivatives at time now-1=n-1  (To be used in the corrector part)

         ntemp=nt1; nt1=nt2; nt2=ntemp

! ------------------------------------------------------------
! LEAP-FROG   LF-AM timestepping 
! We  get the first guess of the future value u*(t_n+1) (prediction by Leap-Frog)
! using the u(t_n) (now)  and u(t_n-1) (old)
! ------------------------------------------------------------

! Compute the derivatives F(u_n,t_n) and it is put in nt1

!  PREDICTOR
!   predicted  fields u*(n+1) are written in old ; uses tendencies in ut(nt1) which was computed above

#if defined(IMPLICITDIFF) || defined(TIMESPLIT)
	   call computetaub(FIELDSUVH,taub,taubu,taubv,tautu,tautv,hhi,hhui,hhvi)   ! pre-compute using u(now)
#endif

   		call tendencies(
     &       EXTFIELDS,
     &       TENDFIELDS1  )    ! TENDFIELDS1 is F1 (w/o diff)

#ifdef TIMESPLIT

! %%%%%%%%%%%%%%% FAST LOOP %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
! Compute rslow using ut(t) (stored in utold(nt1))
 	    call computerslow(utold(0,NYMIN,0,nt1),vtold(0,NYMIN,0,nt1),umtmean,vmtmean,rslowu,rslowv,taubu,taubv,tautu,tautv)

! ************** BEGIN OF FAST LOOP **********************
! Initial ordering of the fields (values 1 or 2 are not important)
!      nowfast=1
!      oldfast=2  
! Initial ordering of tendencies (time derivative fields) (values 1 or 2 are not important)
!      nt1fast=1
!      nt2fast=2
          
          hmean(:,:)=0.d0
          umean(:,:)=0.d0
c2d      vmean(:,:)=0.d0
		 umtmean(:,:)=0.d0	  	
c2d      vmtmean(:,:)=0.d0	  	

    	  htold(:,:,:)=0.
	      umtold(:,:,:)=0.
c2d	      vmtold(:,:,:)=0.

    	  feulerfast=.true.
	
         do ntimefast=1,mstar
!utold(nt2) is changed and will contain the derivatives at time now-1=n-1  (To be used in the corrector part)
         ntemp=nt1fast; nt1fast=nt2fast; nt2fast=ntemp
! ------------------------------------------------------------
! LEAP-FROG   LF-AM timestepping 
! We  get the first guess of the future value u*(t_n+1) (prediction by Leap-Frog)
! using the u(t_n) (now)  and u(t_n-1) (old)
! ------------------------------------------------------------
! Compute the derivatives F(u_n,t_n) and it is put in nt1
#ifndef TEMPERATURE 
	 call tendenciesfast(
     &       FIELDSFAST,
     &       TENDFIELDS1FAST, rslowu,rslowv )
#else
	 call tendenciesfast(
     &       FIELDSFAST,
     &       TENDFIELDS1FAST, rslowu,rslowv,bm3du,bm3dv,bm3dx,bm3dy )
#endif
    	  umtmean(:,:)=umtmean(:,:)+umtold(:,:,nt1fast)*aw(ntimefast)
c2d	      vmtmean(:,:)=vmtmean(:,:)+vmtold(:,:,nt1fast)*aw(ntimefast)

!  PREDICTOR
!   predicted  fields u*(n+1) are written in old ; uses tendencies in ut(nt1) which was computed above
        call timestepLFAM1fast(
     &    FIELDSFAST,
     &    OLDFLDFAST,
     &    TENDFIELDS1FAST )
 
!   predicted  fields u*(t_n+1) are moved to now
!   fields u(t_n) are moved to old

        ntemp=oldfast; oldfast=nowfast; nowfast=ntemp

! Fix the Boundary conditions for u*(t+1) Not NEEDED!
        call boundariesfast(FIELDSFAST)

! -------------------------------------------------------------
! LEAP-FROG   LF-AM timestepping (Correction by Adams-Moulton)
! WE want to get the value u(t_n+1) 
! -------------------------------------------------------------
!  Compute the derivative at t_n+1 with the predicted  fields u*(t_n+1)  
!   ut,... now contains the derivatives at t_n+1  computed  with the predicted  fields u*(t_n+1)  

#ifndef TEMPERATURE 
	  call tendenciesfast(
     &       FIELDSFAST,
     &       TENDFIELDSFAST, rslowu,rslowv )
#else
	  call tendenciesfast(
     &       FIELDSFAST,
     &       TENDFIELDSFAST, rslowu,rslowv,bm3du,bm3dv,bm3dx,bm3dy )
#endif
	
! CORRECTOR
!   'now' will contain the computed fields u(t+1)  */
!   ut is now+1, ut1 is now and ut2 is now-1
  
        call timestepLFAM2fast(
     &    FIELDSFAST,
     &    OLDFLDFAST,
     &    TENDFIELDSFAST,
     &    TENDFIELDS1FAST,
     &    TENDFIELDS2FAST)

        call boundariesfast(FIELDSFAST)

! weighted averages of the barotropic variables
	    hmean(:,:)=hmean(:,:)+hm(:,:,nowfast)*aw(ntimefast)
	    umean(:,:)=umean(:,:)+um(:,:,nowfast)*aw(ntimefast)
c2d	    vmean(:,:)=vmean(:,:)+vm(:,:,nowfast)*aw(ntimefast)

 	   enddo
		
 ! ************** END OF FAST LOOP **********************

  		hm(:,:,nowfast)=hmean(:,:)
		um(:,:,nowfast)=umean(:,:)
c2d 	  vm(:,:,nowfast)=vmean(:,:)

! h is taken from the barotropic integration 
         h(:,:,old)=hmean(:,:) ! hmmmm Should this not be h(:,:,old) ????? JJN

! %%%%%%%%%%%%%%% END FAST LOOP%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

#endif ! TIMESPLIT

        call timestepLFAM1(
     &    FIELDS,
     &    OLDFLD,
     &    TENDFIELDS1 )     ! In case of IMPLICITDIFFZ this is now u*1=u(t-1)+2*dT*F1(u;t)

!   predicted  fields u*(t_n+1) are moved to now
!   fields u(t_n) are moved to old

        ntemp=old; old=now; now=ntemp

#ifdef IMPLICITDIFF
       call solveimplicit(EXTFIELDS,taub,hhi,hhui,hhvi,2.d0)   ! solves u*-F0(u*)=u*1 for u*  ! Still using old hhi JJN!
#endif

! Fix the Boundary conditions for u*(t+1)
        call boundaries(FIELDS)

#ifdef IMPLICITDIFF 
	   call computetaub(FIELDSUVH,taub,taubu,taubv,tautu,tautv,hhi,hhui,hhvi)   ! pre-compute using u*
#endif

! -------------------------------------------------------------
! LEAP-FROG   LF-AM timestepping (Correction by Adams-Moulton)
! WE want to get the value u(t_n+1) 
! -------------------------------------------------------------

!  Compute the derivative at t_n+1 with the predicted  fields u*(t_n+1)  
!   ut,... now contains the derivatives at t_n+1  computed  with the predicted  fields u*(t_n+1)  

! If implicit timestep, do not include vertical diffusion terms in tendencies

!#ifdef IMPLICITDIFFZ
        call tendencies(EXTFIELDS,      ! We only need F1(u*) 
     &                  TENDFIELDS )

! CORRECTOR
!   'now' will contain the computed fields u(t+1)  */
!   ut is now+1, ut1 is now and ut2 is now-1
! h is not touched if TIMESPLIT is defined
  
        call timestepLFAM2(
     &    FIELDS,
     &    OLDFLD,
     &    TENDFIELDS,
     &    TENDFIELDS1,
     &    TENDFIELDS2)

#ifdef IMPLICITDIFF
		call solveimplicit(EXTFIELDS,taub,hhi,hhui,hhvi,1.d0)   ! taub and hhi have been precomputed from u*
#endif


!ifdef TIMESPLIT
! JJN seems redundant!!!
! h is taken from the barotropic integration 
!          h(:,:,now)=hmean(:,:)
!#endif

! -------------------------------------------------------------
! Convective adjustement
! -------------------------------------------------------------

#ifdef TEMPERATURE
#ifdef CONVADJ
	if(mod(ntime,CONVADJ).eq.0) then
	    call convadj(t(1,1,1,now),c(1,1,1,1,now)) 
	    feuler=.true. ! next step should be Euler
	endif
#endif
#endif

#ifdef TIMESPLIT
! Make sure that the vertical average of the u is the same as the barotropic um

#ifdef TEMPERATURE
	    call averagevert(u(1,1,1,now),v(1,1,1,now),t(1,1,1,now),um3d,vm3d,bm3du,bm3dv,bm3dx,bm3dy)
#else
	    call averagevert(u(1,1,1,now),v(1,1,1,now),um3d,vm3d)
#endif
		 do k=2,NNZ
		            u(:,:,k,now)=u(:,:,k,now)-um3d(:,:)+umean(:,:)
c2d            v(:,:,k,now)=v(:,:,k,now)-vm3d(:,:)+vmean(:,:)
		 enddo		
#endif


        call boundaries(FIELDS)

! ---------------------------------
        ntime=ntime+1
    
    	if(mod(ntime,ncfl).eq.0) then 
           call checkCFL(EXTFIELDS)
        endif

        call saveoutput(fout,EXTFIELDS)
! ----------------------------------

!	 write(66,*) ntime, h(2,1,now)

 	enddo


 ! ************** END OF MAIN LOOP **********************
 
      close(fout)
      end
!%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%% 
!  ------------------------ END  of Main ------------------ 
!%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%% 


#ifdef TIMESPLIT

      subroutine makeweightsbaro(m,a,mstar)
!     Computes weights for external mode integration à la ROMS 
!     The weights are such that sum(a)=sum(a*t)=sum(a*t^2)=1
!     with t=0...taustar
!     taustar=mstar/m  is longer than 1 (about 1.31)....

	  implicit none
	  integer n,it,p,q,mstar,m
	  double tau,r,f,dt,a(2*m),t,ap

	  r=0.284578345839306
	  f=0.948986344071042
	  p=2
	  q=4
	  tau=(p+2.d0)*(p+q+2.d0)/(p+1.d0)/(p+q+1.d0)*f;
	  dt=1.d0/m
	  ap=0.
	  do it=1,2*m
	    t=it*dt
	    a(it)=(t/tau)**p*(1-(t/tau)**q)-r*t/tau;
	    if(ap*a(it).lt.0) mstar=it-1
	    ap=a(it)
	  enddo
	  a(:)=a(:)/sum(a(1:mstar))
	  end
! %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%% 
! -----------------------------------------

#ifndef TEMPERATURE
      subroutine tendenciesfast(um,vm,h,umt,vmt,ht,rslowu,rslowv)
#else
	  subroutine tendenciesfast(um,vm,h,umt,vmt,ht,rslowu,rslowv,bm3du,bm3dv,bm3dx,bm3dy )
#endif
! -----------------------------------------
! FAST timesteps for the shallow water eqs.
! Computes  time derivatives  of um,vm,h
! -----------------------------------------
    
      implicit none
#include "commonlac.h"
     
! ---------------------------------
       double h(RANGE2D)
       double um(RANGE2D)
       double vm(RANGE2D)
       double rslowu(RANGE2D)
       double rslowv(RANGE2D)
	   double bm3du(RANGE2D)  ! int(b s ds) interpolated to u 
	   double bm3dv(RANGE2D)  ! int(b s ds) interpolated to u
	   double bm3dx(RANGE2D) ! DX(bm3d)
	   double bm3dy(RANGE2D) ! DY(bm3d)
       double umt(RANGE2D)
       double vmt(RANGE2D)
       double ht(RANGE2D)

       double uhi(RANGE2D)
       double vhi(RANGE2D)      
 
! ---------------------------------    
       double hh(RANGE2D)		! H in h
       double hhi(RANGE2D)               ! 1/H in h
       double hhu(RANGE2D)		! H in a u point 
       double hhui(RANGE2D)		! 1/H in u  
c2d      double hhv(RANGE2D)            ! H in a  v point
c2d      double hhvi(RANGE2D)           ! 1/H in  v
       double hhoz(RANGE2D)		! H in oz    
! ----------------------------------------------------------------------------
! TBC

         double worku(RANGE2D)  
         double workv(RANGE2D)
         double workuhp(RANGE2D) 
         double workvhp(RANGE2D) 
         double workuhi(RANGE2D) 
         double workvozi(RANGE2D) 
         double workuoz(RANGE2D)
 
         double workuh(RANGE2D) 
         double workvh(RANGE2D) 


         double workvoz(RANGE2D) 

         double workuozi(RANGE2D)
         double workvhi(RANGE2D) 

! ----------------------------------------------------------------------------
      integer i,j,k,l

! ----------------------------------------------------------------------------

      j=1  ! Necessary for the 2D CASE 

! ***************************************************************
!  A) Compute some help fields
! ***************************************************************

! ----------------------------------------------
!   1) Compute H=h+z in h,  u and v grid points
! ----------------------------------------------

! Computes total thickness  at h points
         hh(:,:)=h(:,:)+z(:,:)

! Fix land u points to zmin;  Needed to allow computation of oz points 
c2d      hh(:,:)=hh(:,:)*llm(:,:,HGRID)+(1-llm(:,:,HGRID))*zmin

! Computes total thickness hh  at u points : hhu 
! At the boundaries  u ,  hh is multiplied by u =0 in the equation so no exact values are needed for BC

         BLOCKIX2(hh,hhu,0) 

! Computes total thickness hh  at v points : hhv
! At the boundaries v, hh is mulitplied by v =0 in the equation so no BC needed
        
c2d      BLOCKIY2(hh,hhv,0)

! Computes total thickness hh  at oz points : hhoz  
           
c2d      BLOCKIY2(hhu,hhoz,0) ! 

!  From now on we will also need 1/H (divisions are expensive to compute!)

         hhui(0:NX,:)=1.d0/hhu(0:NX,:)
c2d      hhvi(:,0:NY)=1.d0/hhv(:,0:NY)
     	 hhi(:,:)=1.d0/hh(:,:) 

!  ********************************************************************
!  B) Compute dh/dt and w from the depth-integrated continuity equation 
!  ********************************************************************

! -----------------------------------------------------------------
!   1) Compute uH, vH and other help fields
! -----------------------------------------------------------------

	   worku(:,:)=um(:,:)*hhu(:,:)
c2d	      workv(:,:)=vm(:,:)*hhv(:,:)
   
! -----------------------	
  
	 BLOCKDX2(worku,workuhp,-1)      ! D(uh)/DX in c  (no BC needed, as it does not live on boundary) 
c2d	        BLOCKDY2(workv,workvhp,-1)       ! D(vh)/DY in c       
c2d         workuhp(:,:)=workuhp(:,:)+workvhp(:,:)       ! workuhp is now DX(uh)+DY(vh) in c

! -----------------------	
#ifdef IMMBOUND
    	 FIXBC2D(workuhp,HGRID)
#endif

! ------------------------------------------------------------------
!    3) Compute dh/dt: D_t(h)=-D_x(Int(uH;dz))-D_y(Int(vH;dz))
! ------------------------------------------------------------------
        ht(:,:)=-workuhp(:,:)

!  ****************************************************************
!  E)  Momentum equation for U
!  ****************************************************************

#ifdef FIXED_VEL_PROF
       ut=0
c2d       vt=0
	return  
#endif

!  ------------------------------------------------------------------ 
!   1) Advective part:  ADV = - u D_x(u) - v D_y(u)  -  at point u 
!  -------------------------------------------------------------- 

! Velocities interpolated to help points     
         BLOCKIX2(um,workuhi,-1)        ! u velocity interpolated to point h, no BC for u*c_x since Boundary not on h 
c2d         BLOCKIX2(vm,workvozi,0)  ! IX(v) in oz    

! -----------------------------------------------

! Derivatives of u
         BLOCKDX2(um,workuh,-1)           ! DX(u) in h    no BC
c2d      BLOCKDY2(um,workuoz,0) 		    ! DY(u) in o  
c2d      FIXBC2D(workuoz,OGRID)         ! Dy(u)=0 on boundaries - free slip


! -----------------------------------------------
! Now multiply velocities and derivatives
         workuhi(:,:) =workuh(:,:) *workuhi(:,:)   ! u*Dx(u) in h
c2d      workvozi(:,:)=workuoz(:,:)*workvozi(:,:)  ! v*Dy(u) in oz
! -----------------------------------------------

		 BLOCKIX2(workuhi,worku,0)
c2d      BLOCKIY2(workvozi,workv,-1)
		 umt(:,:) = -worku(:,:) + rslowu(:,:)
c2d		 umt(:,:)=umt(:,:)-workv(:,:)
! -----------------------------------------------

! -----------------------------------------------------------------------------------
!   2) Add pressure balance + horizontal and vertical diffusion to u tendency
!  D_t(u) = ADV - g D_x(h) + nu 1/H D_x ( H D_x (u) ) + nu 1/H D_y (H D_y(u) ) )
! -----------------------------------------------------------------------------------

! workuh -->  H du/dx
! workuoz --> H du/dy

         workuh(:,:)=workuh(:,:)*hh(:,:)*znuh(:,:)
c2d      workuoz(:,:)=workuoz(:,:)*hhoz(:,:)*znuoz(:,:)

         BLOCKDX2(h,vhi,0) ! Dx(h)

c		print *,'UMT0:',umt(:,1)
c2d        do j=1,NY
           do i=1,NX
              umt(i,j)= umt(i,j)
     &              -GG*vhi(i,j)                            ! h_x is computed here, but used for u, which is 0 on the boundary anyway
     &              +hhui(i,j)*( DX2(workuh,i,j)
c2d     &           +            DY2(workuoz,i,j-1)
     &    )
#ifdef TEMPERATURE
     & 		    + 2*vhi(i,j)*bm3du(i,j)+ hhu(i,j)*bm3dx(i,j)	  
#endif

           enddo
c2d      enddo
	
#ifdef c2d

!  ****************************************************************
!  F)  Momentum equation for V
!  ****************************************************************

!  ------------------------------------------------------------------ 
!   1) Advective part:  ADV = - u D_x(v) - v D_y(v)  - w/H D_w(v) 
!  ------------------------------------------------------------------ 

! Velocities interpolated to help points
 
         BLOCKIY2(um,workuozi,0)             ! IY(u) in oz    
c2d       BLOCKIY2(vm,workvhi,-1) ! Velocity interpolated to  point h,   no BC for v*c_y since Boundary not on h 

!  ------------------------------------------------------------------ 
! Derivatives of v
       BLOCKDY2(vm,workvh,-1)              ! DY(v) in h, no BC
       BLOCKDX2(vm,workvoz,0)              ! DX(v) in o 
       FIXBC2D(workvoz,OGRID)            ! impose free slip BC    

! -----------------------------------------------
! Now multiply velocities and derivatives 
      workuozi(:,:)=workvoz(:,:)*workuozi(:,:)   ! u*Dx(v) in oz
      workvhi(:,:) =workvh(:,:) *workvhi(:,:)    ! v*Dy(v) in h
! -----------------------------------------------
		
		 BLOCKIX2(workuozi,worku,-1)
c2d       BLOCKIY2(workvhi,workv,0)
		 vmt(:,:)=-worku(:,:)-workv(:,:)+rslowv(:,:)
					     
! -----------------------------------------------------------------------------------
!   2) Add pressure balance + horizontal and vertical diffusion to v tendency
!  D_t(v) = ADV - g D_y(h) + nu 1/H D_x H D_x (v) + nu 1/H D_y H D_y (v)  + 1/H D_w( Kv/H D_w(v) ) 
! -----------------------------------------------------------------------------------


! workvh  --> H dv/dy
! workvoz --> H dv/dx
         workvh(:,:)=workvh(:,:)*hh(:,:)*znuh(:,:)
     	 workvoz(:,:)=workvoz(:,:)*hhoz(:,:)*znuoz(:,:)

       BLOCKDY2(h,vhi,0) ! Dx(h)

        do j=1,NY
          do i=1,NX
            vmt(i,j)=vmt(i,j)
     &           -GG*vhi(i,j) !     h_y is computed but used for v, so OK
     &           +hhvi(i,j)*( DX2(workvoz,i-1,j)
     &           +            DY2(workvh,i,j) )
#ifdef TEMPERATURE
     & 		    + 2*vhi(i,j)*bm3dv(i,j)+ hhv(i,j)*bm3dy(i,j)	  
#endif
          enddo
        enddo

! -----------------------------------------------
#endif

      end


#endif ! TIMESPLIT


! %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%% 
! -----------------------------------------
#ifdef TEMPERATURE
      subroutine tendencies(u,v,w,h,c,t,ut,vt,ht,ct,tt)
#else
      subroutine tendencies(u,v,w,h,c,ut,vt,ht,ct)
#endif
! -----------------------------------------

! -----------------------------------------
! Computes  time derivatives  of u,v,h, c and t
! and derives the current w field
! -----------------------------------------
    
!ON THE  SIDES walls: Free slip  for velocity
!ON THE  SIDES walls: Tracers Zero flux OK
!On TOP surface :  SHEAR AT THE SURFACE for velocity: use  BCs   for du/dz,dv/dz on the points u and v respectively (shear)
!On TOP surface :  Tracers Zero flux

!     Staggered - Arakawa C Grid
!     
!     +           |                     |
!   u(0,2) --- h(1,2) --- u(1,2) --- h(2,2) --
!     +           |                     |
!     +           |                     |
!     +        v(1,1)     oz(1,1)    v(2,1) 
!     +           |                     |
!     +           |                     |
!   u(0,1) --- h(1,1) --- u(1,1) --- h(2,1) --
!     +           |                     |
!     +           |                     |
!     ++++++++ v(1,0) ++++++++++++++ v(2,0) ++
!     (+++ = domain boundary)
!
!   All the coordinates are expressed in the h grid reference frame 
!   Using the workxx arrays we are sacrificing memory for clarity 
!
!  In the vertical, the w fields live 1/2 step above the h fields
!
!                 |
!              w(*,*,1)
!                 | 
!              c(*,*,1)
!                 |
!  ++++++++++  w(*,*,0) ++++++ (boundary) ++++++++++
!                 |               
!              c(*,*,0)

!
!  Explanation of calls to BIX and BDX
!  call BIX(u,workuh,-1)
!  computes u interpolated to the h positions   i.e. u(x-dx/2)     
!  interpolation uses    u(x-dx)  and  u(x) 

!  call BIY(u,workuoz,0)
!  computes u interpolated to the cell center 'oz' positions      
!  i.e. u(x,y,z) interpolated with u(x,y+dy,z)

!  the parameter  -1 or 0 is an offset, 
! ------------------------------------------------------------------

      implicit none
#include "commonlac.h"
#include "commoneco.h"
     
! ---------------------------------
      double h(RANGE2D)
      double u(RANGE)
      double v(RANGE)
      double w(RANGE)
      double c(RANGE,NSCAL)

#ifdef TIMESPLIT
      double htemp(RANGE2D) !!JFAST
#endif

      double ut(RANGE)
      double vt(RANGE)
      double ht(RANGE2D)
      double ct(RANGE,NSCAL)
! ---------------------------------
#ifdef TEMPERATURE
      double t(RANGE)
      double tt(RANGE)
#endif
! ---------------------------------

      double phisurf(RANGE2D),phibulk

      double uhi(RANGE2D)
      double vhi(RANGE2D)      
 
! ---------------------------------    
      double hh(RANGE2D)		! H in h
      double hhi(RANGE2D)               ! 1/H in h
      double hhu(RANGE2D)		! H in a u point 
      double hhui(RANGE2D)		! 1/H in u  
c2d      double hhv(RANGE2D)            ! H in a  v point
c2d      double hhvi(RANGE2D)           ! 1/H in  v
      double hhoz(RANGE2D)		! H in oz
! -------------------------------------

#ifndef IMPLICITDIFF
      double dutop(RANGE2D)
c2d      double dvtop(RANGE2D) 
#endif

      double duhi(RANGE2D)
c2d      double dvhi(RANGE2D) 
    
! ----------------------------------------------------------------------------
c	There are only 6 large fields, all others are equivalences

         double worku(RANGE)  
         double workv(RANGE)  
         double workuh(RANGE) 
         double workvh(RANGE) 
         double workuhp(RANGE) 
         double workvhp(RANGE) 

         double workuoz(RANGE)
         double workuoy(RANGE)

         double workwoyi(RANGE)
         double workwoxi(RANGE)  
         double workuozi(RANGE)
         double workvozi(RANGE)

         double workvoz(RANGE)
         double workvox(RANGE)

         double workuhi(RANGE) 
         double workvhi(RANGE) 

         double workcui(RANGE)  
         double workcvi(RANGE)           
         double workcu(RANGE)  
         double workcv(RANGE)  
         double workcw(RANGE)  
         double workcc(RANGE)  
	     double fctgamma(RANGE)  

         double worktuu(RANGE)  
         double worktvv(RANGE)
         double worktww(RANGE)           
         double worktu(RANGE)  
         double worktv(RANGE)  
         double worktw(RANGE)  
        
	 double workbuoyi(RANGE)

! ----------------------------------------------------------------------------
#ifdef TEMPERATURE
	 double workbuoy(RANGE)    ! Not equivalenced, needed everywhere
#endif
! ----------------------------------------------------------------------------
#if defined(RADIATION) || defined(ECOSYSTEM)
	 double light(RANGE),lightz(RANGE)
#endif
! ----------------------------------------------------------------------------
 
! ----------------------------------------------------------------------------
! We use equivalences to use less memory but maintain clarity
! ----------------------------------------------------------------------------
! These are a bit delicate. 
! Here is a scheme of what is used in which section and of equivalences
! EQUIVALENCE ALONG LINES  ; "-" Means still in use
! +----------------------------------------+
! | A/B |  C  |  D  |     E    |  F        |    
! +----------------------------------------+
! | u   | cu  | tu  | woyi     | uozi      |    
! | v   | cv  | tv  |vozi/buoyi| woxi/buoyi|    
! |     | cc  | tuu | uhi      | vh        |    * also fctgamma is aliased with these (used at the end of C after workcc) 
! |     | cw  | tw  | uh       | vhi       |    
! | vhp | cvi | tvv | uoy      | vox       |    
! | uhp | cui | tww | uoz      | voz       |   
! +----------------------------------------+

	 equivalence(worku,workcu,workwoyi,workuozi,worktu)
	 equivalence(workv,workcv,workwoxi,workvozi,worktv,workbuoyi)
	 equivalence(workvh,workuhi,workcc,fctgamma,worktuu)
	 equivalence(workcw,workuh,workvhi,worktw)
	 equivalence(workuoy,workvox,workcvi,workvhp,worktvv)
	 equivalence(workuoz,workvoz,workcui,workuhp,worktww)

! ----------------------------------------------------------------------------
      integer i,j,k,l,fstop
      double cff,sum,cwind2,swind2,ev
      double am,ap,a0,dp,dm,sp,sm
      
! ----------------------------------------------------------------------------
      
      j=1  ! Necessary for the 2D CASE 

! ***************************************************************
!  A) Compute some help fields
! ***************************************************************

! ----------------------------------------------
!   1) Compute H=h+z in h,  u and v grid points
! ----------------------------------------------

! Computes total thickness  at h points
         hh(:,:)=h(:,:)+z(:,:)

! Fix land u points to zmin;  Needed to allow computation of oz points 
c2d      hh(:,:)=hh(:,:)*llm(:,:,HGRID)+(1-llm(:,:,HGRID))*zmin

! Computes total thickness hh  at u points : hhu 
! At the boundaries  u ,  hh is multiplied by u =0 in the equation so no exact values are needed for BC

         BLOCKIX2(hh,hhu,0) 

! Computes total thickness hh  at v points : hhv
! At the boundaries v, hh is mulitplied by v =0 in the equation so no BC needed
        
c2d      BLOCKIY2(hh,hhv,0)

! Computes total thickness hh  at oz points : hhoz  
           
c2d      BLOCKIY2(hhu,hhoz,0) ! 

!  From now on we will also need 1/H (divisions are expensive to compute!)

         hhui(0:NX,:)=1.d0/hhu(0:NX,:)
c2d      hhvi(:,0:NY)=1.d0/hhv(:,0:NY)
	 hhi(:,:)=1.d0/hh(:,:) 

! ----------------------------------------------------------------------
!  2) Compute top BCs for du/dz,dv/dz on the points u and v respectively 
! ----------------------------------------------------------------------

#ifndef IMPLICITDIFF
      cwind2=cos(WINDRAD) ! 
      swind2=sin(WINDRAD) !  

! Compute   du/dsigma on the points u at TOP SURFACE  sigma =1 
#ifdef LMIXING
         dutop(:,:)=(TAUS/RHOWAT)/kvmin*hhu(:,:)*cwind2	   
#else
#ifdef KVPAR 
         dutop(:,:)=(TAUS/RHOWAT)/Kv(NZ)*cwind2	    ! in this case Kv is already Kv'=Kv/H
#else
	     dutop(:,:)=(TAUS/RHOWAT)/KVERT*hhu(:,:)*cwind2
#endif
#endif

! Compute top BCs for dv/dsigma on the points v at TOP SURFACE   sigma =1  
#ifdef c2d
#ifdef LMIXING
         dvtop(:,:)=(TAUS/RHOWAT)/kvmin*hhv(:,:)*swind2      
#else
#ifdef KVPAR 
         dvtop(:,:)=(TAUS/RHOWAT)/Kv(NZ)*swind2    ! Kv is already Kv'=Kv/H
#else
	 dvtop(:,:)=(TAUS/RHOWAT)/KVERT*hhv(:,:)*swind2
#endif
#endif
#endif	

#endif ! IMPLICITDIFF

! -----------------------------------
! Compute buoyancy from temperature
! -----------------------------------
#ifdef TEMPERATURE
      workbuoy(:,:,:)=GG*(1.-RHOEOS(t(:,:,:))/RHOWAT)  ! buoyancy
#endif

! -------------------------------------------------------------------------------------------
! Compute vertical diffusivity using mixing length theory and Pacanowski and Philander 1981
! In the other cases (Tsanis or Constant) already done in    initvertdiffusion()
!   kvhmix stands for Kv/H 
! -------------------------------------------------------------------------------------------
 
#ifdef LMIXING
         BLOCKIX3(u,workuhi,-1)    ! u velocity interpolated to point h, 
c2d      BLOCKIY3(v,workvhi,-1)  ! v velocity interpolated to point h, 
         BLOCKDZC3(workuhi,worktu,0) ! Du/Dsigma   in w
c2d      BLOCKDZC3(workvhi,worktv,0) ! Dv/Dsigma   in w

! ------------------
! Part with the shear
         worktu(:,:,:)=abs(worktu(:,:,:))
c2d      worktu(:,:,:)=sqrt(worktu(:,:,:)**2+worktv(:,:,:)**2)

#ifndef IMMBOUND
       do k=0,NZ-1
!               kvhmix(:,:,k)=(LAMBDA*(zlw(k)+ZBH)*(1-zlw(k)+ZSH))**2*worktu(:,:,k) 
                kvhmix(:,:,k)=(KARMAN*zlw(k)*(1-zlw(k)))**2*worktu(:,:,k) ! Kv/H  JJIMM
       enddo  
#else
#define SIGMABOT zlw(izoff(:,:,HGRID))
       do k=0,NZ-1
            kvhmix(:,:,k)=(KARMAN*(zlw(k)-SIGMABOT)/(1-SIGMABOT)
    &         *(1-zlw(k)))**2*worktu(:,:,k) ! Kv/H  
       enddo  
#endif

! ------------------
#if defined(TEMPERATURE) && defined(PACANOWSKI)
! Compute N^2=d( buoyancy )/dz  = 1/H d( buoyancy )/dsigma

         BLOCKDZC3(workbuoy,worktv,0) ! H N^2 in w 
         worktv=max(0.d0,worktv)
	 where(worktu(:,:,0:(NZ-1)).gt.0.d0)  worktv(:,:,0:(NZ-1))=worktv(:,:,0:(NZ-1))/(worktu(:,:,0:(NZ-1))**2) ! local gradient Richardson number
         do k=0,NZ-1  ! all derivatives are by sigma, Ri needs still to be multiblied by H
		worktv(:,:,k)= worktv(:,:,k)*hh(:,:)
	 enddo
         kvhmix(:,:,0:(NZ-1))=kvhmix(:,:,0:(NZ-1))*(1+ALPHAPAC*worktv(:,:,0:(NZ-1)))**(NPAC)	! This is Kv/H
#endif

         do k=0,NZ-1
               kvhmix(:,:,k)=kvhmix(:,:,k)+kvmin*hhi(:,:)
         enddo  

#ifdef IMPLICITDIFF
		kvhmix(:,:,NZ)=kvmin*hhi(:,:)
! These will be used later in the solveimplicit part
	 	BLOCKIX3(kvhmix,kvhmixu,0)
c2d				 BLOCKIY3(kvhmix,kvhmixv,0)
#endif

#endif

! -------------------------------------------------------------------------------------------
! Compute light and the  derivative

#if defined(RADIATION) || defined(ECOSYSTEM)
      call complight(c,light,lightz,hhi,hh)
#endif

#ifndef IMPLICITDIFF
#ifdef TEMPERATURE
      phisurf(:,:)=hsurf*(RHOAIR*CAIR)*(Tair-t(:,:,NZ))  ! Top boundary condition [ Flux=-h(Tair-T) ] 
                                                         ! Using top level (error is less than 2%)
#endif
#endif

!  ********************************************************************
!  B) Compute dh/dt and w from the depth-integrated continuity equation 
!  ********************************************************************

! -----------------------------------------------------------------
!   1) Compute uH, vH and other help fields
! -----------------------------------------------------------------

         do k=0,NZ
	   worku(:,:,k)=u(:,:,k)*hhu(:,:)
c2d	   workv(:,:,k)=v(:,:,k)*hhv(:,:)
         enddo	

! -----------------------	
  
	 BLOCKDX3(worku,workuhp,-1)      ! D(uh)/DX in c  (no BC needed, as it does not live on boundary) 
c2d	 BLOCKDY3(workv,workvhp,-1)       ! D(vh)/DY in c
         
c2d      workuhp(:,:,:)=workuhp(:,:,:)+workvhp(:,:,:)       ! workuhp is now DX(uh)+DY(vh) in c

! -----------------------	
#ifdef IMMBOUND
	 FIXBC3D(workuhp,HGRID) 
#endif

! -----------------------------------------------------------------
!   2) Compute vertically integrated Dx(uH)+Dy(vH)
! -----------------------------------------------------------------

	w(:,:,0)=0.d0   ! Integral is 0 at bottom


! --------------------------------

! Regular, square integration
        do k=1,NZ
            w(:,:,k)=w(:,:,k-1)+ workuhp(:,:,k)*zwdif(k) 
        enddo

! ---------------------------------

! ------------------------------------------------------------------
!    3) Compute dh/dt: D_t(h)=-D_x(Int(uH;dz))-D_y(Int(vH;dz))
! ------------------------------------------------------------------

#ifdef TIMESPLIT
        htemp(:,:)=-w(:,:,NZ)
#else
		ht(:,:)=-w(:,:,NZ)
#endif

! ------------------------------------------------------------------
!  4) Obtain "vertical velocity" w=(dsigma/dt*H) = -dh/dt*sigma-Int( Dx(uH)+Dy(vH) dsigma,0,sigma) 
! ------------------------------------------------------------------
	  do k=1,NZ
#ifdef TIMESPLIT
	   w(:,:,k)=-htemp(:,:)*zlw(k)  -w(:,:,k)
#else
	   w(:,:,k)=-ht(:,:)*zlw(k)  -w(:,:,k)
#endif
	  enddo

!  **************************************************************************************
!  C) Scalar advection and diffusion evaluated at point c or h 
!  **************************************************************************************
!   Please remember that c(:,:,:,l) is actually c*H !
!   Advection+Diff.:   D_t(cH) = - D_x(Jx*) - D_y(Jy*)  - D_sigma(Js*)   at point h or c 
!   Fluxes:    Jx*=H(-nu*D_x(c))+uHc;  Jy*=H(-nu*D_y(c))+vHc; Js*=(w+ws)*c-Kc/HDz(c)
! --------------------------------------------------------------------

!%%%%%%%
      do l=1,NSCAL        !  Loop over all scalar species
! ------------------ 
!	 True c (concentration) field
	 do k=0,NZ
	    workcc(:,:,k)=c(:,:,k,l)*hhi(:,:)
	 enddo
! ------------------

! Derivatives
         BLOCKDX3(workcc,workcu,0)  ! Derivative_x of c, DX(c) living in u 
         FIXBC3D(workcu,UGRID)      ! Fix the Boundary condition of  no flux at the side boundary  DX(c) =0    (J_x(c)=0) 

c2d      BLOCKDY3(workcc,workcv,0)  ! derivative Velocity  DY(c) which lives at  point v      
c2d      FIXBC3D(workcv,VGRID)   ! Fix the Boundary condition of  no flux at the side boundary  DY(c) =0    (J_y(c)=0)                  

! ------------------

#ifndef IMPLICITDIFF
        BLOCKDZC3(workcc,workcw,0) !     D(c)/Dsigma which lives at point w 

! ------------------
!  Define     workcw  =  K/H  D(c)/Dsigma

        do k=0,NZ-1
#ifdef LMIXING
             workcw(:,:,k)=workcw(:,:,k)*(kvhmix(:,:,k)+sdiff(l)*hhi(:,:)) !  kvhmix stands for Kv/H **  kvmin  is included AGAIN
#else
#ifdef KVPAR
 	     workcw(:,:,k)=workcw(:,:,k)*(Kv(k)+sdiff(l)*hhi(:,:))  ! In the parabolic case Kv is already Kv/H 
#else
 	     workcw(:,:,k)=workcw(:,:,k)*(KVERT+sdiff(l))*hhi(:,:)
#endif
#endif
        enddo
! -----------------
#endif


! ----------------------------------------------------------------------------
! Add settling velocity for sinking species and the convective flux wc
! not at the top, we have a BC there (already imposed above, workcw(:,:,NZ)=0)
! ----------------------------------------------------------------------------

! workcw -> J_sigma*=   +(w+ws)*c
! - Kc/H dc/d sigma  is now done in the implcit part

      do k=0,NZ-1   
#ifdef	IMPLICITDIFF
	     workcw(:,:,k)=IZC3(workcc,:,:,k)*(wsett(l)+w(:,:,k))
#else
   		 workcw(:,:,k)=-workcw(:,:,k)+IZC3(workcc,:,:,k)*(wsett(l)+w(:,:,k))
#endif
      enddo

#ifndef	IMPLICITDIFF
! -----------------
!  Top BC  J_sigma=0 at the surface 
	  if(l.eq.SCAL_NUTRIENT) then
#ifdef NUTCYCLE
            workcw(:,:,NZ)=ctopnut*NUTCYCLE(ntime,ddt)
#else
            workcw(:,:,NZ)=ctopnut
#endif
	  else
  	    workcw(:,:,NZ)=0.    
	  endif
#else
  		workcw(:,:,NZ)=0.   ! JJN
#endif

! -----------------
! workcu --> Jx*= u H c - nu H dc/dx    at point u
! workcv --> Jy*= v H c - nu H dc/dy    at point v    
!  cH at u and v positions
         BLOCKIX4(c,workcui,0,l)
c2d      BLOCKIY4(c,workcvi,0,l)
! -----------------
      do k=0,NZ
         workcu(:,:,k)=workcui(:,:,k)*u(:,:,k)-workcu(:,:,k)*hhu(:,:)*znu
c2d	 workcv(:,:,k)=workcvi(:,:,k)*v(:,:,k)-workcv(:,:,k)*hhv(:,:)*znu
      enddo


! --------------------------------------------------------------
! Flux corrected transport scheme for scalars
! Guarantees positive definite fields by limiting outgoing flux
! --------------------------------------------------------------
! -----------------
#ifdef FCT_SCAL

      ev=1e-30

! -------------------------------------------------------------
      do k=1,NZ

! -----------------
#ifdef TANPAR
	 cff=1.d0/zwdif(k)
#else
	 cff=cffsigma
#endif
! -----------------

c2d      do j=1,NY
           do i=1,NX

! Compute contribution to tendencies of outgoing flux
  	    fctgamma(i,j,k)=
     &      + (max(workcu(i,j,k),0.d0)-min(workcu(i-1,j,k),0.d0))*cffx
c2d     &   + (max(workcv(i,j,k),0.d0)-min(workcv(i,j-1,k),0.d0))*cffy
#ifndef IMPLICITDIFF
     &      + (max(workcw(i,j,k),0.d0)-min(workcw(i,j,k-1),0.d0))*cff
#endif

! Compute FCT factor = scalar / what should be removed by outgoing flux
            fctgamma(i,j,k)= min(1.d0,c(i,j,k,l)/(ddt*fctgamma(i,j,k)+ev))

           enddo
c2d      enddo
      enddo
! -------------------------------------------------------------     	
! Now correct the outgoing fluxes (external domain boundaries excluded

      do i=1,NX-1
  	    workcu(i,:,:)= fctgamma(i,:,:)*  max(workcu(i,:,:),0.d0) 
     &                    +fctgamma(i+1,:,:)*min(workcu(i,:,:),0.d0)
      enddo
! ---------      
c2d      do j=1,NY-1
c2d  	    workcv(:,j,:)= fctgamma(:,j,:)*  max(workcv(:,j,:),0.d0) 
c2d     &                 +fctgamma(:,j+1,:)*min(workcv(:,j,:),0.d0)
c2d      enddo
! ---------
#ifndef IMPLICITDIFF
      do k=1,NZ-1
  	    workcw(:,:,k)= fctgamma(:,:,k)*  max(workcw(:,:,k),0.d0) 
     &                    +fctgamma(:,:,k+1)*min(workcw(:,:,k),0.d0)
      enddo
#endif

#endif
 
! ---------  
        do k=1,NZ
c2d      do j=1,NY         
           do i=1,NX
              ct(i,j,k,l)=   
     &	        -           DX3(workcu,i-1,j,k)
c2d     &       -           DY3(workcv,i,j-1,k)
#ifndef IMPLICITDIFF
     &	        -           DZW3(workcw,i,j,k-1)
#endif
           enddo
c2d      enddo
        enddo
! -------------------------------------------------------------     	
!%%%%%%%%%%
      enddo                !  End scalar species loop
!%%%%%%%%%%

! --------------------------------------------
! Add ecological reaction terms to tendencies  ***PROBLEM * MUST BE  H*F  ?
! We provide to it five overwritable work arrays 
! --------------------------------------------

#ifdef ECOSYSTEM
#ifdef IMPLICITDIFF
! The rest is defined in ecotendpatlhs (Patankard method)
       call ecotendpatrhs(c,ct,hhi)
#else
  	   call ecotend(c,ct,light,workcc,workcw,workcu,workcv,workcui,hhi)
#endif
#endif
 
!  **************************************************************************************
!  D) Temperature advection and diffusion evaluated at point h 
!  **************************************************************************************
!------------------------------------------------------------------------
#ifdef TEMPERATURE
 
! -------------------------------------------------------------------------------
!  1) Advection:   ADV = - u D_x(T) - v D_y(T)  - w/H D_sigma(T)   at point h 
! -------------------------------------------------------------------------------
! Derivatives
         BLOCKDX3(t,worktu,0)      ! Derivative  DT/DX  living in u
         FIXBC3D(worktu,UGRID)   ! Fix the Boundary condition :  no flux at the side boundary   DT/DX =0    

c2d      BLOCKDY3(t,worktv,0)  ! derivative Velocity  DT/DY which lives at  point v      
c2d      FIXBC3D(worktv,VGRID)   ! Fix the Boundary condition :  no flux at the side boundary  DT/DY =0                    

         BLOCKDZC3(t,worktw,0) !  derivative   DT/Dsigma which lives at  point w 
                
! for now the value at boundaries is not important since (w is zero at surface and bottom) but  worktw  will be the -J_sigma in 2) 
! and  bottom  BC  J_sigma=0  

	 worktw(:,:,NZ)=0.    
	 worktw(:,:,0)=0.    

! Multiply velocities by derivatives 
         worktww(:,:,:)=worktw(:,:,:)*w(:,:,:)      ! w*Dz(t) in w; the BC for dT/dz not needed since w =0 on the top
         worktuu(:,:,:)=worktu(:,:,:)*u(:,:,:)        ! u*Dx(t) 
c2d      worktvv(:,:,:)=worktv(:,:,:)*v(:,:,:)      ! v*Dy(t)

! ------------------- 
        do k=1,NZ
c2d      do j=1,NY
           do i=1,NX
            tt(i,j,k)= 
     &        - IX3(worktuu,i-1,j,k) 
     &        - hhi(i,j)*IZW3(worktww,i,j,k-1)
c2d     &     - IY3(worktvv,i,j-1,k)
           enddo
c2d      enddo
        enddo
! ------------------- 

! ------------------------------------------------------------------------------------
!  2) adding Diffusion: !  D_t(T) = ADV + nu 1/H D_x (H D_x(T)) + nu 1/H D_y (H D_y(T))  + 1/H D_w( Kv/H D_sigma(T) ) at point h or c 
! ------------------------------------------------------------------------------------

#ifndef IMPLICITDIFF
!  Define   worktw = Kv/H D_sigma(T)    (the negative vertical flux of T  )

! ------Condition at the surface on Kv/H D_sigma(T)  ------------- 
      worktw(:,:,NZ)=phisurf(:,:)*(1./(RHOWAT*CWAT))  

! ---------Other points---------- 
        do k=0,NZ-1
#ifdef LMIXING
             worktw(:,:,k)=worktw(:,:,k)*kvhmix(:,:,k)     ! kvhmix is Kv/H 
#else

#ifdef KVPAR
 	     worktw(:,:,k)=worktw(:,:,k)*Kv(k)  ! In the parabolic case Kv is already Kv/H 
#else
 	     worktw(:,:,k)=worktw(:,:,k)*KVERT*hhi(:,:)
#endif

#endif
        enddo
#endif   ! IMPLCITDIFF

! -------------------     
!  J_sigma=0 has been defined in 1)   

! worktu --> H dc/dx at point u
! worktv --> H dc/dy at point v 

      do k=1,NZ
         worktu(:,:,k)=worktu(:,:,k)*hhu(:,:)
c2d	     worktv(:,:,k)=worktv(:,:,k)*hhv(:,:)
      enddo

! ---------------
        do k=NZ,1,-1
c2d      do j=1,NY         
           do i=1,NX
              tt(i,j,k)= tt(i,j,k)  
     &	        + hhi(i,j)*(znu*DX3(worktu,i-1,j,k)
c2d     &       +           znu*DY3(worktv,i,j-1,k) 
#ifndef IMPLICITDIFF
     &	        +               DZW3(worktw,i,j,k-1))
#else
     &	        )
#endif
#ifdef RADIATION
     &          +  lightz(i,j,k) 
#endif

           enddo
c2d      enddo
        enddo
! ----------------
! Make sure temperature does not change on land!
! Please do not remove this. The radiation routines will evolve land points else, leading to ultimate disaster.
        FIXBC3D(tt,HGRID)

#endif


!  ****************************************************************
!  E)  Momentum equation for U
!  ****************************************************************

#ifdef FIXED_VEL_PROF
       ut=0
c2d       vt=0
	return  
#endif

!  ------------------------------------------------------------------ 
!   1) Advective part:  ADV = - u D_x(u) - v D_y(u)  - w/H D_w(u)  at point u 
!  -------------------------------------------------------------- 

! Velocities interpolated to help points     
         BLOCKIX3(u,workuhi,-1)        ! u velocity interpolated to point h, no BC for u*c_x since Boundary not on h 
c2d         BLOCKIX3(v,workvozi,0)  ! IX(v) in oz    
         BLOCKIX3(w,workwoyi,0)	! IX(w) in oy       
! -----------------------------------------------

! Derivatives of u
         BLOCKDX3(u,workuh,-1)           ! DX(u) in h    no BC

c2d      BLOCKDY3(u,workuoz,0) 		! DY(u) in o  
c2d      FIXBC3D(workuoz,OGRID)         ! Dy(u)=0 on boundaries - free slip
 
         BLOCKDZC3(u,workuoy,0) 		! DZ(u) oy 
#ifndef IMPLICITDIFF   
         workuoy(:,:,NZ)=dutop(:,:)              ! TOP BC for du/dz
#endif
! -----------------------------------------------
! Now multiply velocities and derivatives
         workuhi(:,:,:) =workuh(:,:,:) *workuhi(:,:,:)        ! u*Dx(u) in h
c2d      workvozi(:,:,:)=workuoz(:,:,:)*workvozi(:,:,:)  ! v*Dy(u) in oz
         workwoyi(:,:,:)=workuoy(:,:,:)*workwoyi(:,:,:)   ! w*Dz(u) in oy, BC is automatically imposed by w=0

! -----------------------------------------------
      do k=1,NZ
c2d      do j=1,NY
           do i=1,NX
            ut(i,j,k)= 
     &                 -IX3(workuhi,i,j,k)  !  DX(u) u  no need to fix, multiplied by u
     &       -hhui(i,j)*IZW3(workwoyi,i,j,k-1)
c2d     &              -IY3(workvozi,i,j-1,k)   
           enddo
c2d      enddo
      enddo
		
! -----------------------------------------------------------------------------------
!   2) Add pressure balance + horizontal and vertical diffusion to u tendency
!  D_t(u) = ADV - g D_x(h) + nu 1/H D_x ( H D_x (u) ) + nu 1/H D_y (H D_y(u) ) + D_w( Kv/H D_w(u) )
! -----------------------------------------------------------------------------------
#ifndef IMPLICITDIFF

#ifdef LMIXING
             BLOCKIX3(kvhmix,workuhi,0)
	         do k=0,NZ-1
                 workuoy(:,:,k)=workuoy(:,:,k)*workuhi(:,:,k)
             enddo  
             workuoy(:,:,NZ)=workuoy(:,:,NZ)*kvmin*hhui(:,:)
#else
!   multiply vertical help fields by Kv/H
      do k=0,NZ
#ifdef KVPAR
             workuoy(:,:,k)=workuoy(:,:,k)*Kv(k)                     ! in the parabolic case  Kv=diffusivity/H
#else
             workuoy(:,:,k)=workuoy(:,:,k)*KVERT*hhui(:,:)    ! hhui is 1/H
#endif
      enddo
#endif

! -----------------------------------------------
#endif ! IMPLICITDIFF
! -----------------------------------------------


! workuh -->  H du/dx
! workuoz --> H du/dy
      do k=1,NZ
         workuh(:,:,k)=workuh(:,:,k)*hh(:,:)*znuh(:,:)
c2d      workuoz(:,:,k)=workuoz(:,:,k)*hhoz(:,:)*znuoz(:,:)
      enddo
! -----------------------------------------------

! Compute density anomaly field and its derivatives
#ifdef TEMPERATURE
      BLOCKIX3(workbuoy,workbuoyi,0) ! Ix(b) in u
      BLOCKDX2(hh,uhi,0) ! Dx(H)
      duhi(:,:)=0.
#endif

      BLOCKDX2(h,vhi,0) ! Dx(h)

#ifdef TEMPERATURE
         duhi(:,:)=hh(:,:)*workbuoy(:,:,NZ)*(1-zlc(NZ))
#endif
! -----------------------------------------------
      do k=NZ,1,-1
c2d      do j=1,NY
           do i=1,NX
              ut(i,j,k)= ut(i,j,k)
     &              -GG*vhi(i,j)                            ! h_x is computed here, but used for u, which is 0 on the boundary anyway
     &              +hhui(i,j)*( DX3(workuh,i,j,k)
c2d     &           +            DY3(workuoz,i,j-1,k)
#ifndef IMPLICITDIFF
     &              + DZW3(workuoy,i,j,k-1) )
#else
     &               )
#endif
#ifdef TEMPERATURE
     & 		    + workbuoyi(i,j,k)*((zlc(k)-1)*uhi(i,j)+vhi(i,j))+DX2(duhi,i,j)		  
#endif
		  
           enddo
c2d      enddo

#ifdef TEMPERATURE
         duhi(:,:)=duhi(:,:)+hh(:,:)*workbuoy(:,:,k-1)*zcdif(k-1)
#endif
       enddo

!       print *,'Check',u(2,1,izoff(2,1)),u(2,1,izoff(2,1)+1)

#ifdef c2d

!  ****************************************************************
!  F)  Momentum equation for V
!  ****************************************************************

!  ------------------------------------------------------------------ 
!   1) Advective part:  ADV = - u D_x(v) - v D_y(v)  - w/H D_w(v) 
!  ------------------------------------------------------------------ 

! Velocities interpolated to help points
      BLOCKIY3(w,workwoxi,0)             ! IY(w) in ox    
      BLOCKIY3(u,workuozi,0)             ! IY(u) in oz    
c2d      BLOCKIY3(v,workvhi,-1) ! Velocity interpolated to  point h,   no BC for v*c_y since Boundary not on h 

!  ------------------------------------------------------------------ 
! Derivatives of v
      BLOCKDY3(v,workvh,-1)              ! DY(v) in h, no BC

      BLOCKDX3(v,workvoz,0)              ! DX(v) in o 
      FIXBC3D(workvoz,OGRID)            ! impose free slip BC    
    
      BLOCKDZC3(v,workvox,0)             ! DZ(v) ox  Bottom BC is guaranteed by ghost point
#ifndef IMPLICITDIFF
      workvox(:,:,NZ)=dvtop(:,:)        ! TOP BC for dv/dz
#endif
! -----------------------------------------------
! Now multiply velocities and derivatives 
      workuozi(:,:,:)=workvoz(:,:,:)*workuozi(:,:,:)   ! u*Dx(v) in oz
      workwoxi(:,:,:)=workvox(:,:,:)*workwoxi(:,:,:)   ! w*Dz(v) in ox, BC is automatically imposed by w
      workvhi(:,:,:) =workvh(:,:,:) *workvhi(:,:,:)    ! v*Dy(v) in h
! -----------------------------------------------
					     
      do k=1,NZ
         do j=1,NY
            do i=1,NX
             vt(i,j,k)= 
     &                 -IX3(workuozi,i-1,j,k)       
     &                 -IY3(workvhi,i,j,k)           
     &       -hhvi(i,j)*IZW3(workwoxi,i,j,k-1)                 
            enddo
         enddo
      enddo
! -----------------------------------------------------------------------------------
!   2) Add pressure balance + horizontal and vertical diffusion to v tendency
!  D_t(v) = ADV - g D_y(h) + nu 1/H D_x H D_x (v) + nu 1/H D_y H D_y (v)  + 1/H D_w( Kv/H D_w(v) ) 
! -----------------------------------------------------------------------------------

#ifndef IMPLICITDIFF

#ifdef LMIXING
             BLOCKIY3(kvhmix,workvhi,0)
	     do k=0,NZ-1
                 workvox(:,:,k)=workvox(:,:,k)*workvhi(:,:,k)
             enddo  
             workvox(:,:,NZ)=workvox(:,:,NZ)*kvmin*hhvi(:,:) 
#else
      do k=0,NZ
#ifdef KVPAR
             workvox(:,:,k)=workvox(:,:,k)*Kv(k) ! In the parabolic case Kv is already viscosity/H
#else
             workvox(:,:,k)=workvox(:,:,k)*KVERT*hhvi(:,:)
#endif
       enddo
#endif

#endif ! IMPLICITDIFF

! -----------------------------------------------
!    DZ(v) at the top is specified by the BCs already computed 

! workvh  --> H dv/dy
! workvoz --> H dv/dx
      do k=1,NZ
         workvh(:,:,k)=workvh(:,:,k)*hh(:,:)*znuh(:,:)
	 workvoz(:,:,k)=workvoz(:,:,k)*hhoz(:,:)*znuoz(:,:)
      enddo
! -----------------------------------------------

! Compute density anomaly field and its derivatives
#ifdef TEMPERATURE
      BLOCKIY3(workbuoy,workbuoyi,0) ! Iy(b) in u
      BLOCKDY2(hh,uhi,0) ! Dy(H)
      duhi(:,:)=0.
#endif

      BLOCKDY2(h,vhi,0) ! Dx(h)

#ifdef TEMPERATURE
         duhi(:,:)=hh(:,:)*workbuoy(:,:,NZ)*(1-zlc(NZ))
#endif
! -----------------------------------------------

      do k=NZ,1,-1
        do j=1,NY
          do i=1,NX
            vt(i,j,k)=vt(i,j,k)
     &           -GG*vhi(i,j) !     h_y is computed but used for v, so OK
     &           +hhvi(i,j)*( DX3(workvoz,i-1,j,k)
     &           +            DY3(workvh,i,j,k)
#ifndef IMPLICITDIFF
     &           +   DZW3(workvox,i,j,k-1)  )
#else
     &           )
#endif
#ifdef TEMPERATURE		  
     & 		    + workbuoyi(i,j,k)*((zlc(k)-1)*uhi(i,j)+vhi(i,j))+DY2(duhi,i,j)		  
#endif
          enddo
        enddo

#ifdef TEMPERATURE
         duhi(:,:)=duhi(:,:)+hh(:,:)*workbuoy(:,:,k-1)*zcdif(k-1)
#endif

      enddo

#endif



      end
! %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

#ifdef TIMESPLIT

#ifndef TEMPERATURE
      subroutine averagevert(u,v,um,vm)
#else
	  subroutine averagevert(u,v,t,um,vm,bmu,bmv,bmx,bmy)
#endif

!  Inputs: 
!            u,v    3d velocities
!            t temperature
!  Outputs:
!           um,vm   vertical averages of u and v
!           bm = int ( b sigma dsigma)

	  implicit none	
#include "commonlac.h"

   	  double um(RANGE2D)
      double vm(RANGE2D)
      double bmu(RANGE2D)
      double bmv(RANGE2D)
      double bmx(RANGE2D)
      double bmy(RANGE2D)
	  double bm(RANGE2D)

      double u(RANGE)
      double v(RANGE)
      double t(RANGE)
      double workbuoy(RANGE)

      integer k

! -----------------------------------
! Compute buoyancy from temperature
! -----------------------------------
#ifdef TEMPERATURE
      workbuoy(:,:,:)=GG*(1.-RHOEOS(t(:,:,:))/RHOWAT)  ! buoyancy
  
      bm(:,:)=0.
! Int(b*s ds)
      do k=1,NZ
		 bm(:,:)=bm(:,:)+workbuoy(:,:,k)*zlc(k)*zwdif(k)	  
      enddo

      BLOCKDX2(bm,bmx,0)
c2d      BLOCKDY2(bm,bmy,0)
	  BLOCKIX2(bm,bmu,0)
c2d      BLOCKIY2(bm,bmv,0)	
#endif

! Obtain barotropic um and vm

	  um(:,:)=0.
	  vm(:,:)=0.
	
	  do k=1,NZ
         um(:,:)=um(:,:)+u(:,:,k)*zwdif(k)
c2d         vm(:,:)=vm(:,:)+v(:,:,k)*zwdif(k)
	  enddo
		  		
      end

       subroutine computerslow(ut,vt,umt,vmt,rslowu,rslowv,taubu,taubv,tautu,tautv)

!  Inputs: 
!           ut,vt   3d vel  tendencies
!           umt,vmt the barotropic rhs
!  Outputs:
!         rslowum,rslowv the slow tendencies for the barotropic integration

	  implicit none	
#include "commonlac.h"

      double rslowu(RANGE2D)
      double rslowv(RANGE2D)
   	  double umt(RANGE2D)
      double vmt(RANGE2D)
   	  double workumt(RANGE2D)
      double workvmt(RANGE2D)
      double ut(RANGE)
      double vt(RANGE)
      double u(RANGE)
      double v(RANGE)
	  double h(RANGE2D)
	  double taubu(RANGE2D),taubv(RANGE2D)   ! Derivatives at bottom: Kv/H du/dz
	  double tautu(RANGE2D),tautv(RANGE2D)   ! Derivatives at top	
	
      integer k,i

! Obtain barotropic um and vm,  and the vertical integrals of ut and vt
	  workumt(:,:)=0.
c2d       workvmt(:,:)=0.  
	  do k=1,NZ
         workumt(:,:)=workumt(:,:)+ut(:,:,k)*zwdif(k)
c2d         workvmt(:,:)=workvmt(:,:)+vt(:,:,k)*zwdif(k)
	  enddo

!  the old rslows are already in the umts so we remove them first	
		  rslowu(:,:)=workumt(:,:) - ( umt(:,:) -rslowu(:,:) ) 
c2d  	  rslowv(:,:)=workvmt(:,:) - ( vmt(:,:) -rslowv(:,:) ) 
	
! Add Kv/H^2 du/ds at top - Kv/H^2 du/ds at bottom which are missing from ut, vt if we are using implicit
!			print *,taubut(32,1),taubu(32,1),rslowu(32,1)  ! JJNXX
#ifdef IMPLICITDIFF
		  rslowu(:,:)=rslowu(:,:) + ( tautu(:,:) - taubu(:,:) )
c2d  	  rslowv(:,:)=rslowv(:,:) + ( tautv(:,:) - taubv(:,:) )	
#endif

      end



! %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
! --------------------------------------------------------------

      subroutine timestepLFAM1fast(um,vm,h,umold,vmold,hold,umt,vmt,ht)
! --------------------------------------------------------------
 
! --------------------------------------------------------------
!  LF-AM scheme, first part, PREDICTOR: 
!  Entering: unow=u(n) uold=u(n-1) 
!  LEAP-FROG    u*(n+1)=2*ddt*u'(n)+u(n-1)  
!  We dont need the u(n-1) fields any more in old  
!  so that predicted  fields u*(n+1) were written in old


      implicit none
#include "commonlac.h"
! --------------------------------------------------------------
      double h(RANGE2D)
      double um(RANGE2D)
      double vm(RANGE2D)
      double hold(RANGE2D)
      double umold(RANGE2D)
      double vmold(RANGE2D)
      double umt(RANGE2D)
      double vmt(RANGE2D)
      double ht(RANGE2D)
! --------------------------------------------------------------
      integer i,j,k,l
      
     
      umold(:,:)=umold(:,:)+2*ddtfast*umt(:,:)
c2d      vmold(:,:)=vmold(:,:)+2*ddtfast*vmt(:,:)
      hold(:,:)=hold(:,:)+2*ddtfast*ht(:,:)
      
! feuler is a global flag in commonlac.h
      if(feulerfast) then            ! OOPS! We are just starting? Euler-AM then  
c     	 print *,'Euler step!'
         umold(:,:)=um(:,:)+ddtfast*umt(:,:)
c2d      vmold(:,:)=vm(:,:)+ddtfast*vmt(:,:)
         hold(:,:)=h(:,:)+ddtfast*ht(:,:)
      endif
! --------------------------------------------------------------

      feulerfast=.false.

      end


!%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
      subroutine timestepLFAM2fast(um,vm,h,umold,vmold,hold,umt,vmt,ht,umt1,vmt1,ht1,umt2,vmt2,ht2)
! --------------------------------------------------------------
!   LF-AM scheme, second part  CORRECTOR ADAMS-MOULTON
! --------------------------------------------------------------

      implicit none
#include "commonlac.h"

      double h(RANGE2D)
      double um(RANGE2D)
      double vm(RANGE2D)

      double hold(RANGE2D)
      double umold(RANGE2D)
      double vmold(RANGE2D)
      double ht(RANGE2D)
      double umt(RANGE2D)
      double vmt(RANGE2D)
      double ht1(RANGE2D)
      double umt1(RANGE2D)
      double vmt1(RANGE2D)
      double ht2(RANGE2D)
      double umt2(RANGE2D)
      double vmt2(RANGE2D)

      um(:,:)=umold(:,:)+ddtfast*(5*umt(:,:)+8*umt1(:,:)-umt2(:,:))/12
c2d      vm(:,:)=vmold(:,:)+ddtfast*(5*vmt(:,:)+8*vmt1(:,:)-vmt2(:,:))/12
      h(:,:)=hold(:,:)+ddtfast*(5*ht(:,:)+8*ht1(:,:)-ht2(:,:))/12
      end

#endif  ! TIMESPLIT

! %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%








! %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
! --------------------------------------------------------------
#ifdef TEMPERATURE
      subroutine timestepLFAM1(u,v,h,c,t,uold,vold,hold,cold,told,ut,vt,ht,ct,tt)
#else
      subroutine timestepLFAM1(u,v,h,c,uold,vold,hold,cold,ut,vt,ht,ct)
#endif
! --------------------------------------------------------------
 
! --------------------------------------------------------------
!  LF-AM scheme, first part, PREDICTOR: 
!  Entering: unow=u(n) uold=u(n-1) 
!  LEAP-FROG    u*(n+1)=2*ddt*u'(n)+u(n-1)  
!  We dont need the u(n-1) fields any more in old  
!  so that predicted  fields u*(n+1) were written in old


      implicit none
#include "commonlac.h"
! --------------------------------------------------------------
      double h(RANGE2D)
      double u(RANGE)
      double v(RANGE)
      double c(RANGE,NSCAL)

      double hold(RANGE2D)
      double uold(RANGE)
      double vold(RANGE)
      double cold(RANGE,NSCAL)

      double ut(RANGE)
      double vt(RANGE)
      double ht(RANGE2D)
      double ct(RANGE,NSCAL)

#ifdef TEMPERATURE
      double t(RANGE)
      double told(RANGE)
      double tt(RANGE)
#endif
! --------------------------------------------------------------
      integer i,j,k,l
 
      uold(:,:,:)=uold(:,:,:)+2*ddt*ut(:,:,:)
c2d      vold(:,:,:)=vold(:,:,:)+2*ddt*vt(:,:,:)

#ifndef TIMESPLIT
      hold(:,:)=hold(:,:)+2*ddt*ht(:,:)
#endif
      cold(:,:,:,:)=cold(:,:,:,:)+2*ddt*ct(:,:,:,:)

! feuler is a global flag in commonlac.h
      if(feuler) then            ! OOPS! We are just starting? Euler-AM then  
	 print *,'Euler step!'
         uold(:,:,:)=u(:,:,:)+ddt*ut(:,:,:)
c2d      vold(:,:,:)=v(:,:,:)+ddt*vt(:,:,:)
#ifndef TIMESPLIT
         hold(:,:)=h(:,:)+ddt*ht(:,:)
#endif
         cold(:,:,:,:)=c(:,:,:,:)+ddt*ct(:,:,:,:)
      endif
! --------------------------------------------------------------

#ifdef TEMPERATURE
      told(:,:,:)=told(:,:,:)+2*ddt*tt(:,:,:)
      if(feuler) then            ! OOPS! We are just starting? Euler-AM then  
         told(:,:,:)=t(:,:,:)+ddt*tt(:,:,:)
      endif
#endif

      feuler=.false.

      end


!%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
#ifdef TEMPERATURE
      subroutine timestepLFAM2(u,v,h,c,t,uold,vold,hold,cold,told,ut,vt,ht,ct,tt,ut1,vt1,ht1,ct1,tt1,ut2,vt2,ht2,ct2,tt2)
#else
      subroutine timestepLFAM2(u,v,h,c,uold,vold,hold,cold,ut,vt,ht,ct,ut1,vt1,ht1,ct1,ut2,vt2,ht2,ct2)
#endif
 
! --------------------------------------------------------------
!   LF-AM scheme, second part  CORRECTOR ADAMS-MOULTON
! --------------------------------------------------------------

      implicit none
#include "commonlac.h"

      double h(RANGE2D)
      double u(RANGE)
      double v(RANGE)
      double c(RANGE,NSCAL)

      double hold(RANGE2D)
      double uold(RANGE)
      double vold(RANGE)
      double cold(RANGE,NSCAL)

      double ht(RANGE2D)
      double ut(RANGE)
      double vt(RANGE)
      double ct(RANGE,NSCAL)
      double ht1(RANGE2D)
      double ut1(RANGE)
      double vt1(RANGE)
      double ct1(RANGE,NSCAL)
      double ht2(RANGE2D)
      double ut2(RANGE)
      double vt2(RANGE)
      double ct2(RANGE,NSCAL)

#ifdef TEMPERATURE
      double tt(RANGE)
      double tt1(RANGE)
      double tt2(RANGE)
      double t(RANGE)
      double told(RANGE)
#endif

		integer i,k

      u(:,:,:)=uold(:,:,:)+ddt*(5*ut(:,:,:)+8*ut1(:,:,:)-ut2(:,:,:))/12
c2d      v(:,:,:)=vold(:,:,:)+ddt*(5*vt(:,:,:)+8*vt1(:,:,:)-vt2(:,:,:))/12

#ifndef TIMESPLIT
      h(:,:)=hold(:,:)+ddt*(5*ht(:,:)+8*ht1(:,:)-ht2(:,:))/12
#endif
      c(:,:,:,:)=cold(:,:,:,:)+ddt*(5*ct(:,:,:,:)+8*ct1(:,:,:,:)-ct2(:,:,:,:))/12
#ifdef TEMPERATURE
      t(:,:,:)=told(:,:,:)+ddt*(5*tt(:,:,:)+8*tt1(:,:,:)-tt2(:,:,:))/12
#endif

      end

#if defined(TIMESPLIT) || defined(IMPLICITDIFF)
      subroutine computetaub(u,v,h,taub,taubu,taubv,tautu,tautv,hhi,hhui,hhvi)
! --------------------------------------------------------------
!   Computes help fields needed for implicit and timesplit
! --------------------------------------------------------------

      implicit none
#include "commonlac.h"

      double h(RANGE2D)
      double u(RANGE)
      double v(RANGE)
      
	  double hh(RANGE2D)		! H in h
	  double hhi(RANGE2D)       ! 1/H in h
	  double hhu(RANGE2D)		! H in a u point 
	  double hhui(RANGE2D)		! 1/H in u  
      double hhv(RANGE2D)            ! H in a  v point
      double hhvi(RANGE2D)           ! 1/H in  v

	  double taub(RANGE2D)		! bottom stresses for scalars
	  double taubu(RANGE2D)
	  double taubv(RANGE2D)
	  double tautu(RANGE2D)
	  double tautv(RANGE2D)
	
	  double coeff
	  integer i,k,l,j

! Computes total thickness  at h points
	     hh(:,:)=h(:,:)+z(:,:)
! Fix land u points to zmin;  Needed to allow computation of oz points 
c2d      hh(:,:)=hh(:,:)*llm(:,:,HGRID)+(1-llm(:,:,HGRID))*zmin

! Computes total thickness hh  at u points : hhu 
! Computes total thickness hh  at v points : hhv
	     BLOCKIX2(hh,hhu,0) 
c2d      BLOCKIY2(hh,hhv,0)

!  From now on we will also need 1/H (divisions are expensive to compute!)
	     hhui(0:NX,:)=1.d0/hhu(0:NX,:)
c2d      hhvi(:,0:NY)=1.d0/hhv(:,0:NY)
		 hhi(:,:)=1.d0/hh(:,:) 
		
! ------------- Compute bottom stress taub at point c for bottom BCs (used for scalars)  ------------------
! compute  du/dsigma in oy (plane of w)
! compute  dv/dsigma in ox (plane of w)

#ifdef IMMBOUND
	     j=1
c2d      do j=NYMIN,NYMAX
			do i=0,(NX+1)
				 taubu(i,j)=u(i,j,1+izoff(i,j,UGRID))/(zlc(izoff(i,j,UGRID)+1)-zlw(izoff(i,j,UGRID))) 
c2d              taubv(i,j)=v(i,j,1+izoff(i,j,VGRID))/(zlc(izoff(i,j,VGRID)+1)-zlw(izoff(i,j,VGRID)))  
	    	enddo
c2d     	enddo
#else
          taubu(:,:)=u(:,:,1)/zlc(1) 
c2d       taubv(:,:)=v(:,:,1)/zlc(1) 
#endif

! compute  the total shear amplitude in the the w-points 
		 j=1  ! necessary for the 2D case 
#ifndef FIXED_VEL_PROF
c2d      do j=1,NY
			 do i=1,NX
			    taub(i,j)=IX2(taubu,i-1,j)**2
c2d     &  +IY2(taubv,i,j-1)**2
			 enddo
c2d      enddo
#else  
		 taub(1,1)=taubu(1,1)**2 
#endif

! -----------------------------------------  
! compute  the stress in the  w-points 

#ifdef LMIXING
		      taub(:,:)=RHOWAT*kvmin*sqrt(taub(:,:))*hhi(:,:) !tau at bottom  Pacanovski gives kvmin
			  taubu(:,:)=kvmin*taubu(:,:)*hhi(:,:)*hhi(:,:) !taubu is  Kv/H du/dz
c2d		      taubv(:,:)=kvmin*taubv(:,:)*hhi(:,:)*hhi(:,:) !taubu is  Kv/H du/dz
#else
#ifdef KVPAR 
		      taub(:,:)=RHOWAT*Kv(0)*sqrt(taub(:,:))   ! PARABOLIC TSANIS 
		      taubu(:,:)=Kv(0)*(taubu(:,:))*hhi(:,:)   ! PARABOLIC TSANIS 
c2d		      taubv(:,:)=Kv(0)*(taubv(:,:))*hhi(:,:)   ! PARABOLIC TSANIS 

#else
		      taub(:,:)=RHOWAT*KVERT*sqrt(taub(:,:))*hhi(:,:)  ! CONSTANT DIFFUSION 
			  taubu(:,:)=KVERT*taubu(:,:)*hhi(:,:)*hhi(:,:)  !taubu is  Kv/H du/dz
c2d 	      taubv(:,:)=KVERT*taubv(:,:)*hhi(:,:)*hhi(:,:) 
#endif
#endif
			  tautu(:,:)=(TAUS/RHOWAT)*hhui(:,:)*cos(WINDRAD) 
c2d 	      tautv(:,:)=(TAUS/RHOWAT)*hhvi(:,:)*sin(WINDRAD)

		end

! ----------
#endif 


#ifdef IMPLICITDIFF

#define IVARU 1
#define IVARV 2
#define IVART 3

!%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
#ifdef TEMPERATURE
      subroutine solveimplicit(u,v,w,h,c,t,taub,hhi,hhui,hhvi,coeff)
#else
      subroutine solveimplicit(u,v,w,h,c,taub,hhi,hhui,hhvi,coeff)
#endif
 
! --------------------------------------------------------------
!   ADD IMPLICIT DIFFUSION
! --------------------------------------------------------------
! coeff is 5./12. for the corrector and 2 for the predictor

      implicit none
#include "commonlac.h"

      double h(RANGE2D)
      double u(RANGE)
      double v(RANGE)
      double w(RANGE)
      double c(RANGE,NSCAL)
#ifdef TEMPERATURE
	      double t(RANGE)
#endif

      double aa(RANGE)
      double bb(RANGE)
      double cc(RANGE)
	  double rr(RANGE)
	  double hhui(RANGE2D)       ! 1/H in h
	  double hhvi(RANGE2D)       ! 1/H in h      
	  double hhi(RANGE2D)       ! 1/H in h
	  double taub(RANGE2D)		! bottom stresses for scalars
#ifdef ECOSYSTEM
	  double rhat(RANGE,NSCAL)
#endif
	
	  double coeff
	  integer i,k,l,j
	
#ifdef TEMPERATURE
#ifdef ECOSYSTEM
#define MYARGS aa,bb,cc,rr,w,t,rhat,taub
#else
#define MYARGS aa,bb,cc,rr,w,t,taub
#endif
#else
#ifdef ECOSYSTEM
#define MYARGS aa,bb,cc,rr,w,rhat,taub
#else
#define MYARGS aa,bb,cc,rr,w,taub
#endif
#endif
	
#ifdef ECOSYSTEM	
	   call ecotendpatlhs(c,rhat,hhi) !
#endif

#ifdef TEMPERATURE
		call buildlhs(MYARGS,hhi,IVART,coeff)
		call solve_tridiag(aa,bb,cc,rr,t)	
#endif

	 	call buildlhs(MYARGS,hhui,IVARU,coeff)
		call solve_tridiag(aa,bb,cc,rr,u)
c2d	       call buildlhs(MYARGS,hhvi,IVARV,coeff)
c2d  	   call solve_tridiag(aa,bb,cc,rr,v)

		do l=1,NSCAL
			call buildlhs(MYARGS,hhi,IVART+l,coeff)
			call solve_tridiag(aa,bb,cc,rr,c(0,NYMIN,0,l))
	    enddo


      end

! ----------




#ifdef TEMPERATURE
#ifdef ECOSYSTEM
      subroutine buildlhs(aa,bb,cc,rr,w,t,rhat,taub,hhi,ivar,coeff)
#else
      subroutine buildlhs(aa,bb,cc,rr,w,t,taub,hhi,ivar,coeff)
#endif
#else
#ifdef ECOSYSTEM
      subroutine buildlhs(aa,bb,cc,rr,w,rhat,taub,hhi,ivar,coeff)
#else
      subroutine buildlhs(aa,bb,cc,rr,w,taub,hhi,ivar,coeff)
#endif
#endif

! Build LHS matrix for a variable (u,v) which lives on c levels and has no-stress bottom BCs
#include "commonlac.h"
#ifdef TEMPERATURE
      double t(RANGE)
#endif

	  double aa(RANGE)  ! lower tridiag
	  double bb(RANGE)  ! central tridiag
	  double cc(RANGE)  ! upper tridiag
	  double rr(RANGE)  ! additional rhs terms 
	  double w(RANGE)  
	  double dd(RANGE2D)  
	  double ee(RANGE2D)     ! used by scalars
	  double taub(RANGE2D)   ! used by scalars
	  double g6(RANGE2D)	 ! used by deposition for scalars
	  double g7(RANGE2D)
	
#ifdef ECOSYSTEM
	  double rhat(RANGE,NSCAL)   ! Rhat reaction terms Rhat*c+Rstar 
#endif
	
 	  double hhi(RANGE2D)
	  double phisurf(RANGE2D)
	  double cff,cwind2,coeff
	
	  integer k,ivar ! Which var are we considering ivar=1,2 (u,v) ivar=3 (t)

! ---------------------------

#ifdef IMMBOUND
! Used to shorten notation and improve readability. KKH is the new 'ghost' level (0 if no IMMBOUND)
#define KKH izoff(i,j,HGRID)
#define KKU izoff(i,j,UGRID)
#define KKV izoff(i,j,VGRID)
#define KVRANGE i,j
#else
#define KKH 0
#define KKU 0
#define KKV 0
#define KVRANGE :,:
#endif 

#ifdef LMIXING
!   kvhmix stands for Kv/H 
#define MYKVERTT(k) (kvhmix(KVRANGE,k)*cfDZC(k,2)) 
#define MYKVERTU(k) (kvhmixu(KVRANGE,k)*cfDZC(k,2))
#define MYKVERTV(k) (kvhmixv(KVRANGE,k)*cfDZC(k,2))

#ifdef IMMBOUND
#define MYKVERTCBOUND(k) (kvhmix(KVRANGE,k)/(zlc(k+1)-zlw(k))) 
#define MYKVERTUBOUND(k) (kvhmixu(KVRANGE,k)/(zlc(k+1)-zlw(k)))  
#define MYKVERTVBOUND(k) (kvhmixv(KVRANGE,k)/(zlc(k+1)-zlw(k)) )
#else 
#define MYKVERTCBOUND(k) (kvhmix(KVRANGE,k)*cfDZC(k,2)*2) 
#define MYKVERTUBOUND(k) (kvhmixu(KVRANGE,k)*cfDZC(k,2)*2)       
#define MYKVERTVBOUND(k) (kvhmixv(KVRANGE,k)*cfDZC(k,2)*2)        
#endif

#else
#ifdef KVPAR

#define MYKVERTT(k) (Kv(k)*cfDZC(k,2) ) 
#define MYKVERTU(k) (Kv(k)*cfDZC(k,2) )
#define MYKVERTV(k) (Kv(k)*cfDZC(k,2) )

#ifdef IMMBOUND 
#define MYKVERTCBOUND(k) (Kv(k)/(zlc(k+1)-zlw(k))) 
#define MYKVERTUBOUND(k) (Kv(k)/(zlc(k+1)-zlw(k)))  
#define MYKVERTVBOUND(k) (Kv(k)/(zlc(k+1)-zlw(k))) 
#else
#define MYKVERTCBOUND(k) (Kv(k)*cfDZC(k,2)*2)
#define MYKVERTUBOUND(k) (Kv(k)*cfDZC(k,2)*2)       
#define MYKVERTVBOUND(k) (Kv(k)*cfDZC(k,2)*2)        
#endif

#else

#define MYKVERTT(k) (KVERT*hhi(KVRANGE)*cfDZC(k,2))    
#define MYKVERTU(k) (KVERT*hhi(KVRANGE)*cfDZC(k,2)) 
#define MYKVERTV(k) (KVERT*hhi(KVRANGE)*cfDZC(k,2))

#ifdef IMMBOUND
#define MYKVERTCBOUND(k) (KVERT*hhi(KVRANGE)/(zlc(k+1)-zlw(k)))
#define MYKVERTUBOUND(k) (KVERT*hhi(KVRANGE)/(zlc(k+1)-zlw(k)))  
#define MYKVERTVBOUND(k) (KVERT*hhi(KVRANGE)/(zlc(k+1)-zlw(k))) 
#else
#define MYKVERTCBOUND(k) (KVERT*hhi(KVRANGE)*cfDZC(k,2)*2) 
#define MYKVERTUBOUND(k) (KVERT*hhi(KVRANGE)*cfDZC(k,2)*2)       
#define MYKVERTVBOUND(k) (KVERT*hhi(KVRANGE)*cfDZC(k,2)*2)        
#endif

#endif
#endif

#define MYCFF(k) cff*cfDZW(k,2)*hhi(KVRANGE)
		
! --------------------------
	  cff=-coeff*ddt  ! coeff is 5./12 for corrector and 2 for predictor (NEGATIVE!)
	  rr(:,:,:)=0.
! ---------------------------------
	  if((ivar==IVART)) then
	
#ifdef TEMPERATURE
#ifdef RADIATION
	     call radiationsurf(t(0,NYMIN,NZ),phisurf) 
#else	
         phisurf(:,:)=hsurf*(RHOAIR*CAIR)*(Tair-t(:,:,NZ))  ! Top boundary condition [ Flux=-h(Tair-T) ] 
	                                                       ! Using top level (error is less than 2%)   
#endif
#endif	
		 
		 j=1
		 do k=2,NZ-1 
#ifdef IMMBOUND	
c2d      do j=NYMIN,NYMAX
			do i=0,(NX+1)
#endif
		  	cc(KVRANGE,k)=MYCFF(k-1)*MYKVERTT(k)
		  	bb(KVRANGE,k)=1+MYCFF(k-1)*(-MYKVERTT(k)-MYKVERTT(k-1))
		  	aa(KVRANGE,k)=MYCFF(k-1)*MYKVERTT(k-1)
#ifdef IMMBOUND
			enddo
c2d 	 enddo
#endif		
		 enddo

#ifdef IMMBOUND	
c2d      do j=NYMIN,NYMAX
         do i=0,(NX+1)
#endif
	     cc(KVRANGE,KKH+1)=MYCFF(KKH)*MYKVERTT(KKH+1)
	     bb(KVRANGE,KKH+1)=1+MYCFF(KKH)*(-MYKVERTT(KKH+1))
         aa(KVRANGE,KKH+1)=0.
 		 cc(KVRANGE,0:KKH)=0.
		 bb(KVRANGE,0:KKH)=1.
		 aa(KVRANGE,0:KKH)=0.

	  	 cc(KVRANGE,NZ)=0.0
		 bb(KVRANGE,NZ)=1+MYCFF(NZ-1)*(-MYKVERTT(NZ-1))
  		 aa(KVRANGE,NZ)=MYCFF(NZ-1)*MYKVERTT(NZ-1)
  	
 		 rr(KVRANGE,NZ)=-MYCFF(NZ-1)*phisurf(KVRANGE)*(1./(RHOWAT*CWAT))
#ifdef IMMBOUND
		 enddo
c2d		 enddo
#endif
       endif

! ---------------------------------
	  if(ivar==IVARU) then

         cwind2=cos(WINDRAD) !  WINDRAD is defined above
	   
         j=1
		 do k=2,NZ-1 
#ifdef IMMBOUND	
c2d      do j=NYMIN,NYMAX
			do i=0,(NX+1)
#endif
		  	cc(KVRANGE,k)=MYCFF(k-1)*MYKVERTU(k)
		  	bb(KVRANGE,k)=1+MYCFF(k-1)*(-MYKVERTU(k)-MYKVERTU(k-1))
		  	aa(KVRANGE,k)=MYCFF(k-1)*MYKVERTU(k-1)
#ifdef IMMBOUND
		 enddo
c2d 	 enddo
#endif
	     enddo

#ifdef IMMBOUND	
c2d      do j=NYMIN,NYMAX
         do i=0,(NX+1)
#endif
!#ifndef DEBUGDIFF2		
 	     cc(KVRANGE,KKU+1)=MYCFF(KKU)*MYKVERTU(KKU+1)
 		 bb(KVRANGE,KKU+1)=1+MYCFF(KKU)*(-MYKVERTU(KKU+1)-MYKVERTUBOUND(KKU))  ! JJN3
 		 cc(KVRANGE,0:KKU)=0.
		 bb(KVRANGE,0:KKU)=1.
		
		 aa(KVRANGE,0:KKU)=0.
				 
		! zlc(KKU+1)-zlw(KKU)
!#else
!     Alternative, 2nd order formulation
! 	     cc(KVRANGE,KKU+1)=cff*cfDZW(KKU,2)*(MYKVERTU(KKU+1)+MYKVERTU(KKU)*zlc(1)/zlc(2))
! 		 bb(KVRANGE,KKU+1)=1+cff*cfDZW(KKU,2)*(-MYKVERTU(KKU+1)-MYKVERTU(KKU)*zlc(2)/zlc(1)) 
!#endif
 		 aa(KVRANGE,KKU+1)=0.
	
	     cc(KVRANGE,NZ)=0.0	
	     bb(KVRANGE,NZ)=1+MYCFF(NZ-1)*(-MYKVERTU(NZ-1))
		 aa(KVRANGE,NZ)=MYCFF(NZ-1)*MYKVERTU(NZ-1)	
		
		 rr(KVRANGE,NZ)=-MYCFF(NZ-1)*(TAUS/RHOWAT)*cwind2

#ifdef IMMBOUND
		 enddo
c2d		 enddo
#endif

  	  endif
! ---------------------------------

#ifdef C2D
	  if(ivar==IVARV) then	

	    cwind2=sin(WINDRAD) !  WINDRAD is defined above
		
		j=1
	    do k=2,NZ-1 
#ifdef IMMBOUND	
c2d      do j=NYMIN,NYMAX
	     do i=0,(NX+1)
#endif		
		  	cc(KVRANGE,k)=MYCFF(k-1)*MYKVERTV(k)
		  	bb(KVRANGE,k)=1+MYCFF(k-1)*(-MYKVERTV(k)-MYKVERTV(k-1))
		  	aa(KVRANGE,k)=MYCFF(k-1)*MYKVERTV(k-1)
#ifdef IMMBOUND
  	     enddo
c2d		 enddo
#endif
   	    enddo

#ifdef IMMBOUND
c2d      do j=NYMIN,NYMAX
         do i=0,(NX+1)
#endif
!#ifndef DEBUGDIFF2	
		 cc(KVRANGE,KKV+1)=MYCFF(KKV)*MYKVERTV(KKV+1)
		 bb(KVRANGE,KKV+1)=1+MYCFF(KKV)*(-MYKVERTV(KKV+1)-MYKVERTVBOUND(KKV))
		 cc(KVRANGE,0:KKV)=0.
		 bb(KVRANGE,0:KKV)=1.
		 aa(KVRANGE,0:KKV)=0.
!#else
!     Alternative, 2nd order formulation   JJN3
!		 cc(i,j,KKV+1)=cff*cfDZW(KKV,2)*(MYKVERTV(KKV+1)+MYKVERTV(KKV)*zlc(1)/zlc(2))
!		 bb(i,j,KKV+1)=1+cff*cfDZW(KKV,2)*(-MYKVERTV(KKV+1)-MYKVERTV(KKV)*zlc(2)/zlc(1)) 
!#endif	
		 aa(i,j,KKV+1)=0.
		
	  	 cc(KVRANGE,NZ)=0.0
		 bb(KVRANGE,NZ)=1+MYCFF(NZ-1)*(-MYKVERTV(NZ-1))
  		 aa(KVRANGE,NZ)=MYCFF(NZ-1)*MYKVERTV(NZ-1)
                 rr(KVRANGE,NZ)=-MYCFF(NZ-1)*(TAUS/RHOWAT)*cwind2
#ifdef IMMBOUND
		 enddo
c2d		 enddo
#endif

  	  endif
#endif

!  SCALARS


	  if(ivar>IVART) then
		
		l=ivar-IVART	! this is the scalar type number (1...)	
				
		! %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%% SCALARS %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

		! ---------------------------  
		!      J_sigma=  dd*c+ee 
		!   ee should contain all stuff which is not multiplied by c
		!   dd is everything which is multiplied by c
		!-------------------------------------------------------------
!         print *,taub(16,1),tauce
		! ---------   Erosion   ------------------  
		    	ee(:,:)=betac(l)*erom*max(taub(:,:)/tauce-1,0.d0) 
c	    if(l.eq.SCAL_SEDIMENT) then
c			print *,'EROM',l,betac(l),erom,taub(250,1),ee(250,1)
c	    endif

		! -----------------  Deposition   --------------- 
		! Bottom boundary conditions for scalars

			 select case (l)

		!-------------------------------------------------------------
			   case (SCAL_SEDIMENT,SCAL_PHYTO1,SCAL_PHYTO2,SCAL_DETRITUS)
		! -------- J_s=E+D   with  D= ws *c(0)*(1-|taub|/taucd-1) if taub>taucd here D=c(1/2)*dd
		
			      dd(:,:)=min(wsett(l),0.d0)*max(1-taub(:,:)/taucd(l),0.d0)
!			 print *,'dd=',dd(16,1),wsett(l),l
!			 print *,'taub',taub(16,1),taucd(l)
		!-------------------------------------------------------------
			   case (SCAL_NUTRIENT)
		! -------- J_s= cin +E 
		! Just a flux now
			      dd(:,:)=0.d0	    
			      ee(:,:)=cinnut+ee(:,:)

		! Old form J_s=h(cin-c)+E 
		!	      dd(:,:)=-hnut
		!	      ee(:,:)=ee(:,:)+hnut*cinnut
		!-------------------------------------------------------------
			   case (SCAL_ZOO)
		! -------- J_s=ee(:,:)

			      dd(:,:)=0.d0
		!-------------------------------------------------------------
			   case (SCAL_NEUTRAL)
		! -------- Js=0 a pure passive scalar

			      dd(:,:)=0.d0
			      ee(:,:)=0.d0

			    case default
		   	      print *,'Unknown scalar type in BCs!'
			      stop

			 end select

			g6(:,:) = 1./(1.+(wsett(l)-dd(:,:))/MYKVERTCBOUND(KKH))
			g7(:,:) = ee(:,:)/(MYKVERTCBOUND(KKH)+(wsett(l)-dd(:,:)))

		j=1
		do k=2,NZ-1 
#ifdef IMMBOUND	
c2d      do j=NYMIN,NYMAX
			do i=0,(NX+1)
#endif
		 	cc(KVRANGE,k)=MYCFF(k-1)*(MYKVERTT(k) - (w(KVRANGE,k)+wsett(l))*cfIZC(k,2))	
		 
		 	bb(KVRANGE,k)=1+MYCFF(k-1)*(-MYKVERTT(k)-MYKVERTT(k-1)
     &         -(w(KVRANGE,k)+wsett(l))*cfIZC(k,1)+(w(KVRANGE,k-1)+wsett(l))*cfIZC(k-1,2) )
#ifdef ECOSYSTEM
     &         + cff*rhat(KVRANGE,k,l)
#endif
		  	aa(KVRANGE,k)=MYCFF(k-1)*(MYKVERTT(k-1)  +	(w(KVRANGE,k-1)+wsett(l))*cfIZC(k-1,1))
! cfIZC helps to interpolate from c to w
#ifdef IMMBOUND
			enddo
c2d 	 enddo
#endif		
		enddo

! JJNXXXXX

#ifdef IMMBOUND	
c2d      do j=NYMIN,NYMAX
	        do i=0,(NX+1)
#endif

		    cc(KVRANGE,KKH+1)=MYCFF(KKH)*(MYKVERTT(KKH+1) - (w(KVRANGE,KKH+1)+wsett(l))*cfIZC(KKH+1,2))	

		    bb(KVRANGE,KKH+1)=1+MYCFF(KKH)*(-MYKVERTT(KKH+1)
     &         -(w(KVRANGE,KKH+1) + wsett(l))*cfIZC(KKH+1,1) !JJN
     &         +dd(KVRANGE)*g6(KVRANGE) )  
#ifdef ECOSYSTEM
     &         + cff*rhat(KVRANGE,KKH+1,l)
#endif
		aa(KVRANGE,KKH+1)=0.
	        cc(KVRANGE,0:KKH)=0.
		bb(KVRANGE,0:KKH)=1.
#ifdef ECOSYSTEM
     &                + cff*rhat(KVRANGE,(KKH+1):(KKH+1),l)   ! XXXXXJvH GIUSTO?????
#endif
		    aa(KVRANGE,0:KKH)=0.  ! JJN mancano D ed E
c	        rr(KVRANGE,KKH+1)=-MYCFF(KKH)*(ee(KVRANGE)+dd(KVRANGE)*g7(KVRANGE))  !JJ13 This was the mistake
	        rr(KVRANGE,KKH+1)=-cff*cfDZW(KKH,2)*(ee(KVRANGE)+dd(KVRANGE)*g7(KVRANGE))  !JJ13

			cc(KVRANGE,NZ)=0.0
		    bb(KVRANGE,NZ)=1+MYCFF(NZ-1)*(-MYKVERTT(NZ-1)
     &         + (w(KVRANGE,NZ-1)+wsett(l))*cfIZC(NZ-1,2) )    ! JJNXX CHECK SIGN of w
#ifdef ECOSYSTEM
     &                         + cff*rhat(KVRANGE,NZ,l)
#endif
		    aa(KVRANGE,NZ)=MYCFF(NZ-1)*(MYKVERTT(NZ-1)+(w(KVRANGE,NZ-1)+wsett(l))*cfIZC(NZ-1,1))  !JJN

		if (l==SCAL_NUTRIENT) then
#ifdef NUTCYCLE
			rr(KVRANGE,NZ)=cff*cfDZW(NZ-1,2)*ctopnut*NUTCYCLE(ntime,ddt)  ! JJN check sign
#else
			rr(KVRANGE,NZ)=cff*cfDZW(NZ-1,2)*ctopnut  ! JJN check sign
#endif
        else
			rr(KVRANGE,NZ)=0.
		endif
#ifdef IMMBOUND
		    enddo
c2d		 enddo
#endif

   		endif ! l>IVART

  	  end

#endif ! IMPLICITDIFF


!%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
      subroutine checkCFL(u,v,w,h,c)

! --------------------------------------------------------------
! Checks CFL and prints out some diagnostics
! --------------------------------------------------------------

      implicit none
#include "commonlac.h"

      double u(RANGE)
      double v(RANGE)
      double w(RANGE)
      double h(RANGE2D)
      double c(RANGE,NSCAL)

      double cct(NSCAL),vtot
      integer i,j,k,l
      double cflx,cfly,cflz,umax,vmax,wmax,cff
      double uu,vv,ww,rrx,rry,rrz,cc

      vmax=0.0
      umax=0.0
      wmax=0.0

      j=1
      cflz=0.
      do k=1,NZ
c2d         do j=1,NY
      do i=1,NX
	    uu=u(i,j,k)*u(i,j,k)
c2d	    vv=v(i,j,k)*v(i,j,k)
	    ww=w(i,j,k)*w(i,j,k)/(zcdif(k)**2)
            if (uu.gt.umax) umax=uu
c2d	    if (vv.gt.vmax) vmax=vv
            if (ww.gt.cflz) cflz=ww
      enddo
c2d         enddo
      enddo
c2d      vmax=sqrt(vmax)
      umax=sqrt(umax)
      cflz=sqrt(cflz)*ddt
      cflx = ddt * umax / ddx 		! T * V / L
c2d      cfly = ddt * vmax / ddy 
    
      rrx = ddt * znu / (ddx*ddx)
      rry = ddt * znu / (ddy*ddy)
      rrz=0.0
      do k=0,NZ
	  cff=(ddt * KVERT / (zcdif(k)*zcdif(k))) /DEPTH/DEPTH
	  if(rrz.lt.cff) then	    
	    rrz = cff
	  endif
      enddo
	
#ifdef c2d
      write(*,"(i10,e16.7,a6,e16.7,a6,e16.7,a6,e16.7,a6,f8.5,a6,
     &          f8.5,a6,f8.5,a6,f8.5,a6,f8.5)")
c     &       ntime,ntime*ddt," CFLX=",cflx," CFLY=",cfly,
c     &       " Umax=",umax," Vmax=",vmax 
     &       ntime,ntime*ddt," CFLX=",cflx," CFLY=",cfly,
     &       " CFLZ=",cflz,
     &       " Umax=",umax," Vmax",vmax, 
     &       " RRX=",rrx,  " RRY=",rry, " RRZ=",rrz
#else
          write(*,"(i10,e16.7,a6,e16.7,a6,f8.5,a6,f8.5,a6,f8.5,a6,f8.5)") 
     &       ntime,ntime*ddt," CFLX=",cflx," RRX= ",rrx," Umax=",umax ,
     &        " CFLZ=",cflz," RRZ= ",rrz
#endif

#ifdef CHECKSCAL
      open(unit=57,file='scalstats.dat',access='append')
      vtot=sum(llm(:,:,HGRID)*z(:,:))*ddx*ddy
	vtot=1 ! JJ13
      do l=1,NSCAL
	cct(l)=0.d0
        do k=1,NZ
	  cc=sum(   c(:,:,k,l)*llm(:,:,HGRID)   )
#ifdef TANPAR
	  cct(l)=cct(l)+cc*zwdif(k)
#else
	  cct(l)=cct(l)+cc/NZ
#endif
        enddo
	cct(l)=cct(l)*ddx*ddy
      enddo
      write(57,'(i8,7e16.7)') ntime,ntime*ddt,(cct(l)/vtot,l=1,NSCAL)
!      write(57,*) ntime,ntime*ddt,(cct(l)/vtot,l=1,NSCAL)
      close(57)
#endif
      end


!%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
#ifdef TEMPERATURE
      subroutine saveoutput(fout,u,v,w,h,c,t)
#else
      subroutine saveoutput(fout,u,v,w,h,c)
#endif

! ------------------------------------------------------------------
!  Writes to output file outfile
! ------------------------------------------------------------------

      implicit none
#include "commonlac.h"

      double h(RANGE2D)
      double u(RANGE)
      double v(RANGE)
      double w(RANGE)
      double c(RANGE,NSCAL)
#ifdef TEMPERATURE
      double t(RANGE)
#endif
      double hh(RANGE2D)

      double uum(RANGE)
c2d      double vvm(RANGE)  !   

      double hui(RANGE2D)               ! DX H in a u point 
c2d      double hvi(RANGE2D)            ! DY H in a  v point
      double zui(RANGE2D)               ! DX Z in a u point 
c2d      double zvi(RANGE2D)            ! DY Z in a  v point

      double uumi(RANGE)               ! work array in a u point 
c2d      double vvmi(RANGE)            ! work array in a  v point

      double zzu,zzw,uth,hhu,usum,uth0,uth1
      double q1,q2,p1,p2,aa,bb
      integer i,j,k,fout,isec,l,jsec
      double um,vm,cm(NSCAL),err(0:NX),wm,wmreal

      j=1

      if(mod(ntime,nsave).eq.0.0) then 

! ###### Added to have correct w output

         BLOCKDX2(h,hui,0)
c2d      BLOCKDY2(h,hvi,0)

         BLOCKDX2(z,zui,0)
c2d      BLOCKDY2(z,zvi,0)

  
! Compute vertical average of additional terms (eq 8) and interpolate back       
    
         do k=1+izoff(i,j,HGRID),NZ
            uum(:,:,k)=u(:,:,k)*((zui(:,:)+hui(:,:))*(zlc(k)-1)+hui(:,:))
c2d         vvm(:,:,k)=v(:,:,k)*((zvi(:,:)+hvi(:,:))*(zlc(k)-1)+hvi(:,:))
         enddo

c	print *,">-",uum(128,1,26)

         BLOCKIX3(uum,uumi,-1)
c2d        BLOCKIY3(vvm,vvmi,-1)



c	print *,">>",u(128,1,26),uum(128,1,26),uumi(128,1,26),zui(128,1)
         print *, "Saving ",ntime
#ifdef TEMPERATURE
         write(fout) h,u,v,w,c,t
#else
         write(fout) h,u,v,w,c
#endif

	open(unit=50,file='restart.unf',form='unformatted')
	write(50) ntime
#ifdef TEMPERATURE
         write(50) h,u,v,c,t
#else
         write(50) h,u,v,c
#endif
	close(50)

	print *,'written'

! corrected !!!!

      p1=LAMBDA*ZSH
      p2=LAMBDA*(ZBH+1)
      q1=(1+ZBH)*log(1.d0+1.d0/ZBH)-1
      q2=ZSH*log(1.d0+1.d0/ZSH)-1
  
      aa=q2/(p1*q2+q1*p2)
      bb=-q1/(p1*q2+q1*p2)

      hh(:,:)=h(:,:)+z(:,:)
#ifndef c2d
      hhu=IX2(hh,NX/2,1) 
#else
      hhu=IY2(hh,NX/2+1,NY/2) 
#endif   

#ifdef FIXED_VEL_PROF
      hhu=0
#endif
!
!      print *,RHOWAT*USTARS**2*abs(LAMBDA*(aa*(1+zsh)-bb*zbh)),USTARS,uwind
!      stop
!      print *,hhu,aa,bb 
      usum=0

	isec=NX/2
	jsec=NY/2+1

      open(unit=30,file='uprof.dat')
      do k=0,NZ
	  zzu=zlc(k)*hhu
	  
#ifndef KVPAR
	  uth=(TAUS/(4.d0*RHOWAT*KVERT*DEPTH))*zzu*(3*zzu-2*DEPTH)
#else
	  uth=((aa*log(1+zzu/ZBH/hhu)+
     &     bb*log(1.-(zzu/(ZSH+1)/hhu))))*USTARS
#endif 
! Fix Tsanis velocity profiles
#ifdef FIXED_VEL_PROF
	  u(:,:,k)=uth
c2d	  v(:,:,k)=0.d0
#endif
	  usum=usum+u(isec,jsec,k)
	  write(30,'(7e16.7)') zlw(k),w(isec,jsec,k), zlc(k),
     &     u(isec,jsec,k)/USTARS,
     &     uth/USTARS,
     &     v(isec,jsec,k)/USTARS,USTARS

!	  write(70,'(2e16.7)')  zlc(k),  c(1,1,k,1)

      enddo
	write(30,*) 
      close(30)
      close(70)

	open(unit=33,file='errprof.dat')
	do i=0,NX
	err(i)=0
      do k=1,NZ
	  zzu=zlc(k)*hhu

#ifndef KVPAR
	  uth=(TAUS/(4.d0*RHOWAT*KVERT*DEPTH))*zzu*(3*zzu-2*DEPTH)
#else
	  uth=((aa*log(1+zzu/ZBH/hhu)+
     &     bb*log(1-(zzu/(ZSH+1)/hhu))))*USTARS
#endif
	  err(i)=err(i)+(uth-u(i,NY/2+1,k))**2*zwdif(k)
        enddo
        
	write(33,*) ddx*i,sqrt(err(i))
 	enddo
       close(33)

	open(unit=33,file='usurf_noint.dat')
c2d      do j=1,NY
	do i=1,NX
	    write(33,*) u(i,j,NZ),v(i,j,NZ)
	enddo
c2d	 enddo
        close(33)

	vm=0.d0

        open(unit=33,file='cprof.dat')
	do k=1,NZ
#ifdef TEMPERATURE
                write(33,'(9e16.7)') zlc(k),(c(1,1,k,l)/(h(1,1)+z(1,1)),l=1,NSCAL),t(1,1,k)
#else
                write(33,'(8e16.7)') zlc(k),(c(1,1,k,l)/(h(1,1)+z(1,1)),l=1,NSCAL)
#endif
           enddo
        close(33)

        open(unit=33,file='u_int.dat')
	do k=1,NZ
c2d      do j=1,NY
           do i=1,NX
                um=IX3(u,i-1,j,k)
		wm=IZW3(w,i,j,k-1)
c2d                vm=IY3(v,i,j-1,k)
                wmreal=wm+uumi(i,j,k)
c2d                wmreal=wmreal+vvmi(i,j,k)
#ifdef TEMPERATURE
                write(33,'(14e16.7)')  ddx*(i-0.5),zlc(k),um,vm,wm,wmreal,
     &            (c(i,j,k,l)/(h(i,j)+z(i,j)),l=1,NSCAL),t(i,j,k)
#else
                write(33,'(13e16.7)')  ddx*(i-0.5),zlc(k),um,vm,wm,wmreal,
     &             (c(i,j,k,l)/(h(i,j)+z(i,j)),l=1,NSCAL)
#endif
           enddo
c2d      enddo
	enddo
        close(33)

        open(unit=33,file='u_noint.dat')
        do k=1,NZ
c2d      do j=0,NY
           do i=0,NX
               
                write(33,'(9e16.7)') 
     &           u(i,j,k),v(i,j,k),w(i,j,k),(c(i,j,k,l)/(h(i,j)+z(i,j)),l=1,NSCAL)
           enddo
c2d      enddo
        enddo
        close(33)

	open(unit=31,file='uave_int.dat')
c2d      do j=1,NY
	do i=1,NX
	    um=0
	    vm=0
	    wm=0
            wmreal=0
	    do k=izoff(i,j,HGRID)+1,NZ
		um=um+IX3(u,i-1,j,k)*zwdif(k)
c2d		vm=vm+IY3(v,i,j-1,k)*zwdif(k)
		wm=wm+w(i,j,k)*zcdif(k)

                wmreal=wmreal+uumi(i,j,k)*zwdif(k)
c2d             wmreal=wmreal+vvmi(i,j,k)*zwdif(k)
	    enddo
!         	print *,wm,wmreal,uumi(i,j,16),u(i,j,16)
            wmreal=wmreal+wm
	    write(31,'(6e16.7)') ddx*(i-0.5),um,vm,wm,h(i,j),wmreal
	enddo
c2d	 enddo
        close(31)

        open(unit=34,file='uave_nonint.dat')
c2d      do j=1,NY
	do i=1,NX
	    um=0
	    vm=0
	    do k=1,NZ
		um=um+u(i,j,k)*zwdif(k)
		vm=vm+v(i,j,k)*zwdif(k)
	    enddo
	    write(34,'(3e16.7)') um,vm,h(i,j)
	enddo
c2d	 enddo
        close(34)


	open(unit=31,file='cave.dat')
c2d      do j=1,NY
	do i=1,NX
	  do l=1,NSCAL
	    cm(l)=0
	    do k=1,NZ
		cm(l)=cm(l)+c(i,j,k,l)/(h(i,j)+z(i,j))*zwdif(k)
	    enddo
          enddo
	  write(31,'(6e16.7)') (cm(l),l=1,NSCAL)
	enddo
c2d	 enddo
        close(31)

	endif
      end



!%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%% 
! ------------------------------------------------------------------------
          subroutine initparameters()
! ------------------------------------------------------------------------
! ------------------------------------------------------------------------
! Defines grids, vertical distribution of viscosity and other parameters
! ------------------------------------------------------------------------

      implicit none
#include "commonlac.h"     
      integer k
      double zz,zzu,hhu,d0,d1,d2,d3

! --------------------------
!   Numerical parameters
! --------------------------
!rescaling of ddt with NNT
      nsave=nsave*NNT
      numtime=numtime*NNT
      ncfl=ncfl*NNT
      ddt=ddt/real(NNT)
#ifdef TIMESPLIT
      ddtfast=ddt/NFASTSTEPS
#else 
	  ddtfast=ddt
#endif
! --------------------------
      ddx=DOMAINX/real(NX)
      ddy=DOMAINY/real(NY)
      ddsigma=1.D0/real(NZ)
! -----------------------------
! Coefficients for derivatives 
! -----------------------------
         cffx=1./ddx
         cffsigma=1./ddsigma
c2d      cffy=1./ddy
! ------------------------------------------------------------------------------
!   zlc(k) vertical location of  c,u,v  fields
!   k=0    ghost point below ground level
!   k=1    lives starting 1/2 step from the bottom
!   k=NZ   just below the surface 
! --------------------------------------------------------------------------------
!   zlw(k)  vertical location of  w fields     
!   k=0     starts at   the bottom
!   k=NZ    surface 
! --------------------------------------------------------------------------------



#ifdef TANPAR

! Non-regular vertical grids   
      do k=1,NZ	
      zlc(k)=(tanh(((k-0.5)*ddsigma-0.5)/TANPAR)/tanh(0.5/TANPAR)+1)*0.5
      enddo
      zlc(0)=-zlc(1) ! Built to be symmetric 
! ------------------------------------------------------------------------------
      do k=0,NZ	
      zlw(k)=(tanh((k*ddsigma-0.5)/TANPAR)/tanh(0.5/TANPAR)+1)*0.5
      enddo
!   To make sure for precision 
      zlw(NZ)=1.d0   

#ifdef DEBUGTANPAR
		do k=1,NZ	
  				zlc(k)=(zlw(k)+zlw(k-1))/2
  	    enddo
        zlc(0)=-zlc(1)
#endif

#else

! Prepare regular vertical grid
      do k=0,NZ	
          zlc(k)=(k-0.5)*ddsigma
          zlw(k)=k*ddsigma
      enddo
#endif
 
! cfDZW and cfDZC are also used by implicit, so better compute them anyway!
! -----------------------------------------------------------------------------------
! Coefficients to differentiate from the w(k) and w(k+1) levels on the c(k+1) level 
! ----------------------------------------------------------------------------------- 
      do k=0,NZ-1	
          d1=zlc(k+1)-zlw(k)
	  d2=zlw(k+1)-zlc(k+1)

          cfDZW(k,1)=-1/(d1+d2)
	  cfDZW(k,2)=1/(d1+d2)
      enddo
! -----------------------------------------------------------------------------------
! Coefficients to differentiate from  the   c(k) and c(k+1) levels on the w(k) level  
! -----------------------------------------------------------------------------------
      do k=0,NZ-1
          d1=zlw(k)-zlc(k)
	  d2=zlc(k+1)-zlw(k)

          cfDZC(k,1)=-1/(d1+d2)
	  cfDZC(k,2)=1/(d1+d2)
      enddo
! -------------------------------------------------------------------------
! Linear interpolation from the w(k) and w(k+1) levels on the c(k+1) level   
! -------------------------------------------------------------------------
      do k=0,NZ-1
          d1=zlc(k+1)-zlw(k)
	  d2=zlw(k+1)-zlc(k+1)

          cfIZW(k,1)=d2/(d1+d2)
	  cfIZW(k,2)=d1/(d1+d2)
      enddo    
! ----------------------------------------------------------------------------------
! Coefficients for interpolating from the   c(k) and c(k+1) levels on the w(k) level  
! ----------------------------------------------------------------------------------
      do k=0,NZ-1	 
          d1=zlw(k)-zlc(k)
	  d2=zlc(k+1)-zlw(k)

          cfIZC(k,1)=d2/(d1+d2)
	  cfIZC(k,2)=d1/(d1+d2)
      enddo


! ---------------------------------- 
! Other grid point locations
! ----------------------------------
      	do k=1,NZ
	 zwdif(k)=zlw(k)-zlw(k-1)
      	enddo

      	do k=0,NZ-1
	 zcdif(k)=zlc(k+1)-zlc(k)
      	enddo
     	zcdif(NZ)=-9999999.
  
	do k=1,NZ
	zdhp(k)=zlw(k)-zlc(k)
	zdhm(k)=zlc(k)-zlw(k-1)	 
     	enddo
! -------------------------------------
! FILE zgrid.dat OF GRID POINTS

      open(unit=40,file='zgrid.dat')
      do k=0,NZ
          write(40,'(6e16.7)') zlc(k),zwdif(k),zlw(k),zcdif(k),zdhm(k),zdhp(k)
      enddo
      close(40)

      end





!%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
! ------------------------------------------------------------------------
      subroutine inithordiffusion()
! ------------------------------------------------------------------------
 
! -----------------------------------------------------------------
! Defines horizontal  diffusion using masks for sponge layer 
! znuh(:,:) horizontal  diffusion  in h points
! znuoz     horizontal  diffusion  in o points
! -----------------------------------------------------------------

      implicit none
#include "commonlac.h"
      
      double hcount(RANGE2D)
      integer i,j,k
         
! -----------------------------------------------------------------
! Finds h-points having close land points (0 = it has land points close to it)

      hcount(:,:) = 0.d0 
      hcount(1:NX,1:NY)=llm(1:NX,1:NY,UGRID)*llm(0:(NX-1),1:NY,UGRID)
c2d     &                 *llm(1:NX,1:NY,VGRID)*llm(1:NX,0:(NY-1),VGRID)

      znuh(:,:)=znu*spongefac*(1-hcount(:,:))+znu*hcount(:,:)

! -----------------------------------------------------------------
! finds o-points close to land

#ifdef c2d
      hcount(:,:)=0.d0
      hcount(1:NX,1:NY)=llm(1:NX,1:NY,HGRID)*llm(2:(NX+1),1:NY,HGRID)
     &                  *llm(2:(NX+1),2:(NY+1),HGRID)*llm(1:NX,2:(NY+1),HGRID)

      znuoz(:,:)=znu*spongefac*(1-hcount(:,:))+znu*hcount(:,:)
#endif


! -----------------------------------------------------------------
/*
      do j=1,NY
	do i=1,NX
	  write(89,'(2e16.7,2i5)') znuh(i,j),znuoz(i,j),llm(i,j,HGRID),llm(i,j,OGRID)
	enddo
      enddo
      stop
*/    
      end

!%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
! ------------------------------------------------------------------------
      subroutine initvertdiffusion()
! ------------------------------------------------------------------------
! ------------------------------------------------------------------------
! Defines  vertical diffusion
!FOR PARABOLIC DIFFUSION TSANIS this coeff is K/H where K is the diffusion at point w in the vertical
!FOR CONSTANT  DIFFUSION  this coeff is K where K is the diffusion at point w in the vertical
!For LMIXING  it is not used
! ------------------------------------------------------------------------
      implicit none
#include "commonlac.h"      
      integer i,j,k

! -----------------------------------------------------------------
! ------------  VERTICAL DIFFUSION (dependency on H) ---------------
! -----------------------------------------------------------------

      do k=0,NZ
#ifdef KVPAR
!IF TSANIS this coeff is K/H where K is the diffusion at point w in the vertical
	  Kv(k)=LAMBDA*USTARS*(zlw(k)+ZBH)*(1-zlw(k)+ZSH)  
#else
! IF CONSTANT it is the true coefficient 
	  Kv(k)=KVERT
#endif
      enddo
! -----------------------------------------------------------------
!      print *,LAMBDA,USTARS,ZBH,ZSH

      if(mod(ntime,nsave).eq.0) then
      open(unit=40,file='vertdiff.dat')
      do k=0,NZ
	  write(40,*) zlw(k),Kv(k),KVERT
      enddo
      close(40)
      endif

      end



!%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
! --------------------------------------------------------
      subroutine readimmbound()
! --------------------------------------------------------

! --------------------------------------------------------
!     Reads immersed boundary from file
! --------------------------------------------------------

      implicit none
#include "commonlac.h"
      integer i,j,kk,kmin
      double s,r,rr,xx,yy,zz,ss,dmin

! --------------------------------------------------------
      print *,'Reading immersed boundary from file ',immboundfile
      open(unit=10,file=immboundfile,form='formatted')
      j=1

c2d      do j=1,NY
	 do i=1,NX
	    read(10,*) zimm(i,j)

! Position of immersed boundary
            zz=z(i,j)-zimm(i,j)
	    ss=-zz/z(i,j)+1
	    dmin=10
	    kmin=0
	    do kk=0,NZ
		 if(abs(ss-zlw(kk))<dmin) then
		    kmin=kk
		    dmin=abs(ss-zlw(kk))
		 endif
	    enddo
            izoff(i,j,HGRID)=kmin
	    print *,ss,kmin,zimm(i,j)
         enddo
c2d      enddo
      close(10)

! -----------------------

      end


!%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
! --------------------------------------------------------
      subroutine readbathymetry()
! --------------------------------------------------------

! --------------------------------------------------------
!     Reads bathymetry from file

!     Computes  mask for  h (or c or w points)  from bathymetry 
!     1=lake 0 = land 
! --------------------------------------------------------

      implicit none
#include "commonlac.h"
      integer i,j
      double s,r,rr,xx,yy

! --------------------------------------------------------

      open(unit=10,file=bathyfile,form='formatted')
      j=1

! If bathymetry is less than zmin these are land points for h, c or w   on the horizontal surface
c2d      do j=1,NY
	 do i=1,NX
	    read(10,*) z(i,j)
         enddo
c2d      enddo
	z(:,:)=z(:,:)*depth
	print *,depth
       depth=0.0

#ifdef MODBATHY
! Adds a small island at the center
      s=ISLANDBANKS ! width of side banks
      r=ISLANDRADIUS ! Island radius

c2d      do j=1,NY
  	 do i=1,NX
      
      xx=(i-NX/2.)*ddx
      yy=(j-2*NY/3.)*ddy
      rr=sqrt(xx*xx+yy*yy)
      z(i,j)=z(i,j)*(1-0.5*(tanh((rr+s*3+r)/s)*tanh((-rr+s*3+r)/s)+1))
      enddo
c2d      enddo
#endif

c2d      do j=1,NY
	 do i=1,NX  
	    if(z(i,j).lt.zmin) then
	      	z(i,j)=0
	      	llm(i,j,HGRID)=0
	    else 
 		llm(i,j,HGRID)=1
	    endif
	  
	    if(z(i,j).gt.depth) depth=z(i,j)

	 enddo
c2d      enddo

      close(10)

      print *,"Reading bathimetry. Maximum depth=",depth

! Make sure there are no openings on the boundaries

         llm(0,:,HGRID)=0
         llm(NX+1,:,HGRID)=0 
c2d      llm(:,0,HGRID)=0
c2d      llm(:,NY+1,HGRID)=0

! Not really needed but probably avoids divisions by 0

         z(0,:)=1;  z(NX+1,:)=1
c2d      z(:,0)=1; z(:,NY+1)=1


      open(unit=112,file='mybathy.dat')
c2d	do j=1,NY
	  do i=1,NX
	    write(112,*) z(i,j)
	  enddo
c2d	enddo
      close(112)
! -------------------------
      open(unit=92,file='mymask.dat')
      j=1
c2d      do j=1,NY
         do i=1,NX
	    write(92,*) llm(i,j,HGRID),llm(i,j,UGRID),llm(i,j,VGRID),llm(i,j,OGRID)
         enddo
c2d      enddo
      close(92)
! -----------------------

      end


!%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
! --------------------------------------------------------
      subroutine initbathymetry()
! --------------------------------------------------------
! --------------------------------------------------------
!     Initializes bathymetry for h (equivalently  c or  w points)
! --------------------------------------------------------

      implicit none
#include "commonlac.h"
      integer i,j
! --------------------------------------------------------
/*
! a 3D test
	do j=1,NX

	   if((j-0.5).lt.NX/4.) then
	      do i=1,NY
!                 z(i,j)=DEPTH*(1-0.5*(cos((j-0.5)*2*pi/NX.)+1)/2)	            			
		  z(j,i)=DEPTH*(1-0.66*((NX/4.-j+0.5)/NX))
	      enddo

	   elseif ((j-0.5).le.NX*3/4.) then

	      do i=1,NY
		z(j,i)=DEPTH
	      enddo

	   else
	      do i=1,NY
		  z(j,i)=DEPTH*(1-0.66*((j-0.5-NX*3./4.)/NX))
	     enddo
	   endif
	    
	  enddo
*/
! ---------------------------------
 
	do j=1,NY
	  do i=1,NX
	    z(i,j)=DEPTH
	  enddo
	enddo

      open(unit=112,file='mybathy.dat')
	do j=1,NY
	  do i=1,NX
	    write(112,*) z(i,j)
	  enddo
	enddo
      close(112)
! ---------------------------------

! Not really needed but probably avoids divisions by 0

         z(0,:)=DEPTH;  z(NX+1,:)=DEPTH
c2d      z(:,0)=DEPTH; z(:,NY+1)=DEPTH

      end


! %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
! -----------------------------------------------------
      subroutine readllmask()
! -----------------------------------------------------
! -----------------------------------------------------
!     Reads land-lake mask for h (equuivalently c or w  points)  from file llmfile; 1=lake 0 = land 
! -----------------------------------------------------

      implicit none
#include "commonlac.h"

      integer i,j

      open(unit=10,file=llmfile,form='formatted')
   
      j=1
c2d      do j=1,NY
	 do i=1,NX
	    read(10,*) llm(i,j,HGRID)
	 enddo
c2d      enddo

      close(10)

! ------------------------

! Make sure there are no openings on the boundaries

      llm(0,:,HGRID)=0
      llm(NX+1,:,HGRID)=0
c2d      llm(:,0,HGRID)=0
c2d      llm(:,NY+1,HGRID)=0

      end


!%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
! -----------------------------------------------------
      subroutine initmask(llm)
! -----------------------------------------------------


!     Makes simple 'rectangular swimming pool' land-lake mask 1=lake
      implicit none
      integer llm(0:NX+1,NYMIN:NYMAX,0:3)

      llm(:,:,HGRID)=1
      llm(0,:,HGRID)=0
      llm(NX+1,:,HGRID)=0
c2d      llm(:,0,HGRID)=0
c2d      llm(:,NY+1,HGRID)=0

 
      end


!%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%% 
! -----------------------------------------------------
       subroutine makemask3d(llm,llm3d,izoff,zlc)
! -----------------------------------------------------

      implicit none
      integer i,j,k
      integer izoff(0:(NX+1),NYMIN:NYMAX,0:3)
      integer llm(0:(NX+1),NYMIN:NYMAX,0:3)
      integer llmloc(0:(NX+1),NYMIN:NYMAX,0:3)
      integer llm3d(0:(NX+1),NYMIN:NYMAX,0:NZ,0:3)
      double zlc(0:NZ)

      llm3d(:,:,:,:)=0

      j=1
c2d      do j=NYMIN,NYMAX
        do i=0,(NX+1)
          do k=izoff(i,j,HGRID)+1,NZ
             llm3d(i,j,k,HGRID)=llm(i,j,HGRID)
	  enddo
	enddo
c2d      enddo 

      llm3d(0:NX,:,:,UGRID)=2-llm3d(0:NX,:,:,HGRID)-llm3d(1:(NX+1),:,:,HGRID)
!      llm3d(NX+1,:,:,UGRID)=2
c2d      llm3d(:,NYMIN:(NYMAX-1),:,VGRID)=2-llm3d(:,NYMIN:(NYMAX-1),:,HGRID)-llm3d(:,(NYMIN+1):NYMAX,:,HGRID)
!c2d      llm3d(:,NYMAX,:,VGRID)=2
   
      llm3d(:,:,:,UGRID)=llm3d(:,:,:,UGRID)/2
c2d      llm3d(:,:,:,VGRID)=llm3d(:,:,:,VGRID)/2

      izoff(:,:,UGRID)=SUM(llm3d(:,:,:,UGRID),DIM=3)-1
c2d        izoff(:,:,VGRID)=SUM(llm3d(:,:,:,VGRID),DIM=3)-1
      izoff(NX+1,:,UGRID)=0
c2d        izoff(:,NYMAX,VGRID)=0
 
      do k=0,NZ
           llmloc(:,:,HGRID)=llm3d(:,:,k,HGRID)        
           call makevecmask(llmloc) 
           llm3d(:,:,k,:)=llmloc(:,:,:)
      enddo

! XXX
      do k=0,NZ
c2d	    do j=NYMIN,NYMAX
          do i=0,(NX+1)
	     write(160,*) (i-0.5),zlc(k),llm3d(i,j,k,HGRID),llm3d(i,j,k,UGRID)
          enddo
c2d	enddo
      enddo

! XXX
      
c2d	    do j=NYMIN,NYMAX
          do i=0,(NX+1)
	     write(161,*) (i-0.5),izoff(i,j,HGRID),izoff(i,j,UGRID)
          enddo
c2d	enddo

      end

!%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%% 
! -----------------------------------------------------
       subroutine makevecmask(llm)
! -----------------------------------------------------
! -------------------------------------------
!  Extends land-lake mask to vector fields : u v o points 
! 'Rounded corners' for o
! -------------------------------------------

      implicit none
      integer i,j,k
      integer llm(0:(NX+1),NYMIN:NYMAX,0:3)



! if h(i,j) is land then    u(i,j), v(i,j)   are land
         llm(:,:,UGRID)    =llm(:,:,HGRID)
c2d      llm(:,:,VGRID)    =llm(:,:,HGRID)

! for o we count the adjoining land points first
!  o(i,j) is land if h(i,j) is land
c2d      llm(:,:,OGRID)    =1-llm(:,:,HGRID) ! HERE SPECIAL done to count the number of h points nearby


! u(i-1,j),o(i-1,j) and v(i,j-1),o(i,j-1) are land too if h(i,j) is land

      where (llm(1:(NX+1),:,HGRID)==0) 

              llm(0:NX,:,UGRID)=0
c2d           llm(0:NX,:,OGRID)=llm(0:NX,:,OGRID)+1

      end where

c2d      where (llm(:,1:(NY+1),HGRID)==0) 
c2d           llm(:,0:NY,VGRID)=0
c2d           llm(:,0:NY,OGRID)=llm(:,0:NY,OGRID)+1
c2d      end where

c2d      where (llm(1:(NX+1),1:(NY+1),HGRID)==0) 
c2d           llm(0:NX,0:NY,OGRID)=llm(0:NX,0:NY,OGRID)+1
c2d      end where

! Now set to 0 o-points with 2 or more adjoining land h-points

c2d      where(llm(:,:,OGRID).lt.2) llm(:,:,OGRID)=0
c2d      where(llm(:,:,OGRID).ge.2) llm(:,:,OGRID)=1
c2d      llm(:,:,OGRID)=1-llm(:,:,OGRID) 

! -------------------------
      open(unit=92,file='mymask.dat')

      j=1
c2d      do j=1,NY
         do i=1,NX
	    write(92,*) llm(i,j,HGRID),llm(i,j,UGRID),llm(i,j,VGRID),llm(i,j,OGRID)
         enddo
c2d      enddo

      close(92)
! -----------------------

      end


!%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

! -----------------------------------------------------
#ifdef TEMPERATURE
      subroutine initfields(u,v,h,c,t)
#else
      subroutine initfields(u,v,h,c)
#endif
! -----------------------------------------------------

! --------------------
! INITIAL CONDITIONS on (u,v,h,c,t) at t1 and t2
! --------------------

      implicit none
#include "commonlac.h"

      double h(RANGE2D,2)
      double u(RANGE,2)
      double v(RANGE,2)
      double c(RANGE,NSCAL,2)
#ifdef TEMPERATURE
      double t(RANGE,2)
#endif

      double dum,ran2,myz
      integer i,j,k,iseed,l

      iseed=-9371
      dum=ran2(iseed)	

      j=1 ! Necessary for the case  2D
       
! -------------- h VALUES for the waves ----------------
c2d      do j=1,NY
           do i=1,NX
               h(i,j,1)=(2*ran2(iseed)-1)*0.0000
               h(i,j,2)=h(i,j,1)
           enddo
c2d      enddo

! -----------------------------------------------------
	 do k=1,NZ
c2d      do j=0,NY+1
           do i=0,NX+1

            u(i,j,k,1)=0.0
            u(i,j,k,2)=0.0

c2d 	    v(i,j,k,1)=0.0
c2d 	    v(i,j,k,2)=0.0

            c(i,j,k,SCAL_NEUTRAL,1)=cinit(SCAL_NEUTRAL)*exp(-((i-NX/2)**2+(j-NY/2)**2+(k-NZ)**2)/(2.*50))*(z(i,j)+h(i,j,1))
#ifdef TEMPTANH  
		myz=(zlc(k)-1)*z(i,j)/z(NX/2,NY/2+1)
		c(i,j,k,SCAL_NEUTRAL,1)=(tanh((myz+0.333)*10)+10)*z(i,j)
#endif 
#ifdef DEBUGCONSTSCAL
				 
!                c(i,j,k,SCAL_NEUTRAL,1)=cinit(SCAL_NEUTRAL)
#endif
	    c(i,j,k,SCAL_SEDIMENT,1)=cinit(SCAL_SEDIMENT)*(z(i,j)+h(i,j,1)) ! 0
            c(i,j,k,SCAL_SEDIMENT,1)=cinit(SCAL_SEDIMENT)*exp(-((i-NX/2)**2+(j-NY/2)**2+(k-NZ)**2)/(2.*150))*(z(i,j)+h(i,j,1))
            c(i,j,k,SCAL_SEDIMENT,1)=cinit(SCAL_SEDIMENT)*(z(i,j)+h(i,j,1))
#ifdef ECOSYSTEM	   
	    c(i,j,k,SCAL_DETRITUS,1)=cinit(SCAL_DETRITUS)*(z(i,j)+h(i,j,1)) ! 9e-3
	   
! Oligotrophic lake initially [mg/m^3 and mgP/m^3]
	    c(i,j,k,SCAL_NUTRIENT,1)=cinit(SCAL_NUTRIENT)*(z(i,j)+h(i,j,1)) ! 8.

	    c(i,j,k,SCAL_PHYTO1,1)=cinit(SCAL_PHYTO1)*(z(i,j)+h(i,j,1))  ! 14

	    c(i,j,k,SCAL_ZOO,1)=cinit(SCAL_ZOO)*(z(i,j)+h(i,j,1))  ! 500
#ifdef PHYTO2
	    c(i,j,k,SCAL_PHYTO2,1)=cinit(SCAL_PHYTO2)*(z(i,j)+h(i,j,1))  ! 14
#endif
#endif

#ifdef TEMPERATURE
	    t(i,j,k,1)=tinit  
        t(i,j,k,2)=tinit  
#ifdef TEMPTANH  
		myz=(zlc(k)-1)*z(i,j)/z(NX/2,NY/2+1)
		t(i,j,k,1)=(tanh((myz+0.333)*10)+10)
		t(i,j,k,2)=t(i,j,k,1)
#endif 

#endif   
 	   enddo
c2d	 enddo
       enddo
! -----------------------------------------------------
       
       c(:,:,:,:,2)=c(:,:,:,:,1)

      end

!%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
! -----------------------------------------------------
#ifdef TEMPERATURE
      subroutine restart(u,v,h,c,t)
#else
      subroutine restart(u,v,h,c)
#endif
! -----------------------------------------------------
 
! --------------------
! INITIAL CONDITIONS WHEN RESTARTING FROM PREVIOUS FILE
! --------------------
 
      implicit none
#include "commonlac.h"

      double h(RANGE2D)
      double u(RANGE)
      double v(RANGE)
      double c(RANGE,NSCAL)
#ifdef TEMPERATURE
      double t(RANGE)
#endif
      open(unit=50,file='restart.unf',form='unformatted')
	 read(50) ntime
#ifdef TEMPERATURE
         read(50) h,u,v,c,t
#else
         read(50) h,u,v,c
#endif
      close(50)
      print *,'Restarted from step',ntime,'t=',ntime*ddt

      end




!%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%% 
! -------------------------------------------------------------------------------------------------------
#ifdef TEMPERATURE
      subroutine boundaries(u,v,h,c,t)
#else
      subroutine boundaries(u,v,h,c)
#endif
! ---------------------------------------------------------------------------- 
 
! -------------------------------------------------------------------------- 
! Fixes boundaries for all five variables u,v,h,c,t
! EXCEPT FREE SLIP ON THE SIDE USING FIX  directly in tendencies 
!  it takes care of all ghost points (on land and in the bottom layer)

!In subroutine Boundaries
!ON THE  SIDES walls: IMPERMEABILITY  at the land-lake boundaries U =0 or V=0 
!AT THE BOTTOM   NO SLIP  U =0 AND V=0   
!AT THE BOTTOM   No FLUX for T 
!at the bottom   BOUNDARY CONDITIONS   FOR THE FLUX OF SCALARS    


!IN subroutine TENDENCIES
!ON THE  SIDES walls: Free slip  for velocity
!ON THE  SIDES walls: Tracers Zero flux
!On TOP surface :  SHEAR AT THE SURFACE for velocity 
!On TOP surface :  Tracers Zero flux

! ---------------------------------------------------------------------------

      implicit none
#include "commonlac.h"
! ---------------------------------------------------------------------------
      double h(RANGE2D)
      double u(RANGE)
      double v(RANGE)
      double c(RANGE,NSCAL)
#ifdef TEMPERATURE
      double t(RANGE)
#endif

#ifndef IMPLICITDIFF
! ---------------------------------------------------------------------------
      double taub(RANGE2D)		! Bottom stress at w-points
      double taubu(RANGE2D)		! Bottom stress at oy-points (orthogonal plane to   y axis. center  of  h points) 
      double taubv(RANGE2D)		! Bottom stress at ox-points
      double ee(RANGE2D)		! Erosion
      double dd(RANGE2D)		! Deposition
      double hhi(RANGE2D)		! Deposition
#endif

      integer i,j,l,k

! ------------- IMPERMEABILITY ON THE SIDES at the land-lake boundaries  ( U =0 or V=0 according to the points u or v)---------------

#ifndef FIXED_VEL_PROF
      FIXBC3D(u,UGRID)
      FIXBC3D(v,VGRID)

#endif

! ------------- AT THE BOTTOM NON SLIP  U =0 AND V=0  -------------
! We use the fact that the vertical grid is symmetric respect to 0

#ifndef IMPLICITDIFF 
#ifdef IMMBOUND
        j=1
c2d      do j=NYMIN,NYMAX
        do i=0,(NX+1)
      u(i,j,izoff(i,j,UGRID))=-u(i,j,izoff(i,j,UGRID)+1)
c2d      v(i,j,izoff(i,j,VGRID))=-v(i,j,izoff(i,j,VGRID)+1)
	enddo
c2d      enddo
#else  
      u(:,:,0)=-u(:,:,1)  ! IF implicit this is not needed
      v(:,:,0)=-v(:,:,1)
#endif
#endif
! ---------------------------------------------------------------------------

!-------------  No FLUX for T at the bottom
! We use the fact that the vertical grid is symmetric respect to 0
#ifdef TEMPERATURE
#ifndef IMPLICITDIFF
      t(:,:,0)=t(:,:,1)     ! NOT NEEDED anymore since we use an implcit solver!
#endif
#endif
 
! -----------------------------------------
! These are just to check for bugs  
      where(llm(:,:,HGRID)==0) h(:,:)=-99999999.0
      do l=1,NSCAL
         do k=1,NZ 
            where(llm(:,:,HGRID)==0) c(:,:,k,l)=-999.0 
         enddo
      enddo

#ifndef IMPLICITDIFF

! ---------------------------------------------------------------------------
! ------------- Compute bottom stress taub at point c for bottom BCs ------------------

! compute  du/dsigma in oy (plane of w)

	 taubu(:,:)=DZC3(u,:,:,0) 
! -----------------------------------------

! compute  dv/dsigma in ox (plane of w)

c2d      taubv(:,:)=DZC3(v,:,:,0)  

! -----------------------------------------

! compute  the total shear amplitude in the the w-points 
      j=1  ! necessary for the 2D case 

#ifndef FIXED_VEL_PROF
c2d      do j=1,NY
	    do i=1,NX
	       taub(i,j)=IX2(taubu,i-1,j)**2
c2d     &  +IY2(taubv,i,j-1)**2
	    enddo
c2d      enddo

#else  
!        PROBLEM*** this is a test!
         taub(1,1)=taubu(1,1)**2 !*(1.56073254346215 **2)    	
#endif
 
! -----------------------------------------  
! compute  the stress in the the w-points 
        
#define INRANGE 1:NX,1:NY

#ifdef LMIXING
      taub(INRANGE)=RHOWAT*kvmin*sqrt(taub(INRANGE))/(h(INRANGE)+z(INRANGE)) !at bottom  Pacanovski gives kvmin
#else
#ifdef KVPAR 
      taub(INRANGE)=RHOWAT*Kv(0)*sqrt(taub(INRANGE))   ! PARABOLIC TSANIS 
#else
      taub(INRANGE)=RHOWAT*KVERT*sqrt(taub(INRANGE))/(h(INRANGE)+z(INRANGE))  ! CONSTANT DIFFUSION 
#endif

#endif

! ------------------------------
       if(mod(ntime,nsave).eq.0) then
       open(unit=54,file='taub.dat')
      do j=1,NY
	    do i=1,NX
	      write(54,*) ntime,ntime*ddt,taub(i,j)
	    enddo
      enddo 
!      close(54)
      endif
! -------------------------------

! %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%% SCALARS %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
 
        do l=1,NSCAL         
 
! ---------------------------  
!      J_sigma=  dd*c+ee 
!   ee should contain all stuff which is not multiplied by c
!   dd is everything which is multiplied by c
!-------------------------------------------------------------

! ---------   Erosion   ------------------  
	ee(:,:)=betac(l)*erom*max(taub(:,:)/tauce-1,0.d0) 
/*
!	if(l.eq.SCAL_SEDIMENT) print *,ee(1,1)

      open(54+l)
      do j=1,NY
      do i=1,NX
	   write(54+l,*) ee(i,j)
      enddo
      enddo 
      close(54+l)
*/

! -----------------  Deposition   --------------- 
! Bottom boundary conditions for scalars

	 select case (l)

!-------------------------------------------------------------
	   case (SCAL_SEDIMENT,SCAL_PHYTO1,SCAL_PHYTO2,SCAL_DETRITUS)
! -------- J_s=E+D   with  D= ws *c(0)*(1-|taub|/taucd-1) if taub>taucd here D=c(1/2)*dd

	      dd(:,:)=min(wsett(l),0.d0)*max(1-taub(:,:)/taucd(l),0.d0)
!-------------------------------------------------------------
	   case (SCAL_NUTRIENT)
! -------- J_s= cin +E 
! Just a flux now
	      dd(:,:)=0.d0
	    
	      ee(:,:)=cinnut+ee(:,:)
	    
! Old form J_s=h(cin-c)+E 
!	      dd(:,:)=-hnut
!	      ee(:,:)=ee(:,:)+hnut*cinnut
!-------------------------------------------------------------
	   case (SCAL_ZOO)
! -------- J_s=ee(:,:)

	      dd(:,:)=0.d0
!-------------------------------------------------------------
	   case (SCAL_NEUTRAL)
! -------- Js=0 a pure passive scalar

	      dd(:,:)=0.d0
	      ee(:,:)=0.d0

	    case default
   	      print *,'Unknown scalar type in BCs!'
	      stop

	 end select

! ----------    BOUNDARY CONDITIONS   FOR THE SCALARS   -----------------
!     The BC is written as J_sigma=wsett*c-K_v/H dc/dsigma = dd*c+ee
! -----------------------------------------------------------------------

! One defines hhi as  Kv(0)/(Hdsigma)

#ifdef LMIXING
	 hhi(:,:)=(kvmin+sdiff(l))/zcdif(0)/(h(:,:)+z(:,:)) 
#else
#ifdef KVPAR        
         hhi(:,:)=(Kv(0)+sdiff(l)/(h(:,:)+z(:,:)))/zcdif(0) 
#else
	 hhi(:,:)=(KVERT+sdiff(l))/zcdif(0)/(h(:,:)+z(:,:))
#endif
#endif

! -----------------------------------------------------------------------
 
!  (c-)*H =[(c+)*H*((-ws+D)/2+Kv/(Hdz)) + E H]/ ((ws-D)/2+Kv/(Hdz)) 
!In the code  the variable is  C*H (conservative form for tracers except temperature) 
! Checked three times!

         c(:,:,0,l)=(-c(:,:,1,l)*(    (wsett(l)-dd(:,:))*0.5 - hhi(:,:)  )+ee(:,:)*(h(:,:)+z(:,:)) )/
     &                          (    (wsett(l)-dd(:,:))*0.5 + hhi(:,:)  ) 
			
      enddo ! END DO LOOP SCALARS 
		
#endif ! not IMPLICITDIFF
      end

!%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%% 

#ifdef TIMESPLIT

      subroutine boundariesfast(um,vm,h)

! -------------------------------------------------------------------------- 
! Fixes boundaries for all five variables u,v,h
! EXCEPT FREE SLIP ON THE SIDE USING FIX  directly in tendencies 
!  it takes care of all ghost points (on land and in the bottom layer)

!In subroutine Boundaries
!ON THE  SIDES walls: IMPERMEABILITY  at the land-lake boundaries U =0 or V=0 

!IN subroutine TENDENCIES
!ON THE  SIDES walls: Free slip  for velocity

! ---------------------------------------------------------------------------

      implicit none
#include "commonlac.h"
! ---------------------------------------------------------------------------
      double h(RANGE2D)
      double um(RANGE2D)
      double vm(RANGE2D)

! ------------- IMPERMEABILITY ON THE SIDES at the land-lake boundaries  ( U =0 or V=0 according to the points u or v)---------------

#ifndef FIXED_VEL_PROF
      FIXBC2D(um,UGRID)
      FIXBC2D(vm,VGRID)
#endif
 
! -----------------------------------------
! These are just to check for bugs  
      where(llm(:,:,HGRID)==0) h(:,:)=-99999999.0

      end
#endif !TIMESPLIT

!%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%


      FUNCTION ran2(idum)

      INTEGER idum,IM1,IM2,IMM1,IA1,IA2,IQ1,IQ2,IR1,IR2,NTAB,NDIV
      double ran2,AM,EPS,RNMX
      PARAMETER (IM1=2147483563,IM2=2147483399,AM=1./IM1,IMM1=IM1-1,
     *IA1=40014,IA2=40692,IQ1=53668,IQ2=52774,IR1=12211,IR2=3791,
     *NTAB=32,NDIV=1+IMM1/NTAB,EPS=1.2e-7,RNMX=1.-EPS)
      INTEGER idum2,j,k,iv(NTAB),iy
      SAVE iv,iy,idum2
      DATA idum2/123456789/, iv/NTAB*0/, iy/0/
      if (idum.le.0) then
        idum=max(-idum,1)
        idum2=idum
        do 11 j=NTAB+8,1,-1
          k=idum/IQ1
          idum=IA1*(idum-k*IQ1)-k*IR1
          if (idum.lt.0) idum=idum+IM1
          if (j.le.NTAB) iv(j)=idum
11      continue
        iy=iv(1)
      endif
      k=idum/IQ1
      idum=IA1*(idum-k*IQ1)-k*IR1
      if (idum.lt.0) idum=idum+IM1
      k=idum2/IQ2
      idum2=IA2*(idum2-k*IQ2)-k*IR2
      if (idum2.lt.0) idum2=idum2+IM2
      j=1+iy/NDIV
      iy=iv(j)-idum2
      iv(j)=idum
      if(iy.lt.1)iy=iy+IMM1
      ran2=min(AM*iy,RNMX)
      return
      END
C  (C) Copr. 1986-92 Numerical Recipes Software 


!%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
! --------------------------------------------------
      subroutine readparameters
! --------------------------------------------------

! --------------------------------------------------
!     Reads parameters from external parameter file
! --------------------------------------------------
	
      implicit none
#include "commonlac.h"

      integer nskip,nr,l
      open(unit=10, file='params.dat', form='formatted')
 !	Skip comment lines
	nskip=0
	do while (nskip.lt.50)
	   read(10,*,err=5320) nr
	   exit
 5320        continue
           nskip=nskip+1
	enddo	
	if(nskip.ge.50) then 
	   write(*,*) 'No input data found in params.dat'
	   stop
	endif

! -------------------------------
!	Read parameters

	numtime=nr   ! Already read in nr
	read(10,*) nsave
	read(10,*) ncfl
	read(10,*) ddt
	read(10,*) domainx
	read(10,*) domainy
	read(10,*) depth
	read(10,*) windangle ! If Two dimensional must be 0
	read(10,*) uwindamp
	read(10,*) znu
        read(10,*) kvmin
        read(10,*) frestart
	read(10,*) freadbathy
	read(10,*) freadmask
	read(10,*) freadimmbound
	read(10,*) outfile
	read(10,*) bathyfile
	read(10,*) llmfile
	read(10,*) immboundfile
	read(10,*) zmin	
	read(10,*) tauce
	read(10,*) erom
	read(10,*) (betac(l),l=1,NSCAL)
	read(10,*) (wsett(l),l=1,NSCAL)
	read(10,*) (taucd(l),l=1,NSCAL)
	read(10,*) (sdiff(l),l=1,NSCAL) ! Minimum diffusivity of tracers  
	read(10,*) ctopnut
	read(10,*) cinnut
	read(10,*) spongefac
	read(10,*) (cinit(l),l=1,NSCAL)
	read(10,*) tinit
	read(10,*) tair
	read(10,*) hsurf
	read(10,*) phi_top
	close(10)

! Fix time units to seconds from days
	do l=1,NSCAL
	    wsett(l)=wsett(l)/86400
	enddo
	erom=erom/86400
        cinnut=cinnut/86400
	ctopnut=ctopnut/86400
!
	uwind=uwindamp

#ifndef c2d
	windangle=0
#endif

	end


! %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
! Ecosystem routines
! %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
! --------------------------------------------------
      subroutine readparameterseco
! --------------------------------------------------

! --------------------------------------------------
!     Reads ECOSYSTEM PARAMETERS from external parameter file
! --------------------------------------------------
 
      implicit none
#include "commoneco.h"

      integer nskip,l
      double nr

      open(unit=10, file='paramseco.dat', form='formatted')

! -------------------------------
!	Skip comment lines

	nskip=0

	do while (nskip.lt.50)
	   read(10,*,err=5200) nr
	   exit
 5200        continue
           nskip=nskip+1
	enddo
 	
	if(nskip.ge.50) then 
	   write(*,*) 'No input data found in paramseco.dat'
	   stop
	endif
      
! -------------------------------
!	Read parameters

	k_bg=nr   ! Already read in nr  
	read(10,*) k_S    
	read(10,*) k_D    
	read(10,*) mN1_D      
	read(10,*) k_P1   
	read(10,*) k_P2   
	read(10,*) mN1_P1   
	read(10,*) mN1_P2  
	read(10,*) l_P1   
	read(10,*) l_P2    
	read(10,*) mN1_Z  
	read(10,*) l_Z   
	read(10,*)  e_z  
	read(10,*) time_D  
	read(10,*)  mu_Z  
	read(10,*)  P_c  
	read(10,*)  mu_P1      
	read(10,*)  N_P1_c  
	read(10,*)  I_P1_c  
	read(10,*)  mu_P2  
	read(10,*)  N_P2_c  
	read(10,*)  I_P2_c 
	read(10,*)  alpha     
	read(10,*)  I_P2_c  
	read(10,*) gamma_c     
	read(10,*) Z_c      
	read(10,*) Fish        
       close(10)

! Fix time units to seconds from days
        mu_P1=mu_P1/86400
	mu_P2=mu_P2/86400
	l_P1=l_P1/86400
	l_P2=l_P2/86400
	mu_Z=mu_Z/86400
	L_z=L_z/86400
	gamma_c=gamma_c/86400
	time_D=time_D*86400

	end


!%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%% 
! ------------------------------------------------------------------------
      subroutine printparameters()
! ------------------------------------------------------------------------
      implicit none
#include "commonlac.h"
      double cflswx,cflswy

      print *
      print *, "---------------------------"
      print *, "Numerical model parameters:"
      print *, "---------------------------"
      print *
#ifdef IMPLICITDIFF 
      print *, "* Implicit vertical diffusion"
#else
      print *, "* Explicit vertical diffusion"
#endif

#ifdef TIMESPLIT
      print *, "* Split time stepping"
#else
      print *, "* Non-split time stepping"
#endif

#ifdef CONVADJ
      print *, "* Convective adjustment every ", CONVADJ, " steps"
#else
      print *, "* No convective adjustment"
#endif


      print *
      print *, "X resolution =",NX
      print *, "Y resolution =",NY
      print *, "Z resolution =",NZ
      print *, "Total num time steps",numtime     
      print *, "Save every",nsave," time steps"     
      print *, "Time \step : dt = ",ddt
      print *, "Horizontal  grid: dx = ",ddx, "dy = ",ddy
      cflswx=sqrt(GG*depth)*ddtfast/ddx
      cflswy=sqrt(GG*depth)*ddtfast/ddy

      print *, "CFL_sw=sqrt(gH)dtfast/dx = " , cflswx,cflswy,ddtfast

#ifndef c2d
      if((cflswx.ge.0.5).or.(cflswy.ge.0.5)) then
#else
      if((cflswx.ge.0.5)) then
#endif
	  print *,'Shallow water CFL violation. Please increase ddx or reduce (ddtfast, NFASTSTEPS or depth).'
	  stop
      endif

! ----------------------------------- 
#ifdef TANPAR
      print *, "Irregular vert. grid: TANPAR = ",TANPAR
#else
      print *, 'Regular vert. grid: dsigma = ',ddsigma
#endif
! ----------------------------------- 
      print *,'Domain size is ',domainx,'x',domainy,'m'
! ----------------------------------- 
      if(freadbathy) then
	print *,'Bathymetry read from ',bathyfile
      else
        print *,'Domain constant depth =',depth
      endif 
! ----------------------------------- 
      if(freadmask) then
	print *,'Land-lake mask read from ',llmfile
      endif 
! ----------------------------------- 

	print *,'Output written to ',outfile

      print *
      print *, "-----------------"
      print *, "Model parameters:"
      print *, "-----------------"
      print *
      print *,"Surface wind angle =",windangle,' degrees'
      print *,"Surface wind speed =",uwindamp,' m'
! ----------------------------------- 
      print *,"Horizontal diffusivity nu =",znu
! ----------------------------------- 
#ifndef LMIXING

#ifndef KVPAR
      print *,"Constant vertical diffusivity:"
      print *," Vert. diffusivity =",KVERT
#else
      print *,"Wu and Tsanis vertical diffusivity:"
      print *," lambda   =",LAMBDA
      print *," zsh   =",ZSH
      print *," zbh   =",ZBH
#endif

#else
      print *,"Minimum vertical diffusivity kvmin=",kvmin
#ifndef  PACANOWSKI
      print *,"Mixing Length diffusivity without density effects"
#else
      print *,"Mixing Length diffusivity with PACANOWSKI" 
#endif

#endif 
! ----------------------------------- 

      end


! %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
! --------------------------------------------------
      subroutine printparameterseco()
! --------------------------------------------------
 
      implicit none
#include "commoneco.h"

      print *, ""
      print *, "---------------------------"   
      print *, " Ecological Parameters"
      print *, "---------------------------"
     
        print *," k_bg   =", k_bg  
	print *," k_S",   k_S  
	print *," k_D ",  k_D    
	print *," mN1_D",mN1_D      
	print *," k_P1",k_P1   
	print *," k_P2 ",k_P2  
	print *," mN1_P1",mN1_P1   
	print *," mN1_P2 ",mN1_P2  
	print *," l_P1",l_P1   
	print *," l_P2",l_P2    
	print *," mN1_Z",mN1_Z  
	print *," l_Z ", l_Z   
	print *," e_z",e_z  
	print *," time_D ",time_D  
	print *," mu_Z ", mu_Z  
	print *," P_c", P_c  
	print *," mu_P1",mu_P1      
	print *," N_P1_c ",N_P1_c  
	print *," I_P1_c",I_P1_c  
	print *," mu_P2",mu_P2  
	print *," N_P2_c",N_P2_c  
	print *," I_P2_c",I_P2_c 
	print *," alpha ",alpha     
	print *," I_P2_c ",I_P2_c  
	print *," gamma_c  ",gamma_c     
	print *," Z_c",Z_c      
	print *," Fish ",Fish        

      end



! %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
#ifdef ECOSYSTEM
      subroutine ecotendpatrhs(c,ct,hhi)

! -------------------------------------------------------------------------------------------------------
! Ecosystem reaction contribution to tendencies
! Notice that since c(:,:,:,:) is actually c*H, there are some 1/H floating around in the nonlinear terms 
! -------------------------------------------------------------------------------------------------------

      implicit none
#include "commoneco.h"
#include "commonlac.h"
      double hhi(RANGE2D)

      double c(RANGE,NSCAL)
      double ct(RANGE,NSCAL) 

      integer k
      double  cfftime_D

      cfftime_D=1./time_D
     
#ifdef PHYTO2
#define MYFZ(k) mu_Z*hhi(:,:)*(PHY1(:,:,k)+PHY2(:,:,k))/(P_c +(PHY1(:,:,k)+PHY2(:,:,k))*hhi(:,:))  
#else
#define MYFZ(k) mu_Z*PHY1(:,:,k)*hhi(:,:)/(P_c +PHY1(:,:,k)*hhi(:,:))   
#endif

! XXXXXJvH CHECK ALSO THE STUFF ABOVE  in particular the use of hhi

! --------------------------------------------------------------------
!  REACTION PART:   D_t(c) = D_t(c) + Reaction   at point h or c 
! --------------------------------------------------------------------

      do k=1,NZ
              ct(:,:,k,SCAL_DETRITUS)= ct(:,:,k,SCAL_DETRITUS)   
     &        + mN1_P1*l_P1*PHY1(:,:,k)  
#ifdef PHYTO2
     &        + mN1_P2*l_P2*PHY2(:,:,k)  
#endif
     &        + mN1_Z*ZOO(:,:,k)*(l_Z +(1.-e_z)
     &        * MYFZ(k))

              ct(:,:,k,SCAL_NUTRIENT)= ct(:,:,k,SCAL_NUTRIENT)   
     &        +DET(:,:,k)*cfftime_D  
!#ifdef NUTCYCLE
!     &        +ctopnut*NUTCYCLE(ntime,ddt)
!#else
!	 & 	      +ctopnut
!#endif
  	  enddo

      end

#endif


! %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
#ifdef ECOSYSTEM
      subroutine ecotendpatlhs(c,ct,hhi)

! -------------------------------------------------------------------------------------------------------
! Ecosystem reaction contribution to tendencies
! Notice that since c(:,:,:,:) is actually c*H, there are some 1/H floating around in the nonlinear terms 
! -------------------------------------------------------------------------------------------------------

      implicit none
#include "commoneco.h"
#include "commonlac.h"
      double hhi(RANGE2D),hh(RANGE2D)

      double c(RANGE,NSCAL)
      double ct(RANGE,NSCAL) 

      double light(RANGE)
      double workFZ(RANGE) 
      double workFP1(RANGE) 
      double workFP2(RANGE) 
      double workPrefer(RANGE) 
      double workhhi(RANGE) 


      integer k,l,i,j
      double  cfftime_D

      cfftime_D=1./time_D

#if defined(RADIATION) || defined(ECOSYSTEM)
      call complight(c,light,workhhi,hhi,hh) ! lightz is not used, workhhi passed as dummy
#endif

! Needed to convert to actual concentrations
      do k=1,NZ
          workhhi(:,:,k)=hhi(:,:)
      enddo

! ---------------------------------------------------------------------
!  Compute some functions that are needed several times for reactions
! ---------------------------------------------------------------------

! WorkPrefer and workFZ are used simply as support/help arrays here
!   they will be defined below
	    workFZ(:,:,:) = light(:,:,:)/(I_P1_c+ light(:,:,:))	
	    workPrefer(:,:,:) = NUT(:,:,:)*workhhi(:,:,:)/( N_P1_c + NUT(:,:,:)*workhhi(:,:,:))
        workFP1(:,:,:)= mu_P1* min(workPrefer,workFZ)            

! Definition of workPrefer and workFZ     
#ifdef PHYTO2
	    workFZ(:,:,:) = light(:,:,:)/(I_P2_c+ light(:,:,:))	
	    workPrefer(:,:,:) = NUT(:,:,:)*workhhi(:,:,:)/( N_P2_c + NUT(:,:,:)*workhhi(:,:,:))
        workFP2(:,:,:)= mu_P2* min(workPrefer,workFZ)            

        workFZ(:,:,:)=mu_Z*workhhi(:,:,:)*(PHY1(:,:,:)+PHY2(:,:,:))/(P_c +(PHY1(:,:,:)+PHY2(:,:,:))*workhhi(:,:,:))  
#define PREFER1 (alpha)/(alpha*PHY1(:,:,:)+(1-alpha)*PHY2(:,:,:))            
#define PREFER2 (1-alpha)/(alpha*PHY1(:,:,:)+(1-alpha)*PHY2(:,:,:)) 
#else
        workFZ(:,:,:)=mu_Z*PHY1(:,:,:)*workhhi(:,:,:)/(P_c +PHY1(:,:,:)*workhhi(:,:,:))   
#define PREFER1 1.           
#define PREFER2 1.           
#endif

! --------------------------------------------------------------------
!  REACTION PART:   D_t(c) = D_t(c) + Reaction   at point h or c 
! --------------------------------------------------------------------

      do l=1,NSCAL

	select case (l)

! -----------------------------------
     	   case (SCAL_DETRITUS)
 
              ct(:,:,:,l)=     - cfftime_D 
 
! -------------------------------------
      	   case (SCAL_NUTRIENT)
 
              ct(:,:,:,l)=    
     &        - mN1_P1*workFP1(:,:,:)*PHY1(:,:,:)
#ifdef PHYTO2
     &        - mN1_P2*workFP2(:,:,:)*PHY2(:,:,:)
#endif
    
! -------------------------------------
 	   case (SCAL_PHYTO1)

              ct(:,:,:,l)= 
     &           (workFP1(:,:,:)-l_P1)
     &         -   workFZ(:,:,:)*PREFER1*ZOO(:,:,:)

#ifdef PHYTO2
! -------------------------------------
 	   case (SCAL_PHYTO2)

              ct(:,:,:,l)= 
     &           (workFP2(:,:,:)-l_P2)
     &         -   workFZ(:,:,:)*PREFER2*ZOO(:,:,:)
#endif
! -------------------------------------
 	   case (SCAL_ZOO)
 
              ct(:,:,:,l)= 
     &          (e_z*workFZ(:,:,:) -  l_Z)
     &         - gamma_c*ZOO(:,:,:)*workhhi(:,:,:)/(Z_c*Z_c+ZOO(:,:,:)*ZOO(:,:,:)*workhhi(:,:,:)**2)*Fish

! -------------------------------------
	   case (SCAL_NEUTRAL)
! -------------------------------------
	    
! -------------------------------------
	   case (SCAL_SEDIMENT)
! -------------------------------------

! -------------------------------------
	case default
! -------------------------------------

   	 print *,'Unknown scalar type in ecotend!'
	 stop

	 end select

      enddo ! tracer loop

      end

#endif



! %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
#ifdef ECOSYSTEM
      subroutine ecotend(c,ct,light,workFZ,workFP1,workFP2,workPrefer,workhhi,hhi)

! -------------------------------------------------------------------------------------------------------
! Ecosystem reaction contribution to tendencies
! Please provide also 6 work arrays
! Notice that since c(:,:,:,:) is actually c*H, there are some 1/H floating around in the nonlinear terms 
! -------------------------------------------------------------------------------------------------------

      implicit none
#include "commoneco.h"
#include "commonlac.h"
      double hhi(RANGE2D)

      double c(RANGE,NSCAL)
      double ct(RANGE,NSCAL) 

      double light(RANGE)
      double workFZ(RANGE) 
      double workFP1(RANGE) 
      double workFP2(RANGE) 
      double workPrefer(RANGE) 
      double workhhi(RANGE) 


      integer k,l,i,j
      double  cfftime_D

      cfftime_D=1./time_D

! Needed to convert to actual concentrations
      do k=1,NZ
          workhhi(:,:,k)=hhi(:,:)
      enddo

! ---------------------------------------------------------------------
!  Compute some functions that are needed several times for reactions
! ---------------------------------------------------------------------

! WorkPrefer and workFZ are used simply as support/help arrays here
!   they will be defined below
	workFZ(:,:,:) = light(:,:,:)/(I_P1_c+ light(:,:,:))	
	workPrefer(:,:,:) = NUT(:,:,:)*workhhi(:,:,:)/( N_P1_c + NUT(:,:,:)*workhhi(:,:,:))
        workFP1(:,:,:)= mu_P1* min(workPrefer,workFZ)            

!	do k=1,NZ
!	    write(38,'(4e16.7)') zlc(k),workFP1(1,1,k),workPrefer(1,1,k),workFZ(1,1,k)
!	enddo
!	stop
! Definition of workPrefer and workFZ     
#ifdef PHYTO2
	workFZ(:,:,:) = light(:,:,:)/(I_P2_c+ light(:,:,:))	
	workPrefer(:,:,:) = NUT(:,:,:)*workhhi(:,:,:)/( N_P2_c + NUT(:,:,:)*workhhi(:,:,:))
        workFP2(:,:,:)= mu_P2* min(workPrefer,workFZ)            

        workFZ(:,:,:)=mu_Z*workhhi(:,:,:)*(PHY1(:,:,:)+PHY2(:,:,:))/(P_c +(PHY1(:,:,:)+PHY2(:,:,:))*workhhi(:,:,:))  
	workPrefer(:,:,:)= (alpha*PHY1(:,:,:))/(alpha*PHY1(:,:,:)+(1-alpha)*PHY2(:,:,:))            
#else
        workFZ(:,:,:)=mu_Z*PHY1(:,:,:)*workhhi(:,:,:)/(P_c +PHY1(:,:,:)*workhhi(:,:,:))   
        workPrefer(:,:,:)= 1.           
#endif
! --------------------------------------------------------------------
!  REACTION PART:   D_t(c) = D_t(c) + Reaction   at point h or c 
! --------------------------------------------------------------------

      do l=1,NSCAL

	select case (l)

! -----------------------------------
     	   case (SCAL_DETRITUS)
 
              ct(:,:,:,l)= ct(:,:,:,l)   
     &        + mN1_P1*l_P1*PHY1(:,:,:)  
#ifdef PHYTO2
     &        + mN1_P2*l_P2*PHY2(:,:,:)  
#endif
     &        + mN1_Z*ZOO(:,:,:)*(l_Z +(1.-e_z)*workFZ(:,:,:))
     &        - DET(:,:,:)*cfftime_D  

! -------------------------------------
      	   case (SCAL_NUTRIENT)
 
              ct(:,:,:,l)= ct(:,:,:,l)   
     &        - mN1_P1*workFP1(:,:,:)*PHY1(:,:,:)
#ifdef PHYTO2
     &        - mN1_P2*workFP2(:,:,:)*PHY2(:,:,:)
#endif
     &        +DET(:,:,:)*cfftime_D  

! -------------------------------------
 	   case (SCAL_PHYTO1)

              ct(:,:,:,l)= ct(:,:,:,l)  
     &         +  (workFP1(:,:,:)-l_P1)*PHY1(:,:,:)
     &         -   workFZ(:,:,:)*workPrefer(:,:,:)*ZOO(:,:,:)

#ifdef PHYTO2
! -------------------------------------
 	   case (SCAL_PHYTO2)

              ct(:,:,:,l)= ct(:,:,:,l)  
     &         +  (workFP2(:,:,:)-l_P2)*PHY2(:,:,:)
     &         -   workFZ(:,:,:)*(1.-workPrefer(:,:,:))*ZOO(:,:,:)
#endif
! -------------------------------------
 	   case (SCAL_ZOO)
 
              ct(:,:,:,l)= ct(:,:,:,l)  
     &         + (e_z*workFZ(:,:,:) -  l_Z)*ZOO(:,:,:)
     &         - gamma_c*ZOO(:,:,:)*ZOO(:,:,:)*workhhi(:,:,:)/(Z_c*Z_c+ZOO(:,:,:)*ZOO(:,:,:)*workhhi(:,:,:)**2)*Fish

! -------------------------------------
	   case (SCAL_NEUTRAL)
! -------------------------------------
	    
! -------------------------------------
	   case (SCAL_SEDIMENT)
! -------------------------------------

! -------------------------------------
	case default
! -------------------------------------

   	 print *,'Unknown scalar type in ecotend!'
	 stop

	 end select

      enddo ! tracer loop

      end

#endif


! %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
! -------------------------------------------------------------------------------------------------------
      subroutine complight(c,light,lightz,hhi,hh)
! -------------------------------------------------------------------------------------------------------

! -------------------------------------------------------------------------------------------------------
! Computes light with a global amplitude unity.
! Amplitude is applied outside by radiation or Itop
! HERE C is C*H
! -------------------------------------------------------------------------------------------------------

      implicit none
#include "commoneco.h"
#include "commonlac.h"
      double hhi(RANGE2D),hh(RANGE2D)
      double c(RANGE,NSCAL)
      double light(RANGE),lightz(RANGE),phibulk
      integer k,l,i,j

#ifdef ECOSYSTEM
      double k_Det
      k_Det=k_D/mN1_D
#endif

! -------------------------------------------------------------------------------------------------------
! ONE COMPUTES Dln(I)/DZ and NOT Dln(I)/Dsigma

      do k=1,NZ
         lightz(:,:,k)=(
     &        k_bg + hhi(:,:)*(  
     &        k_S* SED(:,:,k) 
#ifdef ECOSYSTEM 
     &        +k_Det* DET(:,:,k) 
     &        +k_P1 * PHY1(:,:,k) 
#ifdef PHYTO2
     &        + k_P2 * PHY2(:,:,k)  
#endif  
#endif
     &        ))
      enddo
! ------------------------------------
!ONE COMPUTES  ln(I/Itop) with ln(I/Itop)=0 at the top surface

      light(:,:,NZ)= - ( zlw(NZ)- zlc(NZ) )*hh(:,:)*lightz(:,:,NZ)

! Regular, square integration

      do k=NZ-1,1,-1
          light(:,:,k)=light(:,:,k+1)-  zcdif(k)*hh(:,:)*lightz(:,:,k)
      enddo
! ------------------------------------
!ONE COMPUTES I/Itop

      light(:,:,:)= exp(light(:,:,:))

! ------------------------------------

!ONE COMPUTES D(I)/DZ 

       lightz(:,:,:)=lightz(:,:,:)*light(:,:,:)

! We introduce the amplitude of light at the top
#ifdef RADIATION
! light expressed in 'ecological units'  muMol ...
      call radiationbulk(phibulk) 

      light(:,:,:)=light(:,:,:)*phibulk/MUMOL2JOULE  ! 0.218 J = 1 muMol;  light is in mMol/m^2/sec; used for ecosystem 
#ifdef TEMPERATURE
      lightz(:,:,:)=lightz(:,:,:)*phibulk*(1./(RHOWAT*CWAT))  ! lightz is in Watt/m^2/rho/C* because it is used for temperature 
#endif   

#else

      phibulk=0.45*(1-ALBEDO)*phi_top
      light(:,:,:)=light(:,:,:)* phibulk/MUMOL2JOULE   ! lightz is not used 
#ifdef DIURNALCYCLE
!      light(:,:,:)=light(:,:,:)*(1+tanh(6*sin(2*pi*ntime*ddt/86400.d0-pi/2)))/2  
      light(:,:,:)=light(:,:,:)*LIGHTCYCLE(ntime,ddt)
#endif
#endif

! ------------------------------------
        if(mod(ntime,nsave).eq.0.0) then
         open(unit=37,file='light.dat')
	 do k=1,NZ
	  do j=1,NY
	    do i=1,NX
		write(37,*) zlc(k),light(i,j,k),lightz(i,j,k)
	    enddo
	  enddo
	 enddo
	 close(37)
	endif
! ------------------------------------

      end

! %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

! -------------------------------------------------------------------------------------------------------
      subroutine radiationbulk(phibulk)
! -------------------------------------------------------------------------------------------------------
! Radiation BCs on the surface 
! Uses Tair
! -------------------------------------------------------------------------------------------------------
      implicit none
#include "commonlac.h"
#include "commoneco.h"

      real*8 day,sp,sunpos
      real*8 phisol,phibulk

 ! --------------------------- 
       day=STARTDAY+ntime*ddt/86400.
 
! --------------------------- 
#ifdef DIURNALCYCLE
       day=STARTDAY+ntime*ddt/86400.
       sp=sunpos(LONGITUDE,LATITUDE,LTM,day)
#endif
! ---------------------------
! Short wave radiation from the sun
#ifdef DIURNALCYCLE
       phisol=AASKY*PHI0*sin(sp)**BBSKY*(1-0.65*CLOUDCOV**2) 
#else
       phisol=phi_top
#endif

! Short wave radiation from the sun penetrating in the bulk (light)
       phibulk=0.45*(1-ALBEDO)*phisol

      end

! -------------------------------------------------------------------------------------------------------

#ifdef TEMPERATURE
      subroutine radiationsurf(temp,phisurf)
! -------------------------------------------------------------------------------------------------------
! Radiation BCs on the surface 
! Uses Tair
! -------------------------------------------------------------------------------------------------------
      implicit none
#include "commonlac.h"
#include "commoneco.h"

      real*8 day,sp,sunpos
      real*8 phisol
      real*8 temp(RANGE2D),humair
      real*8 phiwat(RANGE2D)
      real*8 phiatm,phiconvec(RANGE2D),philatent(RANGE2D),phisurf(RANGE2D)
      integer k,l,i,j
 ! --------------------------- 
       day=STARTDAY+ntime*ddt/86400.
 
! --------------------------- 
#ifdef DIURNALCYCLE
       day=STARTDAY+ntime*ddt/86400.
       sp=sunpos(LONGITUDE,LATITUDE,LTM,day)
#endif
! ---------------------------
! Short wave radiation from the sun
#ifdef DIURNALCYCLE
       phisol=AASKY*PHI0*sin(sp)**BBSKY*(1-0.65*CLOUDCOV**2) 
#else
       phisol=phi_top
#endif

! ---------------------------       
! Short wave radiation from the sun stopped in the first centimeters
       phisol=0.55*(1-ALBEDO)*phisol
! ---------------------------       
! Long wave radiation flux from water (:,:)
       phiwat=EPSWAT*STEFAN*(temp+273.15)**4
! ---------------------------
! Long wave radiation flux from sky
       phiatm=0.937e-5* STEFAN*(TAIR+273.15)**6*(1+0.15*CLOUDCOV**2)
! ---------------------------
! Sensible heat flux (:,:)
       phiconvec=RHOAIR*CAIR*(AACONVEC+BBCONVEC*UWIND)*(temp-TAIR)
! ---------------------------
! Latent heat flux
       humair=exp(2.3026*(7.5*TAIR)/(TAIR+273.15)+0.7858)
       philatent=0.622/PAIR*CEVAP*RHOAIR*UWIND*2475e3*(exp(2.3026*(7.5*temp)/(temp+273.15)+0.7858)-humair)
! ---------------------------
       phisurf=phisol+phiatm-phiconvec-philatent-phiwat

      end
#endif

! %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
 
      function sunpos(lon,lat,ltm,day)

! based on a routine by N.E. Zimmermann (niklaus.zimmermann@wsl.ch)

      real*8      lon,lat
      integer     ltm         ! local time meridian (in hours, NOT degrees. Add 1 for summer time)

      REAL*8      lsn         ! local solar noon (standard time)
      REAL*8      day        ! Julian Day
      INTEGER     hr          ! hour (local standard time)
      INTEGER     min         ! minute (local standard time)
      INTEGER     hangle      ! hourangle (local standard time)   
      integer     jday

      REAL*8        pi          ! PI
      REAL*8        d2r         ! degree to radians conversion
      REAL*8        r2d         ! radians to degree conversion
      REAL*8        decrad      ! declination of the sun
      REAL*8        latrad      ! latitude of site in radians
      REAL*8        sunpos     ! solar altitude in radians
! 
! local variables
!
      REAL*8        ha          ! hour (floating from hr & min)
      REAL*8        harad       ! hourangle in radians

! -------------------------------------------------------------------- *

      pi      = 3.14159265358979
      d2r     = pi/180.0
      r2d     = 1/d2r
      lsn     = 12.0+(float(ltm)-lon/15.0)
      latrad  = lat*d2r
      jday = int(day)
      ha=(day-jday)*24.
      decrad  = 23.45*d2r*sin(d2r*360.*(284.+FLOAT(jday))/365.)
      harad  = (lsn-ha)/24.*2*pi           !solrad is given for the hour preceeding the time given
                       
!      print *,latrad,decrad,harad               
      sunpos = asin((sin(latrad)*sin(decrad))+(cos(latrad)*cos(decrad)*cos(harad)))

      if (sunpos.LT.0.0 .OR. sunpos.GT.pi) sunpos = 0.0 ! sun is below horizon
      
      END


! %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
! -------------------------------------------------------------------------------------------------------
      subroutine convadj(t,c)
! -------------------------------------------------------------------------------------------------------

! -------------------------------------------------------------------------------------------------------
! Trivial convective adjustment works only when TEMPERATURE IS ACTIVE
! -------------------------------------------------------------------------------------------------------

      implicit none
#include "commonlac.h"

      real*8 t(RANGE),c(RANGE,NSCAL)
      integer flips(RANGE2D)
      
      integer k,i,j  !,imask(RANGE2D)

      print *,'*** CONVECTIVE ADJUSTMENT ***'

! Where cannot contain loops, makes the following clumsy structure necessary 
      flips(:,:)=0 
      do k=NZ,2,-1
	 where(RHOEOS(t(:,:,k)).gt.RHOEOS(t(:,:,k-1))) 

            t(:,:,k)=(t(:,:,k-1)*zwdif(k-1)+t(:,:,k)*zwdif(k))/(zwdif(k-1)+zwdif(k))
            c(:,:,k,1)=(c(:,:,k-1,1)*zwdif(k-1)+c(:,:,k,1)*zwdif(k))/(zwdif(k-1)+zwdif(k))

            t(:,:,k-1)=t(:,:,k)
            c(:,:,k-1,1)=c(:,:,k,1)
#if NSCAL>1
	    c(:,:,k,2)=(c(:,:,k-1,2)*zwdif(k-1)+c(:,:,k,2)*zwdif(k))/(zwdif(k-1)+zwdif(k))
            c(:,:,k-1,2)=c(:,:,k,2)
#endif
#if NSCAL>2
	    c(:,:,k,3)=(c(:,:,k-1,3)*zwdif(k-1)+c(:,:,k,3)*zwdif(k))/(zwdif(k-1)+zwdif(k))
            c(:,:,k-1,3)=c(:,:,k,3)
#endif
#if NSCAL>3
	    c(:,:,k,4)=(c(:,:,k-1,4)*zwdif(k-1)+c(:,:,k,4)*zwdif(k))/(zwdif(k-1)+zwdif(k))
            c(:,:,k-1,4)=c(:,:,k,4)
#endif
#if NSCAL>4
	    c(:,:,k,5)=(c(:,:,k-1,5)*zwdif(k-1)+c(:,:,k,5)*zwdif(k))/(zwdif(k-1)+zwdif(k))
            c(:,:,k-1,5)=c(:,:,k,5)
#endif
#if NSCAL>5
	    c(:,:,k,6)=(c(:,:,k-1,6)*zwdif(k-1)+c(:,:,k,6)*zwdif(k))/(zwdif(k-1)+zwdif(k))
            c(:,:,k-1,6)=c(:,:,k,6)
#endif
#if NSCAL>6
	    c(:,:,k,7)=(c(:,:,k-1,7)*zwdif(k-1)+c(:,:,k,7)*zwdif(k))/(zwdif(k-1)+zwdif(k))
            c(:,:,k-1,7)=c(:,:,k,7)
#endif
            flips(:,:)=flips(:,:)+1
	 end where
      enddo
      print *,'Flips: ',sum(flips)
      end

      subroutine solve_tridiag(a,b,c,d,r)
      implicit none
!        a - sub-diagonal (means it is the diagonal below the main diagonal)
!        b - the main diagonal
!        c - sup-diagonal (means it is the diagonal above the main diagonal)
!        d - first right part
! 	     r - INPUT : additional rhs part (RHS=d+r)  OUTPUT: solution
!        n - number of equations -->> NZ

        double a(RANGE)
 		double b(RANGE)
		double c(RANGE)
		double d(RANGE)
		double r(RANGE)
        double m(RANGE2D)
        integer i
 
! initialize c-prime and d-prime
        c(:,:,1) = c(:,:,1)/b(:,:,1)
        d(:,:,1) = (d(:,:,1)+r(:,:,1))/b(:,:,1)
! solve for vectors c-prime and d-prime
         do i = 2,NZ
           m(:,:) = b(:,:,i)-c(:,:,i-1)*a(:,:,i)
           c(:,:,i) = c(:,:,i)/m(:,:)
           d(:,:,i) = ((d(:,:,i)+r(:,:,i))-d(:,:,i-1)*a(:,:,i))/m(:,:)
         enddo
! initialize x
         r(:,:,NZ) = d(:,:,NZ)
! solve for x from the vectors c-prime and d-prime
        do i = NZ-1, 1, -1
          r(:,:,i) = d(:,:,i)-c(:,:,i)*r(:,:,i+1)
        end do
 
       end subroutine solve_tridiag

