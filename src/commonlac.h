! Common file for lac.F 

      double z(RANGE2D)
      double zimm(RANGE2D)
      
      double pi
      double  ddt,ddx,ddy,ddsigma,ddtfast
      integer ntime,numtime,nsave,ncfl,ntimefast

      double znu,Kv(0:NZ),kvmin
      double znuoz(RANGE2D),znuh(RANGE2D)

      double cffx,cffy,cffsigma 
      double zlc(0:NZ),zlw(0:NZ)
      double zwdif(0:NZ),zcdif(0:NZ)
      double zdhm(0:NZ),zdhp(0:NZ)
      double cfIZW(0:NZ,0:3),cfDZW(0:NZ,0:3)
      double cfIZC(0:NZ,0:3),cfDZC(0:NZ,0:3)

#ifdef LMIXING
      double kvhmix(RANGE)
	  double kvhmixu(RANGE)
c2d	  double kvhmixv(RANGE)
#endif  
    
      integer llm(0:NX+1,NYMIN:NYMAX,0:3),izoff(0:NX+1,NYMIN:NYMAX,0:3)
      integer llm3d(0:(NX+1),NYMIN:NYMAX,0:NZ,0:3)
	    
      double domainx,domainy,depth,windangle,uwind,zmin,uwindamp
      double lambda,zsh,zbh,tauce,erom,cinnut,spongefac,ctopnut
      double betac(NSCAL),wsett(NSCAL),taucd(NSCAL),cinit(NSCAL),sdiff(NSCAL)
      double hsurf,tinit,tair,phi_top

      logical freadbathy,freadmask,frestart,feuler,freadimmbound,feulerfast
      character*30,outfile,bathyfile,llmfile,immboundfile
      parameter(pi=3.14159265358979)
      

      common /WUTS/ lambda,zsh,zbh 
      common /MODEL/  domainx,domainy,depth,windangle,uwind,uwindamp
      common /NSYS/ ntime,ntimefast,numtime,nsave,ncfl,feuler,feulerfast
      common /XSYS/ ddt,ddx,ddy,ddsigma,ddtfast
      common /PHYS/ znu,Kv,spongefac,znuh,znuoz,kvmin
      common /FILENAMES/ outfile,bathyfile,llmfile,immboundfile
      common /FREESURF/ z,zimm
      common /CFF/ cffx, cffy, cffsigma
      common /VERTGRID/ zlc,zlw,cfIZC,cfDZC,cfIZW,cfDZW,zwdif,zcdif,zdhm,zdhp
      common /MASK/  llm,llm3d,freadmask,freadbathy,zmin,frestart,freadimmbound,izoff
      common /SEDIMENT/ tauce,taucd,erom,betac,wsett,cinnut,cinit,sdiff,ctopnut
      common /AIR/ hsurf,tinit,tair,phi_top

#ifdef LMIXING
      common /MIX/ kvhmix,kvhmixu
c2d     &             ,kvhmixv
#endif
