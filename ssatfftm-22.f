      PROGRAM NCORPOS
      IMPLICIT REAL *8(a-h,o-z)
      DIMENSION X(36),V(36)	
      dimension a(6),e(6),di(6),w(6),om(6),am(6),dm(6),ami(6),dmimi(6),
     |xpla(6,6),f(6),bx(6,1),by(6,1),bz(6,1)
      REAL*8 XAFFT(1048576),XEFFT(1048576),
     |XIFFT(1048576),XE2FFT(1048576)
c	XA1FFT(1048576),XI1FFT(1048576),
c     >,XA3FFT(1048576),XI3FFT(1048576),
c     >XA4FFT(1048576),XA5FFT(1048576),XE5FFT(1048576),XA6FFT(1048576)
      REAL*8 NLINHAS,DELTAT,ttf,tint
	INTEGER nespec01,nespec05,nespec1,nespec3,
     >nespec5,nespec10,nespec15,dn01,dn05,dn1,dn3,dn5,dn10,dn15
      COMMON/MASSAS/dm0,G,dj2,dj4
      COMMON/IN/ikepler,inicio,pi,pi05,pi15,pi2
      COMMON/CONST/conv
      COMMON/NUMERO/N,neq      
      common/tempo/tint,ttf,Nlinhas
	common/novo/aa,ee
	ioi=79

      OPEN(123,FILE='emax.dat',STATUS='UNKNOWN')
      OPEN(1234,FILE='amax.dat',STATUS='UNKNOWN')
      OPEN(12345,FILE='imax.dat',STATUS='UNKNOWN')

      OPEN(1239,FILE='emax.grd',STATUS='UNKNOWN')
      OPEN(12349,FILE='amax.grd',STATUS='UNKNOWN')
      OPEN(23459,FILE='imax.grd',STATUS='UNKNOWN')

      OPEN(923,FILE='emax-t.dat',STATUS='UNKNOWN')
      OPEN(9234,FILE='amax-t.dat',STATUS='UNKNOWN')
      OPEN(92345,FILE='imax-t.dat',STATUS='UNKNOWN')

	write(1239,'(4hDSAA)')
	write(1239,*)'51 51'
	write(1239,*)'198010.514 198209.3984'
	write(1239,*)'0 +0.01'
   	write(1239,*)'1 100'

	write(12349,'(4hDSAA)')
	write(12349,*)'51 51'
	write(12349,*)'198010.514 198209.3984'
	write(12349,*)'0 +0.01'
   	write(12349,*)'1 100'

	write(23459,'(4hDSAA)')
	write(23459,*)'51 51'
	write(23459,*)'198010.514 198209.3984'
	write(23459,*)'0 +0.01'
   	write(23459,*)'1 100'

c199	Format(F14.8,1x,F14.8,1x,I5)

c      open (522,file='t.DAT', status='unknown') 	

c	Numero de linhas=potencia de 2	
c	2**21
   	Nlinhas=1048576.0d0
 
 	ttf=0.9d0/5.0d0
 	tf=ttf

c    ~516.75  anos
     					  
	tint=Nlinhas*ttf-ttf
c	write(*,*)tint/365.25d0
	

c	ARQUIVOS DE FFT DO SEMI-EIXO		
      OPEN(2,FILE='dn01a.grd',STATUS='UNKNOWN')
      OPEN(22,FILE='e1e2dn01a.dat',STATUS='UNKNOWN')

      OPEN(3,FILE='dn05a.grd',STATUS='UNKNOWN')
      OPEN(33,FILE='e1e2dn05a.dat',STATUS='UNKNOWN')

      OPEN(4,FILE='dn1a.grd',STATUS='UNKNOWN')
      OPEN(44,FILE='e1e2dn1a.dat',STATUS='UNKNOWN')

      OPEN(41,FILE='dn3a.grd',STATUS='UNKNOWN')
      OPEN(441,FILE='e1e2dn3a.dat',STATUS='UNKNOWN')

      OPEN(5,FILE='dn5a.grd',STATUS='UNKNOWN')
      OPEN(55,FILE='e1e2dn5a.dat',STATUS='UNKNOWN')

      OPEN(622,FILE='dn10a.grd',STATUS='UNKNOWN')
      OPEN(662,FILE='e1e2dn10a.dat',STATUS='UNKNOWN')

      OPEN(7,FILE='dn15a.grd',STATUS='UNKNOWN')
      OPEN(77,FILE='e1e2dn15a.dat',STATUS='UNKNOWN')


	write(2,'(4hDSAA)')
	write(2,*)'51 71'
	write(2,*)'198010.514 198209.3984'
	write(2,*)'0 +0.01'
   	write(2,*)'1 100'

	write(3,'(4hDSAA)')
	write(3,*)'51 71'
	write(3,*)'198010.514 198209.3984'
	write(3,*)'0 +0.01'
   	write(3,*)'1 100'

	write(4,'(4hDSAA)')
	write(4,*)'51 71'
	write(4,*)'198010.514 198209.3984'
	write(4,*)'0 +0.01'
   	write(4,*)'1 100'

	write(41,'(4hDSAA)')
	write(41,*)'51 71'
	write(41,*)'198010.514 198209.3984'
	write(41,*)'0 +0.01'
   	write(41,*)'1 100'

	write(5,'(4hDSAA)')
	write(5,*)'51 71'
	write(5,*)'198010.514 198209.3984'
	write(5,*)'0 +0.01'
   	write(5,*)'1 100'

	write(622,'(4hDSAA)')
 	write(622,*)'51 71'
	write(622,*)'198010.514 198209.3984'
	write(622,*)'0 +0.01'
   	write(622,*)'1 100'

	write(7,'(4hDSAA)')
	write(7,*)'51 71'
	write(7,*)'198010.514 198209.3984'
	write(7,*)'0 +0.01'
   	write(7,*)'1 100'

C	ARQUIVOS DE FFT da excentricidade

      OPEN(29,FILE='dn01.grd',STATUS='UNKNOWN')
      OPEN(229,FILE='e1e2dn01.dat',STATUS='UNKNOWN')

      OPEN(39,FILE='dn05.grd',STATUS='UNKNOWN')
      OPEN(339,FILE='e1e2dn05.dat',STATUS='UNKNOWN')

      OPEN(49,FILE='dn1.grd',STATUS='UNKNOWN')
      OPEN(449,FILE='e1e2dn1.dat',STATUS='UNKNOWN')

      OPEN(491,FILE='dn3.grd',STATUS='UNKNOWN')
      OPEN(4491,FILE='e1e2dn3.dat',STATUS='UNKNOWN')

      OPEN(59,FILE='dn5.grd',STATUS='UNKNOWN')
      OPEN(559,FILE='e1e2dn5.dat',STATUS='UNKNOWN')

      OPEN(6229,FILE='dn10.grd',STATUS='UNKNOWN')
      OPEN(6629,FILE='e1e2dn10.dat',STATUS='UNKNOWN')

      OPEN(79,FILE='dn15.grd',STATUS='UNKNOWN')
      OPEN(779,FILE='e1e2dn15.dat',STATUS='UNKNOWN')


	write(29,'(4hDSAA)')
	write(29,*)'51 71'
	write(29,*)'198010.514 198209.3984'
	write(29,*)'0 +0.01'
   	write(29,*)'1 100'

	write(39,'(4hDSAA)')
	write(39,*)'51 71'
	write(39,*)'198010.514 198209.3984'
	write(39,*)'0 +0.01'
   	write(39,*)'1 100'

	write(49,'(4hDSAA)')
	write(49,*)'51 71'
	write(49,*)'198010.514 198209.3984'
	write(49,*)'0 +0.01'
   	write(49,*)'1 100'

	write(491,'(4hDSAA)')
	write(491,*)'51 71'
	write(491,*)'198010.514 198209.3984'
	write(491,*)'0 +0.01'
   	write(491,*)'1 100'

	write(59,'(4hDSAA)')
	write(59,*)'51 71'
	write(59,*)'198010.514 198209.3984'
	write(59,*)'0 +0.01'
   	write(59,*)'1 100'

	write(6229,'(4hDSAA)')
 	write(6229,*)'51 71'
	write(6229,*)'198010.514 198209.3984'
	write(6229,*)'0 +0.01'
   	write(6229,*)'1 100'

	write(79,'(4hDSAA)')
	write(79,*)'51 71'
	write(79,*)'198010.514 198209.3984'
	write(79,*)'0 +0.01'
   	write(79,*)'1 100'

C	ARQUIVOS DE FFT da inclinacao

      OPEN(290,FILE='dn01i.grd',STATUS='UNKNOWN')
      OPEN(2290,FILE='e1e2dn01i.dat',STATUS='UNKNOWN')

      OPEN(390,FILE='dn05i.grd',STATUS='UNKNOWN')
      OPEN(3390,FILE='e1e2dn05i.dat',STATUS='UNKNOWN')

      OPEN(490,FILE='dn1i.grd',STATUS='UNKNOWN')
      OPEN(4490,FILE='e1e2dn1i.dat',STATUS='UNKNOWN')

      OPEN(4901,FILE='dn3i.grd',STATUS='UNKNOWN')
      OPEN(44901,FILE='e1e2dn3i.dat',STATUS='UNKNOWN')

      OPEN(590,FILE='dn5i.grd',STATUS='UNKNOWN')
      OPEN(5590,FILE='e1e2dn5i.dat',STATUS='UNKNOWN')

      OPEN(62290,FILE='dn10i.grd',STATUS='UNKNOWN')
      OPEN(66290,FILE='e1e2dn10i.dat',STATUS='UNKNOWN')

      OPEN(790,FILE='dn15i.grd',STATUS='UNKNOWN')
      OPEN(7790,FILE='e1e2dn15i.dat',STATUS='UNKNOWN')

	write(290,'(4hDSAA)')
	write(290,*)'51 71'
	write(290,*)'198010.514 198209.3984'
	write(290,*)'0 +0.01'
   	write(290,*)'1 100'

	write(390,'(4hDSAA)')
	write(390,*)'51 71'
	write(390,*)'198010.514 198209.3984'
	write(390,*)'0 +0.01'
   	write(390,*)'1 100'

	write(490,'(4hDSAA)')
	write(490,*)'51 71'
	write(490,*)'198010.514 198209.3984'
	write(490,*)'0 +0.01'
   	write(490,*)'1 100'

	write(4901,'(4hDSAA)')
	write(4901,*)'51 71'
	write(4901,*)'198010.514 198209.3984'
	write(4901,*)'0 +0.01'
   	write(4901,*)'1 100'

	write(590,'(4hDSAA)')
	write(590,*)'51 71'
	write(590,*)'198010.514 198209.3984'
	write(590,*)'0 +0.01'
   	write(590,*)'1 100'

	write(62290,'(4hDSAA)')
 	write(62290,*)'51 71'
	write(62290,*)'198010.514 198209.3984'
	write(62290,*)'0 +0.01'
   	write(62290,*)'1 100'

	write(790,'(4hDSAA)')
	write(790,*)'51 71'
	write(790,*)'198010.514 198209.3984'
	write(790,*)'0 +0.01'
   	write(790,*)'1 100'


C	ARQUIVOS DE FFT ecosvarpi

      OPEN(4290,FILE='dn01evarpi.grd',STATUS='UNKNOWN')
      OPEN(42290,FILE='e1e2dn01evarpi.dat',STATUS='UNKNOWN')

      OPEN(4390,FILE='dn05evarpi.grd',STATUS='UNKNOWN')
      OPEN(43390,FILE='e1e2dn05evarpi.dat',STATUS='UNKNOWN')

      OPEN(4490,FILE='dn1evarpi.grd',STATUS='UNKNOWN')
      OPEN(44490,FILE='e1e2dn1evarpi.dat',STATUS='UNKNOWN')

      OPEN(44901,FILE='dn3evarpi.grd',STATUS='UNKNOWN')
      OPEN(444901,FILE='e1e2dn3evarpi.dat',STATUS='UNKNOWN')

      OPEN(4590,FILE='dn5evarpi.grd',STATUS='UNKNOWN')
      OPEN(45590,FILE='e1e2dn5evarpi.dat',STATUS='UNKNOWN')

      OPEN(462290,FILE='dn10evarpi.grd',STATUS='UNKNOWN')
      OPEN(466290,FILE='e1e2dn10i.dat',STATUS='UNKNOWN')

      OPEN(4790,FILE='dn15evarpi.grd',STATUS='UNKNOWN')
      OPEN(47790,FILE='e1e2dn15evarpi.dat',STATUS='UNKNOWN')

	write(4290,'(4hDSAA)')
	write(4290,*)'51 71'
	write(4290,*)'198010.514 198209.3984'
	write(4290,*)'0 +0.01'
   	write(4290,*)'1 100'

	write(4390,'(4hDSAA)')
	write(4390,*)'51 71'
	write(4390,*)'198010.514 198209.3984'
	write(4390,*)'0 +0.01'
   	write(4390,*)'1 100'

	write(4490,'(4hDSAA)')
	write(4490,*)'51 71'
	write(4490,*)'198010.514 198209.3984'
	write(4490,*)'0 +0.01'
   	write(4490,*)'1 100'

	write(44901,'(4hDSAA)')
	write(44901,*)'51 71'
	write(44901,*)'198010.514 198209.3984'
	write(44901,*)'0 +0.01'
   	write(44901,*)'1 100'

	write(4590,'(4hDSAA)')
	write(4590,*)'51 71'
	write(4590,*)'198010.514 198209.3984'
	write(4590,*)'0 +0.01'
   	write(4590,*)'1 100'

	write(462290,'(4hDSAA)')
 	write(462290,*)'51 71'
	write(462290,*)'198010.514 198209.3984'
	write(462290,*)'0 +0.01'
   	write(462290,*)'1 100'

	write(4790,'(4hDSAA)')
	write(4790,*)'51 71'
	write(4790,*)'198010.514 198209.3984'
	write(4790,*)'0 +0.01'
   	write(4790,*)'1 100'

c     constantes numericas
      ikepler=1
      inicio=1       
      kk=1
      pi=4.0d0*datan(1.0d0)      
      pi2=pi+pi
      pi05=0.5d0*pi
      pi15=1.5d0*pi
      conv=pi/180.0d0
      
c**** DADOS DE ENTRADA                                           
c      write(*,*)'Numero de planetas:'
c      read(*,*) N
      N=5
      neq=N+N+N+N+N+N

c	LOOPS

	ei=0.01d0
	ef=0.05d0
	de=(ef-ei)/50.0d0

c	valor medio atual=238040/60268=3.94969137850932
c	valor do artigo

	semi0=198100.d0/60268d0
	semif=198200.d0/60268d0

c	diferen� de 0.030 R_S, ou ~1808 km

	dsemi=(semif-semi0)/50.0d0

c	iu=1	
	do exc=ei,ef+0.000001,de
	do semi=semi0,semif+0.0000001,dsemi
c	write(*,*)exc,semi

c	pause
c	write(*,*)iu
c	iu=iu+1
c	pause

	aa=semi
	ee=exc
      CALL ENTRE(aa,ee,dm,ami,dmimi,a,e,di,w,om,am)

c	COME�  com  emax=excentricidade inicial, amax=..., imax=...
      emax=exc
	amax=semi
	dimax=di(2)

c**** TRANFORMACAO DE ORBITAIS P/ CARTESIANOS DAS CONDICOES INICIAIS      
      CALL ORBXYZ(N,a,e,di,w,om,am,ami,xpla,f)
c      write(*,*)'call orbxyz'                 
c      pause

cccccccccccccccccccccccccccccc     Definicoes: 

c     COORDENADAS a serem utilizadas na subroutine FORCE (VER orbxyz): 
      i=0   
      j=0   
      k=1
      do 2 i=1,N
      do 1 j=4,6

      x(k)=xpla(i,j)      
c      write(2,*)'x coordenadas',k,x(k)
      k=k+1
1     continue
      k=k+3
2     continue
	
c     Velocidades relativas, que serao utilizadas em gaussj para 
c     calcular os momenta (ver orbxyz)
c     Vel. eixo x
      i=0
      k=1
      do 4 i=1,N
      bx(k,1)=xpla(i,1)
c      write(*,*)'vel x',k,bx(k,1)
      k=k+1                       
4     continue
       
      CALL gaussj(N,bx,dmimi)        
	
c      write(*,*)'call gaussj X'                 
c      pause

c     Vel. eixo y             
      i=0
      k=1
      do 41 i=1,N
      by(k,1)=xpla(i,2)
c      write(*,*)'vel y',k,by(k,1)
      k=k+1                       
41     continue
      
      CALL gaussj(N,by,dmimi)
c      write(*,*)'call gaussj Y'                 
c      pause

c     Vel. eixo z             
      i=0
      k=1
      do 42 i=1,N
      bz(k,1)=xpla(i,3)
c      write(2,*)'vel z',k,bz(k,1)
      k=k+1                       
42     continue
      
      CALL gaussj(N,bz,dmimi)
c      write(*,*)'call gaussj Z'                 
c      pause


c     MOMENTA canonicos (VER FORCE):      
      k=4
      l=1
5     continue      

      x(k)=bx(l,1)    
C      write(2,*)'x momenta',k,bx(l,1)
      x(k+1)=by(l,1)
C      write(2,*)'x momenta',k+1,by(l,1)      
      x(k+2)=bz(l,1)
C      write(2,*)'x momenta',k+2,bz(l,1)

      k=k+6      
      l=l+1
      if(l.gt.N) go to 55
      go to 5

55    continue

c     Devemos redefinir a matriz xpla(8,6):
      idef=1
      jdef=1
      j=1
      do 25 idef=1,N
      do 24 jdef=1,3

C      xpla(idef,jdef)=F1(j)
C     modifiquei agora: quero os elementos com as velocidades absolutas
c     para comparar com a Tatiana.
c     Agora temos

      xpla(idef,jdef)=x(j+3)/dm(idef)
c      WRITE(2,*)'XPLA(I,J)',IDEF,JDEF,xpla(idef,jdef)
      xpla(idef,jdef+3)=x(j)
c      WRITE(2,*)'XPLA(I,J)',IDEF,JDEF+3,xpla(idef,jdef+3)
      j=j+1
24     continue
      j=j+3
25    continue

cccccccccccccccccccccccccccccc  FIM   Definicoes: 
      CALL XYZORB(N,xpla,a,e,di,w,om,am,ami,lei,kk) 

c	XA1FFT(1)=a(1)
c	XI1FFT(1)=di(1)

	XAFFT(1)=a(2)
	XEFFT(1)=e(2)
	XIFFT(1)=di(2)
	XE2FFT(1)=e(2)*dcos(w(2)+om(2))

c	XA3FFT(1)=a(3)
c	XI3FFT(1)=di(3)

c	XA4FFT(1)=a(4)

c	XA5FFT(1)=a(5)
c	XE5FFT(1)=e(5)

c	XA6FFT(1)=a(6)


c	CALL ESCREVE1(a(1),e(1),di(1),w(1),om(1),am(1),a(2),e(2),
c	<di(2),w(2),om(2),am(2),ami(1),ami(2),t)


c	CALL ESCREVE1(a(2),e(2),di(2),w(2),om(2),am(2),a(4),e(4),
c	<di(4),w(4),om(4),am(4),ami(2),ami(4),t)

c	CALL ESCREVE2(a(1),e(1),di(1),w(1),om(1),am(1),a(3),e(3),
c	<di(3),w(3),om(3),am(3),ami(1),ami(3),t)

c	CALL ESCREVE3(a(5),e(5),di(5),w(5),om(5),am(5),a(6),e(6),
c	<di(6),w(6),om(6),am(6),ami(6),ami(6),t)

                  
c      CALL ENERGIA(H,X,dmimi,dm)
c      write(*,*)'energia 0'                 
c      pause      
c      H=H*100000.0D0      
c      write(5,*)t,H
c      write(*,*)H
c      write(*,*)
	      
c**** INTEGRACAO NUMERICA E SAIDA DOS DADOS
	jjj=2
      do while (t.lt.tint)      
      CALL RA15(TM,X,V,TF,10,XL,neq,1)      
C      write(*,*)'call ra15'                 
C      pause

      t=t+tf   
c     Neste ponto, ja' tenho as novas x(1),...,x(24) (coordenadas e momenta)

c     Deverei transformar agora as novas COORDENADAS e as novas VELOCIDADES     
C     RELATIVAS em elem. orbitais. As velocidades relativas sao fornecidas na
c     propria sub. FORCE (VER FORCE).
                               
c       CALL FORCE(TM,X,V,F1) 
c      write(*,*)'call force'                 
C      pause

  
c     Devemos redefinir a matriz xpla(8,6):
      idef=1
      jdef=1
      j=1
      do 10 idef=1,N
      do 9  jdef=1,3

C      xpla(idef,jdef)=F1(j)
C     modifiquei agora: quero os elementos com as velocidades absolut
c     para comparar com a Tatiana.
c     Agora temos

      xpla(idef,jdef)=x(j+3)/dm(idef)	
c      WRITE(2,*)'XPLA(I,J)',IDEF,JDEF,xpla(idef,jdef)            
      xpla(idef,jdef+3)=x(j)
c      WRITE(2,*)'XPLA(I,J)',IDEF,JDEF+3,xpla(idef,jdef+3)         
      j=j+1
9     continue
      j=j+3
10    continue


      CALL XYZORB(N,xpla,a,e,di,w,om,am,ami,lei,kk)

      IF(e(2).GT.EMAX) then
	EMAX=e(2)
	TTT=T
      ENDIF 

      IF(a(2).GT.AMAX) then
	AMAX=a(2)
	TTT=T
      ENDIF 

      IF(di(2).GT.DIMAX) then
	DIMAX=di(2)
	TTT=T
      ENDIF 
   

C	vetor de dados para FFT

c	XA1FFT(jjj)=a(1)
c	XI1FFT(jjj)=di(1)

	XAFFT(jjj)=a(2)
	XEFFT(jjj)=e(2)
	XIFFT(jjj)=di(2)
	XE2FFT(jjj)=e(2)*dcos(w(2)+om(2))

c	XA3FFT(jjj)=a(3)
c	XI3FFT(jjj)=di(3)

c	XA4FFT(jjj)=a(4)

c	XA5FFT(jjj)=a(5)
c	XE5FFT(jjj)=e(5)

c	XA6FFT(jjj)=a(6)
	jjj=jjj+1


c	CALL ESCREVE1(a(1),e(1),di(1),w(1),om(1),am(1),a(2),e(2),
c	<di(2),w(2),om(2),am(2),ami(1),ami(2),t)


c	CALL ESCREVE1(a(2),e(2),di(2),w(2),om(2),am(2),a(4),e(4),
c	<di(4),w(4),om(4),am(4),ami(2),ami(4),t)

c	CALL ESCREVE2(a(1),e(1),di(1),w(1),om(1),am(1),a(3),e(3),
c	<di(3),w(3),om(3),am(3),ami(1),ami(3),t)

c	CALL ESCREVE3(a(5),e(5),di(5),w(5),om(5),am(5),a(6),e(6),
c	<di(6),w(6),om(6),am(6),ami(6),ami(6),t)

c      CALL ENERGIA(H,X,dmimi,dm)
c      write(*,*)'energia'                 
c      pause      
c      H=H*100000.0D0      
c      write(5,*)t,H
c      write(*,*)H
c      WRITE(*,*)
c	write(522,*)t,a(2),a(1)

      end do      

c99	Format(F11.3,1x,F12.5,1x,F9.6,1x,F12.8,1x,F10.8)


	write(923,*)ttt,ttt/365.25d0,semi,exc,emax
	write(9234,*)ttt,ttt/365.25d0,semi,exc,amax
	write(92345,*)ttt,ttt/365.25d0,semi,exc,dimax/conv

	write(123,*)semi,exc,emax
	write(1234,*)semi,exc,amax
	write(12345,*)semi,exc,dimax/conv

	write(1239,*)emax
	write(12349,*)amax
	write(23459,*)dimax/conv


	T=0.0D0

c	F  F  T

c     CALL  FFT(XA1FFT,XI1FFT,XA2FFT,XE2FFT,XA3FFT,
c     >XI3FFT,XA4FFT,XA5FFT,XE5FFT,XA6FFT)   
	iii=1
	CALL FFT(XAFFT,nespec01,nespec05,nespec1,nespec3,
     >nespec5,nespec10,nespec15,exc,semi,iii,ioi)

	write(2,*)nespec01
	write(22,*)exc,semi,nespec01

	write(3,*)nespec05
	write(33,*)exc,semi,nespec05

	write(4,*)nespec1
	write(44,*)exc,semi,nespec1

	write(41,*)nespec3
	write(441,*)exc,semi,nespec3

	write(5,*)nespec5
	write(55,*)exc,semi,nespec5

	write(622,*)nespec10
	write(662,*)exc,semi,nespec10

	write(7,*)nespec15
	write(77,*)exc,semi,nespec15
	ioi=0

	iii=2
	CALL FFT(XEFFT,nespec01,nespec05,nespec1,nespec3,
     >nespec5,nespec10,nespec15,exc,semi,iii,ioi)

	write(29,*)nespec01
	write(229,*)exc,semi,nespec01

	write(39,*)nespec05
	write(339,*)exc,semi,nespec05

	write(49,*)nespec1
	write(449,*)exc,semi,nespec1

	write(491,*)nespec3
	write(4491,*)exc,semi,nespec3

	write(59,*)nespec5
	write(559,*)exc,semi,nespec5

	write(6229,*)nespec10
	write(6629,*)exc,semi,nespec10

	write(79,*)nespec15
	write(779,*)exc,semi,nespec15

	iii=3
	CALL FFT(XIFFT,nespec01,nespec05,nespec1,nespec3,
     >nespec5,nespec10,nespec15,exc,semi,iii,ioi)

	write(290,*)nespec01
	write(2290,*)exc,semi,nespec01

	write(390,*)nespec05
	write(3390,*)exc,semi,nespec05

	write(490,*)nespec1
	write(4490,*)exc,semi,nespec1

	write(4901,*)nespec3
	write(44901,*)exc,semi,nespec3

	write(590,*)nespec5
	write(5590,*)exc,semi,nespec5

	write(62290,*)nespec10
	write(66290,*)exc,semi,nespec10

	write(790,*)nespec15
	write(7790,*)exc,semi,nespec15


	iii=4
	CALL FFT(XE2FFT,nespec01,nespec05,nespec1,nespec3,
     >nespec5,nespec10,nespec15,exc,semi,iii,ioi)

	write(4290,*)nespec01
	write(42290,*)exc,semi,nespec01

	write(4390,*)nespec05
	write(43390,*)exc,semi,nespec05

	write(4490,*)nespec1
	write(44490,*)exc,semi,nespec1

	write(44901,*)nespec3
	write(444901,*)exc,semi,nespec3

	write(4590,*)nespec5
	write(45590,*)exc,semi,nespec5

	write(462290,*)nespec10
	write(466290,*)exc,semi,nespec10

	write(4790,*)nespec15
	write(47790,*)exc,semi,nespec15



c	FIM LOOPS
	
	enddo
	enddo



      
      END                               
cccccccccccccccccccccccccccccccccccccccccccccccccc
      SUBROUTINE FFT(XE2FFT,nespec01,nespec05,nespec1,nespec3,
     >nespec5,nespec10,nespec15,exc,semi,iii,ioi)
C	(XA1FFT,XI1FFT,XA2FFT,XE2FFT,XA3FFT,
C     >XI3FFT,XA4FFT,XA5FFT,XE5FFT,XA6FFT)   
      INTEGER n,nn,i,ii,j,kk,isign,iii,ioi
      REAL*8 data(2097152)

      REAL*8 XE2FFT(1048576)
C	XA1FFT(1048576),XI1FFT(1048576),XA2FFT(1048576),
C     >XE2FFT(1048576),XA3FFT(1048576),XI3FFT(1048576),
C     >XA4FFT(1048576),XA5FFT(1048576),XE5FFT(1048576),XA6FFT(1048576)

      REAL*8 FRE(1048576),AMP(1048576),FFRE(1048576)
      REAL*8 AAMP(1048576),ff(1048576),gg(1048576)
      REAL*8 NLINHAS,DELTAT,ttf,tint
      REAL*8 tf,pi,conv,amedio,emedio,dimedio,semi,exc
	INTEGER nespec01,nespec05,nespec1,nespec3,
     >nespec5,nespec10,nespec15,dn01,dn05,dn1,dn3,dn5,dn10,dn15
      REAL*8 a,a1,a2,ampl,ampmax,amplim,peri
      common/tempo/tint,ttf,Nlinhas


      open (100,file='Amedio.DAT', status='unknown') 	
      open (200,file='Emedio.DAT', status='unknown') 	
      open (300,file='Imedio.DAT', status='unknown') 			
      open (400,file='EcosVmedio.DAT', status='unknown') 			

      open (171,file='Amedio.grd', status='unknown') 	
      open (272,file='Emedio.grd', status='unknown') 	
      open (373,file='Imedio.grd', status='unknown') 			
      open (473,file='EcosVmedio.grd', status='unknown') 			

      open (400,file='AmedioEmedioEVmedioImedio.DAT', status='unknown') 	


	if(ioi.eq.79) then

	write(171,'(4hDSAA)')
	write(171,*)'51 71'
	write(171,*)'198010.514 198209.3984'
	write(171,*)'0 +0.01'
   	write(171,*)'1 100'

	write(272,'(4hDSAA)')
	write(272,*)'51 71'
	write(272,*)'198010.514 198209.3984'
	write(272,*)'0 +0.01'
   	write(272,*)'1 100'

	write(373,'(4hDSAA)')
	write(373,*)'51 71'
	write(373,*)'198010.514 198209.3984'
	write(373,*)'0 +0.01'
   	write(373,*)'1 100'

	write(473,'(4hDSAA)')
	write(473,*)'51 71'
	write(473,*)'198010.514 198209.3984'
	write(473,*)'0 +0.01'
   	write(473,*)'1 100'
	endif


     	deltat=ttf
                 
      nn=IDINT(Nlinhas)  
c      write(*,*)nn
         
c      dnn=2*nn


c******************************************************

c	FFT DE E2.
	KK=1
	j=0      

      i=1
      n=0
	do while(n.le.nn-1)

      data(i)=XE2FFT(KK)
	data(i+1)=0.0d0
	
      if(n.eq.0) gg(n+1)=0.0d0
	if(n.ge.1. and. n.le.nn/2-1) gg(n+1)=n/(nn*deltat)
	if(n.eq.nn/2) gg(n+1)=1.0d0/(2.0d0*deltat)
	if(n.ge.nn/2+1. and. n.le.nn-1) then
	ff(n+1)=-gg(n-1-2*j)
	j=j+1
	endif 

	n=n+1      
	i=i+2     
	KK=KK+1
	end do
        
      isign=1                               
      CALL FOUR1(data,nn,isign) 

	ii=1
      n=0
	i=1
	do while(n.le.nn-1)
	tf=2.0d0*((data(i)**2+data(i+1)**2)**0.5d0)/nn

      pi=4.0d0*datan(1.0d0)      
      conv=pi/180.0d0
 
cccccccccccccccccccccccc
      if(n.eq.0) then
	
	if(iii.eq.1) then
	write(100,*) semi,exc,tf/2.0d0
	write(171,*) tf/2.0d0	
	amedio=tf/2.0d0
	endif
	
	if(iii.eq.2) then
	write(200,*) semi,exc,tf/2.0d0
	write(272,*) tf/2.0d0	
	emedio=tf/2.0d0
	endif

	if(iii.eq.3) then 
	write(300,*) semi,exc,(tf/2.0d0)/conv
	write(373,*) (tf/2.0d0)/conv
	dimedio=tf/2.0d0
	endif

	if(iii.eq.4) then 
	write(400,*) semi,exc,(tf/2.0d0)/conv
	write(473,*) (tf/2.0d0)/conv
	dEvmedio=tf/2.0d0
	endif


c199	Format(F11.8,1x,F10.8,1x,F10.8,1x,F12.8,1x,F9.6,1x,F12.8)
	write(400,*)amedio,emedio,dEvmedio,dimedio/conv,semi,exc

	endif
ccccccccccccccccccccccccc


	if(n.ge.1.and.n.le.nn/2-1) THEN
	FRE(ii)=gg(n+1)
	AMP(ii)=tf
c	write(25,*)fre(ii),1.0D0/fre(ii),amp(ii)
	ENDIF

	if(n.eq.nn/2) THEN
	FRE(ii)=gg(n+1)
	AMP(ii)=tf
c	write(25,*)fre(ii),1.0D0/fre(ii),amp(ii)
	ENDIF


      n=n+1
      i=i+2
	ii=ii+1


	end do

cccccccccc PROCURA PELOS PICOS DE FREQUENCIA 
	J=1
	do jj=2,ii
	
	a1=AMP(jj-1)
	a =AMP(jj)
	a2=AMP(jj+1)

	if(a.gt.a1. and. a2.lt.a) then
	AAMP(J)=AMP(jj)
	FFRE(J)=FRE(jj)
C	write(*,*)AAMP(J),FFRE(J)

c	no ultimo loop, o numero de picos e' determinado aqui (J);
c	note que aamp(j) e' dado pelo numero de picos (importante), 
C	E NAO O NUMERO TOTAL DE FREQS.

	J=J+1
	endif

	enddo
C	WRITE(*,*)J

c	devo fazer J=J-1 pois o numero de picos FOI SOMADO DE UM ACIMA 

	J=J-1

cccccccccc PROCURA PELA amplitude MAXIMA

	JJ=1
	ampmax=-1.

	do while(JJ.le.J)

c	FILTRO
	peri=1.0d0/FFRE(JJ)

cccccccccccccccccccc
	ampl=AAMP(JJ)
	if(peri.gt.15.0d0. and. ampl.GT.ampmax) then
        ampmax=ampl
	endif
cccccccccccccccccccc

	JJ=JJ+1

	enddo

	amplim=ampmax/100.0d0
C	write(*,*)ampmax

C	WRITE(*,*)JJ

c	o  mesmo que acima (fazer jj=jj-1)

	JJ=JJ-1

cccccccccccccc 

C	JJ e' o numero total de picos

      dn01=0
      dn05=0
      dn1=0
	dn3=0
      dn5=0
      dn10=0
      dn15=0
      
	J=1
	do while(J.le.JJ)

	ampl=AAMP(J)

	peri=1.0d0/FFRE(J)

	if(peri.gt.15.0d0.and.  ampl.gt.( 0.1d0*amplim)  ) dn01=dn01+1
	if(peri.gt.15.0d0.and.  ampl.gt.( 0.5d0*amplim)  ) dn05=dn05+1
	if(peri.gt.15.0d0.and.  ampl.gt.( 1.0d0*amplim)  ) dn1 =dn1 +1
	if(peri.gt.15.0d0.and.  ampl.gt.( 3.0d0*amplim)  ) dn3 =dn3 +1
	if(peri.gt.15.0d0.and.  ampl.gt.( 5.0d0*amplim)  ) dn5 =dn5 +1
	if(peri.gt.15.0d0.and.  ampl.gt.(10.0d0*amplim)  ) dn10=dn10+1
	if(peri.gt.15.0d0.and.  ampl.gt.(15.0d0*amplim)  ) dn15=dn15+1

	J=J+1

c	Numeros espectrais: definidos abaixo (USO O ULTIMO no artigo)


      nespec01=dn01
      nespec05=dn05
      nespec1 =dn1
      nespec3 =dn3
      nespec5 =dn5
      nespec10=dn10
      nespec15=dn15


	enddo


c	FIM FFT EM E2.

c******************************************************


	RETURN
	END
cccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccc                 
      SUBROUTINE FORCE (TM,X,V,F1)
      IMPLICIT REAL *8(a-h,o-z)
      DIMENSION X(36),V(36),F1(36)
      DIMENSION DDJ2(3),DDJ4(3)
	INTEGER II
      dimension a(6),e(6),di(6),w(6),om(6),am(6),dm(6),ami(6),dmimi(6)
      COMMON/MASSAS/dm0,G,dj2,dj4
      COMMON/NUMERO/N,neq
	common/novo/aa,ee

c*****EQUACOES DE MOVIMENTO CANONICAS PARA N PLANETAS*****
c          EM VARIAVEIS HELIOCENTRICAS CANONICAS

c     Coordenadas (relativas ao Sol)
c     Exemplo: Mercurio e Venus
c      x(1)=xmer, x(2)=ymer, x(3)=zmer
c      x(7)=xven, x(8)=yven, x(9)=zven

c     Momenta (C.M.):
c      x(4)=pxmer, x(5)=pymer, x(6)=pzmer
c      x(10)=pxven, x(11)=pyven, x(12)=pzven


      CALL ENTRE(aa,ee,dm,ami,dmimi,a,e,di,w,om,am)
C       write(*,*)'call entre force'                 

C      pause
      l=1
      i=1
      
1     continue
      if(l+5.gt.neq. and .i.gt.N) go to 2                      
C	WRITE(*,*)I
C	PAUSE
      F1(l)=dmimi(i)*x(l+3)
      
      F1(l+1)=dmimi(i)*x(l+4)
      
      F1(l+2)=dmimi(i)*x(l+5)

      j=1     
      smx=0.0d0
      smy=0.0d0
      smz=0.0d0
      do while(j+5.le.neq)
      if (j.eq.l) go to 21
      smx=smx+x(j+3)
      smy=smy+x(j+4)
      smz=smz+x(j+5)
21    continue
      j=j+6
      enddo  

      F1(l)=F1(l)+smx/dm0
            
      F1(l+1)=F1(l+1)+smy/dm0
      
      F1(l+2)=F1(l+2)+smz/dm0

 
      ri3=(x(l)**2+x(l+1)**2+x(l+2)**2)**(1.5d0)      
c      write(*,*)'ri3',ri3
      
      F1(l+3)=-G*dm(i)*dm0*x(l)/ri3

      F1(l+4)=-G*dm(i)*dm0*x(l+1)/ri3
      
      F1(l+5)=-G*dm(i)*dm0*x(l+2)/ri3      


      j=1
      do while(j+5.le.neq)
    
      if(j.eq.l) go to 12
      if (j.eq.1)  m=1
      if (j.eq.7)  m=2
      if (j.eq.13) m=3
      if (j.eq.19) m=4
      if (j.eq.25) m=5
            
      dij3=( (x(l)-x(j))**2 + (x(l+1)-x(j+1))**2 
     | + (x(l+2)-x(j+2))**2 )**(1.5d0)          
c      write(*,*)'dij3',dij3      
          
      F1(l+3)=F1(l+3)-G*dm(i)*dm(m)*(x(l)-x(j))/dij3
      
      F1(l+4)=F1(l+4)-G*dm(i)*dm(m)*(x(l+1)-x(j+1))/dij3

      F1(l+5)=F1(l+5)-G*dm(i)*dm(m)*(x(l+2)-x(j+2))/dij3

12    continue
      j=j+6
      enddo      
 
	II=i

	CALL ACHAT(dm,X,DDJ2,DDJ4,II)

	dj2x=DDJ2(1)
	dj2y=DDJ2(2)
	dj2z=DDJ2(3)

	dj4x=DDJ4(1)
	dj4y=DDJ4(2)
	dj4z=DDJ4(3)

	F1(l+3)=F1(l+3)-dj2x-dj4x
      F1(l+4)=F1(l+4)-dj2y-dj4y
      F1(l+5)=F1(l+5)-dj2z-dj4z
C	write(*,*)I,F1(l+3)

      l=l+6
      i=i+1
      go to 1                           
2     continue      

      
      return
      end      
     
cccccccccccccccccccccccccccccccccccccc                 
      SUBROUTINE ACHAT(dm,X,DDJ2,DDJ4,II)
      IMPLICIT REAL *8(a-h,o-z)
      DIMENSION dm(6),X(36),DDJ2(3),DDJ4(3)
	INTEGER II
      COMMON/MASSAS/dm0,G,dj2,dj4

	if(II.eq.1) then
      x1=x(1)
      y1=x(2)
      z1=x(3)
      r1=dsqrt(x1**2.0d0+y1**2.0d0+z1**2.0d0)
	dMm=dm(1)
c	write(*,*)x1,y1,z1
c	pause
	endif

	if(II.eq.2) then
      x1=x(7)
      y1=x(8)
      z1=x(9)
      r1=dsqrt(x1**2.0d0+y1**2.0d0+z1**2.0d0)
	dMm=dm(2)
c	write(*,*)dmm
c	pause
	endif

	if(II.eq.3) then
      x1=x(13)
      y1=x(14)
      z1=x(15)
      r1=dsqrt(x1**2.0d0+y1**2.0d0+z1**2.0d0)
	dMm=dm(3)
	endif

	if(II.eq.4) then
      x1=x(19)
      y1=x(20)
      z1=x(21)
      r1=dsqrt(x1**2.0d0+y1**2.0d0+z1**2.0d0)
	dMm=dm(4)
	endif

	if(II.eq.5) then
      x1=x(25)
      y1=x(26)
      z1=x(27)
      r1=dsqrt(x1**2.0d0+y1**2.0d0+z1**2.0d0)
	dMm=dm(5)
	endif

	if(II.eq.6) then
      x1=x(31)
      y1=x(32)
      z1=x(33)
      r1=dsqrt(x1**2.0d0+y1**2.0d0+z1**2.0d0)
	dMm=dm(6)
	endif

c     Derivadas que entram no achatamento (PLANETA 1 (I-ESIMO PLANETA))

c     coordenada x

      r17=r1**7.0d0
      DPX1=-5.0d0*x1/r17
      
      r15=r1**5.0d0
      DPX2=-3.0d0*x1/r15

      r11=r1**11.0d0
      DPX3=-9.0d0*x1/r11

      r19=r1**9.0d0
      DPX4=-7.0d0*x1/r19
      
c     coordenada y

      DPY1=-5.0d0*y1/r17

      DPY2=-3.0d0*y1/r15

      DPY3=-9.0d0*y1/r11

      DPY4=-7.0d0*y1/r19

c     coordenada z

      DPZ1=-5.0d0*z1/r17

      DPZ2=-3.0d0*z1/r15

      DPZ3=-9.0d0*z1/r11

      DPZ4=-7.0d0*z1/r19
c	write(*,*)i,x1,y1,z1
c	pause

c     Termos que aparecem do achatamento (PLANETA 1 (I-ESIMO PLANETA)	 )

      cj2=G*dm0*dMm*dj2/2.0d0
      cj4=G*dm0*dMm*dj4/8.0d0

c 	write(*,*)ii,g,dm0,dj2,dj4
c	pause

      dj2x1=3.0d0*z1*z1*DPX1-DPX2
	DDJ2(1)=cj2*dj2x1
c	write(*,*)ii,ddj2(1)
c	pause

      dj2y1=3.0d0*z1*z1*DPY1-DPY2
	DDJ2(2)=cj2*dj2y1

      dj2z1=3.0d0*( z1*z1*DPZ1+2.0d0*z1/(r1**5.0d0) )-DPZ2
	DDJ2(3)=cj2*dj2z1

cccccccccccc

      dj4x1=35.0d0*(z1**4.0d0)*DPX3-
     -      30.0D0*(z1**2.0d0)*DPX4+
     +      3.0d0*DPX1

	DDJ4(1)=cj4*dj4x1

cccccccccccc

      dj4y1=35.0d0*(z1**4.0d0)*DPY3-
     -      30.0D0*(z1**2.0d0)*DPY4+
     +      3.0d0*DPY1

	DDJ4(2)=cj4*dj4y1

cccccccccccc

      dj4z1=35.0d0*( (z1**4.0d0)*DPZ3+4.0d0*(z1**3.0d0)/(r1**9.0d0) )-
     ~      30.0D0*( (z1**2.0d0)*DPZ4+2.0d0*z1/(r1**7.0d0) )+
     ~      3.0d0*DPZ1

	DDJ4(3)=cj4*dj4z1


      return
      end      

ccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccc      

ccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccc      
c      CALL ENTRE      (dm,ami,dmimi,a,e,di,w,om,am)
      SUBROUTINE ENTRE(aa,ee,dm,ami,dmimi,a,e,di,w,om,am)
      implicit real *8(a-h,o-z)
      dimension a(6),e(6),di(6),w(6),om(6),am(6),dm(6),ami(6),dmimi(6)
      COMMON/CONST/conv
      COMMON/MASSAS/dm0,G,dj2,dj4
c:    Esta rotina informa as constantes de cada planeta
c:    Este grav em UA, Massa Solar e Ano
      
ccccc      G=39.476926421373020d0      
ccccc      dm0=1.0000000000000002d0

	dm0=1.0d0
      pi=4.0d0*datan(1.0d0)      
	grav=4.0d0*pi*pi
c	g=grav

c	HORIZONS DATA SYSTEM PARA DATA 2007-JAN 01- 00:00
      ua=1.49597870691d11                                                   metros

      requat=6.0268d7
      esc=ua/requat

c	em	quilos
c	dmsatt=1.9891d30/5.68462313752d26
c      write(*,*)dmsat
c	=3499.08859722190

c	JACOBSON 2006
	dmsun=3497.9018d0*37940585.2d0/6.6742d-23
c	write(*,*)dmsun
c	pause

c	Esta quantidade e' 1/M
	dmsatt=dmsun/(37931207.7d0/6.6742d-23)
c      write(*,*)dmsatt,i
c	pause

      grav=grav*esc*esc*esc
	G=grav/(dmsatt*365.25d0*365.25d0)

c	JACOBSON 2005
c      dj2= 0.0162906d0
c	dj4=-0.000936d0

C	J2, J4
C	JACOBSON 2006
      dj2= 0.01629071d0
	dj4=-0.0009358d0


c     grav=G      
c     ua=1.4959787d11  metros        
 
c:    GRAV acima em UA**3/(Massa Solar*Ano**2)
c:    UA acima em metros  e Raio equatorial=requat sera em METROS (aqui)




c Mimas
      dm(1)=3.75d22/5.68462313752d29

      ami(1)=G*(dm0+dm(1))
      dmimi(1)=(dm0+dm(1))/(dm0*dm(1))
      a(1)=1.860256305913010d+08/requat
      a(1)=1.860256453660664d+08/requat
      e(1)=2.002855413732700d-02
      e(1)=2.002833386117072d-02
      di(1)=1.568081505522592d+00*conv
      di(1)=1.567961379446476d+00*conv
      w(1)=1.496697121023191d+02*conv
      w(1)=1.497709076854479d+02*conv
      om(1)=9.268827186791158d+01*conv
      om(1)=9.264771224778193d+01*conv
      am(1)=7.298297803108471d+01*conv
      am(1)=7.291866912551491d+01*conv


c Anthe
      dm(2)= dm(1)*(1.5d0/200.0d0)**3.0d0
c10.805d22/5.68462313752d29
      ami(2)=G*(dm0+dm(2))
      dmimi(2)=(dm0+dm(2))/(dm0*dm(2))
      a(2)=aa
c            dperiodo=2.0d0*pi*dsqrt( a(2)**3.0d0/ami(2) )
c 	write(*,*)dperiodo
c	pause
      e(2)=ee

      di(2)=2.116034272544871d-02*conv
      om(2)= 2.844366472947400d+02*conv
      w(2)=7.773816459991164d+01*conv
      am(2)=3.407353235645896d+02*conv


c Enceladus
      dm(3)=10.805d22/5.68462313752d29
      ami(3)=G*(dm0+dm(3))
      dmimi(3)=(dm0+dm(3))/(dm0*dm(3))
      a(3)= 2.384109587128550d+08/requat
      a(3)= 2.384109510175851d+08/requat
c     0.001593665160605842d0*ua/requat
      e(3)=5.409346486942390d-03
      e(3)=5.409346486942390d-03
      di(3)=6.318310116077915d-03*conv
      di(3)=6.363503564759852d-03*conv
      w(3)=2.530723574939369d+02*conv
      w(3)=2.551978022439982d+02*conv
      om(3)=1.086311106169690d+02*conv
      om(3)=1.065048182831290d+02*conv
      am(3)=5.419323397767205d+01*conv
      am(3)=5.419404329277588d+01*conv

c Tethys

      dm(4)= 61.76d22/5.68462313752d29
      ami(4)=G*(dm0+dm(4))
      dmimi(4)=(dm0+dm(4))/(dm0*dm(4))
      a(4)=2.949751551802501d+08/requat
      a(4)=2.949751474610816d+08/requat

      e(4)= 1.070242549617775d-03
      e(4)= 1.064627114234757d-03

      di(4)=1.093915959700871d+00*conv
      di(4)=1.093729066055329d+00*conv

      w(4)=1.554163514461071d+02*conv
      w(4)=1.553196746160572d+02*conv

      om(4)=1.837802356190735d+02*conv
      om(4)=1.837818434564813d+02*conv

      am(4)=3.538882223978986d+02*conv
      am(4)=3.539836991737204d+02*conv


c Dione
       dm(5)=109.572d22/5.68462313752d29


      ami(5)=G*(dm0+dm(5))
      dmimi(5)=(dm0+dm(5))/(dm0*dm(5))
      a(5)=3.776514622766157d+08/requat
      a(5)=3.776514606708938d+08/requat

c      0.002524414137313330d0*ua/requat
      e(5)=1.723691030913839d-03
      e(5)=1.723734703752217d-03

      di(5)=2.983993394557515d-02*conv
      di(5)=2.980714121136082d-02*conv

      w(5)=1.622056727265864d+02*conv
      w(5)=1.620327193810711d+02*conv

      om(5)=1.708605061184919d+02 *conv
      om(5)=1.710361740602899d+02 *conv

      am(5)=2.282374048754474d+02*conv
      am(5)=2.282347444217630d+02*conv
      return
	end       

cccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccc
c      SUBROUTINE ENERGIA(H,X,dmimi,dm)
c      IMPLICIT REAL *8(a-h,o-z)
c      DIMENSION X(24),dmimi(4),dm(4)
c      COMMON/MASSAS/dm0,G,dj2,dj4
c      COMMON/NUMERO/N,neq
      
c      l=1  
c      U0=0.0d0
c      T0=0.0d0        
c      do while(l+5.le.neq)

c      if (l.eq.1)  m=1
c      if (l.eq.7)  m=2
c      if (l.eq.13) m=3
c      if (l.eq.19) m=4

c      ri=dsqrt(x(l)**2+x(l+1)**2+x(l+2)**2)
c      dmom2=x(l+3)**2+x(l+4)**2+x(l+5)**2 
c      U0=U0-dm(m)/ri                               
cc      write(4,*)'u0',u0
c      T0=T0+(dmom2*dmimi(m))
cc      write(4,*)'t0',t0
c      l=l+6

c      enddo
c      H0=G*dm0*U0+T0/2.0d0      
cc      write(4,*)'h0',h0
c      T1=0.0d0
cc      U1=0.0d0
c      j=7
c      do while(j+5.le.neq)

c      if (j.eq.7)  mj=2
c      if (j.eq.13) mj=3
c      if (j.eq.19) mj=4       
c
c      d1j=dsqrt( (x(1)-x(j))**2 + (x(2)-x(j+1))**2 
c     | + (x(3)-x(j+2))**2 )          
c      U1=U1-dm(1)*dm(MJ)/d1j
cc      write(4,*)'u1',u1
c      T1=T1+x(4)*x(j+3)+x(5)*x(j+4)+x(6)*x(j+5)
cc      write(4,*)'t1',t1      
c
c      j=j+6
c      enddo
c c
c      
c      H1=G*U1+T1/dm0      
cc      write(4,*)'h1',h1
c      H=H0+H1 
cc      write(4,*)'h',h
c      return
c      end
cccccccccccccccccccccccccccccccccccccccccccccccccccc
      SUBROUTINE four1(data,nn,isign)
      INTEGER isign,nn
c      REAL data(2*nn)
	REAL*8 data(2*nn)
      INTEGER i,istep,j,m,mmax,n
c      REAL*8 tempi,tempr
      DOUBLE PRECISION theta,wi,wpi,wpr,wr,wtemp,tempi,tempr
      n=2*nn
      j=1
      do 11 i=1,n,2
        if(j.gt.i)then
          tempr=data(j)
          tempi=data(j+1)
          data(j)=data(i)
          data(j+1)=data(i+1)
          data(i)=tempr
          data(i+1)=tempi
        endif
        m=n/2
1       if ((m.ge.2).and.(j.gt.m)) then
          j=j-m
          m=m/2
        goto 1
        endif
        j=j+m
11    continue
      mmax=2

2     if (n.gt.mmax) then
        istep=2*mmax
        theta=6.28318530717959d0/(isign*mmax)
        wpr=-2.d0*sin(0.5d0*theta)**2
        wpi=sin(theta)
        wr=1.d0
        wi=0.d0
        do 13 m=1,mmax,2
          do 12 i=m,n,istep
            j=i+mmax
            tempr=sngl(wr)*data(j)-sngl(wi)*data(j+1)
            tempi=sngl(wr)*data(j+1)+sngl(wi)*data(j)
            data(j)=data(i)-tempr
            data(j+1)=data(i+1)-tempi
            data(i)=data(i)+tempr
            data(i+1)=data(i+1)+tempi
12        continue
          wtemp=wr
          wr=wr*wpr-wi*wpi+wr
          wi=wi*wpr+wtemp*wpi+wi
13      continue
        mmax=istep

      goto 2
      endif
      return
      END


cccccccccccccccccccccccccccccccccccccccccccc     


ccccccccccccccccccccccccccccccccccccccccccccccccccccccccccc

ccccccccccccccccccccccccccccccccccccccccccccccccccccccccccc
	SUBROUTINE ESCREVE1(au,eu,diu,wu,omu,amu,an,en,din,wn,omn,amn,
     >amiura,aminet,t)
      IMPLICIT REAL *8(a-h,o-z)
      COMMON/IN/ikepler,inicio,pi,pi05,pi15,pi2
      COMMON/CONST/conv

c      OPEN(1,FILE='pluexato.dat',STATUS='UNKNOWN')
c      OPEN(2,FILE='plnexato.dat',STATUS='UNKNOWN')   

      OPEN(19,FILE= 'a1e1ED.dat',STATUS='UNKNOWN')
      OPEN(29,FILE= 'a2e2ED.dat',STATUS='UNKNOWN')
      OPEN(39,FILE= 'i1i2ED.dat',STATUS='UNKNOWN')   
	OPEN(49,FILE= 'om1om2ED.dat',STATUS='UNKNOWN')
      OPEN(59,FILE= 'w1w2ED.dat',STATUS='UNKNOWN')
      OPEN(69,FILE= 'lam12ED.dat',STATUS='UNKNOWN')   
      OPEN(669,FILE= 'l1-l2ED.dat',STATUS='UNKNOWN')   
      OPEN(79,FILE= 'l1l2ED.dat',STATUS='UNKNOWN')
      OPEN(89,FILE= 'sig12ED.dat',STATUS='UNKNOWN')   
      OPEN(99,FILE= 'dpiED.dat',STATUS='UNKNOWN')
      OPEN(109,FILE='e1sig1ED.dat',STATUS='UNKNOWN')
      OPEN(119,FILE='e2sig2ED.dat',STATUS='UNKNOWN')
      OPEN(129,FILE='e1vpi1ED.dat',STATUS='UNKNOWN')
      OPEN(139,FILE='e2vpi2ED.dat',STATUS='UNKNOWN')
      OPEN(149,FILE='i1om1ED.dat',STATUS='UNKNOWN')
      OPEN(159,FILE='i2om2ED.dat',STATUS='UNKNOWN')
      OPEN(179,FILE='PERIODOED.dat',STATUS='UNKNOWN')


c      OPEN(3,FILE='testencano.dat',STATUS='UNKNOWN')

c     reducao de wu e wn para 0-2pi
c      if (lei.eq.1) go to 3
c     se lei=1 => singularidade em eu, en (ver xyzorb), 
c     e wu, wn recebem 9.d10 ficticios, e portanto nao existe reducao

C	NODO
      iomu=omu/pi2
      omu=omu-iomu*pi2
      if(omu.lt.0.0d0) omu=omu+pi2                   

      iomn=omn/pi2
      omn=omn-iomn*pi2
      if(omn.lt.0.0d0) omn=omn+pi2                   

C	PERI�IO
      wu=wu-aint(wu/pi2)*pi2
      if (wu.lt.0.0d0) wu=wu+pi2

      wn=wn-aint(wn/pi2)*pi2
      if (wn.lt.0.0d0) wn=wn+pi2

C	ANOMALIA M�IA
      amu=amu-aint(amu/pi2)*pi2
      if (amu.lt.0.0d0) amu=amu+pi2

      amn=amn-aint(amn/pi2)*pi2
      if (amn.lt.0.0d0) amn=amn+pi2

C	LONGITUDE M�IA
	dl1=amu+wu+omu
      idl1=dl1/pi2
      dl1=dl1-idl1*pi2
      if(dl1.lt.0.0d0) dl1=dl1+pi2                   

	dl2=amn+wn+omn
      idl2=dl2/pi2
      dl2=dl2-idl2*pi2
      if(dl2.lt.0.0d0) dl2=dl2+pi2                   

C	�GULO CR�ICO
	sig1=2.0d0*dl2-dl1-wu-omu
	isig1=sig1/pi2
      sig1=sig1-isig1*pi2
      if(sig1.lt.0.0d0) sig1=sig1+pi2                   

	sig2=2.0d0*dl2-dl1-wn-omn
      isig2=sig2/pi2
      sig2=sig2-isig2*pi2
      if(sig2.lt.0.0d0) sig2=sig2+pi2                   

C	DELTA PI
	dpi=(wu+omu)-(wn+omn)
      idpi=dpi/pi2
      dpi=dpi-idpi*pi2
      if(dpi.lt.0.0d0) dpi=dpi+pi2                   

	AA=365.25d0

      write(19,*)t/AA,au,eu
      write(29,*)t/AA,an,en
	write(39,*)t/AA,diu/conv,din/conv

	WRITE(49,*)t/AA,omu/conv,omn/conv
	write(59,*)t/AA,wu/conv,wn/conv
	WRITE(69,*)t/AA,dl1/conv,dl2/conv
	WRITE(669,*)t/AA,(dl1-dl2)/conv
	WRITE(79,*)t/AA,amu/conv,amn/conv
 
 	WRITE(89,*)t/AA,sig1/conv,sig2/conv   
	WRITE(99,*)t/AA,dpi/conv
      write(109,*)eu*dcos(sig1),eu*dsin(sig1),t/AA
      write(119,*)en*dcos(sig2),en*dsin(sig2),t/AA
      write(129,*)eu*dcos(wu+omu),eu*dsin(wu+omu),t/AA
      write(139,*)en*dcos(wn+omn),en*dsin(wn+omn),t/AA
      write(149,*)diu*dcos(omu),diu*dsin(omu),t/AA
      write(159,*)din*dcos(omn),din*dsin(omn),t/AA
      write(179,*)t/AA,2.0d0*pi*dsqrt( (au**3.0d0)/amiura),
     ^2.0d0*pi*dsqrt( (an**3.0d0)/aminet)

	RETURN
	END


cccccccccccccccccccccccccccccccccccccccccccccccccccc
ccccccccccccccccccccccccccccccccccccccccccccccccccccccccccc
	SUBROUTINE ESCREVE2(au,eu,diu,wu,omu,amu,an,en,din,wn,omn,amn,
     >amiura,aminet,t)
      IMPLICIT REAL *8(a-h,o-z)
      COMMON/IN/ikepler,inicio,pi,pi05,pi15,pi2
      COMMON/CONST/conv

c      OPEN(1,FILE='pluexato.dat',STATUS='UNKNOWN')
c      OPEN(2,FILE='plnexato.dat',STATUS='UNKNOWN')   

      OPEN(18,FILE= 'a1e1MT.dat',STATUS='UNKNOWN')
      OPEN(28,FILE= 'a2e2MT.dat',STATUS='UNKNOWN')
      OPEN(38,FILE= 'i1i2MT.dat',STATUS='UNKNOWN')   
	OPEN(48,FILE= 'om1om2MT.dat',STATUS='UNKNOWN')
      OPEN(58,FILE= 'w1w2MT.dat',STATUS='UNKNOWN')
      OPEN(68,FILE= 'lam12MT.dat',STATUS='UNKNOWN')   
      OPEN(668,FILE= '1-l2MT.dat',STATUS='UNKNOWN')  
      OPEN(78,FILE= 'l1l2MT.dat',STATUS='UNKNOWN')
      OPEN(88,FILE= 'sig12MT.dat',STATUS='UNKNOWN')   
      OPEN(98,FILE= 'dpiMT.dat',STATUS='UNKNOWN')
      OPEN(108,FILE='e1sig1MT.dat',STATUS='UNKNOWN')
      OPEN(118,FILE='e2sig2MT.dat',STATUS='UNKNOWN')
      OPEN(128,FILE='e1vpi1MT.dat',STATUS='UNKNOWN')
      OPEN(138,FILE='e2vpi2MT.dat',STATUS='UNKNOWN')
      OPEN(148,FILE='i1om1MT.dat',STATUS='UNKNOWN')
      OPEN(158,FILE='i2om2MT.dat',STATUS='UNKNOWN')
      OPEN(178,FILE='PERIODOMT.dat',STATUS='UNKNOWN')


c      OPEN(3,FILE='testencano.dat',STATUS='UNKNOWN')

c     reducao de wu e wn para 0-2pi
c      if (lei.eq.1) go to 3
c     se lei=1 => singularidade em eu, en (ver xyzorb), 
c     e wu, wn recebem 9.d10 ficticios, e portanto nao existe reducao

C	NODO
      iomu=omu/pi2
      omu=omu-iomu*pi2
      if(omu.lt.0.0d0) omu=omu+pi2                   

      iomn=omn/pi2
      omn=omn-iomn*pi2
      if(omn.lt.0.0d0) omn=omn+pi2                   

C	PERI�IO
      wu=wu-aint(wu/pi2)*pi2
      if (wu.lt.0.0d0) wu=wu+pi2

      wn=wn-aint(wn/pi2)*pi2
      if (wn.lt.0.0d0) wn=wn+pi2

C	ANOMALIA M�IA
      amu=amu-aint(amu/pi2)*pi2
      if (amu.lt.0.0d0) amu=amu+pi2

      amn=amn-aint(amn/pi2)*pi2
      if (amn.lt.0.0d0) amn=amn+pi2

C	LONGITUDE M�IA
	dl1=amu+wu+omu
      idl1=dl1/pi2
      dl1=dl1-idl1*pi2
      if(dl1.lt.0.0d0) dl1=dl1+pi2                   

	dl2=amn+wn+omn
      idl2=dl2/pi2
      dl2=dl2-idl2*pi2
      if(dl2.lt.0.0d0) dl2=dl2+pi2                   

C	�GULO CR�ICO
	sig1=2.0d0*dl2-dl1-wu-omu
	isig1=sig1/pi2
      sig1=sig1-isig1*pi2
      if(sig1.lt.0.0d0) sig1=sig1+pi2                   

	sig2=2.0d0*dl2-dl1-wn-omn
      isig2=sig2/pi2
      sig2=sig2-isig2*pi2
      if(sig2.lt.0.0d0) sig2=sig2+pi2                   

C	DELTA PI
	dpi=(wu+omu)-(wn+omn)
      idpi=dpi/pi2
      dpi=dpi-idpi*pi2
      if(dpi.lt.0.0d0) dpi=dpi+pi2                   

	AA=365.25d0

      write(18,*)t/AA,au,eu
      write(28,*)t/AA,an,en
	write(38,*)t/AA,diu/conv,din/conv

	WRITE(48,*)t/AA,omu/conv,omn/conv
	write(58,*)t/AA,wu/conv,wn/conv
	WRITE(68,*)t/AA,dl1/conv,dl2/conv
	WRITE(668,*)t/AA,(dl1-dl2)/conv
	WRITE(78,*)t/AA,amu/conv,amn/conv
 	WRITE(88,*)t/AA,sig1/conv,sig2/conv   
	WRITE(98,*)t/AA,dpi/conv
      write(108,*)eu*dcos(sig1),eu*dsin(sig1),t/AA
      write(118,*)en*dcos(sig2),en*dsin(sig2),t/AA
      write(128,*)eu*dcos(wu+omu),eu*dsin(wu+omu),t/AA
      write(138,*)en*dcos(wn+omn),en*dsin(wn+omn),t/AA
      write(148,*)diu*dcos(omu),diu*dsin(omu),t/AA
      write(158,*)din*dcos(omn),din*dsin(omn),t/AA
      write(178,*)t/AA,2.0d0*pi*dsqrt( (au**3.0d0)/amiura),
     ^2.0d0*pi*dsqrt( (an**3.0d0)/aminet)

	RETURN
	END


cccccccccccccccccccccccccccccccccccccccccccccccccccc

ccccccccccccccccccccccccccccccccccccccccccccccccccccccccccc
	SUBROUTINE ESCREVE3(au,eu,diu,wu,omu,amu,an,en,din,wn,omn,amn,
     >amiura,aminet,t)
      IMPLICIT REAL *8(a-h,o-z)
      COMMON/IN/ikepler,inicio,pi,pi05,pi15,pi2
      COMMON/CONST/conv

c      OPEN(1,FILE='pluexato.dat',STATUS='UNKNOWN')
c      OPEN(2,FILE='plnexato.dat',STATUS='UNKNOWN')   

      OPEN(17,FILE= 'a1e1DT.dat',STATUS='UNKNOWN')
      OPEN(27,FILE= 'a2e2DT.dat',STATUS='UNKNOWN')
      OPEN(37,FILE= 'i1i2DT.dat',STATUS='UNKNOWN')   
	OPEN(47,FILE= 'om1om2DT.dat',STATUS='UNKNOWN')
      OPEN(57,FILE= 'w1w2DT.dat',STATUS='UNKNOWN')
      OPEN(67,FILE= 'lam12DT.dat',STATUS='UNKNOWN')   
      OPEN(77,FILE= 'l1l2DT.dat',STATUS='UNKNOWN')
      OPEN(87,FILE= 'sig12DT.dat',STATUS='UNKNOWN')   
      OPEN(97,FILE= 'dpiDT.dat',STATUS='UNKNOWN')
      OPEN(107,FILE='e1sig1DT.dat',STATUS='UNKNOWN')
      OPEN(117,FILE='e2sig2DT.dat',STATUS='UNKNOWN')
      OPEN(127,FILE='e1vpi1DT.dat',STATUS='UNKNOWN')
      OPEN(137,FILE='e2vpi2DT.dat',STATUS='UNKNOWN')
      OPEN(147,FILE='i1om1DT.dat',STATUS='UNKNOWN')
      OPEN(157,FILE='i2om2DT.dat',STATUS='UNKNOWN')
      OPEN(177,FILE='PERIODODT.dat',STATUS='UNKNOWN')


c      OPEN(3,FILE='testencano.dat',STATUS='UNKNOWN')

c     reducao de wu e wn para 0-2pi
c      if (lei.eq.1) go to 3
c     se lei=1 => singularidade em eu, en (ver xyzorb), 
c     e wu, wn recebem 9.d10 ficticios, e portanto nao existe reducao

C	NODO
      iomu=omu/pi2
      omu=omu-iomu*pi2
      if(omu.lt.0.0d0) omu=omu+pi2                   

      iomn=omn/pi2
      omn=omn-iomn*pi2
      if(omn.lt.0.0d0) omn=omn+pi2                   

C	PERI�IO
      wu=wu-aint(wu/pi2)*pi2
      if (wu.lt.0.0d0) wu=wu+pi2

      wn=wn-aint(wn/pi2)*pi2
      if (wn.lt.0.0d0) wn=wn+pi2

C	ANOMALIA M�IA
      amu=amu-aint(amu/pi2)*pi2
      if (amu.lt.0.0d0) amu=amu+pi2

      amn=amn-aint(amn/pi2)*pi2
      if (amn.lt.0.0d0) amn=amn+pi2

C	LONGITUDE M�IA
	dl1=amu+wu+omu
      idl1=dl1/pi2
      dl1=dl1-idl1*pi2
      if(dl1.lt.0.0d0) dl1=dl1+pi2                   

	dl2=amn+wn+omn
      idl2=dl2/pi2
      dl2=dl2-idl2*pi2
      if(dl2.lt.0.0d0) dl2=dl2+pi2                   

C	�GULO CR�ICO
	sig1=2.0d0*dl2-dl1-wu-omu
	isig1=sig1/pi2
      sig1=sig1-isig1*pi2
      if(sig1.lt.0.0d0) sig1=sig1+pi2                   

	sig2=2.0d0*dl2-dl1-wn-omn
      isig2=sig2/pi2
      sig2=sig2-isig2*pi2
      if(sig2.lt.0.0d0) sig2=sig2+pi2                   

C	DELTA PI
	dpi=(wu+omu)-(wn+omn)
      idpi=dpi/pi2
      dpi=dpi-idpi*pi2
      if(dpi.lt.0.0d0) dpi=dpi+pi2                   

	AA=365.25d0

      write(17,*)t/AA,au,eu
      write(27,*)t/AA,an,en
	write(37,*)t/AA,diu/conv,din/conv

	WRITE(47,*)t/AA,omu/conv,omn/conv
	write(57,*)t/AA,wu/conv,wn/conv
	WRITE(67,*)t/AA,dl1/conv,dl2/conv
	WRITE(77,*)t/AA,amu/conv,amn/conv
 
 	WRITE(87,*)t/AA,sig1/conv,sig2/conv   
	WRITE(97,*)t/AA,dpi/conv
      write(107,*)eu*dcos(sig1),eu*dsin(sig1),t/AA
      write(117,*)en*dcos(sig2),en*dsin(sig2),t/AA
      write(127,*)eu*dcos(wu+omu),eu*dsin(wu+omu),t/AA
      write(137,*)en*dcos(wn+omn),en*dsin(wn+omn),t/AA
      write(147,*)diu*dcos(omu),diu*dsin(omu),t/AA
      write(157,*)din*dcos(omn),din*dsin(omn),t/AA
      write(177,*)t/AA,2.0d0*pi*dsqrt( (au**3.0d0)/amiura),
     ^2.0d0*pi*dsqrt( (an**3.0d0)/aminet)

	RETURN
	END


cccccccccccccccccccccccccccccccccccccccccccccccccccc


ccccccccccccccccccccccccccccccccccccccccccccccccccccccccccc
      SUBROUTINE gaussj(N,b,dmimi)
      implicit real *8(a-h,o-z)
      INTEGER NMAX
      dimension dmimi(6),a(6,6),b(6,1)
      PARAMETER (NMAX=50)
      INTEGER m,n,i,icol,irow,j,k,l,ll,indxc(NMAX),
     +indxr(NMAX),ipiv(NMAX)
      REAL big,dum,pivinv


C          LINEAR EQUATION SOLUTION BY GAUSS-JORDAN ELIMINATION
C          *****************NUMERICAL RECIPES******************
      
C     a(1:n,1:n) e' uma matriz de entrada de dimensao np X np
C     b(1:n,1:m) e' uma matriz de entrada de dimensao np X mp

c     Saidas:  a(1:n,1:n) e' substituida por sua inversa
c              b(1:n,1:m) e' substituida pela solucao 

      i=0
      j=0
	k=1
	do 2 i=1,n
	do 1 j=1,n

	if(i.eq.j) then
	a(i,j)=dmimi(k)
	k=k+1
	else
	a(i,j)=1.0d0
	endif

1     continue
2     continue               
      i=0
      j=0
      k=0       
cccccccccccc
      m=1
CCCCCCCCCCCC

      do j=1,n
          ipiv(j)=0
      enddo
      
      
      do i=1,n
      
          big=0.
          do j=1,n
              if(ipiv(j).ne.1) then
                  do k=1,n
                      if (ipiv(k).eq.0) then
                          if (abs(a(j,k)).ge.big) then
                              big=abs(a(j,k))
                              irow=j
                              icol=k
                          endif
                      else if (ipiv(k).gt.1) then
c                          write(*,*) 'singular matriz in gaussj'
c                          pause
                      endif
                  enddo
               endif
           enddo
     
      ipiv(icol)=ipiv(icol)+1
      if(irow.ne.icol) then
          do l=1,n
              dum=a(irow,l)
              a(irow,l)=a(icol,l)
              a(icol,l)=dum
          enddo
          do l=1,m
              dum=b(irow,l) 
              b(irow,l)=b(icol,l)
              b(icol,l)=dum
          enddo
      endif
      indxr(i)=irow
      indxc(j)=icol
      if  (a(icol,icol).eq.0) then
c                              write(*,*)'singular matriz in gaussj II'
c                              pause
	endif
      pivinv=1./a(icol,icol)
      a(icol,icol)=1.
      do l=1,n 
          a(icol,l)=a(icol,l)*pivinv
      enddo
      do l=1,m
          b(icol,l)=b(icol,l)*pivinv
      enddo
      do ll=1,n
          if(ll.ne.icol) then
              dum=a(ll,icol)
              a(ll,icol)=0.
              do l=1,n
                  a(ll,l)=a(ll,l)-a(icol,l)*dum
              enddo
              do l=1,m
                  b(ll,l)=b(ll,l)-b(icol,l)*dum
              enddo
          endif
      enddo
      
      enddo
c      do l=n,1,-1
c	write(*,*) 'caca final'
c	write(*,*) n, l, indxr(l), indxc(l)
c          if (indxr(l).ne.indxc(l)) then
c              do k=1,n
c                  dum=a(k,indxr(l))   
c                  a(k,indxr(l))=a(k,indxc(l))
c                  a(k,indxc(l))=dum   
c              enddo
c          endif
c      enddo
      return
      END    

cccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccc
      subroutine orbxyz(N,a,e,dincl,w,om,am,ami,x,f)
      implicit real *8(a-h,o-z)
      COMMON/IN/ikepler,inicio,pi,pi05,pi15,pi2
      dimension a(6),e(6),w(6),om(6),am(6),ami(6),
     |f(6),x(6,6),u(6),dincl(6)
 
c: rotina que passa de elementos orbitais para cartesianos
c: preciso tambem das derivadas no tempo das cartesianas(velocidades)
c: referencia: Fitzpatrick (em velocidades parece que nao esta ok)

c: entrada : elementos orbitais , saida : conversao para cartesianos
c: u=anom.excent. (sera achado via rotina Kepler)
c: am=anom. media
c: om=longitude do nodo (indefinido se I=0,mas tomo zero para que as formulas
c: do caso espacial fiquem validas no caso plano)
c: w=argumento do pericentro (se I=0 este w representara longitude do pericen.)
c: alat=latitude=w+f (representara wtil+f no caso I=0) (wtil=wpi)

c: obtencao de anom. excent. (u) e verdadeira (f) via Kepler

      call kepler(f,u,am,e,N)

      do 50 i=1,N
      inicio=0
      alat=w(i)+f(i)
      r=a(i)*(1.d0-e(i)*dcos(u(i)))
      rp=dsqrt(ami(i)*a(i))*e(i)*dsin(u(i))/r
      rfp=dsqrt(ami(i)*a(i)*(1.d0-e(i)*e(i)))/r
c: velocidades (x(1),x(2),x(3) no caso tridimensional
c: Nesta rotina se I = zero ou e= 0 nao ha nenhum problema, embora
c: para entrar com angulos eu entre com Wtil (no lugar de w)
c: e em lugar de om=nodo=indefinido tomo=0
      cfw=dcos(alat)
      sfw=dsin(alat)
      coom=dcos(om(i))
      siom=dsin(om(i))
      coi=dcos(dincl(i))
      somcoi=siom*coi
      comcoi=coom*coi
      sini=dsin(dincl(i))
      cosi=dcos(dincl(i))
c: no caso plano basta tomar I=dincl=0 e definir om=0.d0
c: tambem w daqui deve ser interpretado como sendo o wpi=long. pericentro
c: tambem alat deve ser pensado como alat=wpi+f

c: velocidades

      x(i,1)=rp*cfw*coom-rfp*sfw*coom-rp*sfw*somcoi-rfp*cfw*somcoi
      
      x(i,2)=rp*cfw*siom-rfp*sfw*siom+rp*sfw*comcoi+rfp*cfw*comcoi
      
      x(i,3)=rp*sfw*sini+rfp*cfw*sini

c: coordenadas

      x(i,4)=r*(cfw*coom-sfw*cosi*siom)    
     
      x(i,5)=r*(cfw*siom+sfw*cosi*coom)
     
      x(i,6)=r*sini*sfw

50    continue      

      return
      end
ccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccc      
C SUBROTINA PARA RESOLVER EQCAO DE KEPLER.  ADO AM ACHO U DEPOIS F
C FIZ TODOS OS TESTES : ESTA RODANDO REDONDO , U1 ,F (anomalias)
      SUBROUTINE KEPLER (f,u1,AM,E0,N)
      IMPLICIT REAL *8(A-H,O-Z) 
      COMMON/IN/ikepler,inicio,pi,pi05,pi15,pi2                       
      dimension f(6),fa(6),am(6),e0(6),u1(6)     

C      READ(*,*) AM,E0
C: metodo de Newton Raphson  
c:ccccccccccccccccccccccccccccccccccccccccccccccccccccccccccc 
c: este bloco eu usava p/ garantir grande precisao na solucao da eq.Kepler
c: apenas p/ converter no inicio os elementos orbitais em x,y,z (qdo tomo kepl=1)
c: voce pode deixar sempre kepl=1, e garante sempre grande precisao, mas fica
c: mais lento.Fica a seu criterio tirar, mas tem que dar o LITERA e PREC)
cccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccc
      
      DO 80 j=1,N
      
      litera=15
      prec=1.d-10
      if(ikepler.eq.1) then
                           litera=35
                           prec=1.d-12
                           end if
ccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccc                           
      U0=AM(j)+E0(j)*DSIN(AM(j))
      DO 10 I=1,litera
      U1(j)=U0-((U0-E0(j)*DSIN(U0)-
     /AM(j)))/(1.D0-E0(j)*DCOS(U0))
      TEST1=DABS(U1(j)-U0)
      DSU1=DSIN(U1(j))
      TEST2=DABS(U1(j)-E0(j)*DSU1-AM(j))
      IF(TEST1.LE.prec.AND.TEST2.LE.prec) GO TO 11
      U0=U1(j)
 10   CONTINUE
c**      WRITE(*,*)' NAO CONVERGENCIA EM KEPLER'
c**      write(*,*)'test2    excent', test2, e0,u1
 11   CONTINUE
      RCF=DCOS(U1(j))-E0(j)
c: rcf=0 ou quase se f=anom. verdad.= proximo de pi/2,pi*3/2.Se rcf<=1.d-11
c: suporemos f =pi/2 ou pi*3/2, mas se isto ocorrer ja no inicio (primeira)
c: chamada de Kepler, posso (opcionalmente) parar a integracao entrando com
c: nova condicao inicial para o Planeta (esperando que nao ocorra o mesmo)
c: (se porem o parametro INICIO e'=0,apenas considero f= pi/2 ou pi*3/2.Esta 
c: aproximacao so' tem problema se estamos usando Kepler para o Planeta,
c: pois no caso do satelite, esta transformacao e' apenas um dado de saida.
      if(dabs(rcf).gt.1.d-11) go to 20
      if(inicio.eq.1)then
      open (10,file='keple',status='unknown')
c      write(10,*)'melhor mudar AM ou E0 inic do Planeta,para que f seja
c     +melhor determinado'
                         stop 
                         end if
      if(dsu1.gt.0.d0) f(j)=pi05
      if(dsu1.lt.0.d0) f(j)=pi15
      return
  20  continue

      RSF=DSQRT(1.D0-E0(j)*E0(j))*DSU1
      FA(j)=DATAN(RSF/RCF)
      F(j)=FA(j)
      
      IF(f(j).ge.0.D0.AND.rcf.LT.0.D0) F(j)=FA(j)+PI
      IF(f(j).lt.0.D0.AND.rcf.lt.0.D0) F(j)=FA(j)+PI
      if(f(j).lt.0.d0.and.rcf.gt.0.d0) F(j)=FA(j)+PI2
      
80    continue                                        

c: apos orbxyz chamar Kepler, ele redefine inicio=0
      RETURN
      END 
      
c: SUBROTINA XYZORB(....) 

c: NELSON, NELSONNNNNNNNNNNNNNNNNNN (13/01/98)

c: Pelo que me lembro, o programa deve funcinar p/ I=0, mas eu quase nao testei pois
c: o nosso problema era espacial. Mas se seu I=0 (sempre plano) e'facilimo mexer
c: neste aqui p/ que funcione sem galho algum. QQr coisa podemos falar depois. 
c: No primeiro caso qdo I=quase zero, zero, e nao zero o programa sempre avisa
c: onde esta a indeterminacao. Idem se excent=zero ou quase. O aviso vem com
c: lei=1 (se lei=zero, esta tudo ok). Sugiro que sempre teste. 
c: Aqui ainda existem alguns PAUSE, mas voce pode tira-los (num processo REMOTO)
c: pois, o aviso de que deu problema voce sabe com o parametro LEI (que voce
c: sempre  deve imprimir , eu sugiro).
cccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccc

c      subroutine xyzorb(x,a,e,dincl,w,om,am,ami,lei,kk)
c     /wtif,f,rcosf,rsinf)      
      subroutine xyzorb(N,x,a,e,dincl,w,om,am,ami,lei,kk)
c      ,wtif,alat,f,rcosf,ww)
      implicit real *8(a-h,o-z)
      dimension x(6,6),am(6),e(6),a(6),dincl(6),w(6),om(6),ami(6),
     |f(6),u(6) 
      COMMON/IN/ikepler,inicio,pi,pi05,pi15,pi2   

c: referencias : Fitzpatrick e Brouwer  
c: entrada: r,v (x(1)....x(6)) , saida : a,e,incl,w,om,am)
c: se kk=1 : converto (acho) tudo. Se kk=0, entao so acho Excent. e Inclin. 
c: ami=G*(Mj+msat), msat=0. ou ami=Amilu=G*(Mj+dml),etc.  depende do corpo    
c: extrv=produto escalar r.v
c: lei=1 se excent. ou inclinacao sao nulos ou quase nulos e lei entao
c: pode ser usado para controlar a impressao de dados de saida.Nlei=0 se
c: nenhum problema de singularidade.

      do 70 i=1,N

      lei=0
      extrv=x(i,1)*x(i,4)+x(i,2)*x(i,5)+x(i,3)*x(i,6)
      r2=x(i,4)*x(i,4)+x(i,5)*x(i,5)+x(i,6)*x(i,6)
      v2=x(i,1)*x(i,1)+x(i,2)*x(i,2)+x(i,3)*x(i,3)
  20  r=dsqrt(r2)
      v=dsqrt(v2)
      amido=ami(i)+ami(i)
      a(i)=ami(i)/(amido/r-v2)
c: achar angulo teta entre r e v (Fitzpatrick p.71)
c: teta serve tambem p/ definir quadrante de u=anom. excentrica (ver p.71)
c: O processo de Fitzpatrick, para achar u e excent. esta' bom (abaixo),mas,
c: prefiro usar o processo de Brouwer para achar excen/ e u (Brouwer p.48)

ccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccc
c: esta parte (acha u , e ) esta ok,bom mas prefiro Brouwer
c:       costeta=extrv/(r*v)
c:       sinteta2=1.d0-costeta*costeta
c:       efitz=dsqrt(1.d0-r*v2*(2.d0-r*v2/ami)*sinteta2/ami)
c:       afitz=r/(2.d0-r*v2/ami)
c:       cosu=(a-r)/(a*e)
c:       write(*,*)'a,r,e,cosu',a,r,e,cosu
c:       pause
c:       ufitz=dacos(cosu)
c:       if(costeta.lt.0.d0) ufitz=-ufitz
c:       u=ufitz
c:       e=efitz
cccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccc
c:              abaixo , Brouwer :   (parece ser melhor)

      esinu=extrv/dsqrt(ami(i)*a(i))
      ecosu=r*v2/ami(i)-1.d0
      eb=dsqrt(esinu*esinu+ecosu*ecosu)
      e(i)=eb
c      open(85,file='xyzprob',status='unknown')        
      if(kk.eq.0) go to 30
      if(e(i).le.1.d-11. or. ecosu.eq.0.d0) then  
                                        w(i)=9.d10
                                        f(i)=9.d10
                                        am(i)=9.d10
                                        lei=1           
c                                        write(85,*)'e=0 ou quase'
c                                        write(*,*)'e=0 ou quase'
c                                        pause
                                        go to 30
                                        end if
      ub=datan(esinu/ecosu)
c: em principio,parece ser esta (abaixo) a melhor forma de definir o quadrante

c*****************************************************************************
C*******MODIFIQUEI A DEFINICAO DO QUADRANTE, FAZENDO DE 0-360 GRAUS***********
C*****************************************************************************
      
      if(ub.ge.0.d0 . and . ecosu.lt.0.d0) ub=ub+pi
      if(ub.lt.0.d0 . and . ecosu.lt.0.d0) ub=ub+pi
      if(ub.lt.0.d0 . and . ecosu.gt.0.d0) ub=ub+pi2          
c:    write(*,*)'ub,ufitz,eb,efitz',ub,ufitz,eb,efitz
      u(i)=ub
      rsinf=dsqrt(1.d0-e(i)*e(i))*dsin(u(i))
      rcosf=dcos(u(i))-e(i)
      if(dabs(rcosf).gt.1.d-11) go to 31
      if(rsinf.gt.0.d0) f(i)=pi05
      if(rsinf.lt.0.d0) f(i)=pi15
      go to 29
  31  continue
      f(i)=datan(rsinf/rcosf)
      if(f(i).ge.0.d0 . and . rcosf.lt.0.d0) f(i)=f(i)+pi
      if(f(i).lt.0.d0 . and . rcosf.lt.0.d0) f(i)=f(i)+pi
      if(f(i).lt.0.d0 . and . rcosf.gt.0.d0) f(i)=f(i)+pi2      
C      sinf=dsin(f)
C      cosf=dcos(f)
C      ff=datan(sinf/cosf)
   29 continue 



      am(i)=u(i)-e(i)*dsin(u(i))
      


  30  continue
c: determinacao da inclinacao
      h=dsqrt(ami(i)*a(i)*(1.d0-e(i)*e(i)))
      cosi=(x(i,4)*x(i,2)-x(i,5)*x(i,1))/h
      if(kk.eq.0) then
                      dincl(i)=cosi
                      return
      end if   

C     Daqui para baixo, modifiquei praticamente toda a subrotina por
c     causa do caso PLANO.      
      difcosi=cosi-1.d0
ccccccccccccccccccc
      IF (difcosi.gt.0) THEN
c     Se o cosi>>>1, entao ha' problemas serios. Testes mostraram
c     que mesmo nos casos corretos, APESAR DE I=0, difcosi pode atingir valores 
C     da ordem de 10d-16, talvez por "imprecisao" do fortran (??). No entanto, aqui, 
c     apesar de absurdo, aproximo para i=0 QUANDO:      
                          IF((dabs(difcosi)*1.0e+16).lt.9 ) then                      
                                       dincl(i)=0.0d0
                                       om(i)=0.0d0
                                       lei=5
                          ELSE                   
C     se difcosi>9.d-16, ai considero cosi>>1:                          
c                          WRITE(*,*)'problemas serios com cosi=',cosi
c                          WRITE(85,*)'problemas serios com cosi=',cosi                          
c                          pause
                          end if
      
c      write(85,*)'cosi>1',cosi,difcosi
      else                                              


                          DINCL(i)=DACOS(COSI)      



c     Para o caso plano:
C     Do mesmo modo: Se cosi<1: alguns testes mostram que qdo cosi<1, mas bem
c     proximo de 1, a diferenca  difcosi=cosi-1=-1.0e-16,-2.0e-16,... . 
c     Assim, quando:
C                  IF( (dabs(difcosi)*1.0e+16).gt.1 ) then       
c  aproximo para i=0 nestes casos tambem.      
c: pequenos erros (ao fazer I=0) nao afetam nada, pois este erro nao
c: sera carregado adiante, e' apenas saida de dado.
C                      lei=6
C                      dincl(i)=0.0d0
C                      om(i)=0.d0
C                  END IF

      end if
cccccccccccccc
      if(e(i).le.5.d-11) then
      w(i)=0.0D0
      lei=7
      return 
      end if
cccccccccccccc
      if(dincl(i).eq.0.0d0) then
c     como i=0, w agora e' a longitude do periastro
      
      if(DABS(X(i,4)).gt.1.D-11) go to 34 
c      WRITE(*,*)'PROBLEMAS COM X(i,4)=',X(i,4)
c      WRITE(85,*)'PROBLEMAS COM X(i,4)=',X(i,4)                               
      lei=8          
      w(i)=88
      return
  34  continue                                   

                      wtif=datan(x(i,5)/x(i,4))

                    if(wtif.ge.0.d0. and . x(i,4).lt.0.d0) wtif=wtif+pi
                    if(wtif.lt.0.d0. and . x(i,4).lt.0.d0) wtif=wtif+pi
                    if(wtif.lt.0.d0. and . x(i,4).gt.0.d0) wtif=wtif+pi2
                      
                      ww=wtif-f(i)               
c     A reducao abaixo e' para os casos onde wtif avanca em relacao a f (ou 
c     vice-versa.                      

C 1)    Reduzindo w aos "senos e cossenos" abaixo, posso evitar problemas nos 
c     casos onde |wtif|<|f| e -360<w<-270 (por exemplo, w=-340 graus).
c     Por exemplo, apos a 1a. volta de wtif: se w=wtif-f=10-350=-340 graus.
c     Este w=-340 e', na verdade, igual a +20, pois o wtif ja' deu uma volta.
c     Para resolver isto, basta fazer: 
                      sinw=dsin(ww)
                      cosw=dcos(ww)
                      w(i)=datan(sinw/cosw)                       
c     Ou seja, sinw=sin(-340)=sin(+20) e cosw=cos(-340)=cos(+20), e portanto
c     o w=arctan(sin20/cos20)=+20>0, como queriamos.

C  2)   Outro problema que esta reducao conserta: se wtif=352 e f=355, w=-2. Se
c     f agora avanca 1 volta �frente de wtif (p.ex, wtif vai p/ 354 e f vai
c     para +4, teremos: w=354-4=350. Com a reducao: sin(350)=-sin(10)<0,
c     cos(350)=cos(10)>0, ou seja, w=arctan(-sin10/cos10)<0=-10 graus, 
c     como queriamos.  

c  3)   Se w<0, mas agora, 0<|w|<90, nao ha' problema, o angulo continua a 
c     sendo w<0, com essa reducao.
                      

c: aqui, em que I=0 ,este w de saida representa  w+om,ie, w e'long. pericentro
c: tambem om nao esta definido (nao existe,mas tomarei zero)
      endif


c      if(x(4).gt.0.0d0.and.x(5).gt.0.0d0) iqwtif=1
c      if(x(4).lt.0.0d0.and.x(5).gt.0.0d0) iqwtif=2
c      if(x(4).lt.0.0d0.and.x(5).lt.0.0d0) iqwtif=3      
c      if(x(4).gt.0.0d0.and.x(5).lt.0.0d0) iqwtif=4
      
c      if(rcosf.gt.0.0d0.and.rsinf.gt.0.0d0) iqf=11
c      if(rcosf.lt.0.0d0.and.rsinf.gt.0.0d0) iqf=22
c      if(rcosf.lt.0.0d0.and.rsinf.lt.0.0d0) iqf=33      
c      if(rcosf.gt.0.0d0.and.rsinf.lt.0.0d0) iqf=44
          
cccccccccccccccccccccccc

c:Aqui Inclin.nao e'zero nem quase,logo,pelo menos omega=om posso achar.
      sinomi=(x(i,5)*x(i,3)-x(i,6)*x(i,2))                                   
      cosomi=-(x(i,6)*x(i,1)-x(i,4)*x(i,3))
      if(dabs(cosomi).gt.1.d-11) go to 45
      if(sinomi.gt.0.d0) om(i)=pi05
      if(sinomi.lt.0.d0) om(i)=pi15
      go to 41
  45  continue
      om(i)=datan(sinomi/cosomi)
      if(om(i).ge.0.d0. and . cosomi.lt.0.d0) om(i)=om(i)+pi
      if(om(i).lt.0.d0. and . cosomi.lt.0.d0) om(i)=om(i)+pi
      if(om(i).lt.0.d0. and . cosomi.gt.0.d0) om(i)=om(i)+pi2      
  41  continue    
      if(eb.le.5.d-11) return 
C: Aqui excent. nao e' nula, posso achar w tambem 

      coslat=x(i,4)*dcos(om(i))+x(i,5)*dsin(om(i))
      sinlat=(x(i,5)*dcos(om(i))-x(i,4)*dsin(om(i)))*dcos(dincl(i))+
     |x(i,6)*dsin(dincl(i))
      if(dabs(coslat).gt.5.d-11) go to 46
      if(sinlat.gt.0.d0) alat=pi05
      if(sinlat.lt.0.d0) alat=pi15
      go to 47
  46  continue
      alat=datan(sinlat/coslat)
      if(alat.ge.0.d0. and . coslat.lt.0.d0) alat=alat+pi
      if(alat.lt.0.d0. and . coslat.lt.0.d0) alat=alat+pi
      if(alat.lt.0.d0. and . coslat.gt.0.d0) alat=alat+pi2      
   47 continue    
      w(i)=alat-f(i)


70    continue      
c: a formula do sinlat vem do McCuskey (posso deduzir facilmente das eqcoes
c: x,y,z  que relacionam elementos orbitais, idem coslat)
      return
      end







CCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCC      
      SUBROUTINE RA15(TM,X,V,TF,LL,XL,NV,NCLASS)

C  Integrador RADAU de Everhart - Version de orden 15
C  Integra ecuaciones de segundo orden:
C    y' = F(y,t)    es NCLASS = 1
C    y" = F(y',y,t) es NCLASS = 2  
C    y" = F(y,t)    es NCLASS = -2
C  TF es t(final) - t(initial). Debe ser negativo para 
C    integracion hacia atras
C  NV es el numero de ecuaciones diferenciales simultaneas
C  LL controla el tamanio de paso, monitoreando el error en los
C    terminos en base a la tolerancia SS = 10**(-LL). Un valor
C    tipico para LL esta entre 6 y 12
C  XL es el tamanio de paso constante si LL < 0. Si no, es el
C    tamanio de paso inicial
C  TM, X y V son tiempo, coordenadas y velocidades iniciales en
C    la entrada. En la salida son devueltas como tiempo, 
C    coordenadas y velocidades finales.
      
      IMPLICIT REAL*8 (A-H,O-Z)
      DIMENSION X(60),V(60),F1(60),FJ(60),C(21),D(21),R(21),
     >     Y(60),Z(60),B(7,60),G(7,60),E(7,60),BD(7,60),H(8),
     >     W(7),U(7),NW(8)
      LOGICAL NSF,NPER,NPQ,NCL,NES
      DATA NW/0,0,1,3,6,10,15,21/
      DATA ZERO, HALF, ONE, SR/0.0D0, 0.5D0, 1.0D0, 1.4D0/
      DATA H/         0.D0, .05626256053692215D0, .18024069173689236D0,
     >.35262471711316964D0, .54715362633055538D0, .73421017721541053D0,
     >.88532094683909577D0, .97752061356128750D0/
C  Estos valores de H son los espaciadores de Gauss - Radau,
C  escalados en el intervalo [0,1]

      NPER=.FALSE.
      NSF=.FALSE.
      NCL=NCLASS.EQ.1
      NPQ=NCLASS.LT.2
      NES=LL.LT.0
C  NCLASS =  1   ==>   NCL = .TRUE.    NPQ = .TRUE.
C  NCLASS = -2   ==>   NCL = .FALSE.   NPQ = .TRUE.
C  NCLASS =  2   ==>   NCL = .FALSE.   NPQ = .FALSE.
C  NPER es .TRUE. solo en la ultima secuencia de integracion.
C  NSF es .FALSE. solo en la secuencia inicial
C  NES es .TRUE. solo si LL es negativo. Entonces el tamanio de
C  paso es XL.
      
      IF (TF.LT.ZERO) THEN
       DIR=-ONE
      ELSE
       DIR=ONE
      END IF
      
      XL=DIR*DABS(XL)
      PW=1./9.
      
      DO N=2,8
       WW=N+N*N
       IF (NCL) WW=N
       W(N-1)=ONE/WW
       WW=N
       U(N-1)=ONE/WW
      END DO
      
      DO K=1,NV
       IF (NCL) V(K)=ZERO
       DO L=1,7
	BD(L,K)=ZERO
	B(L,K)=ZERO
       END DO
      END DO
      
      W1=HALF
      IF (NCL) W1=ONE
      C(1)=-H(2)
      D(1)=H(2)
      R(1)=ONE/(H(3)-H(2))
      LA=1
      LC=1
      
      DO K=3,7
       LB=LA
       LA=LC+1
       LC=NW(K+1)
       C(LA)=-H(K)*C(LB)
       C(LC)=C(LA-1)-H(K)
       D(LA)=H(2)*D(LB)
       D(LC)=-C(LC)
       R(LA)=ONE/(H(K+1)-H(2))
       R(LC)=ONE/(H(K+1)-H(K))
       
       IF (K.NE.3) THEN
	DO L=4,K
	 LD=LA+L-3
	 LE=LB+L-4
	 C(LD)=C(LE)-H(K)*C(LE+1)
	 D(LD)=D(LE)+H(L-1)*D(LE+1)
	 R(LD)=ONE/(H(K+1)-H(L-1))
	END DO
       END IF
      
      END DO
      
      SS=10.**(-LL)
C  Las instrucciones anteriores son calculadas solo una vez
C  en una integracion para inicializar las constantes

C  A continuacion se inicializa el tamanio de paso inicial TP
      IF (NES) THEN
       TP=XL
      ELSE IF (XL.NE.ZERO) THEN
       TP=XL
      ELSE
       TP=0.1D0*DIR
      END IF
      
      IF (TP/TF.GT.HALF) TP=HALF*TF
      NCOUNT=0

C  La linea 1000 es el comienzo de la primera secuencia.
C  NS es el numero de secuencia. 
C  NF es el numero de llamados a la subrrutina FORCE. 
C  NI es el numero de iteraciones en cada secuencia
1000  NS=0
      NF=0
      NI=6
      TM=ZERO
      CALL FORCE (TM, X, V, F1)
      NF=NF+1

C  La linea 2000 es el comienzo de cada secuencia despues de
C  la primera
2000  DO K=1,NV
       G(1,K)=B(1,K)+D( 1)*B(2,K)+D(2)*B(3,K)+
     >   D(4)*B(4,K)+D( 7)*B(5,K)+D(11)*B(6,K)+D(16)*B(7,K)
       G(2,K)=             B(2,K)+D(3)*B(3,K)+
     >   D(5)*B(4,K)+D( 8)*B(5,K)+D(12)*B(6,K)+D(17)*B(7,K)
       G(3,K)=                         B(3,K)+
     >   D(6)*B(4,K)+D( 9)*B(5,K)+D(13)*B(6,K)+D(18)*B(7,K)
       G(4,K)=B(4,K)+D(10)*B(5,K)+D(14)*B(6,K)+D(19)*B(7,K)
       G(5,K)=             B(5,K)+D(15)*B(6,K)+D(20)*B(7,K)
       G(6,K)=                          B(6,K)+D(21)*B(7,K)
       G(7,K)=                                       B(7,K)
      END DO
      
      T=TP
      T2=T*T
      IF (NCL) T2=T
      TVAL=DABS(T)

C  Comienzo de las iteraciones para cada secuencia. Realiza 6
C  iteraciones en la primera secuencia (NI=6)
      DO M=1,NI
       DO J=2,8
	JD=J-1
	JDM=J-2
	S=H(J)
	Q=S
	IF (NCL) Q=ONE
	
C  Aqui calcula los predictores de posicion y velocidad para
C  cada subsecuencia        
	DO K=1,NV
	 A=W(3)*B(3,K)+S*(W(4)*B(4,K)+S*(W(5)*B(5,K)+S*(W(6)*
     >     B(6,K)+S*W(7)*B(7,K))))
	 Y(K)=X(K)+Q*(T*V(K)+T2*S*(F1(K)*W1+S*(W(1)*B(1,K)+
     >        S*(W(2)*B(2,K)+S*A))))
	 
	 IF (.NOT.NPQ) THEN
C         Esta parte solo es calculada cuando NCLASS = 2          
	  A=U(3)*B(3,K)+S*(U(4)*B(4,K)+S*(U(5)*B(5,K)+
     >      S*(U(6)*B(6,K)+S*U(7)*B(7,K))))
	  Z(K)=V(K)+S*T*(F1(K)+S*(U(1)*B(1,K)+S*(U(2)*B(2,K)+
     >      S*A)))
	 END IF
	
	END DO
	
C  Calcula la fuerza FJ al final de cada subsecuencia
	CALL FORCE(TM+S*T, Y, Z, FJ)
	NF=NF+1
	
C  Calcula los nuevos valores de G y de B para la fuerza FJ        
	DO K=1,NV
	 TEMP=G(JD,K)
	 GK=(FJ(K)-F1(K))/S
	 
	 GO TO (102,102,103,104,105,106,107,108),J
 102        G(1,K)=GK
	    GO TO 100
 103        G(2,K)=(GK-G(1,K))*R(1)
	    GO TO 100
 104        G(3,K)=((GK-G(1,K))*R(2)-G(2,K))*R(3)
	    GO TO 100
 105        G(4,K)=(((GK-G(1,K))*R(4)-G(2,K))*R(5)-G(3,K))*
     >             R(6)
	    GO TO 100
 106        G(5,K)=((((GK-G(1,K))*R(7)-G(2,K))*R(8)-G(3,K))*
     >             R(9)-G(4,K))*R(10)
	    GO TO 100
 107        G(6,K)=(((((GK-G(1,K))*R(11)-G(2,K))*R(12)-
     >             G(3,K))*R(13)-G(4,K))*R(14)-G(5,K))*R(15)
	    GO TO 100
 108        G(7,K)=((((((GK-G(1,K))*R(16)-G(2,K))*R(17)-
     >             G(3,K))*R(18)-G(4,K))*R(19)-G(5,K))*R(20)-
     >             G(6,K))*R(21)
 100        CONTINUE
	 
	 TEMP=G(JD,K)-TEMP
	 B(JD,K)=B(JD,K)+TEMP
	 
	 GO TO (200,200,203,204,205,206,207,208),J
 203        B(1,K)=B(1,K)+C(1)*TEMP
	    GO TO 200
 204        B(1,K)=B(1,K)+C(2)*TEMP
	    B(2,K)=B(2,K)+C(3)*TEMP
	    GO TO 200
 205        B(1,K)=B(1,K)+C(4)*TEMP
	    B(2,K)=B(2,K)+C(5)*TEMP
	    B(3,K)=B(3,K)+C(6)*TEMP
	    GO TO 200
 206        B(1,K)=B(1,K)+C(7)*TEMP
	    B(2,K)=B(2,K)+C(8)*TEMP
	    B(3,K)=B(3,K)+C(9)*TEMP
	    B(4,K)=B(4,K)+C(10)*TEMP
	    GO TO 200
 207        B(1,K)=B(1,K)+C(11)*TEMP
	    B(2,K)=B(2,K)+C(12)*TEMP
	    B(3,K)=B(3,K)+C(13)*TEMP
	    B(4,K)=B(4,K)+C(14)*TEMP
	    B(5,K)=B(5,K)+C(15)*TEMP
	    GO TO 200
 208        B(1,K)=B(1,K)+C(16)*TEMP
	    B(2,K)=B(2,K)+C(17)*TEMP
	    B(3,K)=B(3,K)+C(18)*TEMP
	    B(4,K)=B(4,K)+C(19)*TEMP
	    B(5,K)=B(5,K)+C(20)*TEMP
	    B(6,K)=B(6,K)+C(21)*TEMP
 200        CONTINUE
	
	END DO
       END DO
       
C  Final de la secuencia. Calculo del control (HV) del tamanio
C  de paso       
       IF (.NOT.NES.OR.M.GE.NI) THEN
	HV=ZERO
	DO K=1,NV
	 HV=DMAX1(HV,DABS(B(7,K)))
	END DO
	TVAL1=TVAL*TVAL*TVAL*TVAL*TVAL*TVAL*TVAL
	HV=HV*W(7)/TVAL1
       END IF
      
      END DO
      
C  Calcula el nuevo tamanio de paso (TP) y lo compara con HV.
C  Si es menor que lo previsto, reinicia la secuencia anterior
C  usando un paso 0.8 veces menor
      IF (.NOT.NSF) THEN
       
       IF (NES) THEN
	TP=XL
       ELSE
	TP=(SS**PW)/(HV**PW)*DIR
	IF (TP/T.LE.ONE) THEN
	 TP=.8D0*TP
	 NCOUNT=NCOUNT+1
	 IF (NCOUNT.GT.10) RETURN
	 GO TO 1000
	END IF
       END IF
       
       NSF=.TRUE.
      END IF

C  Coordenadas y velocidades y tiempo al final de la secuencia
      DO K=1,NV
       X(K)=X(K)+V(K)*T+T2*(F1(K)*W1+B(1,K)*W(1)+B(2,K)*W(2)+
     >      B(3,K)*W(3)+B(4,K)*W(4)+B(5,K)*W(5)+B(6,K)*W(6)+
     >      B(7,K)*W(7))
       IF (.NOT.NCL) THEN
	V(K)=V(K)+T*(F1(K)+B(1,K)*U(1)+B(2,K)*U(2)+B(3,K)*U(3)
     >       +B(4,K)*U(4)+B(5,K)*U(5)+B(6,K)*U(6)+B(7,K)*U(7))
       END IF
      END DO
      
      TM=TM+T
      NS=NS+1
      
C  Si termino retorna. Caso contrario controla el tamanio de
C  la siguiente secuencia y la ajusta para cubrir exactamente
C  el intervalo de total de integracion (TF)
      IF (NPER) RETURN
      CALL FORCE (TM, X, V, F1)
      NF=NF+1
      
      IF (NES) THEN
       TP=XL
      ELSE
       TP=DIR*(SS**PW)/(HV**PW)
       IF (TP/T.GT.SR) TP=T*SR
      END IF
      
      IF (DIR*(TM+TP).GE.DIR*TF-1.D-8) THEN
       TP=TF-TM
       NPER=.TRUE.
      END IF
  
C  Ahora predice los nuevos valores de B para utilizarlos en
C  la proxima secuencia      
      Q=TP/T
      
      DO K=1,NV
       
       IF (NS.NE.1) THEN
	DO J=1,7
	 BD(J,K)=B(J,K)-E(J,K)  
	END DO
       END IF
       
       E(1,K)=     Q*(B(1,K)+ 2.D0*B(2,K)+ 3.D0*B(3,K)+
     >           4.D0*B(4,K)+ 5.D0*B(5,K)+ 6.D0*B(6,K)+ 
     >           7.D0*B(7,K))
       E(2,K)=               Q**2*(B(2,K)+ 3.D0*B(3,K)+
     >           6.D0*B(4,K)+10.D0*B(5,K)+15.D0*B(6,K)+
     >          21.D0*B(7,K))
       E(3,K)=                            Q**3*(B(3,K)+
     >           4.D0*B(4,K)+10.D0*B(5,K)+20.D0*B(6,K)+
     >          35.D0*B(7,K))
       E(4,K)=  Q**4*(B(4,K)+ 5.D0*B(5,K)+15.D0*B(6,K)+
     >          35.D0*B(7,K))
       E(5,K)=               Q**5*(B(5,K)+ 6.D0*B(6,K)+
     >          21.D0*B(7,K))
       E(6,K)=                            Q**6*(B(6,K)+
     >           7.D0*B(7,K))
       E(7,K)=   Q**7*B(7,K)
       
       DO L=1,7
	B(L,K)=E(L,K)+BD(L,K)
       END DO
      
      END DO

C  A partir de la segunda secuencia, solo realiza dos 
C  iteraciones por cada secuencia
      NI=2
      GO TO 2000

      END
	


c13    FORMAT(1F8.1,1X,1F9.6,1X,1F10.6,1X,1F8.6,1X,1F12.6,
c     /1X,1F12.6,1X,1F12.6)
C14    FORMAT(1F8.1,1X,1F9.6,1X,1F10.6,1X,1F8.6,1X,1F12.6,
C     /1X,1F12.6,1X,1F10.6,1X,i1)
c     ,1X,1F10.6,1X,i1)
c      write(*,13)t,au,am0u/conv,eu,wu/conv,wu/conv,diu/conv
c      ,omu/conv,lei       
c      write(1,13)t,au,am0u/conv,eu,wu/conv,wu/conv,diu/conv
c      ,omu/conv,lei 
c      write(1,*)
c      write(*,12)t,au,am0u/conv,eu,wu/conv,diu/conv,omu/conv,lei
c      write(1,12)t,au,am0u/conv,eu,wu/conv,diu/conv,omu/conv,lei 
c      write(*,*)     
c      write(*,12)t,an,am0n/conv,en,wn/conv,din/conv,omn/conv,lei        
c      write(2,12)t,an,am0n/conv,en,wn/conv,din/conv,omn/conv,lei  
C      write(*,13)t,an,am0n/conv,en,wn/conv,wn/conv,din/conv
c      ,omn/conv,lei        
C      write(2,13)t,an,am0n/conv,en,wn/conv,wn/conv,din/conv
c      ,omn/conv,lei  
c      write(2,*)           
c      write(*,*)   

C      write(*,13)t,au,amu/conv,eu,wu/conv,wwu/conv,diu/conv
C      ,omu/conv,lei       
C      write(1,13)t,au,amu/conv,eu,wu/conv,wwu/conv,diu/conv
C      ,omu/conv,lei
c      write(*,12)t,au,amu/conv,eu,wu/conv,diu/conv,omu/conv,lei       
c      write(1,12)t,au,amu/conv,eu,wu/conv,diu/conv,omu/conv,lei 
c      write(1,*)t,wtifu/conv,fu/conv,xura(4),rcosu,alatu/conv      
c      write(1,*)t,wtifu/conv,xura(5),xura(4),xura(5)/xura(4),
c     ?fu/conv,rcosfu,rsinfu 
c      write(1,*)             
c      write(*,*)     
C      write(*,13)t,an,amn/conv,en,wn/conv,wwn/conv,din/conv
c      ,omn/conv,lei        
C      write(2,13)t,an,amn/conv,en,wn/conv,wwn/conv,din/conv
c      ,omn/conv,lei
c      write(*,12)t,an,amn/conv,en,wn/conv,din/conv,omn/conv,lei        
c      write(2,12)t,an,amn/conv,en,wn/conv,din/conv,omn/conv,lei
c      write(2,*)t,wtifn/conv,fn/conv,xnet(4),rcosn,alatn/conv       
c      write(2,*)t,wtifn/conv,xnet(5),xnet(4),xnet(5)/xnet(4),
c     ?fn/conv,rcosfn,rsinfn 
c      write(2,*)          
c      write(*,*)




ccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccc      

c      SUBROUTINE ENERGIA(H,X)
c      IMPLICIT REAL *8(a-h,o-z)
c      DIMENSION X(12)
c      COMMON/MASSAS/dm0,G,dj2,dj4


c      xr=x(1)-x(4)
c      yr=x(2)-x(5)
c      zr=x(3)-x(6) 
c      d=sqrt(xr**2+yr**2+zr**2)    
c      ru=sqrt(x(1)**2+x(2)**2+x(3)**2)
c      rn=sqrt(x(4)**2+x(5)**2+x(6)**2)                                         

      
c      U0=-cmu/ru-cmn/rn

c      pura2=x(7)**2+x(8)**2+x(9)**2
c      pnet2=x(10)**2+x(11)**2+x(12)**2            
c      T0=0.5d0*dmiu*pura2+0.5d0*dmin*pnet2

c      H0=U0+T0

c      U1=-G*dkm/d        
      
c      T1=x(7)*x(10)+x(8)*x(11)+x(9)*x(12)
      
c      H1=U1+T1      


c      H=H0+H1 
      
c      return
c      end

cccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccc

C	SUBROUTINE SISTEMA(b,ami,N)
C	implicit real *8(a-h,o-z)
C      DIMENSION a(24,24),b(24,1),ami(8)

c	write(*,*)'ordem da matriz, n='
c	read(*,*)n
	
C	m=1

C	k=1
C	do 2 i=1,n
C	do 1 j=1,n

C	if(i.eq.j) then
C	a(i,j)=ami(k)
C	k=k+1
C	else
C	a(i,j)=1.0d0
C	endif

C1     continue
C2     continue

C	call  gaussj(a,n,b,m)

c     Sistema qualquer
     
c     do 4 i=1,n
c	do 3 j=1,n
c	write(*,*)'vetor a(i,j)',i,j,a(i,j)
c	pause
c	read(*,*)a(i,j)
c3	continue
c	write(*,*)'vetor b(i,1)',i,1
c	read(*,*)b(i,1)
c4	continue
c     conferir valores de a(i,j)
c      do 4 l=1,n
c	do 3 ll=1,n
c	write(*,*)'vetor a(l,ll)',l,ll,a(l,ll)
c	pause
c3	continue
c4     continue
	
C	SAIDA
	     
c      do 5 i=1,n

c     matriz inversa
      
c	do 4 j=1,n

c	write(*,*)'i,j,a(i,j)',i,j,a(i,j)	
c	pause
c4	continue
c	write(*,*)'i,b(i,1)',i,b(i,1)
c	pause

c5	continue
C      return
C	END

cccccccccccccccccccccccccccccccccccccccccccccccc

c     reducao de wu e wn para 0-2pi
c      if (lei.eq.1) go to 3
c     se lei=1 => singularidade em eu, en (ver xyzorb), 
c     e wu, wn recebem 9.d10 ficticios, e portanto nao existe reducao
c        inwu=wu/pi2
c        wu=wu-inwu*pi2
c        if(wu.lt.0.0d0) wu=wu+pi2                   
c        inwn=wn/pi2
c        wn=wn-inwn*pi2
c        if(wn.lt.0.0d0) wn=wn+pi2                 
c3     continue          

cccccccccccccccccccccccccccccccccccccccccccccccccc

c      CALL ENERGIA(H,X)       
c      H=H*100000.0D0
c7     format(1f8.1,1x,1f20.15)      
c      write(3,*)t,H          
c      write(*,*)t,H
c      write(*,*)

C      dhu=eu*dsin(wu)
C      dku=eu*dcos(wu)
C      write(*,*)dku,dhu
C      write(4,*)dku,dhu      
C      dhn=en*dsin(wn)
C      dkn=en*dcos(wn)
C      write(*,*)dkn,dhn
C      write(5,*)dkn,dhn      
C      write(*,*) 



C      dhu=eu*dsin(wu)
C      dku=eu*dcos(wu)
C      write(*,*)dku,dhu
C      write(4,*)dku,dhu      
C      dhn=en*dsin(wn)
C      dkn=en*dcos(wn)
C      write(*,*)dkn,dhn
C      write(5,*)dkn,dhn      
C      write(*,*)


cccccccccccccccccccccccccccccccccccccccccccccccc	

c     Conferir os valores de xyz     
c      do 77 linha=1,N
c      write(*,*)'f(i)',f(linha),linha
c      do 66 lcol=1,6
c      write(*,*)'xpla(i,j),linha,coluna',xpla(linha,lcol),linha,lcol
c66    continue
c      write(*,*)                            
c      pause
c77    continue
c      pause


ccccccccccccccccccccccccccccccccccccccccccccccccc
C      CHARACTER NOME1*20  
C      CHARACTER NOME2*20      
C      CHARACTER NOME3*20      
c      CHARACTER NOME4*20  
c      CHARACTER NOME5*20      
C      WRITE(*,*) 'Arq. Plan.:'
C      READ(*,48)NOME1
C      WRITE(*,*) 'Arq. U:'
C      READ(*,48)NOME2      
C      WRITE(*,*) 'Arq. N:'
C      READ(*,48)NOME3     
C48     FORMAT(A)   
cccccccccc

c      do 7 i=1,5	
c      write(*,12)t,a(i),am(i)/conv,e(i),w(i)/conv,di(i)/conv,
c     |om(i)/conv,lei
c      write(1,12)t,a(i),am(i)/conv,e(i),w(i)/conv,di(i)/conv,
c     |om(i)/conv,lei 
c7     continue   
c      write(*,*)
c      write(1,*)
cccccccccccccccccccccccccccc


  
c: Mercurio
c      requat=2.4397d6
c      dj2=0.027d-3
C      write(*,*)'entre dj2 de Mercurio,sabendo que o de Venus=0.027d-3'
C      read(*,*) dj2
c      dm(1)=1.d0/6.023600d6  
c      ami(1)=G*(dm0+dm(1))
c      dmimi(1)=(dm0+dm(1))/(dm0*dm(1))

c      a(1)=0.387098d0
c      e(1)=0.205633d0      
c      w(1)=8.d0*conv
c      om(1)=0.d0*conv 
c      am(1)=1.d0*conv
c      di(1)=0.0d0*conv
c: atual    epsil=0.0d0 (pag. E87-Nautical Almanac 96)        
C        write(*,*)'entre epsil de Mercurio , sendo o atual=0.d0'
C        read(*,*)epsil

  
c: Venus
c      requat=6.0518d6
c:c:     dj2=0.027d-3     c: c: c:
C      write(*,*)' entre J2 de Venus c/ rotacao (hoje 0.027d-3) '
c      , sendo atualmente=0.027d-3'
C      read(*,*)dj2
c      dm(2)=1.d0/408523.5d0
c      ami(2)=G*(dm0+dm(2))
c      dmimi(2)=(dm0+dm(2))/(dm0*dm(2))      

c      a(2)=0.7233d0
c      e(2)=0.0067d0      
c      w(2)=7.d0*conv
c      om(3)=0.d0*conv 
c      am(3)=2.d0*conv
c      di(2)=0.0d0*conv
c: atual    epsil=177.36d0         
C      write(*,*)'entre epsil de Venus , sendo o atual=177.36'
C      read(*,*)epsil


c: Terra                                                                
c      requat=0.637814d7                                                 Naut
c      dj2=1.08263d-3                                                    Naut
                                                    
c: P/Terra:tomarei: do Connaissaince:  semi-eixo, massa,J2,Requat,epsil
c: para aj vamos tomar aj=1.000001d0 (Connaissance 92) em UA
c: esol: tirado Nautical Almanac 96 (~=6 casas) 
c      dm(3)=1.d0/332946.d0
c      ami(3)=G*(dm0+dm(3))
c      dmimi(3)=(dm0+dm(3))/(dm0*dm(3))                                                                         Conn
C      epsil=23.43929115d0                                                Conn
C      write(*,*)'entre epsil da Terra, sendo o atual 23.44 graus'
C      read(*,*) epsil
c      a(3)=1.000001d0     
c      e(3)=0.01672d0                                                      Naut
c      w(3)=6.d0*conv
c      om(3)=0.d0*conv        
c      am(3)=3.d0*conv
c      di(3)=0.0d0*conv


c: Marte
c      requat=3397.d3                                                    Naut
c      dj2=0.001964d0                                                    Naut
c      dm(4)=1.d0/3098710.d0
c      ami(4)=G*(dm0+dm(4))                                               Conn
c      dmimi(4)=(dm0+dm(4))/(dm0*dm(4))

                                                                      
c      a(4)=1.5236d0                                                       Naut
c      e(4)=0.0933d0                                                       Naut
c      w(4)=5.d0*conv
c      om(4)=0.d0*conv 
c      am(4)=4.d0*conv
c      di(4)=0.0d0*conv
c: atual        epsil=25.19d0                                           Naut
C       write(*,*)'entre epsil de Marte, sendo  o atual= 25.19'
C       read (*,*) epsil  
c: se tomo epsil=Incl=29.2 e hnod=w=0 , ha um grande salto em excent.    
c: abaixo dados usados por Adrian
c:    requat=3.3972d6
c:    dj2=0.001964d0
c:    ej=0.0934183
c:    dmj=1.d0/3098710.d0
c:    aj=1.52236636d0
c:    esc=ua/requat
c:    grav=grav*esc*esc*esc
c:    write(*,*)' entre epsil de Marte'
c:    read(*,*) epsil
