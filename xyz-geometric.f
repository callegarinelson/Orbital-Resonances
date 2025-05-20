c:  27-09-2019 :    o correto e'tomar um Grav*Msaturno como UMA so constante, ao inves de 2 constantes separadas

c: 08-set-2019: Atencao: as unidades daqui tem que serem as mesmas do integrador P-SICARDY-08-09-19.for (UA e dias)

c: 29-ago-2019: Hoje usarei Gmp do R-Sicardy (ao inves de Grav, dm0 separadamente). Unidades ficam ainda: UA e dias
              
c---------------------------------------------------------------

c: 23-agosto de 2019: retomei a versao do Giupone de 2018 e agora calculo a=semieixo via Momento angular Hz- apenas isso
c: 23-agosto-2019 :   o novo semi-eixo eh calculado via function Amelhor(....)

c: 12-fev-2018: ajeito o anterior (ELEM-SICARDY-15-10-2017.for) para mandar para Cristian Giupone

c: OUTUBRO-dia 15: O mais novo-IGNORAR os outros anteriores 
c: 14-outubro-2017- contem ate J6 (isto é, apenas p/ SATURNO
c:                  Atencao: o calc. de a=semi-eixo posso fazer pela eqc 42. formalismo usual.
c:                  DESISTI de mexer com aquilo mencionado na pag.245 

c: 11-outubro de 2017: trouxe o ELEM-SICARDY.for do SICARDY-GEOMETRICO-adapto para Sat. de SATURNO-Anthe

c: 22-Setembro-2017: faco SIMPLESMENTE o que esta escrito no Renner-Sicardy: ver pag.244:
c: dado (a,e,I, lambda,wtil,nodo)0, acho (rc,Lc,zc,rptc,Lptc,zptc)0 e frequencias (n,kapa,ni,eta,qui)0,
c: obter  novo (a,e,I, lambda,wtil,nodo)1, obter (rc,Lc,zc,rptc,Lptc,zptc)1 e frequencias (n,kapa,ni,eta,qui)1
c: obter  novo (a,e,I, lambda,wtil,nodo)2, obter (rc,Lc,zc,rptc,Lptc,zptc)2 e frequencias (n,kapa,ni,eta,qui)2
c: NOTE CUIDADOSAMENTE que: (rc,Lc,zc,rptc,Lptc,zptc)s e frequencias (n,kapa,ni,eta,qui)s sao obtidas usando
c:                          (a,e,I, lambda,wtil,nodo)s, de indice s - isso e o que esta no Renner Sicardy,pag 244.
c: NOTE: (r,L,z,rpt,Lpt,zpt) nao sao atualizados, sempre os mesmos.


      implicit none
      real*8 grav, ua, Mp,Mp0,Rp,Gmp, amelhor,r0,b0, BBmelhor,testeBB
      real *8 x(6),t,r,L,rpt,Lpt,rc,Lc,z,zc,rptc,Lptc,zptc,n,ni,a0,a,
     +      a1,e,I, kapa,kapa2,eta2,qui2,c1,c2,J2, alfa1,alfa2,alfaqu,
     +      zpt, J4,adn,adkapa,adni,adeta2,adqui2
      real*8 dt1,dt2,dt3,dt4,dt5,dt6,dt7,dt8,dt9,dt10,dt11,dt12
      real*8 lambda,wtil,om,w,am
      real*8 pi, pi2, conv, ang(5), tmax,parar, sinx,cosx, arco
      real*8 J6 , arcoteste  ! 11-10-2017 para Anthe
      real*8 Hz , r0c,rR2,n02  ! 22-08-2019 para achar semi-eixo (elem. geom.) no fim do programa
      common/barra1/Rp,Gmp,J2,J4,J6
      integer icont, it,itmax, iang(5),mir,k  

      data itmax/100/    ! em geral itmax=50 deve ser suficiente

c: Alguns dados básicos de entrada:
c: Criar um arquivo de sete colunas x(1),x(2),x(3),x(4),x(5),x(6),t= [x,y,z,vx,vy,vz,t]-coord cartesianas e t
c: Por exemplo,aqui : ANTHE-YAHOO-LIMPO-2017.dat", status="old"
c: Criar outro arq. para salvar os elementos geometricos relativos ao arq acima.
c: Por exemplo aqui: open (55,file='GEO3-ANTHE-YAHOO-2017.dat', status='unknown')

   
c: O calculo dos elementos geometricos eh iterativo
     
c             t=1.0d0
c             x(1)=0.0012620112193601d0
c             x(2)=-0.0003904786169965d0
c             x(3)=4.1541288959938d-007
c             x(4)=0.0023612513802838d0
c             x(5)=0.0076557505032543d0
c             x(6)=1.5494141103581d-006

             
             
      pi=dacos(-1.d0)
      pi2=pi+pi
      
c:      grav= 0.2958301934D-3   ! em unidades: UA**3  Msol**(-1) Dia**(-2) antigo era: 0.2959129080D-3

      ua=1.4959787d11   ! metros
      conv=pi/180.d0
c: o x deve ser dado em UA. E vx em UA/dia, portanto verificar Grav, Mp se estao corretos
5     continue

      mir=6
cccccccccccccccccccccccccc      write(*,*)'Entre mir=6 (Sat SATURNO) e mir= 0 (sat Urano) MIR=?'
ccccccccccccccc cccccccccccccccc     read(*,*) mir

      if(mir.ne.6) then   ! Aqui sats de Urano
                   write(*,*)'Desativei sat. de Urano, redefinir mir' 
                   goto 5     

                          else  

c:      Mp0=5.68319d26 ! massa de Saturno (kg), (dsol=1.988544d30 kg Horizons) !29-08-2019 -agora nao preciso Mp nem Mp0
c:      Mp=Mp0/1.988544d30  ! massa de Saturno em unidades de massa Solar ! nao preciso, usarei Gmp do Cooper ou Sicardy
cc*      Rp=60330.d3/ua      ! Raio Equatorial de Saturno  em Ua -  preciso dele sim mais

c:      Rp=60330.d3 /ua  ! atencao em UA  pois usarei Gmp do R-Sicardy ~ Cooper (sao =~s)   convertido em UA e Dias
      Rp=60268.d3/ua  ! Horizons em 28/09/2019 ! Jacobson idem (lembrar: ua eu dei em metros, por isso o D+03 no numerador)
      J2=0.01629071d0  ! Jacobson 2006=Horizons 04-outubro-2019
      J4=-0.00093583d0  ! Jacobson 2006 idem= Horizons
      J6=0.00008614d0  !   ! Jacobson 2006 idem= Horizons 04-Out-2019


c:      tmax=10.05d0  ! Anthe durante 10 dias apenas  - atualmente nao estou usando tmax

                      end if
c:*****************************************************************************************************************************************
c:     GMp=3.7931273d7  ! este eh  G*dmSaturno do Cooper (~Sicardy)- entao nao usarei Grav e dmSaturno separadamente
      
c:     GMp=3.7931273d7*(24.d0*3600.d0)**2/(1.4959787d8 **3) ! (converto para Ua, Dias) eh o GMp integral do Cooper
       GMp=3.79312078d7*(24.d0*3600.d0)**2/(1.4959787d8 **3)! do Horizons no dia 28/09/2019 ! = Jacobson-2006 -OK
c: ****************************************************************************************
c****************************************************************************************
c: Estes valores abaixo, usando Grav e  Massa de Saturno, separadamente, isto eh, considerando 2 constantes, produzem resultados  pobres,
c: corrompidos ( surge uma oscilacao falsa nos elementos geometricos

c:       Mp0=5.6834d26 ! massa de Saturno do Horizons de hoje: 27/09/2019
c:       Mp=Mp0/1.988544d30 
c:       GMp=Grav*Mp   ! usando o grav antigo junto com Mp de hoje    !!  NAO TOMAR  este GMp assim calculado. Usar o GMp observado
c*****************************************************************************************************************************
c*****************************************************************************************************************************
     
c:      open(10,file="SS1-08-09-19.dat", status="old")  ! aqui hoje 
c:      open(55,file='EG1-08-09-19.dat', status='unknown') ! aqui-hoje  
c-----------------------------------------------------------------------
c: 15-Set-2019- processar Anthe e MIMAS para gerar o Phi10 do Nelson

c:      open(10,file="ANTHE-280919-LIMPO.dat", status="old")  ! aqui hoje  29-09-2019   - ultima forma  em ~23-09-2019
c:      open(10,file="MIMAS-280919-LIMPO.dat", status="old")  ! aqui hoje  29-09-2019  - faco os dois: Anthe e Mimas hoje
c:04-10-2019---------------------------------------------------------------------------------------------------------
c:      open(10,file="anthe-300919-LIMPO.dat", status="old")  !  hoje  30-09-2019  - faco : anthe - 3 dias apenas passo 5 minutos



      
c       open(10,file="SS1-mimas.dat", status="unknown") !  hje  04-10-2019- vem de freq-nove4-2.-otubro.for: MIMAS- 152d+3 dias saida 3dias

       open(10,file="SS1.dat", status="unknown")






c:04-10-2019---------------------------------------------------------------------------------------------------------
c:      open(55,file='ANTHE-250919-EG-B.dat', status='unknown') ! hoje 250919 tomo Gmp do Cooper

c:      open(55,file='ANTHE-250919-EG-C.dat', status='unknown') !27-09-2019 : tomo o GM do HORIZON de hoje
c:      open(55,file='ANTHE-250919-EG-D.dat', status='unknown') !27-09-2019 : tomo o GM do HORIZON de hoje e desativo Amelhor
c:
c:      open(10,file="ANTHE2-CURTA.dat", status="old")  ! aqui hoje  23-09-2019 saida 0.5 anos so para obter o semi-eixo em KM
c:      open(55,file='ANTHE2-CURTA-NELSON.dat', status='unknown') ! hoje 23-09-19  saida 0.5anos e obter o semi-eixo em Km

cccc:      open(55,file='ANTHE-250919-EG-E.dat', status='unknown') !27-09-2019 : tomo o Grav antigo e desativo Amelhor
c: o ultimo que usei hoje (27-09-2019)
c:&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&
          
c:        open(55,file='ANTHE-290919-EG.dat', status='unknown') !29-09-2019 : DAR aqui o arquivo que recebera os xyz reduzidos em ELEM. GEO.
c:        open(55,file='MIMAS-290919-EG.dat', status='unknown') !29-09-2019 Faco Anthe e Mimas hoje 29-09-2019
c: ********************************************************************************************************************
c:        open(55,file='anthe-300919-EG.dat', status='new') !30-09-2019 : Faco Anthe -so 3 dias, passo =5minutos
       
       
       
       
       
       
       
       
       
c        open(55,file='EG-Mimas.dat', status='unknown') !04-10-2019 : Faco MIMAS2016-D.dat- 152mil dias saida 3 dias

       open(55,file='EG.dat', status='unknown')








c:**********************************************************************************************************************         
c:           open(55,file='SS1-2018GGG.dat', status='unknown') ! p/guardar elem geomet (RS, rs de RSicardy)
cc*           open(55,file='EG1-AGO-2019-PV.dat', status='unknown') ! guardo elem geomet (RS, rs de RSicardy)-versao de 2019
            

20    continue ! Ler  [x(1)...x(6)] = (x,y,z, vx,vy,vz)= coordenadas e velocidades do satelite de Saturno

c:**********************************************************************************************************************************

c: CUIDADO: em alguns read abaixo, precisei ler o tempo apos o x(6)-- CUIDADO - CUIDADO - veja 

c: 04-OUTUBRO-2019: quando xyz vem do freq-novej4-2.for-meu integrador, entao o TEMPO t vem no fim , apos o x(6)
c:                  dai usar o read (10,*...) adequado     
      
c:      read(10,*,err=100, end=100) t, x(1),x(2),x(3),x(4),x(5),x(6) ! x(k)=x,y,z,vx,vy,vz (k=1,2...6) e t=tempo, vem na col UM
      read(10,*,err=100, end=100) x(1),x(2),x(3),x(4),x(5),x(6),t
c       x(k)=x,y,z,vx,vy,vz (k=1,2...6) e t=tempo, vem na col SETE
     
c: Atencao: aqui, lendo Anthe, satelite de Saturno: os x, vx tem que vir em
c:  UA, UA/dia. Verificar: Gmp deve estar em UA , dia e massa solar 

c:**********************************************************************************************************************************

c-----------------------------------------------------------------------
c: Comeca aqui o processo do Renner-Sicarddy

c: eqcoes 22-25 do Renner-Sicardy
      z=x(3)
      zpt=x(6)   ! velocidade do z=x(6)
      r=dsqrt(x(1)*x(1)+x(2)*x(2))
      hz=x(1)*x(5)-x(2)*x(4)  ! componente vertical do momento angular do satelite (supondo satelite isolado, so 2 corpos)

c------------------------------------------------
      cosx=x(1)   ! na verdade cosx eh imaginado como cosx= r cos(L) e sinx= r sin(L)y e  obter L tal que x=rcos(L), y=rsin(L) 
      sinx=x(2)
      if(cosx.eq.0.d0 . and. sinx.gt.0.d0) then
                                         L=pi/2.d0
                                         go to 70
                                         end if
      if(cosx.eq.0.d0 . and. sinx.lt.0.d0) then
                                         L=pi*1.5d0
                                         go to 70
                                         end if

      arco=datan(sinx/cosx)
      arcoteste=datan2(sinx,cosx)  ! testando datan2

      if(sinx.ge.0.d0 .and .cosx.gt.0.d0) then 
                                           L=arco
      elseif(sinx.ge.0.d0. and. cosx.lt.0.d0) then 
                                           L=arco+pi
      elseif(sinx.le.0.d0. and. cosx.lt.0.d0) then 
                                           L=arco+pi
      elseif(sinx.le.0.d0. and.cosx.gt.0.d0) then 
                                           L=arco+pi2
      end if
c----------------------------------------------------   
70    continue
      rpt=x(4)*dcos(L)+x(5)*dsin(L)   !rpt=velocidade do escalar r, x(4)=vx e x(5)=vy
      Lpt=(-x(4)*dsin(L)+x(5)*dcos(L))/r  !velocidade de L= Lpt
      
c:      print*,'x(4),x(5), L', x(4),x(5),L/conv, arcoteste/conv
c:      pause
c-------------------------------------------------------------------------
c:   forneco valores iniciais para comecar as iteracoes:
      a=r       ! inicial, so uso uma vez   (*)
c      a=150000.497d0/ua  ! inicio com o a sugerido na pagina 246

      b0=r      ! uso para entrada da function AMELHORA
      e=0.d0    ! inicial, so uso uma vez   (*)
      I=0.d0    ! inicial, so uso uma vez   (*)
      rc=0.d0   !              idem         (*)
      Lc=0.d0   !              idem         (*)
      zc=0.d0   !              idem         (*)
      rptc=0.d0               !idem         (*)
      Lptc=0.d0               !idem         (*)
      zptc=0.d0               !idem         (*)

      a0=a  ! salvar para controlar a parada das iteracoes (ver dabs(a0-a)) ! no momento nao estou usando. 
      it=1  ! contar o numero de iteracoes para achar o elemento geometrico

c:   Agora posso comecar as iteracoes. Como dito na pag.244 do Renner, neste ponto, tomando
c:   eqs 22-25  e os a,e,..Lc,..Zptc dados acima, obtenho de 14-21, as frequencias (f1). Com estas (f1)
c:   e usando os dados acima (a,e,I,rc...zptc) resolvo 42-47, obtendo
c:   o primeiro conjuto de elementos geomet (elge1). Com estes novos elge1, acho também novos 
c:   rc,..rptc .. zptc (correc1), usando eqs 36-41 (os iniciais acima (*) nao usarei mais.
c 
c:   Entao: primeiro ciclo completo. Preparo o 2-ciclo: com elge1, acho novas frequencias (f2) e
c:   usando (correc1), uso 42-47, obtenho elge2. Novamente com este elge2, acho correc2 atraves de 36-41.
c:   Assim, fechei ciclo2, i.e, tenho elge2 e correc2. Posso repetir o ciclo 3 , de forma similar.
  

30    continue

cccccccccccccccccccccccccccccc      write(*,*)'oi'
c: eqcoes 14-21 - frequencias - estas eh que entram nas correcoes 36-41
      
      c1=dsqrt(GMp/a**3)
      c2=(Rp/a)**2
      
      n=c1*(1.0d0+0.75d0*c2*J2-9.d0/32.d0*c2**2*J2**2+ 
     +      27.d0/128.d0*c2**3*J2**3+3.d0*c2*J2*e*e-
     +    12.d0*c2*J2*I*I)
c: adicionando termos de ordem J4,J6 em n
      adn=-15.d0*c2*c2*J4/16.d0+45.d0*c2**3*J2*J4/64.d0+
     +     35.d0*c2**3*J6/32.d0
      n=n+c1*adn
      
c------------------------------------------------------------------------
      kapa=c1*(1.0d0-0.75d0*c2*J2-9.d0/32.d0*c2*c2*J2*J2-
     +     27.d0/128.d0*c2**3*J2**3-9.d0*c2*J2*I*I) 
c: adicionando termos de ordem J4 em kapa
      adkapa=45.d0*c2*c2*J4/16.d0+135.d0*c2**3*J2*J4/64.d0-
     +       175.d0*c2**3*J6/32.d0
      kapa=kapa+c1*adkapa
c------------------------------------------------------------------------
      ni=c1*(1.0+2.25d0*c2*J2-81.d0/32.d0*c2*c2*J2*J2+
     +     729.d0/128.d0*c2**3*J2**3+6.d0*c2*J2*e*e-
     +     51.d0/4.d0*c2*J2*I*I)
c: adicionando termos de ordem J4 em ni
      adni=-75.d0*c2*c2*J4/16.d0+675.d0*c2**3*J2*J4/64.d0+
     +      245.d0*c2**3*J6/32.d0
      ni=ni+c1*adni
c----------------------------------------------------------------
      eta2=c1*c1*(1.d0-2.d0*c2*J2)
c: adicionando termos de ordem J4 em eta2
      adeta2=75.d0*c2*c2*J4/8.d0-175.d0*c2**3*J6/8.d0
      eta2=eta2+c1*c1*adeta2
c---------------------------------------------
      qui2=c1*c1*(1.d0+7.5d0*c2*J2)
c: adicionando termos de ordem J4 em qui2
      adqui2=-175.d0*c2*c2*J4/8.d0+735.d0*c2**3*J6/16.d0
      qui2=qui2+c1*c1*adqui2
c------------------------------------------------------
      alfa1=0.333333333333333333d0*(2.d0 *ni+kapa)
      alfa2=ni+ni-kapa
      alfaqu=alfa1*alfa2   ! cuidado que eh alfaqu (sem o a final)
      kapa2=kapa*kapa   ! uso na pg. 243 do Renner
c------------------------------------------------------

c--------------------------------------------------------------------------
c: eqcoes 42-47 : obtendo os elementos geometricos


      dt1=r-rc
      dt2=(Lpt-Lptc-n)/(n+n)
       

      a=dt1/(1.d0-dt2)   ! este eh o formalismo da pg. 244 para achar a. 
     
c------------------------------------------
c: 17-fev-2018 : nao vou usar procedimento da Hz, pag.245 - Parece INUTIL- DESISTO

      dt3=(rpt-rptc)/(a*kapa)
      e=dsqrt(dt2**2+dt3**2)

      dt4=(z-zc)/a
      dt5=(zpt-zptc)/(a*ni)

      I=dsqrt(dt4**2+dt5**2)

      lambda=L-Lc-2.d0*n*dt3/kapa

      dt6=a*kapa*(1.d0-(r-rc)/a)
      dt7=rpt-rptc
c:--------------------------------------------------------------------
      sinx=dt7
      cosx=dt6
      if(cosx.eq.0.d0 . and. sinx.gt.0.d0) then
                                         dt8=pi/2.d0
                                         go to 200
                                         end if
      if(cosx.eq.0.d0 . and. sinx.lt.0.d0) then
                                         dt8=pi*1.5d0
                                         go to 200
                                         end if

      arco=datan(sinx/cosx)   ! tan(lambda-wtil)=sinx/cosx   e dt8=lambda-wtil   - eq.46
      if(sinx.ge.0.d0 .and .cosx.gt.0.d0) then 
                                           dt8=arco
      elseif(sinx.ge.0.d0. and. cosx.lt.0.d0) then 
                                           dt8=arco+pi
      elseif(sinx.le.0.d0. and. cosx.lt.0.d0) then 
                                           dt8=arco+pi
      elseif(sinx.le.0.d0. and.cosx.gt.0.d0) then 
                                           dt8=arco+pi+pi
      end if
c-------------------------------------------------------
200   continue
      wtil=lambda-dt8
      dt10=ni*(z-zc)
      dt9=zpt-zptc
c---------------------------------------------------------------------
      sinx=dt10
      cosx=dt9   ! tan(lambda-nodo)=sinx/cosx
      if(cosx.eq.0.d0 . and. sinx.gt.0.d0) then
                                         dt11=pi/2.d0
                                         go to 300
                                         end if
      if(cosx.eq.0.d0 . and. sinx.lt.0.d0) then
                                         dt11=pi*1.5d0
                                         go to 300
                                         end if

      arco=datan(sinx/cosx)   !dt11=lambda-nodo
      if(sinx.ge.0.d0 .and .cosx.gt.0.d0) then 
                                           dt11=arco
      elseif(sinx.ge.0.d0. and. cosx.lt.0.d0) then 
                                           dt11=arco+pi
      elseif(sinx.le.0.d0. and. cosx.lt.0.d0) then 
                                           dt11=arco+pi
      elseif(sinx.le.0.d0. and.cosx.gt.0.d0) then 
                                           dt11=arco+pi2
      end if
c-------------------------------------------------------------------
300   continue      
      om=lambda-dt11
         
      ang(1)=lambda
      ang(2)=wtil
      ang(3)=om
      ang(4)=ang(2)-ang(3)    ! ang(4) sera o w : pericentro
      ang(5)=lambda-wtil      ! ang(5) sera o am
c------------------------------------------
c: reduzir para 0-pi2
      do icont=1,5
      iang(icont)=ang(icont)/pi2
      ang(icont)=ang(icont)-iang(icont)*pi2
      if(ang(icont).lt.0.d0) ang(icont)=ang(icont)+pi2
      end do  

      lambda=ang(1)
      wtil=ang(2)
      om=ang(3)  
      w=ang(4)
      am=ang(5)        
c-------------------------------------------------------------------------------
      parar=dabs(a-a0)
c      write(*,*)parar
c      pause
      a0=a
      it=it+1
      if(parar.le.1.2d-12)         then
c:      if(it.gt. itmax)               then  ! qualquer dos criterios parece ok (itmax e parar)

c:    salvar os elem geom. achados, pelas eqcoes 42-47, mas posso dar mais uma atualizada no semi-eixo a,
c:   pois tenho novos e,I, entao posso chamar AMELHOR(....)
c:      r0=b0  ! este b0 e'o primeiro r (=a) inicial qdo foi lido o xyz,vxvyvz
c:      Aqui,teoricamente, este a=semi-eixo, via AMELHOR (...) seria melhor do que o a de 42  *** levar este bloco la em cima,eq 42

c: No entanto, testes mostram, este AMELHOR da o mesmo resultado de antes, nao mudou nada.

c:      a=AMELHOR(b0,x(1),x(4),x(2),x(5),e,I) ! ver function AMELHOR, fim do arquivo ! Desisti (Fev-2018)  , veja abaixo, RETOMEI

c:      a=Amelhor(r,Hz,Rp,Gmp,J2,J4,J6,e,I) ! desativei o calculo de a via momento angular

c: Se nao quiser usar o AMELHOR, basta COMENTAR a linha acima, como mostrado acima
cc*      a1=a*ua/1000.d0  ! transformei em Km como fez R-Sicardy
cc*      a1=a1-1.49400365d5

      a1=a  *1.4959787d8  ! apenas para saida transformo a em Km  

      write(55,*)t,lambda/conv, wtil/conv,om/conv,a1,e,I/conv,w/conv,
     +            am/conv    ! gravando elementos geometricos  no arq  55
60    format(e18.8, 3f13.4,e18.8,1x,e16.6,1x,f16.6, 2f12.3)   ! 

61    format(e18.8, 3f13.4,e17.8,1x,e14.5,1x,f14.3, 2f13.3)   !  dependendo precisao de alguma variavel usar 60/61

c-------------------------------------------------------------------
      go to 20  ! vou ler novos xyz vx vy vz

c                                                 stop
                                     end if  
c-----------------------------------------------------------------------------------

c------------------------------------------------------------------------------------
c: achar novas correcoes (rc,Lc,..) para fazer a proxima iteracao
c: eqcoes 36-41 - correcoes do tipo : rc, Lc, *c

84    continue        
      rc=a*e*e*(1.5d0*eta2/kapa2-1.d0-0.5d0*eta2*dcos(2.d0*lambda-
     +   2.d0*wtil)/kapa2)+ 
     +  a*I*I*(0.75d0*qui2/kapa2-1.d0+0.25d0*qui2*dcos(2.d0*lambda-
     +   2.d0*om)/alfaqu) 
      
      Lc=e*e*(0.75d0+0.5*eta2/kapa2)*n/kapa*dsin(2.d0*lambda-2.d0*wtil)
     +   -I*I*0.25d0*qui2*n*dsin(2.d0*lambda-2.d0*om)/(alfaqu*ni)

      zc=a*I*e*(0.5d0*qui2*dsin(2.d0*lambda-wtil-om)/(kapa*alfa1)-
     +   1.5d0*qui2*dsin(wtil-om)/(kapa*alfa2)) 

      rptc=a*e*e*eta2*dsin(2.d0*lambda-2.d0*wtil)/kapa-   !!!!!!!!! eh kapa**1, correto
     +   a*I*I*0.5d0*qui2*ni*dsin(2.d0*lambda-2.d0*om)/alfaqu

      Lptc=e*e*n*(3.5d0-3.d0*eta2/kapa2-0.5d0*kapa2/(n*n)+
     +     (1.5d0+eta2/kapa2)*dcos(2.d0*lambda-2.d0*wtil))+
     +     I*I*n*(2.d0-0.5d0*kapa2/(n*n)-1.5d0*qui2/kapa2-
     +     0.5d0*qui2/alfaqu*dcos(2.d0*lambda-2.d0*om))

      zptc=a*I*e*(0.5d0*qui2*(kapa+ni)*dcos(2.d0*lambda-wtil-om)/
     +     (kapa*alfa1)+1.5d0*qui2*(kapa-ni)*dcos(wtil-om)/(kapa*alfa2))

      go to 30  ! continuo iterando as frequencias e paro quando atinjo  dabs(a-a0)< 1.d-8,etc, ou it=itmax
      
100   continue
      write(*,*) " fim, tempo final=",t
      end
c-----------------------------------------------------------

c----------------------------------------------------------------------------
c: 14-10-2017 : Monto aqui a subrotina que acha o r0 adequado para achar o semi-eixo obedecendo a integral das areas.
c: Ver pg.245 do Renner-Sicardy

      function Amelhor(r,Hz,Rp,Gmp,J2,J4,J6,e,I)
      implicit none
      real*8 r,Hz,Gmp,Rp,J2,J4,J6
      real*8 r0,r0c,rR2,n0,n02,a,e,I, Amelhor
c: BLA, BLA, BLA,   desisto.
c: Como que Hz se conserva quando existem N- satelites ?
c: Obter c3  ???
c:      return
c:      end 
c-----------------------------------------------------------
c 21-agosto -2019 
c: Vou hoje implementar o AMELHOR  baseado no ERNESTO
      r0=r
      r0c=0.d0  ! isso eh imposto so para comecar as iteracoes que calcularao o r0
700   continue
      if(dabs(r0-r0c).ge.1.2d-15) then
                                 r0c=r0
                                 rR2=(Rp/r0)**2
                                 n02=Gmp/r0**3*(1.d0+1.5d0*rR2*J2-
     +                           15.d0/8.d0*rR2**2*J4+
     +                           35.d0/16.d0*rR2**3*J6)  ! esta eh a formula 49 do Renner-Sicardy
                                 n0=dsqrt(n02)
                                 r0=dsqrt(Hz/n0)  ! hz eh a componente z do momento angular do satelite (supondo que seja constante ?)
                                 go to 700
                                 end if
                                 Amelhor= r0*(1.d0+e*e+I*I)
      return
      end







