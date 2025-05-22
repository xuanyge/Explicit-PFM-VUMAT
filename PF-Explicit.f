C VUMAT subroutine for implementing the explicit phase field model
      
C The manuscript is "Ge,X.,Zhou,L.,Bagherifard,S.,Guagliano,M.,A simple and efficient
C implementation of explicit phase field method in ABAQUS to address complex three-dimensional 
C fracture problems,Engineering Fracture Mechanics,doi.org/10.1016/j.engfracmech.2025.111222"

C*****************************************************************

      subroutine vumat(
C Read only -
     1  nblock, ndir, nshr, nstatev, nfieldv, nprops, jInfoArray,
     2  stepTime, totalTime, dtArray, cmname, coordMp, charLength,
     3  props, density, strainInc, relSpinInc,
     4  tempOld, stretchOld, defgradOld, fieldOld,
     3  stressOld, stateOld, enerInternOld, enerInelasOld,
     6  tempNew, stretchNew, defgradNew, fieldNew,
C Write only -
     5  stressNew, stateNew, enerInternNew, enerInelasNew )
C
      include 'vaba_param.inc'
C
      dimension props(nprops), density(nblock),
     1  coordMp(nblock,*),
     2  charLength(*), dtArray(*), strainInc(nblock,ndir+nshr),
     3  relSpinInc(*), tempOld(*),
     4  stretchOld(*), defgradOld(*),
     5  fieldOld(*),  stressOld(nblock,ndir+nshr),
     6  stateOld(nblock,nstatev), enerInternOld(nblock),
     7  enerInelasOld(nblock), tempNew(*),
     8  stretchNew(*), defgradNew(*), fieldNew(*),
     9  stressNew(nblock,ndir+nshr), stateNew(nblock,nstatev),
     2  enerInternNew(nblock), enerInelasNew(nblock), jInfoArray(*),
     3  stran(nblock,ndir+nshr),eigVal(nblock,3),alpha(nblock,3)
C             
      real*8 I1,I2,I3,AS,DS,IS,CM,fai
      parameter (pi=3.141592653589)
C
      character*80 cmname
C      
      E=props(1)   ! Young's modulus
      xnu=props(2) ! Poisson's ratio
      xl=props(3)  ! Regularization length
      Gc=props(4)  ! Fracture toughness
      
      twomu=E/(1.d0+xnu)
      eg=E/(1.d0+xnu)/2.d0
      sixmu=3.d0*twomu
      alamda=twomu*(E-twomu)/(sixmu-2.d0*E)
      
      do i=1,nblock
          
C      Recover the phase field and the history field 
      phi=tempOld(i)
      psit=stateOld(i,7)
      g=(1.d0-phi)**2.d0+1.d-7
      
C      Update stresses      
      trace  = strainInc(i,1) + strainInc(i,2) + strainInc(i,3)
      stressNew(i,1)=stateOld(i,1)+alamda*trace+twomu*strainInc(i,1)
      stressNew(i,2)=stateOld(i,2)+alamda*trace+twomu*strainInc(i,2)
      stressNew(i,3)=stateOld(i,3)+alamda*trace+twomu*strainInc(i,3)
      stressNew(i,4)=stateOld(i,4)+twomu*strainInc(i,4)
      stressNew(i,5)=stateOld(i,5)+twomu*strainInc(i,5)
      stressNew(i,6)=stateOld(i,6)+twomu*strainInc(i,6)
        
C      Update strain
      stran(i,1:3)=stateOld(i,8:10)+strainInc(i,1:3)
      stran(i,4:6)=stateOld(i,11:13)+2.d0*strainInc(i,4:6)
      stateNew(i,8:13)=stran(i,1:6)
        
C      Calculate the strain invariants
      I1=stran(i,1)+stran(i,2)+stran(i,3)
      
      I2=stran(i,1)*stran(i,2)+stran(i,2)*stran(i,3)+
     1     stran(i,1)*stran(i,3)-
     2    (stran(i,4)**2+stran(i,5)**2+stran(i,6)**2)/4.d0
      
      I3=stran(i,1)*stran(i,2)*stran(i,3)
     1    -(stran(i,1)*stran(i,5)**2+stran(i,2)
     2    *stran(i,6)**2+stran(i,3)*stran(i,4)**2)/4.d0
     3    +stran(i,4)*stran(i,5)*stran(i,6)/4.d0
      
C      Calculate the principal strain

      AS=I1/3.0
      IS=I1**2-3.0*I2
      if(IS.gt.1e-16)then
         DS=sqrt(IS)
      else
         DS=0.0
      endif
      
      if(DS.gt.0)then
         CM=(2.0*I1**3-9.0*I1*I2+27.0*I3)/DS**3/2.0
         if(CM.gt.1.0)then
            CM=1.0
         elseif(CM.lt.-1.0)then
            CM=-1.0
         endif
         fai=acos(CM)/3.0
         eigVal(i,1)=AS+2.0*DS*cos(fai)/3.0
         eigVal(i,2)=AS+2.0*DS*cos(fai+2.0*pi/3.0)/3.0
         eigVal(i,3)=AS+2.0*DS*cos(fai+4.0*pi/3.0)/3.0
      else
         eigVal(i,1)=AS
         eigVal(i,2)=AS
         eigVal(i,3)=AS
      endif

C      Spectral tension¨Ccompression split
      
      trp1=(eigVal(i,1)+eigVal(i,2)+eigVal(i,3)+
     1 abs(eigVal(i,1)+eigVal(i,2)+eigVal(i,3)))/2.d0     
       
      trp2=0.d0
      do j=1,3
      trp2=trp2+(eigVal(i,j)+abs(eigVal(i,j)))**2.d0/4.d0
      end do      
      psip=xnu*eg/(1d0-2d0*xnu)*trp1**2d0+eg*trp2

      H=max(psit,psip)
      
C      Update the stress and the history field
      
      stateNew(i,1:6)=stressNew(i,:)
      stressNew(i,:)=stressNew(i,:)*g        
      stateNew(i,7)=H
      
C      Update the inelastic dissipated energy and restrict the phase field rate
      
      if(tempNew(i) .le. tempOld(i)) tempNew(i)=tempOld(i)
        
      if(tempNew(i) .ge. 1.d0) then
      enerInelasNew(i)=-1.d0/xl**2*
     1  dtArray(1)/density(i)+enerInelasOld(i)
      tempNew(i)=1.d0
      else
        enerInelasNew(i)=-(phi/xl**2-2.d0*(1.d0-phi)*H/(Gc*xl))*
     1  dtArray(1)/density(i)+enerInelasOld(i)      
      endif

C       Update the specific internal energy
      stressPower = 0.5 * (
     1    ( stressOld(i,1)+stressNew(i,1) )*strainInc(i,1)
     1    +     ( stressOld(i,2)+stressNew(i,2) )*strainInc(i,2)
     1    +     ( stressOld(i,3)+stressNew(i,3) )*strainInc(i,3)
     1    + 2.d0 *( stressOld(i,4)+stressNew(i,4) )*strainInc(i,4)
     1    + 2.d0 *( stressOld(i,5)+stressNew(i,5) )*strainInc(i,5)
     1    + 2.d0 *( stressOld(i,6)+stressNew(i,6) )*strainInc(i,6) )
C
      enerInternNew(i) = enerInternOld(i) + stressPower / density(i)  
        
      enddo
C
      return
      end
