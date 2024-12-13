
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
      character*80 cmname
C      
      E=props(1) ! Young's modulus
      xnu=props(2) ! Poisson's ratio
      xl=props(3) ! Characteristic Length
      Gc=props(4) ! Fracture toughness
      
      twomu=E/(1.d0+xnu)
      sixmu=3.d0*twomu
      alamda=twomu*(E-twomu)/(sixmu-2.d0*E)
      
      do i=1,nblock
      phi=tempOld(i)
      psit=stateOld(i,1)
      g=(1.d0-phi)**2.d0+1.d-7
      
C      Update stresses      
        trace  = strainInc(i,1) + strainInc(i,2) + strainInc(i,3)
        stressNew(i,1)=stateOld(i,2)+alamda*trace+twomu*strainInc(i,1)
        stressNew(i,2)=stateOld(i,3)+alamda*trace+twomu*strainInc(i,2)
        stressNew(i,3)=stateOld(i,4)+alamda*trace+twomu*strainInc(i,3)
        stressNew(i,4)=stateOld(i,5)+twomu*strainInc(i,4)
        stressNew(i,5)=stateOld(i,6)+twomu*strainInc(i,5)
        stressNew(i,6)=stateOld(i,7)+twomu*strainInc(i,6)
        
C      Update strain
        stran(i,1:3)=stateOld(i,8:10)+strainInc(i,1:3)
        stran(i,4:6)=stateOld(i,11:13)+2.d0*strainInc(i,4:6)
        stateNew(i,8:13)=stran(i,1:6)
        
C      Calculate strain energy
        psi=0.d0
        do j=1,6
        psi=psi+0.5*stressNew(i,j)*stran(i,j)
        enddo  
        
C      Irreversible damage
       H=max(psit,psi)

       stateNew(i,1)=H
       stateNew(i,2:7)=stressNew(i,:)
      
       stressNew(i,:)=stressNew(i,:)*g        
       
C      Inelastic energy      
        enerInelasNew(i)=-(phi/xl**2-2.d0*(1.d0-phi)*H/(Gc*xl))*
     1  dtArray(1)/density(i)+enerInelasOld(i)
        
C      Phase field rate=0      
      if(tempNew(i) .ge. 1.d0) tempNew(i)=1.d0
      if(tempNew(i) .le. tempOld(i)) tempNew(i)=tempOld(i)
      
C      Update the specific internal energy
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
