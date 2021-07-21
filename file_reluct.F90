module ReluctInterp
 contains 

  function LinInterpSimple(x0) RESULT(y0)
	USE Types 
	IMPLICIT None
	real(kind=dp) :: x0, y0

        y0 = 769d0 + 4.7516d0*x0*1d-6

  end function LinInterpSimple


function Matrix3x3Inv(Ri) RESULT(Ro)
	USE Types 
	IMPLICIT None
	real(kind=dp), Dimension (3,3) :: Ri, Ro
	real(kind=dp) :: Det, s

  Det = Ri(1,1)*Ri(2,2)*Ri(3,3)-Ri(1,1)*Ri(2,3)*Ri(3,2)- & 
        Ri(1,2)*Ri(2,1)*Ri(3,3)+Ri(1,2)*Ri(2,3)*Ri(3,1)+ & 
        Ri(1,3)*Ri(2,1)*Ri(3,2)-Ri(1,3)*Ri(2,2)*Ri(3,1)		! direct calculatios version

  if (Det .NE. 0) then
     s = 1d0 / Det
  else
     s = 1d25
  end if
    
  Ro(1,1) =  s *   (Ri(2,2)*Ri(3,3) - Ri(3,2)*Ri(2,3))
  Ro(2,1) = -1*s * (Ri(2,1)*Ri(3,3) - Ri(3,1)*Ri(2,3))
  Ro(3,1) =  s *   (Ri(2,1)*Ri(3,2) - Ri(3,1)*Ri(2,2))
    
  Ro(1,2) = -1*s * (Ri(1,2)*Ri(3,3) - Ri(3,2)*Ri(1,3))
  Ro(2,2) =  s *   (Ri(1,1)*Ri(3,3) - Ri(3,1)*Ri(1,3))
  Ro(3,2) = -1*s * (Ri(1,1)*Ri(3,2) - Ri(3,1)*Ri(1,2))

  Ro(1,3) =  s *   (Ri(1,2)*Ri(2,3) - Ri(2,2)*Ri(1,3))
  Ro(2,3) = -1*s * (Ri(1,1)*Ri(2,3) - Ri(2,1)*Ri(1,3))
  Ro(3,3) =  s *   (Ri(1,1)*Ri(2,2) - Ri(2,1)*Ri(1,2))

  end function Matrix3x3Inv



end module ReluctInterp

!-------------------------------------------------------------------------------

SUBROUTINE reluct_func(Model, n, X, Y)

  USE DefUtils
  USE ReluctInterp

  IMPLICIT NONE
  TYPE(Model_t) :: Model
  INTEGER :: n
  REAL(KIND=dp) :: X(*)
  REAL(KIND=dp), POINTER CONTIG :: Y(:,:)
  
  REAL(KIND=dp) :: PSx, PSy, PSz, mi0, v
  REAL(KIND=dp) :: ssx, ssy, ssz
  REAL(KIND=dp) :: mirx, miry, mirz, mix, miy, miz
  REAL(KIND=dp) :: relx, rely, relz

  Real(KIND=dp), dimension (3,3) :: reluct, reluctT, PA


  PSx = X(1)
  PSy = X(2)
  PSz = X(3)			! receive principal stresses


  PA(1,1) = COS(X(4))
  PA(2,1) = COS(X(5))
  PA(3,1) = COS(X(6))
		
  PA(1,2) = COS(X(7))
  PA(2,2) = COS(X(8))		
  PA(3,2) = COS(X(9))

  PA(1,3) = COS(X(10))		
  PA(2,3) = COS(X(11))
  PA(3,3) = COS(X(12))		!  receive angles


! -- end of variable import --

  mi0 = 1.256637061436d-06   	! mi0=4*3.14159265359*1e-7
  v = 3d-1			! Poison ration

  ssx = PSx-v*PSy-v*PSz
  ssy = PSy-v*PSx-v*PSz
  ssz = PSz-v*PSx-v*PSy		! effective stresses calculated from principal stresses 

  mirx = LinInterpSimple(ssx)
  miry = LinInterpSimple(ssy)
  mirz = LinInterpSimple(ssz)

!  mirx = LinInterp(ssx)
!  miry = LinInterp(ssy)
!  mirz = LinInterp(ssz)	! values of diagonal tensor of relative permeability
				! determined from linear interpolation in function LinInterp

  if (mirx .LT. 20) then 
     mirx = 20
  end if
  if (miry .LT. 20) then 
     miry = 20
  end if
  if (mirz .LT. 20) then 
     mirz = 20
  end if			! limit lower value of relative permeabiltity to 20


  mix = mi0*mirx
  miy = mi0*miry
  miz = mi0*mirz		! calculate the permeability


  relx = 1/mix
  rely = 1/miy
  relz = 1/miz			! calculate diagonal values of reluctivity tensor


  reluctT(1,1) = relx
  reluctT(2,2) = rely
  reluctT(3,3) = relz		! reluctivity tensor matrix
  
  reluctT(1,2) = 0d0
  reluctT(1,3) = 0d0
  reluctT(2,1) = 0d0
  reluctT(2,3) = 0d0
  reluctT(3,1) = 0d0
  reluctT(3,2) = 0d0		! set zeros to other

 
 ! rotation matrix R = PA

  reluct=MatMul(PA,MatMul(reluctT,Matrix3x3Inv(PA)))	 ! rotation of tensor reluctT by rotation matrix R
							 ! reluct = R * teluctT * R'

! -- export variables --

  Y(1,1)=reluct(1,1)
  Y(1,2)=reluct(1,2)
  Y(1,3)=reluct(1,3)

  Y(2,1)=reluct(2,1)
  Y(2,2)=reluct(2,2)
  Y(2,3)=reluct(2,3)

  Y(3,1)=reluct(3,1)
  Y(3,2)=reluct(3,2)
  Y(3,3)=reluct(3,3)

END SUBROUTINE reluct_func

