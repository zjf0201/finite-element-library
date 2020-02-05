 MODULE net
    IMPLICIT NONE
    INTEGER num_unit,num_point,dof,nnz,ennum,nx,ny,nz
    REAL(8),ALLOCATABLE::RHO(:),X(:),Y(:),Z(:)
    real(8) sgmp
    integer info
    INTEGER,ALLOCATABLE::I8(:,:)
    ENDMODULE net

    MODULE LINKED_list
    IMPLICIT NONE
    TYPE NODE
        INTEGER j
        real(8) c
        complex(8) cc
        TYPE(NODE),POINTER::next,prev
    ENDTYPE NODE
    TYPE IndexCeo
        TYPE(NODE),POINTER::head,tail
    ENDTYPE IndexCeo
    TYPE MatLink
        TYPE(IndexCeo),ALLOCATABLE::ht(:) !Õ∑Œ≤÷∏’Î
        INTEGER,ALLOCATABLE::GAi(:)
    ENDTYPE MatLink
    TYPE(NODE),POINTER::Q,P
    TYPE(MatLink) LK,LQ,LQP,Lp,lk1,Lk2,LQ1,LQP1,Lp1,Lgt,le,le1,le2,lkk,lkk1,lkk2,Lg
    ENDMODULE LINKED_list

    MODULE SparseMatrix
    IMPLICIT NONE
    real(8),ALLOCATABLE::A(:),M(:),C(:)
    INTEGER,ALLOCATABLE::IA(:),JA(:),IM(:),JM(:),jc(:),ic(:)
    
    complex(8),ALLOCATABLE::G(:,:),GG(:)
    INTEGER,ALLOCATABLE::Ig(:),Jg(:)
    integer num
    ENDMODULE SparseMatrix

    MODULE Space_Domain_Field
    COMPLEX(8),ALLOCATABLE::Exyz(:),Hy(:),Ex_z(:),Ez_x(:),hz(:),Ex_y(:),Ey_x(:),Ez_y(:),Ey_z(:),Hx(:)
    real(8),ALLOCATABLE::X_right(:),right_b(:),j_right(:)
    ENDMODULE Space_Domain_Field

    
    MODULE GLOBAL_CONST
    IMPLICIT NONE
    REAL(8) pi,miu,epsilon
    PARAMETER(pi=3.1415926535897932384626433832795,miu=1.2566370614359172953850573533118E-6)
    PARAMETER(epsilon=8.85E-12)
    ENDMODULE GLOBAL_CONST
    
    
      MODULE GaussIntegral
    IMPLICIT NONE
    INTEGER n_point
    PARAMETER(n_point=10)
    
    REAL(8),PARAMETER::XIk(n_point)=(/-0.9739,-0.8651,-0.6794,-0.4334,-0.1489,0.1489,0.4334,0.6794,0.8651,0.9739/)
    REAL(8),PARAMETER::Ak(n_point)=(/0.0667,0.1495,0.2191,0.2693,0.2955,0.2955,0.2693,0.2191,0.1495,0.0667/)
   
    ENDMODULE GaussIntegral
    
   