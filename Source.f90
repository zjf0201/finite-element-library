 MODULE net
    IMPLICIT NONE
    INTEGER num_ms,num_vector,num_unit,num_point,dof,nnz,ennum,flag
    integer WE,NS,SX
    REAL(8),ALLOCATABLE::pp(:),XYZ(:,:),edge(:,:),m_pos(:,:),m_n(:),DX(:),DY(:),DZ(:),P(:)
    INTEGER,ALLOCATABLE::I3(:,:),N3(:,:),Ni(:,:)
    character(20) file_flag,file_xyz,file_freq,file_cedian
    ENDMODULE net
      
      module sub
    implicit none
    integer mi  !µü´ú´ÎÊý
    real(8),allocatable::V(:,:),Vt(:,:),A(:,:),T(:,:),vtmt(:,:),vtm(:,:)
    real(8) norm
    real(8) fx
    end module
    
    MODULE FREQUENCY
    IMPLICIT NONE
    INTEGER num_freq,N_t,ncmin,ncmax,ncmin1,ncmax1
      real(8) del,omgmin,omgmax,dw
    real(8),allocatable::Freq(:),Time(:),Hfim(:,:),OMG_min_max(:),Exim(:,:)
    ENDMODULE FREQUENCY
    
    