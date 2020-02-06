 module file
    implicit none
    character*200 Cmdfile,Netfile,Net_num,file_freq
    character*200 file1,file2,file3,file4
 endmodule
       
 MODULE ParaofElecDipole
    IMPLICIT NONE
    REAL(8)::Idl(3)=(/1q0,0q0,0q0/),tau=1,Posion(3)=(/0.0,0.0,0.0/),Posion1(3)=(/1.0,0.0,0.0/)
    CONTAINS
    FUNCTION DELTA(xx)
    IMPLICIT NONE
    REAL(8) DELTA
    REAL(8) xx,xx1
    xx1=ABS(xx)
    IF(xx1>=2*tau) THEN
        DELTA=0.0
    ELSEIF(xx1>=tau) THEN
        DELTA=(xx1-2*tau)**2/(4*tau**3)
    ELSEIF(xx1<tau) THEN
        DELTA=(2*tau*tau-xx1*xx1)/(4*tau**3)
    ENDIF
    ENDFUNCTION DELTA
      function DELTAx(xx)
    implicit none
    real(8) deltax
    real(8) xx
     IF(xx>=2*tau) THEN
        DELTAx=0.0
    ELSEIF(xx>=tau) THEN
        DELTAx=(xx+2*tau)/(tau*tau)-4.0/tau
    ELSEIF(xx>=-tau) THEN
        DELTAx=-(xx+2*tau)/(tau*tau)+2.0/tau
    else if(xx>-2*tau)then
        deltax=(xx+2*tau)/(tau*tau)
    else
        deltax=0.0
    end if
    end function
 ENDMODULE ParaofElecDipole

    
    MODULE FREQUENCY
    IMPLICIT NONE
    
    INTEGER num_freq
    real,allocatable::Freq(:)
    ENDMODULE FREQUENCY
    
    
    module sub
    implicit none
    integer mi  !µü´ú´ÎÊý
    real(8),allocatable::V(:,:),Vt(:,:),A(:,:),T(:,:),vtmt(:,:),vtm(:,:)
    real(8) norm
    real(8) fx
    end module
    
