 PROGRAM VFEM_Modelreduction
    IMPLICIT NONE
    INTEGER NS,WE,SX,I,J,N,int1,int2
    call readfile();
 
    call read_flag_XYZ();
    call allocate1();
    CALL GET_PP()
    PRINT*,'�����ʷ֡�����'
    CALL GRID()
    call measure()
    PRINT*,'��ʼ���㡣����'
    call system_clock(int1)
    CALL CACULATE_EM()
     
   call system_clock(int2)
    
     print*,(int2-int1)/10000
    ENDPROGRAM