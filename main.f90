
    PROGRAM LJR_CSEM
    !Ö÷³ÌÐò
     use  lapack95
     IMPLICIT NONE
     INTEGER i,j,int1,int2
     real(8) t1,t2,y0
     CHARACTER*200 Netfile,Net_num
     !****read**********************************
 
    
  call system_clock(int1)  
          CALL READNET()
          CALL Allocate1() 
          CALL CAL_SPACE_FILED()
  
   call system_clock(int2)
    
     print*,(int2-int1)/10000
    call system('pause')
    end program
