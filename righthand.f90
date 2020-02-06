subroutine right()
   !  右端项，赋予边界值
   !  上边界Ex给1，下边界剖分足够远自动满足
   !   Liu Jiren and zhang Jifeng
    use Space_Domain_Field,only:right_b
    USE NET,ONLY:nx,ny,info,num_point
    implicit none
     
    real x1,y1

    right_b=0.0
    if(info==0) then
    right_b(1:nx*ny)=1.0                        !TM
    else
     right_b(1+num_point:num_point+nx*ny)=1.0   !TE
    end if
    
    end subroutine

    

