    
    subroutine right_B_air()
     use GaussIntegral
     USE NET,only:xyz,NI,N3,I3,PP,NS,WE,SX,dz,flag,num_vector
    USE Space_Domain_Field,ONLY:r_b
    implicit none
    integer i
    R_b=0.0
    !Ey-TEģʽ
    !Ex-TMģʽ
        !����������ߣ��ҵ���������һ������
        do i=1,num_vector !�Ƚ�����е��z���꣬�Ƿ����dz(1)
           
            if((xyz(NI(i,1),3)+xyz(NI(i,2),3))/2==dz(1)) then!����dz��1��˵���ڶ��㡣
               if(flag==1)then !Ey  TEģʽ
                    if((xyz(NI(i,1),1)==xyz(NI(i,2),1)).and.(xyz(NI(i,1),2)/=xyz(NI(i,2),2))) then
                       R_b(i)=DBLE(1.0)
                    end if
               else        !Ex  TMģʽ
                   if((xyz(NI(i,1),2)==xyz(NI(i,2),2)).and.(xyz(NI(i,1),1)/=xyz(NI(i,2),1))) then
                       R_b(i)=DBLE(1.0)
                    end if  
              end if
            end if  
            
            end do
    end subroutine
    
    
    