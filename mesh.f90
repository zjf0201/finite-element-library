  SUBROUTINE GRID()
    use net
    IMPLICIT NONE
    INTEGER I,J,K,N,M,L,K1,KK
    !NΪ��������ĵ�Ԫ����MΪ�������LΪ�ڵ���
    !INTEGER I3((NS-1)*(WE-1)*(SX-1),8),N3((NS-1)*(WE-1)*(SX-1),12),NI((WE-1)*NS*SX+WE*(NS-1)*SX+WE*NS*(SX-1),2)
     REAL(8) SUM
     ALLOCATE(DX(WE),DY(NS),DZ(SX))
    N=(NS-1)*(WE-1)*(SX-1)
    L=NS*WE*SX
    M=(WE-1)*NS*SX+WE*(NS-1)*SX+WE*NS*(SX-1)
    !�õ��ڵ���
      !DIAN
      I3(1,1)=1
      I3(1,2)=2
      I3(1,3)=WE+2
      I3(1,4)=WE+1      
      I3(1,5)=I3(1,1)+WE*NS
      I3(1,6)=I3(1,2)+WE*NS
      I3(1,7)=I3(1,3)+WE*NS
      I3(1,8)=I3(1,4)+WE*NS
      !XIAN
      DO I=2,WE-1
      DO J=1,8
      I3(I,J)=I3(I-1,J)+1
      ENDDO
      ENDDO
      !MIAN
      K=WE
      DO I=2,NS-1
      DO K1=1,WE-1
      DO J=1,8
      I3(K,J)=I3(K-WE+1,J)+WE
      ENDDO
      K=K+1
      ENDDO
      ENDDO
      !TI
      DO I=2,(SX-1)
      DO K1=1,(WE-1)*(NS-1)
      DO J=1,8
      I3(K,J)=I3(K-(WE-1)*(NS-1),J)+WE*NS
      ENDDO
      K=K+1
      ENDDO
      ENDDO 
      !�õ���߱��
      !DIAN
      N3(1,1)=1
      N3(1,2)=WE-1+1
      N3(1,3)=N3(1,1)+(WE-1)*NS
      N3(1,4)=N3(1,2)+(WE-1)*NS
      N3(1,5)=1+(WE-1)*SX*NS
      N3(1,6)=N3(1,5)+(NS-1)*WE
      N3(1,7)=N3(1,5)+1
      N3(1,8)=N3(1,7)+(NS-1)*WE
      N3(1,9)=1+(WE-1)*NS*SX+(NS-1)*WE*SX
      N3(1,10)=N3(1,9)+1
      N3(1,11)=N3(1,9)+WE
      N3(1,12)=N3(1,11)+1
      !XIAN
      DO I=2,WE-1
      DO J=1,12
      N3(I,J)=N3(I-1,J)+1
      ENDDO
      ENDDO
      !MIAN
      K=WE
      DO I=2,NS-1
      DO K1=1,WE-1
      DO J=1,4
      N3(K,J)=N3(K-WE+1,J)+WE-1
      ENDDO
      DO J=5,12
      N3(K,J)=N3(K-WE+1,J)+WE
      ENDDO
      K=K+1
      ENDDO
      ENDDO
      !TI
      DO I=2,(SX-1)
      DO K1=1,(WE-1)*(NS-1)
      DO J=1,4
      N3(K,J)=N3(K-(WE-1)*(NS-1),J)+(WE-1)*NS
      ENDDO
      DO J=5,8
      N3(K,J)=N3(K-(WE-1)*(NS-1),J)+WE*(NS-1)
      ENDDO
      DO J=9,12
      N3(K,J)=N3(K-(WE-1)*(NS-1),J)+WE*NS
      ENDDO
      K=K+1
      ENDDO
      ENDDO 
      !�õ���߶�Ӧ�Ľڵ�
      !X
      DO I=1,WE-1
      NI(I,1)=I
      NI(I,2)=I+1
      ENDDO
      DO I=WE,(WE-1)*NS
      NI(I,1)=NI(I-(WE-1),1)+WE
      NI(I,2)=NI(I-(WE-1),2)+WE
      ENDDO
      DO I=(WE-1)*NS+1,(WE-1)*NS*SX
      NI(I,1)=NI(I-(WE-1)*NS,1)+WE*NS
      NI(I,2)=NI(I-(WE-1)*NS,2)+WE*NS
      ENDDO 
      !Y
      K1=1
      DO I=(WE-1)*NS*SX+1,(WE-1)*NS*SX+WE
      NI(I,1)=K1
      NI(I,2)=K1+WE
      K1=K1+1
      ENDDO
      DO I=(WE-1)*NS*SX+WE+1,(WE-1)*NS*SX+(NS-1)*WE
      NI(I,1)=NI(I-WE,1)+WE
      NI(I,2)=NI(I-WE,2)+WE
      ENDDO
      DO I=(WE-1)*NS*SX+(NS-1)*WE+1,(WE-1)*NS*SX+(NS-1)*WE*SX
      NI(I,1)=NI(I-WE*(NS-1),1)+WE*NS
      NI(I,2)=NI(I-WE*(NS-1),2)+WE*NS
      ENDDO
      !Z
      K1=1
      DO I=(WE-1)*NS*SX+(NS-1)*WE*SX+1,(WE-1)*NS*SX+(NS-1)*WE*SX+WE
      NI(I,1)=K1
      NI(I,2)=K1+WE*NS
      K1=K1+1
      ENDDO
      DO I=(WE-1)*NS*SX+(NS-1)*WE*SX+WE+1,(WE-1)*NS*SX+(NS-1)*WE*SX+WE*NS
      NI(I,1)=NI(I-WE,1)+WE
      NI(I,2)=NI(I-WE,2)+WE
      ENDDO
      DO I=(WE-1)*NS*SX+(NS-1)*WE*SX+WE*NS+1,(WE-1)*NS*SX+(NS-1)*WE*SX+WE*NS*(SX-1)
      NI(I,1)=NI(I-WE*NS,1)+WE*NS
      NI(I,2)=NI(I-WE*NS,2)+WE*NS
      ENDDO
 
      !�õ��ڵ�����
       OPEN(100,FILE=file_xyz)
       read(100,*)
       READ(100,*)
       READ(100,*)
       READ(100,*)(DX(I),I=1,WE)
       !write(*,*)dx
       READ(100,*)
       READ(100,*)(DY(I),I=1,NS)
       READ(100,*)
       READ(100,*)(DZ(I),I=1,SX)
       CLOSE(100)

       DO I=1,WE
       DO J=1,NS
       DO K=1,SX
       XYZ((K-1)*WE*NS+(J-1)*WE+I,1)=DX(I)
       XYZ((K-1)*WE*NS+(J-1)*WE+I,2)=DY(J)
       XYZ((K-1)*WE*NS+(J-1)*WE+I,3)=DZ(K)
       ENDDO
       ENDDO
       ENDDO
 
    ENDSUBROUTINE
    
    
     SUBROUTINE  GET_PP()
   use net,only:pp,NS,SX,WE,P,file_flag,file_xyz,file_freq,file_cedian,flag
   IMPLICIT NONE
   INTEGER I,J,K,N,L,NX,NY,NZ
   CHARACTER(80) FILENAME
   real,allocatable::DX(:),DY(:),DZ(:),xmin(:),xmax(:),ymin(:),ymax(:),zmin(:),zmax(:)
   
    allocate(DX(WE),DY(NS),DZ(SX))
      close(100)
       OPEN(100,FILE=file_xyz)
       read(100,*)
       READ(100,*)
       READ(100,*)
       READ(100,*)(DX(I),I=1,WE)
       !write(*,*)dx
       READ(100,*)
       READ(100,*)(DY(I),I=1,NS)
       READ(100,*)
       READ(100,*)(DZ(I),I=1,SX)
       CLOSE(100)
   close(11)
   open(11,file=file_flag)
   read(11,*)
   read(11,*)
   read(11,*)
   read(11,*)
   read(11,*)
   read(11,*)N
   allocate(xmin(n),xmax(n),ymin(n),ymax(n),zmin(n),zmax(n),P(n))
   do i=1,n
       NX=0;NY=0;NZ=0
       read(11,*) xmin(i),xmax(i),ymin(i),ymax(i),zmin(i),zmax(i),P(i)
    
       do j=1,we
           if(Dx(j)==xmin(i).and.nx==0)then
               xmin(i)=j;NX=NX+1
           elseif(Dx(j)==xmax(i))then
               xmax(i)=j;NX=NX+1
           end if
            if(nx==2)goto 101
       end do
      
101    IF(NX/=2) THEN
          
         print*,"û���ҵ���",i,"�����壡��,x����������"
         pause
       END IF
     
        do j=1,ns
           if(Dy(j)==ymin(i).and.ny==0)then
               ymin(i)=j;ny=ny+1
            
           elseif(Dy(j)==ymax(i))then
               ymax(i)=j;ny=ny+1
    
           end if
           if(ny==2)goto 102
       end do
102           IF(Ny/=2) THEN
               
         print*,"û���ҵ���",i,"�����壡��,y����������"
         pause
              END IF
               
        do j=1,sx
           if(Dz(j)==zmin(i).and.nz==0)then
               zmin(i)=j;nz=nz+1
           
           elseif(Dz(j)==zmax(i))then
               zmax(i)=j;nz=nz+1
             
           end if
            if(nz==2)goto 103
        end do
103         IF(Nz/=2) THEN
         print*,"û���ҵ���",i,"�����壡��,z����������"
         pause
            END IF
             ! print*,zmin(i),zmax(i)
   end do  
   
do L=1,n
    if(P(L)<=0.0000001) P(L)=0.000099  !0.00009

   DO I=zmin(L),zmax(L)-1   !z����
   DO J=ymin(L),YMAX(L)-1    !y����
    do k=xmin(L),XMAX(L)-1   !x����
    PP((I-1)*(WE-1)*(NS-1)+(J-1)*(WE-1)+K)=P(L)
    end do
   ENDDO
   ENDDO
  
end do
   
    ENDSUBROUTINE
   !******************************
  
subroutine measure()    
    use net
    implicit none
    integer i,k,c
   ! real(8) xyz(num_point,3)
    OPEN(11,FILE=file_cedian,ACTION='READ')
    read(11,*) num_ms
    allocate(m_pos(3,num_ms),m_n(num_ms))
    do i=1,num_ms
        read(11,*) c,m_pos(:,i)
    end do
  !  print*,num_ms
    !�ҵ����ĵ�Ԫλ�ã�
    m_n=0
    do i=1,num_ms
    do k=1,num_point
   
        if(Xyz(k,1)==m_pos(1,i))then
           if(xyz(k,2)==m_pos(2,i))then
             if(xyz(k,3)==m_pos(3,i))then
                m_n(i)=k
                goto 101
             end if
           end if   
        end if
        
101 end do
           if(m_n(i)==0)then
               print*,"û���ҵ����λ�ã���"
               pause
           else
               print*,"���ҵ����",i,"��λ��"
           end if
    end do 
     close(11)
    end 
    