subroutine readfile()
    use net,only: file_flag,file_xyz,file_freq,file_cedian
    implicit none
    open(11,file="cmd.txt")
    read(11,*)
    read(11,*) file_flag
    read(11,*)
    read(11,*) file_xyz
    read(11,*)
    read(11,*) file_freq
    read(11,*)
    read(11,*) file_cedian
close(11)
end subroutine
subroutine read_flag_XYZ()
USE SUB,ONLY:MI
USE net,ONLY:we,NS,SX,flag,file_flag,file_xyz
 implicit none
 close(11)
open(11,file=file_flag)
read(11,*)
read(11,*)FLAG
read(11,*)
read(11,*)mi !��������
close(11)  
open(11,file=file_xyz)
read(11,*)we,NS,SX
close(11)

end     
    
    
    
    
    subroutine allocate1()
    use net
    use  Space_Domain_Field
    use sub,only:mi,V,VT,VTMT,T,VTM
    implicit none
    num_point=NS*WE*SX;num_unit=(NS-1)*(WE-1)*(SX-1);num_vector=(WE-1)*NS*SX+WE*(NS-1)*SX+WE*NS*(SX-1)
    allocate(I3((NS-1)*(WE-1)*(SX-1),8),N3((NS-1)*(WE-1)*(SX-1),12),NI((WE-1)*NS*SX+WE*(NS-1)*SX+WE*NS*(SX-1),2))
    allocate(XYZ(NS*WE*SX,3))
    ALLOCATE(PP((NS-1)*(WE-1)*(SX-1))) 
    ALLOCATE(EH(num_vector,1),X_right(num_vector),R_b(num_vector),j_right(num_vector))
   !ģ�ͽ�����ؾ���
   allocate(V(num_vector,mi+1),Vt(mi+1,num_vector),T(mi+1,mi+1))
   allocate(vtmt(num_vector,mi+1),vtm(mi+1,num_vector))
    end subroutine
    
       SUBROUTINE CACULATE_EM()
   use net,only:xyz,Ni,N3,I3,PP,NS,WE,SX,num_ms
   USE SUB,only:fx,mi
   use FREQUENCY
    IMPLICIT NONE
    INTEGER I,J,MM,K1,K,int1,int2
    REAL(8) S
    REAL(8) PI,U0
    !real(8),allocatable::omg(:)
    PI=3.14159265358979323846
    U0=4.E-7*PI
!***********************************
    !************************CSAMTƵ��***************************************8
    open(400,file="freq.txt")
    read(400,*) num_freq  !Ƶ��ĸ���
    allocate(freq(num_freq),Hfim(num_ms,num_freq))
    do i=1,num_freq
        read(400,*)freq(i)
    end do
    !******************************��һ������װAA��M��ע�⣬AA��M��Ϊʵ������*********************************
     fx=-dble(2.0)*pi*sqrt(freq(1)*freq(num_freq))

      CALL HECHENG()
    !******************************�ڶ�������װright_b,J����**************************************************
   CALL right_B_air() 
   close(11)
   open(11,file="time.txt")
     call system_clock(int1)
      print*,'����ģ�ͽ���...'
     !****************************����������M-1J��x_right*******************************************************
   
    !*******���Ĳ�����J_right��x,VTMT,VTM,��C����*******
    
    !*******���岽����T����*******
   call model_reduction()
     call system_clock(int2)
     write(11,*) "MR",(int2-int1)/10000
     close(11)
    print*,'ģ�ͽ�����ϣ�'
    print* 
    print*,'���ڼ���糡ֵ...'
  call subspace()
  
    ENDSUBROUTINE
!********************************
!***************����ϳɣ�������װAA��M����**************
    SUBROUTINE HECHENG() 
         use SparseMatrix,only:AA,IA,JA,M,IM,JM
         use net,only:XYZ,N3,NI,I3,NS,SX,WE,PP
         use Space_Domain_Field,only:EH
      IMPLICIT NONE
      INTEGER L,M1,M2,I,J,K1,K2,ND,N1((WE-1)*NS*SX+WE*(NS-1)*SX+WE*NS*(SX-1))
      REAL(8) S
      real(8) EA(12,12),Ma(12,12)
   
      TYPE NODE
        real(8) A
        INTEGER KI,KJ
        TYPE(NODE),POINTER:: NEXT
      END TYPE NODE
      
      TYPE NODE1
        TYPE(NODE),POINTER:: HEAD
      END TYPE NODE1   
         
      TYPE(NODE1),ALLOCATABLE:: H(:),Hm(:)
      TYPE(NODE),POINTER:: Q,P,P1,Q1
      allocate(IA((WE-1)*NS*SX+WE*(NS-1)*SX+WE*NS*(SX-1)+1),IM((WE-1)*NS*SX+WE*(NS-1)*SX+WE*NS*(SX-1)+1))
      ALLOCATE(H((WE-1)*NS*SX+WE*(NS-1)*SX+WE*NS*(SX-1)),Hm((WE-1)*NS*SX+WE*(NS-1)*SX+WE*NS*(SX-1)))
      DO I=1,(WE-1)*NS*SX+WE*(NS-1)*SX+WE*NS*(SX-1)
       NULLIFY(H(I).HEAD)   ! ʹָ�������ͷָ��ָ���ָ�� ����AA����
       NULLIFY(Hm(I).HEAD)   ! ʹָ�������ͷָ��ָ���ָ�� ����M����
      ENDDO
      !******************************************��ʼ��תAA����****************************************
      M1=0;M2=0
      DO L=1,(NS-1)*(WE-1)*(SX-1)
       CALL UKE1(L,XYZ,N3,NI,I3,S,NS,SX,WE,EA,PP,Ma)
        DO I=1,12
        DO J=I,12
            IF(ABS(EA(I,J))>0.0) THEN
                IF(N3(L,I)<=N3(L,J))THEN
                K1=N3(L,I)
                K2=N3(L,J)
                ELSE
                K1=N3(L,J)
                K2=N3(L,I)
                ENDIF   ! ��֤�洢�����ϰ벿��
       ! �ж��Ƿ�Ϊ�ձ�ͷ----------------------------------------------
                IF(.NOT.ASSOCIATED(H(K1).HEAD)) THEN  !���ָ��Ϊ��ָ��
                    ALLOCATE(Q)    !�����ַ�ռ�
                    Q.A=EA(I,J)    !��Nodeָ�����ݸ�ֵ
                    Q.KI=K1
                    Q.KJ=K2
                    NULLIFY(Q.NEXT) ! ��Q��Nextָ��ָ���
                  H(K1).HEAD=>Q     ! H��k1��.ͷָ��ָ��Q
                  M1=M1+1
                  N1(K1)=1
                ELSE
                    P1=>H(K1).HEAD
                    NULLIFY(P)
                  !Ѱ��������Ԫ�ص���ͬλ��
                 WW: DO WHILE(ASSOCIATED(P1)) !���P1��Ϊ��
                       IF(K2==P1.KJ)THEN
                         P=>P1
                        EXIT WW   !�ҵ���ͬλ��Ԫ���˳����Ҳ���ָ���ƶ����¸�����
                       ELSE
                         P1=>P1.NEXT     !����ĺ��˷���
                       ENDIF
                     ENDDO WW
                  !������ظ��������ۼӣ����û�д浽��һ���ڵ�
                     IF(ASSOCIATED(P))THEN    !���P��ΪΪ�գ���Ϊ�գ���
                       P.A=P.A+EA(I,J)
                     ELSE
                          ALLOCATE(Q)
                           Q.A=EA(I,J)
                           Q.KI=K1
                           Q.KJ=K2
                           NULLIFY(Q.NEXT)
                          IF(Q.KJ<=H(K1).HEAD.KJ)THEN
                            Q.NEXT=>H(K1).HEAD
                            H(K1).HEAD=>Q
                            M1=M1+1
                            N1(K1)=N1(K1)+1    
                          ELSE
                            NULLIFY(P1) 
                            P1=>H(K1).HEAD
                  WW2:     DO WHILE(ASSOCIATED(P1).AND.Q.KJ>P1.KJ)
                             Q1=>P1
                             P1=>P1.NEXT  
                             EXIT WW2             
                           ENDDO WW2
                           IF(ASSOCIATED(P1))THEN 
                           M1=M1+1
                           N1(K1)=N1(K1)+1
                            Q.NEXT=>Q1.NEXT
                            Q1.NEXT=>Q
                           ELSE
                           M1=M1+1
                           N1(K1)=N1(K1)+1
                           Q1.NEXT=>Q
                           ENDIF
                         ENDIF
                     ENDIF
                ENDIF  
             ENDIF 
           ENDDO
        ENDDO
      ENDDO  
    ! ת��Ϊ��̬����  !�ͷ�����--------------------------------
       ALLOCATE(JA(M1),AA(M1))
  
       IA(1)=1
       J=1
       DO I=1,(WE-1)*NS*SX+WE*(NS-1)*SX+WE*NS*(SX-1)
          P=>H(I).HEAD
          DO WHILE(ASSOCIATED(P))
             AA(J)=P.A
             JA(J)=P.KJ
             Q1=>P.NEXT
             DEALLOCATE(P)
             P=>Q1
             J=J+1 
          ENDDO
             IA(I+1)=J
             IF(N1(I)/=(IA(I+1)-IA(I)).OR.N1(I)==0) THEN
             PRINT*,'YOUCUOWU'
             PAUSE
             ENDIF
       ENDDO
       DEALLOCATE(H)
       !***********************************************************��װM����**********************************************************
        M1=0;N1=0
      DO L=1,(NS-1)*(WE-1)*(SX-1)
       CALL UKE1(L,XYZ,N3,NI,I3,S,NS,SX,WE,EA,PP,Ma)
        DO I=1,12
        DO J=I,12
            IF(ABS(MA(I,J))>0.0) THEN
                IF(N3(L,I)<=N3(L,J))THEN
                K1=N3(L,I)
                K2=N3(L,J)
                ELSE
                K1=N3(L,J)
                K2=N3(L,I)
                ENDIF   ! ��֤�洢�����ϰ벿��
       ! �ж��Ƿ�Ϊ�ձ�ͷ----------------------------------------------
                IF(.NOT.ASSOCIATED(Hm(K1).HEAD)) THEN  !���ָ��Ϊ��ָ��
                    ALLOCATE(Q)    !�����ַ�ռ�
                    Q.A=MA(I,J)    !��Nodeָ�����ݸ�ֵ
                    Q.KI=K1
                    Q.KJ=K2
                    NULLIFY(Q.NEXT) ! ��Q��Nextָ��ָ���
                  Hm(K1).HEAD=>Q     ! H��k1��.ͷָ��ָ��Q
                  M1=M1+1
                  N1(K1)=1
                ELSE
                    P1=>Hm(K1).HEAD
                    NULLIFY(P)
                  !Ѱ��������Ԫ�ص���ͬλ��
                 WWW: DO WHILE(ASSOCIATED(P1)) !���P1��Ϊ��
                       IF(K2==P1.KJ)THEN
                         P=>P1
                        EXIT WWW   !�ҵ���ͬλ��Ԫ���˳����Ҳ���ָ���ƶ����¸�����
                       ELSE
                         P1=>P1.NEXT     !����ĺ��˷���
                       ENDIF
                     ENDDO WWW
                  !������ظ��������ۼӣ����û�д浽��һ���ڵ�
                     IF(ASSOCIATED(P))THEN    !���P��ΪΪ�գ���Ϊ�գ���
                       P.A=P.A+MA(I,J)
                     ELSE
                          ALLOCATE(Q)
                           Q.A=MA(I,J)
                           Q.KI=K1
                           Q.KJ=K2
                           NULLIFY(Q.NEXT)
                          IF(Q.KJ<=Hm(K1).HEAD.KJ)THEN
                            Q.NEXT=>Hm(K1).HEAD
                            Hm(K1).HEAD=>Q
                            M1=M1+1
                            N1(K1)=N1(K1)+1    
                          ELSE
                            NULLIFY(P1) 
                            P1=>Hm(K1).HEAD
                  WWW2:     DO WHILE(ASSOCIATED(P1).AND.Q.KJ>P1.KJ)
                             Q1=>P1
                             P1=>P1.NEXT  
                             EXIT WWW2             
                           ENDDO WWW2
                           IF(ASSOCIATED(P1))THEN 
                           M1=M1+1
                           N1(K1)=N1(K1)+1
                            Q.NEXT=>Q1.NEXT
                            Q1.NEXT=>Q
                           ELSE
                           M1=M1+1
                           N1(K1)=N1(K1)+1
                           Q1.NEXT=>Q
                           ENDIF
                         ENDIF
                     ENDIF
                ENDIF  
             ENDIF 
           ENDDO
        ENDDO
      ENDDO  
    ! ת��Ϊ��̬����  !�ͷ�����--------------------------------
       ALLOCATE(JM(M1),M(M1))
       IM(1)=1
       J=1
       DO I=1,(WE-1)*NS*SX+WE*(NS-1)*SX+WE*NS*(SX-1)
          P=>Hm(I).HEAD
          DO WHILE(ASSOCIATED(P))
             M(J)=P.A
             JM(J)=P.KJ
             Q1=>P.NEXT
             DEALLOCATE(P)
             P=>Q1
             J=J+1 
          ENDDO
             IM(I+1)=J
             IF(N1(I)/=(IM(I+1)-IM(I)).OR.N1(I)==0) THEN
             PRINT*,'YOUCUOWU'
             PAUSE
             ENDIF
       ENDDO
       DEALLOCATE(Hm)
       
      ENDSUBROUTINE
      !***************************************
      SUBROUTINE UKE1(K,XYZ,N3,NI,I3,S,NS,SX,WE,EA,PP,Ma)
      USE SUB,only:fx
      IMPLICIT NONE
      INTEGER I,J,K,NS,SX,WE
      INTEGER N3((NS-1)*(WE-1)*(SX-1),12),NI((WE-1)*NS*SX+WE*(NS-1)*SX+WE*NS*(SX-1),2),I3((NS-1)*(WE-1)*(SX-1),8)
      REAL(8) A1(4,4),A2(4,4),A3(4,4),A4(4,4),PP((NS-1)*(WE-1)*(SX-1)),A31(4,4),LX,LY,LZ,PI,U0,S,XYZ(NS*WE*SX,3)
      real(8) EA(12,12),EA1(12,12),EA2(12,12),II
      real(8) Ma(12,12)
      DATA A1/2.,-2.,1.,-1.,-2.,2.,-1.,1.,1.,-1.,2.,-2.,-1.,1.,-2.,2./
      DATA A2/2.,1.,-2.,-1.,1.,2.,-1.,-2.,-2.,-1.,2.,1.,-1.,-2.,1.,2./
      DATA A3/2.,-2.,1.,-1.,1.,-1.,2.,-2.,-2.,2.,-1.,1.,-1.,1.,-2.,2./
      DATA A4/4.,2.,2.,1.,2.,4.,1.,2.,2.,1.,4.,2.,1.,2.,2.,4./
      PI=3.14159265358979323846
      U0=4.E-7*PI
      II=CMPLX(0,1.0)
      A31=TRANSPOSE(A3)
      LX=ABS(XYZ(I3(K,1),1)-XYZ(I3(K,2),1))
      LY=ABS(XYZ(I3(K,1),2)-XYZ(I3(K,4),2))
      LZ=ABS(XYZ(I3(K,1),3)-XYZ(I3(K,5),3))
      DO I=1,4
      DO J=1,4
      EA1(I,J)=LX*LZ/6./LY*A1(I,J)+LX*LY/6./LZ*A2(I,J)
      EA1(I+4,J+4)=LX*LY/6./LZ*A1(I,J)+LZ*LY/6./LX*A2(I,J)
      EA1(I+8,J+8)=LY*LZ/6./LX*A1(I,J)+LX*LZ/6./LY*A2(I,J)
      EA1(I+4,J)=-LZ/6.*A31(I,J)
    
      EA1(I+8,J)=-LY/6.*A3(I,J)
      EA1(I+8,J+4)=-LX/6.*A31(I,J)
      
      EA1(I,J+4)=-LZ/6.*A3(I,J)   
  
     EA1(I,J+8)=-LY/6.*A31(I,J)
      EA1(I+4,J+8)=-LX/6.*A3(I,J)
      
      ENDDO
      ENDDO
      EA2=0.0
      DO I=1,4
      DO J=1,4
      EA2(I,J)=U0*PP(K)*LX*LY*LZ/36.0*A4(I,J)
      EA2(I+4,J+4)=U0*PP(K)*LX*LY*LZ/36.0*A4(I,J)
      EA2(I+8,J+8)=U0*PP(K)*LX*LY*LZ/36.0*A4(I,J)
      ENDDO
      ENDDO 
      DO I=1,12
      DO J=1,12
      EA(I,J)=EA1(I,J)-fx*EA2(i,j)
      MA(i,j)=EA2(i,j)
      ENDDO
      ENDDO
      ENDSUBROUTINE
     !******����������ģ�ͽ��׵��Ҷ���*****
subroutine mj_inv()
!    ����M��*right_b=x_right
!    �൱�����M*��M^-1 *right_b��=right_b��
!    ģ�ͽ��׽ⷽ��Ϊ  (C-fx*M)E=M^-1 right_b
! Liu Jiren and Zhang Jifeng

    use Space_Domain_Field,only:r_b,X_right
    use SparseMatrix,only:M,IM,JM
    use net,only:num_vector
     IMPLICIT NONE
    INTEGER lda, ldb, maxfct, mnum, mtype, phase, nrhs, error, msglvl
    INTEGER  iparm(64),n,i,ii
    INTEGER*8 pt(64)
    integer,allocatable::perm(:)
    integer,external:: mkl_get_max_threads
    REAL*8 ddum
    lda=Im(num_vector+1)-1;ldb=num_vector
    allocate( perm(ldb) )
    pt=0 !ָ���ʼ��
    maxfct=1;mnum=1;mtype=2;phase=13;n=ldb;perm=0;nrhs=1;x_right=dble(0.0)
    iparm(1) = 0
    error = 0 ! ������Ϣ��ֵ��0 
    msglvl = 0 ! ����ʾͳ����Ϣ 


    call pardiso (pt, maxfct, mnum, mtype, phase, n, m, im, jm, perm, nrhs, iparm, msglvl, r_b, x_right, error)
    if(error/=0) print*,'ʧ��'
    phase=-1
    call pardiso (pt, maxfct, mnum, mtype, phase, n, m, im, jm, perm, nrhs, iparm, msglvl, r_b, x_right, error)

    deallocate( perm )

    end 
    
    !*****************��ʼģ�ͽ���**********************
subroutine model_reduction()
!    ���õ����ظ�����fx
!     m�λش�
!    ����������AD���̵ı�����ʽ����
!    ģ�ͽ��׽ⷽ��Ϊ  (C-fx*M)E=M^-1 b =>(C-fx*M)�洢��A��IA,JA��
!    Liu Jiren and Zhang Jifeng
    use Space_Domain_Field,only:X_right,J_right,R_b
    use SparseMatrix,only:M,IM,JM,AA,IA,JA,C,JC,IC
    use sub,only:mi,v,vtm,norm,T,vt
    use net,only:num_vector
     IMPLICIT NONE
    INTEGER lda, ldb, maxfct, mnum, mtype, phase, nrhs, error, msglvl
     real(8) vi(num_vector),vj(num_vector),y(num_vector),vtcvI(mi+1),mvj(num_vector),mvj1(num_vector),mvj2(num_vector)
     real(8) y1(num_vector),y2(num_vector)
    INTEGER  iparm(64),n,step
    INTEGER*8 pt(64),ierr
    integer,allocatable::perm(:)
    integer,external:: mkl_get_max_threads
    REAL*8 ddum,x(num_vector)
    real normm,dot
    integer,allocatable:: iw(:)
    integer i,j,p,ii
     call MJ_inv()

  
    lda=IA(num_vector+1)-1;ldb=num_vector
     allocate( perm(ldb),iw(ldb) )
    pt=0 !ָ���ʼ��
    maxfct=1;mnum=1;mtype=-2;phase=13;n=ldb;perm=0;nrhs=1;x=dble(0.0)
    iparm(1) = 0
    error = 0 ! ������Ϣ��ֵ��0 
    msglvl = 0 ! ����ʾͳ����Ϣ 
    mvj1=dble(0.0);mvj2=dble(0.0);y1=dble(0.0);y2=dble(0.0)
  
    norm=normm(x_right,num_vector)               ! ��||b||m�������ı�����ʽ
    ! ������
    x_right(:)=x_right(:)/norm                    ! ��V1
     V(:,1)=x_right(:)                            ! ����v1
     call atmux ( n, x_right, mvj1, m, jm, im) 
     call amux ( n, x_right, mvj2,m, jm, im )      ! ���һ���Ҷ��� M*V1������mvj
     mvj=mvj1+mvj2
     step=mi/10 
     do j=1,mi,1 
       if(mod(j,step)==0) print 20,int(real(j)/mi*95)
20      format(i3,'% Completed...')
        if(j==1)then
                                                  ! pardiso ��⣬�洢��X��
           call pardiso (pt, maxfct, mnum, mtype, phase, n, Aa, ia, ja, perm, nrhs, iparm, msglvl, mvj, x, error)
           if(error/=0) print*,'ʧ��'
        else                                      ! m�λش�
             call pardiso (pt, maxfct, mnum, mtype,33, n, Aa, ia, ja, perm, nrhs, iparm, msglvl, mvj, x, error)
             if(error/=0) print*,'ʧ��'
        end if

        do i=1,j
            vi(:)=v(:,i)
            x(:)=x(:)-dot(x,vi,num_vector)*Vi(:)  ! ��� x=x-(x,vi)*vi
        end do
           V(:,j+1)=x(:)/normm(x,num_vector)      ! ��� x/||x||m
            vi(:)= V(:,j+1)
           mvj1=0.0;mvj2=0.0
         call atmux ( n, vi, mvj1, m, jm, im)       !������������˻�
         call amux ( n, vi, mvj2,m, jm, im )        !������������˻�
         mvj=mvj1+mvj2                              ! ����M*Vj+1
       end do
      phase=-1                                      ! �ͷ�
         call pardiso (pt, maxfct, mnum, mtype, phase, n, Aa, ia, ja, perm, nrhs, iparm, msglvl,mVj, x, error)
   print*,"���ڼ���VTM"

    vt=transpose(V)  
    call get_c();                                   ! �õ�c����C,IC��JC
     do i=1,mi+1
        vi(:)=v(:,i)    
      
        call atmux(num_vector, vi, y1, c, jc, ic)
        call amux (num_vector, vi, y2, c, jc, ic)
        y(:)=y1(:)+y2(:) 
        vtcvi=matmul(vt,y)                           !�õ�T�ĵ�i������          
      
        !
        T(:,i)=vtcVI(:)                              !�ô�T����
    
        
     end do
     
    
      deallocate( perm,vtm,AA,IA,JA,M,JM,IM,C,IC,JC,X_right,R_b,j_right)
    end subroutine
    
  
!*******************�õ�C����************************
subroutine get_c()
    use sub,only:fx
    use SparseMatrix,only:M,IM,JM,AA,IA,JA,C,JC,IC
    use net,only:num_vector
    implicit none
    integer nrow,ncol,job,nzmax,ierr,i,nm
    integer,allocatable:: iw(:)
    ierr=0
    nzmax=Ia(num_vector+1)-1;nm=Im(num_vector+1)-1
    nrow=num_vector;ncol=num_vector
   allocate(jc(nzmax),ic(nrow+1),c(nzmax),iw(ncol))

   M(:)=fx*M(:)
    call aplb(nrow,ncol,1,Aa,ja,ia,m,jm,im,C,jc,ic,nzmax,iw,ierr)
    if(ierr>0)then
        print*,ierr,"ת��C����ʧ��"
    end if
    M(:)=M(:)/fx     !�ָ�M

    end subroutine