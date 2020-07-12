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
read(11,*)mi !迭代次数
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
   !模型降阶相关矩阵
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
    !************************CSAMT频点***************************************8
    open(400,file="freq.txt")
    read(400,*) num_freq  !频点的个数
    allocate(freq(num_freq),Hfim(num_ms,num_freq))
    do i=1,num_freq
        read(400,*)freq(i)
    end do
    !******************************第一步，组装AA和M，注意，AA和M均为实数矩阵*********************************
     fx=-dble(2.0)*pi*sqrt(freq(1)*freq(num_freq))

      CALL HECHENG()
    !******************************第二步、组装right_b,J矩阵**************************************************
   CALL right_B_air() 
   close(11)
   open(11,file="time.txt")
     call system_clock(int1)
      print*,'正在模型降阶...'
     !****************************第三步、求M-1J，x_right*******************************************************
   
    !*******第四步、求J_right，x,VTMT,VTM,求C矩阵*******
    
    !*******第五步、求T矩阵*******
   call model_reduction()
     call system_clock(int2)
     write(11,*) "MR",(int2-int1)/10000
     close(11)
    print*,'模型降阶完毕！'
    print* 
    print*,'正在计算电场值...'
  call subspace()
  
    ENDSUBROUTINE
!********************************
!***************矩阵合成，用于组装AA和M矩阵**************
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
       NULLIFY(H(I).HEAD)   ! 使指针数组的头指针指向空指针 用于AA矩阵
       NULLIFY(Hm(I).HEAD)   ! 使指针数组的头指针指向空指针 用于M矩阵
      ENDDO
      !******************************************开始组转AA矩阵****************************************
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
                ENDIF   ! 保证存储矩阵上半部分
       ! 判断是否为空表头----------------------------------------------
                IF(.NOT.ASSOCIATED(H(K1).HEAD)) THEN  !如果指针为空指针
                    ALLOCATE(Q)    !分配地址空间
                    Q.A=EA(I,J)    !给Node指针数据赋值
                    Q.KI=K1
                    Q.KJ=K2
                    NULLIFY(Q.NEXT) ! 让Q的Next指针指向空
                  H(K1).HEAD=>Q     ! H（k1）.头指针指向Q
                  M1=M1+1
                  N1(K1)=1
                ELSE
                    P1=>H(K1).HEAD
                    NULLIFY(P)
                  !寻找数组中元素的相同位置
                 WW: DO WHILE(ASSOCIATED(P1)) !如果P1不为空
                       IF(K2==P1.KJ)THEN
                         P=>P1
                        EXIT WW   !找到相同位置元素退出，找不到指针移动到下个数据
                       ELSE
                         P1=>P1.NEXT     !链表的后退方向
                       ENDIF
                     ENDDO WW
                  !如果有重复，进行累加，如果没有存到后一个节点
                     IF(ASSOCIATED(P))THEN    !如果P不为为空，若为空？？
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
    ! 转化为静态数组  !释放链表--------------------------------
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
       !***********************************************************组装M矩阵**********************************************************
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
                ENDIF   ! 保证存储矩阵上半部分
       ! 判断是否为空表头----------------------------------------------
                IF(.NOT.ASSOCIATED(Hm(K1).HEAD)) THEN  !如果指针为空指针
                    ALLOCATE(Q)    !分配地址空间
                    Q.A=MA(I,J)    !给Node指针数据赋值
                    Q.KI=K1
                    Q.KJ=K2
                    NULLIFY(Q.NEXT) ! 让Q的Next指针指向空
                  Hm(K1).HEAD=>Q     ! H（k1）.头指针指向Q
                  M1=M1+1
                  N1(K1)=1
                ELSE
                    P1=>Hm(K1).HEAD
                    NULLIFY(P)
                  !寻找数组中元素的相同位置
                 WWW: DO WHILE(ASSOCIATED(P1)) !如果P1不为空
                       IF(K2==P1.KJ)THEN
                         P=>P1
                        EXIT WWW   !找到相同位置元素退出，找不到指针移动到下个数据
                       ELSE
                         P1=>P1.NEXT     !链表的后退方向
                       ENDIF
                     ENDDO WWW
                  !如果有重复，进行累加，如果没有存到后一个节点
                     IF(ASSOCIATED(P))THEN    !如果P不为为空，若为空？？
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
    ! 转化为静态数组  !释放链表--------------------------------
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
     !******第三步、求模型降阶的右端项*****
subroutine mj_inv()
!    计算M逆*right_b=x_right
!    相当于求解M*（M^-1 *right_b）=right_b；
!    模型降阶解方程为  (C-fx*M)E=M^-1 right_b
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
    pt=0 !指针初始化
    maxfct=1;mnum=1;mtype=2;phase=13;n=ldb;perm=0;nrhs=1;x_right=dble(0.0)
    iparm(1) = 0
    error = 0 ! 错误信息初值置0 
    msglvl = 0 ! 不显示统计信息 


    call pardiso (pt, maxfct, mnum, mtype, phase, n, m, im, jm, perm, nrhs, iparm, msglvl, r_b, x_right, error)
    if(error/=0) print*,'失败'
    phase=-1
    call pardiso (pt, maxfct, mnum, mtype, phase, n, m, im, jm, perm, nrhs, iparm, msglvl, r_b, x_right, error)

    deallocate( perm )

    end 
    
    !*****************开始模型降阶**********************
subroutine model_reduction()
!    利用单个重复极点fx
!     m次回代
!    采用正交化AD过程的变体形式构建
!    模型降阶解方程为  (C-fx*M)E=M^-1 b =>(C-fx*M)存储在A，IA,JA中
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
    pt=0 !指针初始化
    maxfct=1;mnum=1;mtype=-2;phase=13;n=ldb;perm=0;nrhs=1;x=dble(0.0)
    iparm(1) = 0
    error = 0 ! 错误信息初值置0 
    msglvl = 0 ! 不显示统计信息 
    mvj1=dble(0.0);mvj2=dble(0.0);y1=dble(0.0);y2=dble(0.0)
  
    norm=normm(x_right,num_vector)               ! 求||b||m，范数的变体形式
    ! 正交化
    x_right(:)=x_right(:)/norm                    ! 求V1
     V(:,1)=x_right(:)                            ! 赋给v1
     call atmux ( n, x_right, mvj1, m, jm, im) 
     call amux ( n, x_right, mvj2,m, jm, im )      ! 求第一个右端项 M*V1，赋给mvj
     mvj=mvj1+mvj2
     step=mi/10 
     do j=1,mi,1 
       if(mod(j,step)==0) print 20,int(real(j)/mi*95)
20      format(i3,'% Completed...')
        if(j==1)then
                                                  ! pardiso 求解，存储在X中
           call pardiso (pt, maxfct, mnum, mtype, phase, n, Aa, ia, ja, perm, nrhs, iparm, msglvl, mvj, x, error)
           if(error/=0) print*,'失败'
        else                                      ! m次回代
             call pardiso (pt, maxfct, mnum, mtype,33, n, Aa, ia, ja, perm, nrhs, iparm, msglvl, mvj, x, error)
             if(error/=0) print*,'失败'
        end if

        do i=1,j
            vi(:)=v(:,i)
            x(:)=x(:)-dot(x,vi,num_vector)*Vi(:)  ! 求解 x=x-(x,vi)*vi
        end do
           V(:,j+1)=x(:)/normm(x,num_vector)      ! 求解 x/||x||m
            vi(:)= V(:,j+1)
           mvj1=0.0;mvj2=0.0
         call atmux ( n, vi, mvj1, m, jm, im)       !求下三角与其乘积
         call amux ( n, vi, mvj2,m, jm, im )        !求上三角与其乘积
         mvj=mvj1+mvj2                              ! 计算M*Vj+1
       end do
      phase=-1                                      ! 释放
         call pardiso (pt, maxfct, mnum, mtype, phase, n, Aa, ia, ja, perm, nrhs, iparm, msglvl,mVj, x, error)
   print*,"正在计算VTM"

    vt=transpose(V)  
    call get_c();                                   ! 得到c矩阵，C,IC，JC
     do i=1,mi+1
        vi(:)=v(:,i)    
      
        call atmux(num_vector, vi, y1, c, jc, ic)
        call amux (num_vector, vi, y2, c, jc, ic)
        y(:)=y1(:)+y2(:) 
        vtcvi=matmul(vt,y)                           !得到T的第i列数据          
      
        !
        T(:,i)=vtcVI(:)                              !得打T矩阵
    
        
     end do
     
    
      deallocate( perm,vtm,AA,IA,JA,M,JM,IM,C,IC,JC,X_right,R_b,j_right)
    end subroutine
    
  
!*******************得到C矩阵************************
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
        print*,ierr,"转换C矩阵失败"
    end if
    M(:)=M(:)/fx     !恢复M

    end subroutine