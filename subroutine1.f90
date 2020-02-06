SUBROUTINE ReadNet()
    !读入网格参数
    !网格剖分数量
    !读入频点个数
     USE NET
     use FREQUENCY  
     use GLOBAL_CONST,only:pi
     use sub,only:fx,mi
     use file
IMPLICIT NONE
    INTEGER i,c
    character*80 cc
    
!*************读入文件名*********** 
   ! read(*,*) cmdfile
    cmdfile='cmdfile.txt'
     OPEN(11,FILE=cmdfile)
    read(11,*) cc,Netfile
    read(11,*) cc,Net_num
    read(11,*) cc,file_freq
    read(11,*) cc,file1
    read(11,*) cc,file2
    read(11,*) cc,file3
    read(11,*) cc,file4
    CLOSE(11)  
!*************读入单元参数************   
    OPEN(11,FILE=Netfile)
    READ(11,*) num_point
 
    ALLOCATE(X(num_point),Y(num_point),Z(num_point)) !分别为节点坐标数组
    DO I=1,num_point
       READ(11,*) c,X(i),y(i),z(i)
    end do
       READ(11,*) num_unit
    ALLOCATE(I8(8,num_unit),RHO(num_unit))
    do i=1,num_unit
       READ(11,*) I8(:,i),RHO(i)
    end do
    CLOSE(11)
    ennum= 8;
	dof = 3 ;
!************读入剖分大小*********  
    OPEN(11,FILE=Net_num)
    read(11,*) cc,nx
    read(11,*) cc,ny
    read(11,*) cc,nz
close(11)

!*************读入频点**************
  OPEN(11,FILE=file_freq)
    read(11,*) num_freq,info,mi
    allocate(Freq(num_freq))
    do i=1,num_freq
    read(11,*)Freq(i)
    end do
close(11)
!*************计算单个频点**********
   fx=-dble(2.0)*pi*dsqrt(freq(1)*freq(num_freq))
    ENDSUBROUTINE
    
    
    
 !****************************开辟动态数组*************************************
SUBROUTINE Allocate1()
    USE NET,ONLY:num_point,dof
    USE LINKED_list,ONLY:LK,LQ,LQP,Lp,lk1,Lk2,le1,le2,le
    USE SparseMatrix,ONLY:IA,IM,ic
    USE Space_Domain_Field,ONLY:Exyz,right_b,x_right,J_right,Ex_z,Ez_x,Hy,hz,Ex_y,Ey_x,Ez_y,Ey_z,HX
    use sub
    use FREQUENCY
    IMPLICIT NONE
    !分别为链表行索引指针、总刚矩阵行元素个数索引分配空间
    ALLOCATE(LK.ht(num_point),LK.GAi(num_point),LQ.ht(num_point),LQ.GAi(num_point),LQP.ht(num_point),LQP.GAi(num_point),LP.ht(num_point),LP.GAi(num_point))
    ALLOCATE(Le.ht(num_point),Le.GAi(num_point),Le1.ht(num_point),Le1.GAi(num_point),Le2.ht(num_point),Le2.GAi(num_point))
    !为行起点索引数组分配空间
    allocate(LK1.ht(num_point),LK1.GAi(num_point),LK2.ht(num_point),LK2.GAi(num_point))

    ALLOCATE(IA(dof*num_point+1))
    ALLOCATE(IM(dof*num_point+1))

    
    ALLOCATE(Exyz(dof*num_point),X_right(dof*num_point),right_b(dof*num_point),j_right(dof*num_point))
   !模型降阶相关矩阵
   allocate(V(3*num_point,mi+1),Vt(mi+1,3*num_point),T(mi+1,mi+1))
   allocate(vtmt(3*num_point,mi+1),vtm(mi+1,3*num_point))
   !未知量
   allocate(Ex_z(num_point),Ez_x(num_point),Hy(num_point),Hz(num_point),Hx(num_point),Ey_z(num_point),Ez_y(num_point))
    allocate(Ex_y(num_point),Ey_x(num_point))
    
    
    ENDSUBROUTINE
    
   
    
  !****************************主体计算*************************************   
SUBROUTINE CAL_SPACE_FILED()
    IMPLICIT NONE
    integer j
    print*
    print*,'正在组装整体矩阵...'
    !********第一步、组装A、M********
    CALL GK() !集成总体刚度矩阵，存储在链表中，链表元素行索引为head
    
    CALL CONVERT()  !将链表中的刚阵转存在指针数组中,并释放链表
 
    print*,'正在组装右端...'
    !********第二步、组装right_b,J矩阵*******
    CALL right()
    print*,'整体矩阵组装完毕！'
    print*
    print*,'正在模型降阶...'
     !******第三步、求M-1J，x_right*****
    call MJ_inv()
    print*,'模型降阶...'
    !*******第四步、求J_right，x,VTMT,VTM,求C矩阵*******
    
    !*******第五步、求T矩阵*******
   call model_reduction()
    print*,'模型降阶完毕！'
    print* 
    print*,'正在计算电场值...'
    call subspace()
  
    ENDSUBROUTINE 
    
    
subroutine  K8(Xx,Yy,Zz,sgm,Ke,ke1,ke2,kpe,qe,qpe,Me,me1,me2)
!    [ ke ,kpe, qe ]  
!    [kpet,ke1, qpe]
!    [qet, qpet,ke2]
!    组装矩阵 ，该矩阵为模型降阶算法左端项
!    组装M矩阵
!    [ Me , 0 , 0]  
!    [ 0 , Me1, 0]
!    [ 0 , 0 ,Me2]
!    (C+i*omg*M)E=b
!    模型降阶解方程为  (C-fx*M)E=M^-1 b
! Liu Jiren and Zhang Jifeng
    use net,only:sgmp,dof,ennum
    use GLOBAL_CONST
    use sub,only:fx
    implicit none
    real(8) KE(8,8),k11,ke1(8,8),ke2(8,8),kpe(8,8),qpe(8,8),qe(8,8),Me(8,8),Me1(8,8),Me2(8,8)
    REAL(8) Xx(8),Zz(8),Yy(8)
    REAL(8) sgm,omg,y0,z0,x0
    real(8) Nx(8,8),Ny(8,8),Nz(8,8),Nij(8,8),Ni(8)
    real(8) Nxy(8,8),Nxz(8,8),Nyz(8,8),k12
    real(8) a,b,c
    integer i,j,kk2,m
    
    a=abs(xx(2)-xx(1))
    b=abs(yy(3)-yy(2))
    c=abs(zz(5)-zz(1))
    ke=0.0;ke1=0.0;ke2=0.0;kpe=0.0;me=0.0;me1=0.0;me2=0.0;qpe=0.0;qe=0.0
   
    data Nx/ 4.0,-4.0,-2.0, 2.0, 2.0,-2.0,-1.0, 1.0,&
            -4.0, 4.0, 2.0,-2.0,-2.0, 2.0, 1.0,-1.0,&
            -2.0, 2.0, 4.0,-4.0,-1.0, 1.0, 2.0,-2.0,&
             2.0,-2.0,-4.0, 4.0, 1.0,-1.0,-2.0, 2.0,&
             2.0,-2.0,-1.0, 1.0, 4.0,-4.0,-2.0, 2.0,&
            -2.0, 2.0, 1.0,-1.0,-4.0, 4.0, 2.0,-2.0,&
            -1.0, 1.0, 2.0,-2.0,-2.0, 2.0, 4.0,-4.0,&
             1.0,-1.0,-2.0, 2.0, 2.0,-2.0,-4.0, 4.0/
    
    data Ny/4.0, 2.0,-2.0,-4.0, 2.0, 1.0,-1.0,-2.0 ,&
            2.0, 4.0,-4.0,-2.0, 1.0, 2.0,-2.0,-1.0 ,&
           -2.0,-4.0, 4.0, 2.0,-1.0,-2.0, 2.0, 1.0 ,&
           -4.0,-2.0, 2.0, 4.0,-2.0,-1.0, 1.0, 2.0 ,&
            2.0, 1.0,-1.0,-2.0, 4.0, 2.0,-2.0,-4.0 ,&
            1.0, 2.0,-2.0,-1.0, 2.0, 4.0,-4.0,-2.0 ,&
           -1.0,-2.0, 2.0, 1.0,-2.0,-4.0, 4.0, 2.0 ,&
           -2.0,-1.0, 1.0, 2.0,-4.0,-2.0, 2.0, 4.0/
    
    data Nz/4.0, 2.0, 1.0, 2.0,-4.0,-2.0,-1.0,-2.0,&
            2.0, 4.0, 2.0, 1.0,-2.0,-4.0,-2.0,-1.0,&
            1.0, 2.0, 4.0, 2.0,-1.0,-2.0,-4.0,-2.0,&
            2.0, 1.0, 2.0, 4.0,-2.0,-1.0,-2.0,-4.0,&
           -4.0,-2.0,-1.0,-2.0, 4.0, 2.0, 1.0, 2.0,&
           -2.0,-4.0,-2.0,-1.0, 2.0, 4.0, 2.0, 1.0,&
           -1.0,-2.0,-4.0,-2.0, 1.0, 2.0, 4.0, 2.0,&
           -2.0,-1.0,-2.0,-4.0, 2.0, 1.0, 2.0, 4.0/
    
    data Nij/8.0, 4.0, 2.0, 4.0, 4.0, 2.0, 1.0, 2.0,&
             4.0, 8.0, 4.0, 2.0, 2.0, 4.0, 2.0, 1.0,&
             2.0, 4.0, 8.0, 4.0, 1.0, 2.0, 4.0, 2.0,&
             4.0, 2.0, 4.0, 8.0, 2.0, 1.0, 2.0, 4.0,&
             4.0, 2.0, 1.0, 2.0, 8.0, 4.0, 2.0, 4.0,&
             2.0, 4.0, 2.0, 1.0, 4.0, 8.0, 4.0, 2.0,&
             1.0, 2.0, 4.0, 2.0, 2.0, 4.0, 8.0, 4.0,&
             2.0, 1.0, 2.0, 4.0, 4.0, 2.0, 4.0, 8.0/
    
    data Nxy/2.0, 2.0,-2.0,-2.0, 1.0, 1.0,-1.0,-1.0,&
            -2.0,-2.0, 2.0, 2.0,-1.0,-1.0, 1.0, 1.0,&
            -2.0,-2.0, 2.0, 2.0,-1.0,-1.0, 1.0, 1.0,&
             2.0, 2.0,-2.0,-2.0, 1.0, 1.0,-1.0,-1.0,&
             1.0, 1.0,-1.0,-1.0, 2.0, 2.0,-2.0,-2.0,&
            -1.0,-1.0, 1.0, 1.0,-2.0,-2.0, 2.0, 2.0,&
            -1.0,-1.0, 1.0, 1.0,-2.0,-2.0, 2.0, 2.0,&
             1.0, 1.0,-1.0,-1.0, 2.0, 2.0,-2.0,-2.0/
    data Nxz/2.0, 2.0, 1.0, 1.0,-2.0,-2.0,-1.0,-1.0,&
            -2.0,-2.0,-1.0,-1.0, 2.0, 2.0, 1.0, 1.0,&
            -1.0,-1.0,-2.0,-2.0, 1.0, 1.0, 2.0, 2.0,&
             1.0, 1.0, 2.0, 2.0,-1.0,-1.0,-2.0,-2.0,&
             2.0, 2.0, 1.0, 1.0,-2.0,-2.0,-1.0,-1.0,&
            -2.0,-2.0,-1.0,-1.0, 2.0, 2.0, 1.0, 1.0,&
            -1.0,-1.0,-2.0,-2.0, 1.0, 1.0, 2.0, 2.0,&
             1.0, 1.0, 2.0, 2.0,-1.0,-1.0,-2.0,-2.0 /
    
    data Nyz/2.0, 1.0, 1.0, 2.0,-2.0,-1.0,-1.0,-2.0,&
             1.0, 2.0, 2.0, 1.0,-1.0,-2.0,-2.0,-1.0,&
            -1.0,-2.0,-2.0,-1.0, 1.0, 2.0, 2.0, 1.0,&
            -2.0,-1.0,-1.0,-2.0, 2.0, 1.0, 1.0, 2.0,&
             2.0, 1.0, 1.0, 2.0,-2.0,-1.0,-1.0,-2.0,&
             1.0, 2.0, 2.0, 1.0,-1.0,-2.0,-2.0,-1.0,&
            -1.0,-2.0,-2.0,-1.0, 1.0, 2.0, 2.0, 1.0,&
            -2.0,-1.0,-1.0,-2.0, 2.0, 1.0, 1.0, 2.0/
    
    kk2=dcmplx(0.0,1.0)
     do i=1,ennum
       do j=1,ennum
  
        k11=fx*dble(miu*sgm*Nij(i,j)*(a*b*c/dble(216.0)))                                                       !fx*k2e..k2e=-miu*sgm*NIJ*(a*b*c/dble(216.0)
        
        KE(i,J)= dble(Nx(i,j)*(b*c/dble(36.0)/a)+Ny(i,j)*(a*c/dble(36.0)/b)+Nz(i,j)*(a*b/dble(36.0)/c))-k11 ! 模型降阶左端所用矩阵对角线矩阵位置 1,1 
        ke1(i,j)=dble(Nx(i,j)*(b*c/dble(36.0)/a)+Ny(i,j)*(a*c/dble(36.0)/b)+Nz(i,j)*(a*b/dble(36.0)/c))-k11 ! Ny(i,j)*(a*c/dble(36.0)/b)         2,2
        ke2(i,j)=dble(Nx(i,j)*(b*c/dble(36.0)/a)+Ny(i,j)*(a*c/dble(36.0)/b)+Nz(i,j)*(a*b/dble(36.0)/c))-k11 !   +Nz(i,j)*(a*b/dble(36.0)/c)      3,3
      
        KPE(I,J)=dble((c/dble(24.0))*(Nxy(i,j)-Nxy(j,i)))  !非对角线位置 1,2
        qe(i,j)=dble((b/dble(24.0))*(Nxz(i,j)-Nxz(j,i)))   !             1,3
        qpe(i,j)=dble((a/dble(24.0))*(Nyz(i,j)-Nyz(j,i)))   !             2,3
        
       end do                                            
     end do
     
     
     do i=1,ennum
         do j=1,ennum
            k11=dble(miu*sgm*Nij(i,j)*(a*b*c/dble(216.0)))
            Me(i,j)=k11                                     
            me1(i,j)=k11
            me2(i,j)=k11
         end do
     end do

    end subroutine
    
    
 SUBROUTINE GK()
 ! 总体合成,双向链表
 ! Liu Jiren and Zhang Jifeng
    USE NET,ONLY:I8,RHO,X,Y,Z,num_unit,num_point,dof,ennum
    USE LINKED_list,ONLY:Q,P,LK,LQ,LQP,Lp,lk1,Lk2,Le,le1,le2
    IMPLICIT NONE
    real(8) KE(8,8),ke1(8,8),ke2(8,8),kpe(8,8),qpe(8,8),qe(8,8),me(8,8),me1(8,8),me2(8,8),kee(8,8),kee2(8,8),kee1(8,8)
    REAL(8) Xx(8),Zz(8),Yy(8)
    REAL(8) sgm,omg
    INTEGER i1,j1,i,j,k,it,step,ii1,jj2
    !初始化
    DO k=1,num_point 
        !指向对角元素
        ALLOCATE(Q);Q.j=k;Q.c=0.0;NULLIFY(Q.next,Q.prev);LK.ht(k).head=>Q;LK.ht(k).tail=>Q
        ALLOCATE(Q);Q.j=k;Q.c=0.0;NULLIFY(Q.next,Q.prev);LK1.ht(k).head=>Q;LK1.ht(k).tail=>Q
        ALLOCATE(Q);Q.j=k;Q.c=0.0;NULLIFY(Q.next,Q.prev);LK2.ht(k).head=>Q;LK2.ht(k).tail=>Q
        !指向非对角矩阵的第一列
        ALLOCATE(Q);Q.j=1;Q.c=0.0;NULLIFY(Q.next,Q.prev);Lq.ht(k).head=>Q;Lq.ht(k).tail=>Q
        ALLOCATE(Q);Q.j=1;Q.c=0.0;NULLIFY(Q.next,Q.prev);LQp.ht(k).head=>Q;LQp.ht(k).tail=>Q
        ALLOCATE(Q);Q.j=1;Q.c=0.0;NULLIFY(Q.next,Q.prev);Lp.ht(k).head=>Q;Lp.ht(k).tail=>Q
     
        !指向对角元素
        ALLOCATE(Q);Q.j=k;Q.c=0.0;NULLIFY(Q.next,Q.prev);Le.ht(k).head=>Q;Le.ht(k).tail=>Q
        ALLOCATE(Q);Q.j=k;Q.c=0.0;NULLIFY(Q.next,Q.prev);Le1.ht(k).head=>Q;Le1.ht(k).tail=>Q
        ALLOCATE(Q);Q.j=k;Q.c=0.0;NULLIFY(Q.next,Q.prev);Le2.ht(k).head=>Q;Le2.ht(k).tail=>Q
 
    ENDDO
 LK.GAi=1;LK1.GAi=1;Lk2.GAi=1;LQp.GAi=1;LP.GAi=1;Lq.GAi=1;le.gai=1;le1.gai=1;le2.gai=1;
    step=num_unit/10

 DO k=1,num_unit
        if(mod(k,step)==0) print 20,int(real(k)/num_unit*95)
20      format(i3,'% Completed...')
        !print*,"xx"
        Xx=(/X(I8(1,k)),X(I8(2,k)),X(I8(3,k)),X(I8(4,k)),X(I8(5,k)),X(I8(6,k)),X(I8(7,k)),X(I8(8,k))/)
        Zz=(/Z(I8(1,k)),Z(I8(2,k)),Z(I8(3,k)),Z(I8(4,k)),Z(I8(5,k)),Z(I8(6,k)),Z(I8(7,k)),Z(I8(8,k))/)
        Yy=(/Y(I8(1,k)),Y(I8(2,k)),Y(I8(3,k)),Y(I8(4,k)),Y(I8(5,k)),Y(I8(6,k)),Y(I8(7,k)),Y(I8(8,k))/)
        sgm=1.0/RHO(k)
     
CALL K8(Xx,Yy,Zz,sgm,KE,ke1,ke2,kpe,qe,qpe,me,me1,me2)
   
          DO i1=1,8
                 DO j1=i1,8
                i=I8(i1,k);j=I8(j1,k)
               IF(i>j) THEN
                    it=i;i=j;j=it
              ENDIF
                !**********************************************************
                IF(( ABS(KE(i1,j1)) )>1E-45) THEN
                    !为链表LK添加元素
                    P=>LK.ht(i).tail
                    DO 
                        IF(j==P.j) THEN
                            P.c=P.c+KE(i1,j1);EXIT
                        ELSEIF(j<P.j) THEN
                            P=>P.prev
                        ELSEIF(j>P.j) THEN !后方插入节点
                            ALLOCATE(Q);Q.j=j;Q.c=KE(i1,j1);LK.GAi(i)=LK.GAi(i)+1
                            Q.prev=>P;Q.next=>P.next
                            IF(ASSOCIATED(P.next)) THEN
                                P.next.prev=>Q
                            ELSE
                                LK.ht(i).tail=>Q
                            ENDIF
                            P.next=>Q;EXIT
                        ENDIF
                    ENDDO
                endif
               !**********************************************************
                IF(( ABS(KE1(i1,j1)) )>1E-45) THEN
                       P=>LK1.ht(i).tail
                    DO 
                        !print*,ke
                        IF(j==P.j) THEN
                            P.c=P.c+KE1(i1,j1);EXIT
                        ELSEIF(j<P.j) THEN
                            P=>P.prev
                        ELSEIF(j>P.j) THEN !后方插入节点
                            ALLOCATE(Q);Q.j=j;Q.c=KE1(i1,j1);LK1.GAi(i)=LK1.GAi(i)+1
                            Q.prev=>P;Q.next=>P.next
                            IF(ASSOCIATED(P.next)) THEN
                                P.next.prev=>Q
                            ELSE
                                LK1.ht(i).tail=>Q
                            ENDIF
                            P.next=>Q;EXIT
                        ENDIF
                    ENDDO
                endif
                  !**********************************************************
                   IF(( ABS(KE2(i1,j1)) )>1E-45) THEN
                        P=>LK2.ht(i).tail
                      
                    DO 
                        IF(j==P.j) THEN
                            P.c=P.c+KE2(i1,j1);EXIT
                        ELSEIF(j<P.j) THEN
                            P=>P.prev;
                        ELSEIF(j>P.j) THEN !后方插入节点
                            ALLOCATE(Q);Q.j=j;Q.c=KE2(i1,j1);LK2.GAi(i)=LK2.GAi(i)+1
                            Q.prev=>P;Q.next=>P.next
                            IF(ASSOCIATED(P.next)) THEN
                                P.next.prev=>Q
                            ELSE
                                LK2.ht(i).tail=>Q
                            ENDIF
                            P.next=>Q;EXIT
                        ENDIF
                    ENDDO
                   end if
                     !**********************************************************
                    IF(( ABS(me(i1,j1)) )>1E-45) THEN
                        P=>Le.ht(i).tail
                      
                    DO 
                        IF(j==P.j) THEN
                            P.c=P.c+me(i1,j1);EXIT
                        ELSEIF(j<P.j) THEN
                            P=>P.prev;
                        ELSEIF(j>P.j) THEN !后方插入节点
                            ALLOCATE(Q);Q.j=j;Q.c=me(i1,j1);Le.GAi(i)=Le.GAi(i)+1
                            Q.prev=>P;Q.next=>P.next
                            IF(ASSOCIATED(P.next)) THEN
                                P.next.prev=>Q
                            ELSE
                                Le.ht(i).tail=>Q
                            ENDIF
                            P.next=>Q;EXIT
                        ENDIF
                    ENDDO
                    end if
                      !**********************************************************
                      IF(( ABS(me1(i1,j1)) )>1E-45) THEN
                        P=>Le1.ht(i).tail
                      
                    DO 
                        IF(j==P.j) THEN
                            P.c=P.c+me1(i1,j1);EXIT
                        ELSEIF(j<P.j) THEN
                            P=>P.prev;
                        ELSEIF(j>P.j) THEN !后方插入节点
                            ALLOCATE(Q);Q.j=j;Q.c=me1(i1,j1);Le1.GAi(i)=Le1.GAi(i)+1
                            Q.prev=>P;Q.next=>P.next
                            IF(ASSOCIATED(P.next)) THEN
                                P.next.prev=>Q
                            ELSE
                                Le1.ht(i).tail=>Q
                            ENDIF
                            P.next=>Q;EXIT
                        ENDIF
                    ENDDO
                      end if
                      !**********************************************************
                         IF(( ABS(me2(i1,j1)) )>1E-45) THEN
                        P=>Le2.ht(i).tail
                      
                    DO 
                        IF(j==P.j) THEN
                            P.c=P.c+me2(i1,j1);EXIT
                        ELSEIF(j<P.j) THEN
                            P=>P.prev;
                        ELSEIF(j>P.j) THEN !后方插入节点
                            ALLOCATE(Q);Q.j=j;Q.c=me2(i1,j1);Le2.GAi(i)=Le2.GAi(i)+1
                            Q.prev=>P;Q.next=>P.next
                            IF(ASSOCIATED(P.next)) THEN
                                P.next.prev=>Q
                            ELSE
                                Le2.ht(i).tail=>Q
                            ENDIF
                            P.next=>Q;EXIT
                        ENDIF
                    ENDDO
                         end if
                 end do
          end do
          
          
           DO i1=1,8
                 DO j1=1,8
                i=I8(i1,k);j=I8(j1,k) 
               !**********************************************************
                  IF(abs( (KPE(i1,j1)) )>1E-45) THEN  
                    P=>Lp.ht(i).tail
                    DO 
                        IF(j==P.j) THEN
                            P.c=P.c+Kpe(i1,j1);EXIT
                        ELSEIF(j<P.j) THEN
                            P=>P.prev
                        ELSEIF(j>P.j) THEN !后方插入节点
                            ALLOCATE(Q);Q.j=j;Q.c=KpE(i1,j1);Lp.GAi(i)=Lp.GAi(i)+1
                            Q.prev=>P;Q.next=>P.next
                            IF(ASSOCIATED(P.next)) THEN
                                P.next.prev=>Q
                            ELSE
                                Lp.ht(i).tail=>Q
                            ENDIF
                            P.next=>Q;EXIT
                        ENDIF
                    ENDDO
                  end if
             !********************************************************** 
               IF(abs( (qe(i1,j1)) )>1E-45) THEN  
                    P=>Lq.ht(i).tail
                    DO 
                        IF(j==P.j) THEN
                            P.c=P.c+qe(i1,j1);EXIT
                        ELSEIF(j<P.j) THEN
                            P=>P.prev
                        ELSEIF(j>P.j) THEN !后方插入节点
                            ALLOCATE(Q);Q.j=j;Q.c=qe(i1,j1);Lq.GAi(i)=Lq.GAi(i)+1
                            Q.prev=>P;Q.next=>P.next
                            IF(ASSOCIATED(P.next)) THEN
                                P.next.prev=>Q
                            ELSE
                                Lq.ht(i).tail=>Q
                            ENDIF
                            P.next=>Q;EXIT
                        ENDIF
                            
                    ENDDO
               end if
                 !**********************************************************
                     IF(abs( (qpe(i1,j1)) )>1E-45) THEN  
                    P=>Lqp.ht(i).tail
                    DO 
                        IF(j==P.j) THEN
                            P.c=P.c+qpe(i1,j1);EXIT
                        ELSEIF(j<P.j) THEN
                            P=>P.prev
                        ELSEIF(j>P.j) THEN !后方插入节点
                            ALLOCATE(Q);Q.j=j;Q.c=qpe(i1,j1);Lqp.GAi(i)=Lqp.GAi(i)+1
                            Q.prev=>P;Q.next=>P.next
                            IF(ASSOCIATED(P.next)) THEN
                                P.next.prev=>Q
                            ELSE
                                Lqp.ht(i).tail=>Q
                            ENDIF
                            P.next=>Q;EXIT
                        ENDIF
                    ENDDO
                     endif
                     
             ENDDO 
        ENDDO 
          end do
    ENDSUBROUTINE
    
SUBROUTINE CONVERT()
! 总体合成,释放链表
! Liu Jiren and Zhang Jifeng
    USE NET,ONLY:num_point
    USE LINKED_list,ONLY:Q,P,LK,LQ,LQP,Lp,lk1,Lk2,le1,le2,le
    USE SparseMatrix,ONLY:A,IA,JA,M,IM,JM
    IMPLICIT NONE
    INTEGER i,j,ii,iii
    IA(1)=1
    IM(1)=1
   
     DO i=1,num_point
        IA(i+1)=IA(i)+LK.GAi(i)+LP.Gai(i)+LQ.GAI(I)
        IM(i+1)=Im(i)+LE.GAi(i)
     ENDDO
     
     DO i=num_point+1,2*NUM_POINT
        IA(i+1)=IA(i)+LK1.GAi(i-NUM_POINT)+LQP.Gai(i-NUM_POINT)
        IM(i+1)=IM(i)+Le1.GAi(i-num_point)
     ENDDO
     
      DO i=2*NUM_POINT+1,3*NUM_POINT
        IA(i+1)=IA(i)+LK2.GAi(i-2*NUM_POINT)
        IM(i+1)=IM(i)+LE2.GAi(i-2*num_point)
      ENDDO
 
    i=IA(3*num_point+1)-1
    ii=IM(3*num_point+1)-1
   
    ALLOCATE(A(i),JA(i))
    allocate(M(ii),JM(ii))
 
    DO i=1,num_point
        P=>LK.ht(i).head
        DO j=IA(i),IA(i)+LK.GAi(i)-1
            JA(j)=P.j;A(J)=P.c;Q=>P.next;deALLOCATE(P);P=>Q
        ENDDO
        
        P=>LP.ht(i).head
        DO j=IA(i)+LK.GAi(i),IA(i)+LK.GAi(i)+LP.GAI(i)-1
            JA(j)=P.j+num_point;A(J)=P.c;Q=>P.next;deALLOCATE(P);P=>Q
        ENDDO
        
          P=>LQ.ht(i).head
        DO j=IA(i)+LK.GAi(i)+LP.GAI(i),IA(I+1)-1
            JA(j)=P.j+2*num_point;A(J)=P.c;Q=>P.next;deALLOCATE(P);P=>Q
        ENDDO
        
         P=>LK1.ht(i).head
        DO j=IA(i+num_point),IA(i+NUM_POINT)+LK1.GAi(i)-1
            JA(j)=P.j+num_point;A(J)=P.c;Q=>P.next;deALLOCATE(P);P=>Q
        ENDDO
        
          P=>LQP.ht(i).head
        DO j=IA(i+NUM_POINT)+LK1.GAi(i),IA(i+num_point+1)-1
            JA(j)=P.j+2*num_point;A(J)=P.c;Q=>P.next;deALLOCATE(P);P=>Q
        ENDDO
        
            P=>LK2.ht(i).head
        DO j=IA(i+2*NUM_POINT),IA(i+2*num_point+1)-1
            JA(j)=P.j+2*num_point;A(J)=P.c;Q=>P.next;deALLOCATE(P);P=>Q
        ENDDO
    
         P=>Le.ht(i).head
        DO j=Im(i),Im(i)+Le.GAi(i)-1
            Jm(j)=P.j;m(j)=P.c;Q=>P.next;deALLOCATE(P);P=>Q
        ENDDO
        
         P=>Le1.ht(i).head
        DO j=Im(i+num_point),Im(i+NUM_POINT)+Le1.GAi(i)-1
            Jm(j)=P.j+num_point;m(j)=P.c;Q=>P.next;deALLOCATE(P);P=>Q
        ENDDO
        
         P=>Le2.ht(i).head
        DO j=Im(i+2*NUM_POINT),Im(i+2*num_point+1)-1
            Jm(j)=P.j+2*num_point;m(j)=P.c;Q=>P.next;deALLOCATE(P);P=>Q
        ENDDO
    ENDDO

    deallocate(LK.ht,LQ.ht,LQP.ht,Lp.ht,lk1.ht,Lk2.ht,le1.ht,le2.ht,le.ht)
    deallocate(LK.gai,LQ.gai,LQP.gai,Lp.gai,lk1.gai,Lk2.gai,le1.gai,le2.gai,le.gai)
   
    ENDSUBROUTINE
  
   !******第三步、求模型降阶的右端项*****
subroutine mj_inv()
!    计算M逆*right_b=x_right
!    相当于求解M*（M^-1 *right_b）=right_b；
!    模型降阶解方程为  (C-fx*M)E=M^-1 right_b
! Liu Jiren and Zhang Jifeng

    use Space_Domain_Field,only:right_b,X_right
    use SparseMatrix,only:M,IM,JM
    use net,only:num_point
     IMPLICIT NONE
    INTEGER lda, ldb, maxfct, mnum, mtype, phase, nrhs, error, msglvl
    INTEGER  iparm(64),n,i
    INTEGER*8 pt(64)
    integer,allocatable::perm(:)
    integer,external:: mkl_get_max_threads
    REAL*8 ddum
    lda=Im(3*num_point+1)-1;ldb=3*num_point
    allocate( perm(ldb) )
    pt=0 !指针初始化
    maxfct=1;mnum=1;mtype=-2;phase=13;n=ldb;perm=0;nrhs=1;x_right=dble(0.0)
    iparm(1) = 0
    error = 0 ! 错误信息初值置0 
    msglvl = 0 ! 不显示统计信息 
    
    call pardiso (pt, maxfct, mnum, mtype, phase, n, m, im, jm, perm, nrhs, iparm, msglvl, right_b, x_right, error)
    if(error/=0) print*,'失败'
    phase=-1
    call pardiso (pt, maxfct, mnum, mtype, phase, n, m, im, jm, perm, nrhs, iparm, msglvl, right_b, x_right, error)

    deallocate( perm )
    end 
    
!*****************开始模型降阶**********************
subroutine model_reduction()
!    利用单个重复极点fx
!     m次回代
!    采用正交化AD过程的变体形式构建
!    模型降阶解方程为  (C-fx*M)E=M^-1 b =>(C-fx*M)存储在A，IA,JA中
!    Liu Jiren and Zhang Jifeng
    use Space_Domain_Field,only:X_right,J_right,right_b
    use SparseMatrix,only:M,IM,JM,A,IA,JA,C,JC,IC
    use sub,only:mi,v,vtm,norm,T,vt
    use net,only:num_point
     IMPLICIT NONE
    INTEGER lda, ldb, maxfct, mnum, mtype, phase, nrhs, error, msglvl
     real(8) vi(3*num_point),vj(3*num_point),y(3*num_point),vtcvI(mi+1),mvj(3*num_point),mvj1(3*num_point),mvj2(3*num_point)
     real(8) y1(3*num_point),y2(3*num_point)
    INTEGER  iparm(64),n
    INTEGER*8 pt(64),ierr
    integer,allocatable::perm(:)
    integer,external:: mkl_get_max_threads
    REAL*8 ddum,x(3*num_point)
    real normm,dot
    integer,allocatable:: iw(:)
    integer i,j,p
   
    lda=IA(3*num_point+1)-1;ldb=3*num_point
     allocate( perm(ldb),iw(ldb) )
    pt=0 !指针初始化
    maxfct=1;mnum=1;mtype=-2;phase=13;n=ldb;perm=0;nrhs=1;x=dble(0.0)
    iparm(1) = 0
    error = 0 ! 错误信息初值置0 
    msglvl = 0 ! 不显示统计信息 
    mvj1=dble(0.0);mvj2=dble(0.0);y1=dble(0.0);y2=dble(0.0)
  
    norm=normm(x_right,3*num_point)               ! 求||b||m，范数的变体形式
    ! 正交化
    x_right(:)=x_right(:)/norm                    ! 求V1
     V(:,1)=x_right(:)                            ! 赋给v1
     call atmux ( n, x_right, mvj1, m, jm, im) 
     call amux ( n, x_right, mvj2,m, jm, im )      ! 求第一个右端项 M*V1，赋给mvj
     mvj=mvj1+mvj2
     
     do j=1,mi,1  
        if(j==1)then
                                                  ! pardiso 求解，存储在X中
           call pardiso (pt, maxfct, mnum, mtype, phase, n, a, ia, ja, perm, nrhs, iparm, msglvl, mvj, x, error)
           if(error/=0) print*,'失败'
        else                                      ! m次回代
             call pardiso (pt, maxfct, mnum, mtype,33, n, a, ia, ja, perm, nrhs, iparm, msglvl, mvj, x, error)
             if(error/=0) print*,'失败'
        end if

        do i=1,j
            vi(:)=v(:,i)
            x(:)=x(:)-dot(x,vi,3*num_point)*Vi(:)  ! 求解 x=x-(x,vi)*vi
        end do
           V(:,j+1)=x(:)/normm(x,3*num_point)      ! 求解 x/||x||m

            vi(:)= V(:,j+1)
           mvj1=0.0;mvj2=0.0
         call atmux ( n, vi, mvj1, m, jm, im)       !求下三角与其乘积
         call amux ( n, vi, mvj2,m, jm, im )        !求上三角与其乘积
         mvj=mvj1+mvj2                              ! 计算M*Vj+1
      
       end do
      phase=-1                                      ! 释放
         call pardiso (pt, maxfct, mnum, mtype, phase, n, a, ia, ja, perm, nrhs, iparm, msglvl,mVj, x, error)
   print*,"正在计算VTM"

    vt=transpose(V)  
    call get_c();                                   ! 得到c矩阵，C,IC，JC
     do i=1,mi+1
        vi(:)=v(:,i)    
       ! print*, vi!v的第i列       
        call atmux(3*num_point, vi, y1, c, jc, ic)
        call amux (3*num_point, vi, y2, c, jc, ic)
        y(:)=y1(:)+y2(:)

        vtcvi=matmul(vt,y)                           !得到T的第i列数据      
        T(:,i)=vtcVI(:)                              !得打T矩阵
     end do
     
      deallocate( perm,vtm,A,IA,JA,M,JM,IM,C,IC,JC)
    end subroutine
    
  
!*******************得到C矩阵************************
subroutine get_c()
    use sub,only:fx
    use SparseMatrix,only:M,IM,JM,A,IA,JA,C,JC,IC
    use net,only:num_point
    implicit none
    integer nrow,ncol,job,nzmax,ierr,i,nm
    integer,allocatable:: iw(:)
    ierr=0
    nzmax=Ia(3*num_point+1)-1;nm=Im(3*num_point+1)-1
    nrow=3*num_point;ncol=3*num_point
   allocate(jc(nzmax),ic(nrow+1),c(nzmax),iw(ncol))

   M(:)=fx*M(:)
    call aplb(nrow,ncol,1,a,ja,ia,m,jm,im,C,jc,ic,nzmax,iw,ierr)
    if(ierr>0)then
        print*,ierr,"转换C矩阵失败"
    end if
    M(:)=M(:)/fx     !恢复M

    end subroutine
    
    
   !*******************计算H场************************ 
      SUBROUTINE Cal_Assist_Component()
    USE Space_Domain_Field,ONLY:Ex_z,Ez_x,Exyz,Hy,Ex_y,Ey_x,hz,Hx,Ez_y,Ey_z
    USE GLOBAL_CONST,ONLY:epsilon,miu,pi
    USE NET,ONLY:num_point,num_unit,I8,X,RHO,Y,Z
    use FREQUENCY
    IMPLICIT NONE
    REAL(8) omg,sgm
    REAL(8) xi_x,gama_z,eta_y!全局坐标对局部坐标的导数
    REAL(8),PARAMETER::xi(8)=(/-1.0,1.0,1.0,-1.0,-1.0,1.0,1.0,-1.0/),eta(8)=(/-1.0,-1.0,1.0,1.0,-1.0,-1.0,1.0,1.0/)
    REAL(8),PARAMETER::gama(8)=(/-1.0,-1.0,-1.0,-1.0,1.0,1.0,1.0,1.0/)
    REAL(8) Xx(8),Zz(8),yy(8),a,b,c
    REAL(8) temp1,temp2,temp3,temp4,temp5,temp6
    REAL(8) N_xi(8),N_eta(8),N_gama(8) !基函数对自然坐标的导数
    COMPLEX(8) zc,yc,k2,ke2,iky,ky2
    COMPLEX(8) ctemp1,ctemp2,ctemp3,ctemp4,ctemp5,ctemp6
    INTEGER(2),ALLOCATABLE::NU(:)
    INTEGER::i,j,k,ii,jj,numf
    !计算Ex、Ez沿z、x的导数
    ALLOCATE(NU(num_point))
 
    
   do numf=1,num_freq
       omg=2*pi*freq(numf)
       Ez_x=0.0;Ex_z=0.0;Ex_y=0.0;Ey_x=0.0;Ez_y=0.0;Ey_z=0.0;nu=0.0
    DO k=1,num_unit
     
        Xx=(/X(I8(1,k)),X(I8(2,k)),X(I8(3,k)),X(I8(4,k)),X(I8(5,k)),X(I8(6,k)),X(I8(7,k)),X(I8(8,k))/)
        Zz=(/Z(I8(1,k)),Z(I8(2,k)),Z(I8(3,k)),Z(I8(4,k)),Z(I8(5,k)),Z(I8(6,k)),Z(I8(7,k)),Z(I8(8,k))/)
        Yy=(/Y(I8(1,k)),Y(I8(2,k)),Y(I8(3,k)),Y(I8(4,k)),Y(I8(5,k)),Y(I8(6,k)),Y(I8(7,k)),Y(I8(8,k))/)
        
        sgm=1/RHO(k)
    if(zz(1)<0)then
    a=abs(xx(2)-xx(1))
    b=abs(yy(3)-yy(2))
    c=abs(zz(5)-zz(1))
    
        DO j=1,8
           temp1=1.0-xi(j);temp2=1.0-eta(j);temp3=1.0+xi(j);temp4=1.0+eta(j)
           temp5=1.0-gama(j);temp6=1.0+gama(j);
                 N_ETA(1)=1.0/(8)*(-1)*temp1*temp5
                N_ETA(2)=1.0/(8)*(-1)*temp3*temp5
                N_ETA(3)=1.0/(8)*(1)*temp3*temp5
                N_ETA(4)=1.0/(8)*(1)*temp1*temp5
                N_ETA(5)=1.0/(8)*(-1)*temp1*temp6
                N_ETA(6)=1.0/(8)*(-1)*temp3*temp6
                N_ETA(7)=1.0/(8)*(1)*temp3*temp6
                N_ETA(8)=1.0/(8)*(1)*temp1*temp6
                
                 N_XI(1)=1.0/(8)*(-1)*temp2*temp5
                N_XI(2)=1.0/(8)*(1)*temp2*temp5
                N_XI(3)=1.0/(8)*(1)*temp4*temp5
                N_XI(4)=1.0/(8)*(-1)*temp4*temp5
                N_XI(5)=1.0/(8)*(-1)*temp2*temp6
                N_XI(6)=1.0/(8)*(1)*temp2*temp6
                N_XI(7)=1.0/(8)*(1)*temp4*temp6
                N_XI(8)=1.0/(8)*(-1)*temp4*temp6
                
                N_GAMA(1)=1.0/(8)*(-1)*temp1*temp2
                N_gama(2)=1.0/(8)*(-1)*temp3*temp2
                N_gama(3)=1.0/(8)*(-1)*temp3*temp4
                N_gama(4)=1.0/(8)*(-1)*temp1*temp4
                N_gama(5)=1.0/(8)*(1)*temp1*temp2
                N_gama(6)=1.0/(8)*(1)*temp3*temp2
                N_gama(7)=1.0/(8)*(1)*temp3*temp4
                N_gama(8)=1.0/(8)*(1)*temp1*temp4
      
                xi_x=2.0/a;gama_z=2.0/c;eta_y=2.0/b
         
            ctemp1=0q0;ctemp2=0q0;ctemp3=0q0;ctemp4=0q0;CTEMP5=0q0;ctemp6=0q0
            DO i=1,8
                ii=I8(i,k)
          !      ctemp1=ctemp1+N_xi(i)*xi_x *(Exyz(numf,ii+2*num_point))   !Ez to x
           !     ctemp2=ctemp2+N_gama(i)*gama_z*(Exyz(numf,ii));           !ex to z
           !     ctemp3=ctemp3+N_xi(i)*xi_x *(Exyz(numf,ii+num_point))     !eY  to x
            !    ctemp4=ctemp4+N_ETA(i)*ETA_Y *(Exyz(numf,ii))             !Ex to y
                 
            !    ctemp5=ctemp5+N_eta(i)*eta_y *(Exyz(numf,ii+2*num_point))                 !ez  to y
             !   ctemp6=ctemp6+N_gama(i)*gama_z *(Exyz(numf,ii+num_point))                 !Ey to z
            ENDDO
            
            jj=I8(j,k)
            Ez_x(jj)=Ez_x(jj)+ctemp1
            Ex_z(jj)=Ex_z(jj)+ctemp2
            
            Ex_y(jj)=Ex_y(jj)+ctemp4
            Ey_x(jj)=Ey_x(jj)+ctemp3
            
            Ez_y(jj)=Ez_y(jj)+ctemp5
            Ey_z(jj)=Ey_z(jj)+ctemp6
            
            NU(jj)=NU(jj)+1
        ENDDO
    end if
    ENDDO
    DO k=1,num_point
           Ex_z(k)=Ex_z(k)/NU(k)!取平均
           Ez_x(k)=Ez_x(k)/NU(k)
       !    Hy(numf,k)=-dcmplx(0.0,1.0)/(omg*miu)*(Ex_z(k)-Ez_x(k))
           
              
           Ex_y(k)=Ex_y(k)/NU(k)!取平均
           Ey_x(k)=Ey_x(k)/NU(k)
       !     Hz(k)=-dcmplx(0,1)/(omg*miu)*(Ey_x(k)-Ex_y(k))
            
             Ez_y(k)=Ez_y(k)/NU(k)!取平均
             Ey_z(k)=Ey_z(k)/NU(k)
        !    Hx(numf,k)=-dcmplx(0.0,1.0)/(omg*miu)*(Ez_y(k)-Ey_z(k))
    ENDDO
    end do
    deALLOCATE(NU)
    ENDSUBROUTINE
   
subroutine CSr_store(Gt,n)
use SparseMatrix,only:G,Ig,Jg,GG,num
use LINKED_list,only:Lg,Q,P
use FREQUENCY,only: num_freq
implicit none
 integer i,j,k,n,i1,j1
  complex(8) Gt(n,n)
  ALLOCATE(Lg.ht(n),Lg.GAi(n))

  
    DO k=1,n!指向第一列元素
    ALLOCATE(Q);Q.j=1;Q.cc=0.0;NULLIFY(Q.next,Q.prev);Lg.ht(k).head=>Q;Lg.ht(k).tail=>Q
    end do
    Lg.Gai=1
    do i1=1,n
        do j1=1,n
            i=i1;j=j1
            IF(( ABS(Gt(i1,j1)) )>1E-45) THEN
                    !为链表LK添加元素
                    P=>Lg.ht(i).tail
                    DO 
                        !print*,ke
                        IF(j==P.j) THEN
                            P.cc=P.cc+Gt(i1,j1);EXIT
                        ELSEIF(j<P.j) THEN
                            P=>P.prev
                        ELSEIF(j>P.j) THEN !后方插入节点
                            ALLOCATE(Q);Q.j=j;Q.cc=Gt(i1,j1);Lg.GAi(i)=Lg.GAi(i)+1
                            Q.prev=>P;Q.next=>P.next
                            IF(ASSOCIATED(P.next)) THEN
                                P.next.prev=>Q
                            ELSE
                                Lg.ht(i).tail=>Q
                            ENDIF
                            P.next=>Q;EXIT
                        ENDIF
                    ENDDO
            endif  
        end do
    end do
    Ig(1)=1
     DO i=1,n
        Ig(i+1)=Ig(i)+Lg.GAi(i)
     ENDDO
         i=Ig(n+1)-1
         num=i
    ALLOCATE(g(num_freq,i),Jg(i))
    
    DO i=1,n
        P=>Lg.ht(i).head
        DO j=Ig(i),Ig(i)+Lg.GAi(i)-1
            Jg(j)=P.j;g(1:num_freq,j)=P.cc;Q=>P.next;deALLOCATE(P);P=>Q
        ENDDO
    end do
    
deallocate(LG.ht,lg.gai)
    end subroutine
    

