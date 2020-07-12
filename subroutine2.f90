  subroutine aplb ( nrow, ncol, job, a, ja, ia, b, jb, ib, c, jc, ic, nzmax, &
  iw, ierr )

!*****************************************************************************80
!
!! APLB performs the CSR matrix sum C = A + B.
!
!  Modified:
!
!    07 January 2004
!
!  Author:
!
!    Youcef Saad
!
!  Parameters:
!
!    Input, integer ( kind = 4 ) NROW, the row dimension of A and B.
!
!    Input, integer ( kind = 4 ) NCOL, the column dimension of A and B.
!
!    Input, integer ( kind = 4 ) JOB.  When JOB = 0, only the structure
!    (i.e. the arrays jc, ic) is computed and the
!    real values are ignored.
!
!    Input, real A(*), integer ( kind = 4 ) JA(*), IA(NROW+1), the matrix in CSR
!    Compressed Sparse Row format.
!
! b,
! jb,
! ib      =  Matrix B in compressed sparse row format.
!
! nzmax      = integer ( kind = 4 ). The  length of the arrays c and jc.
!         amub will stop if the result matrix C  has a number
!         of elements that exceeds exceeds nzmax. See ierr.
!
! on return:
!
! c,
! jc,
! ic      = resulting matrix C in compressed sparse row sparse format.
!
! ierr      = integer ( kind = 4 ). serving as error message.
!         ierr = 0 means normal return,
!         ierr > 0 means that amub stopped while computing the
!         i-th row  of C with i = ierr, because the number
!         of elements in C exceeds nzmax.
!
! work arrays:
!
! iw      = integer ( kind = 4 ) work array of length equal to the number of
!         columns in A.
!
  implicit none

  integer ( kind = 4 ) ncol
  integer ( kind = 4 ) nrow

  real ( kind = 8 ) a(*)
  REAL ( kind = 8 ) b(*)
  REAL ( kind = 8 ) c(*)
  integer ( kind = 4 ) ia(nrow+1)
  integer ( kind = 4 ) ib(nrow+1)
  integer ( kind = 4 ) ic(nrow+1)
  integer ( kind = 4 ) ierr
  integer ( kind = 4 ) ii
  integer ( kind = 4 ) iw(ncol)
  integer ( kind = 4 ) ja(*)
  integer ( kind = 4 ) jb(*)
  integer ( kind = 4 ) jc(*)
  integer ( kind = 4 ) jcol
  integer ( kind = 4 ) job
  integer ( kind = 4 ) jpos
  integer ( kind = 4 ) k
  integer ( kind = 4 ) ka
  integer ( kind = 4 ) kb
  integer ( kind = 4 ) len
  integer ( kind = 4 ) nzmax
  logical values

  values = ( job /= 0 )
  ierr = 0
  len = 0
  ic(1) = 1
  iw(1:ncol) = 0

  do ii = 1, nrow
!
!  Row I.
!
     do ka = ia(ii), ia(ii+1)-1

        len = len + 1
        jcol = ja(ka)

        if ( nzmax < len ) then
          ierr = ii
          return
        end if

        jc(len) = jcol
        if ( values ) then
          c(len) = a(ka)
        end if
        iw(jcol) = len
     end do

     do kb = ib(ii), ib(ii+1)-1

        jcol = jb(kb)
        jpos = iw(jcol)

        if ( jpos == 0 ) then

           len = len + 1

           if ( nzmax < len ) then
             ierr = ii
             return
           end if

           jc(len) = jcol
           if ( values ) then
             c(len) = b(kb)
           end if
           iw(jcol)= len
        else
           if ( values ) then
             c(jpos) = c(jpos) + b(kb)
           end if
        end if

     end do

     do k = ic(ii), len
       iw(jc(k)) = 0
     end do

     ic(ii+1) = len+1
  end do

  return
    end
   

    
   
 
    
     subroutine amux ( n, x, y, a, ja, ia )

!*****************************************************************************80
!
!! AMUX multiplies a CSR matrix A times a vector.
!
!  Discussion:
!
!    This routine multiplies a matrix by a vector using the dot product form.
!    Matrix A is stored in compressed sparse row storage.
!
!  Modified:
!
!    07 January 2004
!
!  Author:
!
!    Youcef Saad
!
!  Parameters:
!
!    Input, integer ( kind = 4 ) N, the row dimension of the matrix.
!
!    Input, real X(*), and array of length equal to the column dimension 
!    of A.
!
!    Input, real A(*), integer ( kind = 4 ) JA(*), IA(NROW+1), the matrix in CSR
!    Compressed Sparse Row format.
!
!    Output, real Y(N), the product A * X.
!
  implicit none

  integer ( kind = 4 ) n

  real ( kind = 8 ) a(*)
  integer ( kind = 4 ) i,j
  integer ( kind = 4 ) ia(*)
  integer ( kind = 4 ) ja(*)
  integer ( kind = 4 ) k
  real ( kind = 8 ) t
  real ( kind = 8 ) x(*)
  real ( kind = 8 ) y(n)

  do i = 1, n
!
!  Compute the inner product of row I with vector X.
!
    t = 0.0D+00
  
    do k = ia(i), ia(i+1)-1
       ! if(i==ja(k))then
          !  a(k)=a(k)-1256.63
           !  t = t + a(k) * x(ja(k))
           !  else
   if(ja(k)==i)then
       t=t
       else
               
            t = t + a(k) * x(ja(k))
       end if
       
      
    end do
   
    y(i) = t

  end do

  return
    end

   
   
   subroutine atmux ( n, x, y, a, ja, ia )

!*****************************************************************************80
!
!! ATMUX computes A' * x for a CSR matrix A.
!
!  Discussion:
!
!    This routine multiplies the transpose of a matrix by a vector when the
!    original matrix is stored in compressed sparse row storage. Can also be
!    viewed as the product of a matrix by a vector when the original
!    matrix is stored in the compressed sparse column format.
!
!  Modified:
!
!    07 January 2004
!
!  Author:
!
!    Youcef Saad
!
!  Parameters:
!
!    Input, integer ( kind = 4 ) N, the row dimension of the matrix.
!
!    Input, real X(*), an array whose length is equal to the
!    column dimension of A.
!
!    Output, real Y(N), the product A' * X.
!
!    Input, real A(*), integer ( kind = 4 ) JA(*), IA(N+1), the matrix in CSR
!    Compressed Sparse Row format.
!
  implicit none

  integer ( kind = 4 ) n

  real ( kind = 8 ) a(*)
  integer ( kind = 4 ) i
  integer ( kind = 4 ) ia(*)
  integer ( kind = 4 ) ja(*)
  integer ( kind = 4 ) k
  real ( kind = 8 ) x(*)
  real ( kind = 8 ) y(n)

  y(1:n) = 0.0D+00

  do i = 1, n
    do k = ia(i), ia(i+1)-1
      y(ja(k)) = y(ja(k)) + x(i) * a(k)
    end do
  end do

  return
    end 
    
  
!    
!!*******************求逆************************          
subroutine CSR_inv(Gt,n,Gt_inv)
   !use LINKED_list,only:Lgt,Q,P
    implicit none
    integer i,j,k,n,i1,j1
    complex(8) Gt(n,n),GT_inv(n,n)
     INTEGER,ALLOCATABLE::IAg(:),JAg(:)
     complex(8),ALLOCATABLE::Ag(:)
      TYPE NODE
        INTEGER j
        real(8) c
        complex(8) cc
        TYPE(NODE),POINTER::next,prev
    ENDTYPE NODE
    TYPE IndexCeo
        TYPE(NODE),POINTER::head,tail
    ENDTYPE IndexCeo
    TYPE MatLink
        TYPE(IndexCeo),ALLOCATABLE::ht(:) !头尾指针
        INTEGER,ALLOCATABLE::GAi(:)
    ENDTYPE MatLink
    TYPE(NODE),POINTER::Q,P
    TYPE(MatLink) Lgt
     
     
  ALLOCATE(Lgt.ht(n),Lgt.GAi(n))
  allocate(IAg(n+1))
  !print*,"求逆"
    DO k=1,n!指向第一列元素
    ALLOCATE(Q);Q.j=1;Q.cc=0.0;NULLIFY(Q.next,Q.prev);Lgt.ht(k).head=>Q;Lgt.ht(k).tail=>Q
    end do
    Lgt.Gai=1
    do i1=1,n
        do j1=1,n
            i=i1;j=j1
            IF(( ABS(Gt(i1,j1)) )>1E-45) THEN
                    !为链表LK添加元素
                    P=>Lgt.ht(i).tail
                    DO 
                        !print*,ke
                        IF(j==P.j) THEN
                            P.cc=P.cc+Gt(i1,j1);EXIT
                        ELSEIF(j<P.j) THEN
                            P=>P.prev
                        ELSEIF(j>P.j) THEN !后方插入节点
                            ALLOCATE(Q);Q.j=j;Q.cc=Gt(i1,j1);Lgt.GAi(i)=Lgt.GAi(i)+1
                            Q.prev=>P;Q.next=>P.next
                            IF(ASSOCIATED(P.next)) THEN
                                P.next.prev=>Q
                            ELSE
                                Lgt.ht(i).tail=>Q
                            ENDIF
                            P.next=>Q;EXIT
                        ENDIF
                    ENDDO
            endif  
        end do
    end do
    IAg(1)=1
     DO i=1,n
        IAg(i+1)=IAg(i)+Lgt.GAi(i)
     ENDDO
         i=IAg(n+1)-1
    ALLOCATE(Ag(i),JAg(i))
     DO i=1,n
        P=>Lgt.ht(i).head
        DO j=IAg(i),IAg(i)+Lgt.GAi(i)-1
            JAg(j)=P.j;Ag(j)=P.cc;Q=>P.next;deALLOCATE(P);P=>Q
        ENDDO
     end do
      i=IAg(n+1)-1

    call pard(Ag,JAg,i,IAg,n,gt_inv)
 
    deallocate(iag,jag,ag,lgt.ht,lgt.gai)
  
    end 
    
    subroutine pard(A,JA,lda,IA,ldb,x)
    IMPLICIT NONE
    INTEGER lda, ldb, maxfct, mnum, mtype, phase, n, nrhs, error, msglvl
    INTEGER ja(lda),ia(ldb+1), iparm(64) 
    INTEGER*8 pt(64)
    integer,allocatable::perm(:)
    integer,external:: mkl_get_max_threads
    REAL*8 ddum
    COMPLEX(8) a(lda),b(ldb,ldb),x(ldb,ldb),xx(ldb),bb(ldb)
    integer i,j
    allocate( perm(ldb) )
   b=dcmplx(0.0,0.0)
    do i=1,ldb
        b(i,i)=dcmplx(1.0,0.0)
    end do
 
    
    pt=0 !指针初始化
    maxfct=1;mnum=1;mtype=13;phase=13;n=ldb;perm=0;nrhs=ldb;x=dble(0.0)
    iparm(1) = 0
    !iparm(1) = 1 ! iparm不使用默认值
    !iparm(2:64)=0; iparm(2)=3; iparm(4) =32; iparm(10) =8; iparm(18) =-1; iparm(24) =1; 
    error = 0 ! 错误信息初值置0 
    msglvl = 0 ! 不显示统计信息 

    call pardiso (pt, maxfct, mnum, mtype, phase, n, a, ia, ja, perm, nrhs, iparm, msglvl, b, x, error)
    if(error/=0) print*,'失败'

    phase=-1
    call pardiso (pt, maxfct, mnum, mtype, phase, n, a, ia, ja, perm, nrhs, iparm, msglvl, b, x, error)
    deallocate( perm )

    end   
  
    
    
!*******************M—范数************************  
   real function normm(x,n)
   implicit none
    integer n,i
    real(8) x(n)
    real dot
    normm=sqrt(dot(x,x,n))
    end function
    
    
!*******************M-内积************************    
   real function dot(x,y,n)
    use SparseMatrix,only:M,IM,JM
    implicit none
    integer n,i
    real(8) x(n),mx(n),y(n),t,mx1(n),mx2(n)
    call atmux(n,x,mx1,m,jm,im)
    call amux ( n, x, mx2,m, jm, im )
    mx=mx1+mx2
    t=0.0
    do i=1,n
        t=t+y(I)*mx(i)
    end do
    dot=t
    end function  
    
  