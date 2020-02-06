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
   
    
subroutine aplb1 ( nrow, ncol, job, a, ja, ia, b, jb, ib, c, jc, ic, &
  nzmax, ierr )
    

!*****************************************************************************80
!
!! APLB1 performs the sum C = A + B for sorted CSR matrices.
!
!  Discussion:
!
!    The difference between this routine and APLB is that here the 
!    resulting matrix is such that the elements of each row are sorted,
!    with increasing column indices in each row, provided the original
!    matrices are sorted in the same way.
!
!    This routine will not work if either of the two input matrices is 
!    not sorted.
!
!  Modified:
!
!    11 January 2004
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
!    Compressed Sparse Row format with entries sorted.
!
! b,
! jb,
! ib      =  Matrix B in compressed sparse row format with entries sorted
!        ascendly in each row
!
! nzmax      = integer ( kind = 4 ). The  length of the arrays c and jc.
!         amub will stop if the result matrix C  has a number
!         of elements that exceeds exceeds nzmax. See ierr.
!
! on return:
!
! c,
! jc,
! ic      = resulting matrix C in compressed sparse row sparse format
!         with entries sorted ascendly in each row.
!
! ierr      = integer ( kind = 4 ). serving as error message.
!         ierr = 0 means normal return,
!         ierr > 0 means that amub stopped while computing the
!         i-th row  of C with i = ierr, because the number
!         of elements in C exceeds nzmax.
!
  implicit none

  integer ( kind = 4 ) nrow

  real ( kind = 8 ) a(*)
  real ( kind = 8 ) b(*)
  real ( kind = 8 ) c(*)
  integer ( kind = 4 ) i
  integer ( kind = 4 ) ia(nrow+1)
  integer ( kind = 4 ) ib(nrow+1)
  integer ( kind = 4 ) ic(nrow+1)
  integer ( kind = 4 ) ierr
  integer ( kind = 4 ) j1
  integer ( kind = 4 ) j2
  integer ( kind = 4 ) ja(*)
  integer ( kind = 4 ) jb(*)
  integer ( kind = 4 ) jc(*)
  integer ( kind = 4 ) job
  integer ( kind = 4 ) ka
  integer ( kind = 4 ) kamax
  integer ( kind = 4 ) kb
  integer ( kind = 4 ) kbmax
  integer ( kind = 4 ) kc
  integer ( kind = 4 ) ncol
  integer ( kind = 4 ) nzmax
  logical values

  values = ( job /= 0 )
  ierr = 0
  kc = 1
  ic(1) = kc

  do i = 1, nrow

    ka = ia(i)
    kb = ib(i)
    kamax = ia(i+1) - 1
    kbmax = ib(i+1) - 1

    do

      if ( ka <= kamax ) then
        j1 = ja(ka)
      else
        j1 = ncol + 1
      end if

      if ( kb <= kbmax ) then
        j2 = jb(kb)
      else
        j2 = ncol + 1
      end if
!
!  Three cases
!
      if ( j1 == j2 ) then
        if ( values ) then
          c(kc) = a(ka) + b(kb)
        end if
        jc(kc) = j1
        ka = ka + 1
        kb = kb + 1
        kc = kc + 1
      else if ( j1 < j2 ) then
        jc(kc) = j1
        if ( values ) then
          c(kc) = a(ka)
        end if
        ka = ka + 1
        kc = kc + 1
      else if ( j2 < j1 ) then
        jc(kc) = j2
        if ( values ) then
          c(kc) = b(kb)
        end if
        kb = kb + 1
        kc = kc + 1
      end if

      if ( nzmax < kc ) then
        ierr = i
        return
      end if

      if ( kamax < ka .and. kbmax < kb ) then
        exit
      end if

     end do

     ic(i+1) = kc

  end do

  return
    end
    
   
    subroutine amub ( nrow, ncol, job, a, ja, ia, b, jb, ib, c, jc, ic, nzmax, &
  iw, ierr )

!*****************************************************************************80
!
!! AMUB performs the matrix product C = A * B.
!
!  Discussion:
!
!    The column dimension of B is not needed.
!
!  Modified:
!
!    08 January 2004
!
!  Author:
!
!    Youcef Saad
!
!  Parameters:
!
!    Input, integer ( kind = 4 ) NROW, the row dimension of the matrix.
!
!    Input, integer ( kind = 4 ) NCOL, the column dimension of the matrix.
!
!    Input, integer ( kind = 4 ) JOB, job indicator.  When JOB = 0, only the
!    structure is computed, that is, the arrays JC and IC, but the real values
!    are ignored.
!
!    Input, real A(*), integer ( kind = 4 ) JA(*), IA(NROW+1), the matrix in CSR
!    Compressed Sparse Row format.
!
!    Input, b, jb, ib, matrix B in compressed sparse row format.
!
!    Input, integer ( kind = 4 ) NZMAX, the length of the arrays c and jc.
!    The routine will stop if the result matrix C  has a number
!    of elements that exceeds exceeds NZMAX.
!
! on return:
!
! c,
! jc,
! ic    = resulting matrix C in compressed sparse row sparse format.
!
! ierr      = integer ( kind = 4 ). serving as error message.
!         ierr = 0 means normal return,
!         ierr > 0 means that amub stopped while computing the
!         i-th row  of C with i = ierr, because the number
!         of elements in C exceeds nzmax.
!
! work arrays:
!
!  iw      = integer ( kind = 4 ) work array of length equal to the number of
!         columns in A.
!
  implicit none

  integer ( kind = 4 ) ncol
  integer ( kind = 4 ) nrow
  integer ( kind = 4 ) nzmax

  real ( kind = 8 ) a(*)
  real ( kind = 8 ) b(*)
  real ( kind = 8 ) c(nzmax)
  integer ( kind = 4 ) ia(nrow+1)
  integer ( kind = 4 ) ib(ncol+1)
  integer ( kind = 4 ) ic(ncol+1)
  integer ( kind = 4 ) ierr
  integer ( kind = 4 ) ii
  integer ( kind = 4 ) iw(ncol)
  integer ( kind = 4 ) ja(*)
  integer ( kind = 4 ) jb(*)
  integer ( kind = 4 ) jc(nzmax)
  integer ( kind = 4 ) jcol
  integer ( kind = 4 ) jj
  integer ( kind = 4 ) job
  integer ( kind = 4 ) jpos
  integer ( kind = 4 ) k
  integer ( kind = 4 ) ka
  integer ( kind = 4 ) kb
  integer ( kind = 4 ) len
  real ( kind = 8 ) scal
  logical values

  values = ( job /= 0 )
  len = 0
  ic(1) = 1
  ierr = 0
!
!  Initialize IW.
!
  iw(1:ncol) = 0

  do ii = 1, nrow
!
!  Row I.
!
    do ka = ia(ii), ia(ii+1)-1

      if ( values ) then
        scal = a(ka)
      end if

      jj = ja(ka)

      do kb = ib(jj), ib(jj+1)-1

           jcol = jb(kb)
           jpos = iw(jcol)

           if ( jpos == 0 ) then
              len = len + 1
              if ( nzmax < len ) then
                 ierr = ii
                 return
              end if
              jc(len) = jcol
              iw(jcol)= len
              if ( values ) then
                c(len) = scal * b(kb)
              end if
           else
              if ( values ) then
                c(jpos) = c(jpos) + scal * b(kb)
              end if
           end if

         end do

    end do

    do k = ic(ii), len
      iw(jc(k)) = 0
    end do

    ic(ii+1) = len + 1

  end do

  return
    end

    
    
  subroutine aplsca ( nrow, a, ja, ia, scal, iw )

!*****************************************************************************80
!
!! APLSCA adds a scalar to the diagonal entries of a sparse matrix A :=A + s I
!
!  Discussion:
!
!    The column dimension of A is not needed.
!
!    important: the matrix A may be expanded slightly to allow for
!    additions of nonzero elements to previously nonexisting diagonals.
!    The is no checking as to whether there is enough space appended
!    to the arrays a and ja. if not sure allow for n additional
!    elements.
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
!    Input, integer ( kind = 4 ) NROW, the row dimension of the matrix.
!
!    Input, real A(*), integer ( kind = 4 ) JA(*), IA(NROW+1), the matrix in CSR
!    Compressed Sparse Row format.
!
!    Input, real SCAL, a scalar to be added to the diagonal entries.
!
! on return:
!
!
! a,
! ja,
! ia      = matrix A with diagonal elements shifted (or created).
!
! iw    = integer ( kind = 4 ) work array of length n. On return iw will
!         contain  the positions of the diagonal entries in the
!         output matrix. (i.e., a(iw(k)), ja(iw(k)), k = 1,...n,
!         are the values/column indices of the diagonal elements
!         of the output matrix. ).
!
  implicit none

  integer ( kind = 4 ) nrow

  real ( kind = 8 ) a(*)
  integer ( kind = 4 ) ia(nrow+1)
  integer ( kind = 4 ) icount
  integer ( kind = 4 ) ii
  integer ( kind = 4 ) iw(*)
  integer ( kind = 4 ) j
  integer ( kind = 4 ) ja(*)
  integer ( kind = 4 ) k
  integer ( kind = 4 ) k1
  integer ( kind = 4 ) k2
  integer ( kind = 4 ) ko
  real ( kind = 8 ) scal
  logical test

  call diapos ( nrow, ja, ia, iw )
  icount = 0

  do j = 1, nrow

     if ( iw(j) == 0 ) then
        icount = icount + 1
     else
        a(iw(j)) = a(iw(j)) + scal
     end if

  end do
!
!  If no diagonal elements to insert in data structure, return.
!
  if ( icount == 0 ) then
    return
  end if
!
!  Shift the nonzero elements if needed, to allow for created
!  diagonal elements.
!
  ko = ia(nrow+1) + icount
!
!  Copy rows backward.
!
  do ii = nrow, 1, -1
!
!  Go through row II.
!
     k1 = ia(ii)
     k2 = ia(ii+1) - 1
     ia(ii+1) = ko
     test = ( iw(ii) == 0 )

     do k = k2, k1, -1

        j = ja(k)

        if ( test .and. j < ii ) then
           test = .false.
           ko = ko - 1
           a(ko) = scal
           ja(ko) = ii
           iw(ii) = ko
        end if

        ko = ko - 1
        a(ko) = a(k)
        ja(ko) = j

    end do
!
!  The diagonal element has not been added yet.
!
     if ( test ) then
        ko = ko - 1
        a(ko) = scal
        ja(ko) = ii
        iw(ii) = ko
     end if

  end do

  ia(1) = ko

  return
    end
    
      subroutine aplsca_complex ( nrow, a, ja, ia, scal )

!*****************************************************************************80
!
!! APLSCA adds a scalar to the diagonal entries of a sparse matrix A :=A + s I
!
!  Discussion:
!
!    The column dimension of A is not needed.
!
!    important: the matrix A may be expanded slightly to allow for
!    additions of nonzero elements to previously nonexisting diagonals.
!    The is no checking as to whether there is enough space appended
!    to the arrays a and ja. if not sure allow for n additional
!    elements.
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
!    Input, integer ( kind = 4 ) NROW, the row dimension of the matrix.
!
!    Input, real A(*), integer ( kind = 4 ) JA(*), IA(NROW+1), the matrix in CSR
!    Compressed Sparse Row format.
!
!    Input, real SCAL, a scalar to be added to the diagonal entries.
!
! on return:
!
!
! a,
! ja,
! ia      = matrix A with diagonal elements shifted (or created).
!
! iw    = integer ( kind = 4 ) work array of length n. On return iw will
!         contain  the positions of the diagonal entries in the
!         output matrix. (i.e., a(iw(k)), ja(iw(k)), k = 1,...n,
!         are the values/column indices of the diagonal elements
!         of the output matrix. ).
!
  implicit none

  integer ( kind = 4 ) nrow

  complex( kind = 8 ) a(*)
  integer ( kind = 4 ) ia(nrow+1)
  integer ( kind = 4 ) icount
  integer ( kind = 4 ) ii
  integer ( kind = 4 ) iw(nrow)
  integer ( kind = 4 ) j
  integer ( kind = 4 ) ja(*)
  integer ( kind = 4 ) k
  integer ( kind = 4 ) k1
  integer ( kind = 4 ) k2
  integer ( kind = 4 ) ko
  complex ( kind = 8 ) scal
  logical test

  call diapos ( nrow, ja, ia, iw )
  icount = 0

  do j = 1, nrow

     if ( iw(j) == 0 ) then
        icount = icount + 1
     else
        a(iw(j)) = a(iw(j)) + scal
     end if

  end do
!
!  If no diagonal elements to insert in data structure, return.
!
  if ( icount == 0 ) then
    return
  end if
!
!  Shift the nonzero elements if needed, to allow for created
!  diagonal elements.
!
  ko = ia(nrow+1) + icount
!
!  Copy rows backward.
!
  do ii = nrow, 1, -1
!
!  Go through row II.
!
     k1 = ia(ii)
     k2 = ia(ii+1) - 1
     ia(ii+1) = ko
     test = ( iw(ii) == 0 )

     do k = k2, k1, -1

        j = ja(k)

        if ( test .and. j < ii ) then
           test = .false.
           ko = ko - 1
           a(ko) = scal
           ja(ko) = ii
           iw(ii) = ko
        end if

        ko = ko - 1
        a(ko) = a(k)
        ja(ko) = j

    end do
!
!  The diagonal element has not been added yet.
!
     if ( test ) then
        ko = ko - 1
        a(ko) = scal
        ja(ko) = ii
        iw(ii) = ko
     end if

  end do

  ia(1) = ko

  return
    end
    
    
    
    subroutine diapos ( n, ja, ia, idiag )

!*****************************************************************************80
!
!! DIAPOS returns the positions of the diagonal elements of a sparse matrix.
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
!    Input, JA(*), IA(N+1), the matrix information, (but no values) 
!    in CSR Compressed Sparse Row format.
!
!    Output, integer ( kind = 4 ) IDIAG(N); the I-th entry of IDIAG points to the 
!    diagonal element A(I,I) in the arrays A and JA.  That is,
!    A(IDIAG(I)) = element A(I,I) of matrix A.  If no diagonal element 
!    is found, the entry is set to 0.
!
  implicit none

  integer ( kind = 4 ) n

  integer ( kind = 4 ) i
  integer ( kind = 4 ) ia(n+1)
  integer ( kind = 4 ) idiag(n)
  integer ( kind = 4 ) ja(*)
  integer ( kind = 4 ) k

  idiag(1:n) = 0
!
!  Sweep through the data structure.
!
  do i = 1, n
    do k = ia(i), ia(i+1) -1
      if ( ja(k) == i ) then
        idiag(i) = k
      end if
    end do
  end do

  return
    end
    
    
subroutine diamua ( nrow, job, a, ja, ia, diag, b )

!*****************************************************************************80
!
!! DIAMUA performs the matrix by matrix product B = Diag * A.
!
!  Discussion:
!
!    The column dimension of A is not needed.
!
!    The algorithm is in-place; that is, B can take the place of A.
!    in this case use job=0.
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
!    Input, integer ( kind = 4 ) NROW, the row dimension of the matrix.
!
!    Input, integer ( kind = 4 ) JOB, indicates the job to be done.
!    0, means get array B only;
!    1, means get B, and the integer ( kind = 4 ) arrays IB and JB.
!
!    Input, real A(*), integer ( kind = 4 ) JA(*), IA(NROW+1), the matrix in CSR
!    Compressed Sparse Row format.
!
!    Input, real DIAG(N), a diagonal matrix stored as a vector.
!
!    Output, real B(*), integer ( kind = 4 ) JB(*), 
!    integer ( kind = 4 ) IB(NROW+1), the resulting 
!    matrix B in compressed sparse row sparse format.
!
  implicit none

  integer ( kind = 4 ) nrow

  real ( kind = 8 ) a(*)
  real ( kind = 8 ) b(*)
  real ( kind = 8 ) diag(nrow)
  integer ( kind = 4 ) ia(nrow+1)
!  integer ( kind = 4 ) ib(nrow+1)
  integer ( kind = 4 ) ii
  integer ( kind = 4 ) ja(*)
!  integer ( kind = 4 ) jb(*)
  integer ( kind = 4 ) job
  integer ( kind = 4 ) k
  integer ( kind = 4 ) k1
  integer ( kind = 4 ) k2
  real ( kind = 8 ) scal

  do ii = 1, nrow
!
!  Normalize each row.
!
     k1 = ia(ii)
     k2 = ia(ii+1) - 1
     scal = diag(ii)
     b(k1:k2) = a(k1:k2) * scal

  end do

  if ( job == 0 ) then
    return
  end if

 ! ib(1) = ia(1)

  do ii = 1, nrow
  ! ib(ii) = ia(ii)
    do k = ia(ii), ia(ii+1)-1
   !   jb(k) = ja(k)
    end do
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
    
    
      MODULE HANKLE_COE
    IMPLICIT NONE
    REAL::HANKLE(0:1,-53:58)=&
    &(/0.1154026E-05, 0.2930619E-11, 0.1452833E-05, 0.4522038E-11,&
    & 0.1829008E-05, 0.7361382E-11, 0.2302585E-05, 0.1135885E-11,&
    & 0.2898783E-05, 0.1849096E-10, 0.3649351E-05, 0.2853213E-10,&
    & 0.4594261E-05, 0.4644718E-10, 0.5783832E-05, 0.7166948E-10,&
    & 0.7281413E-05, 0.1166700E-09, 0.9166757E-05, 0.1800256E-09,&
    & 0.1154026E-04, 0.2930619E-09, 0.1452833E-04, 0.4522038E-09,&
    & 0.1829008E-04, 0.7361382E-09, 0.2302585E-04, 0.1135885E-08,&
    & 0.2898783E-04, 0.1849096E-08, 0.3649352E-04, 0.2853213E-08,&
    & 0.4594261E-04, 0.4644718E-08, 0.5783833E-04, 0.7166948E-08,&
    & 0.7281411E-04, 0.1166700E-07, 0.9166760E-04, 0.1800256E-07,&
    & 0.1154026E-03, 0.2930619E-07, 0.1452834E-03, 0.4522038E-07,&
    & 0.1829007E-03, 0.7361381E-07, 0.2302586E-03, 0.1135885E-06,&
    & 0.2898779E-03, 0.1849095E-06, 0.3649354E-03, 0.2853212E-06,&
    & 0.4594250E-03, 0.4644716E-06, 0.5783834E-03, 0.7166942E-06,&
    & 0.7281377E-03, 0.1166699E-05, 0.9166749E-03, 0.1800252E-05,&
    & 0.1154015E-02, 0.2930610E-05, 0.1452826E-02, 0.4522015E-05,&
    & 0.1828968E-02, 0.7361325E-05, 0.2302545E-02, 0.1135870E-04,&
    & 0.2898640E-02, 0.1849059E-04, 0.3649167E-02, 0.2853122E-04,&
    & 0.4593733E-02, 0.4644490E-04, 0.5783032E-02, 0.7166375E-04,&
    & 0.7279415E-02, 0.1166557E-03, 0.9163407E-02, 0.1799894E-03,&
    & 0.1153257E-01, 0.2929711E-03, 0.1451458E-01, 0.4519758E-03,&
    & 0.1826012E-01, 0.7355654E-03, 0.2297010E-01, 0.1134446E-02,&
    & 0.2887026E-01, 0.1845483E-02, 0.3626918E-01, 0.2844143E-02,&
    & 0.4547940E-01, 0.4621947E-02, 0.5694082E-01, 0.7109806E-02,&
    & 0.7098731E-01, 0.1152369E-01, 0.8809954E-01, 0.1764345E-01,&
    & 0.1082239E+00, 0.2840762E-01, 0.1312505E+00, 0.4297706E-01,&
    & 0.1550557E+00, 0.6803326E-01, 0.1763715E+00, 0.9978460E-01,&
    & 0.1856277E+00, 0.1510705E+00, 0.1697780E+00, 0.2035406E+00,&
    & 0.1034052E+00, 0.2712354E+00,-0.3025832E-01, 0.2760739E+00,&
    &-0.2275744E+00, 0.2166920E+00,-0.3621732E+00,-0.7837237E-01,&
    &-0.2055004E+00,-0.3406756E+00, 0.3373949E+00,-0.3606937E+00,&
    & 0.3176899E+00, 0.5130245E+00,-0.5137622E+00,-0.5947247E-01,&
    & 0.3091303E+00,-0.1951171E+00,-0.1267576E+00, 0.1992356E+00,&
    & 0.4619679E-01,-0.1385216E+00,-0.1809687E-01, 0.8793210E-01,&
    & 0.8354260E-02,-0.5506971E-01,-0.4473683E-02, 0.3456378E-01,&
    & 0.2619748E-02,-0.2175272E-01,-0.1601714E-02, 0.1371003E-01,&
    & 0.9977179E-03,-0.8646564E-02,-0.6262758E-03, 0.5454628E-02,&
    & 0.3943388E-03,-0.3441389E-02,-0.2486064E-03, 0.2171307E-02,&
    & 0.1568086E-03,-0.1369986E-02,-0.9892663E-04, 0.8643991E-03,&
    & 0.6241524E-04,-0.5453979E-03,-0.3938054E-04, 0.3441225E-03,&
    & 0.2484724E-04,-0.2171226E-03,-0.1567749E-04, 0.1369976E-03,&
    & 0.9891818E-05,-0.8643964E-04,-0.6241312E-05, 0.5453972E-04,&
    & 0.3938001E-05,-0.3441224E-04,-0.2484710E-05, 0.2171265E-04,&
    & 0.1567746E-05,-0.1369976E-04,-0.9891809E-06, 0.8643963E-05,&
    & 0.6241309E-06,-0.5453972E-05,-0.3938000E-06, 0.3441224E-05,&
    & 0.2484710E-06,-0.2171265E-05,-0.1567746E-06, 0.1369976E-05,&
    & 0.0000000E+00,-0.8643964E-06, 0.0000000E+00, 0.5453972E-06,&
    & 0.0000000E+00,-0.3441224E-06, 0.0000000E+00, 0.2171265E-06,&
    & 0.0000000E+00,-0.1369976E-06, 0.0000000E+00, 0.8643963E-07,&
    & 0.0000000E+00,-0.5453972E-07, 0.0000000E+00, 0.3441224E-07,&
    & 0.0000000E+00,-0.2171265E-07, 0.0000000E+00, 0.1369976E-07,&
    & 0.0000000E+00,-0.8643963E-08, 0.0000000E+00, 0.5453972E-08,&
    & 0.0000000E+00,-0.3441224E-08, 0.0000000E+00, 0.2171265E-08,&
    & 0.0000000E+00,-0.1369976E-08, 0.0000000E+00, 0.8643963E-09,&
    & 0.0000000E+00,-0.5453972E-09, 0.0000000E+00, 0.3441224E-09,&
    & 0.0000000E+00,-0.2171265E-09, 0.0000000E+00, 0.1369976E-09/)
    ENDMODULE HANKLE_COE

    
    SUBROUTINE Fwd1D(omg,PP,n,x,y,Ex,Hy)
    USE HANKLE_COE
    IMPLICIT NONE
    INTEGER layer,i,j,n
    REAL pe,x,y
    REAL u0,omg,PP(n),P((n+1)/2),TH((n+1)/2-1)
    COMPLEX ex,hy,rr1,r1,ctemp,ctemp1,ctemp2,CCOTH,CACOTH,f1,f2,f3,f4,hi3,hi4,hi5,hi6,hi7,hi8,hi9,er,ep,hr,hp
    REAL pi,delta,sinp,cosp,r,wu,temp,t1,t2
    REAL M(-53:58)
    COMPLEX KK1((n+1)/2),M1((n+1)/2)
    PARAMETER(pi=3.141592653,u0=1.2566370612E-6,pe=1.0*0.5d0/pi)
    layer=(n+1)/2;delta=log(10.)/10.;r=SQRT(x**2+y**2)
    cosp=x/r;sinp=y/r;wu=omg*u0
    P=PP(1:layer);TH=PP(layer+1:n)
    DO i=1,layer
        KK1(i)=CMPLX(0.,omg*u0/P(i))
    ENDDO
    DO i=-53,58
        M(i)=EXP(i*delta)/r
    ENDDO
    hi3=0.;hi4=0.;hi5=0.;hi6=0.;hi7=0.;hi8=0.;hi9=0.
    DO i=-53,58
        DO j=1,Layer
            M1(j)=CSQRT(M(i)**2-KK1(j))
        ENDDO
        rr1=1.;r1=1.;ctemp=0.;ctemp1=0.;ctemp2=0.
        DO j=layer,2,-1
            ctemp=M1(j-1)/M1(j)
            ctemp=ctemp*rr1
            t1=REAL(ctemp)
            t2=AIMAG(ctemp)
            IF(ABS(t2).LE..1E-35) t2=0.
            IF(ABS(t1).LE..1E-35) t1=0.
            ctemp=CMPLX(t1,t2)
            ctemp1=M1(j-1)*TH(j-1)
            ctemp1=CEXP(-2.*ctemp1)
            ctemp1=(-ctemp1+1.)/(ctemp1+1.)
            t1=REAL(ctemp1)
            t2=AIMAG(ctemp1)
            IF(ABS(t2).LE..1E-35) t2=0.
            IF(ABS(t1).LE..1E-35) t1=0.
            ctemp1=CMPLX(t1,t2)
            rr1=(ctemp+ctemp1)/(1.+ctemp*ctemp1)
            ctemp2=P(j-1)/P(j)*M1(j-1)/M1(j)*r1
            t1=REAL(ctemp2)
            t2=AIMAG(ctemp2)
            IF(ABS(t2).LE..1E-35) t2=0.
            IF(ABS(t1).LE..1E-35) t1=0.
            ctemp2=CMPLX(t1,t2)
            r1=(ctemp2+ctemp1)/(1.+ctemp2*ctemp1)
        ENDDO
        ctemp=M(i)+M1(1)/rr1
        f1=1/ctemp-1/(M(i)+M1(1))
        f2=(1/r1-1)*M1(1)
        f3=M1(1)/rr1/ctemp
        f4=M(i)/ctemp
        hi4=(f4-0.5)*HANKLE(1,i)+hi4
        hi5=(f3-0.5)*M(i)*HANKLE(0,i)+hi5
        hi6=f1*HANKLE(1,i)+hi6
        hi7=f2*M(i)*HANKLE(0,i)+hi7
        hi8=f2*HANKLE(1,i)+hi8
        hi9=f1*M(i)*HANKLE(0,i)+hi9
    ENDDO
    hi3=hi3/r;hi4=hi4/r;hi5=hi5/r;hi6=hi6/r
    hi7=hi7/r;hi8=hi8/r;hi9=hi9/r
    ctemp=CMPLX(0.,wu)
    er=pe*cosp*(ctemp*hi6/r-P(1)*hi7+P(1)*hi8/r)
    ep=pe*sinp*(ctemp*hi6/r-ctemp*hi9+P(1)*hi8/r)
    temp=SQRT(wu/(2*P(1)))*r
    ctemp=CMPLX(-temp,temp)
    er=er+pe*P(1)*cosp*(1+CEXP(ctemp)*(1-ctemp))/r**3
    ep=ep+pe*P(1)*sinp*(2-CEXP(ctemp)*(1-ctemp))/r**3
    ex=er*cosp-ep*sinp
    temp=pe/r/r/2
    hr=-pe*sinp*(hi4/r+hi5)-temp*sinp
    hp=pe*cosp*hi4/r+temp*cosp
    hy=hr*sinp+hp*cosp
    ENDSUBROUTINE 
    
    
!*******************求逆************************          
subroutine CSR_inv(Gt,n,Gt_inv)
    use LINKED_list,only:Lgt,Q,P
    implicit none
    integer i,j,k,n,i1,j1
    complex(8) Gt(n,n),GT_inv(n,n)
     INTEGER,ALLOCATABLE::IAg(:),JAg(:)
     complex(8),ALLOCATABLE::Ag(:)
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
       ! print*,b(i,i),i
    end do
    !print*,b
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