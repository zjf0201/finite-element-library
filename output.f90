subroutine subspace()
    !模型降阶后进行频点计算,测深或者测网
    !求解逆矩阵
    !构造传递函数（与陆地类似）
    !求解每个频点的Hy
    !Liu JiRen and Zhang JiFeng
    use file
    use sub,only:V,Vt,T,mi
    use Space_Domain_Field,only:right_b,EXYZ,Hy,hz,hx
    use net,only:num_point
    use FREQUENCY,only:num_freq,freq
    use sub,only:norm,mi
    use net,only:x,y,z,info
    use GLOBAL_CONST
    use ParaofElecDipole
    use SparseMatrix,only:G,Ig,jg,gg,num
   implicit none 
   complex(8) Gt(mi+1,mi+1),Gt_inv(mi+1,mi+1)
   COMPLEX(8)  hyy,ex1,ex2,ex3,ex4,ez2,ez1,c,ez3,ez4,hyy1,ey,ez,hzz,ey1,hxx,scal
    real(8)    E1(mi+1),omg,dx
    integer k,i,j,m
    allocate(Ig(mi+1+1))
    open(11,file=file1)
    open(12,file=file2)
    open(13,file=file3)
    open(14,file=file4)
    Gt=dcmplx(0.0,0.0)
    E1=dble(0.0)
    E1(1)=dble(1.0)
    GT(:,:)=T(:,:)
    do i=1,mi+1
            if(Gt(i,i)==0)print*,i
    end do
     E1(:)=norm*E1(:)
 !call csr_store(Gt,mi+1)   
!******************************频点计算******************************        
!!$OMP PARALLEL DO 
   do k=1,num_freq
       print*,"正在计算频点",k
           omg=dble(2.0)*pi*FREQ(k) 
            do i=1,mi+1,1
        Gt(i,i)=Gt(i,i)+dcmplx(0.0,omg)      !Gt+iwI
            end do

        call CSR_inv(Gt,mi+1,Gt_inv)         !求逆
        Gt_inv=Gt_inv*-1!*(dcmplx(0.0,-omg)) !传递函数
        Exyz=matmul(matmul(V,gt_inv),e1)
        !求平均Hy
        call Cal_Assist_Component1(omg)       !求Hy hx
!!$OMP END  PARALLEL DO
 ! call Cal_Assist_Component() 
!do k=1,num_freq  
  do j=1,1
      do m=1,1       
   do i=1,num_point
       if(y(i)==0+(j-1)*dx.and.z(i)==0.0.and.x(i)==0+(m-1)*dx)then
           ex1=exyz(i)
           ey1=exyz(i+num_point)
           hyy=Hy(i)
           hxx=Hx(i)
       end if            
   end do
    if(info==0)then
          write(11,*)  1.0/(omg*miu)*abs(ex1/hyy)**2     !TM
          write(12,*)  180/pi*atan(imag(ex1/hyy)/real(ex1/hyy))
    else
          write(13,*) 1.0/(omg*miu)*abs(ey1/hxx)**2      !TE
          write(14,*) 180/pi*atan(imag(ey1/hxx)/real(ey1/hxx))
           ! write(13,*) FREQ(k),abs(hxx)
            !write(14,*) FREQ(k),abs(hyy)
    end if
    
      end do 
  end do
        GT(:,:)=T(:,:)
   end do
!
 100 format(1X,3F12.3)    
    end subroutine
  
  !*******************计算H场************************ 
      SUBROUTINE Cal_Assist_Component1(omg)
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
                ctemp1=ctemp1+N_xi(i)*xi_x *(Exyz(ii+2*num_point))   !Ez to x
                ctemp2=ctemp2+N_gama(i)*gama_z*(Exyz(ii));           !ex to z
                ctemp3=ctemp3+N_xi(i)*xi_x *(Exyz(ii+num_point))     !eY  to x
                ctemp4=ctemp4+N_ETA(i)*ETA_Y *(Exyz(ii))             !Ex to y
                 
                ctemp5=ctemp5+N_eta(i)*eta_y *(Exyz(ii+2*num_point))                 !ez  to y
                ctemp6=ctemp6+N_gama(i)*gama_z *(Exyz(ii+num_point))                 !Ey to z
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
           Hy(k)=-dcmplx(0.0,1.0)/(omg*miu)*(Ex_z(k)-Ez_x(k))
           
              
           Ex_y(k)=Ex_y(k)/NU(k)!取平均
           Ey_x(k)=Ey_x(k)/NU(k)
       !     Hz(k)=-dcmplx(0,1)/(omg*miu)*(Ey_x(k)-Ex_y(k))
            
             Ez_y(k)=Ez_y(k)/NU(k)!取平均
             Ey_z(k)=Ey_z(k)/NU(k)
            Hx(k)=-dcmplx(0.0,1.0)/(omg*miu)*(Ez_y(k)-Ey_z(k))
    ENDDO

    deALLOCATE(NU)
    ENDSUBROUTINE