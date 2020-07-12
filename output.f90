subroutine subspace()
    !模型降阶后进行频点计算,测深或者测网
    !求解逆矩阵
    !构造传递函数（与陆地类似）
    !求解每个频点的Hy
    !Liu JiRen and Zhang JiFeng
    use blas95
    use sub,only:V,Vt,T,mi
    use Space_Domain_Field,only:EH!,Hy,hz,hx
    use net,only:num_vector,XYZ
    use FREQUENCY,only:num_freq,freq
    use sub,only:norm,mi
    use GLOBAL_CONST
    use SparseMatrix,only:G,Ig,jg,gg,num
   implicit none 
   complex(8) Gt(mi+1,mi+1),Gt_inv(mi+1,mi+1),E1(mi+1,1)
   complex(8),allocatable::GD(:,:)
   COMPLEX(8)  hyy,ex1,ex2,ex3,ex4,ez2,ez1,c,ez3,ez4,hyy1,ey,ez,hzz,ey1,hxx!,scal
    real(8)   omg,dx
    integer k,i,j,m,i1,j1,k1
    integer int1,int2
    allocate(Ig(mi+1+1))
    allocate(GD(num_vector,mi+1))
    Gt=dcmplx(0.0,0.0)
    E1=dcmplx(0.0,0.0)
    E1(1,1)=dcmplx(1.0,0.0)
    GT(:,:)=T(:,:)
    do i=1,mi+1
            if(Gt(i,i)==0)print*,i
    end do
     E1(:,1)=norm*E1(:,1)
 !call csr_store(Gt,mi+1)   
!******************************频点计算******************************        

close(101)
open(101,file="time.txt",access='append')
!!$OMP PARALLEL DO 
   do k=1,num_freq
       call system_clock(int1)
       print*,"正在计算频点",k
           omg=dble(2.0)*pi*FREQ(k) 
            do i=1,mi+1,1
        Gt(i,i)=Gt(i,i)+dcmplx(0.0,omg)      !Gt+iwI
            end do
        ! print*,Gt
        call CSR_inv(Gt,mi+1,Gt_inv)         !求逆
        Gt_inv=Gt_inv                       !传递函数
        GD=dcmplx(0.0,0.0)
        do i1=1,mi+1
           do j1=1,num_vector
               GD(j1,i1)=dcmplx(0.0,0.0)
               do k1=1,mi+1
                   GD(j1,i1)=GD(j1,i1)+V(j1,k1)*gt_inv(k1,i1)
               end do
           end do
       end do
         do i1=1,1
           do j1=1,num_vector
               EH(j1,i1)=dcmplx(0.0,0.0)
               do k1=1,mi+1
                   EH(j1,i1)=EH(j1,i1)+GD(j1,k1)*E1(k1,i1)
               end do
           end do
       end do
    
        !求平均Hy
        call Cal_Assist_Component1(omg,k)       !求Ex,Hz并输出

        GT(:,:)=T(:,:)
        
        call system_clock(int2)
        write(101,*) FREQ(k),(int2-int1)/10000
   end do
!!$OMP END  PARALLEL DO
   deallocate(GD,T,V,Vt,EH)
 100 format(1X,3F12.3)    
    end subroutine
    

    
    
  !*******************计算H场************************ 
      SUBROUTINE Cal_Assist_Component1(omg,nf)
     USE NET,ONLY:M_POS,M_N,NUM_MS,flag,num_vector
    ! use file
    USE Space_Domain_Field
    USE GLOBAL_CONST
    use FREQUENCY
    IMPLICIT NONE
    INTEGER i,j,k,ii,jj,m,nf
    REAL i1, r1, i2, r2,dx
    REAL(8)::omg,y0,res
    real(8) xx(4),yy(4),zz(4),Zm,xm,ym,x0,z0,a,b,x12(2,3),Ni(4),bi(4),ci(4),di(4),Ve
    COMPLEX(8) Ex0,Hy0,hyy,ex1,ex2,ex3,ex4,ez2,ez1,Exx,c,ez3,ez4,Ey1,e
    !CHARACTER*80 Outfile,cmdfile
    complex(8)Ex_z,Ex_y,Ez_x,Ey_x,Tzy,Tzx
    complex(8) Hy,Hz,hx
    real(8) hz1,ex
 
    !计算Ex、Ez沿z、x的导数
  
   do k=1,num_ms
       x12(1,:)=m_pos(:,k);x12(2,:)=m_pos(:,k)
       m=m_n(k)
     call qiu_E_H(ex1,ey1,hy,hz,hx,m,x12,omg)
if(flag==1)then !TE模式 Ey/hx
   close(11);close(12)
         open(11,file="TE_PHASE_RES.txt",access='append')
         IF(NF==1.AND.K==1)WRITE(11,"(6a15)") "X","Y","Z","FREQ","PHASE","RES"
         write(11,"(6F15.6)") x12(1,1),x12(1,2),X12(1,3),omg/2/pi,180/pi*atan(imag(ey1/hx)/real(ey1/hx)), dble(1.0)/OMG/MIU*abs(ey1/hx)**2 
else   !TM模式Ex/hy
     close(11)
         open(11,file="TM_PHASE_RES.txt",access='append')
         IF(NF==1.AND.K==1)WRITE(11,"(6a15)") "X","Y","Z","FREQ","PHASE","RES"
          write(11,"(6F15.6)") x12(1,1),x12(1,2),X12(1,3),omg/2/pi,180/pi*atan(imag(ex1/hy)/real(ex1/hy)),dble(1.0)/OMG/MIU*abs(ex1/hy)**2  

end if
  end do

    100 format(1X,4F12.3)  
    ENDSUBROUTINE


  subroutine qiu_E_H(Ex1,EY1,hy,hz,hx,m,x12,s)
    use net,ONLY:XYZ,N3,i3,pp,NUM_VECTOR,NUM_UNIT,NUM_POINT
    USE Space_Domain_Field,ONLY:EH
    implicit none
     REAL(8) S
    complex(8) Hy,Hz,hx
    integer i,j,k,n,m,pos
    real(8) x12(2,3),xx(8),yy(8),zz(8)
    complex(8) Ex_z,Ex_y,Ez_x,Ey_x,Ex1,Ey1,Ez1,e,Ez_y,Ey_z
     real(8) Ve,omg
     real(8) dmax,w
     real(8) Nxi(4),Nxi_z(4),Nzi_x(4),Nxi_y(4),Nzi_y(4),Nyi_z(4),Nyi_x(4),NYI(4)
     real(8) xc,yc,zc,x0,y0,z0,pi,u0,a,b,c
    PI=3.14159265358979323846
    U0=4.E-7*PI
    pos=m;
         Ex1=DCMPLX(0.0,0.0);Ey1=DCMPLX(0.0,0.0);Ez1=DCMPLX(0.0,0.0);e=DCMPLX(0.0,0.0);ex_y=DCMPLX(0.0,0.0);ex_z=DCMPLX(0.0,0.0);ez_x=DCMPLX(0.0,0.0);ey_x=DCMPLX(0.0,0.0);n=0;ez_y=DCMPLX(0.0,0.0);Ey_z=DCMPLX(0.0,0.0)
    do k=1,num_unit
      if(dble(1.0)/pp(k)>10000)then
        do j=1,8
           if( I3(k,j)==pos) goto 101    !定点
           if(j==8) goto 102
        end do       
101     Xx=(/Xyz(I3(k,1),1),Xyz(I3(k,2),1),Xyz(I3(k,3),1),Xyz(I3(k,4),1),Xyz(I3(k,5),1),Xyz(I3(k,6),1),Xyz(I3(k,7),1),Xyz(I3(k,8),1)/)
        Yy=(/Xyz(I3(k,1),2),Xyz(I3(k,2),2),Xyz(I3(k,3),2),Xyz(I3(k,4),2),Xyz(I3(k,5),2),Xyz(I3(k,6),2),Xyz(I3(k,7),2),Xyz(I3(k,8),2)/)
        Zz=(/Xyz(I3(k,1),3),Xyz(I3(k,2),3),Xyz(I3(k,3),3),Xyz(I3(k,4),3),Xyz(I3(k,5),3),Xyz(I3(k,6),3),Xyz(I3(k,7),3),Xyz(I3(k,8),3)/)
        a=abs(xx(2)-xx(1)) !边长
        b=abs(yy(3)-yy(2))
        c=abs(zz(5)-zz(1))
        xc=(xx(1)+xx(2))/dble(2.0) !中心坐标
        yc=(yy(1)+yy(3))/dble(2.0)
        zc=(zz(1)+zz(5))/dble(2.0)
        x0=x12(1,1);y0=x12(1,2);z0=x12(1,3)
        
        Nxi(1)=1.0/(b*c)*(yc+b/dble(2.0)-y0)*(zc+c/dble(2.0)-z0)
        Nxi(2)=1.0/(b*c)*(y0+b/dble(2.0)-yc)*(zc+c/dble(2.0)-z0)
        Nxi(3)=1.0/(b*c)*(yc+b/dble(2.0)-y0)*(z0+c/dble(2.0)-zc)
        Nxi(4)=1.0/(b*c)*(y0+b/dble(2.0)-yc)*(z0+c/dble(2.0)-zc)
          !Nyi
        Nyi(1)=1.0/(a*c)*(xc+a/dble(2.0)-x0)*(zc+c/dble(2.0)-z0)
        Nyi(2)=1.0/(a*c)*(xc+a/dble(2.0)-x0)*(z0+c/dble(2.0)-zc)
        Nyi(3)=1.0/(a*c)*(x0+a/dble(2.0)-xc)*(zc+c/dble(2.0)-z0)
        Nyi(4)=1.0/(a*c)*(x0+a/dble(2.0)-xc)*(z0+c/dble(2.0)-zc)
        
     !Nx_z   
        Nxi_z(1)=1.0/(b*c)*(yc+b/dble(2.0)-y0)*(-1)
        Nxi_z(2)=1.0/(b*c)*(y0+b/dble(2.0)-yc)*(-1)
        Nxi_z(3)=1.0/(b*c)*(yc+b/dble(2.0)-y0)*(1)
        Nxi_z(4)=1.0/(b*c)*(y0+b/dble(2.0)-yc)*(1)
     !Nz_x   
        Nzi_x(1)=1.0/(b*a)*(yc+b/dble(2.0)-y0)*(-1)
        Nzi_x(2)=1.0/(b*a)*(yc+b/dble(2.0)-y0)*(1)
        Nzi_x(3)=1.0/(b*a)*(y0+b/dble(2.0)-yc)*(-1)
        Nzi_x(4)=1.0/(b*a)*(y0+b/dble(2.0)-yc)*(1)
    !Nz_y    
        Nzi_y(1)=1.0/(b*a)*(xc+a/dble(2.0)-x0)*(-1)
        Nzi_y(2)=1.0/(b*a)*(x0+a/dble(2.0)-xc)*(-1)
        Nzi_y(3)=1.0/(b*a)*(xc+a/dble(2.0)-x0)*(1)
        Nzi_y(4)=1.0/(b*a)*(x0+a/dble(2.0)-xc)*(1)
     !Nx_y   
        Nxi_y(1)=1.0/(b*c)*(zc+c/dble(2.0)-z0)*(-1)
        Nxi_y(2)=1.0/(b*c)*(zc+c/dble(2.0)-z0)*(1)
        Nxi_y(3)=1.0/(b*c)*(z0+c/dble(2.0)-zc)*(-1)
        Nxi_y(4)=1.0/(b*c)*(z0+c/dble(2.0)-zc)*(1)
      !Ny_x
        Nyi_x(1)=1.0/(a*c)*(zc+c/dble(2.0)-z0)*(-1)
        Nyi_x(2)=1.0/(a*c)*(z0+c/dble(2.0)-zc)*(-1)
        Nyi_x(3)=1.0/(a*c)*(zc+c/dble(2.0)-z0)*(1)
        Nyi_x(4)=1.0/(a*c)*(z0+c/dble(2.0)-zc)*(1)
      !Ny_z
        Nyi_z(1)=1.0/(a*c)*(xc+a/dble(2.0)-x0)*(-1)
        Nyi_z(2)=1.0/(a*c)*(xc+a/dble(2.0)-x0)*(1)
        Nyi_z(3)=1.0/(a*c)*(x0+a/dble(2.0)-xc)*(-1)
        Nyi_z(4)=1.0/(a*c)*(x0+a/dble(2.0)-xc)*(1)
        do i=1,4
         ex1=ex1+ Nxi(i)*Eh(N3(k,i),1)   
         ey1=ey1+ Nyi(i)*Eh(N3(k,i+4),1)
         ex_z=ex_z+Nxi_z(i)*Eh(N3(k,i),1)
         Ez_x=ez_x+Nzi_x(i)*Eh(N3(k,i+8),1)
         Ey_x=Ey_x+Nyi_x(i)*Eh(N3(k,i+4),1)
         Ex_y=Ex_y+Nxi_y(i)*Eh(N3(k,i),1)
         Ez_y=Ez_y+Nzi_y(i)*Eh(N3(k,i+8),1)
         Ey_z=Ey_z+Nyi_z(i)*Eh(N3(k,i+4),1)
        end do
        n=n+1
102     end if 
        end do
        ex1=ex1/n
        ey1=ey1/n
        Hy=-dcmplx(0.0,1.0)*(Ex_z/n-Ez_x/n)/s/u0!
        Hz=dcmplx(-(imag(Ey_x)/n-imag(Ex_y)/n),-(real(Ey_x)/n-real(Ex_y)/n))/s/u0!
        hx=-dcmplx(0.0,1.0)*(Ez_y/n-Ey_z/n)/s/u0
    end 