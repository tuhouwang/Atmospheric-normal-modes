module param_mod
    implicit none
    integer,        parameter :: rkind            = 8
    integer,        parameter :: MAX_FILENAME_LEN = 200

    real(rkind),    parameter :: pi = 4.0_rkind * atan(1.0_rkind)   ! 双精度pi值
    complex(rkind), parameter :: ci = cmplx(0.0, 1.0, kind=rkind)   ! 复数
end module

module cheb_mod
    use param_mod
    implicit none

    interface Cheb
        module procedure ChebReal
        module procedure ChebComplex
    end interface

    interface Convolution
        module procedure ConvolutionReal
        module procedure ConvolutionComplex
    end interface

contains

    !---------------------------------------------------------------------------  
    ! DESCRIPTION: 
    !> @brief Chebyshev正变换。已知原函数在CGL节点上的取值，求其Chebyshev变换后的展开系数.
    !> 给定原函数在n个Chebyshev-Gauss-Lobbatto节点z(1:n)上的取值fx(1:n)，
    !> 求n阶Chebyshev变换后的展开系数fk(0:n)。
    !
    !> @param[in]  z(1:n)   n个Chebyshev-Gauss-Lobbatto节点
    !> @param[in] fx(1:n)   原函数在给定 z(1:n) 处的值
    !> @param[out]  fk(0:n) 各阶 Chebyshev 多项式对应的展开系数
    !> @return none
    !---------------------------------------------------------------------------  
    subroutine ChebReal(fk, fx, z)
        implicit none
        real(rkind),	intent(in)    :: z(:)
        real(rkind),    intent(inout) :: fk(:), fx(:)
        real(rkind),	allocatable   :: T(:,:)
        integer                       :: m, i
        m = size(fx,1)
        allocate(T(m, m))             

        T = 0.0_rkind
        fk= 0.0_rkind
        
        do i=0, m-1
            T(:,i+1) = cos(i * acos(z))	  
        enddo   

        fx(1) = fx(1) * 0.5_rkind
        fx(m) = fx(m) * 0.5_rkind

        fk = matmul( transpose(T), fx ) * 2.0_rkind / (m - 1)
        
        fk(1) = fk(1) * 0.5_rkind
        fk(m) = fk(m) * 0.5_rkind
        
    end subroutine ChebReal

    subroutine ChebComplex(fk, fx, z)
        implicit none
        real(rkind),    intent(in)    :: z(:)
        complex(rkind), intent(inout) :: fk(:), fx(:)
        real(rkind),    allocatable   :: T(:, :)
        integer                       :: m, i

        m = size(fx, 1)
        allocate ( T(m, m) )

        T  = 0.0_rkind
        fk = 0.0_rkind

        do i = 0, m - 1
            T(:, i + 1) = cos( i * acos(z) )
        end do

        fx(1) = fx(1) * 0.5_rkind
        fx(m) = fx(m) * 0.5_rkind

        fk = matmul( transpose(T), fx ) * 2.0_rkind / (m - 1)

        fk(1) = fk(1) * 0.5_rkind
        fk(m) = fk(m) * 0.5_rkind

    end subroutine ChebComplex

   !---------------------------------------------------------------------------  
    ! DESCRIPTION: 
    !> @brief Chebyshev反变换。已知一组函数的Chebyshev变换系数，求取这组函数在指定位置处的函数近似值.
    !> 给定d个函数关于m阶Chebyshev变换的展开系数fk(0:m, 1:d)，以及给定的自变量离散点z(1:n)，
    !> 求取原函数在这些离散点上的函数近似值fx(1:n, 1:d)。
    !
    !> @param[in]  z(1:n)   表示原函数要在这些离散点上求取近似值
    !> @param[in]  fk(0:m, 1:d)   fk(:, j) 表示各阶 Chebyshev 多项式对应的展开系数（列向量）
    !> @param[out] fx(1:n, 1:d)   fx(:, j) 表示原函数在给定 z(1:n) 处的近似值
    !> @return none
    !---------------------------------------------------------------------------  
    subroutine InvChebMatrix(fx, fk, z)
        implicit none
        real(rkind),    intent(in)  :: z(:)
        complex(rkind), intent(in)  :: fk(:, :)
        complex(rkind), intent(out) :: fx(:, :)
        integer                     :: i
        real(rkind)                 :: T(size(z), size(fk, 1))

        T  = 0.0_rkind
        fx = 0.0_rkind

        do i = 0, size(fk, 1) - 1
            T(:, i+1) = cos( i * acos( z(:) ) )
        end do

        fx = matmul(T, fk)
    end subroutine InvChebMatrix

    subroutine DerivationMatrix(D, m)
        implicit none
        integer,      intent(in) :: m
        real(rkind),  intent(out):: D(m, m)
        integer                  :: i, j
    
        D = 0.0_rkind
        do i = 1, m
        do j = i + 1, m
            if ( mod((i+j), 2) == 1 ) D(i, j) = 2.0_rkind * j - 2.0_rkind
        end do
        end do
        D(1, :) = D(1, :) * 0.5_rkind
    
    end subroutine DerivationMatrix
    
    subroutine ConvolutionReal(Co, v)
        implicit none
        real(rkind), intent(in)   :: v  (:)
        real(rkind), intent(out)  :: Co (:, :)
        integer                   :: i, j, k, n
    
        Co = 0.0_rkind
    
		n  = size(v)
		do  i = 1, n
		do  k = 1, n
		     j= k-i+1
		     if (1 <= j .and. j <= n) then
			    Co(k,i) = Co(k,i) + v(j) 
		     endif
		     j = i-k+1
		     if (j <=n .and. j >= 1) then
				Co(k,i) = Co(k,i) + v(j) 
		     endif
		     if (k >1) then 
				j = i+k-1   
			    if (1 <= j .and. j <= n) then
					Co(k,i) = Co(k,i) + v(j)
				endif 
		     endif   
		     Co(k,i) = Co(k,i) * 0.5_rkind
		enddo
		enddo
    
    end subroutine ConvolutionReal
    
    subroutine ConvolutionComplex(Co,v)
        implicit none
        complex(rkind), intent(in)   ::v(:)
        complex(rkind), intent(out)  ::Co(:,:)
        integer                      ::i, j, k, n  
      
        Co = 0.0_rkind + ci * 0.0_rkind
    
		n  = size(v)
		do  i = 1, n
		do  k = 1, n
		     j= k-i+1
		     if (1 <= j .and. j <= n) then
			    Co(k,i) = Co(k,i) + v(j) 
		     endif
		     j = i-k+1
		     if (j <=n .and. j >= 1) then
				Co(k,i) = Co(k,i) + v(j) 
		     endif
		     if (k >1) then 
				j = i+k-1   
			    if (1 <= j .and. j <= n) then
					Co(k,i) = Co(k,i) + v(j)
				endif 
		     endif   
		     Co(k,i) = Co(k,i) * 0.5_rkind
		enddo
		enddo
    
    end subroutine ConvolutionComplex

end module cheb_mod

module linsys_mod
    use param_mod
    implicit none

contains

    subroutine MatrcInv(A, IA)
        implicit none
        complex(rkind)              :: A(:, :), IA(:, :)
        complex(rkind), allocatable :: B(:, :)
        integer                     :: m, i

        m = size(A, 1)
        allocate ( B(m, m) )
        IA = 0.0_rkind
        do i = 1, m
            IA(i, i) = 1.0_rkind
        end do
        B = A
        call MatrcUpper(B, IA)
        call MatrcLower(B, IA)
        do i = 1, m
            IA(i, :) = IA(i, :) / B(i, i)
        end do
        deallocate(B)
    end subroutine MatrcInv

    subroutine MatrcUpper(P, S)
        implicit none
        integer        :: m, i, j
        complex(rkind) :: P(:, :), S(:, :), E
        m = size(P, 1)
        do i = 1, m-1
        do j = i+1, m
            E = P(j, i) / P(i, i)
            P(j, i:m) = P(j, i:m) - P(i, i:m) * E
            S(j, :) = S(j, :) - S(i, :) * E
        end do
        end do
    end subroutine MatrcUpper

    subroutine MatrcLower(P, S)
        implicit none
        integer        :: m, i, j
        complex(rkind) :: P(:, :), S(:, :), E

        m = size(P, 1)
        do i = m, 2, -1
        do j = i-1, 1, -1
            E = P(j, i) / P(i, i)
            P(j, 1:m) = P(j, 1:m) - P(i, 1:m) * E
            S(j, :) = S(j, :) - S(i, :) * E
        end do
        end do
    end subroutine MatrcLower

end module linsys_mod

module util_mod
    use param_mod
    implicit none

contains

    subroutine assert (cond, msg)
        implicit none
        logical,          intent(in) :: cond
        character(len=*), intent(in) :: msg

        if (.not. cond) then
            write(*, *) "ERROR : ", msg
            stop
        end if
    end subroutine assert
end module util_mod

module CTSM_mod
    use param_mod
    use cheb_mod
    use linsys_mod
    implicit none

contains

    subroutine ReadEnvParameter(casename, N, cpmax, dr, zs, zr, dz, rmax, freq, H, Zground,&
        tlmin, tlmax, alpha, c, data_file)
        use util_mod
        implicit none
        character(len=MAX_FILENAME_LEN), intent(out) :: casename
		character(len=MAX_FILENAME_LEN), intent(in)  :: data_file
        integer,                         intent(out) :: N
        real(rkind),                     intent(out) :: cpmax
        real(rkind),                     intent(out) :: freq		
        real(rkind),                     intent(out) :: zs
        real(rkind),                     intent(out) :: zr
        real(rkind),                     intent(out) :: dz		
        real(rkind),                     intent(out) :: rmax
        real(rkind),                     intent(out) :: dr		
        real(rkind),                     intent(out) :: H
        real(rkind)                     			 :: Zg_re
        real(rkind)                     			 :: Zg_im		
        real(rkind),                     intent(out) :: tlmin
        real(rkind),                     intent(out) :: tlmax
        real(rkind), allocatable,        intent(out) :: alpha(:), c(:)
        real(rkind), allocatable, dimension(:)       :: temp_alpha, dep, temp_c
        integer										 :: nw,i
		complex(rkind),                  intent(out) :: Zground
                                        
        open(unit=1, status='unknown', file=data_file) 
    
        read (1, *) casename
        read (1, *) N
        read (1, *) cpmax
        read (1, *) freq
        read (1, *) zs
        read (1, *) zr
        read (1, *) dz		
        read (1, *) rmax
        read (1, *) dr
        read (1, *) H
        read (1, *) Zg_re
        read (1, *) Zg_im		
        read (1, *) tlmin
        read (1, *) tlmax
        read (1, *) nw
        
        !read the param_mod of ocean
        allocate(dep(nw), temp_c(nw), temp_alpha(nw))

        if (H > 0.0) then
            do i = 1, nw
                read(1, *) dep(i), temp_c(i), temp_alpha(i)
            end do
        else
            call assert(.false., "H must greater than 0!")
        end if
        close(1)

		Zground = cmplx(Zg_re, Zg_im, kind=rkind)

        call assert( N > 2, 'N must greater than 2!' )
        call assert( dep(1) == 0.0 .and. dep(nw) == H,  &
                     'input sound profile is unsuitable!' )
        call assert( H / dz == floor(H / dz), &
                     'The input dz unsuitable!' )
					 
        call assert( zs > 0 .and. zs < H .and. zr > 0 .and. zr < H,'bad zs!')
		
        call assert( rmax / dr == floor(rmax / dr), 'Please reinput the dr and rmax!' )

        !Interplating to the CGL points
        allocate(c(N+1), alpha(N+1))

        call Interpolation(dep,temp_c,c,temp_alpha,alpha,N,0.0_rkind,H)	

        deallocate(dep, temp_c, temp_alpha)

    end subroutine ReadEnvParameter

    subroutine Interpolation(dep,b1,b2,d1,d2,N,s,t)
        implicit none
        integer,       intent(in)   ::N
        real(rkind),   intent(in)   ::dep(:), b1(:), d1(:), s, t
        real(rkind),   intent(out)  ::b2(:), d2(:)
        real(rkind)                 ::x(N+1), z(N+1)
        integer                     ::i, j, m
        
        m = size(dep)
        do i = 1, N + 1
            x(i) = cos((i - 1) * pi / N)
            z(i) = ((t + s) / (t - s) - x(i)) * (t - s) / 2.0
        end do

        do i=1,N+1	
            do j=1,m-1
                if((z(i) >= dep(j)).and.(z(i) <= dep(j+1))) then
                    b2(i) = (z(i) - dep(j)) / (dep(j+1) - dep(j)) * b1(j+1)&
                            +(dep(j+1)-z(i))/ (dep(j+1) - dep(j)) * b1(j)
                    
                    d2(i) = (z(i) - dep(j)) / (dep(j+1) - dep(j)) * d1(j+1)&
                            +(dep(j+1)-z(i))/ (dep(j+1) - dep(j)) * d1(j)				
                    
                endif				
            enddo 
            
        if(z(i) > dep(m)) then
            b2(i) = b1(m)
            d2(i) = d1(m)							
        end if
            
        if(z(i) < dep(1)) then
            b2(i) = b1(1)
            d2(i) = d1(1)					
            endif
        enddo	
        
    end subroutine Interpolation

    subroutine Initialization(N, freq, rmax, dr, c, alpha, H, nr, r, k)
        implicit none
        integer,                     intent(in)  :: N
        integer,                     intent(out) :: nr
        real(rkind),                 intent(in)  :: freq
        real(rkind),                 intent(in)  :: rmax
        real(rkind),                 intent(in)  :: dr
        real(rkind),                 intent(in)  :: H
        real(rkind), dimension(N+1), intent(in)  :: c, alpha
        real(rkind),    allocatable, intent(out) :: r(:)
        complex(rkind), allocatable, intent(out) :: k(:)
        real(rkind), dimension(N+1)              :: x
        real(rkind)                              :: w
        integer                                  :: i

        w  = 2.0_rkind * pi * freq
        nr = int(rmax / dr)

        allocate ( r(nr) )
		
        do i = 1, nr
            r(i) = i * dr
        end do

        allocate (k(N+1))
        k = w / c * (1.0_rkind + ci * alpha / (40.0_rkind * pi * log10(exp(1.0_rkind))))
        
    end subroutine Initialization

    subroutine EigenValueVector(N, H, k, Zground, kr, eigvector)
        implicit none
        integer,                     intent(in)    :: N
        real(rkind),                 intent(in)    :: H
        complex(rkind),              intent(inout) :: k(N+1)
        complex(rkind), allocatable, intent(out)   :: kr(:)
        complex(rkind), allocatable, intent(out)   :: eigvector(:, :)
        complex(rkind), dimension(N+1)    	       :: Pg, Pa
        real(rkind),    dimension(N+1, N+1)        :: D		
        real(rkind)                                :: x(N+1), ub(N+1), mb(N+1)		
        complex(rkind)                             :: A(N+1, N+1), tem1(N+1), tem2(N+1, N+1), Zground	
        complex(rkind)                             :: L11(N-1, N-1)
        complex(rkind)                             :: L12(N-1,   2)
        complex(rkind)                             :: L21(2,   N-1)
        complex(rkind)                             :: L22(2,     2)
        complex(rkind)                             :: L  (N-1, N-1)
        complex(rkind)                             :: IL22(2, 2)
        complex(rkind)                             :: v2(2, N-1)
        complex(rkind)                             :: VL(N-1)
        complex(rkind)                             :: VR(N-1, N-1)
        complex(rkind)                             :: WORK (2*(N-1))
        real(rkind)                                :: RWORK(2*(N-1))
        integer                                    :: i, info, j(1)

        allocate(kr(N-1), eigvector(N+1, N-1))

        call DerivationMatrix(D, N+1)

        do i = 1, N + 1
            x(i) = cos( (i-1) * pi / N )
			ub(i)= 1.0_rkind
			mb(i)= (-1.0_rkind) ** (i-1)
        end do
	
		Pg = -2.0_rkind / H * matmul(ub, D) + ci * k(1) / Zground * ub
		Pa = -2.0_rkind / H * matmul(mb, D) - ci * k(N+1) * mb
		
		k = k ** 2
		A = 4.0_rkind / H / H * matmul(D, D)
        call Cheb(tem1, k, x)
        call Convolution(tem2, tem1)
        A = A + tem2
		
		A(N,  :) = Pg
		A(N+1,:) = Pa

        L11 = A(1:N-1, 1:N-1)
        L12 = A(1:N-1, N:N+1)
        L21 = A(N:N+1, 1:N-1)
        L22 = A(N:N+1, N:N+1)

        call MatrcInv(L22, IL22)
        L = L11 - matmul(L12, matmul(IL22, L21))

        call zgeev('N', 'V', N-1, L, N-1, kr, VL, 1, VR, N-1, WORK, 2*(N-1), RWORK, INFO)

        v2 = - matmul(matmul(IL22, L21), VR)

        eigvector = 0.0_rkind
        eigvector(1:N-1, :) = VR
        eigvector(N:N+1, :) = v2

        kr = sqrt(kr)
        !L stores the sorted eigenvectors. VL stores the sorted eigenvalues
        do i = 1, N-1
            j = maxloc(real(kr))
            VL(i)    = kr(j(1))
            A(:, i)  = eigvector(:, j(1))
            kr(j(1)) = - 1.0_rkind
        end do

        kr = VL
        eigvector = A(:,1:N-1)

    end subroutine EigenValueVector

    subroutine NumOfModes(freq, kr, nmodes, cpmax)
        use util_mod
        implicit none
        integer,        intent(out):: nmodes
        real(rkind),    intent(in) :: freq
        real(rkind),    intent(in) :: cpmax
        complex(rkind), intent(in) :: kr(:)
        real(rkind), allocatable   :: cp(:)
        integer                    :: i

        allocate ( cp(size(kr)) )
        cp = 2.0_rkind * pi * freq / real(kr, kind=rkind)
        nmodes = 0
        do i = 1, size(kr)
            if (cp(i) <= cpmax) nmodes = i
        end do

        call assert( nmodes /= 0, 'Incorrect maximum phase speed input!')

        deallocate(cp)
    end subroutine NumofModes

    subroutine GenerateModes(nmodes, dz, H, eigvector, psi, z)
        implicit none
        integer,                         intent(in)  :: nmodes
        real(rkind),                     intent(in)  :: dz, H
        complex(rkind), dimension(:, :), intent(in)  :: eigvector
        real(rkind),    allocatable,     intent(out) :: z(:)
        complex(rkind), allocatable,     intent(out) :: psi(:, :)
        real(rkind),    allocatable      		     :: zt(:)
        integer                                      :: i

        allocate ( z(int(H / dz)), zt(int(H / dz)) )

        do i = 1, int(H / dz) + 1
            z (i) = (i - 1) * dz
			zt(i) = - 2.0_rkind / H * z(i) + 1.0_rkind
        end do

        allocate ( psi(size(z), nmodes) )
        call InvChebMatrix(psi, eigvector(:, 1:nmodes), zt)
		
        deallocate(zt)
    end subroutine GenerateModes

    subroutine Normalization(eigvector, H, N, nmodes, psi)
        implicit none
        complex(rkind), dimension(:, :), intent(in)  :: eigvector
        complex(rkind), dimension(:, :), intent(out) :: psi

        integer            :: i, N, nmodes,k
        real(rkind)        :: H, x(N+1), P(N+1)
        complex(rkind)     :: Co1(N+1, N+1), f(N+1), norm,temp(N+1,nmodes)

        do i = 1, N+1
            x(i) = cos( (i-1) * pi / N )
        end do
!---------------------------------------------------------------------------
!        P = 0.0_rkind
!        do i = 0, N, 2
!            P(i+1) = - 2.0_rkind / (i * i - 1.0_rkind)
!        end do

!        do i = 1, nmodes
!            call Convolution(Co1, eigvector(:, i))
!            f = matmul(Co1, eigvector(:, i))
!            norm = sqrt(dot_product(P, f) * H * 0.5_rkind)
!			 psi(:, i) = psi(:, i) / norm
!        end do

!----------------------------------------------------------------------------
        call InvChebMatrix(temp, eigvector(:, 1:nmodes), x)
		do i = 1, nmodes		
			f = temp(:, i)
			forall(k=1:N+1) f(k) = f(k) * f(k)			
			norm = ChebGaussQuadrature(f) * H * 0.5_rkind			
			psi(:,i) = psi(:,i) / sqrt(norm)	
		enddo

    end subroutine Normalization
	
	function ChebGaussQuadrature(f)
		implicit none
		complex(rkind),intent(in)::f(:)
		complex(rkind) ChebGaussQuadrature
		real(rkind),allocatable ::w(:),x(:)
		integer(rkind) n,k
		
		n  = size(f);
		allocate(w(n),x(n))
		w=pi/(n-1.0_rkind)
		w(1)=w(1)*0.5_rkind
		w(n)=w(n)*0.5_rkind
		forall(k=1:n)  x(k)=cos((k-1)*pi/(n-1))
		
		ChebGaussQuadrature = 0.0_rkind;

		do k = 1, n
			 ChebGaussQuadrature = ChebGaussQuadrature + w(k) * f(k) * (sqrt(1 - x(k)**2))
		enddo

		deallocate(w,x)

	end function ChebGaussQuadrature	
	
    subroutine SynthesizeSoundField(nmodes, nr, r, kr, zs, z, dz, psi, tl)
        implicit none
        integer,                  intent(in)    :: nmodes, nr
        real(rkind),              intent(in)    :: r(nr), zs, z(:), dz
        real(rkind), allocatable, intent(inout) :: tl(:, :)
        complex(rkind),           intent(in)    :: kr(:)
        complex(rkind),           intent(inout) :: psi(:, :)

        complex(rkind), allocatable             :: p(:, :)
        complex(rkind)                          :: bessel(nmodes, nr)
        complex(rkind)                          :: psizs(nmodes, nmodes)
        integer                                 :: i, k

        do k = 1, nmodes
        do i = 1, nr
            bessel(k, i) = exp(ci * r(i) * kr(k)) / ( sqrt( r(i) * kr(k) ) )
        end do
        end do
            
        allocate(p(size(psi, 1), nr), tl(size(psi, 1), nr))

        psizs(:, :)=0.0_rkind
        do k = 1, nmodes
            call Interpolation_zs(z,zs,psi(:,k),psizs(k,k))
        end do
        
        psi=matmul(psi, psizs)
        p = matmul(psi, bessel)
        p = p * sqrt( 2 * pi)
        
        tl = - 20.0_rkind * log10(abs(p))
    end subroutine SynthesizeSoundField

	subroutine Interpolation_zs(z,zs,v,p)
		implicit none
		real(rkind),   intent(in) ::z(:),zs
		complex(rkind),intent(in) ::v(:)
		complex(rkind),intent(out)::p
		integer :: j

		do j=1,size(z)-1
			if (zs>=z(j).and.zs<=z(j+1)) then
				p=(zs-z(j))/(z(j+1)-z(j))*v(j+1) +(z(j+1)-zs)/(z(j+1)-z(j))*v(j)
			endif				
		enddo 

	end subroutine

    subroutine SaveSoundField(filename, tlmin, tlmax, r, z, tl)
        implicit none
        character(len=MAX_FILENAME_LEN), intent(in) :: filename
        real(rkind),                     intent(in) :: tlmin, tlmax
        real(rkind), dimension(:),       intent(in) :: r, z
        real(rkind), dimension(:, :),    intent(in) :: tl

        open(unit=20, status='unknown', file=filename, access='stream', form='unformatted')
        write(20)  size(z), size(r), tlmin, tlmax, z, r, tl
        close(20)
    end subroutine SaveSoundField	
	

end module CTSM_mod
