module param_mod
    implicit none
    integer,        parameter :: rkind            = 8
    integer,        parameter :: MAX_FILENAME_LEN = 200

    real(rkind),    parameter :: pi = 4.0_rkind * atan(1.0_rkind)   ! 双精度pi值
    complex(rkind), parameter :: ci = cmplx(0.0, 1.0, kind=rkind)   ! 复数
end module

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

module CCSM_mod
    use param_mod
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
            x(i) = - cos((i - 1) * pi / N)
            z(i) = ((t + s) / (t - s) + x(i)) * (t - s) / 2.0
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

    subroutine Initialization(N, freq, rmax, dr, c, alpha, z, H, nr, r, k)
        implicit none
        integer,                     intent(in)  :: N
        integer,                     intent(out) :: nr
        real(rkind),                 intent(in)  :: freq
        real(rkind),                 intent(in)  :: rmax
        real(rkind),                 intent(in)  :: dr
        real(rkind),                 intent(in)  :: H
		real(rkind),    allocatable, intent(out) :: z(:)		
        real(rkind), dimension(N+1), intent(in)  :: c, alpha
        real(rkind),    allocatable, intent(out) :: r(:)
        complex(rkind), allocatable, intent(out) :: k(:)
        real(rkind), dimension(N+1)              :: x
        real(rkind)                              :: w
        integer                                  :: i

        w  = 2.0_rkind * pi * freq
        nr = int(rmax / dr)

        allocate ( r(nr) , z(N+1))
		
        do i = 1, nr
            r(i) = i * dr
        end do
		
        do i = 1, N+1
            x(i) = - cos( (i-1) * pi / N )
			z(i) = H / 2 * x(i) + H / 2
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
        real(rkind),    dimension(N+1, N+1)        :: D		
        real(rkind)                                :: x(N+1)	
        complex(rkind)                             :: A(N+1, N+1), B(N+1,N+1), Zground	
        complex(rkind)                             :: L11(N-1, N-1)
        complex(rkind)                             :: L12(N-1, 2)
        complex(rkind)                             :: L21(2, N-1)
        complex(rkind)                             :: L22(2, 2)
        complex(rkind)                             :: L(N-1, N-1)
        complex(rkind)                             :: IL22(2, 2)
        complex(rkind)                             :: v2(2, N-1)
        complex(rkind)                             :: VL(N-1)
        complex(rkind)                             :: VR(N-1, N-1)
        complex(rkind)                             :: WORK (2*(N-1))
        real(rkind)                                :: RWORK(2*(N-1))
        integer                                    :: i, info, j(1)

        allocate(kr(N-1), eigvector(N+1, N-1))

        do i = 1, N + 1
            x(i) = - cos( (i-1) * pi / N )
        end do

        call DerivationMatrix(D, x)

        A = 4.0_rkind / H / H * matmul(D, D)
        do i = 1, N + 1
			A(i, i) = A(i, i) + k(i) ** 2
        end do

		A(1,  :)    = 2.0_rkind / H * D(1, :)
		A(1,  1)    = A(1, 1) + ci * k(1) / Zground
		A(N+1,:)    = 2.0_rkind / H * D(N + 1, :)
		A(N+1, N+1) = A(N + 1, N + 1) - ci * k(N + 1);
		
		B(1:N-1,:) = A(2:N, :)
		B(N    ,:) = A(1,   :)
		B(N+1,  :) = A(N+1, :)		
		A(:,1:N-1) = B(:, 2:N)
		A(:,    N) = B(:,   1)
		A(:,  N+1) = B(:, N+1)
		

        L11 = A(1:N-1, 1:N-1)
        L12 = A(1:N-1, N:N+1)
        L21 = A(N:N+1, 1:N-1)
        L22 = A(N:N+1, N:N+1)

        call MatrcInv(L22, IL22)
        L = L11 - matmul(L12, matmul(IL22, L21))

        call zgeev('N', 'V', N-1, L, N-1, kr, VL, 1, VR, N-1, WORK, 2*(N-1), RWORK, INFO)

        v2 = - matmul(matmul(IL22, L21), VR)

        eigvector = 0.0_rkind
        eigvector(2:N, :)  = VR
		eigvector(1,   :)  = v2(1,:)
        eigvector(N+1, :)  = v2(2,:)

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
	
	subroutine DerivationMatrix(D, x)
		implicit none
		real(rkind),intent(out)::D(:,:),x(:)
		real(rkind),allocatable::c(:)
		real(rkind) :: v,summ
		integer :: i,j,m
		m=size(x)
		allocate(c(m))
		D=0.0_rkind
		c=1.0_rkind
		c(1)=2.0_rkind
		c(m)=2.0_rkind
		c(2:m:2)=-1.0_rkind

		do i=1,m
			summ=0.0_rkind
			do j=1,m
				if(i/=j) then
				   v=c(i)/c(j)/(x(i)-x(j))
				   D(i,j)=v
				   summ=summ+v
				endif
			enddo
			D(i,i)=-summ
		enddo

	end subroutine

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

	subroutine Normalization(nmodes,N,eigvector,psi,z,H,zs,psizs)
		implicit none
		integer,intent(in)::nmodes, N
		real(rkind),   intent(in)::zs, H, z(N+1)
		complex(rkind),intent(out),allocatable::psi(:,:),psizs(:)
		complex(rkind),intent(in)::eigvector(N+1,N-1)
		complex(rkind) :: f(N+1),norm
		integer:: i, k
		
		allocate(psi(N+1,nmodes), psizs(nmodes))

		psi = eigvector(:, 1:nmodes)

		do i = 1, nmodes		
			f = eigvector(:, i)
			forall(k=1:N+1) f(k) = f(k) * f(k)			
			norm = ChebGaussQuadrature(f) * H * 0.5_rkind			
			psi(:,i) = psi(:,i) / sqrt(norm)	
			call Interpolation_zs(z,zs,psi(:,i),psizs(i))
		enddo

	end	

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

    subroutine SynthesizeSoundField(nmodes, nr, r, kr, zs, dz, psi, psizs, tl)
        implicit none
        integer,                  intent(in)    :: nmodes, nr
        real(rkind),              intent(in)    :: r(nr), zs, dz
        real(rkind), allocatable, intent(inout) :: tl(:, :)
        complex(rkind),           intent(in)    :: kr(:)
        complex(rkind),           intent(inout) :: psi(:, :)

        complex(rkind), allocatable             :: p(:, :)
        complex(rkind)                          :: bessel(nmodes, nr)
        complex(rkind)                          :: psizs(nmodes),psizs2(nmodes,nmodes)
        integer                                 :: i, k

		psizs2=0.0_rkind
        do k = 1, nmodes
			do i = 1, nr
				bessel(k, i) = exp(ci * r(i) * kr(k)) / ( sqrt( r(i) * kr(k) ) )
			end do
			psizs2(k, k) = psizs(k)
        end do
            
        allocate(p(size(psi,1), nr), tl(size(psi,1), nr))
        
        psi=matmul(psi, psizs2)
        p = matmul(psi, bessel)
        p = p * sqrt( 2 * pi)
        
        tl = - 20.0_rkind * log10(abs(p))
    end subroutine SynthesizeSoundField

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
	

end module CCSM_mod
