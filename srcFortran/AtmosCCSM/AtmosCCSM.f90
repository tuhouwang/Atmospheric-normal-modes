
program AtmosCCSM
    use param_mod
    use CCSM_mod
    implicit none
    
    external zgeev
    !---------------------------------------------------------------------------------
    ! Declare the variable needed later.
    character(len=MAX_FILENAME_LEN)             :: casename
    character(len=MAX_FILENAME_LEN)             :: data_file = "input_CAA.txt"
    character(len=MAX_FILENAME_LEN)             :: filename  = "tl.bin"
    integer                                     :: N
    integer                                     :: nr
    integer                                     :: nmodes
    real(rkind)                                 :: cpmax
    real(rkind)                                 :: dr
    real(rkind)                                 :: zs
    real(rkind)                                 :: zr
    real(rkind)                                 :: rmax
    real(rkind)                                 :: freq
    real(rkind)                                 :: H
    real(rkind)                                 :: dz
    real(rkind)                                 :: tlmin
    real(rkind)                                 :: tlmax
    real(rkind),   allocatable, dimension(:)    :: alpha, c
    real(rkind),   allocatable, dimension(:)    :: r
    complex(rkind),allocatable, dimension(:)    :: k, kr, psizs
    complex(rkind),allocatable, dimension(:, :) :: eigvector, psi
    real(rkind),   allocatable, dimension(:)    :: z
    real(rkind),   allocatable, dimension(:, :) :: tl
    complex(rkind)  						    :: Zground
	
    call ReadEnvParameter(casename, N, cpmax, dr, zs, zr, dz, rmax, freq, H, Zground,&
        tlmin, tlmax, alpha, c, data_file)		
    
	call Initialization(N, freq, rmax, dr, c, alpha, z, H, nr, r, k)
    
    call EigenValueVector(N, H, k, Zground, kr, eigvector)

    call NumOfModes(freq,kr,nmodes,cpmax)

    call Normalization(nmodes,N,eigvector,psi,z,H,zs,psizs)

    call SynthesizeSoundField(nmodes, nr, r, kr, zs, dz, psi, psizs, tl)	
	
    call SaveSoundField(filename,tlmin,tlmax,r,z,tl)

    deallocate(alpha,c,r,k,kr,psizs,eigvector,psi,z,tl)	

end program AtmosCCSM
