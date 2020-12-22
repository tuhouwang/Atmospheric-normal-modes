
program AtmosCTSM
    use param_mod
    use CTSM_mod
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
    complex(rkind),allocatable, dimension(:)    :: k, kr
    complex(rkind),allocatable, dimension(:, :) :: eigvector, psi
    real(rkind),   allocatable, dimension(:)    :: z
    real(rkind),   allocatable, dimension(:, :) :: tl
    complex(rkind)  						    :: Zground
	
    call ReadEnvParameter(casename, N, cpmax, dr, zs, zr, dz, rmax, freq, H, Zground,&
        tlmin, tlmax, alpha, c, data_file)		
    
	call Initialization(N, freq, rmax, dr, c, alpha, H, nr, r, k)
    
    call EigenValueVector(N, H, k, Zground, kr, eigvector)

    call NumOfModes(freq,kr,nmodes,cpmax)
    
    call GenerateModes(nmodes, dz, H, eigvector, psi, z)

    call Normalization(eigvector, H, N, nmodes, psi)
		
    call SynthesizeSoundField(nmodes, nr, r, kr, zs, z, dz, psi, tl)	
	
    call SaveSoundField(filename,tlmin,tlmax,r,z,tl)

    deallocate(alpha,c,r,k,kr,eigvector,psi,z,tl)	

end program AtmosCTSM
