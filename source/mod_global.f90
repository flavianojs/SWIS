!Program to calculate the spinwaves disperion of non-collinear magnetic structure. It allows to account for the symmetric exchange interaction the DMI, external mag. fields and anisotropy.
!Flaviano J. dos Santos, April 2016
!Forschungszentrum Jülich, Germany
module mod_global
   implicit none
   !!!!!!!! Globals, constants and non-ajustables !!!!!!!!!!!!!!!!!!!!!!
   integer, parameter :: pc=8 !Precision of the calculations
   complex(kind=pc), parameter :: ii=(0.d0,1.d0), czero=(0.d0,0.d0), cone=(1.d0,0.d0), ctwo=(2.d0,0.d0)
   real(kind=pc), parameter :: pi=2.d0*asin(1.d0)
   !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
   

   !!!!!!!!!!!!!!!!! Simulation control board (Parameters) !!!!!!!!!!!!!
   character(len=100) :: inputcardname
   character(len=100) :: basisname, latticefile
   integer :: dims(3) !Determines the dimension of the problem. Shall be filled by 1 or 0
   integer :: n_basis_cells(3) = [1.d0,1.d0,1.d0] !To expand the unit cell to form the simulation box
   integer :: npt, naucell, ncellpdim=2, npath, ncpdim ! n. of k-point in the disperion; n. of atom in the unit cell; number of real space unit cell for the Fourier transformation
   integer :: twonaucell, effnkpt, nptomega=0, n_nndists=0
   integer :: max_num_kpt_path=50, max_num_ani=100
   integer, allocatable :: nkpt_npath(:)
   real(kind=pc) :: latcons
   real(kind=pc) :: a1(3), big_a1(3), b1(3)
   real(kind=pc) :: a2(3), big_a2(3), b2(3)
   real(kind=pc) :: a3(3), big_a3(3), b3(3)
   real(kind=pc) :: polarization1(2) = [0.d0,0.d0], polarization2(2) = [0.d0,0.d0], polarization3(2) = [0.d0,0.d0]

   !The Interactions
   real(kind=pc) :: Jnn=0.0d0, all_Jnn(20)=0.d0, all_nndist(20)=0.d0, Dnn=0.d0, kanis=0.d0, kanis2=0.d0, hm0=0.d0, muB=1.d0, interactions_cutoff_radius = 1.d6
   real(kind=pc) :: hmagunitvec(3), kaniunitvec(3) = [0.d0, 0.d0, 1.d0], kaniunitvec2(3) = [0.d0, 0.d0, 1.d0], J_D_scaling=0.123456789d0
   real(kind=pc) :: mu_s=1.d0, gamma=1.d0 ! mag mom equals to gamma times the spin operator: mu_s = gamma * S
   real(kind=pc) :: mu_s_vec(1000) = 0.123456789d0

   !Calculations to be performed
   logical :: Tneel=.false.
   logical :: spirit=.false.
   logical :: printD=.false.
   logical :: toprint=.false.
   logical :: unfolding=.false.
   logical :: analytics=.false.
   logical :: calc_occup=.false. !Occup. number is calculated on the first k point of the k-path
   logical :: toprint_Rij=.false.
   logical :: kpoint_mesh=.false.
   logical :: internalunit=.false.
   logical :: neutron_factor=.false. !Neutron scattering polarisation factor which ensures only components of spin perpendicular to k are observed
   logical :: spirit_old_method=.false.
   logical :: atoms_distinguisable=.false.
   logical :: constant_energy_plot=.false.
   logical :: circular_cartesian_convertion_factor=.true. !Add the convertion term from circular to catersian spin compenents when computing the neutron scattering cross section
   integer :: mode=1
   integer :: spin_dyn = 1
   !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!


   !!!!!!!!!!!!!!!! Auto determined global variables !!!!!!!!!!!!!!!!!!!
   real(kind=pc) :: Jx, Jy, Jz, maxomega=0.123456789d0, minomega=0.123456789d0, eta=1.d-1 , rescalf=1.1d0, zero_toler=1.d-5
   real(kind=pc) :: Dx, Dy, Dz, hmag(3), bb(3,3), nndist=0.d0
   real(kind=pc), allocatable :: r0j(:,:), anglesphi(:), anglestheta(:), Si(:), Si_aux(:), basis(:,:) , kpoints(:,:), dkpoints(:), symmpts(:,:), syptsMataux(:,:), syptsMat(:,:), Rotmat(:,:,:), Uweight(:), mu_s_vec_aux(:), n_anisotropy_aux(:,:), n_anisotropy(:,:), anisotropy_vecs(:,:)
   integer :: ninteraction_pairs, ncell !Number of unit cells
   integer :: origcell !The index of the origin cell
   integer :: num_anisotropies = 0
   complex(kind=pc), allocatable :: highsymm(:,:), Rotpm(:,:,:), g(:,:)
   complex(kind=pc) :: paulimatrix(3,3,2,2), paulimatrix_cartesian(3,2,2), paulimatrix_cartesian_rotated(3,2,2)
   real(kind=pc) :: beg_cpu_time, positions(100000,3), DJvect(100000,4)
   real(kind=pc) ::  kaniunitvec3(3), kanis3
   integer :: ijda_db_dc(100000,5)
   logical :: maxomegaOK
   !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!

contains

   ! Returns the inverse of a matrix calculated by finding the LU
   ! decomposition.  Depends on LAPACK.
   function inv(A) result(Ainv)
      implicit none
      real(pc), dimension(:,:), intent(in) :: A
      real(pc), dimension(size(A,1),size(A,2)) :: Ainv

      real(pc), dimension(size(A,1)) :: work  ! work array for LAPACK
      integer, dimension(size(A,1)) :: ipiv   ! pivot indices
      integer :: n, info

      ! External procedures defined in LAPACK
      external DGETRF
      external DGETRI

      ! Store A in Ainv to prevent it from being overwritten by LAPACK
      Ainv = A
      n = size(A,1)

      ! DGETRF computes an LU factorization of a general M-by-N matrix A
      ! using partial pivoting with row interchanges.
      call DGETRF(n, n, Ainv, n, ipiv, info)

      if (info /= 0) then
         stop 'Matrix is numerically singular!'
      end if

      ! DGETRI computes the inverse of a matrix using the LU factorization
      ! computed by DGETRF.
      call DGETRI(n, Ainv, n, ipiv, work, n, info)

      if (info /= 0) then
         stop 'Matrix inversion failed!'
      end if
   end function inv

!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
   subroutine invers(matriz)
      implicit none
      integer :: lwork 
      complex(kind=pc), dimension(:,:), intent(in) :: matriz
      complex(kind=pc), dimension(size(matriz,1)*4) :: work  
      integer, dimension(size(matriz,1)) :: ipiv
      integer :: nn,info

      nn = size(matriz,1)
      lwork = 4*nn
      info = 0

      call zgetrf(nn,nn,matriz,nn,ipiv,info)
      if (info/=0) then
         write(*,*)'ifail=',info
      end if
      call zgetri(nn,matriz,nn,ipiv,work,lwork,info)

      if (info/=0) then
         write(*,*)'ifail=',info
         stop 'Fudeu!!! Inversao da NAG foi pro caralho!'
      end if

      return
   end subroutine invers

!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
   subroutine initialization()
      implicit none
      real(kind=pc), parameter :: c0(3) = [0.123456789, 0.123456789, 0.123456789]

      integer :: i,j,k, dim, counter, reading_status, d1, d2, d3, l1, l2, ierr, num_atoms_small_basis, polariz
      real(kind=pc) :: basisaux(3), amat(3,3), bmat(3,3), net_magnetization(3), Sx, Sy, Sz, pos_norm, big_a1_norm, big_a2_norm, big_a3_norm, Rotmat_aux(3,3), polarizations(3,2)
      real(kind=pc), allocatable :: unit_cell_basis(:,:)
      real(kind=pc) :: c1(3)=c0, c2(3)=c0, c3(3)=c0, position(3)
      complex(kind=pc) :: HP(3,3), HPinv(3,3), aux(3,3)
      character(len=1000) :: read_in_data, lines, spirit_input = "", pairfile
      character(len=100) :: word
      logical :: badindexing, flag1, file_exists
      logical :: spirit_input_given = .false., c1c2c3set_given = .false.

      ! c0 = [0.123456789, 0.123456789, 0.123456789]
      ! c1 = c0
      ! c2 = c0
      ! c3 = c0

      ! A list to be read on the input file: inputcardsky.f90
      namelist /input/ &
      !===== crystal strucutre ===========================================================
         dims, a1, a2, a3, latcons, naucell, ncellpdim, c1, c2, c3, internalunit, &
      !===== k path and energy sampling ==================================================
         npt, nptomega, npath, maxomega, minomega, eta, ncpdim, &
         polarization1, polarization2, polarization3, &
      !===== inputing and controlling interactions =======================================
         Jnn, J_D_scaling, Dnn, nndist, all_Jnn, all_nndist, interactions_cutoff_radius, &
         hm0, muB, hmagunitvec, &
         kanis, kanis2, kaniunitvec, kaniunitvec2, &
         mu_s, gamma, mu_s_vec, &
      !===== logical var: turn on and off features =======================================
         spirit, spirit_input, unfolding, analytics, calc_occup, constant_energy_plot, & 
         toprint, printD, toprint_Rij, spirit_old_method, atoms_distinguisable, &
         kpoint_mesh, Tneel, neutron_factor, circular_cartesian_convertion_factor, &
      !===== others ======================================================================
         mode, spin_dyn, zero_toler

      namelist /input2/ Uweight, syptsMataux, basisname, pairfile, latticefile, n_anisotropy_aux
         
      ! Reading the inputs and allocating variables   
      inquire(file=inputcardname, exist=file_exists)
      if(file_exists) then
         open(unit=101, file=inputcardname, status="old")
         print "(a,a)", " Inputcard file: ", inputcardname
      else
         print "(3a)", " File: '", trim(inputcardname), "' does NOT exist. Stopping."
         print *
         stop 
      end if

      read(101, nml=input)
      !write(*,nml=input)

      if(neutron_factor) then
         print "(a)", " Neutron scattering polarization factor is ON."
      end if
      if(circular_cartesian_convertion_factor) then
         print "(a)", " Sum convertion from Circular to Cartesian on inelastic cross section is ON."
      end if

      !Reading the magnetic moment vector
      allocate( mu_s_vec_aux(naucell) )

      ! Search for the location of the first element equal to 0.123456789d0, 
      i = FINDLOC(mu_s_vec, 0.123456789d0, 1)

      if( i < naucell + 1 ) then
         mu_s_vec_aux = mu_s
         
         if( i > 1 ) then
            print "(3a)", " In file: '", trim(inputcardname), "'. ATTENTION: Number of magmetic moment does not match number of sites. Using 'mu_s' for all sites."
         end if
      else
         mu_s_vec_aux = mu_s_vec(1:naucell) 
      end if
      print "(a,1000f8.4)", " Mag moms= ", mu_s_vec_aux

      !If to compute spectrum of a energy cut, we use a kpoint mesh instead of a kpoint path.
      if( constant_energy_plot ) kpoint_mesh = .true.

      !Determining the number of nn interactions inputed
      do i=1, 20
         if( abs(all_nndist(i)) < zero_toler ) then
            n_nndists = i - 1
            exit
         end if
      end do
      if( n_nndists .ne. 0) print *, "Interactions will be inputed from all_Jnn by measuring the n.n distances."
      if( n_nndists .ne. 0) print *, "n_nndists:", n_nndists

      if( spirit_input .ne. "" ) then
         spirit_input_given = .true.
         spirit_old_method = .true.
         print *
         print "(2a)", "Reading SPIRIT inputcard: ", trim(spirit_input)

         inquire(file=spirit_input, exist=file_exists)
         if(file_exists) then
            open(unit=102, file=spirit_input)
         else
            print "(3a)", " File: '", trim(spirit_input), "' does NOT exist. Stopping."
            print *
            stop 
         end if

         flag1 = .false.
         num_atoms_small_basis = -1
         counter = -1
         do 
            read(unit=102, fmt="(a)", iostat=ierr ) lines
            if(ierr /= 0) exit
            if(lines == "" ) cycle

            read(lines, *) word

            !______________________________________________________
            if( word == "n_basis_cells") then
               read(lines, *) word, n_basis_cells
               naucell = n_basis_cells(1)*n_basis_cells(2)*n_basis_cells(3)
               print "(a,3i4)", " n_basis_cells ", n_basis_cells
               print "(a,3i4)", " naucell ", naucell
            end if
            !-------------------------------------------------------

            !______________________________________________________
            if( 0 < counter .and. counter <= num_atoms_small_basis ) then
               read(lines, *) basisaux
               unit_cell_basis(counter,:) = basisaux
               print "(3f12.8)", basisaux 
               counter = counter + 1
            end if

            if( flag1 ) then
               read(lines, *) num_atoms_small_basis
               print "(a)", " basis"  
               print "(a,i0)", " ", num_atoms_small_basis  
               allocate( unit_cell_basis(num_atoms_small_basis,3) )
               flag1 = .false.
               counter = 1
            end if

            if( word == "basis") flag1 = .true.
            !-------------------------------------------------------

         end do
         print "(a)", "END Reading SPIRIT inputcard."
         print *
      end if

      twonaucell = 2*naucell
      allocate( syptsMataux(3,max_num_kpt_path+2), g(twonaucell,twonaucell), Uweight(naucell), n_anisotropy_aux(5,max_num_ani) )

      Uweight = 1.d0 !Default value
      syptsMataux = 0.123456789d0
      n_anisotropy_aux = 0.123456789d0

      read(101, nml=input2)

      do i = 1, max_num_kpt_path+2
         if( syptsMataux(1,i) == 0.123456789d0 ) then 
            npath = i - 2
            print "(a,i0)", " Attributing 'npath' = ", npath
            exit
         end if
         if( i == max_num_kpt_path+2 ) then
            print *, "It is better that you make 'max_num_kpt_path' larger than", max_num_kpt_path, "Stopping."
            stop
         end if
      end do
      allocate( syptsMat(npath+1,3) )

      do i = 1, max_num_ani
         if( n_anisotropy_aux(1,i) == 0.123456789d0 ) then 
            num_anisotropies = i-1
            print "(a,i0,a)", " Number of anisotropy lines = ", num_anisotropies
            exit
         end if
         if( i == max_num_ani ) then
            print *, "It is better that you make 'max_num_ani' larger than", max_num_ani, "Stopping."
            stop
         end if
      end do
      allocate( n_anisotropy(num_anisotropies,5) )
      n_anisotropy = transpose(n_anisotropy_aux(:,1:num_anisotropies))
      print *, "  i      Kx             Ky             Kz             K"
      do j=1, num_anisotropies
         print "(i4,' '4f15.10)", int(n_anisotropy(j,1)), n_anisotropy(j,2:)  
      end do

      if(nndist==0.d0) then
         print *, "No n.n distance inputted! Setting it to its default value: nndist=latcons"
         nndist = latcons
      end if

      !g is a diagonal matrix of first upper half entries equal 1 and second half equal to -1. Diagonalization of the D dynamical matrix and on the Bogoliubov transformation. 
      g(:,:) = 0.d0
      do i=1, naucell
         g(i,i) = 1.d0
         g(i+naucell,i+naucell) = -1.d0
      end do

      syptsMat = transpose(syptsMataux(:,1:npath+1))

      !make sure we have a unitarian vector
      kaniunitvec = kaniunitvec / norm2(kaniunitvec)
      kaniunitvec = kaniunitvec / norm2(kaniunitvec)
      kaniunitvec = kaniunitvec * kanis * gamma**2/(mu_s**2)

      kaniunitvec2 = kaniunitvec2 / norm2(kaniunitvec2)
      kaniunitvec2 = kaniunitvec2 * kanis2 * gamma**2/(mu_s**2)

      allocate(anisotropy_vecs(naucell,3))

      do i=1, naucell
         anisotropy_vecs(i,:) = kaniunitvec + kaniunitvec2
      end do
      ! print *, 'kaniunitvec + kaniunitvec2'
      ! do i=1, naucell
      !    print *, anisotropy_vecs(i,:)
      ! end do
      do i=1, num_anisotropies
         j = int(n_anisotropy(i,1))
         kaniunitvec3 = n_anisotropy(i,2:4)
         kanis3 = n_anisotropy(i,5)

         !make sure we have a unitarian vector
         kaniunitvec3 = kaniunitvec3 / norm2(kaniunitvec3)
         kaniunitvec3 = kaniunitvec3 * kanis3 * gamma**2/(mu_s**2)

         anisotropy_vecs(j+1,:) = anisotropy_vecs(j+1,:) + kaniunitvec3
      end do
      ! print *, 'with individual anisotropies'
      ! do i=1, naucell
      !    print *,  anisotropy_vecs(i,:)
      ! end do

      hmagunitvec = hmagunitvec / norm2(hmagunitvec)
      hmag(:) = hmagunitvec * (muB * hm0) * gamma

      Jx=Jnn; Jy=Jnn; Jz=Jnn
      big_a1 = a1 * latcons * n_basis_cells(1)
      big_a2 = a2 * latcons * n_basis_cells(2)
      big_a3 = a3 * latcons * n_basis_cells(3)

      big_a1_norm = norm2( big_a1 * (ncellpdim-1) )
      big_a2_norm = norm2( big_a2 * (ncellpdim-1) )
      big_a3_norm = norm2( big_a3 * (ncellpdim-1) )

      allocate( anglesphi(naucell), anglestheta(naucell), Si(naucell), Si_aux(naucell), basis(naucell,3), Rotmat(naucell,3,3), Rotpm(naucell,3,3) )
      
      !To calculate the k path where the dispersion will be calculated
      dim = sum(dims)
      amat(1,:) = a1 !/latcons*n_basis_cells(1)
      amat(2,:) = a2 !/latcons*n_basis_cells(2)
      amat(3,:) = a3 !/latcons*n_basis_cells(3)
      bmat = 0.d0 !This is important to initialize the part of the bmat matrix that wont be initialized on the next line
      bmat(1:dim,1:dim) = 2.d0*pi* transpose(inv(amat(1:dim,1:dim)))

      b1 = bmat(1,:)!*norm2(a1)
      b2 = bmat(2,:)!*norm2(a2)
      b3 = bmat(3,:)!*norm2(a3)
      
      call kpointsinitialization()
      !end: To calculate the k path where the dispersion will be calculated

      if(nptomega==0) then
         print *, "Setting 'nptomega' = 'effnkpt' "
         nptomega = effnkpt
      end if

      !Reading in the spin configuration
      inquire(file=basisname, exist=file_exists)
      if(file_exists) then
         open(unit=90, file=basisname, status="old")
      else
         print "(3a)", " File: '", trim(basisname), "' does NOT exist. Stopping."
         print *
         stop 
      end if

      !Skipping comment lines
      do 
         read(unit=90, fmt="(a)") read_in_data
         if( read_in_data(1:1) .NE. "#" ) then
            backspace(90)
            exit
         end if
      end do

      if(.not.spirit) then
         do counter=1, naucell
            read( unit=90, fmt=* ) basis(counter,:), anglesphi(counter), anglestheta(counter), Si(counter)
            ! print *, Si(counter)*sin(anglestheta(counter))*cos(anglesphi(counter)), Si(counter)*sin(anglestheta(counter))*sin(anglesphi(counter)), Si(counter)*cos(anglestheta(counter))
         end do
 
      else !if spirit = .true.
         print *, "'Spirit' mode is on."

         !Reading the spin configuration: spin orientation
         do i = 1, naucell 
            !In the 'spirit' mode, the spin orientation are read in cartesian coordinates.
            read(unit=90, fmt=*) Sx, Sy, Sz

            Sx = Sx * (mu_s/gamma)
            Sy = Sy * (mu_s/gamma) !+ hm0/(2*(kanis+20))               !This is the analytical solution for the spin config. for external field perpendicular to the mag. mom...
            Sz = Sz * (mu_s/gamma) !* cos( asin( hm0/(2*(kanis+20)) ) )
            Si(i) = norm2([Sx, Sy, Sz])
            anglestheta(i) = acos( Sz/Si(i))

            if( Si(i) < zero_toler ) then
               anglesphi(i) = 0.d0
            else
               anglesphi(i) = atan2( Sy, Sx )
            end if

            ! Special vectors with varying magnetic moment to account for the magnetic field
            Sx = Sx * mu_s_vec_aux(i) / mu_s
            Sy = Sy * mu_s_vec_aux(i) / mu_s 
            Sz = Sz * mu_s_vec_aux(i) / mu_s
            Si_aux(i) = norm2([Sx, Sy, Sz])
         end do

         if( spirit_input_given ) then
            counter = 0
            do k = 1, n_basis_cells(3)
               do j = 1, n_basis_cells(2)
                  do i = 1, n_basis_cells(1)
                     counter = counter + 1
                     basis(counter,:) = (i*a1 + j*a2 + k*a3)
                  end do 
               end do
            end do
            print *, "Position of each spin in the simulation box computed."
         else 
            !Reading the atom positions. In the 'spirit' mode, the atom positions are read from a separate file
            inquire(file=latticefile, exist=file_exists)
            if(file_exists) then
               open(unit=91, file=latticefile)
            else
               print "(3a)", " File: '", trim(latticefile), "' does NOT exist. Stopping."
               print *
               stop 
            end if

            !Check for commented line
            read(unit=91, fmt="(a)") read_in_data
            if( read_in_data(1:1) .NE. "#" ) rewind(91)

            do i = 1, naucell 
               read(unit=91, fmt=*) basisaux
               if( internalunit ) then 
                  basis(i,:) = a1*basisaux(1)+ a2*basisaux(2)+a3*basisaux(3)
               else
                  basis(i,:) = basisaux
               end if
            end do
         end if

         !Reading the interaction pairs
         if( n_nndists == 0 ) then ! If the interaction is NOT to be found by only measuring the n.n distances
            inquire(file=pairfile, exist=file_exists)
            if(file_exists) then
               open(unit=92, file=pairfile, status="old")
            else
               print "(3a)", " File: '", trim(pairfile), "' does NOT exist. Stopping."
               print *
               stop 
            end if

            !This factor of 2 is to ajust difference between this hamiltonian here and the spirit code
            if( J_D_scaling == 0.123456789d0 ) J_D_scaling = 2.d0
            print "(a,f7.3,a,f7.3,a)", " J and D are being rescalled by 'J_D_scaling/(mu_s**2)'=", J_D_scaling, "/", mu_s**2, "."
            print "(a)", " When 'Spirit mode' is on, if not given, 'J_D_scaling' is automatically set to 2 ajusting ..."
            print "(a)", "                    ... the differnce between the inner hamiltonian and the 'spirit code' one."

            if( abs(interactions_cutoff_radius - 1.d6) > zero_toler ) print *, "Interaction cuf-off: ", interactions_cutoff_radius

            if( spirit_input_given ) then
               print *, "Interactions pair positions will be calculated from d1*a1, d2*a2, d3*a3."
            else
               if( norm2(c1-c0) < zero_toler .or. norm2(c2-c0) < zero_toler .or. norm2(c3-c0) < zero_toler ) then
                  print "(a/,a)", " Basis index MUST start in 0 in the interactions file 'pairXX.txt'.", " Set of c1, c2, c3 NOT given." 
                  c1c2c3set_given = .false.
               else
                  print *, "Set of c1, c2, c3 GIVEN"
                  c1c2c3set_given = .true.
                  atoms_distinguisable = .false.
               end if
            end if

            counter = 0
            badindexing = .true.
            DJvect = 0.d0
            do 
               read(unit=92, fmt="(a)", iostat=reading_status) read_in_data
               if( reading_status < 0 ) exit
               if( reading_status > 0 ) stop "Sub initialization() error: On reading interaction pair file 'pairfile'. Stopping."
               
               read( read_in_data, *) word
               if( read_in_data(1:1) == "#" .or. trim(word) == "i" ) cycle

               read( read_in_data, fmt=* , iostat=reading_status) l1, l2, d1, d2, d3, Dx, Dy, Dz, Jnn

               if( reading_status .ne. 0 ) cycle
               ! if( reading_status .ne. 0 ) stop "Sub initialization() error: On reading interaction pair file 'pairfile'. Stopping."

               if( l1 == 0 ) badindexing = .false.
               if( l1 == naucell - 1 ) atoms_distinguisable = .true. !Test if the interaction for the last atom was given, then it assumes that the atom are considered distinguisible, and that the set of interaction for each atom is given

               counter = counter + 1
               if(counter>100000) stop "Sub initialization() error: There are more than 10000 entry on the 'pairfile'. 'mod_global' has to be modified."

               if( spirit_input_given ) then
                  positions(counter,:) = d1*a1 + d2*a2 + d3*a3
               else
                  if( c1c2c3set_given ) then
                     positions(counter,:) = d1*c1 + d2*c2 + d3*c3
                  else                        !          Position of second atom                     Position of the first
                     positions(counter,:) = ( basis(l2+1,:) + d1*big_a1 + d2*big_a2 + d3*big_a3 )    -  basis(l1+1,:)  
                     ijda_db_dc(counter,:) = [ l1+1, l2+1, d1, d2, d3 ] !+1 in the first two elements because in this program, the atom indices start counting at 1 and not zero (like in the Spirit-code)
                  end if
               end if

               pos_norm = norm2(positions(counter,:))

               do
                  if( ((pos_norm+1                 > big_a1_norm) .or. (pos_norm+1                 > big_a2_norm)) .and. &
                      ((interactions_cutoff_radius > big_a1_norm) .or. (interactions_cutoff_radius > big_a2_norm)) ) then
                      ! stop "Cluster of interaction larger than simulation box. Please increase 'ncellpdim' in the inputcard_XX.inp. Stopping."
                     ncellpdim = ncellpdim+1
                     big_a1_norm = norm2( big_a1 * (ncellpdim-1) )
                     big_a2_norm = norm2( big_a2 * (ncellpdim-1) )
                     big_a3_norm = norm2( big_a3 * (ncellpdim-1) )
                  else
                     exit
                  end if
               end do

               if( counter == 1 ) print *, "Interactions pair positions"
               if( counter <= 10 .and. counter<= naucell ) then
                  position = positions(counter,:)
                  print "(3(f12.8))", position
               end if  

               DJvect(counter,:) = [Dx, Dy, Dz, Jnn] * J_D_scaling * gamma**2/(mu_s**2)

               if( Dnn .ne. 0.d0 ) then
                  DJvect(counter,:) = [Dx*Dnn/norm2([Dx,Dy,Dz]), Dy*Dnn/norm2([Dx,Dy,Dz]), Dz*Dnn/norm2([Dx,Dy,Dz]), Jnn] * J_D_scaling * gamma**2/(mu_s**2)
               end if
            end do
            close(unit=92)

            print "(a, i0)", " ncellpdim = ", ncellpdim

            if( badindexing ) stop "Most probably your indexing is wrong in the interaction file. Basis atom index should start at zero. Stopping."
            if( atoms_distinguisable ) then
               print "(a)", " Atoms are DIStinguisible, and a set of interaction for each atom should have been given."
            else
               print "(a)", " Atoms are INDIStinguisible, and a set of interaction for the first atom only shoud have been given."
            end if            

            ninteraction_pairs = counter 
            print *, "Interaction parameters read from file. Number of interaction pairs:", ninteraction_pairs

         end if !in not spirit
      end if !If n_nndists == 0
      !end: Reading in the spin configuration

      !This determines the position in real space of each unit cell
      call crystallattice()
      !end: This determines the position in real space of each unit cell

      net_magnetization = 0.d0
      do i = 1, naucell
         net_magnetization = net_magnetization + [ cos(anglesphi(i))*sin(anglestheta(i))*Si(i), sin(anglesphi(i))*sin(anglestheta(i))*Si(i), cos(anglestheta(i))*Si(i) ]
      end do
      net_magnetization = net_magnetization/( naucell * (mu_s/gamma) )
      print "(a,3F22.16)", " Net magnetization / ( naucell * (mu_s/gamma) ) ", net_magnetization

      ! net_magnetization = net_magnetization/naucell
      ! if( mu_s .ne. 1.d0 .and. gamma .eq. 2.d0) then
      !    print "(a,3F22.16)", " Net magnetization / ( naucell * (mu_s/gamma) ) ", (net_magnetization)/norm2(net_magnetization)
      ! else
      !    print "(a,3F22.16)", " Net magnetization / ( naucell * (mu_s/gamma) ) ", net_magnetization * (gamma/mu_s)
      ! end if
      
      !Computing rotation matrixes
      do i = 1, naucell
         Rotmat(i,:,:) = rotationmatrix( anglesphi(i), anglestheta(i) )
      end do
      !Matrix of xyz --> +-z transformation (S+,S-,Sz) = HP (Sx,Sy,Sz)
      HP(1,:) = [ cone ,  ii  , czero ]
      HP(2,:) = [ cone , -ii  , czero ]
      HP(3,:) = [ czero, czero, cone  ]

      HPinv(1,:) = [  0.5d0*cone, 0.5d0*cone, czero ]
      HPinv(2,:) = [ -0.5d0*ii  , 0.5d0*ii  , czero ]
      HPinv(3,:) = [       czero,      czero, cone  ] 

      do i = 1, naucell
         aux = matmul( HP, Rotmat(i,:,:) )
         Rotpm(i,:,:) = matmul( aux, HPinv )

         ! !Rotating with respect to the polarization axis              sigma x                    !for test                              sigma y                                      sigma z
         ! aux(1,:) = [HP(1,1)*Rotmat(i,1,1) + HP(1,2)*Rotmat(i,2,1) + HP(1,3)*Rotmat(i,3,1), &    !for test
         !             HP(1,1)*Rotmat(i,1,2) + HP(1,2)*Rotmat(i,2,2) + HP(1,3)*Rotmat(i,3,2), &    !for test
         !             HP(1,1)*Rotmat(i,1,3) + HP(1,2)*Rotmat(i,2,3) + HP(1,3)*Rotmat(i,3,3)]      !for test
         ! aux(2,:) = [HP(2,1)*Rotmat(i,1,1) + HP(2,2)*Rotmat(i,2,1) + HP(2,3)*Rotmat(i,3,1), &    !for test
         !             HP(2,1)*Rotmat(i,1,2) + HP(2,2)*Rotmat(i,2,2) + HP(2,3)*Rotmat(i,3,2), &    !for test
         !             HP(2,1)*Rotmat(i,1,3) + HP(2,2)*Rotmat(i,2,3) + HP(2,3)*Rotmat(i,3,3)]      !for test
         ! aux(3,:) = [HP(3,1)*Rotmat(i,1,1) + HP(3,2)*Rotmat(i,2,1) + HP(3,3)*Rotmat(i,3,1), &    !for test
         !             HP(3,1)*Rotmat(i,1,2) + HP(3,2)*Rotmat(i,2,2) + HP(3,3)*Rotmat(i,3,2), &    !for test
         !             HP(3,1)*Rotmat(i,1,3) + HP(3,2)*Rotmat(i,2,3) + HP(3,3)*Rotmat(i,3,3)]      !for test

         ! !Rotating with respect to the polarization axis              sigma x                    !for test                             sigma y                                      sigma z
         ! Rotpm(i,1,:) = [aux(1,1)*HPinv(1,1) + aux(1,2)*HPinv(2,1) + aux(1,3)*HPinv(3,1), &      !for test
         !             aux(1,1)*HPinv(1,2) + aux(1,2)*HPinv(2,2) + aux(1,3)*HPinv(3,2), &          !for test
         !             aux(1,1)*HPinv(1,3) + aux(1,2)*HPinv(2,3) + aux(1,3)*HPinv(3,3)]            !for test
         ! Rotpm(i,2,:) = [aux(2,1)*HPinv(1,1) + aux(2,2)*HPinv(2,1) + aux(2,3)*HPinv(3,1), &      !for test
         !             aux(2,1)*HPinv(1,2) + aux(2,2)*HPinv(2,2) + aux(2,3)*HPinv(3,2), &          !for test
         !             aux(2,1)*HPinv(1,3) + aux(2,2)*HPinv(2,3) + aux(2,3)*HPinv(3,3)]            !for test
         ! Rotpm(i,3,:) = [aux(3,1)*HPinv(1,1) + aux(3,2)*HPinv(2,1) + aux(3,3)*HPinv(3,1), &      !for test
         !             aux(3,1)*HPinv(1,2) + aux(3,2)*HPinv(2,2) + aux(3,3)*HPinv(3,2), &          !for test
         !             aux(3,1)*HPinv(1,3) + aux(3,2)*HPinv(2,3) + aux(3,3)*HPinv(3,3)]            !for test


         ! Rotpm(i,:,:) =  Rotmat(i,:,:) !for test
         ! aux = matmul(aux,HPinv)
         ! print "(3('(',f4.1,',',f4.1,') '))", aux(1,:)
         ! print "(3('(',f4.1,',',f4.1,') '))", aux(2,:)
         ! print "(3('(',f4.1,',',f4.1,') '))", aux(3,:); pause
      end do

      !Pauli Matrices
      !Sigma x
      paulimatrix_cartesian(1,1,:)=[0.d0, 1.d0]
      paulimatrix_cartesian(1,2,:)=[1.d0, 0.d0]
      
      !Sigma y
      paulimatrix_cartesian(2,1,:)=[czero,   -ii]
      paulimatrix_cartesian(2,2,:)=[   ii, czero]

      !Sigma z
      paulimatrix_cartesian(3,1,:)=[1.d0, 0.d0]
      paulimatrix_cartesian(3,2,:)=[0.d0,-1.d0]

      polarizations(1,:) = polarization1 
      polarizations(2,:) = polarization2 
      polarizations(3,:) = polarization3   

      !Coverting Pauli matrices 
      do polariz = 1, 3
         Rotmat_aux(:,:) = rotationmatrix( polarizations(polariz,1), polarizations(polariz,2) )

         !Rotating with respect to the polarization axis              sigma x                                          sigma y                                      sigma z
         paulimatrix_cartesian_rotated(1,:,:) = Rotmat_aux(1,1)*paulimatrix_cartesian(1,:,:) + Rotmat_aux(1,2)*paulimatrix_cartesian(2,:,:) + Rotmat_aux(1,3)*paulimatrix_cartesian(3,:,:) 
         paulimatrix_cartesian_rotated(2,:,:) = Rotmat_aux(2,1)*paulimatrix_cartesian(1,:,:) + Rotmat_aux(2,2)*paulimatrix_cartesian(2,:,:) + Rotmat_aux(2,3)*paulimatrix_cartesian(3,:,:) 
         paulimatrix_cartesian_rotated(3,:,:) = Rotmat_aux(3,1)*paulimatrix_cartesian(1,:,:) + Rotmat_aux(3,2)*paulimatrix_cartesian(2,:,:) + Rotmat_aux(3,3)*paulimatrix_cartesian(3,:,:) 

         !Changing from the xyz to the +-z notation
         paulimatrix(polariz,1,:,:) = 0.5d0*(paulimatrix_cartesian_rotated(1,:,:)+ii*paulimatrix_cartesian_rotated(2,:,:))
         paulimatrix(polariz,2,:,:) = 0.5d0*(paulimatrix_cartesian_rotated(1,:,:)-ii*paulimatrix_cartesian_rotated(2,:,:))
         paulimatrix(polariz,3,:,:) = paulimatrix_cartesian_rotated(3,:,:)
         ! paulimatrix(polariz,1,:,:)=paulimatrix_cartesian_rotated(1,:,:) !for test
         ! paulimatrix(polariz,2,:,:)=paulimatrix_cartesian_rotated(2,:,:) !for test
         ! paulimatrix(polariz,3,:,:)=paulimatrix_cartesian_rotated(3,:,:) !for test
      end do
      !end: computing rotation matrixes
   end subroutine initialization

!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
   subroutine crystallattice()
      !Determine the position of each unit cell around a reference point
      !These positions will be used to Fourier transform functions to the reciprocal space
      implicit none
      integer :: i, j, k
      real(kind=pc) :: r0jaux( (2*ncellpdim+1)**3 ,3), r0(3), rj(3)

      ncell = 0
      r0jaux = 0.d0
      r0 = 0*big_a1 + 0*big_a2 + 0*big_a3
      !The lattice can be 1, 2 or 3 dimensional depending on the inputted array "dims"
      do i = -ncellpdim*dims(1), ncellpdim*dims(1)
      do j = -ncellpdim*dims(2), ncellpdim*dims(2)
      do k = -ncellpdim*dims(3), ncellpdim*dims(3)
         ncell = ncell + 1
         if(i==0 .and. j==0 .and. k==0) then
            origcell = ncell
            r0jaux(ncell,:) = r0
         else
            rj = i*big_a1 + j*big_a2 + k*big_a3
            r0jaux(ncell,:) = rj-r0
         end if
         ! print *, "i, j, k, ncell", i, j, k, ncell
      end do
      end do
      end do
      allocate( r0j(ncell,3) )
      r0j = r0jaux(1:ncell,:)
      print *, "Origcell: ", origcell
   end subroutine crystallattice

!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
   subroutine kpointsinitialization()
      !Compute and store k points on the chosen path over which the dispersion will be calculated
      implicit none
      integer :: i, j, counter
      real(kind=pc) :: k(3), prevk(3), pathelement(3), totalpathlength

      allocate( symmpts(npath+1,3), nkpt_npath(npath) )

      do i = 1, npath+1
         symmpts(i,:) = syptsMat(i,1)*b1 + syptsMat(i,2)*b2 + syptsMat(i,3)*b3
      end do

      !Create kpoint mesh (area)
      if( kpoint_mesh ) then

         print "(a,3i4)", " Creating kpoint mesh aroung the first given high symmetry point "

         effnkpt = npt**2
         allocate( kpoints(effnkpt,3), dkpoints(effnkpt) )

         counter = 0
         prevk = symmpts(1,:)
         do i = 1, npt
            do j = 1, npt
               counter = counter + 1

               k = 1.d0*b1*norm2(a1) * ((i-npt/2.d0)/npt) + 1.0d0*b2*norm2(a2) * ((j-npt/2.d0)/npt) + symmpts(1,:)
            
               kpoints(counter, :) = k 
               dkpoints(counter) = norm2(k-prevk)

               prevk = k
            end do
         end do

      else ! To create a kpoint path instead of an area

         totalpathlength = 0.d0
         do j = 1, npath
            k = symmpts(j+1,:)-symmpts(j,:)
            totalpathlength = totalpathlength + sqrt(dot_product(k,k)) 
         end do
    
         effnkpt = 0
         do j = 1, npath
            k = symmpts(j+1,:)-symmpts(j,:)

            ! The mininum amount of point is now 1, also to accept path of length zero, i.e., to calculate on a single point.
            nkpt_npath(j) = max( int( dble(npt) * sqrt(dot_product(k,k)) / totalpathlength ), 1 )
            
            !This if is to accept path of length zero, i.e., to calculate on a single point.
            ! if( abs(totalpathlength) < zero_toler ) then
            !    nkpt_npath(j) = 1
            ! else 
            !    nkpt_npath(j) = max( int( dble(npt) * sqrt(dot_product(k,k)) / totalpathlength ), 1 )
            ! end if

            effnkpt = effnkpt + nkpt_npath(j)
            if(nkpt_npath(j)<1) then
               print "(a,i0,a,i0,a)", " Number of k points in the path ", j, " is too small: ", nkpt_npath(j), ". Increase 'npt'. Stopping."
               stop
            end if
         end do
         effnkpt = effnkpt+1 !One more point to include the last kpoint

         allocate( kpoints(effnkpt,3), dkpoints(effnkpt) )

         counter = 0
         prevk = symmpts(1,:)
         !Each path is a line segment given by two points informed in the inputcard
         !Loop of the paths
         do j = 1, npath
            pathelement = ( symmpts(j+1,:) - symmpts(j,:) ) / dble(nkpt_npath(j))
            do i = 1, nkpt_npath(j)
               counter = counter + 1
               k = symmpts(j,:) + (i-1)*pathelement
               kpoints(counter,:) = k
               dkpoints(counter) = norm2(k-prevk)
               prevk = k
            end do
         end do
         !Including the last kpoint
         counter = counter + 1
         k = symmpts(npath+1,:)
         kpoints(counter,:) = k
         dkpoints(counter) = norm2(k-prevk)

      end if ! if( kpoint_mesh )

   end subroutine kpointsinitialization
 
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
   function rotationmatrix(phi,theta)
      !Return the rotation matrix given an angle phi and a theta
      !rotationmatrix = R_z(phi) . R_y(theta)
      implicit none
      real(kind=pc), intent(in) :: phi, theta
      real(kind=pc) :: rotationmatrix(3,3)

      rotationmatrix(1,:) = [ cos(phi)*cos(theta), -sin(phi), cos(phi)*sin(theta) ]
      rotationmatrix(2,:) = [ sin(phi)*cos(theta),  cos(phi), sin(phi)*sin(theta) ]
      rotationmatrix(3,:) = [     -sin(theta)    ,    0.d0  ,     cos(theta)      ]
      return
   end function rotationmatrix

!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
   function cross(a, b)
      implicit none
      real(kind=pc) :: cross(3)
      real(kind=pc), intent(in) :: a(3), b(3)

      cross(1) = a(2) * b(3) - a(3) * b(2)
      cross(2) = a(3) * b(1) - a(1) * b(3)
      cross(3) = a(1) * b(2) - a(2) * b(1)
   end function cross

!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
   subroutine commutationmatrix(Dmatrix)
      implicit none
      complex(kind=pc), intent(out) :: Dmatrix(effnkpt,twonaucell,twonaucell) !Dmatrix is the one to be diagonalized to provide the spinwave dispersion.
      integer :: m, n, l1, l2, i, j, nnn, sumnnn, ikptn, n_inter_per_atom(naucell)
      real(kind=pc) :: r0jaux(3), Rotmati(3,3), Rotmatj(3,3), norm_r0j, k(3)
      complex(kind=pc) :: HPmatrixi(3,3), HPmatrixj(3,3), HPmatrixj_aux(3,3), HPbase(3,3), Dk(naucell,naucell,2,2), D0j(naucell,ncell,naucell,2,2), A0j_tilde(naucell,ncell,naucell,2,2), J0j(naucell,ncell,naucell,3,3), hmagvec(naucell,3), Jtilde0zz(naucell), Jtemp(3,3), Jtemp2(3,3), sumJxy
      logical :: indices_match

      !A base format of the HP transformation matrix
      HPbase(1,:) = [ cone , cone , czero ]
      HPbase(2,:) = [-ii, ii, czero ]
      HPbase(3,:) = [ czero, czero, czero ]  

      J0j = 0.d0

      !Calculation of the interaction matrix, which will be determined by the n.n distance and the inputted parameters
      if( (.not. spirit) .or. spirit_old_method .or. (n_nndists .ne. 0) ) then
         print "(a)", " Dmatrix: scanning the real space method. (old)"
         sumnnn = 0
         do l1 = 1, naucell
            Rotmati = Rotmat(l1,:,:)

            HPmatrixi = HPbase
            HPmatrixi(3,3) = sqrt(ctwo/Si(l1))
            HPmatrixi = sqrt(Si(l1)/2.d0)*HPmatrixi   

            nnn = 0
            do l2 = 1, naucell
               Rotmatj = Rotmat(l2,:,:)
               
               HPmatrixj = HPbase
               HPmatrixj(3,3) = sqrt(ctwo/Si(l2))
               HPmatrixj = sqrt(Si(l2)/2.d0)*HPmatrixj
               
               HPmatrixj_aux = HPbase
               HPmatrixj_aux(3,3) = sqrt(ctwo/Si_aux(l2))
               HPmatrixj_aux = sqrt(Si_aux(l2)/2.d0)*HPmatrixj_aux

               do j = 1, ncell 
                  !        Position of the second atom         Position of the first atom 
                  r0jaux =( r0j(j,:) + basis(l2,:) )   -  ( r0j(origcell,:) + basis(l1,:) )
                  norm_r0j = norm2(r0jaux)

                  Jtemp = 0.d0
                  if(spirit_old_method .and. n_nndists==0) then
                     if(norm_r0j <= interactions_cutoff_radius) then
                        do i = 1, ninteraction_pairs

                           if( atoms_distinguisable ) then
                              !Checking for bound
                              norm_r0j = norm2( r0jaux-positions(i,:) )
                              ! norm_r0j = norm2( r0j(j,:) - ( big_a1*ijda_db_dc(i,3) + big_a2*ijda_db_dc(i,4) + big_a3*ijda_db_dc(i,5) ) )

                              !Checking for the right indicex
                              indices_match = .false.
                              if( l1 == ijda_db_dc(i,1) .and. l2 == ijda_db_dc(i,2)  ) indices_match = .true.

                           else !if the atoms all the same, that is, if the unit cell are made of atoms of a simulation box
                              !Checking for bound
                              norm_r0j = norm2( r0jaux-positions(i,:) )
                              indices_match = .true. !if the atoms are indistiguisible, we do NOT check the indeces
                           end if

                           if( norm_r0j < zero_toler .and. indices_match ) then
                              nnn=nnn+1

                              Dx = DJvect(i,1)
                              Dy = DJvect(i,2)
                              Dz = DJvect(i,3)
                              Jx = DJvect(i,4)
                              Jy = DJvect(i,4)
                              Jz = DJvect(i,4)

                              Jtemp(1,:) = [ Jx, Dz,-Dy ]
                              Jtemp(2,:) = [-Dz, Jy, Dx ]
                              Jtemp(3,:) = [ Dy,-Dx, Jz ]
                              exit
                           end if
                        end do
                     end if
                  elseif( n_nndists .ne. 0 ) then
                     if( toprint_Rij .and. norm_r0j < 10.d0 ) print *, "|Rij|=", norm_r0j

                     do i=1, n_nndists
                        nndist = all_nndist(i)

                        !Test to see if the current pair of atoms are at n.n distance
                        if( abs(norm_r0j - nndist) < zero_toler .and. norm_r0j > zero_toler ) then
                           Jnn =  all_Jnn(i) * J_D_scaling * gamma**2/(mu_s**2)

                           Jx = Jnn; Jy = Jnn; Jz = Jnn
                           Dx = 0.d0; Dy = 0.d0; Dz = 0.d0

                           Jtemp(1,:) = [ Jx, Dz,-Dy ]
                           Jtemp(2,:) = [-Dz, Jy, Dx ]
                           Jtemp(3,:) = [ Dy,-Dx, Jz ]

                           nnn=nnn+1
                        end if
                     end do
                          !Test to see if the current pair of atoms are at n.n distance
                  elseif( abs(norm_r0j - nndist) < zero_toler .and. norm_r0j > zero_toler ) then
                     !The DMI vector is always rotated 90 degrees in relation to the vector connecting the pair of atoms
                     Dx =-Dnn*r0jaux(2)/norm_r0j
                     Dy = Dnn*r0jaux(1)/norm_r0j
                     Dz = 0.d0
                     ! Dz = Dnn*r0jaux(1)/norm_r0j
                     ! Dy = 0.d0

                     Jtemp(1,:) = [ Jx, Dz,-Dy ]
                     Jtemp(2,:) = [-Dz, Jy, Dx ]
                     Jtemp(3,:) = [ Dy,-Dx, Jz ]

                     nnn=nnn+1
                  end if !if spirit mode

                  if(j==origcell .and. l1==l2) then
                     Jtemp(1,1) = Jtemp(1,1) + 2.0d0*(anisotropy_vecs(l1,1))
                     Jtemp(2,2) = Jtemp(2,2) + 2.0d0*(anisotropy_vecs(l1,2))
                     Jtemp(3,3) = Jtemp(3,3) + 2.0d0*(anisotropy_vecs(l1,3))
                     ! Jtemp(1,1) = Jtemp(1,1) + 2.0d0*(kaniunitvec(1) + kaniunitvec2(1))
                     ! Jtemp(2,2) = Jtemp(2,2) + 2.0d0*(kaniunitvec(2) + kaniunitvec2(2))
                     ! Jtemp(3,3) = Jtemp(3,3) + 2.0d0*(kaniunitvec(3) + kaniunitvec2(3))
                  end if
                  ! if(l1==2 .and. l2==2 .and. j==116) print "(a,8f9.4)", "J0j_aux", Jtemp(1,1), Jtemp(1,2), Jtemp(2,1), Jtemp(2,2) !For test, to be deleted

                  Jtemp2 = matmul(  Jtemp, Rotmatj )
                  Jtemp = matmul( transpose(Rotmati) , Jtemp2 )

                  if(.not. Tneel) Jtemp2 = matmul( Jtemp , HPmatrixj )
                  if(.not. Tneel) Jtemp = matmul( conjg(transpose(HPmatrixi)) , Jtemp2 )

                  if( j==origcell .and. l1==l2 ) then
                     hmagvec(l2,:) = matmul( hmag(:), Rotmatj )
                     if(.not. Tneel) hmagvec(l2,:) = matmul( hmagvec(l2,:), HPmatrixj_aux )
                     ! if(.not. Tneel) hmagvec(l2,:) = matmul( hmagvec(l2,:), HPmatrixj )
                  end if

                  J0j(l1,j,l2,:,:) = Jtemp 
               end do !over unit cells
            end do !over atoms in the unit cell l2
            sumnnn = sumnnn + nnn
         end do !over atoms in the unit cell l1
      end if

      !SELF NOTE: I still have to implement the method below even when the atoms are indistiguisible
      if( spirit .and. (.not. spirit_old_method) .and. n_nndists == 0 ) then
         print "(a)", " Dmatrix: using the direct attribution method for the Spirit mode. (new)"
         n_inter_per_atom = 0
         do i = 1, ninteraction_pairs
            Dx = DJvect(i,1)
            Dy = DJvect(i,2)
            Dz = DJvect(i,3)
            Jx = DJvect(i,4)
            Jy = DJvect(i,4)
            Jz = DJvect(i,4)

            Jtemp = 0.d0

            Jtemp(1,:) = [ Jx, Dz,-Dy ]
            Jtemp(2,:) = [-Dz, Jy, Dx ]
            Jtemp(3,:) = [ Dy,-Dx, Jz ]

            l1 = ijda_db_dc(i,1)
            l2 = ijda_db_dc(i,2)

            if( 1 > l1 > naucell .or. 1 > l2 > naucell  ) stop("Atom indeces failed to test: 1 <= [l1, l2] <= 'naucell'. Stopping.")
            
            do n = 1, ncell
               r0jaux = ijda_db_dc(i,3)*a1 + ijda_db_dc(i,4)*a2 + ijda_db_dc(i,5)*a3 
               if( norm2(r0jaux - r0j(n,:)) < zero_toler ) then
                  j = n
                  exit
               end if
            end do
            ! if(l1==2 .and. l2==2 .and. j==116) print "(a,8f9.4)" , "J0j    ", Jtemp(1,1), Jtemp(1,2), Jtemp(2,1), Jtemp(2,2)  !For test, to be deleted

            if( 1 > j > ncell  ) stop("Unit cell index failed the test: 1 <= j <= 'ncell'. Stopping.")

            n_inter_per_atom(l1) = n_inter_per_atom(l1) + 1

            Rotmati = Rotmat(l1,:,:)

            HPmatrixi = HPbase
            HPmatrixi(3,3) = sqrt(ctwo/Si(l1))
            HPmatrixi = sqrt(Si(l1)/2.d0)*HPmatrixi   

            Rotmatj = Rotmat(l2,:,:)
            
            HPmatrixj = HPbase
            HPmatrixj(3,3) = sqrt(ctwo/Si(l2))
            HPmatrixj = sqrt(Si(l2)/2.d0)*HPmatrixj

            Jtemp2 = matmul(  Jtemp, Rotmatj )
            Jtemp = matmul( transpose(Rotmati) , Jtemp2 )

            if( .not. Tneel ) Jtemp2 = matmul( Jtemp , HPmatrixj )
            if( .not. Tneel ) Jtemp = matmul( conjg(transpose(HPmatrixi)) , Jtemp2 )

            J0j(l1,j,l2,:,:) = Jtemp
         end do

         j = origcell
         do l1=1, naucell
            l2 = l1
            
            Jtemp = 0.d0

            Rotmatj = Rotmat(l2,:,:)
            
            HPmatrixj = HPbase
            HPmatrixj(3,3) = sqrt(ctwo/Si(l2))
            HPmatrixj = sqrt(Si(l2)/2.d0)*HPmatrixj

            HPmatrixj_aux = HPbase
            HPmatrixj_aux(3,3) = sqrt(ctwo/Si_aux(l2))
            HPmatrixj_aux = sqrt(Si_aux(l2)/2.d0)*HPmatrixj_aux
            
            Jtemp(1,1) = 2.0d0*(anisotropy_vecs(l1,1))
            Jtemp(2,2) = 2.0d0*(anisotropy_vecs(l1,2))
            Jtemp(3,3) = 2.0d0*(anisotropy_vecs(l1,3))
            ! Jtemp(1,1) = 2.0d0*(kaniunitvec(1) + kaniunitvec2(1))
            ! Jtemp(2,2) = 2.0d0*(kaniunitvec(2) + kaniunitvec2(2))
            ! Jtemp(3,3) = 2.0d0*(kaniunitvec(3) + kaniunitvec2(3))
            
            hmagvec(l2,:) = matmul( hmag(:), Rotmatj )
            if( .not. Tneel ) hmagvec(l2,:) = matmul( hmagvec(l2,:), HPmatrixj_aux )
            ! if( .not. Tneel ) hmagvec(l2,:) = matmul( hmagvec(l2,:), HPmatrixj )

            J0j(l1,j,l2,:,:) = Jtemp
         end do

         sumnnn = sum( n_inter_per_atom )
      end if

      !This was used to compare the new Dmatrix construction from the pair interaction file agaist the old method. Efentually to be deleted.
      ! do l1=1, naucell
      ! do l2=1, naucell
      ! do j=1, ncell
      !    if( abs( J0j(l1,j,l2,1,1) - J0j_aux(l1,j,l2,1,1) ) > zero_toler ) print "(3(a,i3), 4f14.8)", "MERDA 1 for l1=", l1, " l2=", l2, " j=", j, J0j(l1,j,l2,1,1), J0j_aux(l1,j,l2,1,1)
      !    if( abs( J0j(l1,j,l2,1,2) - J0j_aux(l1,j,l2,1,2) ) > zero_toler ) print "(3(a,i3), 4f14.8)", "MERDA 2 for l1=", l1, " l2=", l2, " j=", j, J0j(l1,j,l2,1,2), J0j_aux(l1,j,l2,1,2)
      !    if( abs( J0j(l1,j,l2,2,1) - J0j_aux(l1,j,l2,2,1) ) > zero_toler ) print "(3(a,i3), 4f14.8)", "MERDA 3 for l1=", l1, " l2=", l2, " j=", j, J0j(l1,j,l2,2,1), J0j_aux(l1,j,l2,2,1)
      !    if( abs( J0j(l1,j,l2,2,2) - J0j_aux(l1,j,l2,2,2) ) > zero_toler ) print "(3(a,i3), 4f14.8)", "MERDA 4 for l1=", l1, " l2=", l2, " j=", j, J0j(l1,j,l2,2,2), J0j_aux(l1,j,l2,2,2)
      ! end do   
      ! end do   
      ! end do   

      !Informs the number of n.n sites
      print "(a,f18.10)", " Average number of atoms in the cluster (n. near neighbours):             ", dble(sumnnn)/dble(naucell)

      !Tests the vanishing of the linear term of the HP transformation
      do l1 = 1, naucell
         sumJxy = hmagvec(l1,2)
         do l2 = 1, naucell
            do i = 1, ncell
               sumJxy = sumJxy + J0j(l1,i,l2,1,3)*Si(l2)
               ! sumJxy = sumJxy + J0j(l1,i,l2,1,3)*Si(l2)
            end do
         end do
         !In the ground state the following elements should be zero
         ! if( abs( sumJxy ) > zero_toler ) print "(a,i4,2f16.8)", " Warning: |Jxz|+|Jyz| finite for l1, sumJxy=", l1, sumJxy
         if( abs( sumJxy ) > zero_toler ) print "(a,i4,2f16.8)", achar(27)//"[31m Warning: |Jxz|+|Jyz| finite for l1, sumJxy="//achar(27)//"[0m", l1, sumJxy
      end do

      !Collects the circular parts (+,-) of the the matrix J0j
      A0j_tilde=0.d0
      do l1 = 1, naucell; do l2 = 1, naucell
         do i = 1, ncell
            A0j_tilde(l1,i,l2,:,:) = J0j(l1,i,l2,1:2,1:2)
         end do
      end do; end do

      !Calculate Jzz(0)
      Jtilde0zz = 0.d0
      do l1 = 1, naucell; do l2 = 1, naucell
         do i = 1, ncell
            Jtilde0zz(l1) = Jtilde0zz(l1) + Si(l2) * J0j(l1,i,l2,3,3)
         end do
      end do; end do
      ! print *,"Constant:", Jtilde0zz; stop


      !Dij construction
      if( Tneel ) then
         D0j = 0.d0 
         D0j(:,:,:,1,1) = A0j_tilde(:,:,:,1,1) 
      else
         j = origcell
         do l1=1, naucell
            A0j_tilde(l1,j,l1,1,1) = A0j_tilde(l1,j,l1,1,1) - ( hmagvec(l1,3) + Jtilde0zz(l1) )
            A0j_tilde(l1,j,l1,2,2) = A0j_tilde(l1,j,l1,2,2) - ( hmagvec(l1,3) + Jtilde0zz(l1) )
         end do

         D0j(:,:,:,1,:) = -A0j_tilde(:,:,:,1,:) 
         D0j(:,:,:,2,:) =  A0j_tilde(:,:,:,2,:)
      end if

      do ikptn = 1, effnkpt
         k = kpoints(ikptn,:)
         !Fourier Transformation
         Dk(:,:,:,:)=0.d0
         do l1 = 1, naucell; do l2 = 1, naucell
            do j = 1, ncell
               Dk(l1,l2,:,:) = Dk(l1,l2,:,:) + exp(  ii*dot_product(k,r0j(j,:))  ) * D0j(l1,j,l2,:,:)
            end do
         end do; end do

         if( printD .and. ikptn == 1) then
            print *, "Printing the dynamical matrix for the 1st k point. Only work for two atoms in the unit cell."
            print "(2f18.10,a,2f18.10)", real(Dk(1,1,1,1)), real(Dk(1,1,1,2))," |", real(Dk(1,2,1,1)), real(Dk(1,2,1,2))
            print "(2f18.10,a,2f18.10)", real(Dk(1,1,2,1)), real(Dk(1,1,2,2))," |", real(Dk(1,2,2,1)), real(Dk(1,2,2,2))
            print "(2f18.10,a,2f18.10)", real(Dk(2,1,1,1)), real(Dk(2,1,1,2))," |", real(Dk(2,2,1,1)), real(Dk(2,2,1,2))
            print "(2f18.10,a,2f18.10)", real(Dk(2,1,2,1)), real(Dk(2,1,2,2))," |", real(Dk(2,2,2,1)), real(Dk(2,2,2,2))
         end if

         !Expanding Dk matrix (4 dim.) to a 2 dim. matrix.
         do i = 1, 2; do j = 1, 2
            do l1 = 1, naucell; do l2 = 1, naucell
               m = (i-1)*naucell+l1
               n = (j-1)*naucell+l2
               Dmatrix(ikptn,m,n) = Dk(l1,l2,i,j)
            end do; end do
         end do; end do
         ! do l1 = 1, naucell; do l2 = 1, naucell
         !    do i = 1, 2; do j = 1, 2
         !       m = 2*(l1-1)+1
         !       n = 2*(l2-1)+1
         !       Dmatrix(ikptn,m:m+1,n:n+1) = Dk(l1,l2,:,:)
         !    end do; end do
         ! end do; end do

         !Transforming D -> D + epsilon . g , where g is a matrix of entries 1 in the first half of the diagonal elements and -1 in the other half. This avoids to have null energies and therefore possible null eigen-vec.
         do j=1, naucell
            Dmatrix(ikptn,j,j)=Dmatrix(ikptn,j,j) + zero_toler*0.d0
            Dmatrix(ikptn,j+naucell,j+naucell)=Dmatrix(ikptn,j+naucell,j+naucell) - zero_toler*0.d0
         end do

! Dmatrix(ikptn,2,:)=Dmatrix(ikptn,2,:)+Dmatrix(ikptn,6,:)
! do m=1, twonaucell
!    print "(8f6.1)", (REAL(Dmatrix(ikptn,m,n)),n=1,twonaucell)
! end do; stop   
      end do !Loop of k

   end subroutine commutationmatrix

!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
   subroutine printinginput()
      !Subroutine to display the input parameters.
      implicit none
      character(len=50) :: formt, suffix
      integer :: i, j, naucellaux, m

      naucellaux=naucell
      if(naucell>5) naucellaux = 5
      write(formt,fmt="(a,i0,a)") "(f12.8,'| ',", naucellaux+1, "f16.8)"
      print *
      print "(a)", "High symmetry point path location:"
      print *, " Path      | Correspondent frequency value"
      do i = 1, npath
         write( unit=*, fmt=formt) real(highsymm(i,1)),  ( real(highsymm(i,naucell+1-j)), j=0, naucellaux)  
      end do

      open(unit=456,file="outputfiles/highsympoints.gnu")
      write(formt,fmt="(a,a,a,i0,a)") "('set xtics (", '""' ,"  0.0',", npath, "(a,f12.8),')')"
      write( unit=456, fmt=formt) (', ""', real(highsymm(i,1)), i=1, npath)  

      print *
      print "(a)", "Input parameters:"
      print *
      print "(a,3i2,a)", " dims =",dims, ","
      print "(4(a,i0),a)", " npt = ",npt, ", naucell = ", naucell, ", ncellpdim = ", ncellpdim, ", npath = ", npath, ","
      print "(a,f4.1,a)", " latcons =", latcons, "d0,"
      print "(3(a,f5.2),a)", " a1 = ", a1(1), "d0   ", a1(2), "d0   ", a1(3), "d0,"
      print "(3(a,f5.2),a)", " a2 = ", a2(1), "d0   ", a2(2), "d0   ", a2(3), "d0,"
      print "(3(a,f5.2),a)", " a3 = ", a3(1), "d0   ", a3(2), "d0   ", a3(3), "d0,"
      print "(3(a,f5.2),a)", " b1 = ", b1(1), "d0   ", b1(2), "d0   ", b1(3), "d0,"
      print "(3(a,f5.2),a)", " b2 = ", b2(1), "d0   ", b2(2), "d0   ", b2(3), "d0,"
      print "(3(a,f5.2),a)", " b3 = ", b3(1), "d0   ", b3(2), "d0   ", b3(3), "d0,"
      print "(2(a,f5.2),a)", " polarization1 = ", polarization1(1), "d0   ", polarization1(2), "d0,"
      print "(2(a,f5.2),a)", " polarization2 = ", polarization2(1), "d0   ", polarization2(2), "d0,"
      print "(2(a,f5.2),a)", " polarization3 = ", polarization3(1), "d0   ", polarization3(2), "d0,"
      print *
      print "(4(a,f4.1),a)", " Jnn =", Jnn, "d0, Dnn =", Dnn, "d0, kanis =", kanis, "d0, hm0 =", hm0, "d0,"
      print "(3(a,f5.2),a)", " hmagunitvec = ",hmagunitvec(1), "d0 ",hmagunitvec(2), "d0 ",hmagunitvec(3), "d0,"
      print *
      print "(a)", " syptsMat = "
      do i = 1, npath+1
         print "(3f7.2)", (syptsMat(i,j),j=1,3)
      end do
      print *
      print *, "Basis:"
      print *, "Positions           |      phi  theta"
      do i = 1, naucellaux
         print "(3f7.2,a,2f6.2)", basis(i,:),"     ", anglesphi(i), anglestheta(i)
         ! print "(3f15.8,a,2f15.8)", basis(i,:),"     ", anglesphi(i)/pi * 180, anglestheta(i)/pi * 180
      end do
      if(naucellaux<naucell) print *, "..."

      i=index(inputcardname,"_")+1
      j=index(inputcardname,".")-1
      suffix = inputcardname(i:j) 

      !Storage file openning
      formt="latticExt_"//trim(suffix)//".dat"
      open( unit=666, file=formt )

      !To print the crystal
      do i = -dims(1)*spin_dyn,dims(1)*spin_dyn; do m = -dims(2)*spin_dyn,dims(2)*spin_dyn
      do j = 1, naucell
         
         !These 4 next lines were used to write a spin spiral lattice in 2d with a 16a pitch
         ! k = -[0.d0,2.d0*pi/15.d0,0.d0]
         ! r0j(i,:) = a1*i + a2*m + basis(j,:)
         ! anglesphi(j) = pi/2.d0
         ! anglestheta(j) = dot_product(k,r0j(i,:))
         
         write(unit=666,fmt='(9f14.8)') big_a1*i + big_a2*m + basis(j,:), anglesphi(j), anglestheta(j), Si(j), big_a1*i + big_a2*m
         
         ! k = r0j(i,:)
         ! write(unit=666,fmt='(6f14.8)') k + basis(j,:), anglesphi(j), anglestheta(j), Si(j)
      end do
      end do; end do

      ! end do
   end subroutine printinginput

!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
   function omega(k)
      implicit none
      complex(kind=pc), dimension(2) :: omega
      real(kind=pc), intent(in) :: k(3)
      real, parameter :: S=1.0d0
      real :: aa1(3), aa2(3)
      integer :: info=0
      real(kind=pc) :: alpha, q(3), c1, c2, s1, s2, cphi, sphi, k_a1, k_a2, phi!, JD, qq
      complex(kind=pc) :: Dpp, Dpm, prefact1pp, prefact2pp, prefact1pm, prefact2pm, J00, mat(2,2), evalues(2), leftevector(2,2), rightevector(2,2)!, Dmp
      complex(kind=pc) :: work(4) !work space for the diagonalization routine.
      double precision :: rwork(4) !work space for the diagonalization routine.     

      alpha=atan2(Dnn,Jnn)

      aa1=latcons*[1.0d0,0.0d0,0.0d0]
      aa2=latcons*[0.0d0,1.0d0,0.0d0]
      
      q=[1.d0,1.d0,0.d0]
      q=q*alpha/latcons
      if(abs(Jnn)<abs(Dnn)) q=(alpha+pi)/latcons
      
      q=[1.d0,1.d0,0.d0]
      q = q*pi/2.d0
      phi = pi/4.d0

      c1 = cos(dot_product(q,aa1))
      c2 = cos(dot_product(q,aa2))
      s1 = sin(dot_product(q,aa1))
      s2 = sin(dot_product(q,aa2))
      cphi = cos(phi)
      sphi = sin(phi)
      k_a1 = dot_product(k,aa1)
      k_a2 = dot_product(k,aa2)

!     Dpp = -S*(  cos(dot_product(k,aa1))*(Jnn*c1+Dnn*s1+Jnn) + cos(dot_product(k,aa2))*(Jnn*c2+Dnn*s2+Jnn) - 2.0*Jnn*(c1+c2) - 2.d0*Dnn*(s1+s2)  )
!     Dmp = -S*(  cos(dot_product(k,aa1))*(Jnn*c1+Dnn*s1-Jnn) + cos(dot_product(k,aa2))*(Jnn*c2+Dnn*s2-Jnn)  )
!     omega=sqrt( Dpp**2 - Dmp**2 )

!     prefact1pp = (1.d0-ii)*Dnn/sqrt(2.d0)
!     prefact2pp = (1.d0+ii)*Dnn/sqrt(2.d0)

!     prefact1pm = (1.d0-ii)*Dnn/sqrt(2.d0)
!     prefact2pm = (1.d0+ii)*Dnn/sqrt(2.d0)

!     J00 = 2.d0*Jnn*(c1+c2) + 2.d0*Dnn*(s1*cphi+s2*sphi)
!     J00 = 2.d0*Jnn*(c1+c2) + sqrt(2.d0)*Dnn*(s1+s2)
!     J00 = 2.d0*sqrt(2.d0)*Dnn


      prefact1pp =      Jnn+Dnn*cos(phi)-ii*Dnn*sin(phi)
      prefact2pp =      Jnn+Dnn*sin(phi)+ii*Dnn*cos(phi)
      prefact1pm =-(   -Dnn+Jnn*cos(phi)+ii*Jnn*sin(phi))*(cos(phi)-ii*sin(phi))
      prefact2pm =-(-ii*Dnn+Jnn*cos(phi)+ii*Jnn*sin(phi))*(cos(phi)-ii*sin(phi))
      J00 = 2.d0*Dnn*(cos(phi)+sin(phi))


      Dpp = cos(k_a1)*prefact1pp + cos(k_a2)*prefact2pp - J00
      Dpm = cos(k_a1)*prefact1pm + cos(k_a2)*prefact2pm 

      mat(1,:) = [-Dpp,-Dpm]  
      mat(2,:) = [conjg(Dpm),conjg(Dpp)]

      !Diagonalizing the D(k) matrix
      call ZGEEV('V','V', 2, mat(:,:), 2, evalues, leftevector, 2, rightevector, 2, work, 2*2, rwork, info)

      omega = evalues

! !Spin spiral
!     alpha=atan2(Dnn,Jnn)
!     a=1.d0
!     aa1=a*[1.0d0,0.0d0,0.0d0]
!     aa2=a*[0.0d0,1.0d0,0.0d0]

!     qq=alpha/a
!     if(abs(Jnn)<abs(Dnn)) qq=(alpha+pi-2.d0*pi)/a
!     JD=Jnn*cos(qq*a)+Dnn*sin(qq*a)

!     omega(1)= 2.0d0*sqrt( JD* (1.0d0-cos(k(1)*a)) * (JD-Jnn*cos(k(1)*a)) )
!     omega(2)=-2.0d0*sqrt( JD* (1.0d0-cos(k(1)*a)) * (JD-Jnn*cos(k(1)*a)) )

! !end spin spiral
   end function omega

!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
   subroutine sort(values,leftevector,rightevector,matrix_size)
      implicit none
      !It sorts the values contained in "values" in a decreasing order
      integer, intent(in) :: matrix_size
      complex(kind=pc), intent(inout) :: values(matrix_size), leftevector(matrix_size,matrix_size), rightevector(matrix_size,matrix_size)
      complex(kind=pc) :: valsaux(matrix_size), leftevector_init(matrix_size,matrix_size), rightevector_init(matrix_size,matrix_size)
      integer :: i, j, k, newpos(matrix_size)

      do i = 1, matrix_size
        do j = 1, i-1
          ! if( dble(values(i))>dble(valsaux(j)) ) then
          ! if( ( dble(values(i))-aimag(values(i)) ) > ( dble(valsaux(j))-aimag(valsaux(j)) ) ) then
          if( ( dble(values(i))+aimag(values(i)) ) > ( dble(valsaux(j))+aimag(valsaux(j)) ) ) then
            do k = i-1, j, -1
              valsaux(k+1)=valsaux(k)
              newpos(k+1)=newpos(k)
            end do
            exit
          end if
        end do
        valsaux(j)=values(i)
        newpos(j)=i
      end do
      
      values=valsaux

      rightevector_init = rightevector
      leftevector_init = leftevector
      do i = 1, matrix_size
         rightevector(:,i) = rightevector_init(:,newpos(i))
         leftevector(:,i) = leftevector_init(:,newpos(i))
      end do
   end subroutine sort

!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
   function spectralfunction(values,evector,k)
      implicit none
      real(kind=pc) :: spectralfunction(nptomega,2)
      complex(kind=pc), intent(in) :: values(naucell), evector(naucell,naucell)
      real(kind=pc), intent(in) :: k(3)
      real(kind=pc) :: omega, diff(3), incrementomega
      complex(kind=pc) ::  sum1, sum2, auxevector(naucell,naucell)
      integer :: n, i, j, p, mu, nasmallucell=4
      ! integer :: n, i, j, p, mu, nasmallucell=64

      ! auxevector=trp_cjg(evector)
      auxevector=evector

      incrementomega = (maxomega - minomega)/dble(nptomega-1)
      do p = 1, nptomega !Number of point in the y direction
         omega = minomega + incrementomega*dble(p-1)
         sum1 = 0.d0
         do n = 1, naucell
            sum2 = 0.d0
            do mu = 0,naucell-nasmallucell,nasmallucell
               do i = 1+mu, nasmallucell+mu; do j=1+mu, nasmallucell+mu
                 diff = basis(j,:) - basis(i,:)
                 sum2 = sum2 + exp( ii * dot_product(k, diff) ) * auxevector(i,n) * conjg(auxevector(j,n))
                 ! print *, "evector: ", evector(i,n); pause
               end do; end do !i, j
            end do !mu
            ! print *, "Sum2:", sum2; pause
            sum1 = sum1 + sum2/(omega-values(n)+ii*eta)
         end do !n
         ! print *,  "Sum1:", sum1; pause
         spectralfunction(p,:) = [ omega, -dimag(sum1)/(dble(naucell)*pi) ]
      end do !p, Number of point in the y direction
   end function spectralfunction

!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
   subroutine occupationnumber(evector,evalues,mode)
      implicit none
      complex(kind=pc), intent(in) :: evector(twonaucell,twonaucell), evalues(twonaucell)
      integer, intent(in) :: mode
      real(kind=pc) :: sumall, occupation, imag_occupation, real_occupation
      integer :: j, kk, pp, modeindex

      open(unit=99,file="occupation.dat",status="old")
      open(unit=1001,file="eigenvector.dat",status="old")

      sumall = 0.d0
      ! do j = 1, twonaucell
      do j = 1, naucell
         occupation = 0.d0
         do modeindex = 1, 64
            occupation = occupation + ABS(evector(modeindex,j+naucell))**2
         end do
         ! occupation = ABS(evector(mode,j))**2 + ABS(evector(mode,j+naucell))**2
         real_occupation = real(evector(naucell+1-mode,j))
         imag_occupation = aimag(evector(naucell+1-mode,j))

         sumall = sumall + occupation
         do kk=-dims(1),dims(1); do pp=-dims(2),dims(2)
            write(unit=99,fmt= "(6f20.12)") basis(j,:)+(big_a1*kk+big_a2*pp), occupation, real_occupation, imag_occupation
            write(unit=1001,fmt= "(6f20.12)") basis(j,:)+(big_a1*kk+big_a2*pp), atan2(imag_occupation, real_occupation ), 1.57079632679490, 1.d0
         end do; end do
      end do
      print *, "Magnon occupation number calculations of mode: ", mode
      print *, "Sum of the occupation number of all sites: ", sumall
      print *, "    k: ", kpoints(1,:)
      print *, "k/2pi: ", kpoints(1,:)/(2.d0*pi)
      print *, "Mode energy: ", evalues(naucell+1-mode)
      stop
   end subroutine occupationnumber

!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
   subroutine printtime(message)
!$    use omp_lib
      implicit none
!$    integer :: nthreads
      real(kind=pc) ::  end_cpu_time, timedif
      character(len=*), intent(in) :: message
            
!$    nthreads = omp_get_num_threads()
      call cpu_time (end_cpu_time)
!$    end_cpu_time = omp_get_wtime()
      timedif=end_cpu_time - beg_cpu_time
!! $    timedif=timedif/nthreads
      print "(a,': running time = ',i0,'h:',i0,'m:',i0,'s')", message, &
                                       int(timedif/3600.0), & 
                                       int(mod(timedif,3600.0)/60.0),&
                                       int(mod(mod(timedif,3600.0),60.0)) 

      return
   end subroutine printtime

!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
   subroutine getriddegeneracy(i,evalues,rightevector,halfevalues,halfrightevector)
      implicit none
      integer, intent(in) :: i
      complex(kind=pc), intent(in)  :: evalues(twonaucell), rightevector(twonaucell,twonaucell)
      complex(kind=pc), intent(out) :: halfevalues(naucell), halfrightevector(naucell,naucell)
      integer :: m, j

      m = 0
      do j = 1, twonaucell
         if( sum(abs(rightevector(1:naucell,j))**2) > sum(abs(rightevector(naucell+1:twonaucell,j))**2) ) then
            m = m + 1
            halfevalues(m) = evalues(j)
            halfrightevector(:,m) = rightevector(1:naucell,j)
         end if
         if( m == naucell ) exit
      end do
      if( m < naucell ) print *, "Warning: Degeneracy processing failed! kpoint i=", i
   end subroutine getriddegeneracy

!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
   ! subroutine subunfolding(i,eigenvalues,eigenvector,spectra)
!$ subroutine subunfolding(i,eigenvalues,eigenvector,spectra,lck)
!$    use omp_lib
      implicit none
      integer, intent(in) :: i
      complex(kind=pc), intent(in) :: eigenvalues(naucell), eigenvector(naucell,naucell)
      real(kind=pc), intent(out) :: spectra(nptomega,2)
!$      integer(kind = OMP_lock_kind) :: lck

      if( i == 1 ) then
         if(maxomega == 0.123456789d0) maxomega = rescalf*dble(eigenvalues(naucell))!The max. omega will be proportional to the highest eigen-energy
         if(minomega == 0.123456789d0) minomega = dble(eigenvalues(1)) - abs(rescalf-1.d0)*dble(eigenvalues(naucell))
         ! eta = eta*maxomega
         maxomegaOK = .true.
      else
         do while(.not.maxomegaOK)
!$          call OMP_set_lock(lck)
!$          call OMP_unset_lock(lck)
         end do
      end if

      spectra(:,:) = spectralfunction(eigenvalues,eigenvector,kpoints(i,:))
   end subroutine subunfolding

!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
!$ subroutine unfoldingnoncol(kpointindex,eigenvalues,eigenvector,k,spectra,lck)
   ! subroutine unfoldingnoncol(kpointindex,eigenvalues,eigenvector,k,spectra)
!$    use omp_lib
      implicit none
      integer, intent(in) :: kpointindex
      complex(kind=pc), intent(in) :: eigenvalues(naucell), eigenvector(2,2,naucell,naucell)
      real(kind=pc), intent(in) :: k(3)
      real(kind=pc), intent(out) :: spectra(nptomega,13)!, expected_valueS(naucell,3)
      integer :: r, p, mu, nu, alpha, beta, sii, ssf, polariz, i, j
      real(kind=pc) :: omega, diff(3), incrementomega, polarizations(3,2)!, tempD(4,8)
      complex(kind=pc) :: N_alp_bet, N_mu_nu, gamma_scat(3,2,2)!, F_mu_r(3), expected_aux(3)
      real(kind=pc) :: C_alpha_beta, delta_alpha_beta
!$    integer(kind = OMP_lock_kind) :: lck
      ! character(len=50) :: suffix, formt
      complex(kind=pc) :: factor(3,3), Matrix(3,3), factor_alpha_beta
      integer :: alpha_, beta_

      ! Calculate factor to account for the dipole interaction between neutron and electron in neutron scattering measurements.
      ! The neutron cross section includes a summation over x,y,z like
      ! sum_{alpha,beta} (delta_{alpha, beta} - k_alpha*k_beta/k^2) S_{alpha, beta}
      ! However, in my calculations, I sum over +, -, z. 
      ! I determinined that sum_{alpha,beta} C_{alpha,beta}*S_alpha*S_beta = sum_{gamma,eta}*S_gamma*S_eta *sum_{alpha,beta} C_{alpha, beta} M_{alpha,gamma}*M_{beta,eta}
      ! where alpha, beta run over x,y,z and gamma, eta over +, -, z. And M = [[1,1,0], [-i, i, 0], [0, 0, 1]]
      ! In the code below, however, we use variable 'alpha', 'beta' for gamma and eta, and 'alpha_', 'beta_' for alpha and beta in the equations above.
      ! M derives from the transformation S^+- = 1/2 (S^x +- i S^y) actually from S^x = S^+ - S^-, S^y = -i(S^+ - S^-)   ==>  vec(Sx) = M vec(S+)

      if(circular_cartesian_convertion_factor) then
         Matrix(1,:) = [ cone,  cone, czero]
         Matrix(2,:) = [-ii  ,  ii  , czero]
         Matrix(3,:) = [czero, czero,  cone]
      
         factor = 0.d0
         do alpha = 1, 3; do beta = 1, 3
            factor_alpha_beta = 0.d0
            do alpha_=1, 3; do beta_=1,3
               
               if( alpha_ == beta_ ) then
                  delta_alpha_beta = 1.0d0
               else
                  delta_alpha_beta = 0.0d0
               end if

               if(neutron_factor) then
                  C_alpha_beta = ( delta_alpha_beta - k(alpha_)*k(beta_)/norm2(k) )
               else
                  C_alpha_beta = 1.d0
               end if

               factor_alpha_beta = factor_alpha_beta + Matrix(alpha_,alpha)*Matrix(beta_,beta) * C_alpha_beta
            end do; end do
            factor(alpha,beta) = factor_alpha_beta
         end do; end do
      else
         factor = 1.d0
      end if

      if( kpointindex == 1 ) then
         if(maxomega == 0.123456789d0) maxomega = rescalf*dble(eigenvalues(naucell))!The max. omega will be proportional to the highest eigen-energy
         if(minomega == 0.123456789d0) minomega = dble(eigenvalues(1)) - abs(rescalf-1.d0)*dble(eigenvalues(naucell))
         ! eta = eta*maxomega
         maxomegaOK = .true.

         ! write(format,fmt="(a,f7.3,a,f7.3,a)") "'s/set yrange.*/set yrange[", minomega, " :", maxomega, " ]/'"
         ! format="sed -i '' "//format(1:50)//" disp_unfol.gnu"
         ! call system(format)
      else
         do while(.not.maxomegaOK)
!$          call OMP_set_lock(lck)
!$          call OMP_unset_lock(lck)
         end do
      end if

      incrementomega = (maxomega - minomega)/dble(nptomega-1)
      do p = 1, nptomega !Number of points in the y direction
         omega = minomega + incrementomega*dble(p-1)

         gamma_scat = 0.d0
         do alpha = 1, 3; do beta = 1, 3
            N_alp_bet = 0.d0
            do mu = 1, naucell; do nu = 1, naucell
               diff = basis(nu,:) - basis(mu,:)
               N_mu_nu = 0.d0
               do r = 1, naucell
                  ! N_mu_nu = N_mu_nu + delta(omega-eigenvalues(r)) * (eigenvector(r,mu)) * conjg(eigenvector(r,nu))
                  ! N_mu_nu = N_mu_nu + delta(omega-eigenvalues(r)) * (eigenvector(1,1,mu,r)) * conjg(eigenvector(1,1,nu,r))

                  N_mu_nu = N_mu_nu + delta(omega-eigenvalues(r), eigenvalues(r))  * ( &
                             Rotpm(mu,alpha,2) * Rotpm(nu,beta, 1) * eigenvector(2,1,nu,r ) * conjg(  eigenvector(2,1,mu,r ) ) + &
                             Rotpm(mu,alpha,2) * Rotpm(nu,beta, 2) * eigenvector(1,1,nu,r ) * conjg(  eigenvector(2,1,mu,r ) ) + &
                             Rotpm(mu,alpha,1) * Rotpm(nu,beta, 1) * eigenvector(2,1,nu,r ) * conjg(  eigenvector(1,1,mu,r ) ) + &
                             Rotpm(mu,alpha,1) * Rotpm(nu,beta, 2) * eigenvector(1,1,nu,r ) * conjg(  eigenvector(1,1,mu,r ) )   &
                           )

                  ! N_mu_nu = N_mu_nu + delta(omega-eigenvalues(r))  * ( &
                  !            Rotpm(mu,alpha,1) * Rotpm(nu,beta, 1)  * eigenvector(2,2,mu,r ) * eigenvector(2,1,nu,r ) + &
                  !            Rotpm(mu,alpha,2) * Rotpm(nu,beta, 1)  * eigenvector(1,2,mu,r ) * eigenvector(2,1,nu,r ) + &
                  !            Rotpm(mu,alpha,1) * Rotpm(nu,beta, 2)  * eigenvector(2,2,mu,r ) * eigenvector(1,1,nu,r ) + &
                  !            Rotpm(mu,alpha,2) * Rotpm(nu,beta, 2)  * eigenvector(1,2,mu,r ) * eigenvector(1,1,nu,r )   &
                  !          )
               end do

               !These commented lines below were for a test for Manuel, about the local or global reference system for the spin wave wave vectores.
               ! if(alpha==1 .and. beta==2) then
               !    N_mu_nu = 2.d0 * sqrt( Si(mu) * Si(nu) ) * N_mu_nu
               ! else
               !    N_mu_nu = 0.d0
               ! end if

               ! N_mu_nu = 2.d0 * sqrt( Si(mu) * Si(nu) ) * Rotpm(mu,alpha,1)*Rotpm(nu,beta,2) * N_mu_nu
               N_mu_nu = 2.d0 * sqrt( Si(mu) * Si(nu) ) * N_mu_nu / naucell !I am trying to normalize the spectra intensity for any number of atoms in the unit cell
               N_alp_bet = N_alp_bet + exp( -ii*dot_product(k, diff) ) * Uweight(mu)*Uweight(nu) * N_mu_nu * factor(alpha,beta)
             end do; end do !mu, nu
            do sii=1,2; do ssf=1,2
               do polariz = 1, 3
                  gamma_scat(polariz,sii,ssf) = gamma_scat(polariz,sii,ssf) + paulimatrix(polariz,alpha,sii,ssf)*paulimatrix(polariz,beta,ssf,sii)*N_alp_bet
               end do
            end do; end do
         end do; end do  !alpha, beta
         gamma_scat = gamma_scat * 2.d0*pi!/dble(naucell)

         spectra(p,:) = [ omega, ( abs(gamma_scat(polariz,1,1)), abs(gamma_scat(polariz,1,2)), abs(gamma_scat(polariz,2,2)), abs(gamma_scat(polariz,2,1)), polariz = 1, 3 ) ]
      end do !p, Number of points in the y direction

      !The next lines prints on the screen the eigenvector and eigenvalues
      !It it useful to understand the structure them. For example, the eigenvector
      !should be block diagonal for FM, and with off-diagonal when D or non-collinearity
      !is present.

      ! print "(a,100f12.8)", "real E=", real(eigenvalues)
      ! print "(a,100f12.8)", "imag E=", aimag(eigenvalues)
      ! tempD(1:2,1:2) =  real(eigenvector(1,1,:,:))
      ! tempD(1:2,3:4) =  real(eigenvector(1,2,:,:))
      ! tempD(3:4,1:2) =  real(eigenvector(2,1,:,:))
      ! tempD(3:4,3:4) =  real(eigenvector(2,2,:,:))

      ! tempD(1:2,5:6) = aimag(eigenvector(1,1,:,:))
      ! tempD(1:2,7:8) = aimag(eigenvector(1,2,:,:))
      ! tempD(3:4,5:6) = aimag(eigenvector(2,1,:,:))
      ! tempD(3:4,7:8) = aimag(eigenvector(2,2,:,:))
      ! call printmatrix( "Real Rmatrix", tempD, "f12.4" )
      ! pause

   end subroutine unfoldingnoncol

!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
!$ subroutine inelastic_intensity(kpointindex,eigenvalues,eigenvector,k,ine_intensities)
!$    use omp_lib
      implicit none
      integer, intent(in) :: kpointindex
      complex(kind=pc), intent(in) :: eigenvalues(naucell), eigenvector(2,2,naucell,naucell)
      real(kind=pc), intent(in) :: k(3)
      real(kind=pc), intent(out) :: ine_intensities(naucell,13)!, expected_valueS(naucell,3)
      integer :: r, p, mu, nu, alpha, beta, sii, ssf, polariz, i, j
      real(kind=pc) :: omega, diff(3), incrementomega!, tempD(4,8)
      complex(kind=pc) :: N_alp_bet, N_mu_nu, gamma_scat(3,2,2)!, F_mu_r(3), expected_aux(3)
      real(kind=pc) :: C_alpha_beta, delta_alpha_beta
      complex(kind=pc) :: factor(3,3), Matrix(3,3), factor_alpha_beta
      integer :: alpha_, beta_

      ! Calculate factor to account for the dipole interaction between neutron and electron in neutron scattering measurements.
      ! The neutron cross section includes a summation over x,y,z like
      ! sum_{alpha,beta} (delta_{alpha, beta} - k_alpha*k_beta/k^2) S_{alpha, beta}
      ! However, in my calculations, I sum over +, -, z. 
      ! I determinined that sum_{alpha,beta} C_{alpha,beta}*S_alpha*S_beta = sum_{gamma,eta}*S_gamma*S_eta *sum_{alpha,beta} C_{alpha, beta} M_{alpha,gamma}*M_{beta,eta}
      ! where alpha, beta run over x,y,z and gamma, eta over +, -, z. And M = [[1,1,0], [-i, i, 0], [0, 0, 1]]
      ! In the code below, however, we use variable 'alpha', 'beta' for gamma and eta, and 'alpha_', 'beta_' for alpha and beta in the equations above.
      ! M derives from the transformation S^+- = 1/2 (S^x +- i S^y) actually from S^x = S^+ - S^-, S^y = -i(S^+ - S^-)   ==>  vec(Sx) = M vec(S+)

      if(circular_cartesian_convertion_factor) then
         Matrix(1,:) = [ cone,  cone, czero]
         Matrix(2,:) = [-ii  ,  ii  , czero]
         Matrix(3,:) = [czero, czero,  cone]
      
         factor = 0.d0
         do alpha = 1, 3; do beta = 1, 3
            factor_alpha_beta = 0.d0
            do alpha_=1, 3; do beta_=1,3
               
               if( alpha_ == beta_ ) then
                  delta_alpha_beta = 1.0d0
               else
                  delta_alpha_beta = 0.0d0
               end if

               if(neutron_factor) then
                  C_alpha_beta = ( delta_alpha_beta - k(alpha_)*k(beta_)/norm2(k) )
               else
                  C_alpha_beta = 1.d0
               end if

               factor_alpha_beta = factor_alpha_beta + Matrix(alpha_,alpha)*Matrix(beta_,beta) * C_alpha_beta
            end do; end do
            factor(alpha,beta) = factor_alpha_beta
         end do; end do
      else
         factor = 1.d0
      end if

      do r = 1, naucell !Number of points in the y direction
         omega = eigenvalues(r)

         gamma_scat = 0.d0
         do alpha = 1, 3; do beta = 1, 3
            N_alp_bet = 0.d0
            do mu = 1, naucell; do nu = 1, naucell
               diff = basis(nu,:) - basis(mu,:)
               N_mu_nu = 0.d0
               ! do r = 1, naucell

                  N_mu_nu = N_mu_nu + delta(omega-eigenvalues(r), eigenvalues(r))  * ( &
                             Rotpm(mu,alpha,2) * Rotpm(nu,beta, 1) * eigenvector(2,1,nu,r ) * conjg(  eigenvector(2,1,mu,r ) ) + &
                             Rotpm(mu,alpha,2) * Rotpm(nu,beta, 2) * eigenvector(1,1,nu,r ) * conjg(  eigenvector(2,1,mu,r ) ) + &
                             Rotpm(mu,alpha,1) * Rotpm(nu,beta, 1) * eigenvector(2,1,nu,r ) * conjg(  eigenvector(1,1,mu,r ) ) + &
                             Rotpm(mu,alpha,1) * Rotpm(nu,beta, 2) * eigenvector(1,1,nu,r ) * conjg(  eigenvector(1,1,mu,r ) )   &
                           )

                  ! N_mu_nu = N_mu_nu + delta(omega-eigenvalues(r))  * ( &
                  !            Rotpm(mu,alpha,1) * Rotpm(nu,beta, 1)  * eigenvector(2,2,mu,r ) * eigenvector(2,1,nu,r ) + &
                  !            Rotpm(mu,alpha,2) * Rotpm(nu,beta, 1)  * eigenvector(1,2,mu,r ) * eigenvector(2,1,nu,r ) + &
                  !            Rotpm(mu,alpha,1) * Rotpm(nu,beta, 2)  * eigenvector(2,2,mu,r ) * eigenvector(1,1,nu,r ) + &
                  !            Rotpm(mu,alpha,2) * Rotpm(nu,beta, 2)  * eigenvector(1,2,mu,r ) * eigenvector(1,1,nu,r )   &
                  !          )
               ! end do

               N_mu_nu = 2.d0 * sqrt( Si(mu) * Si(nu) ) * N_mu_nu / naucell !I am trying to normalize the spectra intensity for any number of atoms in the unit cell
               N_alp_bet = N_alp_bet + exp( -ii*dot_product(k, diff) ) * Uweight(mu)*Uweight(nu) * N_mu_nu * factor(alpha,beta)
            end do; end do !mu, nu
            do sii=1,2; do ssf=1,2
               do polariz = 1, 3
                  gamma_scat(polariz,sii,ssf) = gamma_scat(polariz,sii,ssf) + paulimatrix(polariz,alpha,sii,ssf)*paulimatrix(polariz,beta,ssf,sii)*N_alp_bet
               end do
            end do; end do
         end do; end do  !alpha, beta
         gamma_scat = gamma_scat * 2.d0*pi!/dble(naucell)

         ine_intensities(r,:) = [ omega, ( abs(gamma_scat(polariz,1,1)), abs(gamma_scat(polariz,1,2)), abs(gamma_scat(polariz,2,2)), abs(gamma_scat(polariz,2,1)), polariz = 1, 3 ) ]
      end do !p, Number of points in the y direction

   end subroutine inelastic_intensity

!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
   ! Calculating the spin-wave mode angular momentum expected value of the spin operators
   subroutine angular_momentum(kpointindex,eigenvector,k,expected_valueS)
      implicit none
      integer, intent(in) :: kpointindex
      complex(kind=pc), intent(in) :: eigenvector(2,2,naucell,naucell)
      real(kind=pc), intent(in) :: k(3)
      complex(kind=pc), intent(out) :: expected_valueS(naucell,3)
      integer :: r, mu
      real(kind=pc) :: Rotmat_aux(3,3)
      complex(kind=pc) :: N_mu_nu, F_mu_r(3), expected_aux(3)

      expected_valueS = 0.d0
      do r = 1, naucell
         expected_aux = 0.d0
         do mu = 1, naucell
            ! ! N_mu_nu =  - eigenvector(1,2,mu,r)*eigenvector(2,1,mu,naucell+1-r) - eigenvector(1,1,mu,r)*eigenvector(2,2,mu,naucell+1-r)
            ! N_mu_nu =  - eigenvector(2,1,mu,r)*conjg(eigenvector(2,1,mu,r)) - eigenvector(1,1,mu,r)*conjg(eigenvector(1,1,mu,r))
            ! ! N_mu_nu =  - eigenvector(1,2,mu,r)*conjg(eigenvector(1,2,mu,r)) - eigenvector(1,1,mu,r)*conjg(eigenvector(1,1,mu,r))

            ! N_mu_nu =  - eigenvector(1,1,mu,r)*eigenvector(2,2,mu,r) - eigenvector(1,2,mu,r)*eigenvector(2,1,mu,r)
            N_mu_nu =  - eigenvector(1,1,mu,r)*conjg(eigenvector(1,1,mu,r)) - eigenvector(1,2,mu,r)*conjg(eigenvector(1,2,mu,r))

            Rotmat_aux(:,:) = Rotmat(mu,:,:)

            !local expected value vector
            F_mu_r = [ czero, czero, N_mu_nu ] 
            !Rotation to the global reference frame
            F_mu_r = matmul(Rotmat_aux, F_mu_r)

            expected_aux = expected_aux + F_mu_r
         end do

         expected_valueS(r,:) =  expected_aux
      end do
   end subroutine angular_momentum

   function delta(omega, energy)
      complex(kind=pc), intent(in) :: omega, energy
      real(kind=pc) :: delta
      ! delta = ((eta*energy)/pi) / (omega**2 + (eta*energy)**2) * energy
      delta = (eta/pi) / (omega**2 + eta**2) 
   end function

   subroutine outputdata(disp_matrix,spectra,ine_intensities,expected_valueS,rightevector,evalues)
      !Storage making
      implicit none
      real(kind=pc), intent(in) :: disp_matrix(effnkpt,twonaucell), spectra(effnkpt,nptomega,13), ine_intensities(effnkpt,naucell,13)
      complex(kind=pc), intent(in) :: rightevector(effnkpt,twonaucell,twonaucell), evalues(effnkpt,twonaucell), expected_valueS(effnkpt,naucell,3)
      integer :: i, j, counter, m, n
      real(kind=pc) :: k(3), prevk(3), path, units=2.d0*pi/1.d0, spectra_aux(13), disp_matrix_aux(twonaucell), Re_expected_valueS(3), Im_expected_valueS(3)
      character(len=70) :: formt, formt2, suffix

      i=index(inputcardname,"_")+1
      j=index(inputcardname,".")-1
      suffix = inputcardname(i:j)

      !Storage file openning
      formt="dispersion_"//trim(suffix)//".dat"
      open( unit=90, file=formt )

      formt="disp_analy_"//trim(suffix)//".dat"
      open( unit=91, file=formt )

      formt="disp_unfol_"//trim(suffix)//".dat"
      open( unit=92, file=formt )

      formt="ine_intensities_"//trim(suffix)//".dat"
      open( unit=93, file=formt )
      
      formt="precession_"//trim(suffix)//".dat"
      open( unit=123, file=formt )

      formt="kpath_"//trim(suffix)//".dat"
      open( unit=333, file=formt )

      formt="outputfiles/angular_momentum_"//trim(suffix)//".dat"
      open( unit=987, file=formt )

      do i=-2,2; do j=-2,2
         write(333, fmt=*) (b1*i + b2*j)
      end do; end do
      write(333, fmt=*) "#"

      write(formt,fmt="(a,i0,a)") "(", twonaucell+1, "f16.8)"
      write(formt2,fmt="(a,i0,a)") "(", 2*((2*ncpdim+1)**2+1), "f16.8)"
      
      write( unit=123, fmt="(a)") "#Kveci  Mode    Energies      K-vector x    y             z              Re R++1       Im R++1       Re R--1       Im R--1        Re R++2       Im R++2       Re R--2       Im R--2    "

      if(.not. constant_energy_plot ) then
         path = 0.d0
         counter = 0
         do i = 1, effnkpt
            k = kpoints(i,:)
            write( unit=333, fmt=* ) k
            path = path + dkpoints(i)
            !File "dispersion_**.dat" will contain the result via the numberical solution
            disp_matrix_aux = disp_matrix(i,:)
            write( unit=90, fmt=formt ) path/units, disp_matrix_aux
            
            !File "dispersion_analytics.dat" contains the analytical solution
            if(analytics) write( unit=91, fmt=formt2) path/units, ((real(omega(k+m*b1+n*b2)),m=-ncpdim,ncpdim),n=-ncpdim,ncpdim)
            
            !File "disp_unfol_XX.dat" contains the unfolded spectral function
            do j = 1, nptomega
               spectra_aux =  spectra(i,j,:)
               write( unit=92, fmt= "( 14(f25.15,' ') )" ) path/units, spectra_aux 
            end do
            write( unit=92, fmt= * )    

            !File "ine_intensities_XX.dat" contains the inelastic intensities
            do j = 1, naucell
               spectra_aux =  ine_intensities(i,j,:)
               write( unit=93, fmt= "( 14(f25.15,' ') )" ) path/units, spectra_aux 
            end do
            write( unit=93, fmt= * )      
            
            !Saving the results for the high symmetry points to plot on the screen
            n = 0
            do j=1, npath
               n = n + nkpt_npath(j)
               if( i == n+1 ) then
                  counter = counter + 1
                  highsymm(counter,:) = [path/units, disp_matrix(i,:)]
               end if
            end do
            if( i == effnkpt ) then
               counter = counter + 1
               highsymm(counter,:) = [path/units, disp_matrix(i,:)]
            end if
            
            ! Writing down the eigenvector for the precession visualization
            do mode = 1, naucell
               n = naucell+1 - mode
               write( unit=123, fmt="(2i6, 4f14.8,' ')", advance="no" ) i, mode, dble(evalues(i,n))+aimag(evalues(i,n)), k
               do m = 1, 1 !This was used to repeat the unit cell info for a crystal construction, but this is done on the visualization script now
                  do j = 1, naucell 
                     write( unit=123, fmt="(4f14.8, ' ')", advance="no" ) &
                        real(rightevector(i,j,n)),         aimag(rightevector(i,j,n)), &
                        real(rightevector(i,j+naucell,n)), aimag(rightevector(i,j+naucell,n)) 
                  end do
               end do
               write( unit=123, fmt=* )! spectra(i,j,:)
            end do

            !Angular momentum
            write(unit=987, fmt="(a,i0,a,3f14.6)") "kpointindex: ", i, " K=", k
            do mode = 1, naucell
               Re_expected_valueS = real(  expected_valueS(i,mode,:) )
               Im_expected_valueS = aimag( expected_valueS(i,mode,:) )
               write(unit=987, fmt="(a,i3,a,3f12.8,a,3f12.8,a,1f8.4,a,1f8.4,a,2f12.8)") "   mode= ", mode, "   Re <S_r>", Re_expected_valueS, "   Im <S_r>", Im_expected_valueS, "   Re|<S_r>|", norm2(Re_expected_valueS), "   Im|<S_r>|", norm2(Im_expected_valueS) , "   E_r= ", evalues(i,mode) 
            end do
            write(unit=987, fmt=*) 

         end do !Loop of k points

      else ! If constant_energy_plot = .true.
         formt="outputfiles/constant_energy_"//trim(suffix)//".csv"
         open( unit=92, file=formt )
         write( unit=92, fmt="(a)" ) " x, y, radius, omega, px_ch1, px_ch2, px_ch3, px_ch4, py_ch1, py_ch2, py_ch3, py_ch4, pz_ch1, pz_ch2, pz_ch3, pz_ch4  "

         do i=1, effnkpt
            k = kpoints(i,:)
            spectra_aux =  spectra(i,1,:)

            write( unit=92, fmt="(16(f14.8,','))" ) k(1:2), norm2(k(1:2)), spectra_aux 
         end do
      end if 

      write(333, fmt=*) "# ^ Reciprocal space lattice"

      write(333, fmt=*) 0.0, 0.0, 0.0 
      write(333, fmt=*) b1
      write(333, fmt=*) "# ^ b1"
      write(333, fmt=*) 0.0, 0.0, 0.0 
      write(333, fmt=*) b2 
      write(333, fmt=*) "# ^ b2"

      do i=-1,1; do j=-1,1
         k = (b1*i + b2*j)
         if(norm2(k) <= norm2(b1)*1.2 .and. norm2(k) <= norm2(b2)*1.2 ) then
            prevk = [-k(2), k(1), k(3)]
            write(333, fmt=*) 0.5*k - 0.5*prevk   
            write(333, fmt=*) 0.5*k + 0.5*prevk   
            write(333, fmt=*) "# ^ edge of the 1st BZ"
         end if
      end do; end do
   end subroutine outputdata

   subroutine printmatrix(title,A,format)
      implicit none
      real(kind=pc), intent(in) :: A(:,:)
      character(len=*), intent(in) :: format, title
      character(len=5) :: leftform 
      integer :: lin, col
      integer :: i, j

      lin=size(A(:,1))
      col=size(A(1,:))
      
      write(leftform,fmt="(a,i0)") '(', col

      print *, title
      do i=1, lin
         print leftform // format // ')', ( A(i,j), j=1, col )
      end do
   end subroutine printmatrix

   function trp_cjg(A) result(trp_cjgA)
      implicit none
      complex(pc), dimension(:,:), intent(in) :: A
      complex(pc), dimension(size(A,1),size(A,2)) :: trp_cjgA
      trp_cjgA = transpose(conjg(A))
   end function trp_cjg

   function trimatmul(A,B,C) result(ABC)
      implicit none
      complex(pc), dimension(:,:), intent(in) :: A, B, C
      complex(pc), dimension(size(A,1),size(A,2)) :: ABC 
      ABC = matmul(A, matmul(B,C) )
   end function trimatmul
end module mod_global
