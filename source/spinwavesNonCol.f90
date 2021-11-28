!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
!!!!!!!!!!!!!!!!!!!!!!!!!!!! MAIN PROGRAM !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
program dispersion
   use mod_global
!$ use omp_lib
   implicit none
   external :: bogoliubov_transf
!$ integer :: nthreads,mythread
!$ integer(kind = OMP_lock_kind) :: lck
   integer :: i,j,m, info=0, SDIM 
   real(kind=pc) :: k(3), summation
   complex(kind=pc), allocatable :: Dmatrix(:,:,:), evalues(:), evaluesTOT(:,:), leftevector(:,:), rightevector(:,:), rightevectorTOT(:,:,:), halfevector(:,:), evectorUnfol(:,:,:,:), halfevalues(:), diagonal(:,:), nodiagonal(:,:) !Dmatrix is the one to be diagonalized to provide the spinwave dispersion. evalues store the eigen energies from this diagonalization, left/rightevector store the eigenvectors.
   complex(kind=pc), allocatable :: work(:) !work space for the diagonalization routine.
   double precision, allocatable :: rwork(:), rwork2(:) !work space for the diagonalization routine.
   logical, allocatable :: bwork(:) !work space for the diagonalization routine.
   real(kind=pc), allocatable :: spectra(:,:,:), spectra_aux(:,:), ine_intensities(:,:,:), ine_intensities_aux(:,:), disp_matrix(:,:)
   complex(kind=pc), allocatable :: expected_valueS(:,:,:), expected_valueS_aux(:,:)
   character(len=70) :: formt
   complex(kind=pc) :: number
   logical :: SELECT

   call cpu_time (beg_cpu_time)
!$ beg_cpu_time = omp_get_wtime()

   call getarg(1, inputcardname)
   if(inputcardname=="") stop "No inputcard file name informed. Please, give one as argument to the program. Stopping..."

   !Initializations: Reads inputcard and initializes variables and secondary parameters
   call initialization()
   allocate( Dmatrix(effnkpt,twonaucell,twonaucell) )
   ! allocate( evalues(twonaucell), evaluesTOT(effnkpt,twonaucell), halfevalues(naucell), leftevector(twonaucell,twonaucell), rightevector(twonaucell, twonaucell), rightevectorTOT(effnkpt, twonaucell, twonaucell), halfevector(naucell, naucell), evectorUnfol(2,2,naucell,naucell), work(2*twonaucell), rwork(2*twonaucell), rwork2(twonaucell), spectra(effnkpt,nptomega,5), disp_matrix(effnkpt,naucell) )
   ! allocate( diagonal(twonaucell,twonaucell), nodiagonal(twonaucell,twonaucell), bwork(twonaucell) )
   allocate( rightevectorTOT(effnkpt, twonaucell, twonaucell), evaluesTOT(effnkpt,twonaucell), spectra(effnkpt,nptomega,13), spectra_aux(nptomega,13), ine_intensities(effnkpt,naucell,13), ine_intensities_aux(naucell,13), expected_valueS(effnkpt,naucell,3), expected_valueS_aux(naucell,3), disp_matrix(effnkpt,twonaucell), highsymm(npath,naucell+1) )

!$ call omp_init_lock(lck)
   maxomegaOK = .false.

   !Calculates the D(k) matrix, which eigen-values correspond to the spinwaves dispersion
   call commutationmatrix(Dmatrix)

   !Loop over k points
!$omp parallel default(none) &
!$omp& private(mythread,i,k,nodiagonal,diagonal,evalues,leftevector,rightevector,work,rwork,rwork2,info,m,j,halfevalues,spectra_aux,ine_intensities_aux,expected_valueS_aux,halfevector,evectorUnfol,formt,number,summation,bwork,SDIM) &
!$omp& shared(nthreads,lck,effnkpt,kpoints,Dmatrix,twonaucell,naucell,calc_occup,mode,unfolding,maxomega,rescalf,eta,maxomegaOK,spectra,ine_intensities,expected_valueS,disp_matrix,rightevectorTOT,evaluesTOT,zero_toler,g,SELECT,nptomega,Tneel)
!$ mythread = omp_get_thread_num()
!$ if(mythread.eq.0) then
!$    nthreads = omp_get_num_threads()
!$    write(*,"(' OPENMP Number of threads: ',i0)") nthreads
      print *
      call printtime("End of Initiatization. Starting k-point loop")
      write(formt,fmt="(a,i3,a,i4,a,3f7.3)") " Progress: ",0*100/effnkpt, " %  K-point: ", 1, "      k: ", kpoints(1,:)
      call printtime(formt)
!$ end if
   allocate( evalues(twonaucell), halfevalues(naucell), leftevector(twonaucell,twonaucell), rightevector(twonaucell, twonaucell),  halfevector(naucell, naucell), evectorUnfol(2,2,naucell,naucell), work(2*twonaucell), rwork(2*twonaucell), rwork2(twonaucell) )
   allocate( diagonal(twonaucell,twonaucell), nodiagonal(twonaucell,twonaucell), bwork(twonaucell) )
!$omp do schedule(dynamic)
   do i = 1, effnkpt
      k = kpoints(i,:)
      
      nodiagonal = Dmatrix(i,:,:)
      !Diagonalizing the D(k) matrix
      ! call zgees  ('V','N', SELECT, twonaucell, Dmatrix(i,:,:), twonaucell, SDIM, evalues, rightevector, twonaucell, work, 2*twonaucell, rwork2, bwork, info)
      call zgeev('V','V', twonaucell, Dmatrix(i,:,:), twonaucell, evalues, leftevector, twonaucell, rightevector, twonaucell, work, 2*twonaucell, rwork, info)
      
      !Sorting in a crescent order the energy values
      if( .not. Tneel) call sort(evalues,leftevector,rightevector,twonaucell)

      !Makes the Bogoliubov transformation making from the eigenvectors new operators that obey bosons commutation relation
      ! if(i==2) &
      call bogoliubov_transf(i,nodiagonal,evalues,leftevector,rightevector)

      ! !Getting rid of the opposite energy degeneracy
      ! call getriddegeneracy(i,evalues,rightevector,halfevalues,halfevector)
      halfevalues = dble(evalues(1:naucell)) + aimag(evalues(1:naucell))
      halfevector = leftevector(1:naucell,1:naucell)
      ! halfevector=rightevector(1:naucell,1:naucell)

      evectorUnfol(1,1,:,:) = rightevector(1:naucell,1:naucell)
      evectorUnfol(2,1,:,:) = rightevector(1+naucell:twonaucell,1:naucell)
      evectorUnfol(1,2,:,:) = rightevector(1:naucell,1+naucell:twonaucell)
      evectorUnfol(2,2,:,:) = rightevector(1+naucell:twonaucell,1+naucell:twonaucell)

      !Unfolding. Spectral function calculation
      if(unfolding) &
!$       call unfoldingnoncol(i,halfevalues,evectorUnfol,k,spectra_aux,lck)
         ! call unfoldingnoncol(i,halfevalues,halfevector,k,spectra_aux)

!$       call inelastic_intensity(i,halfevalues,evectorUnfol,k,ine_intensities_aux)

!       if(unfolding) &
! !$       call subunfolding(i,halfevalues,halfevector,spectra_aux,lck)
!          call subunfolding(i,halfevalues,halfevector,spectra_aux)

      spectra(i,:,:) = spectra_aux
      ine_intensities(i,:,:) = ine_intensities_aux

      call angular_momentum(i,evectorUnfol,k,expected_valueS_aux)

      expected_valueS(i,:,:) = expected_valueS_aux

      !Calculation of the magnon occupation number
      if(i==1 .and. calc_occup ) &
         call occupationnumber(leftevector,evalues,mode)

      !Dispersion not unfolded
      disp_matrix(i,:) = dble(evalues) + aimag(evalues)

      !For the precession visualization
      rightevectorTOT(i,:,:) = rightevector
      evaluesTOT(i,:) = evalues

      !Informe on the terminal which points in being calculated (but once every 100 pts.)
      write(formt,fmt="(a,i3,a,i4,a,3f7.3)") " Progress: ",i*100/effnkpt, " %  K-point: ", i, "      k: ", k
      !This if allow to have this k points. For a single k point calculation, for example.
      if(effnkpt<5) then 
         call printtime(formt)
      else
         if(mod(i,effnkpt/5) == 0) call printtime(formt)
      end if

      ! deallocate( evalues, evaluesTOT, halfevalues, leftevector, rightevector, rightevectorTOT, halfevector, evectorUnfol, work, rwork, rwork2, spectra, disp_matrix )
      ! deallocate( diagonal, nodiagonal, bwork )

   end do !do i=1, effnkpt: Loop over k points
!$omp end do
!$omp end parallel

   !Storage making
   call outputdata( disp_matrix, spectra, ine_intensities, expected_valueS, rightevectorTOT, evaluesTOT )

   !Display the input parameters
   call printinginput()
 
!$omp parallel
!$ if(omp_get_thread_num()==0) then
   call printtime("PROGRAM END")
!$ end if
!$omp end parallel

end program dispersion



















