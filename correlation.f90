!==================MAIN PROGRAM===============================!
program linear_correlation
implicit none 
REAL,ALLOCATABLE,DIMENSION(:,:) :: COORD
REAL,ALLOCATABLE,DIMENSION(:,:) :: COVARMAT,TWOBODYCORR,EIGENVEC
REAL,ALLOCATABLE,DIMENSION(:)   :: COVARVEC,EIGENVAL,WORK
INTEGER :: natom,nframe,frame,i,j,k,jj,ii,threenatom,INFO
real :: temp, temp1, temp2, temp3
logical :: ESSENTIAL_DYN,MI  
! Three body variables
logical :: threebody
REAL,ALLOCATABLE,DIMENSION(:) :: MARGINAL
REAL,ALLOCATABLE,DIMENSION(:,:) :: AVERAGEPROD 
REAL,ALLOCATABLE,DIMENSION(:,:,:) :: COVARMAT3B,THREEBODYCORR
!Centrality
logical :: centrality
REAL,ALLOCATABLE,DIMENSION(:) ::  TWOBODYVEC
! RMSF
logical :: dormsf
real    :: rmsf 
real    :: rmtemp,rmtemp0
! MI-damp
logical :: damp
real,allocatable,dimension(:,:) :: distances
! READING VARIABLES
character (3)     :: char1
integer           :: int2
character (2)     :: char3
character (3)     :: char4
integer           :: int5
character (1)     :: char5
real           :: step  
real              :: StrengthI,  StrengthJ,  ConstraintI,  ConstraintJ,  Overlap,  CTE,  time,  mixE,  mixE1,  mixE2,  consE
real, allocatable :: coupling(:)
logical           :: couplingcorr
! animation
integer           :: nfr_animation
real,allocatable, dimension(:)   :: vector, wt
! Energies
integer  :: int3
real, allocatable, dimension(:,:) :: energy
real     :: re1, re2,re3, re4, re5
real     :: avgenerdiff,avgcoupling,maxenerdiff
logical  :: flip,skip
real     :: hbond
!.............................................................!
!.............................................................!

!*******************! 
    natom=   ! number of atoms in the system
    nframe=  ! number of frames in the trajectory
    couplingcorr=.false.
!*******************!
   ALLOCATE(coord(3*natom,nframe+1),COVARVEC(3*natom),coupling(nframe),energy(nframe,2),wt(nframe)) ! The possitions in each frame are stored and in the last possition of coord we store their average
  ALLOCATE(vector(3*natom))  

! Read the coordinates in xyz format
   write(*,*) 'Reading coordinates'
   open(1,file='qm.xyz') ! NAME OF TRAJ FILE
   open(2,file='E0.dat')
   open(3,file='E1.dat')
   do frame=1,nframe
      ii=1
      read(1,*)
      read(1,*)
      read(2,*) energy(frame,1)
      read(3,*) energy(frame,2)
      do i=1,natom
         read(1,*)  char5, coord(ii,frame), coord(ii+1,frame), coord(ii+2,frame)
         ii=ii+3
      enddo
   enddo
   close(1)

! Calculate the average coordinate 
  write(*,*) 'Calculating averages'
  call average(coord,natom,nframe) ! the average matrix is stored in the coord(:,nframe+1)

!! The coordenate vector is transformed into a displacement vector substracting the average from the coordinates i.e.  X = X - <X>
   write(*,*) 'Evaluating displacement vector (X - <X>)'
   do i=1,3*natom
      do j=1,nframe
         coord(i,j)=coord(i,j)-coord(i,nframe+1)
      enddo
   enddo

! Flip the energies
     flip=.false.
     open(unit=99,file='fliped-energies.txt')
     DO frame=1,nframe
        temp=energy(frame,2)-energy(frame,1)
        if(temp.lt.0.000) flip=.true.
        if(flip) temp=-temp
        if((mod(frame-1,500).eq.0).and.(frame.ne.1)) flip=.false.
 !       if(flip) write(99,*) frame
 !       write(99,*) frame,temp, energy(frame,1), energy(frame,2)
        if(flip) then
          re1=energy(frame,1)
          re2=energy(frame,2)
          energy(frame,2)=re1
          energy(frame,1)=re2
        endif
        write(99,*) frame, energy(frame,1), energy(frame,2)
     ENDDO

! compute energy-diff correlation vector

     avgenerdiff=0
     maxenerdiff=0
     DO frame=1,nframe
        temp=energy(frame,1)-energy(frame,2)
        if(maxenerdiff.lt.temp) maxenerdiff=temp
        avgenerdiff=avgenerdiff+temp
     ENDDO
     avgenerdiff=avgenerdiff/nframe
     ii=1
     DO i=1,3*natom
        re1=0.0000
        re2=0.0000
        DO frame=1,nframe
           re1=re1+coord(i,frame)**2
           re2=re2+(energy(frame,1)-energy(frame,2)+avgenerdiff)**2
        ENDDO
        DO frame=1,nframe
           re3=energy(frame,2)-energy(frame,1)+avgenerdiff
           re4=coord(i,frame)
           covarvec(i)=covarvec(i)+((re3*re4)/(sqrt(re1)*sqrt(re2)))
        ENDDO
        covarvec(i)=covarvec(i)/nframe
        if(mod(i,3).eq.0) ii=ii+3
     ENDDO

  open(unit=100,file='coordinate.txt')
  do i=1,3*natom
     write(100,*) covarvec(i)
  enddo

  covarvec=covarvec/norm2(covarvec)

! Write energy-diff mode animation
      nfr_animation=100
      write(*,*) nfr_animation
      call ANIMATE_EM(coord,covarvec,natom,nframe,nfr_animation)
      call system('mv animode-1.xyz animode-energy-diff.xyz')

! Project the coupling zeroth mean trajectory into our coordinate
      open(unit=333,file='projection.txt')
      open(unit=334,file='reactants.txt')
      open(unit=335,file='CI.txt')
      open(unit=336,file='products.txt')
      open(unit=337,file='h-bond-dyn.txt')
      flip=.false.
      skip=.false.
      DO frame=1,nframe-1
         read(337,*) re1,hbond
         DO i=1,3*natom
            vector(i)=coord(i,frame)
         ENDDO
         wt(frame)=dot_product(covarvec,vector)
         skip=.false.
         write(333,*) frame, wt(frame)
         temp=energy(frame,2)-energy(frame,1)
         temp1=energy(frame+1,2)-energy(frame+1,1)
         if(temp.lt.0.000) then
            flip=.true.
         else
            flip=.false.
         endif  
         if((temp1.lt.0.0000).and.(.not.flip)) skip=.true.          
         if(skip) write(335,*) temp, hbond !, wt(frame)
         if((flip).and.(.not. skip)) write(336,*) temp, hbond !, wt(frame)
         if((.not.flip).and.(.not.skip)) write(334,*) temp,hbond !, wt(frame)  
      ENDDO
      call PROJECTION_EM(coord,covarvec,natom,nframe,wt)
      call system('mv projection-md.xyz projection-md-energy-diff.xyz')
contains
!====================AUXILIAR FUNCTIONS======================!
       subroutine  PROJECTION_EM(coord,EM,natom,nframes,wt)
       implicit none
       REAL, PARAMETER      :: PI        = 3.1415926       ! some constants
       REAL,PARAMETER       :: amplitude = 0.1
       integer, INTENT(in)  :: nframes,natom
       REAL, INTENT(in)     :: coord(3*natom,nframes+1), EM(3*natom), wt(nframes)
       REAL                 :: oneovernsteps
       character(3)         :: char1
       integer              :: int2
       character(2)         :: char3
       character(3)         :: char4
       character(1)            :: char5
       character            :: algo
       integer              :: int5
       real                 :: re1,re2,re3,re4,re5,nada
!..............................................................................................................................!
       close(1000)
       close(3000)
       open(unit=1000, file='qm.xyz')
       open(unit=3000, file='projection-md.xyz')
       ii=1
!       read(1000,*)
!       write(3000,*) 'CRYST1    0.000    0.000    0.000  90.00  90.00  90.00 P 1 '
       DO i=1,nframes
          read(1000,*)
          read(1000,*)
          write(3000,*) natom
          write(3000,*) ''
          DO j=1,natom
             read(1000,*)  char1
             char5=char1
             re1=coord(ii,nframes+1) + EM(ii)*wt(i)
             re2=coord(ii+1,nframes+1) + EM(ii+1)*wt(i)
             re3=coord(ii+2,nframes+1) + EM(ii+2)*wt(i)
             ii=ii+3
             write(3000,*) char5,re1,re2,re3
          ENDDO
          ii=1
!          write(3000,*) 'END'
          rewind(1000)
       ENDDO

756  format(A4,2X,I5,2X,A4,A4,I5,4X,3F8.3,2F6.2,10X,2A2)

       return;end subroutine
!%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%!
       subroutine  ANIMATE_EM(coord,EM,natom,nframes,nsteps)
       implicit none
       REAL, PARAMETER      :: PI        = 3.1415926       ! some constants
       REAL,PARAMETER       :: amplitude = 5.0
       integer, INTENT(in)  :: nsteps,nframes,natom
       REAL, INTENT(in)     :: coord(3*natom,nframes+1), EM(3*natom)
       REAL                 :: oneovernsteps
       character(3)         :: char1
       integer              :: int2
       character(2)         :: char3
       character(3)         :: char4
       character(1)            :: char5
       character            :: algo
       integer              :: int5
       real                 :: re1,re2,re3,re4,re5,nada
!..............................................................................................................................!
       oneovernsteps=1.00000000000/nsteps
       open(unit=1000, file='qm.xyz')
       open(unit=3000, file='animode-1.xyz')
       ii=1
!       read(1000,*)
!       write(3000,*) 'CRYST1    0.000    0.000    0.000  90.00  90.00  90.00 P 1 '
       DO i=1,nsteps
          read(1000,*)
          read(1000,*)
!          write(*,*) sin(i*2*PI*oneovernsteps)
          write(3000,*) natom
          write(3000,*) ''
          DO j=1,natom
             read(1000,*)  char1
             char5=char1
!             write(99,*) EM(ii)*amplitude*sin(i*2*PI*oneovernsteps)
             re1=coord(ii,nframes+1) + EM(ii)*amplitude*sin(i*2*PI*oneovernsteps)
             re2=coord(ii+1,nframes+1) + EM(ii+1)*amplitude*sin(i*2*PI*oneovernsteps)
             re3=coord(ii+2,nframes+1) + EM(ii+2)*amplitude*sin(i*2*PI*oneovernsteps)
             ii=ii+3
             write(3000,*) char5,re1,re2,re3
          ENDDO
          ii=1
!          write(3000,*) 'END'
          rewind(1000)
       ENDDO

756  format(A4,2X,I5,2X,A4,A4,I5,4X,3F8.3,2F6.2,10X,2A2)

       return;end subroutine
!%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%!
      function COVAR(COORD,natom,nframe) result(COVARMAT)
      implicit none
!....................................................................................................!
!     Covariance matrix calculation:
!     This is maybe the simplest and most intuitive way to calculate the covariance matrix.   
!     This is currently the most demanding computation in the code and can be improved a lot.
!....................................................................................................!
      real, intent(in)    ::   COORD(:,:)
      integer, intent(in) ::   natom
      integer, intent(in) ::   nframe
      real,allocatable    ::   COVARMAT(:,:)
      real                ::   scratch,oneovernframe
      integer             ::   N,k,i,j 
     
      allocate(COVARMAT(3*natom,3*natom))
      oneovernframe=1.00000/nframe 
      scratch=0.00000
      do i=1,3*natom
         do j=1,i
           do k=1,nframe
               scratch=scratch+coord(i,k)*coord(j,k)
           enddo
           scratch=scratch*oneovernframe
           COVARMAT(i,j)=scratch
           COVARMAT(j,i)=scratch
           scratch=0.00000000
         enddo
       enddo
       return;end function 
!!------------------------------------------------------------!
       subroutine  AVERAGE(coord,natom,nframe) 
       implicit none
!....................................................................................................!
!      This subroutine calculates the average value of the coordinates and stores those values in coord(:,nframe+1)
!....................................................................................................!
       integer,intent(in)    :: natom,nframe
       real, intent(inout)      :: coord(3*natom,nframe+1)
       real                  :: avg,oneovernframe
       integer               :: i,j
       
       oneovernframe=1.0000/nframe
       do i=1,3*natom
          coord(i,nframe+1)=0.00000 
          avg=sum(coord(i,:))
          avg=avg/nframe
          coord(i,nframe+1)=avg
       enddo
     
       return;end subroutine
!!------------------------------------------------------------!
      function DET(COVARMAT,i,j) result(determinant)
       implicit none
!....................................................................................................!
! 3x3 determinant:
!      |a  b  c| 
!  det |d  e  f| = aei + bfg + cdh - ceg - bdi - afh 
!      |g  h  i|
! Carefull must be taken, i and j are indexes that label each atom while the indexes in COVARMAT span over al the coordenates of each atom.
!....................................................................................................!
       integer, intent(in)   ::  i,j
       real, intent(in)      ::  COVARMAT(:,:)
       real                  ::  determinant
       real                  ::  temp
       temp=COVARMAT((i-1)*3+1,(j-1)*3+1)
       temp=temp*COVARMAT((i-1)*3+2,(j-1)*3+2)
       temp=temp*COVARMAT((i-1)*3+3,(j-1)*3+3)
       determinant=temp
       temp=COVARMAT((i-1)*3+1,(j-1)*3+2)
       temp=temp*COVARMAT((i-1)*3+2,(j-1)*3+3)
       temp=temp*COVARMAT((i-1)*3+3,(j-1)*3+1)
       determinant=determinant+temp
       temp=COVARMAT((i-1)*3+1,(j-1)*3+3)
       temp=temp*COVARMAT((i-1)*3+2,(j-1)*3+1)
       temp=temp*COVARMAT((i-1)*3+3,(j-1)*3+2)
       determinant=determinant+temp
       temp=COVARMAT((i-1)*3+1,(j-1)*3+3)
       temp=temp*COVARMAT((i-1)*3+2,(j-1)*3+2)
       temp=temp*COVARMAT((i-1)*3+3,(j-1)*3+1)
       determinant=determinant-temp
       temp=COVARMAT((i-1)*3+1,(j-1)*3+2)
       temp=temp*COVARMAT((i-1)*3+2,(j-1)*3+1)
       temp=temp*COVARMAT((i-1)*3+3,(j-1)*3+3)
       determinant=determinant-temp
       temp=COVARMAT((i-1)*3+1,(j-1)*3+1)
       temp=temp*COVARMAT((i-1)*3+2,(j-1)*3+3)
       temp=temp*COVARMAT((i-1)*3+3,(j-1)*3+2)
       determinant=determinant-temp
      
      end function DET
!------------------------------------------------------------!
    REAL function DETv1(COVARMAT, ii,jj,n)
    IMPLICIT NONE
!.....................................................................................................!
! ONLY WORKS WHEN n=6 or n=3 (i.e. 6x6 or 3x3 matrices) !!!!!!!!!!!!!
! 3x3 or 6x6 determinant
! The subroutine is based on two key points:
! [1] A determinant is unaltered when row operations are performed: Hence, using this principle,
!     row operations (column operations would work as well) are used
!     to convert the matrix into upper traingular form
! [2] The determinant of a triangular matrix is obtained by finding the product of the diagonal elements
! Note: In order to make this function more efficient, it would be better to get rid of the variable matrix or vectorize the loops.
!.....................................................................................................!
    INTEGER, INTENT(IN)  :: n,ii,jj
    REAL, INTENT(IN)     :: COVARMAT(:,:)
    REAL, DIMENSION(n,n) :: matrix
    REAL                 :: m, temp
    INTEGER              :: i, j, k, l
    LOGICAL              :: DetExists = .TRUE.
    l = 1
    if(n.eq.6) then
       do i=1,3
          do j=1,3 
             matrix(i,j)=COVARMAT((ii-1)*3+i,(ii-1)*3+j)
             matrix(i,j+3)=COVARMAT((ii-1)*3+i,(jj-1)*3+j)
             matrix(i+3,j)=COVARMAT((jj-1)*3+i,(ii-1)*3+j)
             matrix(i+3,j+3)=COVARMAT((jj-1)*3+i,(jj-1)*3+j)
          enddo
       enddo
    else if(n.eq.3) then
       do i=1,3
          do j=1,3
             matrix(i,j)=COVARMAT((ii-1)*3+i,(jj-1)*3+j)
          enddo
       enddo
    else 
       stop 'wrong determinant size, the current version of DETv1 only allows determinant of 6x6 or 3x3 matrices'
    endif
    
    !Convert to upper triangular form
    do k = 1, n-1
        if (matrix(k,k) == 0) then
            DetExists = .FALSE.
            do i = k+1, n
                if (matrix(i,k) /= 0) then
                    do j = 1, n
                        temp = matrix(i,j)
                        matrix(i,j)= matrix(k,j)
                        matrix(k,j) = temp
                    enddo
                    DetExists = .TRUE.
                    l=-l
                    exit
                endif
            end do
            if (DetExists .EQV. .FALSE.) then
               DETv1=0.0000
               write(15,*) n,ii,jj
               return
            endif
        endif
        do j = k+1, n
            m = matrix(j,k)/matrix(k,k)
            do i = k+1, n
                matrix(j,i) = matrix(j,i) - m*matrix(k,i)
            enddo
        enddo
    enddo
    
    !Calculate determinant by finding product of diagonal elements
    DETv1 = l
    do i = 1, n
        DETv1 = DETv1 * matrix(i,i)
    enddo
    
end function DETv1

    REAL function DETv2(COVAR,n)
    IMPLICIT NONE
!.....................................................................................................!
! ONLY WORKS WHEN n=6 or n=3 (i.e. 6x6 or 3x3 matrices) !!!!!!!!!!!!!
! 3x3 or 6x6 determinant
! The subroutine is based on two key points:
! [1] A determinant is unaltered when row operations are performed: Hence, using this principle,
!     row operations (column operations would work as well) are used
!     to convert the matrix into upper traingular form
! [2] The determinant of a triangular matrix is obtained by finding the product of the diagonal elements
! Note: In order to make this function more efficient, it would be better to get rid of the variable matrix or vectorize the loops.
!.....................................................................................................!
    INTEGER, INTENT(IN)  :: n
    REAL,intent(in)      :: COVAR(:,:)
    REAL,allocatable     :: matrix(:,:)
    REAL                 :: m, temp
    INTEGER              :: i, j, k, l
    LOGICAL              :: DetExists = .TRUE.
    l = 1
    if(n.eq.6) then
       allocate(matrix(n,n))
       matrix=covar
    elseif(n.eq.3) then
       allocate(matrix(n,n))
       do i=1,3
          do j=1,3
             matrix(i,j)=COVAR(i+3,j+3)
          enddo
       enddo
    else 
       stop 'wrong determnant dimension'
    endif
    !Convert to upper triangular form
    do k = 1, n-1
        if (matrix(k,k) == 0) then
            DetExists = .FALSE.
            do i = k+1, n
                if (matrix(i,k) /= 0) then
                    do j = 1, n
                        temp = matrix(i,j)
                        matrix(i,j)= matrix(k,j)
                        matrix(k,j) = temp
                    enddo
                    DetExists = .TRUE.
                    l=-l
                    exit
                endif
            end do
            if (DetExists .EQV. .FALSE.) then
               DETv2=0.0000
               stop 'Determinants could not be calculated'
            endif
        endif
        do j = k+1, n
            m = matrix(j,k)/matrix(k,k)
            do i = k+1, n
                matrix(j,i) = matrix(j,i) - m*matrix(k,i)
            enddo
        enddo
    enddo

    !Calculate determinant by finding product of diagonal elements
    DETv2 = l
    do i = 1, n
        DETv2 = DETv2 * matrix(i,i)
    enddo

    end function DETv2

!%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%!
    end program
