! Fortran file contains all subroutines needed for FE-calculations
! with the C3D4C1 finite element as presented at 
! https://github.com/Ra-Na/C3D4C1-Abaqus-UEL 
!
! Structure:
! ==========
!
!    UEL      <- Abaqus entry point
!    | A         Performs index-rearrangement from N x N matrices 
!    | |         to higher-index objects 
!    V |
!    UELR     <- Performs actually the element calculations, calls 
!    | A         auxiliary subroutines and UMATR
!    | |         
!    V |         
!    UMATR    <- Defines the mechanical material behavior
!                in terms of first Piola-Kirchhoff-stresses
!                and first Piola-Kirchhoff-like hyperstresses 
!                as a function of the first and second displacement
!                gradient
!    

      SUBROUTINE UEL(RHS,AMATRX,SVARS,ENERGY,NDOFEL,NRHS,NSVARS,
     1     PROPS,NPROPS,COORDS,MCRD,NNODE,U,DU,V,A,JTYPE,TIME,DTIME,
     2     KSTEP,KINC,JELEM,PARAMS,NDLOAD,JDLTYP,ADLMAG,PREDEF,
     3     NPREDF,LFLAGS,MLVARX,DDLMAG,MDLOAD,PNEWDT,JPROPS,NJPROP,
     4     PERIOD)
      IMPLICIT NONE
! Abaqus defs: i,j,k,l,m,n INTEGER, all other DOUBLE PREC
      INTEGER MLVARX,NDOFEL,MDLOAD,MCRD,NNODE,NSVARS,NRHS,NDLOAD
      INTEGER NPROPS,JTYPE,JELEM,KSTEP,KINC,NJPROP,NPREDF
      INTEGER LFLAGS(*),JPROPS(*),JDLTYP(MDLOAD,*)
      DOUBLE PRECISION RHS(MLVARX,*),AMATRX(NDOFEL,NDOFEL),PROPS(*),
     1 SVARS(*),ENERGY(8),COORDS(MCRD,NNODE),U(NDOFEL),
     2 DU(MLVARX,*),V(NDOFEL),A(NDOFEL),TIME(2),PARAMS(*),
     3 ADLMAG(MDLOAD,*),DDLMAG(MDLOAD,*),DTIME,PNEWDT,PERIOD,
     4 PREDEF(2,NPREDF,NNODE)
! Own definitions
      INTEGER nooffields,ndofelemental,IDUMMY,svatpg,unodes
      INTEGER i,j,k,l,m,n
      
      parameter (unodes=4)         ! unique nodes
      parameter (nooffields=3)     ! no. of C1-cont fields
      parameter (ndofelemental=4)  ! no. of elemental DOF per node and field
      
      double precision koordinaten(3,unodes)                             ! 3 x 4
      double precision dof(nooffields,unodes,ndofelemental)              ! 3 x 4 x 4
      double precision force(nooffields,unodes,ndofelemental)            ! 3 x 4 x 4
      double precision dforce_ddof(nooffields,unodes,ndofelemental,      ! 3 x 4 x 4  x  3 x 4 x 4
     &                             nooffields,unodes,ndofelemental)
      integer indexmatrix(3,4)
      integer indexvektor1(12),indexvektor2(12),echo
      
      double precision H1(3,3),H2(3,3),H3(3,3),H4(3,3)
      double precision u1(3),u2(3),u3(3),u4(3),piofo(10)
      double precision rdummy,rd2,rd3,rd4,uav(3),cav(3)
      
      character(len=200) :: filename1,filename2,file1,file2,nostr,
     &                      filename3,file3,dirname,filename4,file4,
     &                      filename5,file5,filename6,file6
      
      
      INTEGER*4 getcwd,status
            
      logical :: exist


!     Set 1 or 0 to turn output on or of
      echo=1
!     Create output files: when jelem=1 a new iteration starts, files
!     get overwritten, the last one remains for examination (last 
!     increment, last iteration)
      if(echo.eq.1)then
!       filename 1: integration point values
!       filename 3: unique nodal values
        dirname="~/data/ecke/sehrfein/a0/"
        filename1="IPVALS20_"
        filename3="UNVALS20_"
        write(nostr,'(I2.2)') kinc
        nostr="0"
        file1=trim(dirname)//trim(filename1)//trim(nostr)//".txt"
        file3=trim(dirname)//trim(filename3)//trim(nostr)//".txt"
        if(jelem.eq.1) then
          inquire(file=file1, exist=exist)
          if (exist) then
            close(28)
            open (28, file=file1, form='unformatted')
            close (28, status='delete') 
          end if
          open(28, file=file1, status="new", action="write")

          inquire(file=file3, exist=exist)
          if (exist) then
            close(30)
            open (30, file=file3, form='unformatted')
            close (30, status='delete') 
          end if
          open(30, file=file3, status="new", action="write")
        end if  
      end if
      
      indexmatrix(1,1)=1
      indexmatrix(2,1)=2
      indexmatrix(3,1)=3
      indexmatrix(1,2)=4
      indexmatrix(1,3)=5
      indexmatrix(1,4)=6
      indexmatrix(2,2)=7
      indexmatrix(2,3)=8
      indexmatrix(2,4)=9
      indexmatrix(3,2)=10
      indexmatrix(3,3)=11
      indexmatrix(3,4)=12
      
!     rearrange indexing
      do i=1,unodes
      do j=1,3
        koordinaten(j,i)=COORDS(j,i)
      end do
      end do

      do i=1,3
      do j=1,4
      do k=1,4
        dof(i,j,k)=u(12*(j-1)+indexmatrix(i,k))
      end do
      end do
      end do

!     extract displacement and displacement gradient
!     at each corner node
      u1(1)=dof(1,1,1)
      u1(2)=dof(2,1,1)
      u1(3)=dof(3,1,1)
      
      u2(1)=dof(1,2,1)
      u2(2)=dof(2,2,1)
      u2(3)=dof(3,2,1)

      u3(1)=dof(1,3,1)
      u3(2)=dof(2,3,1)
      u3(3)=dof(3,3,1)

      u4(1)=dof(1,4,1)
      u4(2)=dof(2,4,1)
      u4(3)=dof(3,4,1)

      H1(1,1)=dof(1,1,2)
      H1(1,2)=dof(1,1,3)
      H1(1,3)=dof(1,1,4)
      H1(2,1)=dof(2,1,2)
      H1(2,2)=dof(2,1,3)
      H1(2,3)=dof(2,1,4)
      H1(3,1)=dof(3,1,2)
      H1(3,2)=dof(3,1,3)
      H1(3,3)=dof(3,1,4)

      H2(1,1)=dof(1,2,2)
      H2(1,2)=dof(1,2,3)
      H2(1,3)=dof(1,2,4)
      H2(2,1)=dof(2,2,2)
      H2(2,2)=dof(2,2,3)
      H2(2,3)=dof(2,2,4)
      H2(3,1)=dof(3,2,2)
      H2(3,2)=dof(3,2,3)
      H2(3,3)=dof(3,2,4)
      
      H3(1,1)=dof(1,3,2)
      H3(1,2)=dof(1,3,3)
      H3(1,3)=dof(1,3,4)
      H3(2,1)=dof(2,3,2)
      H3(2,2)=dof(2,3,3)
      H3(2,3)=dof(2,3,4)
      H3(3,1)=dof(3,3,2)
      H3(3,2)=dof(3,3,3)
      H3(3,3)=dof(3,3,4)
      
      H4(1,1)=dof(1,4,2)
      H4(1,2)=dof(1,4,3)
      H4(1,3)=dof(1,4,4)
      H4(2,1)=dof(2,4,2)
      H4(2,2)=dof(2,4,3)
      H4(2,3)=dof(2,4,4)
      H4(3,1)=dof(3,4,2)
      H4(3,2)=dof(3,4,3)
      H4(3,3)=dof(3,4,4)

!     write some output
      if(echo.eq.1) then 
      write(28,*) jelem, 29, 8
      write(30,*) jelem, 4 , 3
      write(30,'(12f20.10)') u1,H1
      write(30,'(12f20.10)') u2,H2
      write(30,'(12f20.10)') u3,H3
      write(30,'(12f20.10)') u4,H4
      end if
      
!     call actual user element routine
      call uelr(koordinaten,dof,force,dforce_ddof,echo,piofo)

!     rearrange force and linearization for the Abaqus-indexing
     
      do i=1,3
      do j=1,4
      do k=1,4
        RHS(12*(j-1)+indexmatrix(i,k),1)=-force(i,j,k)
      end do
      end do
      end do
      
      do i=1,3
      do j=1,4
      do k=1,4
      do l=1,3
      do m=1,4
      do n=1,4
        AMATRX(12*(j-1)+indexmatrix(i,k),12*(m-1)+indexmatrix(l,n))=
     &         dforce_ddof(i,j,k,l,m,n)
      end do
      end do
      end do
      end do
      end do
      end do
      return
      end


!
!
!
!
!
!
!
!
!
!
!
!


      subroutine uelr(coords,dof,force,dforce_ddof,echo,piofo)
! Input:      reference coordinates       : coords
!             DOF-Matrix                  : dof
! Output:     to DOF conjugate forces     : force
!             linearization               : dforce_ddof 
      implicit none
      double precision coords(3,4),dof(3,4,4),force(3,4,4)
      double precision dforce_ddof(3,4,4,3,4,4)

      double precision matrixi(20,20)
      double precision q(20,20)
      double precision matrix2(20,16)
      double precision matrix3(20,16)
      

      double precision x1,x2,x3,x4,y1,y2,y3,y4,z1,z2,z3,z4
      double precision phi1,phi1x,phi1y,phi1z
      double precision phi2,phi2x,phi2y,phi2z
      double precision phi3,phi3x,phi3y,phi3z
      double precision phi4,phi4x,phi4y,phi4z
      
      integer LDA,NPIVOT,IERR,MX(20),MY(20),i,j,k,l,m,n,gpindex,i2
      integer s1,s2,s3,s4,s5,s6
      
      double precision x,y,z,B0(4,4),B1(4,4,3),B2(4,4,3,3),rdummy
      double precision u(3),h1(3,3),h2(3,3,3),detf,ps1,ps2,rd2,val
      
      double precision GP(50,4),F(3,3),GPakt(50,4),det33,weight
      integer ngp,echo
      
      double precision T(3,3),T3(3,3,3),pass(3)
      double precision DTDH(3,3,3,3)
      double precision DTDH3(3,3,3,3,3)
      double precision DT3DH(3,3,3,3,3)
      double precision DT3DH3(3,3,3,3,3,3)
      
      double precision norm333,piofo(10)
      
      x1=coords(1,1)
      y1=coords(2,1)
      z1=coords(3,1)

      x2=coords(1,2)
      y2=coords(2,2)
      z2=coords(3,2)

      x3=coords(1,3)
      y3=coords(2,3)
      z3=coords(3,3)

      x4=coords(1,4)
      y4=coords(2,4)
      z4=coords(3,4)

! matrix2 maps from the 16 dof to the 20 pseudo-dof 
      matrix2(1,1)=1
      matrix2(1,2)=0
      matrix2(1,3)=0
      matrix2(1,4)=0
      matrix2(1,5)=0
      matrix2(1,6)=0
      matrix2(1,7)=0
      matrix2(1,8)=0
      matrix2(1,9)=0
      matrix2(1,10)=0
      matrix2(1,11)=0
      matrix2(1,12)=0
      matrix2(1,13)=0
      matrix2(1,14)=0
      matrix2(1,15)=0
      matrix2(1,16)=0
      matrix2(2,1)=0
      matrix2(2,2)=1
      matrix2(2,3)=0
      matrix2(2,4)=0
      matrix2(2,5)=0
      matrix2(2,6)=0
      matrix2(2,7)=0
      matrix2(2,8)=0
      matrix2(2,9)=0
      matrix2(2,10)=0
      matrix2(2,11)=0
      matrix2(2,12)=0
      matrix2(2,13)=0
      matrix2(2,14)=0
      matrix2(2,15)=0
      matrix2(2,16)=0
      matrix2(3,1)=0
      matrix2(3,2)=0
      matrix2(3,3)=1
      matrix2(3,4)=0
      matrix2(3,5)=0
      matrix2(3,6)=0
      matrix2(3,7)=0
      matrix2(3,8)=0
      matrix2(3,9)=0
      matrix2(3,10)=0
      matrix2(3,11)=0
      matrix2(3,12)=0
      matrix2(3,13)=0
      matrix2(3,14)=0
      matrix2(3,15)=0
      matrix2(3,16)=0
      matrix2(4,1)=0
      matrix2(4,2)=0
      matrix2(4,3)=0
      matrix2(4,4)=1
      matrix2(4,5)=0
      matrix2(4,6)=0
      matrix2(4,7)=0
      matrix2(4,8)=0
      matrix2(4,9)=0
      matrix2(4,10)=0
      matrix2(4,11)=0
      matrix2(4,12)=0
      matrix2(4,13)=0
      matrix2(4,14)=0
      matrix2(4,15)=0
      matrix2(4,16)=0
      matrix2(5,1)=0
      matrix2(5,2)=0
      matrix2(5,3)=0
      matrix2(5,4)=0
      matrix2(5,5)=1
      matrix2(5,6)=0
      matrix2(5,7)=0
      matrix2(5,8)=0
      matrix2(5,9)=0
      matrix2(5,10)=0
      matrix2(5,11)=0
      matrix2(5,12)=0
      matrix2(5,13)=0
      matrix2(5,14)=0
      matrix2(5,15)=0
      matrix2(5,16)=0
      matrix2(6,1)=0
      matrix2(6,2)=0
      matrix2(6,3)=0
      matrix2(6,4)=0
      matrix2(6,5)=0
      matrix2(6,6)=1
      matrix2(6,7)=0
      matrix2(6,8)=0
      matrix2(6,9)=0
      matrix2(6,10)=0
      matrix2(6,11)=0
      matrix2(6,12)=0
      matrix2(6,13)=0
      matrix2(6,14)=0
      matrix2(6,15)=0
      matrix2(6,16)=0
      matrix2(7,1)=0
      matrix2(7,2)=0
      matrix2(7,3)=0
      matrix2(7,4)=0
      matrix2(7,5)=0
      matrix2(7,6)=0
      matrix2(7,7)=1
      matrix2(7,8)=0
      matrix2(7,9)=0
      matrix2(7,10)=0
      matrix2(7,11)=0
      matrix2(7,12)=0
      matrix2(7,13)=0
      matrix2(7,14)=0
      matrix2(7,15)=0
      matrix2(7,16)=0
      matrix2(8,1)=0
      matrix2(8,2)=0
      matrix2(8,3)=0
      matrix2(8,4)=0
      matrix2(8,5)=0
      matrix2(8,6)=0
      matrix2(8,7)=0
      matrix2(8,8)=1
      matrix2(8,9)=0
      matrix2(8,10)=0
      matrix2(8,11)=0
      matrix2(8,12)=0
      matrix2(8,13)=0
      matrix2(8,14)=0
      matrix2(8,15)=0
      matrix2(8,16)=0
      matrix2(9,1)=0
      matrix2(9,2)=0
      matrix2(9,3)=0
      matrix2(9,4)=0
      matrix2(9,5)=0
      matrix2(9,6)=0
      matrix2(9,7)=0
      matrix2(9,8)=0
      matrix2(9,9)=1
      matrix2(9,10)=0
      matrix2(9,11)=0
      matrix2(9,12)=0
      matrix2(9,13)=0
      matrix2(9,14)=0
      matrix2(9,15)=0
      matrix2(9,16)=0
      matrix2(10,1)=0
      matrix2(10,2)=0
      matrix2(10,3)=0
      matrix2(10,4)=0
      matrix2(10,5)=0
      matrix2(10,6)=0
      matrix2(10,7)=0
      matrix2(10,8)=0
      matrix2(10,9)=0
      matrix2(10,10)=1
      matrix2(10,11)=0
      matrix2(10,12)=0
      matrix2(10,13)=0
      matrix2(10,14)=0
      matrix2(10,15)=0
      matrix2(10,16)=0
      matrix2(11,1)=0
      matrix2(11,2)=0
      matrix2(11,3)=0
      matrix2(11,4)=0
      matrix2(11,5)=0
      matrix2(11,6)=0
      matrix2(11,7)=0
      matrix2(11,8)=0
      matrix2(11,9)=0
      matrix2(11,10)=0
      matrix2(11,11)=1
      matrix2(11,12)=0
      matrix2(11,13)=0
      matrix2(11,14)=0
      matrix2(11,15)=0
      matrix2(11,16)=0
      matrix2(12,1)=0
      matrix2(12,2)=0
      matrix2(12,3)=0
      matrix2(12,4)=0
      matrix2(12,5)=0
      matrix2(12,6)=0
      matrix2(12,7)=0
      matrix2(12,8)=0
      matrix2(12,9)=0
      matrix2(12,10)=0
      matrix2(12,11)=0
      matrix2(12,12)=1
      matrix2(12,13)=0
      matrix2(12,14)=0
      matrix2(12,15)=0
      matrix2(12,16)=0
      matrix2(13,1)=0
      matrix2(13,2)=0
      matrix2(13,3)=0
      matrix2(13,4)=0
      matrix2(13,5)=0
      matrix2(13,6)=0
      matrix2(13,7)=0
      matrix2(13,8)=0
      matrix2(13,9)=0
      matrix2(13,10)=0
      matrix2(13,11)=0
      matrix2(13,12)=0
      matrix2(13,13)=1
      matrix2(13,14)=0
      matrix2(13,15)=0
      matrix2(13,16)=0
      matrix2(14,1)=0
      matrix2(14,2)=0
      matrix2(14,3)=0
      matrix2(14,4)=0
      matrix2(14,5)=0
      matrix2(14,6)=0
      matrix2(14,7)=0
      matrix2(14,8)=0
      matrix2(14,9)=0
      matrix2(14,10)=0
      matrix2(14,11)=0
      matrix2(14,12)=0
      matrix2(14,13)=0
      matrix2(14,14)=1
      matrix2(14,15)=0
      matrix2(14,16)=0
      matrix2(15,1)=0
      matrix2(15,2)=0
      matrix2(15,3)=0
      matrix2(15,4)=0
      matrix2(15,5)=0
      matrix2(15,6)=0
      matrix2(15,7)=0
      matrix2(15,8)=0
      matrix2(15,9)=0
      matrix2(15,10)=0
      matrix2(15,11)=0
      matrix2(15,12)=0
      matrix2(15,13)=0
      matrix2(15,14)=0
      matrix2(15,15)=1
      matrix2(15,16)=0
      matrix2(16,1)=0
      matrix2(16,2)=0
      matrix2(16,3)=0
      matrix2(16,4)=0
      matrix2(16,5)=0
      matrix2(16,6)=0
      matrix2(16,7)=0
      matrix2(16,8)=0
      matrix2(16,9)=0
      matrix2(16,10)=0
      matrix2(16,11)=0
      matrix2(16,12)=0
      matrix2(16,13)=0
      matrix2(16,14)=0
      matrix2(16,15)=0
      matrix2(16,16)=1
      matrix2(17,1)=0.3333333333333333
      matrix2(17,2)=(-2*x1 + x2 + x4)/18.
      matrix2(17,3)=(-2*y1 + y2 + y4)/18.
      matrix2(17,4)=(-2*z1 + z2 + z4)/18.
      matrix2(17,5)=0.3333333333333333
      matrix2(17,6)=(x1 - 2*x2 + x4)/18.
      matrix2(17,7)=(y1 - 2*y2 + y4)/18.
      matrix2(17,8)=(z1 - 2*z2 + z4)/18.
      matrix2(17,9)=0
      matrix2(17,10)=0
      matrix2(17,11)=0
      matrix2(17,12)=0
      matrix2(17,13)=0.3333333333333333
      matrix2(17,14)=(x1 + x2 - 2*x4)/18.
      matrix2(17,15)=(y1 + y2 - 2*y4)/18.
      matrix2(17,16)=(z1 + z2 - 2*z4)/18.
      matrix2(18,1)=0.3333333333333333
      matrix2(18,2)=(-2*x1 + x2 + x3)/18.
      matrix2(18,3)=(-2*y1 + y2 + y3)/18.
      matrix2(18,4)=(-2*z1 + z2 + z3)/18.
      matrix2(18,5)=0.3333333333333333
      matrix2(18,6)=(x1 - 2*x2 + x3)/18.
      matrix2(18,7)=(y1 - 2*y2 + y3)/18.
      matrix2(18,8)=(z1 - 2*z2 + z3)/18.
      matrix2(18,9)=0.3333333333333333
      matrix2(18,10)=(x1 + x2 - 2*x3)/18.
      matrix2(18,11)=(y1 + y2 - 2*y3)/18.
      matrix2(18,12)=(z1 + z2 - 2*z3)/18.
      matrix2(18,13)=0
      matrix2(18,14)=0
      matrix2(18,15)=0
      matrix2(18,16)=0
      matrix2(19,1)=0.3333333333333333
      matrix2(19,2)=(-2*x1 + x3 + x4)/18.
      matrix2(19,3)=(-2*y1 + y3 + y4)/18.
      matrix2(19,4)=(-2*z1 + z3 + z4)/18.
      matrix2(19,5)=0
      matrix2(19,6)=0
      matrix2(19,7)=0
      matrix2(19,8)=0
      matrix2(19,9)=0.3333333333333333
      matrix2(19,10)=(x1 - 2*x3 + x4)/18.
      matrix2(19,11)=(y1 - 2*y3 + y4)/18.
      matrix2(19,12)=(z1 - 2*z3 + z4)/18.
      matrix2(19,13)=0.3333333333333333
      matrix2(19,14)=(x1 + x3 - 2*x4)/18.
      matrix2(19,15)=(y1 + y3 - 2*y4)/18.
      matrix2(19,16)=(z1 + z3 - 2*z4)/18.
      matrix2(20,1)=0
      matrix2(20,2)=0
      matrix2(20,3)=0
      matrix2(20,4)=0
      matrix2(20,5)=0.3333333333333333
      matrix2(20,6)=(-2*x2 + x3 + x4)/18.
      matrix2(20,7)=(-2*y2 + y3 + y4)/18.
      matrix2(20,8)=(-2*z2 + z3 + z4)/18.
      matrix2(20,9)=0.3333333333333333
      matrix2(20,10)=(x2 - 2*x3 + x4)/18.
      matrix2(20,11)=(y2 - 2*y3 + y4)/18.
      matrix2(20,12)=(z2 - 2*z3 + z4)/18.
      matrix2(20,13)=0.3333333333333333
      matrix2(20,14)=(x2 + x3 - 2*x4)/18.
      matrix2(20,15)=(y2 + y3 - 2*y4)/18.
      matrix2(20,16)=(z2 + z3 - 2*z4)/18.
      
! the q-matrix maps from the 20 polynomial coeffs. to the 20 pseudo-dof
      q(1,1)=1
      q(1,2)=z1
      q(1,3)=z1**2
      q(1,4)=z1**3
      q(1,5)=y1
      q(1,6)=y1*z1
      q(1,7)=y1*z1**2
      q(1,8)=y1**2
      q(1,9)=y1**2*z1
      q(1,10)=y1**3
      q(1,11)=x1
      q(1,12)=x1*z1
      q(1,13)=x1*z1**2
      q(1,14)=x1*y1
      q(1,15)=x1*y1*z1
      q(1,16)=x1*y1**2
      q(1,17)=x1**2
      q(1,18)=x1**2*z1
      q(1,19)=x1**2*y1
      q(1,20)=x1**3
      q(2,1)=0
      q(2,2)=0
      q(2,3)=0
      q(2,4)=0
      q(2,5)=0
      q(2,6)=0
      q(2,7)=0
      q(2,8)=0
      q(2,9)=0
      q(2,10)=0
      q(2,11)=1
      q(2,12)=z1
      q(2,13)=z1**2
      q(2,14)=y1
      q(2,15)=y1*z1
      q(2,16)=y1**2
      q(2,17)=2*x1
      q(2,18)=2*x1*z1
      q(2,19)=2*x1*y1
      q(2,20)=3*x1**2
      q(3,1)=0
      q(3,2)=0
      q(3,3)=0
      q(3,4)=0
      q(3,5)=1
      q(3,6)=z1
      q(3,7)=z1**2
      q(3,8)=2*y1
      q(3,9)=2*y1*z1
      q(3,10)=3*y1**2
      q(3,11)=0
      q(3,12)=0
      q(3,13)=0
      q(3,14)=x1
      q(3,15)=x1*z1
      q(3,16)=2*x1*y1
      q(3,17)=0
      q(3,18)=0
      q(3,19)=x1**2
      q(3,20)=0
      q(4,1)=0
      q(4,2)=1
      q(4,3)=2*z1
      q(4,4)=3*z1**2
      q(4,5)=0
      q(4,6)=y1
      q(4,7)=2*y1*z1
      q(4,8)=0
      q(4,9)=y1**2
      q(4,10)=0
      q(4,11)=0
      q(4,12)=x1
      q(4,13)=2*x1*z1
      q(4,14)=0
      q(4,15)=x1*y1
      q(4,16)=0
      q(4,17)=0
      q(4,18)=x1**2
      q(4,19)=0
      q(4,20)=0
      q(5,1)=1
      q(5,2)=z2
      q(5,3)=z2**2
      q(5,4)=z2**3
      q(5,5)=y2
      q(5,6)=y2*z2
      q(5,7)=y2*z2**2
      q(5,8)=y2**2
      q(5,9)=y2**2*z2
      q(5,10)=y2**3
      q(5,11)=x2
      q(5,12)=x2*z2
      q(5,13)=x2*z2**2
      q(5,14)=x2*y2
      q(5,15)=x2*y2*z2
      q(5,16)=x2*y2**2
      q(5,17)=x2**2
      q(5,18)=x2**2*z2
      q(5,19)=x2**2*y2
      q(5,20)=x2**3
      q(6,1)=0
      q(6,2)=0
      q(6,3)=0
      q(6,4)=0
      q(6,5)=0
      q(6,6)=0
      q(6,7)=0
      q(6,8)=0
      q(6,9)=0
      q(6,10)=0
      q(6,11)=1
      q(6,12)=z2
      q(6,13)=z2**2
      q(6,14)=y2
      q(6,15)=y2*z2
      q(6,16)=y2**2
      q(6,17)=2*x2
      q(6,18)=2*x2*z2
      q(6,19)=2*x2*y2
      q(6,20)=3*x2**2
      q(7,1)=0
      q(7,2)=0
      q(7,3)=0
      q(7,4)=0
      q(7,5)=1
      q(7,6)=z2
      q(7,7)=z2**2
      q(7,8)=2*y2
      q(7,9)=2*y2*z2
      q(7,10)=3*y2**2
      q(7,11)=0
      q(7,12)=0
      q(7,13)=0
      q(7,14)=x2
      q(7,15)=x2*z2
      q(7,16)=2*x2*y2
      q(7,17)=0
      q(7,18)=0
      q(7,19)=x2**2
      q(7,20)=0
      q(8,1)=0
      q(8,2)=1
      q(8,3)=2*z2
      q(8,4)=3*z2**2
      q(8,5)=0
      q(8,6)=y2
      q(8,7)=2*y2*z2
      q(8,8)=0
      q(8,9)=y2**2
      q(8,10)=0
      q(8,11)=0
      q(8,12)=x2
      q(8,13)=2*x2*z2
      q(8,14)=0
      q(8,15)=x2*y2
      q(8,16)=0
      q(8,17)=0
      q(8,18)=x2**2
      q(8,19)=0
      q(8,20)=0
      q(9,1)=1
      q(9,2)=z3
      q(9,3)=z3**2
      q(9,4)=z3**3
      q(9,5)=y3
      q(9,6)=y3*z3
      q(9,7)=y3*z3**2
      q(9,8)=y3**2
      q(9,9)=y3**2*z3
      q(9,10)=y3**3
      q(9,11)=x3
      q(9,12)=x3*z3
      q(9,13)=x3*z3**2
      q(9,14)=x3*y3
      q(9,15)=x3*y3*z3
      q(9,16)=x3*y3**2
      q(9,17)=x3**2
      q(9,18)=x3**2*z3
      q(9,19)=x3**2*y3
      q(9,20)=x3**3
      q(10,1)=0
      q(10,2)=0
      q(10,3)=0
      q(10,4)=0
      q(10,5)=0
      q(10,6)=0
      q(10,7)=0
      q(10,8)=0
      q(10,9)=0
      q(10,10)=0
      q(10,11)=1
      q(10,12)=z3
      q(10,13)=z3**2
      q(10,14)=y3
      q(10,15)=y3*z3
      q(10,16)=y3**2
      q(10,17)=2*x3
      q(10,18)=2*x3*z3
      q(10,19)=2*x3*y3
      q(10,20)=3*x3**2
      q(11,1)=0
      q(11,2)=0
      q(11,3)=0
      q(11,4)=0
      q(11,5)=1
      q(11,6)=z3
      q(11,7)=z3**2
      q(11,8)=2*y3
      q(11,9)=2*y3*z3
      q(11,10)=3*y3**2
      q(11,11)=0
      q(11,12)=0
      q(11,13)=0
      q(11,14)=x3
      q(11,15)=x3*z3
      q(11,16)=2*x3*y3
      q(11,17)=0
      q(11,18)=0
      q(11,19)=x3**2
      q(11,20)=0
      q(12,1)=0
      q(12,2)=1
      q(12,3)=2*z3
      q(12,4)=3*z3**2
      q(12,5)=0
      q(12,6)=y3
      q(12,7)=2*y3*z3
      q(12,8)=0
      q(12,9)=y3**2
      q(12,10)=0
      q(12,11)=0
      q(12,12)=x3
      q(12,13)=2*x3*z3
      q(12,14)=0
      q(12,15)=x3*y3
      q(12,16)=0
      q(12,17)=0
      q(12,18)=x3**2
      q(12,19)=0
      q(12,20)=0
      q(13,1)=1
      q(13,2)=z4
      q(13,3)=z4**2
      q(13,4)=z4**3
      q(13,5)=y4
      q(13,6)=y4*z4
      q(13,7)=y4*z4**2
      q(13,8)=y4**2
      q(13,9)=y4**2*z4
      q(13,10)=y4**3
      q(13,11)=x4
      q(13,12)=x4*z4
      q(13,13)=x4*z4**2
      q(13,14)=x4*y4
      q(13,15)=x4*y4*z4
      q(13,16)=x4*y4**2
      q(13,17)=x4**2
      q(13,18)=x4**2*z4
      q(13,19)=x4**2*y4
      q(13,20)=x4**3
      q(14,1)=0
      q(14,2)=0
      q(14,3)=0
      q(14,4)=0
      q(14,5)=0
      q(14,6)=0
      q(14,7)=0
      q(14,8)=0
      q(14,9)=0
      q(14,10)=0
      q(14,11)=1
      q(14,12)=z4
      q(14,13)=z4**2
      q(14,14)=y4
      q(14,15)=y4*z4
      q(14,16)=y4**2
      q(14,17)=2*x4
      q(14,18)=2*x4*z4
      q(14,19)=2*x4*y4
      q(14,20)=3*x4**2
      q(15,1)=0
      q(15,2)=0
      q(15,3)=0
      q(15,4)=0
      q(15,5)=1
      q(15,6)=z4
      q(15,7)=z4**2
      q(15,8)=2*y4
      q(15,9)=2*y4*z4
      q(15,10)=3*y4**2
      q(15,11)=0
      q(15,12)=0
      q(15,13)=0
      q(15,14)=x4
      q(15,15)=x4*z4
      q(15,16)=2*x4*y4
      q(15,17)=0
      q(15,18)=0
      q(15,19)=x4**2
      q(15,20)=0
      q(16,1)=0
      q(16,2)=1
      q(16,3)=2*z4
      q(16,4)=3*z4**2
      q(16,5)=0
      q(16,6)=y4
      q(16,7)=2*y4*z4
      q(16,8)=0
      q(16,9)=y4**2
      q(16,10)=0
      q(16,11)=0
      q(16,12)=x4
      q(16,13)=2*x4*z4
      q(16,14)=0
      q(16,15)=x4*y4
      q(16,16)=0
      q(16,17)=0
      q(16,18)=x4**2
      q(16,19)=0
      q(16,20)=0
      q(17,1)=1
      q(17,2)=(z1 + z2 + z4)/3.
      q(17,3)=(z1 + z2 + z4)**2/9.
      q(17,4)=(z1 + z2 + z4)**3/27.
      q(17,5)=(y1 + y2 + y4)/3.
      q(17,6)=((y1 + y2 + y4)*(z1 + z2 + z4))/9.
      q(17,7)=((y1 + y2 + y4)*(z1 + z2 + z4)**2)/27.
      q(17,8)=(y1 + y2 + y4)**2/9.
      q(17,9)=((y1 + y2 + y4)**2*(z1 + z2 + z4))/27.
      q(17,10)=(y1 + y2 + y4)**3/27.
      q(17,11)=(x1 + x2 + x4)/3.
      q(17,12)=((x1 + x2 + x4)*(z1 + z2 + z4))/9.
      q(17,13)=((x1 + x2 + x4)*(z1 + z2 + z4)**2)/27.
      q(17,14)=((x1 + x2 + x4)*(y1 + y2 + y4))/9.
      q(17,15)=((x1 + x2 + x4)*(y1 + y2 + y4)*(z1 + z2 + z4))/27.
      q(17,16)=((x1 + x2 + x4)*(y1 + y2 + y4)**2)/27.
      q(17,17)=(x1 + x2 + x4)**2/9.
      q(17,18)=((x1 + x2 + x4)**2*(z1 + z2 + z4))/27.
      q(17,19)=((x1 + x2 + x4)**2*(y1 + y2 + y4))/27.
      q(17,20)=(x1 + x2 + x4)**3/27.
      q(18,1)=1
      q(18,2)=(z1 + z2 + z3)/3.
      q(18,3)=(z1 + z2 + z3)**2/9.
      q(18,4)=(z1 + z2 + z3)**3/27.
      q(18,5)=(y1 + y2 + y3)/3.
      q(18,6)=((y1 + y2 + y3)*(z1 + z2 + z3))/9.
      q(18,7)=((y1 + y2 + y3)*(z1 + z2 + z3)**2)/27.
      q(18,8)=(y1 + y2 + y3)**2/9.
      q(18,9)=((y1 + y2 + y3)**2*(z1 + z2 + z3))/27.
      q(18,10)=(y1 + y2 + y3)**3/27.
      q(18,11)=(x1 + x2 + x3)/3.
      q(18,12)=((x1 + x2 + x3)*(z1 + z2 + z3))/9.
      q(18,13)=((x1 + x2 + x3)*(z1 + z2 + z3)**2)/27.
      q(18,14)=((x1 + x2 + x3)*(y1 + y2 + y3))/9.
      q(18,15)=((x1 + x2 + x3)*(y1 + y2 + y3)*(z1 + z2 + z3))/27.
      q(18,16)=((x1 + x2 + x3)*(y1 + y2 + y3)**2)/27.
      q(18,17)=(x1 + x2 + x3)**2/9.
      q(18,18)=((x1 + x2 + x3)**2*(z1 + z2 + z3))/27.
      q(18,19)=((x1 + x2 + x3)**2*(y1 + y2 + y3))/27.
      q(18,20)=(x1 + x2 + x3)**3/27.
      q(19,1)=1
      q(19,2)=(z1 + z3 + z4)/3.
      q(19,3)=(z1 + z3 + z4)**2/9.
      q(19,4)=(z1 + z3 + z4)**3/27.
      q(19,5)=(y1 + y3 + y4)/3.
      q(19,6)=((y1 + y3 + y4)*(z1 + z3 + z4))/9.
      q(19,7)=((y1 + y3 + y4)*(z1 + z3 + z4)**2)/27.
      q(19,8)=(y1 + y3 + y4)**2/9.
      q(19,9)=((y1 + y3 + y4)**2*(z1 + z3 + z4))/27.
      q(19,10)=(y1 + y3 + y4)**3/27.
      q(19,11)=(x1 + x3 + x4)/3.
      q(19,12)=((x1 + x3 + x4)*(z1 + z3 + z4))/9.
      q(19,13)=((x1 + x3 + x4)*(z1 + z3 + z4)**2)/27.
      q(19,14)=((x1 + x3 + x4)*(y1 + y3 + y4))/9.
      q(19,15)=((x1 + x3 + x4)*(y1 + y3 + y4)*(z1 + z3 + z4))/27.
      q(19,16)=((x1 + x3 + x4)*(y1 + y3 + y4)**2)/27.
      q(19,17)=(x1 + x3 + x4)**2/9.
      q(19,18)=((x1 + x3 + x4)**2*(z1 + z3 + z4))/27.
      q(19,19)=((x1 + x3 + x4)**2*(y1 + y3 + y4))/27.
      q(19,20)=(x1 + x3 + x4)**3/27.
      q(20,1)=1
      q(20,2)=(z2 + z3 + z4)/3.
      q(20,3)=(z2 + z3 + z4)**2/9.
      q(20,4)=(z2 + z3 + z4)**3/27.
      q(20,5)=(y2 + y3 + y4)/3.
      q(20,6)=((y2 + y3 + y4)*(z2 + z3 + z4))/9.
      q(20,7)=((y2 + y3 + y4)*(z2 + z3 + z4)**2)/27.
      q(20,8)=(y2 + y3 + y4)**2/9.
      q(20,9)=((y2 + y3 + y4)**2*(z2 + z3 + z4))/27.
      q(20,10)=(y2 + y3 + y4)**3/27.
      q(20,11)=(x2 + x3 + x4)/3.
      q(20,12)=((x2 + x3 + x4)*(z2 + z3 + z4))/9.
      q(20,13)=((x2 + x3 + x4)*(z2 + z3 + z4)**2)/27.
      q(20,14)=((x2 + x3 + x4)*(y2 + y3 + y4))/9.
      q(20,15)=((x2 + x3 + x4)*(y2 + y3 + y4)*(z2 + z3 + z4))/27.
      q(20,16)=((x2 + x3 + x4)*(y2 + y3 + y4)**2)/27.
      q(20,17)=(x2 + x3 + x4)**2/9.
      q(20,18)=((x2 + x3 + x4)**2*(z2 + z3 + z4))/27.
      q(20,19)=((x2 + x3 + x4)**2*(y2 + y3 + y4))/27.
      q(20,20)=(x2 + x3 + x4)**3/27.

!     we nned to invert the 20 x 20 matrix to map from the 20 p-DOF
!     to the 20 polynomial coefficients
      LDA=20
      NPIVOT=20
      call PIVOT(q,LDA,NPIVOT,matrixi,PS1,PS2,IERR,MX,MY,VAL)
      
      if(ierr.ne.1) then
        write(*,*)"FATAL: matrix inversion failed "
!        read(*,*)
      end if

!     now unite the 20 x 20 matrix with the 20 x 16 matrix
!     the results maps from the 16 DOF to the 20 polynomial coefficients      
      do i=1,20
      do j=1,16
        matrix3(i,j)=0.d0
        do k=1,20
          matrix3(i,j)=matrix3(i,j)+matrixi(i,k)*matrix2(k,j)
        end do
      end do
      end do

!     prepare numerical integration
      ngp = 29
      
! Integration points and weights for reference element ((0 0 0)(1 0 0)(0 1 0)(0 0 1))
      GP(1,1)=0.25
      GP(1,2)=0.25
      GP(1,3)=0.25
      GP(1,4)=0.015066881743357949738327730999091
      GP(2,1)=0.057426917317356819579978725140823
      GP(2,2)=0.057426917317356819579978725140823
      GP(2,3)=0.057426917317356819579978725140823
      GP(2,4)=0.0031866390464985314763201441565449
      GP(3,1)=0.057426917317356819579978725140823
      GP(3,2)=0.057426917317356819579978725140823
      GP(3,3)=0.827719248047929541260063824577531
      GP(3,4)=0.0031866390464985314763201441565449
      GP(4,1)=0.057426917317356819579978725140823
      GP(4,2)=0.827719248047929541260063824577531
      GP(4,3)=0.057426917317356819579978725140823
      GP(4,4)=0.0031866390464985314763201441565449
      GP(5,1)=0.827719248047929541260063824577531
      GP(5,2)=0.057426917317356819579978725140823
      GP(5,3)=0.057426917317356819579978725140823
      GP(5,4)=0.0031866390464985314763201441565449
      GP(6,1)=0.23129854365191466342385344099185
      GP(6,2)=0.23129854365191466342385344099185
      GP(6,3)=0.48605102857060727870919871076851
      GP(6,4)=0.0072691564011109382427152201950078
      GP(7,1)=0.48605102857060727870919871076851
      GP(7,2)=0.23129854365191466342385344099185
      GP(7,3)=0.05135188412556339444309440724779
      GP(7,4)=0.0072691564011109382427152201950078
      GP(8,1)=0.48605102857060727870919871076851
      GP(8,2)=0.05135188412556339444309440724779
      GP(8,3)=0.23129854365191466342385344099185
      GP(8,4)=0.0072691564011109382427152201950078
      GP(9,1)=0.05135188412556339444309440724779
      GP(9,2)=0.48605102857060727870919871076851
      GP(9,3)=0.23129854365191466342385344099185
      GP(9,4)=0.0072691564011109382427152201950078
      GP(10,1)=0.486051028570607278709198710768507
      GP(10,2)=0.23129854365191466342385344099185
      GP(10,3)=0.23129854365191466342385344099185
      GP(10,4)=0.0072691564011109382427152201950078
      GP(11,1)=0.23129854365191466342385344099185
      GP(11,2)=0.23129854365191466342385344099185
      GP(11,3)=0.05135188412556339444309440724779
      GP(11,4)=0.0072691564011109382427152201950078
      GP(12,1)=0.23129854365191466342385344099185
      GP(12,2)=0.05135188412556339444309440724779
      GP(12,3)=0.486051028570607278709198710768507
      GP(12,4)=0.0072691564011109382427152201950078
      GP(13,1)=0.05135188412556339444309440724779
      GP(13,2)=0.23129854365191466342385344099185
      GP(13,3)=0.23129854365191466342385344099185
      GP(13,4)=0.0072691564011109382427152201950078
      GP(14,1)=0.23129854365191466342385344099185
      GP(14,2)=0.486051028570607278709198710768507
      GP(14,3)=0.23129854365191466342385344099185
      GP(14,4)=0.0072691564011109382427152201950078
      GP(15,1)=0.23129854365191466342385344099185
      GP(15,2)=0.486051028570607278709198710768507
      GP(15,3)=0.05135188412556339444309440724779
      GP(15,4)=0.0072691564011109382427152201950078
      GP(16,1)=0.23129854365191466342385344099185
      GP(16,2)=0.05135188412556339444309440724779
      GP(16,3)=0.23129854365191466342385344099185
      GP(16,4)=0.0072691564011109382427152201950078
      GP(17,1)=0.05135188412556339444309440724779
      GP(17,2)=0.23129854365191466342385344099185
      GP(17,3)=0.486051028570607278709198710768507
      GP(17,4)=0.0072691564011109382427152201950078
      GP(18,1)=0.047569098814722959646021419203138
      GP(18,2)=0.047569098814722959646021419203138
      GP(18,3)=0.60810798940152808607439022163335
      GP(18,4)=0.0043019459936652776758729763917752
      GP(19,1)=0.60810798940152808607439022163335
      GP(19,2)=0.047569098814722959646021419203138
      GP(19,3)=0.29675381296902599463356693996037
      GP(19,4)=0.0043019459936652776758729763917752
      GP(20,1)=0.60810798940152808607439022163335
      GP(20,2)=0.29675381296902599463356693996037
      GP(20,3)=0.047569098814722959646021419203138
      GP(20,4)=0.0043019459936652776758729763917752
      GP(21,1)=0.29675381296902599463356693996037
      GP(21,2)=0.60810798940152808607439022163335
      GP(21,3)=0.047569098814722959646021419203138
      GP(21,4)=0.0043019459936652776758729763917752
      GP(22,1)=0.60810798940152808607439022163335
      GP(22,2)=0.047569098814722959646021419203138
      GP(22,3)=0.047569098814722959646021419203138
      GP(22,4)=0.0043019459936652776758729763917752
      GP(23,1)=0.047569098814722959646021419203138
      GP(23,2)=0.047569098814722959646021419203138
      GP(23,3)=0.29675381296902599463356693996037
      GP(23,4)=0.0043019459936652776758729763917752
      GP(24,1)=0.047569098814722959646021419203138
      GP(24,2)=0.29675381296902599463356693996037
      GP(24,3)=0.60810798940152808607439022163335
      GP(24,4)=0.0043019459936652776758729763917752
      GP(25,1)=0.29675381296902599463356693996037
      GP(25,2)=0.047569098814722959646021419203138
      GP(25,3)=0.047569098814722959646021419203138
      GP(25,4)=0.0043019459936652776758729763917752
      GP(26,1)=0.047569098814722959646021419203138
      GP(26,2)=0.60810798940152808607439022163335
      GP(26,3)=0.047569098814722959646021419203138
      GP(26,4)=0.0043019459936652776758729763917752
      GP(27,1)=0.047569098814722959646021419203138
      GP(27,2)=0.60810798940152808607439022163335
      GP(27,3)=0.29675381296902599463356693996037
      GP(27,4)=0.0043019459936652776758729763917752
      GP(28,1)=0.047569098814722959646021419203138
      GP(28,2)=0.29675381296902599463356693996037
      GP(28,3)=0.047569098814722959646021419203138
      GP(28,4)=0.0043019459936652776758729763917752
      GP(29,1)=0.29675381296902599463356693996037
      GP(29,2)=0.047569098814722959646021419203138
      GP(29,3)=0.60810798940152808607439022163335
      GP(29,4)=0.0043019459936652776758729763917752
      
!     deform reference element to current element
      F(1,1)=coords(1,2)-coords(1,1)
      F(1,2)=coords(1,3)-coords(1,1)
      F(1,3)=coords(1,4)-coords(1,1)

      F(2,1)=coords(2,2)-coords(2,1)
      F(2,2)=coords(2,3)-coords(2,1)
      F(2,3)=coords(2,4)-coords(2,1)

      F(3,1)=coords(3,2)-coords(3,1)
      F(3,2)=coords(3,3)-coords(3,1)
      F(3,3)=coords(3,4)-coords(3,1)
            
!     reference integration points to current integration points,
!     change weight by jacobian 
      detf=det33(F)
      GPakt=0.d0

      do i=1,ngp
        do j=1,3
          do k=1,3
            GPakt(i,j)=GPakt(i,j)+F(j,k)*GP(i,k)
          end do
          GPakt(i,j)=GPakt(i,j)+coords(j,1)
        end do
        GPakt(i,4)=GP(i,4)*detf
      end do

!     initialise everything (piofo = piola for output)
      piofo=0.d0
      force=0.d0
      dforce_ddof=0.d0

!     cycle over integration points and sum up forces and derivatives
      do gpindex=1,ngp
        x=GPakt(gpindex,1)
        y=GPakt(gpindex,2)
        z=GPakt(gpindex,3)
        weight=GPakt(gpindex,4)

!       Get the matrices that are multiplied by the DOF to get the 
!       function value and the first and second gradient at the 
!       integration points
        call Bmatrices(x,y,z,matrix3,B0,B1,B2)
        do k=1,3
          u(k)=0.d0
          do i=1,4 
          do j=1,4
            u(k)=u(k)+dof(k,i,j)*B0(i,j)
          end do
          end do
          do l=1,3
            H1(k,l)=0.d0
            do i=1,4 
            do j=1,4
              H1(k,l)=H1(k,l)+dof(k,i,j)*B1(i,j,l)
            end do
            end do
            do m=1,3
              H2(k,l,m)=0.d0
              do i=1,4 
              do j=1,4
                H2(k,l,m)=H2(k,l,m)+dof(k,i,j)*B2(i,j,l,m)
              end do
              end do
            end do
          end do
        end do

!       pass function value and first and second gradient to material law 
        call umatr(H1,H2,T,T3,DTDH,DTDH3,DT3DH,DT3DH3,pass)
        
!       output
        if(echo.eq.1) then
          write(28,'(9F20.10)') x,y,z,pass(1),pass(2),pass(3),u(1),u(2)
     &                         ,u(3)
        end if

!       piofo is the cell average, it is getting summed and weighed. 
!       the division with detf is necessary to get a volume-specific 
!       output that does not depend on the element size 
        piofo(1:3)=piofo(1:3)+pass(1:3)*weight/detf 
        piofo(4)=detf

!       Calculate residual forces from stresses
        do i=1,3
        do j=1,4
        do k=1,4
          do l=1,3
!                 3 4 4        3 4 4    3 3     4 4 3 
            force(i,j,k)=force(i,j,k)+T(i,l)*B1(j,k,l)*weight
            do m=1,3
!                   3 4 4        3 4 4     3 3 3     4 4 3 3 
              force(i,j,k)=force(i,j,k)+T3(i,l,m)*B2(j,k,l,m)*weight
            end do  
          end do
        end do
        end do
        end do

!       calculate element linearization from material tangent 
!       first PK stress WRT first displacement gradient
        do i=1,3 
        do j=1,4
        do k=1,4
        do l=1,3
        do m=1,4
        do n=1,4
          do s1=1,3
          do s2=1,3
!                       3 4 4 3 4 4              3 4 4 3 4 4     
            dforce_ddof(i,j,k,l,m,n)=dforce_ddof(i,j,k,l,m,n)+
     &       DTDH(i,s1,l,s2)*B1(j,k,s1)*B1(m,n,s2)*weight
!                 3 3| 3 3      4 4 3      4 4 3
          end do
          end do
        end do
        end do
        end do
        end do
        end do
        end do

!       calculate element linearization from material tangent 
!       first PK stress WRT second displacement gradient
        do i=1,3 
        do j=1,4
        do k=1,4
        do l=1,3
        do m=1,4
        do n=1,4
          do s1=1,3
          do s2=1,3
          do s3=1,3
!                       3 4 4 3 4 4              3 4 4 3 4 4     
            dforce_ddof(i,j,k,l,m,n)=dforce_ddof(i,j,k,l,m,n)+
     &       DTDH3(i,s1,l,s2,s3)*B1(j,k,s1)*B2(m,n,s2,s3)*weight
!                  3 3| 3  3  3     4 4 3      4 4 3 3
          end do
          end do
          end do
        end do
        end do
        end do
        end do
        end do
        end do

!       calculate element linearization from material tangent 
!       first PK hyperstress WRT first displacement gradient
        do i=1,3 
        do j=1,4
        do k=1,4
        do l=1,3
        do m=1,4
        do n=1,4
          do s1=1,3
          do s2=1,3
          do s3=1,3
!                       3 4 4 3 4 4              3 4 4 3 4 4     
            dforce_ddof(i,j,k,l,m,n)=dforce_ddof(i,j,k,l,m,n)+
     &       DT3DH(i,s1,s3,l,s2)*B2(j,k,s1,s3)*B1(m,n,s2)*weight
!                  3 3 3 | 3  3     4 4 3 3       4 4 3 
          end do
          end do
          end do
        end do
        end do
        end do
        end do
        end do
        end do
        
!       calculate element linearization from material tangent 
!       first PK hyperstress WRT second displacement gradient
        do i=1,3 
        do j=1,4
        do k=1,4
        do l=1,3
        do m=1,4
        do n=1,4
          do s1=1,3
          do s2=1,3
          do s3=1,3
          do s4=1,3
!                       3 4 4 3 4 4              3 4 4 3 4 4     
            dforce_ddof(i,j,k,l,m,n)=dforce_ddof(i,j,k,l,m,n)+
     &       DT3DH3(i,s1,s3,l,s2,s4)*B2(j,k,s1,s3)*B2(m,n,s2,s4)*weight
!                   3 3 3  | 3 3 3      4 4 3 3       4 4 3 3
          end do
          end do
          end do
          end do
        end do
        end do
        end do
        end do
        end do
        end do
      end do
      
      return
      end

!
!
!    Auxiliary subroutine that returns 20, 20x3 and 20x3x3 matrices 
!    depending on the location such that the product with the DOF
!    gives the field value and the first and second gradient
!
      subroutine Bmatrices(x,y,z,Mi,B0,B1,B2)
      implicit none
      double precision x,y,z
      double precision Mi(20,20)
      double precision B0(4,4),B1(4,4,3),B2(4,4,3,3)

!                      value,gradient,gradientgradient      
      double precision n(20),gn(20,3),ggn(20,3,3)
      
      integer i,j,k,l,m,imat(4,4)
      
      imat(1,1)=1
      imat(1,2)=2
      imat(1,3)=3
      imat(1,4)=4

      imat(2,1)=5
      imat(2,2)=6
      imat(2,3)=7
      imat(2,4)=8

      imat(3,1)=9
      imat(3,2)=10
      imat(3,3)=11
      imat(3,4)=12

      imat(4,1)=13
      imat(4,2)=14
      imat(4,3)=15
      imat(4,4)=16

! n      
! the shape function monomials 
      n(1)=1
      n(2)=z
      n(3)=z**2
      n(4)=z**3
      n(5)=y
      n(6)=y*z
      n(7)=y*z**2
      n(8)=y**2
      n(9)=y**2*z
      n(10)=y**3
      n(11)=x
      n(12)=x*z
      n(13)=x*z**2
      n(14)=x*y
      n(15)=x*y*z
      n(16)=x*y**2
      n(17)=x**2
      n(18)=x**2*z
      n(19)=x**2*y
      n(20)=x**3

! the shape function monomials first gradients 
      gn(1,1)=0
      gn(1,2)=0
      gn(1,3)=0
      gn(2,1)=0
      gn(2,2)=0
      gn(2,3)=1
      gn(3,1)=0
      gn(3,2)=0
      gn(3,3)=2*z
      gn(4,1)=0
      gn(4,2)=0
      gn(4,3)=3*z**2
      gn(5,1)=0
      gn(5,2)=1
      gn(5,3)=0
      gn(6,1)=0
      gn(6,2)=z
      gn(6,3)=y
      gn(7,1)=0
      gn(7,2)=z**2
      gn(7,3)=2*y*z
      gn(8,1)=0
      gn(8,2)=2*y
      gn(8,3)=0
      gn(9,1)=0
      gn(9,2)=2*y*z
      gn(9,3)=y**2
      gn(10,1)=0
      gn(10,2)=3*y**2
      gn(10,3)=0
      gn(11,1)=1
      gn(11,2)=0
      gn(11,3)=0
      gn(12,1)=z
      gn(12,2)=0
      gn(12,3)=x
      gn(13,1)=z**2
      gn(13,2)=0
      gn(13,3)=2*x*z
      gn(14,1)=y
      gn(14,2)=x
      gn(14,3)=0
      gn(15,1)=y*z
      gn(15,2)=x*z
      gn(15,3)=x*y
      gn(16,1)=y**2
      gn(16,2)=2*x*y
      gn(16,3)=0
      gn(17,1)=2*x
      gn(17,2)=0
      gn(17,3)=0
      gn(18,1)=2*x*z
      gn(18,2)=0
      gn(18,3)=x**2
      gn(19,1)=2*x*y
      gn(19,2)=x**2
      gn(19,3)=0
      gn(20,1)=3*x**2
      gn(20,2)=0
      gn(20,3)=0
! the shape function monomials second gradients 
      ggn(1,1,1)=0
      ggn(1,1,2)=0
      ggn(1,1,3)=0
      ggn(1,2,1)=0
      ggn(1,2,2)=0
      ggn(1,2,3)=0
      ggn(1,3,1)=0
      ggn(1,3,2)=0
      ggn(1,3,3)=0
      ggn(2,1,1)=0
      ggn(2,1,2)=0
      ggn(2,1,3)=0
      ggn(2,2,1)=0
      ggn(2,2,2)=0
      ggn(2,2,3)=0
      ggn(2,3,1)=0
      ggn(2,3,2)=0
      ggn(2,3,3)=0
      ggn(3,1,1)=0
      ggn(3,1,2)=0
      ggn(3,1,3)=0
      ggn(3,2,1)=0
      ggn(3,2,2)=0
      ggn(3,2,3)=0
      ggn(3,3,1)=0
      ggn(3,3,2)=0
      ggn(3,3,3)=2
      ggn(4,1,1)=0
      ggn(4,1,2)=0
      ggn(4,1,3)=0
      ggn(4,2,1)=0
      ggn(4,2,2)=0
      ggn(4,2,3)=0
      ggn(4,3,1)=0
      ggn(4,3,2)=0
      ggn(4,3,3)=6*z
      ggn(5,1,1)=0
      ggn(5,1,2)=0
      ggn(5,1,3)=0
      ggn(5,2,1)=0
      ggn(5,2,2)=0
      ggn(5,2,3)=0
      ggn(5,3,1)=0
      ggn(5,3,2)=0
      ggn(5,3,3)=0
      ggn(6,1,1)=0
      ggn(6,1,2)=0
      ggn(6,1,3)=0
      ggn(6,2,1)=0
      ggn(6,2,2)=0
      ggn(6,2,3)=1
      ggn(6,3,1)=0
      ggn(6,3,2)=1
      ggn(6,3,3)=0
      ggn(7,1,1)=0
      ggn(7,1,2)=0
      ggn(7,1,3)=0
      ggn(7,2,1)=0
      ggn(7,2,2)=0
      ggn(7,2,3)=2*z
      ggn(7,3,1)=0
      ggn(7,3,2)=2*z
      ggn(7,3,3)=2*y
      ggn(8,1,1)=0
      ggn(8,1,2)=0
      ggn(8,1,3)=0
      ggn(8,2,1)=0
      ggn(8,2,2)=2
      ggn(8,2,3)=0
      ggn(8,3,1)=0
      ggn(8,3,2)=0
      ggn(8,3,3)=0
      ggn(9,1,1)=0
      ggn(9,1,2)=0
      ggn(9,1,3)=0
      ggn(9,2,1)=0
      ggn(9,2,2)=2*z
      ggn(9,2,3)=2*y
      ggn(9,3,1)=0
      ggn(9,3,2)=2*y
      ggn(9,3,3)=0
      ggn(10,1,1)=0
      ggn(10,1,2)=0
      ggn(10,1,3)=0
      ggn(10,2,1)=0
      ggn(10,2,2)=6*y
      ggn(10,2,3)=0
      ggn(10,3,1)=0
      ggn(10,3,2)=0
      ggn(10,3,3)=0
      ggn(11,1,1)=0
      ggn(11,1,2)=0
      ggn(11,1,3)=0
      ggn(11,2,1)=0
      ggn(11,2,2)=0
      ggn(11,2,3)=0
      ggn(11,3,1)=0
      ggn(11,3,2)=0
      ggn(11,3,3)=0
      ggn(12,1,1)=0
      ggn(12,1,2)=0
      ggn(12,1,3)=1
      ggn(12,2,1)=0
      ggn(12,2,2)=0
      ggn(12,2,3)=0
      ggn(12,3,1)=1
      ggn(12,3,2)=0
      ggn(12,3,3)=0
      ggn(13,1,1)=0
      ggn(13,1,2)=0
      ggn(13,1,3)=2*z
      ggn(13,2,1)=0
      ggn(13,2,2)=0
      ggn(13,2,3)=0
      ggn(13,3,1)=2*z
      ggn(13,3,2)=0
      ggn(13,3,3)=2*x
      ggn(14,1,1)=0
      ggn(14,1,2)=1
      ggn(14,1,3)=0
      ggn(14,2,1)=1
      ggn(14,2,2)=0
      ggn(14,2,3)=0
      ggn(14,3,1)=0
      ggn(14,3,2)=0
      ggn(14,3,3)=0
      ggn(15,1,1)=0
      ggn(15,1,2)=z
      ggn(15,1,3)=y
      ggn(15,2,1)=z
      ggn(15,2,2)=0
      ggn(15,2,3)=x
      ggn(15,3,1)=y
      ggn(15,3,2)=x
      ggn(15,3,3)=0
      ggn(16,1,1)=0
      ggn(16,1,2)=2*y
      ggn(16,1,3)=0
      ggn(16,2,1)=2*y
      ggn(16,2,2)=2*x
      ggn(16,2,3)=0
      ggn(16,3,1)=0
      ggn(16,3,2)=0
      ggn(16,3,3)=0
      ggn(17,1,1)=2
      ggn(17,1,2)=0
      ggn(17,1,3)=0
      ggn(17,2,1)=0
      ggn(17,2,2)=0
      ggn(17,2,3)=0
      ggn(17,3,1)=0
      ggn(17,3,2)=0
      ggn(17,3,3)=0
      ggn(18,1,1)=2*z
      ggn(18,1,2)=0
      ggn(18,1,3)=2*x
      ggn(18,2,1)=0
      ggn(18,2,2)=0
      ggn(18,2,3)=0
      ggn(18,3,1)=2*x
      ggn(18,3,2)=0
      ggn(18,3,3)=0
      ggn(19,1,1)=2*y
      ggn(19,1,2)=2*x
      ggn(19,1,3)=0
      ggn(19,2,1)=2*x
      ggn(19,2,2)=0
      ggn(19,2,3)=0
      ggn(19,3,1)=0
      ggn(19,3,2)=0
      ggn(19,3,3)=0
      ggn(20,1,1)=6*x
      ggn(20,1,2)=0
      ggn(20,1,3)=0
      ggn(20,2,1)=0
      ggn(20,2,2)=0
      ggn(20,2,3)=0
      ggn(20,3,1)=0
      ggn(20,3,2)=0
      ggn(20,3,3)=0      
! B0 bauen
      do i=1,4
      do j=1,4
        B0(i,j)=0.d0
        do m=1,20
          B0(i,j)=B0(i,j)+n(m)*Mi(m,imat(i,j))
        end do
      end do
      end do      
      
! B1 bauen
      do i=1,4
      do j=1,4
      do k=1,3
        B1(i,j,k)=0.d0
        do m=1,20
          B1(i,j,k)=B1(i,j,k)+gn(m,k)*Mi(m,imat(i,j))
        end do
      end do
      end do      
      end do
      
! B2 bauen
      do i=1,4
      do j=1,4
      do k=1,3
      do l=1,3
        B2(i,j,k,l)=0.d0
        do m=1,20
          B2(i,j,k,l)=B2(i,j,k,l)+ggn(m,k,l)*Mi(m,imat(i,j))
        end do
      end do
      end do      
      end do
      end do

      return
      end

!
!
!   auxiliary subroutine: norm of a 3 x 3 matrix
!
!
      function norm33(f)
      implicit none
      double precision f(3,3),norm33
      norm33=dsqrt(f(1,1)*f(1,1)+f(2,2)*f(2,2)+f(3,3)*f(3,3)+
     &       f(1,2)*f(1,2)+f(1,3)*f(1,3)+f(2,3)*f(2,3)+   
     &       f(2,1)*f(2,1)+f(3,1)*f(3,1)+f(3,2)*f(3,2))
      return
      end

!
!
!   auxiliary subroutine: norm of a 3 x 3 x 3 "matrix"
!
!
      function norm333(f)
      implicit none
      double precision f(3,3,3),norm333
      integer i,j,k
      norm333=0.d0
      do i=1,3
      do j=1,3
      do k=1,3
      norm333=norm333+f(i,j,k)**2
      end do
      end do
      end do
      norm333=dsqrt(norm333)
      return
      end

!
!
!   auxiliary subroutine: determinant of a 3 x 3 matrix
!
!
      function det33(f)
      implicit none
      double precision f(3,3),det33
      det33=f(1,1)*f(2,2)*f(3,3)+
     &      f(1,2)*f(2,3)*f(3,1)+
     &      f(1,3)*f(2,1)*f(3,2)-
     &      f(3,1)*f(2,2)*f(1,3)-
     &      f(3,2)*f(2,3)*f(1,1)-
     &      f(3,3)*f(2,1)*f(1,2)
      return
      end

!
!
!   auxiliary subroutine: print i x j matrix on screen
!
!
      subroutine printmatrix(m,i,j)
      integer i,j,c1,c2
      double precision m(i,j)
      do c1=1,i
      write(*,3333) (m(c1,c2),c2=1,j)
      end do
3333  Format (20(F11.4))
      return
      end


!
!
!   auxiliary subroutine for general N x N matrix inversion
!
! 
      SUBROUTINE PIVOT (A,LDA,N,B,S1,S2,IERR,MX,MY,VAL)

C*****************************************************************
C                                                                *
C  DIESES UNTERPROGRAMM BERECHNET DIE INVERSE EINER REELLEN      *
C  QUADRATISCHEN NXN-MATRIX  A  MIT DEM AUSTAUSCHVERFAHREN       *
C  (AUCH METHODE DES PIVOTISIERENS GENANNT).                     *
C                                                                *
C                                                                *
C  EINGABEPARAMETER:                                             *
C  =================                                             *
C  A        : 2-DIM. FELD (1:LDA,1:N), DAS DIE ZU INVERTIERENDE  *
C             MATRIX A ENTH�T.                                  *
C  LDA      : F�RENDE DIMENSION VON A, WIE IM RUFENDEN          *
C             PROGRAMM VEREINBART.                               *
C  N        : ORDNUNG DER MATRIX A.                              *
C                                                                *
C                                                                *
C  AUSGABEPARAMETER:                                             *
C  =================                                             *
C  B        : 2-DIM. FELD (1:LDA,1:N), DAS DIE INVERSE VON A     *
C             ENTH�T.                                           *
C  S1,S2    : KONTROLLGRENZEN; S1 IST DIE SUMME DER BETR�E      *
C             DER DIAGONALELEMENTE DER MATRIX A*B-E  (E=EINHEITS-*
C             MATRIX),  S2 IST DIE SUMME DER BETR�E ALLER       *
C             �RIGEN ELEMENTE; A*B-E  IST THEORETISCH DIE NULL- *
C             MATRIX.                                            *
C  IERR     : = 1, INVERSE MATRIX VON A GEFUNDEN.                *
C             = 2, MATRIX A IST SINGUL�, DIE INVERSE EXISTIERT  *
C                  NICHT.                                        *
C  VAL      : PIVOTELEMENT, FALLS A ALS SINGUL� ERKANNT WURDE.  *
C                                                                *
C                                                                *
C  HILFSFELDER:                                                  *
C  ============                                                  *
C  MX,MY    : 1-DIM. INTEGER FELDER (1:N).                       *
C                                                                *
C----------------------------------------------------------------*
C                                                                *
C  BEN�IGTE UNTERPROGRAMME: MACHPD                              *
C                                                                *
C*****************************************************************
C                                                                *
C  AUTOR     : GISELA ENGELN-M�LGES                             *
C  DATUM     : 18.05.1987                                        *
C  QUELLCODE : FORTRAN 77                                        *
C                                                                *
C*****************************************************************
C
C  DEKLARATIONEN.
C
      IMPLICIT DOUBLE PRECISION (A-H,O-Z)
      DOUBLE PRECISION A(LDA,N),B(LDA,N)
      INTEGER MX(N),MY(N)
C
C  BERECHNUNG DES MASCHINENRUNDUNGSFEHLERS FMACHP.
C
      FMACHP=1.0D0
    5 FMACHP=0.5D0*FMACHP
      IF(MACHPD(1.0D0+FMACHP) .EQ. 1) GOTO 5
      FMACHP=FMACHP*2.0D0
C
C  UMSPEICHERN DER MATRIX A AUF B.
C


      DO 10 I=1,N
         DO 20 J=1,N
            B(I,J)=A(I,J)
   20    CONTINUE
   10 CONTINUE
C
C  VORBESETZEN DER PIVOTVEKTOREN MX UND MY MIT NULL.
C
      DO 30 I=1,N
         MX(I)=0
         MY(I)=0
   30 CONTINUE
C
C  BESTIMMUNG DES PIVOTELEMENTES.
C
      DO 40 I=1,N
         PIVO=0.0D0
         DO 50 IX=1,N
            IF(MX(IX) .EQ. 0) THEN
               DO 60 IY=1,N
                  IF(MY(IY) .EQ. 0) THEN
                     IF(DABS(B(IX,IY)) .GT. DABS(PIVO)) THEN
                        PIVO=B(IX,IY)
                        NX=IX
                        NY=IY
                     END IF
                  END IF
   60          CONTINUE
            END IF
   50    CONTINUE
C
C  FALLS DAS PIVOTELEMENT NULL IST, IST DIE MATRIX SINGUL�.
C
         IF(DABS(PIVO) .LT. 4.0D0*FMACHP) THEN
            VAL=PIVO
            IERR=2
            RETURN
         END IF
C
C  MERKEN DER INDIZES DES PIVOTELEMENTES.
C
         MX(NX)=NY
         MY(NY)=NX
C
C  BERECHNUNG DER MATRIXELEMENTE GEM� DER RECHENREGELN F�
C  EINEN AUSTAUSCHSCHRITT.
C
         DUMMY=1.0D0/PIVO
         DO 70 J=1,N
            IF(J .NE. NX) THEN
               FACTOR=B(J,NY)*DUMMY
               DO 80 K=1,N
                  B(J,K)=B(J,K)-B(NX,K)*FACTOR
   80          CONTINUE
               B(J,NY)=-FACTOR
            END IF
   70    CONTINUE
         DO 90 K=1,N
            B(NX,K)=B(NX,K)*DUMMY
   90    CONTINUE
         B(NX,NY)=DUMMY
   40 CONTINUE
C
C  ZEILEN- UND SPALTENVERTAUSCHUNGEN R�KG�GIG MACHEN.
C
      DO 100 I=1,N-1
         DO 110 M=I,N
            IF(MX(M) .EQ. I) GOTO 120
  110    CONTINUE
  120    J=M
         IF(J .NE. I) THEN
            DO 130 K=1,N
               H=B(I,K)
               B(I,K)=B(J,K)
               B(J,K)=H
  130       CONTINUE
            MX(J)=MX(I)
            MX(I)=I
         END IF
         DO 140 M=I,N
            IF(MY(M) .EQ. I) GOTO 150
  140    CONTINUE
  150    J=M
         IF(J .NE. I) THEN
            DO 160 K=1,N
               H=B(K,I)
               B(K,I)=B(K,J)
               B(K,J)=H
  160       CONTINUE
            MY(J)=MY(I)
            MY(I)=I
         END IF
  100 CONTINUE
C
C  BILDUNG DER DIFFERENZMATRIX S=(A*B-E), E=EINHEITSMATRIX.
C  BILDUNG DER SUMME S1 DER BETR�E DER DIAGONALELEMENTE VON S
C  UND DER SUMME S2 DER BETR�E ALLER �RIGEN GLIEDER.
C  THEORETISCH M�TEN S1 UND S2 NULL SEIN.
C
      S1=0.0D0
      S2=0.0D0
      DO 170 I=1,N
         DO 180 J=1,N
            H=0.0D0
            DO 190 K=1,N
               H=H+A(I,K)*B(K,J)
  190       CONTINUE
            IF(I .EQ. J) THEN
               S1=S1+DABS(H-1.0D0)
            ELSE
               S2=S2+DABS(H)
            END IF
  180    CONTINUE
  170 CONTINUE
      IERR=1
      RETURN
      END

      INTEGER FUNCTION MACHPD(X)      
      DOUBLE PRECISION X      
      MACHPD=0
      IF (1.0D0 .LT. X) MACHPD=1      
      RETURN      
      END

!
! Material subroutine UMATR
! This material subroutine states the linearizations explicitely. 
!
! INPUT:    H(3 x 3)        first displacement gradient u o Nabla
!           H3(3 x 3 x 3)   second displacement gradient u o Nabla o Nablba
!
! OUTPUT:   T : first Piola-Kirchhoff-stresses (power-conjugate to H)
!           T3 : first Piola-Kirchhoff-hyperstresses (power-conjugate to  H3)
!           DxDy : derivative of x WRT y (right index gradients)
!       
! The material law:  Starting from Ciarlet-Geymonat strain energy, we extend 
!                    a canonical term to account for a strain gradient contribution
!                    to the elastic energy:
!
!        lambda                   mu                   alpha
!   w = -------- (III-ln(III)) + ---- (I - ln(III)) + ------- CCC ... CCC
!           4                     2                      2
!
!   
!           classical strain energy                 |  strain gradient extension
!
!
! with the first and third invariant I,III of C (resp. B, right and left Cauchy-Green
! tensors), Lame's constants lambda und mu and the right Cauchy-Green-hypertensor 
!
!       CCC= F^T . grad(grad(x)),
!
! the norm of which is penalized energetically. Note that unlike the classical
! right Cauchy-Green tensor, CCC vanishes when no deformation occurs.
!
! Of w: the first summand contains the dilatational contribution, the second contains
! the distortion contribution. 
!
! The second displacement gradient H and the second deformation gradient are equal.
! The following holds:
!
! F=H+1
! F3=grad(grad(x)=grad(grad(u))=H3
! dx/dF=dx/dH
! dx/dF3=dx/dH3
!
! CCC=F^T.F3=F^T.H3
!
! The elastic law results from the derivatives of w WRT H and H3.
!
! T=dw/dF=mu*F + (0.5d0*lambda*(detC-1.d0)-mu)*F^(-T)
!
! T3=dw/dF3=alpha*F.CCC=alpha B.F2=alpha B.H3
!
! The tangents are likewise obtained as the second derivative and are
! summarized in the accompanying PDF.

      subroutine umatr(H,F3,T,T3,DTDH,DTDH3,DT3DH,DT3DH3,pass)
      implicit none

! In- and Output
      double precision H(3,3),F3(3,3,3),T(3,3),T3(3,3,3)
      double precision DTDH(3,3,3,3)
      double precision DTDH3(3,3,3,3,3)
      double precision DT3DH(3,3,3,3,3)
      double precision DT3DH3(3,3,3,3,3,3)
      
! Material parameters
      double precision lambda,mu,alpha,pass(3)
      
! Auxiliary quantities 
      double precision B(3,3)
      double precision F(3,3),Finv(3,3)
      double precision C3(3,3,3)
      double precision delta(3,3)
      double precision detC,detF,det33
      double precision tmp1,tmp2,g
      double precision cauchy(3,3),cauchydev(3,3)
      double precision pressure,mises
      
      integer i,j,k,l,m,n,o,p
      
      lambda=200000.d0*0.3d0/(1.3d0*0.4d0)
      mu=    200000.d0/2.6d0
      alpha= 0.d0

! Identity tensor
      do i=1,3
      do j=1,3
        if(i.eq.j)then
          delta(i,j)=1.d0
        else
          delta(i,j)=0.d0
        end if
      end do
      end do

! various auxiliary quantities
      F=H+delta
      B=matmul(F,transpose(F))
      call invert3(F,Finv)
      detC=det33(B)
      do i=1,3
      do j=1,3
      do k=1,3
        C3(i,j,k)=0.d0
        do l=1,3
          C3(i,j,k)=C3(i,j,k)+F(l,i)*F3(l,j,k)
        end do
      end do
      end do
      end do

! strain energy
!        lambda                   mu                   alpha
!   w = -------- (III-ln(III)) + ---- (I - ln(III)) + ------- CCC ... CCC
!           4                     2                      2
    
      g=0.d0
      do i=1,3
      do j=1,3
      do k=1,3
        g=g+C3(i,j,k)*C3(i,j,k)
      end do
      end do
      end do
    
      pass(1)=0.25d0*lambda*(detc-dlog(detc)-1.d0)+
     &        0.5d0*mu*(B(1,1)+B(2,2)+B(3,3)-dlog(detc)-3.d0)
      pass(2)=0.5d0*alpha*g
     
! calculate first PK stresses
      T=mu*F+(0.5d0*lambda*(detC-1.d0)-mu)*transpose(Finv)
      do i=1,3
      do j=1,3
        do k=1,3
        do l=1,3
          T(i,j)=T(i,j)+alpha*F3(i,k,l)*C3(j,k,l)
        end do
        end do
      end do
      end do

! calculate first PK hyperstresses
      do i=1,3
      do j=1,3
      do k=1,3
        T3(i,j,k)=0.d0
        do l=1,3
          T3(i,j,k)=T3(i,j,k)+alpha*F(i,l)*C3(l,j,k)
        end do
      end do
      end do
      end do

! dT/dF 
      tmp2=lambda*detC
      tmp1=0.5d0*tmp2-0.5d0*lambda-mu
     
      do i=1,3
      do j=1,3
      do k=1,3
      do l=1,3
        DTDH(i,j,k,l)=mu*delta(i,k)*delta(j,l)+
     &  -tmp1*Finv(j,k)*Finv(l,i) + 
     &  +tmp2*Finv(j,i)*Finv(l,k)
      end do
      end do
      end do
      end do

      do m=1,3
      do n=1,3
      do o=1,3
      do p=1,3
        do j=1,3
        do l=1,3
          DTDH(m,n,o,p)=DTDH(m,n,o,p)+
     &        alpha*F3(m,j,l)*F3(o,j,l)*delta(n,p)
        end do
        end do
      end do
      end do
      end do
      end do


! dT/dF3 
      do i=1,3
      do j=1,3
      do k=1,3
      do l=1,3
      do m=1,3
        DTDH3(i,j,k,l,m)=alpha*(delta(i,k)*C3(j,l,m)+F(k,j)*F3(i,l,m))
      end do
      end do
      end do
      end do
      end do
      
! dT3/dF 
      do i=1,3
      do j=1,3
      do k=1,3
      do l=1,3
      do m=1,3
        DT3DH(i,j,k,l,m)=alpha*(delta(i,l)*C3(m,j,k)+F(i,m)*F3(l,j,k))
      end do
      end do
      end do
      end do
      end do
      
! dT3/dF3 
      do i=1,3
      do j=1,3
      do k=1,3
      do l=1,3
      do m=1,3
      do n=1,3
        DT3DH3(i,j,k,l,m,n)=alpha*(B(i,l)*delta(j,m)*delta(k,n))
      end do
      end do
      end do
      end do
      end do
      end do

!     Mises equivalent stress
      cauchy=matmul(T,transpose(F))/sqrt(detC)
      pressure=-(cauchy(1,1)+cauchy(2,2)+cauchy(3,3))/3.d0
      cauchydev=cauchy+delta*pressure
      mises=dsqrt(cauchy(1,1)**2+
     &            cauchy(2,2)**2+
     &            cauchy(3,3)**2+
     &      2.d0* cauchy(1,2)**2+
     &      2.d0* cauchy(1,3)**2+
     &      2.d0* cauchy(2,3)**2) * dsqrt(1.5d0)
      
      
!      pass(2)=pressure
      pass(3)=mises 
              
      return
      end


! c**** inverse matrix mi der 3x3 matrix m
      subroutine invert3(m,mi)
      implicit none
      double precision m(3,3),mi(3,3),hi3m,det33
      hi3m=det33(m)
      if(abs(hi3m).gt.1.d-12)then
        mi(1,1)=(m(2,2)*m(3,3)-m(3,2)*m(2,3))/hi3m
        mi(2,1)=(m(3,1)*m(2,3)-m(2,1)*m(3,3))/hi3m
        mi(3,1)=(m(2,1)*m(3,2)-m(3,1)*m(2,2))/hi3m
        mi(1,2)=(m(3,2)*m(1,3)-m(1,2)*m(3,3))/hi3m
        mi(2,2)=(m(1,1)*m(3,3)-m(3,1)*m(1,3))/hi3m
        mi(3,2)=(m(3,1)*m(1,2)-m(1,1)*m(3,2))/hi3m
        mi(1,3)=(m(1,2)*m(2,3)-m(2,2)*m(1,3))/hi3m
        mi(2,3)=(m(2,1)*m(1,3)-m(1,1)*m(2,3))/hi3m
        mi(3,3)=(m(1,1)*m(2,2)-m(2,1)*m(1,2))/hi3m
      end if
      return
      end


! Material subroutine UMATR2
! Requires only the definition of the strain energy, the material tangent
! is determined numerically by finite differences
!
! INPUT:    H(3 x 3)        first displacement gradient u o Nabla
!           H3(3 x 3 x 3)   second displacement gradient u o Nabla o Nablba
!
! OUTPUT:   T : first Piola-Kirchhoff-stresses (power-conjugate to H)
!           T3 : first Piola-Kirchhoff-hyperstresses (power-conjugate to  H3)
!           DxDy : derivative of x WRT y (right index gradients)

      subroutine umatr2(H,F3,T,T3,DTDH,DTDH3,DT3DH,DT3DH3,pass)
      implicit none

! In- and Output
      double precision H(3,3),F3(3,3,3),T(3,3),T3(3,3,3)
      double precision DTDH(3,3,3,3)
      double precision DTDH3(3,3,3,3,3)
      double precision DT3DH(3,3,3,3,3)
      double precision DT3DH3(3,3,3,3,3,3)
      
! Material parameters
      double precision lambda,mu,alpha,pass(3),emo,nu
      
! Auxiliary quantities
      double precision B(3,3)
      double precision F(3,3),Finv(3,3)
      double precision C3(3,3,3)
      double precision delta(3,3)
      double precision detC,detF,det33
      double precision tmp1,tmp2,g
      double precision cauchy(3,3),cauchydev(3,3)
      double precision pressure,mises,diff
      double precision T0(3,3),T30(3,3,3)
      
      integer i,j,k,l,m,n,o,p
      
      emo=200000.d0
      nu=0.3d0
      lambda=emo*nu/(1.d0+nu)/(1.d0-2.d0*nu)
      mu=    emo/2.d0/(1.d0+nu)
c~       alpha= 200000.0d0
      
! finite difference
      diff=1.d-3

! Identity tensor
      do i=1,3
      do j=1,3
        if(i.eq.j)then
          delta(i,j)=1.d0
        else
          delta(i,j)=0.d0
        end if
      end do
      end do

! Get stresses 
      call getstresses(H,F3,T0,T30,delta,lambda,mu,alpha,pass)

! Get linearization by finite differences
      do i=1,3
      do j=1,3
        H(i,j)=H(i,j)+diff
        call getstresses(H,F3,T,T3,delta,lambda,mu,alpha,pass)
        H(i,j)=H(i,j)-diff
        do k=1,3
        do l=1,3
          dTdH(k,l,i,j)=(T(k,l)-T0(k,l))/diff
        end do
        end do
        do k=1,3
        do l=1,3
        do m=1,3
          dT3dH(k,l,m,i,j)=(T3(k,l,m)-T30(k,l,m))/diff
        end do
        end do
        end do
      end do
      end do
      
      do i=1,3
      do j=1,3
      do k=1,3
        F3(i,j,k)=F3(i,j,k)+diff
        call getstresses(H,F3,T,T3,delta,lambda,mu,alpha,pass)
        F3(i,j,k)=F3(i,j,k)-diff
        do l=1,3
        do m=1,3
          dTdH3(l,m,i,j,k)=(T(l,m)-T0(l,m))/diff
        end do
        end do
        do l=1,3
        do m=1,3
        do n=1,3
          dT3dH3(l,m,n,i,j,k)=(T3(l,m,n)-T30(l,m,n))/diff
        end do
        end do
        end do
      end do
      end do
      end do
    
! get strain energy     
      call getenergy(H,F3,T,T3,delta,lambda,mu,alpha,pass)
      
      return
      end


! Getstresses: gets stresses by finite differences dw/dH
      subroutine getstresses(H,F3,T,T3,delta,lambda,mu,alpha,pifo)
      implicit none
      double precision H(3,3),F3(3,3,3),T(3,3),T3(3,3,3)
      double precision alpha,delta(3,3),lambda,mu,pifo(2)
      
      double precision diff,g0,pass(3)
      integer i,j,k
      
      diff=1.d-3
      
      call getenergy(H,F3,T,T3,delta,lambda,mu,alpha,pass)
      pifo(1:2)=pass(1:2)
      g0=pass(3)
      
      do i=1,3
      do j=1,3
        H(i,j)=H(i,j)+diff
        call getenergy(H,F3,T,T3,delta,lambda,mu,alpha,pass)
        T(i,j)=(pass(3)-g0)/diff
        H(i,j)=H(i,j)-diff
      end do
      end do
    
      do i=1,3
      do j=1,3
      do k=1,3
        F3(i,j,k)=F3(i,j,k)+diff
        call getenergy(H,F3,T,T3,delta,lambda,mu,alpha,pass)
        T3(i,j,k)=(pass(3)-g0)/diff
        F3(i,j,k)=F3(i,j,k)-diff
      end do
      end do
      end do
      
      return
      end


! Define strain energy for UMATR2
      subroutine getenergy(H,F3,T,T3,delta,lambda,mu,alpha,pass)
      implicit none
      double precision H(3,3),F3(3,3,3),T(3,3),T3(3,3,3)
      double precision alpha,delta(3,3),lambda,mu,pass(3)
      
      double precision F(3,3),B(3,3),Finv(3,3),detC
      double precision rd,g,det33 

      integer i,j,k,l
      
      double precision dehn2(3,3),dehn3(3,3,3),c1,c2,c3,c4,c5
      double precision v12(3),v3(3)
      
      double precision kappa,sk,ck,l1,l2,l3,l4,l5,sqrt5

!     from H to F      
      F=H+delta
      
!     E_G o Nabla corresponds sym_12(F^T.F3)
      do i=1,3
      do j=1,3
      do k=1,3
           dehn3(i,j,k)=0.5d0*(F3(i,j,k)+F3(j,i,k))
      end do
      end do
      end do

!     E_G=1/2 (F^T F - I)
      dehn2=0.5d0*(matmul(transpose(F),F)-delta)
      
!     classical el. energy in pass(1)
      pass(1)=0.d0
      do i=1,3
      do j=1,3
        pass(1)=pass(1)+dehn2(i,j)**2.d0
      end do
      end do
      pass(1)=pass(1)*2.d0*mu
      pass(1)=pass(1)+lambda*(dehn2(1,1)+dehn2(2,2)+dehn2(3,3))**2.d0
      pass(1)=pass(1)*0.5d0
      
!     non-classical el. energy: v12=E3_ijj=E3_jij
!                                         v3=E3_jji
      v12=0.d0
      v3=0.d0
      do i=1,3
      do j=1,3
        v12(i)=v12(i)+dehn3(i,j,j)
        v3(i)=v3(i)+dehn3(j,j,i)
      end do
      end do

      g=0.d0
      rd=0.d0
      do i=1,3
      do j=1,3
      do k=1,3
        g=g+dehn3(i,j,k)*dehn3(i,k,j)
        rd=rd+dehn3(i,j,k)*dehn3(i,j,k)
      end do
      end do
      end do
      
!     lambda1 = 2 c4 - 2 c3
!     lambda2 = 8 c3 + 2 c4
!     lambda3 = 2 c1+4 c2+c3+2 c4+(3 c5)/2-1/2 Sqrt[4 (24 c1^2+16 c1 c2+16 c2^2+28 c1 c3+16 c2 c3+9 c3^2)+4 (6 c1-8 c2+c3) c5+9 c5^2]
!     lambda4 = 1/2 (4 c1+8 c2+2 c3+4 c4+3 c5+\[Sqrt](4 (24 c1^2+16 c1 c2+16 c2^2+28 c1 c3+16 c2 c3+9 c3^2)+4 (6 c1-8 c2+c3) c5+9 c5^2))

!     c2 -> 1/240 (-60 c1 - 10 l1 - 20 l2 + 15 (l3 + l4) + 
!           Sqrt[5] \[Sqrt](-3600 c1^2 - 100 l1^2 - 16 l2^2 - 24 l2 l3 + 
!           45 l3^2 - 6 (4 l2 + 21 l3) l4 + 45 l4^2 + 
!           120 c1 (10 l1 - 4 l2 - 3 (l3 + l4)) + 
!           20 l1 (4 l2 + 3 (l3 + l4)))), 
!     c3 -> 1/6 (-l1 + l2), 
!     c4 -> 1/6 (2 l1 + l2), 
!     c5 -> 1/90 (-60 c1 - 20 l1 - 10 l2 + 15 (l3 + l4) - 
!           Sqrt[5] \[Sqrt](-3600 c1^2 - 100 l1^2 - 16 l2^2 - 24 l2 l3 + 
!           45 l3^2 - 6 (4 l2 + 21 l3) l4 + 45 l4^2 + 
!           120 c1 (10 l1 - 4 l2 - 3 (l3 + l4)) + 
!           20 l1 (4 l2 + 3 (l3 + l4))))


! Parameter c1 bis c5: EW l1 bis l4

!       0 | 0      | 0     | 0.5  |  33   :  1 | 1 | 100 | 1
!    16.5 |  -8.25 | -16.5 | 33.5 | -33   :  100 | 1 | 1 | 1
!    -6.6 |   -6.6 | 16.5  | 17   | -6.6  :  1 | 100 | 1 | 1
!    -9.9 |  14.85 |    0  | 0.5  |  39,6 :  1 | 1 | 100 | 100


!      c1=1000.d0
!      c2=1000.d0
!      c3=1000.d0
!      c4=2000.d0
!      c5=1000.d0 ! liefert 3x2300, 7x800,8x200 als EW      

!      c1=0.0
!      c2=0.0
!      c3=0.0
!      c4=0.5
!      c5=33.0

      
! Eigenwerte:

! 1 - 5-fach, Projektor fest
! 2 - 7-fach, Projektor fest
! 3 - 3-fach, Projektor hängt vom Winkel kappa ab
! 4 - 3-fach, Projektor hängt vom Winkel kappa ab

      l1=1.d-3
      l2=1.d-3
      l3=1.d-3
      l4=1.d-3 
      
      kappa=3.141d0/4.d0
      sk=dsin(kappa)
      ck=dcos(kappa)
      sqrt5=dsqrt(5.d0)
            
      c1=l1/6.d0-l2/15.d0-(l3+l4)/20.d0+
     &   (l4-l3)*(sk-sqrt5*ck)/20.d0
      c2=(-10.d0*l1-8.d0*l2+9.d0*(l3+l4)-9.d0*(l4-l3)*sk)/120.d0
      c3=(l2-l1)/6.d0
      c4=(2.d0*l1+l2)/6.d0
      c5=(-15.d0*l1-3.d0*l2+9.d0*(l4+l3)+
     &   3.d0*(l4-l3)*(2.d0*sk+sqrt5*ck) )/45.d0

      pass(2)=  4.d0*c1*dot_product(v12,v3)+
     &          4.d0*c2*dot_product(v12,v12)+
     &          4.d0*c3*g+
     &          2.d0*c4*rd+
     &               c5*dot_product(v3,v3)

      pass(2)=10000.d0*(dehn3(1,2,1)**2.d0+
     &                  dehn3(2,1,1)**2.d0)
 
      pass(2)=0.5d0*pass(2)
      
!     sum classical and non-classical
      pass(3)=pass(1)+pass(2)
     
      return
      end







!
!
!     Programm for testing the user element
!
!


      PROGRAM GETAMATRIX
      implicit none
      
      INTEGER MLVARX,NDOFEL,MDLOAD,MCRD,NNODE,NSVARS,NRHS,NDLOAD
      INTEGER NPROPS,JTYPE,JELEM,KSTEP,KINC,NJPROP,NPREDF
      INTEGER JDLTYP(1,1),LFLAGS(1),JPROPS(1)
      DOUBLE PRECISION RHS(48,1),AMATRX(48,48),
     1     SVARS(1),ENERGY(8),PROPS(1),COORDS(3,4),
     2     U(48),DU(48,1),V(48),A(48),TIME(2),
     3     PARAMS(1),ADLMAG(1,1),PNEWDT,DTIME,PERIOD,
     4     DDLMAG(1,1),PREDEF(2,1,4)
     
      double precision rhs0(48,1),amatn(48,48),delta
c~       character(len=100) :: filename
      
      integer i,j

c~       filename="~/dec/uni/Projekte/DOF_extension_in_ABQ/Paper1/FERechnung
c~      &en/ecke/result.txt"
     
c~       write(*,*)trim(filename)
c~       write(*,*)"a"
c~       read(*,*)
     
      MDLOAD=1
      MLVARX=48
      NDOFEL=48
      NPREDF=1
      NNODE=4
      MCRD=3
      NSVARS=1
      
      COORDS(1,1)=0.d0
      COORDS(2,1)=0.d0
      COORDS(3,1)=0.d0
      COORDS(1,2)=1.d0
      COORDS(2,2)=0.d0
      COORDS(3,2)=0.d0
      COORDS(1,3)=0.d0
      COORDS(2,3)=1.d0
      COORDS(3,3)=0.d0
      COORDS(1,4)=0.d0
      COORDS(2,4)=0.d0
      COORDS(3,4)=1.d0
      
      U=0.d0

      CALL UEL(RHS,AMATRX,SVARS,ENERGY,NDOFEL,NRHS,NSVARS,
     1     PROPS,NPROPS,COORDS,MCRD,NNODE,U,DU,V,A,JTYPE,TIME,DTIME,
     2     KSTEP,KINC,JELEM,PARAMS,NDLOAD,JDLTYP,ADLMAG,PREDEF,
     3     NPREDF,LFLAGS,MLVARX,DDLMAG,MDLOAD,PNEWDT,JPROPS,NJPROP,
     4     PERIOD)
      
      open(unit=28,file="ELSTIFFNESS_A.txt",
     &     action="write",status="replace")
      do i=1,48
        write(28,'(48F20.10)') (AMATRX(i,j),j=1,48)
      end do
      close(28)

      rhs0=rhs
      
      delta=0.0001d0
      do i=1,48
        u=0.d0
        u(i)=delta
        CALL UEL(RHS,AMATRX,SVARS,ENERGY,NDOFEL,NRHS,NSVARS,
     1     PROPS,NPROPS,COORDS,MCRD,NNODE,U,DU,V,A,JTYPE,TIME,DTIME,
     2     KSTEP,KINC,JELEM,PARAMS,NDLOAD,JDLTYP,ADLMAG,PREDEF,
     3     NPREDF,LFLAGS,MLVARX,DDLMAG,MDLOAD,PNEWDT,JPROPS,NJPROP,
     4     PERIOD)
        amatn(i,:)=(rhs0(:,1)-rhs(:,1))/delta
      end do
            
      open(unit=28,file="ELSTIFFNESS_N.txt",
     &     action="write",status="replace")
      do i=1,48
        write(28,'(48F20.10)') (amatn(i,j),j=1,48)
      end do
      close(28)
      end program
