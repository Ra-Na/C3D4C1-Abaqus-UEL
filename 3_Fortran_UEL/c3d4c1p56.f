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
!       The Tilde may cause trouble, replace with your home folder
        dirname="~/"
        filename1="IPVALS56_"
        filename3="UNVALS56_"
        write(nostr,'(I2.2)') kinc
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

      double precision matrixi(56,56)
      double precision q(56,56)
      double precision matrix2(56,16)
      double precision matrix3(56,16)
      

      double precision x1,x2,x3,x4,y1,y2,y3,y4,z1,z2,z3,z4
      double precision phi1,phi1x,phi1y,phi1z
      double precision phi2,phi2x,phi2y,phi2z
      double precision phi3,phi3x,phi3y,phi3z
      double precision phi4,phi4x,phi4y,phi4z
      
      integer LDA,NPIVOT,IERR,MX(56),MY(56),i,j,k,l,m,n,gpindex,i2
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

! matrix2 maps from the 16 dof to the 56 pseudo-dof 
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
      matrix2(21,1)=0
      matrix2(21,2)=0.3333333333333333
      matrix2(21,3)=0
      matrix2(21,4)=0
      matrix2(21,5)=0
      matrix2(21,6)=0.3333333333333333
      matrix2(21,7)=0
      matrix2(21,8)=0
      matrix2(21,9)=0
      matrix2(21,10)=0
      matrix2(21,11)=0
      matrix2(21,12)=0
      matrix2(21,13)=0
      matrix2(21,14)=0.3333333333333333
      matrix2(21,15)=0
      matrix2(21,16)=0
      matrix2(22,1)=0
      matrix2(22,2)=0
      matrix2(22,3)=0.3333333333333333
      matrix2(22,4)=0
      matrix2(22,5)=0
      matrix2(22,6)=0
      matrix2(22,7)=0.3333333333333333
      matrix2(22,8)=0
      matrix2(22,9)=0
      matrix2(22,10)=0
      matrix2(22,11)=0
      matrix2(22,12)=0
      matrix2(22,13)=0
      matrix2(22,14)=0
      matrix2(22,15)=0.3333333333333333
      matrix2(22,16)=0
      matrix2(23,1)=0
      matrix2(23,2)=0
      matrix2(23,3)=0
      matrix2(23,4)=0.3333333333333333
      matrix2(23,5)=0
      matrix2(23,6)=0
      matrix2(23,7)=0
      matrix2(23,8)=0.3333333333333333
      matrix2(23,9)=0
      matrix2(23,10)=0
      matrix2(23,11)=0
      matrix2(23,12)=0
      matrix2(23,13)=0
      matrix2(23,14)=0
      matrix2(23,15)=0
      matrix2(23,16)=0.3333333333333333
      matrix2(24,1)=0
      matrix2(24,2)=0.3333333333333333
      matrix2(24,3)=0
      matrix2(24,4)=0
      matrix2(24,5)=0
      matrix2(24,6)=0.3333333333333333
      matrix2(24,7)=0
      matrix2(24,8)=0
      matrix2(24,9)=0
      matrix2(24,10)=0.3333333333333333
      matrix2(24,11)=0
      matrix2(24,12)=0
      matrix2(24,13)=0
      matrix2(24,14)=0
      matrix2(24,15)=0
      matrix2(24,16)=0
      matrix2(25,1)=0
      matrix2(25,2)=0
      matrix2(25,3)=0.3333333333333333
      matrix2(25,4)=0
      matrix2(25,5)=0
      matrix2(25,6)=0
      matrix2(25,7)=0.3333333333333333
      matrix2(25,8)=0
      matrix2(25,9)=0
      matrix2(25,10)=0
      matrix2(25,11)=0.3333333333333333
      matrix2(25,12)=0
      matrix2(25,13)=0
      matrix2(25,14)=0
      matrix2(25,15)=0
      matrix2(25,16)=0
      matrix2(26,1)=0
      matrix2(26,2)=0
      matrix2(26,3)=0
      matrix2(26,4)=0.3333333333333333
      matrix2(26,5)=0
      matrix2(26,6)=0
      matrix2(26,7)=0
      matrix2(26,8)=0.3333333333333333
      matrix2(26,9)=0
      matrix2(26,10)=0
      matrix2(26,11)=0
      matrix2(26,12)=0.3333333333333333
      matrix2(26,13)=0
      matrix2(26,14)=0
      matrix2(26,15)=0
      matrix2(26,16)=0
      matrix2(27,1)=0
      matrix2(27,2)=0
      matrix2(27,3)=0
      matrix2(27,4)=0
      matrix2(27,5)=0
      matrix2(27,6)=0.3333333333333333
      matrix2(27,7)=0
      matrix2(27,8)=0
      matrix2(27,9)=0
      matrix2(27,10)=0.3333333333333333
      matrix2(27,11)=0
      matrix2(27,12)=0
      matrix2(27,13)=0
      matrix2(27,14)=0.3333333333333333
      matrix2(27,15)=0
      matrix2(27,16)=0
      matrix2(28,1)=0
      matrix2(28,2)=0
      matrix2(28,3)=0
      matrix2(28,4)=0
      matrix2(28,5)=0
      matrix2(28,6)=0
      matrix2(28,7)=0.3333333333333333
      matrix2(28,8)=0
      matrix2(28,9)=0
      matrix2(28,10)=0
      matrix2(28,11)=0.3333333333333333
      matrix2(28,12)=0
      matrix2(28,13)=0
      matrix2(28,14)=0
      matrix2(28,15)=0.3333333333333333
      matrix2(28,16)=0
      matrix2(29,1)=0
      matrix2(29,2)=0
      matrix2(29,3)=0
      matrix2(29,4)=0
      matrix2(29,5)=0
      matrix2(29,6)=0
      matrix2(29,7)=0
      matrix2(29,8)=0.3333333333333333
      matrix2(29,9)=0
      matrix2(29,10)=0
      matrix2(29,11)=0
      matrix2(29,12)=0.3333333333333333
      matrix2(29,13)=0
      matrix2(29,14)=0
      matrix2(29,15)=0
      matrix2(29,16)=0.3333333333333333
      matrix2(30,1)=0
      matrix2(30,2)=0
      matrix2(30,3)=0
      matrix2(30,4)=0
      matrix2(30,5)=0
      matrix2(30,6)=0.3333333333333333
      matrix2(30,7)=0
      matrix2(30,8)=0
      matrix2(30,9)=0
      matrix2(30,10)=0.3333333333333333
      matrix2(30,11)=0
      matrix2(30,12)=0
      matrix2(30,13)=0
      matrix2(30,14)=0.3333333333333333
      matrix2(30,15)=0
      matrix2(30,16)=0
      matrix2(31,1)=0
      matrix2(31,2)=0
      matrix2(31,3)=0
      matrix2(31,4)=0
      matrix2(31,5)=0
      matrix2(31,6)=0
      matrix2(31,7)=0.3333333333333333
      matrix2(31,8)=0
      matrix2(31,9)=0
      matrix2(31,10)=0
      matrix2(31,11)=0.3333333333333333
      matrix2(31,12)=0
      matrix2(31,13)=0
      matrix2(31,14)=0
      matrix2(31,15)=0.3333333333333333
      matrix2(31,16)=0
      matrix2(32,1)=0
      matrix2(32,2)=0
      matrix2(32,3)=0
      matrix2(32,4)=0
      matrix2(32,5)=0
      matrix2(32,6)=0
      matrix2(32,7)=0
      matrix2(32,8)=0.3333333333333333
      matrix2(32,9)=0
      matrix2(32,10)=0
      matrix2(32,11)=0
      matrix2(32,12)=0.3333333333333333
      matrix2(32,13)=0
      matrix2(32,14)=0
      matrix2(32,15)=0
      matrix2(32,16)=0.3333333333333333
      matrix2(33,1)=0.5
      matrix2(33,2)=(-x1 + x2)/8.
      matrix2(33,3)=(-y1 + y2)/8.
      matrix2(33,4)=(-z1 + z2)/8.
      matrix2(33,5)=0.5
      matrix2(33,6)=(x1 - x2)/8.
      matrix2(33,7)=(y1 - y2)/8.
      matrix2(33,8)=(z1 - z2)/8.
      matrix2(33,9)=0
      matrix2(33,10)=0
      matrix2(33,11)=0
      matrix2(33,12)=0
      matrix2(33,13)=0
      matrix2(33,14)=0
      matrix2(33,15)=0
      matrix2(33,16)=0
      matrix2(34,1)=0.5
      matrix2(34,2)=(-x1 + x3)/8.
      matrix2(34,3)=(-y1 + y3)/8.
      matrix2(34,4)=(-z1 + z3)/8.
      matrix2(34,5)=0
      matrix2(34,6)=0
      matrix2(34,7)=0
      matrix2(34,8)=0
      matrix2(34,9)=0.5
      matrix2(34,10)=(x1 - x3)/8.
      matrix2(34,11)=(y1 - y3)/8.
      matrix2(34,12)=(z1 - z3)/8.
      matrix2(34,13)=0
      matrix2(34,14)=0
      matrix2(34,15)=0
      matrix2(34,16)=0
      matrix2(35,1)=0.5
      matrix2(35,2)=(-x1 + x4)/8.
      matrix2(35,3)=(-y1 + y4)/8.
      matrix2(35,4)=(-z1 + z4)/8.
      matrix2(35,5)=0
      matrix2(35,6)=0
      matrix2(35,7)=0
      matrix2(35,8)=0
      matrix2(35,9)=0
      matrix2(35,10)=0
      matrix2(35,11)=0
      matrix2(35,12)=0
      matrix2(35,13)=0.5
      matrix2(35,14)=(x1 - x4)/8.
      matrix2(35,15)=(y1 - y4)/8.
      matrix2(35,16)=(z1 - z4)/8.
      matrix2(36,1)=0
      matrix2(36,2)=0
      matrix2(36,3)=0
      matrix2(36,4)=0
      matrix2(36,5)=0.5
      matrix2(36,6)=(-x2 + x3)/8.
      matrix2(36,7)=(-y2 + y3)/8.
      matrix2(36,8)=(-z2 + z3)/8.
      matrix2(36,9)=0.5
      matrix2(36,10)=(x2 - x3)/8.
      matrix2(36,11)=(y2 - y3)/8.
      matrix2(36,12)=(z2 - z3)/8.
      matrix2(36,13)=0
      matrix2(36,14)=0
      matrix2(36,15)=0
      matrix2(36,16)=0
      matrix2(37,1)=0
      matrix2(37,2)=0
      matrix2(37,3)=0
      matrix2(37,4)=0
      matrix2(37,5)=0.5
      matrix2(37,6)=(-x2 + x4)/8.
      matrix2(37,7)=(-y2 + y4)/8.
      matrix2(37,8)=(-z2 + z4)/8.
      matrix2(37,9)=0
      matrix2(37,10)=0
      matrix2(37,11)=0
      matrix2(37,12)=0
      matrix2(37,13)=0.5
      matrix2(37,14)=(x2 - x4)/8.
      matrix2(37,15)=(y2 - y4)/8.
      matrix2(37,16)=(z2 - z4)/8.
      matrix2(38,1)=0
      matrix2(38,2)=0
      matrix2(38,3)=0
      matrix2(38,4)=0
      matrix2(38,5)=0
      matrix2(38,6)=0
      matrix2(38,7)=0
      matrix2(38,8)=0
      matrix2(38,9)=0.5
      matrix2(38,10)=(-x3 + x4)/8.
      matrix2(38,11)=(-y3 + y4)/8.
      matrix2(38,12)=(-z3 + z4)/8.
      matrix2(38,13)=0.5
      matrix2(38,14)=(x3 - x4)/8.
      matrix2(38,15)=(y3 - y4)/8.
      matrix2(38,16)=(z3 - z4)/8.
      matrix2(39,1)=0
      matrix2(39,2)=0.5
      matrix2(39,3)=0
      matrix2(39,4)=0
      matrix2(39,5)=0
      matrix2(39,6)=0.5
      matrix2(39,7)=0
      matrix2(39,8)=0
      matrix2(39,9)=0
      matrix2(39,10)=0
      matrix2(39,11)=0
      matrix2(39,12)=0
      matrix2(39,13)=0
      matrix2(39,14)=0
      matrix2(39,15)=0
      matrix2(39,16)=0
      matrix2(40,1)=0
      matrix2(40,2)=0
      matrix2(40,3)=0.5
      matrix2(40,4)=0
      matrix2(40,5)=0
      matrix2(40,6)=0
      matrix2(40,7)=0.5
      matrix2(40,8)=0
      matrix2(40,9)=0
      matrix2(40,10)=0
      matrix2(40,11)=0
      matrix2(40,12)=0
      matrix2(40,13)=0
      matrix2(40,14)=0
      matrix2(40,15)=0
      matrix2(40,16)=0
      matrix2(41,1)=0
      matrix2(41,2)=0
      matrix2(41,3)=0
      matrix2(41,4)=0.5
      matrix2(41,5)=0
      matrix2(41,6)=0
      matrix2(41,7)=0
      matrix2(41,8)=0.5
      matrix2(41,9)=0
      matrix2(41,10)=0
      matrix2(41,11)=0
      matrix2(41,12)=0
      matrix2(41,13)=0
      matrix2(41,14)=0
      matrix2(41,15)=0
      matrix2(41,16)=0
      matrix2(42,1)=0
      matrix2(42,2)=0.5
      matrix2(42,3)=0
      matrix2(42,4)=0
      matrix2(42,5)=0
      matrix2(42,6)=0
      matrix2(42,7)=0
      matrix2(42,8)=0
      matrix2(42,9)=0
      matrix2(42,10)=0.5
      matrix2(42,11)=0
      matrix2(42,12)=0
      matrix2(42,13)=0
      matrix2(42,14)=0
      matrix2(42,15)=0
      matrix2(42,16)=0
      matrix2(43,1)=0
      matrix2(43,2)=0
      matrix2(43,3)=0.5
      matrix2(43,4)=0
      matrix2(43,5)=0
      matrix2(43,6)=0
      matrix2(43,7)=0
      matrix2(43,8)=0
      matrix2(43,9)=0
      matrix2(43,10)=0
      matrix2(43,11)=0.5
      matrix2(43,12)=0
      matrix2(43,13)=0
      matrix2(43,14)=0
      matrix2(43,15)=0
      matrix2(43,16)=0
      matrix2(44,1)=0
      matrix2(44,2)=0
      matrix2(44,3)=0
      matrix2(44,4)=0.5
      matrix2(44,5)=0
      matrix2(44,6)=0
      matrix2(44,7)=0
      matrix2(44,8)=0
      matrix2(44,9)=0
      matrix2(44,10)=0
      matrix2(44,11)=0
      matrix2(44,12)=0.5
      matrix2(44,13)=0
      matrix2(44,14)=0
      matrix2(44,15)=0
      matrix2(44,16)=0
      matrix2(45,1)=0
      matrix2(45,2)=0.5
      matrix2(45,3)=0
      matrix2(45,4)=0
      matrix2(45,5)=0
      matrix2(45,6)=0
      matrix2(45,7)=0
      matrix2(45,8)=0
      matrix2(45,9)=0
      matrix2(45,10)=0
      matrix2(45,11)=0
      matrix2(45,12)=0
      matrix2(45,13)=0
      matrix2(45,14)=0.5
      matrix2(45,15)=0
      matrix2(45,16)=0
      matrix2(46,1)=0
      matrix2(46,2)=0
      matrix2(46,3)=0.5
      matrix2(46,4)=0
      matrix2(46,5)=0
      matrix2(46,6)=0
      matrix2(46,7)=0
      matrix2(46,8)=0
      matrix2(46,9)=0
      matrix2(46,10)=0
      matrix2(46,11)=0
      matrix2(46,12)=0
      matrix2(46,13)=0
      matrix2(46,14)=0
      matrix2(46,15)=0.5
      matrix2(46,16)=0
      matrix2(47,1)=0
      matrix2(47,2)=0
      matrix2(47,3)=0
      matrix2(47,4)=0.5
      matrix2(47,5)=0
      matrix2(47,6)=0
      matrix2(47,7)=0
      matrix2(47,8)=0
      matrix2(47,9)=0
      matrix2(47,10)=0
      matrix2(47,11)=0
      matrix2(47,12)=0
      matrix2(47,13)=0
      matrix2(47,14)=0
      matrix2(47,15)=0
      matrix2(47,16)=0.5
      matrix2(48,1)=0
      matrix2(48,2)=0
      matrix2(48,3)=0
      matrix2(48,4)=0
      matrix2(48,5)=0
      matrix2(48,6)=0.5
      matrix2(48,7)=0
      matrix2(48,8)=0
      matrix2(48,9)=0
      matrix2(48,10)=0.5
      matrix2(48,11)=0
      matrix2(48,12)=0
      matrix2(48,13)=0
      matrix2(48,14)=0
      matrix2(48,15)=0
      matrix2(48,16)=0
      matrix2(49,1)=0
      matrix2(49,2)=0
      matrix2(49,3)=0
      matrix2(49,4)=0
      matrix2(49,5)=0
      matrix2(49,6)=0
      matrix2(49,7)=0.5
      matrix2(49,8)=0
      matrix2(49,9)=0
      matrix2(49,10)=0
      matrix2(49,11)=0.5
      matrix2(49,12)=0
      matrix2(49,13)=0
      matrix2(49,14)=0
      matrix2(49,15)=0
      matrix2(49,16)=0
      matrix2(50,1)=0
      matrix2(50,2)=0
      matrix2(50,3)=0
      matrix2(50,4)=0
      matrix2(50,5)=0
      matrix2(50,6)=0
      matrix2(50,7)=0
      matrix2(50,8)=0.5
      matrix2(50,9)=0
      matrix2(50,10)=0
      matrix2(50,11)=0
      matrix2(50,12)=0.5
      matrix2(50,13)=0
      matrix2(50,14)=0
      matrix2(50,15)=0
      matrix2(50,16)=0
      matrix2(51,1)=0
      matrix2(51,2)=0
      matrix2(51,3)=0
      matrix2(51,4)=0
      matrix2(51,5)=0
      matrix2(51,6)=0.5
      matrix2(51,7)=0
      matrix2(51,8)=0
      matrix2(51,9)=0
      matrix2(51,10)=0
      matrix2(51,11)=0
      matrix2(51,12)=0
      matrix2(51,13)=0
      matrix2(51,14)=0.5
      matrix2(51,15)=0
      matrix2(51,16)=0
      matrix2(52,1)=0
      matrix2(52,2)=0
      matrix2(52,3)=0
      matrix2(52,4)=0
      matrix2(52,5)=0
      matrix2(52,6)=0
      matrix2(52,7)=0.5
      matrix2(52,8)=0
      matrix2(52,9)=0
      matrix2(52,10)=0
      matrix2(52,11)=0
      matrix2(52,12)=0
      matrix2(52,13)=0
      matrix2(52,14)=0
      matrix2(52,15)=0.5
      matrix2(52,16)=0
      matrix2(53,1)=0
      matrix2(53,2)=0
      matrix2(53,3)=0
      matrix2(53,4)=0
      matrix2(53,5)=0
      matrix2(53,6)=0
      matrix2(53,7)=0
      matrix2(53,8)=0.5
      matrix2(53,9)=0
      matrix2(53,10)=0
      matrix2(53,11)=0
      matrix2(53,12)=0
      matrix2(53,13)=0
      matrix2(53,14)=0
      matrix2(53,15)=0
      matrix2(53,16)=0.5
      matrix2(54,1)=0
      matrix2(54,2)=0
      matrix2(54,3)=0
      matrix2(54,4)=0
      matrix2(54,5)=0
      matrix2(54,6)=0
      matrix2(54,7)=0
      matrix2(54,8)=0
      matrix2(54,9)=0
      matrix2(54,10)=0.5
      matrix2(54,11)=0
      matrix2(54,12)=0
      matrix2(54,13)=0
      matrix2(54,14)=0.5
      matrix2(54,15)=0
      matrix2(54,16)=0
      matrix2(55,1)=0
      matrix2(55,2)=0
      matrix2(55,3)=0
      matrix2(55,4)=0
      matrix2(55,5)=0
      matrix2(55,6)=0
      matrix2(55,7)=0
      matrix2(55,8)=0
      matrix2(55,9)=0
      matrix2(55,10)=0
      matrix2(55,11)=0.5
      matrix2(55,12)=0
      matrix2(55,13)=0
      matrix2(55,14)=0
      matrix2(55,15)=0.5
      matrix2(55,16)=0
      matrix2(56,1)=0
      matrix2(56,2)=0
      matrix2(56,3)=0
      matrix2(56,4)=0
      matrix2(56,5)=0
      matrix2(56,6)=0
      matrix2(56,7)=0
      matrix2(56,8)=0
      matrix2(56,9)=0
      matrix2(56,10)=0
      matrix2(56,11)=0
      matrix2(56,12)=0.5
      matrix2(56,13)=0
      matrix2(56,14)=0
      matrix2(56,15)=0
      matrix2(56,16)=0.5

! the q-matrix maps from the 56 polynomial coeffs. to the 56 pseudo-dof
      q(1,1)=1
      q(1,2)=z1
      q(1,3)=z1**2
      q(1,4)=z1**3
      q(1,5)=z1**4
      q(1,6)=z1**5
      q(1,7)=y1
      q(1,8)=y1*z1
      q(1,9)=y1*z1**2
      q(1,10)=y1*z1**3
      q(1,11)=y1*z1**4
      q(1,12)=y1**2
      q(1,13)=y1**2*z1
      q(1,14)=y1**2*z1**2
      q(1,15)=y1**2*z1**3
      q(1,16)=y1**3
      q(1,17)=y1**3*z1
      q(1,18)=y1**3*z1**2
      q(1,19)=y1**4
      q(1,20)=y1**4*z1
      q(1,21)=y1**5
      q(1,22)=x1
      q(1,23)=x1*z1
      q(1,24)=x1*z1**2
      q(1,25)=x1*z1**3
      q(1,26)=x1*z1**4
      q(1,27)=x1*y1
      q(1,28)=x1*y1*z1
      q(1,29)=x1*y1*z1**2
      q(1,30)=x1*y1*z1**3
      q(1,31)=x1*y1**2
      q(1,32)=x1*y1**2*z1
      q(1,33)=x1*y1**2*z1**2
      q(1,34)=x1*y1**3
      q(1,35)=x1*y1**3*z1
      q(1,36)=x1*y1**4
      q(1,37)=x1**2
      q(1,38)=x1**2*z1
      q(1,39)=x1**2*z1**2
      q(1,40)=x1**2*z1**3
      q(1,41)=x1**2*y1
      q(1,42)=x1**2*y1*z1
      q(1,43)=x1**2*y1*z1**2
      q(1,44)=x1**2*y1**2
      q(1,45)=x1**2*y1**2*z1
      q(1,46)=x1**2*y1**3
      q(1,47)=x1**3
      q(1,48)=x1**3*z1
      q(1,49)=x1**3*z1**2
      q(1,50)=x1**3*y1
      q(1,51)=x1**3*y1*z1
      q(1,52)=x1**3*y1**2
      q(1,53)=x1**4
      q(1,54)=x1**4*z1
      q(1,55)=x1**4*y1
      q(1,56)=x1**5
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
      q(2,11)=0
      q(2,12)=0
      q(2,13)=0
      q(2,14)=0
      q(2,15)=0
      q(2,16)=0
      q(2,17)=0
      q(2,18)=0
      q(2,19)=0
      q(2,20)=0
      q(2,21)=0
      q(2,22)=1
      q(2,23)=z1
      q(2,24)=z1**2
      q(2,25)=z1**3
      q(2,26)=z1**4
      q(2,27)=y1
      q(2,28)=y1*z1
      q(2,29)=y1*z1**2
      q(2,30)=y1*z1**3
      q(2,31)=y1**2
      q(2,32)=y1**2*z1
      q(2,33)=y1**2*z1**2
      q(2,34)=y1**3
      q(2,35)=y1**3*z1
      q(2,36)=y1**4
      q(2,37)=2*x1
      q(2,38)=2*x1*z1
      q(2,39)=2*x1*z1**2
      q(2,40)=2*x1*z1**3
      q(2,41)=2*x1*y1
      q(2,42)=2*x1*y1*z1
      q(2,43)=2*x1*y1*z1**2
      q(2,44)=2*x1*y1**2
      q(2,45)=2*x1*y1**2*z1
      q(2,46)=2*x1*y1**3
      q(2,47)=3*x1**2
      q(2,48)=3*x1**2*z1
      q(2,49)=3*x1**2*z1**2
      q(2,50)=3*x1**2*y1
      q(2,51)=3*x1**2*y1*z1
      q(2,52)=3*x1**2*y1**2
      q(2,53)=4*x1**3
      q(2,54)=4*x1**3*z1
      q(2,55)=4*x1**3*y1
      q(2,56)=5*x1**4
      q(3,1)=0
      q(3,2)=0
      q(3,3)=0
      q(3,4)=0
      q(3,5)=0
      q(3,6)=0
      q(3,7)=1
      q(3,8)=z1
      q(3,9)=z1**2
      q(3,10)=z1**3
      q(3,11)=z1**4
      q(3,12)=2*y1
      q(3,13)=2*y1*z1
      q(3,14)=2*y1*z1**2
      q(3,15)=2*y1*z1**3
      q(3,16)=3*y1**2
      q(3,17)=3*y1**2*z1
      q(3,18)=3*y1**2*z1**2
      q(3,19)=4*y1**3
      q(3,20)=4*y1**3*z1
      q(3,21)=5*y1**4
      q(3,22)=0
      q(3,23)=0
      q(3,24)=0
      q(3,25)=0
      q(3,26)=0
      q(3,27)=x1
      q(3,28)=x1*z1
      q(3,29)=x1*z1**2
      q(3,30)=x1*z1**3
      q(3,31)=2*x1*y1
      q(3,32)=2*x1*y1*z1
      q(3,33)=2*x1*y1*z1**2
      q(3,34)=3*x1*y1**2
      q(3,35)=3*x1*y1**2*z1
      q(3,36)=4*x1*y1**3
      q(3,37)=0
      q(3,38)=0
      q(3,39)=0
      q(3,40)=0
      q(3,41)=x1**2
      q(3,42)=x1**2*z1
      q(3,43)=x1**2*z1**2
      q(3,44)=2*x1**2*y1
      q(3,45)=2*x1**2*y1*z1
      q(3,46)=3*x1**2*y1**2
      q(3,47)=0
      q(3,48)=0
      q(3,49)=0
      q(3,50)=x1**3
      q(3,51)=x1**3*z1
      q(3,52)=2*x1**3*y1
      q(3,53)=0
      q(3,54)=0
      q(3,55)=x1**4
      q(3,56)=0
      q(4,1)=0
      q(4,2)=1
      q(4,3)=2*z1
      q(4,4)=3*z1**2
      q(4,5)=4*z1**3
      q(4,6)=5*z1**4
      q(4,7)=0
      q(4,8)=y1
      q(4,9)=2*y1*z1
      q(4,10)=3*y1*z1**2
      q(4,11)=4*y1*z1**3
      q(4,12)=0
      q(4,13)=y1**2
      q(4,14)=2*y1**2*z1
      q(4,15)=3*y1**2*z1**2
      q(4,16)=0
      q(4,17)=y1**3
      q(4,18)=2*y1**3*z1
      q(4,19)=0
      q(4,20)=y1**4
      q(4,21)=0
      q(4,22)=0
      q(4,23)=x1
      q(4,24)=2*x1*z1
      q(4,25)=3*x1*z1**2
      q(4,26)=4*x1*z1**3
      q(4,27)=0
      q(4,28)=x1*y1
      q(4,29)=2*x1*y1*z1
      q(4,30)=3*x1*y1*z1**2
      q(4,31)=0
      q(4,32)=x1*y1**2
      q(4,33)=2*x1*y1**2*z1
      q(4,34)=0
      q(4,35)=x1*y1**3
      q(4,36)=0
      q(4,37)=0
      q(4,38)=x1**2
      q(4,39)=2*x1**2*z1
      q(4,40)=3*x1**2*z1**2
      q(4,41)=0
      q(4,42)=x1**2*y1
      q(4,43)=2*x1**2*y1*z1
      q(4,44)=0
      q(4,45)=x1**2*y1**2
      q(4,46)=0
      q(4,47)=0
      q(4,48)=x1**3
      q(4,49)=2*x1**3*z1
      q(4,50)=0
      q(4,51)=x1**3*y1
      q(4,52)=0
      q(4,53)=0
      q(4,54)=x1**4
      q(4,55)=0
      q(4,56)=0
      q(5,1)=1
      q(5,2)=z2
      q(5,3)=z2**2
      q(5,4)=z2**3
      q(5,5)=z2**4
      q(5,6)=z2**5
      q(5,7)=y2
      q(5,8)=y2*z2
      q(5,9)=y2*z2**2
      q(5,10)=y2*z2**3
      q(5,11)=y2*z2**4
      q(5,12)=y2**2
      q(5,13)=y2**2*z2
      q(5,14)=y2**2*z2**2
      q(5,15)=y2**2*z2**3
      q(5,16)=y2**3
      q(5,17)=y2**3*z2
      q(5,18)=y2**3*z2**2
      q(5,19)=y2**4
      q(5,20)=y2**4*z2
      q(5,21)=y2**5
      q(5,22)=x2
      q(5,23)=x2*z2
      q(5,24)=x2*z2**2
      q(5,25)=x2*z2**3
      q(5,26)=x2*z2**4
      q(5,27)=x2*y2
      q(5,28)=x2*y2*z2
      q(5,29)=x2*y2*z2**2
      q(5,30)=x2*y2*z2**3
      q(5,31)=x2*y2**2
      q(5,32)=x2*y2**2*z2
      q(5,33)=x2*y2**2*z2**2
      q(5,34)=x2*y2**3
      q(5,35)=x2*y2**3*z2
      q(5,36)=x2*y2**4
      q(5,37)=x2**2
      q(5,38)=x2**2*z2
      q(5,39)=x2**2*z2**2
      q(5,40)=x2**2*z2**3
      q(5,41)=x2**2*y2
      q(5,42)=x2**2*y2*z2
      q(5,43)=x2**2*y2*z2**2
      q(5,44)=x2**2*y2**2
      q(5,45)=x2**2*y2**2*z2
      q(5,46)=x2**2*y2**3
      q(5,47)=x2**3
      q(5,48)=x2**3*z2
      q(5,49)=x2**3*z2**2
      q(5,50)=x2**3*y2
      q(5,51)=x2**3*y2*z2
      q(5,52)=x2**3*y2**2
      q(5,53)=x2**4
      q(5,54)=x2**4*z2
      q(5,55)=x2**4*y2
      q(5,56)=x2**5
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
      q(6,11)=0
      q(6,12)=0
      q(6,13)=0
      q(6,14)=0
      q(6,15)=0
      q(6,16)=0
      q(6,17)=0
      q(6,18)=0
      q(6,19)=0
      q(6,20)=0
      q(6,21)=0
      q(6,22)=1
      q(6,23)=z2
      q(6,24)=z2**2
      q(6,25)=z2**3
      q(6,26)=z2**4
      q(6,27)=y2
      q(6,28)=y2*z2
      q(6,29)=y2*z2**2
      q(6,30)=y2*z2**3
      q(6,31)=y2**2
      q(6,32)=y2**2*z2
      q(6,33)=y2**2*z2**2
      q(6,34)=y2**3
      q(6,35)=y2**3*z2
      q(6,36)=y2**4
      q(6,37)=2*x2
      q(6,38)=2*x2*z2
      q(6,39)=2*x2*z2**2
      q(6,40)=2*x2*z2**3
      q(6,41)=2*x2*y2
      q(6,42)=2*x2*y2*z2
      q(6,43)=2*x2*y2*z2**2
      q(6,44)=2*x2*y2**2
      q(6,45)=2*x2*y2**2*z2
      q(6,46)=2*x2*y2**3
      q(6,47)=3*x2**2
      q(6,48)=3*x2**2*z2
      q(6,49)=3*x2**2*z2**2
      q(6,50)=3*x2**2*y2
      q(6,51)=3*x2**2*y2*z2
      q(6,52)=3*x2**2*y2**2
      q(6,53)=4*x2**3
      q(6,54)=4*x2**3*z2
      q(6,55)=4*x2**3*y2
      q(6,56)=5*x2**4
      q(7,1)=0
      q(7,2)=0
      q(7,3)=0
      q(7,4)=0
      q(7,5)=0
      q(7,6)=0
      q(7,7)=1
      q(7,8)=z2
      q(7,9)=z2**2
      q(7,10)=z2**3
      q(7,11)=z2**4
      q(7,12)=2*y2
      q(7,13)=2*y2*z2
      q(7,14)=2*y2*z2**2
      q(7,15)=2*y2*z2**3
      q(7,16)=3*y2**2
      q(7,17)=3*y2**2*z2
      q(7,18)=3*y2**2*z2**2
      q(7,19)=4*y2**3
      q(7,20)=4*y2**3*z2
      q(7,21)=5*y2**4
      q(7,22)=0
      q(7,23)=0
      q(7,24)=0
      q(7,25)=0
      q(7,26)=0
      q(7,27)=x2
      q(7,28)=x2*z2
      q(7,29)=x2*z2**2
      q(7,30)=x2*z2**3
      q(7,31)=2*x2*y2
      q(7,32)=2*x2*y2*z2
      q(7,33)=2*x2*y2*z2**2
      q(7,34)=3*x2*y2**2
      q(7,35)=3*x2*y2**2*z2
      q(7,36)=4*x2*y2**3
      q(7,37)=0
      q(7,38)=0
      q(7,39)=0
      q(7,40)=0
      q(7,41)=x2**2
      q(7,42)=x2**2*z2
      q(7,43)=x2**2*z2**2
      q(7,44)=2*x2**2*y2
      q(7,45)=2*x2**2*y2*z2
      q(7,46)=3*x2**2*y2**2
      q(7,47)=0
      q(7,48)=0
      q(7,49)=0
      q(7,50)=x2**3
      q(7,51)=x2**3*z2
      q(7,52)=2*x2**3*y2
      q(7,53)=0
      q(7,54)=0
      q(7,55)=x2**4
      q(7,56)=0
      q(8,1)=0
      q(8,2)=1
      q(8,3)=2*z2
      q(8,4)=3*z2**2
      q(8,5)=4*z2**3
      q(8,6)=5*z2**4
      q(8,7)=0
      q(8,8)=y2
      q(8,9)=2*y2*z2
      q(8,10)=3*y2*z2**2
      q(8,11)=4*y2*z2**3
      q(8,12)=0
      q(8,13)=y2**2
      q(8,14)=2*y2**2*z2
      q(8,15)=3*y2**2*z2**2
      q(8,16)=0
      q(8,17)=y2**3
      q(8,18)=2*y2**3*z2
      q(8,19)=0
      q(8,20)=y2**4
      q(8,21)=0
      q(8,22)=0
      q(8,23)=x2
      q(8,24)=2*x2*z2
      q(8,25)=3*x2*z2**2
      q(8,26)=4*x2*z2**3
      q(8,27)=0
      q(8,28)=x2*y2
      q(8,29)=2*x2*y2*z2
      q(8,30)=3*x2*y2*z2**2
      q(8,31)=0
      q(8,32)=x2*y2**2
      q(8,33)=2*x2*y2**2*z2
      q(8,34)=0
      q(8,35)=x2*y2**3
      q(8,36)=0
      q(8,37)=0
      q(8,38)=x2**2
      q(8,39)=2*x2**2*z2
      q(8,40)=3*x2**2*z2**2
      q(8,41)=0
      q(8,42)=x2**2*y2
      q(8,43)=2*x2**2*y2*z2
      q(8,44)=0
      q(8,45)=x2**2*y2**2
      q(8,46)=0
      q(8,47)=0
      q(8,48)=x2**3
      q(8,49)=2*x2**3*z2
      q(8,50)=0
      q(8,51)=x2**3*y2
      q(8,52)=0
      q(8,53)=0
      q(8,54)=x2**4
      q(8,55)=0
      q(8,56)=0
      q(9,1)=1
      q(9,2)=z3
      q(9,3)=z3**2
      q(9,4)=z3**3
      q(9,5)=z3**4
      q(9,6)=z3**5
      q(9,7)=y3
      q(9,8)=y3*z3
      q(9,9)=y3*z3**2
      q(9,10)=y3*z3**3
      q(9,11)=y3*z3**4
      q(9,12)=y3**2
      q(9,13)=y3**2*z3
      q(9,14)=y3**2*z3**2
      q(9,15)=y3**2*z3**3
      q(9,16)=y3**3
      q(9,17)=y3**3*z3
      q(9,18)=y3**3*z3**2
      q(9,19)=y3**4
      q(9,20)=y3**4*z3
      q(9,21)=y3**5
      q(9,22)=x3
      q(9,23)=x3*z3
      q(9,24)=x3*z3**2
      q(9,25)=x3*z3**3
      q(9,26)=x3*z3**4
      q(9,27)=x3*y3
      q(9,28)=x3*y3*z3
      q(9,29)=x3*y3*z3**2
      q(9,30)=x3*y3*z3**3
      q(9,31)=x3*y3**2
      q(9,32)=x3*y3**2*z3
      q(9,33)=x3*y3**2*z3**2
      q(9,34)=x3*y3**3
      q(9,35)=x3*y3**3*z3
      q(9,36)=x3*y3**4
      q(9,37)=x3**2
      q(9,38)=x3**2*z3
      q(9,39)=x3**2*z3**2
      q(9,40)=x3**2*z3**3
      q(9,41)=x3**2*y3
      q(9,42)=x3**2*y3*z3
      q(9,43)=x3**2*y3*z3**2
      q(9,44)=x3**2*y3**2
      q(9,45)=x3**2*y3**2*z3
      q(9,46)=x3**2*y3**3
      q(9,47)=x3**3
      q(9,48)=x3**3*z3
      q(9,49)=x3**3*z3**2
      q(9,50)=x3**3*y3
      q(9,51)=x3**3*y3*z3
      q(9,52)=x3**3*y3**2
      q(9,53)=x3**4
      q(9,54)=x3**4*z3
      q(9,55)=x3**4*y3
      q(9,56)=x3**5
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
      q(10,11)=0
      q(10,12)=0
      q(10,13)=0
      q(10,14)=0
      q(10,15)=0
      q(10,16)=0
      q(10,17)=0
      q(10,18)=0
      q(10,19)=0
      q(10,20)=0
      q(10,21)=0
      q(10,22)=1
      q(10,23)=z3
      q(10,24)=z3**2
      q(10,25)=z3**3
      q(10,26)=z3**4
      q(10,27)=y3
      q(10,28)=y3*z3
      q(10,29)=y3*z3**2
      q(10,30)=y3*z3**3
      q(10,31)=y3**2
      q(10,32)=y3**2*z3
      q(10,33)=y3**2*z3**2
      q(10,34)=y3**3
      q(10,35)=y3**3*z3
      q(10,36)=y3**4
      q(10,37)=2*x3
      q(10,38)=2*x3*z3
      q(10,39)=2*x3*z3**2
      q(10,40)=2*x3*z3**3
      q(10,41)=2*x3*y3
      q(10,42)=2*x3*y3*z3
      q(10,43)=2*x3*y3*z3**2
      q(10,44)=2*x3*y3**2
      q(10,45)=2*x3*y3**2*z3
      q(10,46)=2*x3*y3**3
      q(10,47)=3*x3**2
      q(10,48)=3*x3**2*z3
      q(10,49)=3*x3**2*z3**2
      q(10,50)=3*x3**2*y3
      q(10,51)=3*x3**2*y3*z3
      q(10,52)=3*x3**2*y3**2
      q(10,53)=4*x3**3
      q(10,54)=4*x3**3*z3
      q(10,55)=4*x3**3*y3
      q(10,56)=5*x3**4
      q(11,1)=0
      q(11,2)=0
      q(11,3)=0
      q(11,4)=0
      q(11,5)=0
      q(11,6)=0
      q(11,7)=1
      q(11,8)=z3
      q(11,9)=z3**2
      q(11,10)=z3**3
      q(11,11)=z3**4
      q(11,12)=2*y3
      q(11,13)=2*y3*z3
      q(11,14)=2*y3*z3**2
      q(11,15)=2*y3*z3**3
      q(11,16)=3*y3**2
      q(11,17)=3*y3**2*z3
      q(11,18)=3*y3**2*z3**2
      q(11,19)=4*y3**3
      q(11,20)=4*y3**3*z3
      q(11,21)=5*y3**4
      q(11,22)=0
      q(11,23)=0
      q(11,24)=0
      q(11,25)=0
      q(11,26)=0
      q(11,27)=x3
      q(11,28)=x3*z3
      q(11,29)=x3*z3**2
      q(11,30)=x3*z3**3
      q(11,31)=2*x3*y3
      q(11,32)=2*x3*y3*z3
      q(11,33)=2*x3*y3*z3**2
      q(11,34)=3*x3*y3**2
      q(11,35)=3*x3*y3**2*z3
      q(11,36)=4*x3*y3**3
      q(11,37)=0
      q(11,38)=0
      q(11,39)=0
      q(11,40)=0
      q(11,41)=x3**2
      q(11,42)=x3**2*z3
      q(11,43)=x3**2*z3**2
      q(11,44)=2*x3**2*y3
      q(11,45)=2*x3**2*y3*z3
      q(11,46)=3*x3**2*y3**2
      q(11,47)=0
      q(11,48)=0
      q(11,49)=0
      q(11,50)=x3**3
      q(11,51)=x3**3*z3
      q(11,52)=2*x3**3*y3
      q(11,53)=0
      q(11,54)=0
      q(11,55)=x3**4
      q(11,56)=0
      q(12,1)=0
      q(12,2)=1
      q(12,3)=2*z3
      q(12,4)=3*z3**2
      q(12,5)=4*z3**3
      q(12,6)=5*z3**4
      q(12,7)=0
      q(12,8)=y3
      q(12,9)=2*y3*z3
      q(12,10)=3*y3*z3**2
      q(12,11)=4*y3*z3**3
      q(12,12)=0
      q(12,13)=y3**2
      q(12,14)=2*y3**2*z3
      q(12,15)=3*y3**2*z3**2
      q(12,16)=0
      q(12,17)=y3**3
      q(12,18)=2*y3**3*z3
      q(12,19)=0
      q(12,20)=y3**4
      q(12,21)=0
      q(12,22)=0
      q(12,23)=x3
      q(12,24)=2*x3*z3
      q(12,25)=3*x3*z3**2
      q(12,26)=4*x3*z3**3
      q(12,27)=0
      q(12,28)=x3*y3
      q(12,29)=2*x3*y3*z3
      q(12,30)=3*x3*y3*z3**2
      q(12,31)=0
      q(12,32)=x3*y3**2
      q(12,33)=2*x3*y3**2*z3
      q(12,34)=0
      q(12,35)=x3*y3**3
      q(12,36)=0
      q(12,37)=0
      q(12,38)=x3**2
      q(12,39)=2*x3**2*z3
      q(12,40)=3*x3**2*z3**2
      q(12,41)=0
      q(12,42)=x3**2*y3
      q(12,43)=2*x3**2*y3*z3
      q(12,44)=0
      q(12,45)=x3**2*y3**2
      q(12,46)=0
      q(12,47)=0
      q(12,48)=x3**3
      q(12,49)=2*x3**3*z3
      q(12,50)=0
      q(12,51)=x3**3*y3
      q(12,52)=0
      q(12,53)=0
      q(12,54)=x3**4
      q(12,55)=0
      q(12,56)=0
      q(13,1)=1
      q(13,2)=z4
      q(13,3)=z4**2
      q(13,4)=z4**3
      q(13,5)=z4**4
      q(13,6)=z4**5
      q(13,7)=y4
      q(13,8)=y4*z4
      q(13,9)=y4*z4**2
      q(13,10)=y4*z4**3
      q(13,11)=y4*z4**4
      q(13,12)=y4**2
      q(13,13)=y4**2*z4
      q(13,14)=y4**2*z4**2
      q(13,15)=y4**2*z4**3
      q(13,16)=y4**3
      q(13,17)=y4**3*z4
      q(13,18)=y4**3*z4**2
      q(13,19)=y4**4
      q(13,20)=y4**4*z4
      q(13,21)=y4**5
      q(13,22)=x4
      q(13,23)=x4*z4
      q(13,24)=x4*z4**2
      q(13,25)=x4*z4**3
      q(13,26)=x4*z4**4
      q(13,27)=x4*y4
      q(13,28)=x4*y4*z4
      q(13,29)=x4*y4*z4**2
      q(13,30)=x4*y4*z4**3
      q(13,31)=x4*y4**2
      q(13,32)=x4*y4**2*z4
      q(13,33)=x4*y4**2*z4**2
      q(13,34)=x4*y4**3
      q(13,35)=x4*y4**3*z4
      q(13,36)=x4*y4**4
      q(13,37)=x4**2
      q(13,38)=x4**2*z4
      q(13,39)=x4**2*z4**2
      q(13,40)=x4**2*z4**3
      q(13,41)=x4**2*y4
      q(13,42)=x4**2*y4*z4
      q(13,43)=x4**2*y4*z4**2
      q(13,44)=x4**2*y4**2
      q(13,45)=x4**2*y4**2*z4
      q(13,46)=x4**2*y4**3
      q(13,47)=x4**3
      q(13,48)=x4**3*z4
      q(13,49)=x4**3*z4**2
      q(13,50)=x4**3*y4
      q(13,51)=x4**3*y4*z4
      q(13,52)=x4**3*y4**2
      q(13,53)=x4**4
      q(13,54)=x4**4*z4
      q(13,55)=x4**4*y4
      q(13,56)=x4**5
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
      q(14,11)=0
      q(14,12)=0
      q(14,13)=0
      q(14,14)=0
      q(14,15)=0
      q(14,16)=0
      q(14,17)=0
      q(14,18)=0
      q(14,19)=0
      q(14,20)=0
      q(14,21)=0
      q(14,22)=1
      q(14,23)=z4
      q(14,24)=z4**2
      q(14,25)=z4**3
      q(14,26)=z4**4
      q(14,27)=y4
      q(14,28)=y4*z4
      q(14,29)=y4*z4**2
      q(14,30)=y4*z4**3
      q(14,31)=y4**2
      q(14,32)=y4**2*z4
      q(14,33)=y4**2*z4**2
      q(14,34)=y4**3
      q(14,35)=y4**3*z4
      q(14,36)=y4**4
      q(14,37)=2*x4
      q(14,38)=2*x4*z4
      q(14,39)=2*x4*z4**2
      q(14,40)=2*x4*z4**3
      q(14,41)=2*x4*y4
      q(14,42)=2*x4*y4*z4
      q(14,43)=2*x4*y4*z4**2
      q(14,44)=2*x4*y4**2
      q(14,45)=2*x4*y4**2*z4
      q(14,46)=2*x4*y4**3
      q(14,47)=3*x4**2
      q(14,48)=3*x4**2*z4
      q(14,49)=3*x4**2*z4**2
      q(14,50)=3*x4**2*y4
      q(14,51)=3*x4**2*y4*z4
      q(14,52)=3*x4**2*y4**2
      q(14,53)=4*x4**3
      q(14,54)=4*x4**3*z4
      q(14,55)=4*x4**3*y4
      q(14,56)=5*x4**4
      q(15,1)=0
      q(15,2)=0
      q(15,3)=0
      q(15,4)=0
      q(15,5)=0
      q(15,6)=0
      q(15,7)=1
      q(15,8)=z4
      q(15,9)=z4**2
      q(15,10)=z4**3
      q(15,11)=z4**4
      q(15,12)=2*y4
      q(15,13)=2*y4*z4
      q(15,14)=2*y4*z4**2
      q(15,15)=2*y4*z4**3
      q(15,16)=3*y4**2
      q(15,17)=3*y4**2*z4
      q(15,18)=3*y4**2*z4**2
      q(15,19)=4*y4**3
      q(15,20)=4*y4**3*z4
      q(15,21)=5*y4**4
      q(15,22)=0
      q(15,23)=0
      q(15,24)=0
      q(15,25)=0
      q(15,26)=0
      q(15,27)=x4
      q(15,28)=x4*z4
      q(15,29)=x4*z4**2
      q(15,30)=x4*z4**3
      q(15,31)=2*x4*y4
      q(15,32)=2*x4*y4*z4
      q(15,33)=2*x4*y4*z4**2
      q(15,34)=3*x4*y4**2
      q(15,35)=3*x4*y4**2*z4
      q(15,36)=4*x4*y4**3
      q(15,37)=0
      q(15,38)=0
      q(15,39)=0
      q(15,40)=0
      q(15,41)=x4**2
      q(15,42)=x4**2*z4
      q(15,43)=x4**2*z4**2
      q(15,44)=2*x4**2*y4
      q(15,45)=2*x4**2*y4*z4
      q(15,46)=3*x4**2*y4**2
      q(15,47)=0
      q(15,48)=0
      q(15,49)=0
      q(15,50)=x4**3
      q(15,51)=x4**3*z4
      q(15,52)=2*x4**3*y4
      q(15,53)=0
      q(15,54)=0
      q(15,55)=x4**4
      q(15,56)=0
      q(16,1)=0
      q(16,2)=1
      q(16,3)=2*z4
      q(16,4)=3*z4**2
      q(16,5)=4*z4**3
      q(16,6)=5*z4**4
      q(16,7)=0
      q(16,8)=y4
      q(16,9)=2*y4*z4
      q(16,10)=3*y4*z4**2
      q(16,11)=4*y4*z4**3
      q(16,12)=0
      q(16,13)=y4**2
      q(16,14)=2*y4**2*z4
      q(16,15)=3*y4**2*z4**2
      q(16,16)=0
      q(16,17)=y4**3
      q(16,18)=2*y4**3*z4
      q(16,19)=0
      q(16,20)=y4**4
      q(16,21)=0
      q(16,22)=0
      q(16,23)=x4
      q(16,24)=2*x4*z4
      q(16,25)=3*x4*z4**2
      q(16,26)=4*x4*z4**3
      q(16,27)=0
      q(16,28)=x4*y4
      q(16,29)=2*x4*y4*z4
      q(16,30)=3*x4*y4*z4**2
      q(16,31)=0
      q(16,32)=x4*y4**2
      q(16,33)=2*x4*y4**2*z4
      q(16,34)=0
      q(16,35)=x4*y4**3
      q(16,36)=0
      q(16,37)=0
      q(16,38)=x4**2
      q(16,39)=2*x4**2*z4
      q(16,40)=3*x4**2*z4**2
      q(16,41)=0
      q(16,42)=x4**2*y4
      q(16,43)=2*x4**2*y4*z4
      q(16,44)=0
      q(16,45)=x4**2*y4**2
      q(16,46)=0
      q(16,47)=0
      q(16,48)=x4**3
      q(16,49)=2*x4**3*z4
      q(16,50)=0
      q(16,51)=x4**3*y4
      q(16,52)=0
      q(16,53)=0
      q(16,54)=x4**4
      q(16,55)=0
      q(16,56)=0
      q(17,1)=1
      q(17,2)=(z1 + z2 + z4)/3.
      q(17,3)=(z1 + z2 + z4)**2/9.
      q(17,4)=(z1 + z2 + z4)**3/27.
      q(17,5)=(z1 + z2 + z4)**4/81.
      q(17,6)=(z1 + z2 + z4)**5/243.
      q(17,7)=(y1 + y2 + y4)/3.
      q(17,8)=((y1 + y2 + y4)*(z1 + z2 + z4))/9.
      q(17,9)=((y1 + y2 + y4)*(z1 + z2 + z4)**2)/27.
      q(17,10)=((y1 + y2 + y4)*(z1 + z2 + z4)**3)/81.
      q(17,11)=((y1 + y2 + y4)*(z1 + z2 + z4)**4)/243.
      q(17,12)=(y1 + y2 + y4)**2/9.
      q(17,13)=((y1 + y2 + y4)**2*(z1 + z2 + z4))/27.
      q(17,14)=((y1 + y2 + y4)**2*(z1 + z2 + z4)**2)/81.
      q(17,15)=((y1 + y2 + y4)**2*(z1 + z2 + z4)**3)/243.
      q(17,16)=(y1 + y2 + y4)**3/27.
      q(17,17)=((y1 + y2 + y4)**3*(z1 + z2 + z4))/81.
      q(17,18)=((y1 + y2 + y4)**3*(z1 + z2 + z4)**2)/243.
      q(17,19)=(y1 + y2 + y4)**4/81.
      q(17,20)=((y1 + y2 + y4)**4*(z1 + z2 + z4))/243.
      q(17,21)=(y1 + y2 + y4)**5/243.
      q(17,22)=(x1 + x2 + x4)/3.
      q(17,23)=((x1 + x2 + x4)*(z1 + z2 + z4))/9.
      q(17,24)=((x1 + x2 + x4)*(z1 + z2 + z4)**2)/27.
      q(17,25)=((x1 + x2 + x4)*(z1 + z2 + z4)**3)/81.
      q(17,26)=((x1 + x2 + x4)*(z1 + z2 + z4)**4)/243.
      q(17,27)=((x1 + x2 + x4)*(y1 + y2 + y4))/9.
      q(17,28)=((x1 + x2 + x4)*(y1 + y2 + y4)*(z1 + z2 + z4))/27.
      q(17,29)=((x1 + x2 + x4)*(y1 + y2 + y4)*(z1 + z2 + z4)**2)/81.
      q(17,30)=((x1 + x2 + x4)*(y1 + y2 + y4)*(z1 + z2 + z4)**3)/243.
      q(17,31)=((x1 + x2 + x4)*(y1 + y2 + y4)**2)/27.
      q(17,32)=((x1 + x2 + x4)*(y1 + y2 + y4)**2*(z1 + z2 + z4))/81.
      q(17,33)=((x1 + x2 + x4)*(y1 + y2 + y4)**2*(z1 + z2 + z4)**2)/243.
      q(17,34)=((x1 + x2 + x4)*(y1 + y2 + y4)**3)/81.
      q(17,35)=((x1 + x2 + x4)*(y1 + y2 + y4)**3*(z1 + z2 + z4))/243.
      q(17,36)=((x1 + x2 + x4)*(y1 + y2 + y4)**4)/243.
      q(17,37)=(x1 + x2 + x4)**2/9.
      q(17,38)=((x1 + x2 + x4)**2*(z1 + z2 + z4))/27.
      q(17,39)=((x1 + x2 + x4)**2*(z1 + z2 + z4)**2)/81.
      q(17,40)=((x1 + x2 + x4)**2*(z1 + z2 + z4)**3)/243.
      q(17,41)=((x1 + x2 + x4)**2*(y1 + y2 + y4))/27.
      q(17,42)=((x1 + x2 + x4)**2*(y1 + y2 + y4)*(z1 + z2 + z4))/81.
      q(17,43)=((x1 + x2 + x4)**2*(y1 + y2 + y4)*(z1 + z2 + z4)**2)/243.
      q(17,44)=((x1 + x2 + x4)**2*(y1 + y2 + y4)**2)/81.
      q(17,45)=((x1 + x2 + x4)**2*(y1 + y2 + y4)**2*(z1 + z2 + z4))/243.
      q(17,46)=((x1 + x2 + x4)**2*(y1 + y2 + y4)**3)/243.
      q(17,47)=(x1 + x2 + x4)**3/27.
      q(17,48)=((x1 + x2 + x4)**3*(z1 + z2 + z4))/81.
      q(17,49)=((x1 + x2 + x4)**3*(z1 + z2 + z4)**2)/243.
      q(17,50)=((x1 + x2 + x4)**3*(y1 + y2 + y4))/81.
      q(17,51)=((x1 + x2 + x4)**3*(y1 + y2 + y4)*(z1 + z2 + z4))/243.
      q(17,52)=((x1 + x2 + x4)**3*(y1 + y2 + y4)**2)/243.
      q(17,53)=(x1 + x2 + x4)**4/81.
      q(17,54)=((x1 + x2 + x4)**4*(z1 + z2 + z4))/243.
      q(17,55)=((x1 + x2 + x4)**4*(y1 + y2 + y4))/243.
      q(17,56)=(x1 + x2 + x4)**5/243.
      q(18,1)=1
      q(18,2)=(z1 + z2 + z3)/3.
      q(18,3)=(z1 + z2 + z3)**2/9.
      q(18,4)=(z1 + z2 + z3)**3/27.
      q(18,5)=(z1 + z2 + z3)**4/81.
      q(18,6)=(z1 + z2 + z3)**5/243.
      q(18,7)=(y1 + y2 + y3)/3.
      q(18,8)=((y1 + y2 + y3)*(z1 + z2 + z3))/9.
      q(18,9)=((y1 + y2 + y3)*(z1 + z2 + z3)**2)/27.
      q(18,10)=((y1 + y2 + y3)*(z1 + z2 + z3)**3)/81.
      q(18,11)=((y1 + y2 + y3)*(z1 + z2 + z3)**4)/243.
      q(18,12)=(y1 + y2 + y3)**2/9.
      q(18,13)=((y1 + y2 + y3)**2*(z1 + z2 + z3))/27.
      q(18,14)=((y1 + y2 + y3)**2*(z1 + z2 + z3)**2)/81.
      q(18,15)=((y1 + y2 + y3)**2*(z1 + z2 + z3)**3)/243.
      q(18,16)=(y1 + y2 + y3)**3/27.
      q(18,17)=((y1 + y2 + y3)**3*(z1 + z2 + z3))/81.
      q(18,18)=((y1 + y2 + y3)**3*(z1 + z2 + z3)**2)/243.
      q(18,19)=(y1 + y2 + y3)**4/81.
      q(18,20)=((y1 + y2 + y3)**4*(z1 + z2 + z3))/243.
      q(18,21)=(y1 + y2 + y3)**5/243.
      q(18,22)=(x1 + x2 + x3)/3.
      q(18,23)=((x1 + x2 + x3)*(z1 + z2 + z3))/9.
      q(18,24)=((x1 + x2 + x3)*(z1 + z2 + z3)**2)/27.
      q(18,25)=((x1 + x2 + x3)*(z1 + z2 + z3)**3)/81.
      q(18,26)=((x1 + x2 + x3)*(z1 + z2 + z3)**4)/243.
      q(18,27)=((x1 + x2 + x3)*(y1 + y2 + y3))/9.
      q(18,28)=((x1 + x2 + x3)*(y1 + y2 + y3)*(z1 + z2 + z3))/27.
      q(18,29)=((x1 + x2 + x3)*(y1 + y2 + y3)*(z1 + z2 + z3)**2)/81.
      q(18,30)=((x1 + x2 + x3)*(y1 + y2 + y3)*(z1 + z2 + z3)**3)/243.
      q(18,31)=((x1 + x2 + x3)*(y1 + y2 + y3)**2)/27.
      q(18,32)=((x1 + x2 + x3)*(y1 + y2 + y3)**2*(z1 + z2 + z3))/81.
      q(18,33)=((x1 + x2 + x3)*(y1 + y2 + y3)**2*(z1 + z2 + z3)**2)/243.
      q(18,34)=((x1 + x2 + x3)*(y1 + y2 + y3)**3)/81.
      q(18,35)=((x1 + x2 + x3)*(y1 + y2 + y3)**3*(z1 + z2 + z3))/243.
      q(18,36)=((x1 + x2 + x3)*(y1 + y2 + y3)**4)/243.
      q(18,37)=(x1 + x2 + x3)**2/9.
      q(18,38)=((x1 + x2 + x3)**2*(z1 + z2 + z3))/27.
      q(18,39)=((x1 + x2 + x3)**2*(z1 + z2 + z3)**2)/81.
      q(18,40)=((x1 + x2 + x3)**2*(z1 + z2 + z3)**3)/243.
      q(18,41)=((x1 + x2 + x3)**2*(y1 + y2 + y3))/27.
      q(18,42)=((x1 + x2 + x3)**2*(y1 + y2 + y3)*(z1 + z2 + z3))/81.
      q(18,43)=((x1 + x2 + x3)**2*(y1 + y2 + y3)*(z1 + z2 + z3)**2)/243.
      q(18,44)=((x1 + x2 + x3)**2*(y1 + y2 + y3)**2)/81.
      q(18,45)=((x1 + x2 + x3)**2*(y1 + y2 + y3)**2*(z1 + z2 + z3))/243.
      q(18,46)=((x1 + x2 + x3)**2*(y1 + y2 + y3)**3)/243.
      q(18,47)=(x1 + x2 + x3)**3/27.
      q(18,48)=((x1 + x2 + x3)**3*(z1 + z2 + z3))/81.
      q(18,49)=((x1 + x2 + x3)**3*(z1 + z2 + z3)**2)/243.
      q(18,50)=((x1 + x2 + x3)**3*(y1 + y2 + y3))/81.
      q(18,51)=((x1 + x2 + x3)**3*(y1 + y2 + y3)*(z1 + z2 + z3))/243.
      q(18,52)=((x1 + x2 + x3)**3*(y1 + y2 + y3)**2)/243.
      q(18,53)=(x1 + x2 + x3)**4/81.
      q(18,54)=((x1 + x2 + x3)**4*(z1 + z2 + z3))/243.
      q(18,55)=((x1 + x2 + x3)**4*(y1 + y2 + y3))/243.
      q(18,56)=(x1 + x2 + x3)**5/243.
      q(19,1)=1
      q(19,2)=(z1 + z3 + z4)/3.
      q(19,3)=(z1 + z3 + z4)**2/9.
      q(19,4)=(z1 + z3 + z4)**3/27.
      q(19,5)=(z1 + z3 + z4)**4/81.
      q(19,6)=(z1 + z3 + z4)**5/243.
      q(19,7)=(y1 + y3 + y4)/3.
      q(19,8)=((y1 + y3 + y4)*(z1 + z3 + z4))/9.
      q(19,9)=((y1 + y3 + y4)*(z1 + z3 + z4)**2)/27.
      q(19,10)=((y1 + y3 + y4)*(z1 + z3 + z4)**3)/81.
      q(19,11)=((y1 + y3 + y4)*(z1 + z3 + z4)**4)/243.
      q(19,12)=(y1 + y3 + y4)**2/9.
      q(19,13)=((y1 + y3 + y4)**2*(z1 + z3 + z4))/27.
      q(19,14)=((y1 + y3 + y4)**2*(z1 + z3 + z4)**2)/81.
      q(19,15)=((y1 + y3 + y4)**2*(z1 + z3 + z4)**3)/243.
      q(19,16)=(y1 + y3 + y4)**3/27.
      q(19,17)=((y1 + y3 + y4)**3*(z1 + z3 + z4))/81.
      q(19,18)=((y1 + y3 + y4)**3*(z1 + z3 + z4)**2)/243.
      q(19,19)=(y1 + y3 + y4)**4/81.
      q(19,20)=((y1 + y3 + y4)**4*(z1 + z3 + z4))/243.
      q(19,21)=(y1 + y3 + y4)**5/243.
      q(19,22)=(x1 + x3 + x4)/3.
      q(19,23)=((x1 + x3 + x4)*(z1 + z3 + z4))/9.
      q(19,24)=((x1 + x3 + x4)*(z1 + z3 + z4)**2)/27.
      q(19,25)=((x1 + x3 + x4)*(z1 + z3 + z4)**3)/81.
      q(19,26)=((x1 + x3 + x4)*(z1 + z3 + z4)**4)/243.
      q(19,27)=((x1 + x3 + x4)*(y1 + y3 + y4))/9.
      q(19,28)=((x1 + x3 + x4)*(y1 + y3 + y4)*(z1 + z3 + z4))/27.
      q(19,29)=((x1 + x3 + x4)*(y1 + y3 + y4)*(z1 + z3 + z4)**2)/81.
      q(19,30)=((x1 + x3 + x4)*(y1 + y3 + y4)*(z1 + z3 + z4)**3)/243.
      q(19,31)=((x1 + x3 + x4)*(y1 + y3 + y4)**2)/27.
      q(19,32)=((x1 + x3 + x4)*(y1 + y3 + y4)**2*(z1 + z3 + z4))/81.
      q(19,33)=((x1 + x3 + x4)*(y1 + y3 + y4)**2*(z1 + z3 + z4)**2)/243.
      q(19,34)=((x1 + x3 + x4)*(y1 + y3 + y4)**3)/81.
      q(19,35)=((x1 + x3 + x4)*(y1 + y3 + y4)**3*(z1 + z3 + z4))/243.
      q(19,36)=((x1 + x3 + x4)*(y1 + y3 + y4)**4)/243.
      q(19,37)=(x1 + x3 + x4)**2/9.
      q(19,38)=((x1 + x3 + x4)**2*(z1 + z3 + z4))/27.
      q(19,39)=((x1 + x3 + x4)**2*(z1 + z3 + z4)**2)/81.
      q(19,40)=((x1 + x3 + x4)**2*(z1 + z3 + z4)**3)/243.
      q(19,41)=((x1 + x3 + x4)**2*(y1 + y3 + y4))/27.
      q(19,42)=((x1 + x3 + x4)**2*(y1 + y3 + y4)*(z1 + z3 + z4))/81.
      q(19,43)=((x1 + x3 + x4)**2*(y1 + y3 + y4)*(z1 + z3 + z4)**2)/243.
      q(19,44)=((x1 + x3 + x4)**2*(y1 + y3 + y4)**2)/81.
      q(19,45)=((x1 + x3 + x4)**2*(y1 + y3 + y4)**2*(z1 + z3 + z4))/243.
      q(19,46)=((x1 + x3 + x4)**2*(y1 + y3 + y4)**3)/243.
      q(19,47)=(x1 + x3 + x4)**3/27.
      q(19,48)=((x1 + x3 + x4)**3*(z1 + z3 + z4))/81.
      q(19,49)=((x1 + x3 + x4)**3*(z1 + z3 + z4)**2)/243.
      q(19,50)=((x1 + x3 + x4)**3*(y1 + y3 + y4))/81.
      q(19,51)=((x1 + x3 + x4)**3*(y1 + y3 + y4)*(z1 + z3 + z4))/243.
      q(19,52)=((x1 + x3 + x4)**3*(y1 + y3 + y4)**2)/243.
      q(19,53)=(x1 + x3 + x4)**4/81.
      q(19,54)=((x1 + x3 + x4)**4*(z1 + z3 + z4))/243.
      q(19,55)=((x1 + x3 + x4)**4*(y1 + y3 + y4))/243.
      q(19,56)=(x1 + x3 + x4)**5/243.
      q(20,1)=1
      q(20,2)=(z2 + z3 + z4)/3.
      q(20,3)=(z2 + z3 + z4)**2/9.
      q(20,4)=(z2 + z3 + z4)**3/27.
      q(20,5)=(z2 + z3 + z4)**4/81.
      q(20,6)=(z2 + z3 + z4)**5/243.
      q(20,7)=(y2 + y3 + y4)/3.
      q(20,8)=((y2 + y3 + y4)*(z2 + z3 + z4))/9.
      q(20,9)=((y2 + y3 + y4)*(z2 + z3 + z4)**2)/27.
      q(20,10)=((y2 + y3 + y4)*(z2 + z3 + z4)**3)/81.
      q(20,11)=((y2 + y3 + y4)*(z2 + z3 + z4)**4)/243.
      q(20,12)=(y2 + y3 + y4)**2/9.
      q(20,13)=((y2 + y3 + y4)**2*(z2 + z3 + z4))/27.
      q(20,14)=((y2 + y3 + y4)**2*(z2 + z3 + z4)**2)/81.
      q(20,15)=((y2 + y3 + y4)**2*(z2 + z3 + z4)**3)/243.
      q(20,16)=(y2 + y3 + y4)**3/27.
      q(20,17)=((y2 + y3 + y4)**3*(z2 + z3 + z4))/81.
      q(20,18)=((y2 + y3 + y4)**3*(z2 + z3 + z4)**2)/243.
      q(20,19)=(y2 + y3 + y4)**4/81.
      q(20,20)=((y2 + y3 + y4)**4*(z2 + z3 + z4))/243.
      q(20,21)=(y2 + y3 + y4)**5/243.
      q(20,22)=(x2 + x3 + x4)/3.
      q(20,23)=((x2 + x3 + x4)*(z2 + z3 + z4))/9.
      q(20,24)=((x2 + x3 + x4)*(z2 + z3 + z4)**2)/27.
      q(20,25)=((x2 + x3 + x4)*(z2 + z3 + z4)**3)/81.
      q(20,26)=((x2 + x3 + x4)*(z2 + z3 + z4)**4)/243.
      q(20,27)=((x2 + x3 + x4)*(y2 + y3 + y4))/9.
      q(20,28)=((x2 + x3 + x4)*(y2 + y3 + y4)*(z2 + z3 + z4))/27.
      q(20,29)=((x2 + x3 + x4)*(y2 + y3 + y4)*(z2 + z3 + z4)**2)/81.
      q(20,30)=((x2 + x3 + x4)*(y2 + y3 + y4)*(z2 + z3 + z4)**3)/243.
      q(20,31)=((x2 + x3 + x4)*(y2 + y3 + y4)**2)/27.
      q(20,32)=((x2 + x3 + x4)*(y2 + y3 + y4)**2*(z2 + z3 + z4))/81.
      q(20,33)=((x2 + x3 + x4)*(y2 + y3 + y4)**2*(z2 + z3 + z4)**2)/243.
      q(20,34)=((x2 + x3 + x4)*(y2 + y3 + y4)**3)/81.
      q(20,35)=((x2 + x3 + x4)*(y2 + y3 + y4)**3*(z2 + z3 + z4))/243.
      q(20,36)=((x2 + x3 + x4)*(y2 + y3 + y4)**4)/243.
      q(20,37)=(x2 + x3 + x4)**2/9.
      q(20,38)=((x2 + x3 + x4)**2*(z2 + z3 + z4))/27.
      q(20,39)=((x2 + x3 + x4)**2*(z2 + z3 + z4)**2)/81.
      q(20,40)=((x2 + x3 + x4)**2*(z2 + z3 + z4)**3)/243.
      q(20,41)=((x2 + x3 + x4)**2*(y2 + y3 + y4))/27.
      q(20,42)=((x2 + x3 + x4)**2*(y2 + y3 + y4)*(z2 + z3 + z4))/81.
      q(20,43)=((x2 + x3 + x4)**2*(y2 + y3 + y4)*(z2 + z3 + z4)**2)/243.
      q(20,44)=((x2 + x3 + x4)**2*(y2 + y3 + y4)**2)/81.
      q(20,45)=((x2 + x3 + x4)**2*(y2 + y3 + y4)**2*(z2 + z3 + z4))/243.
      q(20,46)=((x2 + x3 + x4)**2*(y2 + y3 + y4)**3)/243.
      q(20,47)=(x2 + x3 + x4)**3/27.
      q(20,48)=((x2 + x3 + x4)**3*(z2 + z3 + z4))/81.
      q(20,49)=((x2 + x3 + x4)**3*(z2 + z3 + z4)**2)/243.
      q(20,50)=((x2 + x3 + x4)**3*(y2 + y3 + y4))/81.
      q(20,51)=((x2 + x3 + x4)**3*(y2 + y3 + y4)*(z2 + z3 + z4))/243.
      q(20,52)=((x2 + x3 + x4)**3*(y2 + y3 + y4)**2)/243.
      q(20,53)=(x2 + x3 + x4)**4/81.
      q(20,54)=((x2 + x3 + x4)**4*(z2 + z3 + z4))/243.
      q(20,55)=((x2 + x3 + x4)**4*(y2 + y3 + y4))/243.
      q(20,56)=(x2 + x3 + x4)**5/243.
      q(21,1)=0
      q(21,2)=0
      q(21,3)=0
      q(21,4)=0
      q(21,5)=0
      q(21,6)=0
      q(21,7)=0
      q(21,8)=0
      q(21,9)=0
      q(21,10)=0
      q(21,11)=0
      q(21,12)=0
      q(21,13)=0
      q(21,14)=0
      q(21,15)=0
      q(21,16)=0
      q(21,17)=0
      q(21,18)=0
      q(21,19)=0
      q(21,20)=0
      q(21,21)=0
      q(21,22)=1
      q(21,23)=(z1 + z2 + z4)/3.
      q(21,24)=(z1 + z2 + z4)**2/9.
      q(21,25)=(z1 + z2 + z4)**3/27.
      q(21,26)=(z1 + z2 + z4)**4/81.
      q(21,27)=(y1 + y2 + y4)/3.
      q(21,28)=((y1 + y2 + y4)*(z1 + z2 + z4))/9.
      q(21,29)=((y1 + y2 + y4)*(z1 + z2 + z4)**2)/27.
      q(21,30)=((y1 + y2 + y4)*(z1 + z2 + z4)**3)/81.
      q(21,31)=(y1 + y2 + y4)**2/9.
      q(21,32)=((y1 + y2 + y4)**2*(z1 + z2 + z4))/27.
      q(21,33)=((y1 + y2 + y4)**2*(z1 + z2 + z4)**2)/81.
      q(21,34)=(y1 + y2 + y4)**3/27.
      q(21,35)=((y1 + y2 + y4)**3*(z1 + z2 + z4))/81.
      q(21,36)=(y1 + y2 + y4)**4/81.
      q(21,37)=(2*(x1 + x2 + x4))/3.
      q(21,38)=(2*(x1 + x2 + x4)*(z1 + z2 + z4))/9.
      q(21,39)=(2*(x1 + x2 + x4)*(z1 + z2 + z4)**2)/27.
      q(21,40)=(2*(x1 + x2 + x4)*(z1 + z2 + z4)**3)/81.
      q(21,41)=(2*(x1 + x2 + x4)*(y1 + y2 + y4))/9.
      q(21,42)=(2*(x1 + x2 + x4)*(y1 + y2 + y4)*(z1 + z2 + z4))/27.
      q(21,43)=(2*(x1 + x2 + x4)*(y1 + y2 + y4)*(z1 + z2 + z4)**2)/81.
      q(21,44)=(2*(x1 + x2 + x4)*(y1 + y2 + y4)**2)/27.
      q(21,45)=(2*(x1 + x2 + x4)*(y1 + y2 + y4)**2*(z1 + z2 + z4))/81.
      q(21,46)=(2*(x1 + x2 + x4)*(y1 + y2 + y4)**3)/81.
      q(21,47)=(x1 + x2 + x4)**2/3.
      q(21,48)=((x1 + x2 + x4)**2*(z1 + z2 + z4))/9.
      q(21,49)=((x1 + x2 + x4)**2*(z1 + z2 + z4)**2)/27.
      q(21,50)=((x1 + x2 + x4)**2*(y1 + y2 + y4))/9.
      q(21,51)=((x1 + x2 + x4)**2*(y1 + y2 + y4)*(z1 + z2 + z4))/27.
      q(21,52)=((x1 + x2 + x4)**2*(y1 + y2 + y4)**2)/27.
      q(21,53)=(4*(x1 + x2 + x4)**3)/27.
      q(21,54)=(4*(x1 + x2 + x4)**3*(z1 + z2 + z4))/81.
      q(21,55)=(4*(x1 + x2 + x4)**3*(y1 + y2 + y4))/81.
      q(21,56)=(5*(x1 + x2 + x4)**4)/81.
      q(22,1)=0
      q(22,2)=0
      q(22,3)=0
      q(22,4)=0
      q(22,5)=0
      q(22,6)=0
      q(22,7)=1
      q(22,8)=(z1 + z2 + z4)/3.
      q(22,9)=(z1 + z2 + z4)**2/9.
      q(22,10)=(z1 + z2 + z4)**3/27.
      q(22,11)=(z1 + z2 + z4)**4/81.
      q(22,12)=(2*(y1 + y2 + y4))/3.
      q(22,13)=(2*(y1 + y2 + y4)*(z1 + z2 + z4))/9.
      q(22,14)=(2*(y1 + y2 + y4)*(z1 + z2 + z4)**2)/27.
      q(22,15)=(2*(y1 + y2 + y4)*(z1 + z2 + z4)**3)/81.
      q(22,16)=(y1 + y2 + y4)**2/3.
      q(22,17)=((y1 + y2 + y4)**2*(z1 + z2 + z4))/9.
      q(22,18)=((y1 + y2 + y4)**2*(z1 + z2 + z4)**2)/27.
      q(22,19)=(4*(y1 + y2 + y4)**3)/27.
      q(22,20)=(4*(y1 + y2 + y4)**3*(z1 + z2 + z4))/81.
      q(22,21)=(5*(y1 + y2 + y4)**4)/81.
      q(22,22)=0
      q(22,23)=0
      q(22,24)=0
      q(22,25)=0
      q(22,26)=0
      q(22,27)=(x1 + x2 + x4)/3.
      q(22,28)=((x1 + x2 + x4)*(z1 + z2 + z4))/9.
      q(22,29)=((x1 + x2 + x4)*(z1 + z2 + z4)**2)/27.
      q(22,30)=((x1 + x2 + x4)*(z1 + z2 + z4)**3)/81.
      q(22,31)=(2*(x1 + x2 + x4)*(y1 + y2 + y4))/9.
      q(22,32)=(2*(x1 + x2 + x4)*(y1 + y2 + y4)*(z1 + z2 + z4))/27.
      q(22,33)=(2*(x1 + x2 + x4)*(y1 + y2 + y4)*(z1 + z2 + z4)**2)/81.
      q(22,34)=((x1 + x2 + x4)*(y1 + y2 + y4)**2)/9.
      q(22,35)=((x1 + x2 + x4)*(y1 + y2 + y4)**2*(z1 + z2 + z4))/27.
      q(22,36)=(4*(x1 + x2 + x4)*(y1 + y2 + y4)**3)/81.
      q(22,37)=0
      q(22,38)=0
      q(22,39)=0
      q(22,40)=0
      q(22,41)=(x1 + x2 + x4)**2/9.
      q(22,42)=((x1 + x2 + x4)**2*(z1 + z2 + z4))/27.
      q(22,43)=((x1 + x2 + x4)**2*(z1 + z2 + z4)**2)/81.
      q(22,44)=(2*(x1 + x2 + x4)**2*(y1 + y2 + y4))/27.
      q(22,45)=(2*(x1 + x2 + x4)**2*(y1 + y2 + y4)*(z1 + z2 + z4))/81.
      q(22,46)=((x1 + x2 + x4)**2*(y1 + y2 + y4)**2)/27.
      q(22,47)=0
      q(22,48)=0
      q(22,49)=0
      q(22,50)=(x1 + x2 + x4)**3/27.
      q(22,51)=((x1 + x2 + x4)**3*(z1 + z2 + z4))/81.
      q(22,52)=(2*(x1 + x2 + x4)**3*(y1 + y2 + y4))/81.
      q(22,53)=0
      q(22,54)=0
      q(22,55)=(x1 + x2 + x4)**4/81.
      q(22,56)=0
      q(23,1)=0
      q(23,2)=1
      q(23,3)=(2*(z1 + z2 + z4))/3.
      q(23,4)=(z1 + z2 + z4)**2/3.
      q(23,5)=(4*(z1 + z2 + z4)**3)/27.
      q(23,6)=(5*(z1 + z2 + z4)**4)/81.
      q(23,7)=0
      q(23,8)=(y1 + y2 + y4)/3.
      q(23,9)=(2*(y1 + y2 + y4)*(z1 + z2 + z4))/9.
      q(23,10)=((y1 + y2 + y4)*(z1 + z2 + z4)**2)/9.
      q(23,11)=(4*(y1 + y2 + y4)*(z1 + z2 + z4)**3)/81.
      q(23,12)=0
      q(23,13)=(y1 + y2 + y4)**2/9.
      q(23,14)=(2*(y1 + y2 + y4)**2*(z1 + z2 + z4))/27.
      q(23,15)=((y1 + y2 + y4)**2*(z1 + z2 + z4)**2)/27.
      q(23,16)=0
      q(23,17)=(y1 + y2 + y4)**3/27.
      q(23,18)=(2*(y1 + y2 + y4)**3*(z1 + z2 + z4))/81.
      q(23,19)=0
      q(23,20)=(y1 + y2 + y4)**4/81.
      q(23,21)=0
      q(23,22)=0
      q(23,23)=(x1 + x2 + x4)/3.
      q(23,24)=(2*(x1 + x2 + x4)*(z1 + z2 + z4))/9.
      q(23,25)=((x1 + x2 + x4)*(z1 + z2 + z4)**2)/9.
      q(23,26)=(4*(x1 + x2 + x4)*(z1 + z2 + z4)**3)/81.
      q(23,27)=0
      q(23,28)=((x1 + x2 + x4)*(y1 + y2 + y4))/9.
      q(23,29)=(2*(x1 + x2 + x4)*(y1 + y2 + y4)*(z1 + z2 + z4))/27.
      q(23,30)=((x1 + x2 + x4)*(y1 + y2 + y4)*(z1 + z2 + z4)**2)/27.
      q(23,31)=0
      q(23,32)=((x1 + x2 + x4)*(y1 + y2 + y4)**2)/27.
      q(23,33)=(2*(x1 + x2 + x4)*(y1 + y2 + y4)**2*(z1 + z2 + z4))/81.
      q(23,34)=0
      q(23,35)=((x1 + x2 + x4)*(y1 + y2 + y4)**3)/81.
      q(23,36)=0
      q(23,37)=0
      q(23,38)=(x1 + x2 + x4)**2/9.
      q(23,39)=(2*(x1 + x2 + x4)**2*(z1 + z2 + z4))/27.
      q(23,40)=((x1 + x2 + x4)**2*(z1 + z2 + z4)**2)/27.
      q(23,41)=0
      q(23,42)=((x1 + x2 + x4)**2*(y1 + y2 + y4))/27.
      q(23,43)=(2*(x1 + x2 + x4)**2*(y1 + y2 + y4)*(z1 + z2 + z4))/81.
      q(23,44)=0
      q(23,45)=((x1 + x2 + x4)**2*(y1 + y2 + y4)**2)/81.
      q(23,46)=0
      q(23,47)=0
      q(23,48)=(x1 + x2 + x4)**3/27.
      q(23,49)=(2*(x1 + x2 + x4)**3*(z1 + z2 + z4))/81.
      q(23,50)=0
      q(23,51)=((x1 + x2 + x4)**3*(y1 + y2 + y4))/81.
      q(23,52)=0
      q(23,53)=0
      q(23,54)=(x1 + x2 + x4)**4/81.
      q(23,55)=0
      q(23,56)=0
      q(24,1)=0
      q(24,2)=0
      q(24,3)=0
      q(24,4)=0
      q(24,5)=0
      q(24,6)=0
      q(24,7)=0
      q(24,8)=0
      q(24,9)=0
      q(24,10)=0
      q(24,11)=0
      q(24,12)=0
      q(24,13)=0
      q(24,14)=0
      q(24,15)=0
      q(24,16)=0
      q(24,17)=0
      q(24,18)=0
      q(24,19)=0
      q(24,20)=0
      q(24,21)=0
      q(24,22)=1
      q(24,23)=(z1 + z2 + z3)/3.
      q(24,24)=(z1 + z2 + z3)**2/9.
      q(24,25)=(z1 + z2 + z3)**3/27.
      q(24,26)=(z1 + z2 + z3)**4/81.
      q(24,27)=(y1 + y2 + y3)/3.
      q(24,28)=((y1 + y2 + y3)*(z1 + z2 + z3))/9.
      q(24,29)=((y1 + y2 + y3)*(z1 + z2 + z3)**2)/27.
      q(24,30)=((y1 + y2 + y3)*(z1 + z2 + z3)**3)/81.
      q(24,31)=(y1 + y2 + y3)**2/9.
      q(24,32)=((y1 + y2 + y3)**2*(z1 + z2 + z3))/27.
      q(24,33)=((y1 + y2 + y3)**2*(z1 + z2 + z3)**2)/81.
      q(24,34)=(y1 + y2 + y3)**3/27.
      q(24,35)=((y1 + y2 + y3)**3*(z1 + z2 + z3))/81.
      q(24,36)=(y1 + y2 + y3)**4/81.
      q(24,37)=(2*(x1 + x2 + x3))/3.
      q(24,38)=(2*(x1 + x2 + x3)*(z1 + z2 + z3))/9.
      q(24,39)=(2*(x1 + x2 + x3)*(z1 + z2 + z3)**2)/27.
      q(24,40)=(2*(x1 + x2 + x3)*(z1 + z2 + z3)**3)/81.
      q(24,41)=(2*(x1 + x2 + x3)*(y1 + y2 + y3))/9.
      q(24,42)=(2*(x1 + x2 + x3)*(y1 + y2 + y3)*(z1 + z2 + z3))/27.
      q(24,43)=(2*(x1 + x2 + x3)*(y1 + y2 + y3)*(z1 + z2 + z3)**2)/81.
      q(24,44)=(2*(x1 + x2 + x3)*(y1 + y2 + y3)**2)/27.
      q(24,45)=(2*(x1 + x2 + x3)*(y1 + y2 + y3)**2*(z1 + z2 + z3))/81.
      q(24,46)=(2*(x1 + x2 + x3)*(y1 + y2 + y3)**3)/81.
      q(24,47)=(x1 + x2 + x3)**2/3.
      q(24,48)=((x1 + x2 + x3)**2*(z1 + z2 + z3))/9.
      q(24,49)=((x1 + x2 + x3)**2*(z1 + z2 + z3)**2)/27.
      q(24,50)=((x1 + x2 + x3)**2*(y1 + y2 + y3))/9.
      q(24,51)=((x1 + x2 + x3)**2*(y1 + y2 + y3)*(z1 + z2 + z3))/27.
      q(24,52)=((x1 + x2 + x3)**2*(y1 + y2 + y3)**2)/27.
      q(24,53)=(4*(x1 + x2 + x3)**3)/27.
      q(24,54)=(4*(x1 + x2 + x3)**3*(z1 + z2 + z3))/81.
      q(24,55)=(4*(x1 + x2 + x3)**3*(y1 + y2 + y3))/81.
      q(24,56)=(5*(x1 + x2 + x3)**4)/81.
      q(25,1)=0
      q(25,2)=0
      q(25,3)=0
      q(25,4)=0
      q(25,5)=0
      q(25,6)=0
      q(25,7)=1
      q(25,8)=(z1 + z2 + z3)/3.
      q(25,9)=(z1 + z2 + z3)**2/9.
      q(25,10)=(z1 + z2 + z3)**3/27.
      q(25,11)=(z1 + z2 + z3)**4/81.
      q(25,12)=(2*(y1 + y2 + y3))/3.
      q(25,13)=(2*(y1 + y2 + y3)*(z1 + z2 + z3))/9.
      q(25,14)=(2*(y1 + y2 + y3)*(z1 + z2 + z3)**2)/27.
      q(25,15)=(2*(y1 + y2 + y3)*(z1 + z2 + z3)**3)/81.
      q(25,16)=(y1 + y2 + y3)**2/3.
      q(25,17)=((y1 + y2 + y3)**2*(z1 + z2 + z3))/9.
      q(25,18)=((y1 + y2 + y3)**2*(z1 + z2 + z3)**2)/27.
      q(25,19)=(4*(y1 + y2 + y3)**3)/27.
      q(25,20)=(4*(y1 + y2 + y3)**3*(z1 + z2 + z3))/81.
      q(25,21)=(5*(y1 + y2 + y3)**4)/81.
      q(25,22)=0
      q(25,23)=0
      q(25,24)=0
      q(25,25)=0
      q(25,26)=0
      q(25,27)=(x1 + x2 + x3)/3.
      q(25,28)=((x1 + x2 + x3)*(z1 + z2 + z3))/9.
      q(25,29)=((x1 + x2 + x3)*(z1 + z2 + z3)**2)/27.
      q(25,30)=((x1 + x2 + x3)*(z1 + z2 + z3)**3)/81.
      q(25,31)=(2*(x1 + x2 + x3)*(y1 + y2 + y3))/9.
      q(25,32)=(2*(x1 + x2 + x3)*(y1 + y2 + y3)*(z1 + z2 + z3))/27.
      q(25,33)=(2*(x1 + x2 + x3)*(y1 + y2 + y3)*(z1 + z2 + z3)**2)/81.
      q(25,34)=((x1 + x2 + x3)*(y1 + y2 + y3)**2)/9.
      q(25,35)=((x1 + x2 + x3)*(y1 + y2 + y3)**2*(z1 + z2 + z3))/27.
      q(25,36)=(4*(x1 + x2 + x3)*(y1 + y2 + y3)**3)/81.
      q(25,37)=0
      q(25,38)=0
      q(25,39)=0
      q(25,40)=0
      q(25,41)=(x1 + x2 + x3)**2/9.
      q(25,42)=((x1 + x2 + x3)**2*(z1 + z2 + z3))/27.
      q(25,43)=((x1 + x2 + x3)**2*(z1 + z2 + z3)**2)/81.
      q(25,44)=(2*(x1 + x2 + x3)**2*(y1 + y2 + y3))/27.
      q(25,45)=(2*(x1 + x2 + x3)**2*(y1 + y2 + y3)*(z1 + z2 + z3))/81.
      q(25,46)=((x1 + x2 + x3)**2*(y1 + y2 + y3)**2)/27.
      q(25,47)=0
      q(25,48)=0
      q(25,49)=0
      q(25,50)=(x1 + x2 + x3)**3/27.
      q(25,51)=((x1 + x2 + x3)**3*(z1 + z2 + z3))/81.
      q(25,52)=(2*(x1 + x2 + x3)**3*(y1 + y2 + y3))/81.
      q(25,53)=0
      q(25,54)=0
      q(25,55)=(x1 + x2 + x3)**4/81.
      q(25,56)=0
      q(26,1)=0
      q(26,2)=1
      q(26,3)=(2*(z1 + z2 + z3))/3.
      q(26,4)=(z1 + z2 + z3)**2/3.
      q(26,5)=(4*(z1 + z2 + z3)**3)/27.
      q(26,6)=(5*(z1 + z2 + z3)**4)/81.
      q(26,7)=0
      q(26,8)=(y1 + y2 + y3)/3.
      q(26,9)=(2*(y1 + y2 + y3)*(z1 + z2 + z3))/9.
      q(26,10)=((y1 + y2 + y3)*(z1 + z2 + z3)**2)/9.
      q(26,11)=(4*(y1 + y2 + y3)*(z1 + z2 + z3)**3)/81.
      q(26,12)=0
      q(26,13)=(y1 + y2 + y3)**2/9.
      q(26,14)=(2*(y1 + y2 + y3)**2*(z1 + z2 + z3))/27.
      q(26,15)=((y1 + y2 + y3)**2*(z1 + z2 + z3)**2)/27.
      q(26,16)=0
      q(26,17)=(y1 + y2 + y3)**3/27.
      q(26,18)=(2*(y1 + y2 + y3)**3*(z1 + z2 + z3))/81.
      q(26,19)=0
      q(26,20)=(y1 + y2 + y3)**4/81.
      q(26,21)=0
      q(26,22)=0
      q(26,23)=(x1 + x2 + x3)/3.
      q(26,24)=(2*(x1 + x2 + x3)*(z1 + z2 + z3))/9.
      q(26,25)=((x1 + x2 + x3)*(z1 + z2 + z3)**2)/9.
      q(26,26)=(4*(x1 + x2 + x3)*(z1 + z2 + z3)**3)/81.
      q(26,27)=0
      q(26,28)=((x1 + x2 + x3)*(y1 + y2 + y3))/9.
      q(26,29)=(2*(x1 + x2 + x3)*(y1 + y2 + y3)*(z1 + z2 + z3))/27.
      q(26,30)=((x1 + x2 + x3)*(y1 + y2 + y3)*(z1 + z2 + z3)**2)/27.
      q(26,31)=0
      q(26,32)=((x1 + x2 + x3)*(y1 + y2 + y3)**2)/27.
      q(26,33)=(2*(x1 + x2 + x3)*(y1 + y2 + y3)**2*(z1 + z2 + z3))/81.
      q(26,34)=0
      q(26,35)=((x1 + x2 + x3)*(y1 + y2 + y3)**3)/81.
      q(26,36)=0
      q(26,37)=0
      q(26,38)=(x1 + x2 + x3)**2/9.
      q(26,39)=(2*(x1 + x2 + x3)**2*(z1 + z2 + z3))/27.
      q(26,40)=((x1 + x2 + x3)**2*(z1 + z2 + z3)**2)/27.
      q(26,41)=0
      q(26,42)=((x1 + x2 + x3)**2*(y1 + y2 + y3))/27.
      q(26,43)=(2*(x1 + x2 + x3)**2*(y1 + y2 + y3)*(z1 + z2 + z3))/81.
      q(26,44)=0
      q(26,45)=((x1 + x2 + x3)**2*(y1 + y2 + y3)**2)/81.
      q(26,46)=0
      q(26,47)=0
      q(26,48)=(x1 + x2 + x3)**3/27.
      q(26,49)=(2*(x1 + x2 + x3)**3*(z1 + z2 + z3))/81.
      q(26,50)=0
      q(26,51)=((x1 + x2 + x3)**3*(y1 + y2 + y3))/81.
      q(26,52)=0
      q(26,53)=0
      q(26,54)=(x1 + x2 + x3)**4/81.
      q(26,55)=0
      q(26,56)=0
      q(27,1)=0
      q(27,2)=0
      q(27,3)=0
      q(27,4)=0
      q(27,5)=0
      q(27,6)=0
      q(27,7)=0
      q(27,8)=0
      q(27,9)=0
      q(27,10)=0
      q(27,11)=0
      q(27,12)=0
      q(27,13)=0
      q(27,14)=0
      q(27,15)=0
      q(27,16)=0
      q(27,17)=0
      q(27,18)=0
      q(27,19)=0
      q(27,20)=0
      q(27,21)=0
      q(27,22)=1
      q(27,23)=(z2 + z3 + z4)/3.
      q(27,24)=(z2 + z3 + z4)**2/9.
      q(27,25)=(z2 + z3 + z4)**3/27.
      q(27,26)=(z2 + z3 + z4)**4/81.
      q(27,27)=(y2 + y3 + y4)/3.
      q(27,28)=((y2 + y3 + y4)*(z2 + z3 + z4))/9.
      q(27,29)=((y2 + y3 + y4)*(z2 + z3 + z4)**2)/27.
      q(27,30)=((y2 + y3 + y4)*(z2 + z3 + z4)**3)/81.
      q(27,31)=(y2 + y3 + y4)**2/9.
      q(27,32)=((y2 + y3 + y4)**2*(z2 + z3 + z4))/27.
      q(27,33)=((y2 + y3 + y4)**2*(z2 + z3 + z4)**2)/81.
      q(27,34)=(y2 + y3 + y4)**3/27.
      q(27,35)=((y2 + y3 + y4)**3*(z2 + z3 + z4))/81.
      q(27,36)=(y2 + y3 + y4)**4/81.
      q(27,37)=(2*(x2 + x3 + x4))/3.
      q(27,38)=(2*(x2 + x3 + x4)*(z2 + z3 + z4))/9.
      q(27,39)=(2*(x2 + x3 + x4)*(z2 + z3 + z4)**2)/27.
      q(27,40)=(2*(x2 + x3 + x4)*(z2 + z3 + z4)**3)/81.
      q(27,41)=(2*(x2 + x3 + x4)*(y2 + y3 + y4))/9.
      q(27,42)=(2*(x2 + x3 + x4)*(y2 + y3 + y4)*(z2 + z3 + z4))/27.
      q(27,43)=(2*(x2 + x3 + x4)*(y2 + y3 + y4)*(z2 + z3 + z4)**2)/81.
      q(27,44)=(2*(x2 + x3 + x4)*(y2 + y3 + y4)**2)/27.
      q(27,45)=(2*(x2 + x3 + x4)*(y2 + y3 + y4)**2*(z2 + z3 + z4))/81.
      q(27,46)=(2*(x2 + x3 + x4)*(y2 + y3 + y4)**3)/81.
      q(27,47)=(x2 + x3 + x4)**2/3.
      q(27,48)=((x2 + x3 + x4)**2*(z2 + z3 + z4))/9.
      q(27,49)=((x2 + x3 + x4)**2*(z2 + z3 + z4)**2)/27.
      q(27,50)=((x2 + x3 + x4)**2*(y2 + y3 + y4))/9.
      q(27,51)=((x2 + x3 + x4)**2*(y2 + y3 + y4)*(z2 + z3 + z4))/27.
      q(27,52)=((x2 + x3 + x4)**2*(y2 + y3 + y4)**2)/27.
      q(27,53)=(4*(x2 + x3 + x4)**3)/27.
      q(27,54)=(4*(x2 + x3 + x4)**3*(z2 + z3 + z4))/81.
      q(27,55)=(4*(x2 + x3 + x4)**3*(y2 + y3 + y4))/81.
      q(27,56)=(5*(x2 + x3 + x4)**4)/81.
      q(28,1)=0
      q(28,2)=0
      q(28,3)=0
      q(28,4)=0
      q(28,5)=0
      q(28,6)=0
      q(28,7)=1
      q(28,8)=(z2 + z3 + z4)/3.
      q(28,9)=(z2 + z3 + z4)**2/9.
      q(28,10)=(z2 + z3 + z4)**3/27.
      q(28,11)=(z2 + z3 + z4)**4/81.
      q(28,12)=(2*(y2 + y3 + y4))/3.
      q(28,13)=(2*(y2 + y3 + y4)*(z2 + z3 + z4))/9.
      q(28,14)=(2*(y2 + y3 + y4)*(z2 + z3 + z4)**2)/27.
      q(28,15)=(2*(y2 + y3 + y4)*(z2 + z3 + z4)**3)/81.
      q(28,16)=(y2 + y3 + y4)**2/3.
      q(28,17)=((y2 + y3 + y4)**2*(z2 + z3 + z4))/9.
      q(28,18)=((y2 + y3 + y4)**2*(z2 + z3 + z4)**2)/27.
      q(28,19)=(4*(y2 + y3 + y4)**3)/27.
      q(28,20)=(4*(y2 + y3 + y4)**3*(z2 + z3 + z4))/81.
      q(28,21)=(5*(y2 + y3 + y4)**4)/81.
      q(28,22)=0
      q(28,23)=0
      q(28,24)=0
      q(28,25)=0
      q(28,26)=0
      q(28,27)=(x2 + x3 + x4)/3.
      q(28,28)=((x2 + x3 + x4)*(z2 + z3 + z4))/9.
      q(28,29)=((x2 + x3 + x4)*(z2 + z3 + z4)**2)/27.
      q(28,30)=((x2 + x3 + x4)*(z2 + z3 + z4)**3)/81.
      q(28,31)=(2*(x2 + x3 + x4)*(y2 + y3 + y4))/9.
      q(28,32)=(2*(x2 + x3 + x4)*(y2 + y3 + y4)*(z2 + z3 + z4))/27.
      q(28,33)=(2*(x2 + x3 + x4)*(y2 + y3 + y4)*(z2 + z3 + z4)**2)/81.
      q(28,34)=((x2 + x3 + x4)*(y2 + y3 + y4)**2)/9.
      q(28,35)=((x2 + x3 + x4)*(y2 + y3 + y4)**2*(z2 + z3 + z4))/27.
      q(28,36)=(4*(x2 + x3 + x4)*(y2 + y3 + y4)**3)/81.
      q(28,37)=0
      q(28,38)=0
      q(28,39)=0
      q(28,40)=0
      q(28,41)=(x2 + x3 + x4)**2/9.
      q(28,42)=((x2 + x3 + x4)**2*(z2 + z3 + z4))/27.
      q(28,43)=((x2 + x3 + x4)**2*(z2 + z3 + z4)**2)/81.
      q(28,44)=(2*(x2 + x3 + x4)**2*(y2 + y3 + y4))/27.
      q(28,45)=(2*(x2 + x3 + x4)**2*(y2 + y3 + y4)*(z2 + z3 + z4))/81.
      q(28,46)=((x2 + x3 + x4)**2*(y2 + y3 + y4)**2)/27.
      q(28,47)=0
      q(28,48)=0
      q(28,49)=0
      q(28,50)=(x2 + x3 + x4)**3/27.
      q(28,51)=((x2 + x3 + x4)**3*(z2 + z3 + z4))/81.
      q(28,52)=(2*(x2 + x3 + x4)**3*(y2 + y3 + y4))/81.
      q(28,53)=0
      q(28,54)=0
      q(28,55)=(x2 + x3 + x4)**4/81.
      q(28,56)=0
      q(29,1)=0
      q(29,2)=1
      q(29,3)=(2*(z2 + z3 + z4))/3.
      q(29,4)=(z2 + z3 + z4)**2/3.
      q(29,5)=(4*(z2 + z3 + z4)**3)/27.
      q(29,6)=(5*(z2 + z3 + z4)**4)/81.
      q(29,7)=0
      q(29,8)=(y2 + y3 + y4)/3.
      q(29,9)=(2*(y2 + y3 + y4)*(z2 + z3 + z4))/9.
      q(29,10)=((y2 + y3 + y4)*(z2 + z3 + z4)**2)/9.
      q(29,11)=(4*(y2 + y3 + y4)*(z2 + z3 + z4)**3)/81.
      q(29,12)=0
      q(29,13)=(y2 + y3 + y4)**2/9.
      q(29,14)=(2*(y2 + y3 + y4)**2*(z2 + z3 + z4))/27.
      q(29,15)=((y2 + y3 + y4)**2*(z2 + z3 + z4)**2)/27.
      q(29,16)=0
      q(29,17)=(y2 + y3 + y4)**3/27.
      q(29,18)=(2*(y2 + y3 + y4)**3*(z2 + z3 + z4))/81.
      q(29,19)=0
      q(29,20)=(y2 + y3 + y4)**4/81.
      q(29,21)=0
      q(29,22)=0
      q(29,23)=(x2 + x3 + x4)/3.
      q(29,24)=(2*(x2 + x3 + x4)*(z2 + z3 + z4))/9.
      q(29,25)=((x2 + x3 + x4)*(z2 + z3 + z4)**2)/9.
      q(29,26)=(4*(x2 + x3 + x4)*(z2 + z3 + z4)**3)/81.
      q(29,27)=0
      q(29,28)=((x2 + x3 + x4)*(y2 + y3 + y4))/9.
      q(29,29)=(2*(x2 + x3 + x4)*(y2 + y3 + y4)*(z2 + z3 + z4))/27.
      q(29,30)=((x2 + x3 + x4)*(y2 + y3 + y4)*(z2 + z3 + z4)**2)/27.
      q(29,31)=0
      q(29,32)=((x2 + x3 + x4)*(y2 + y3 + y4)**2)/27.
      q(29,33)=(2*(x2 + x3 + x4)*(y2 + y3 + y4)**2*(z2 + z3 + z4))/81.
      q(29,34)=0
      q(29,35)=((x2 + x3 + x4)*(y2 + y3 + y4)**3)/81.
      q(29,36)=0
      q(29,37)=0
      q(29,38)=(x2 + x3 + x4)**2/9.
      q(29,39)=(2*(x2 + x3 + x4)**2*(z2 + z3 + z4))/27.
      q(29,40)=((x2 + x3 + x4)**2*(z2 + z3 + z4)**2)/27.
      q(29,41)=0
      q(29,42)=((x2 + x3 + x4)**2*(y2 + y3 + y4))/27.
      q(29,43)=(2*(x2 + x3 + x4)**2*(y2 + y3 + y4)*(z2 + z3 + z4))/81.
      q(29,44)=0
      q(29,45)=((x2 + x3 + x4)**2*(y2 + y3 + y4)**2)/81.
      q(29,46)=0
      q(29,47)=0
      q(29,48)=(x2 + x3 + x4)**3/27.
      q(29,49)=(2*(x2 + x3 + x4)**3*(z2 + z3 + z4))/81.
      q(29,50)=0
      q(29,51)=((x2 + x3 + x4)**3*(y2 + y3 + y4))/81.
      q(29,52)=0
      q(29,53)=0
      q(29,54)=(x2 + x3 + x4)**4/81.
      q(29,55)=0
      q(29,56)=0
      q(30,1)=0
      q(30,2)=0
      q(30,3)=0
      q(30,4)=0
      q(30,5)=0
      q(30,6)=0
      q(30,7)=0
      q(30,8)=0
      q(30,9)=0
      q(30,10)=0
      q(30,11)=0
      q(30,12)=0
      q(30,13)=0
      q(30,14)=0
      q(30,15)=0
      q(30,16)=0
      q(30,17)=0
      q(30,18)=0
      q(30,19)=0
      q(30,20)=0
      q(30,21)=0
      q(30,22)=1
      q(30,23)=(z1 + z3 + z4)/3.
      q(30,24)=(z1 + z3 + z4)**2/9.
      q(30,25)=(z1 + z3 + z4)**3/27.
      q(30,26)=(z1 + z3 + z4)**4/81.
      q(30,27)=(y1 + y3 + y4)/3.
      q(30,28)=((y1 + y3 + y4)*(z1 + z3 + z4))/9.
      q(30,29)=((y1 + y3 + y4)*(z1 + z3 + z4)**2)/27.
      q(30,30)=((y1 + y3 + y4)*(z1 + z3 + z4)**3)/81.
      q(30,31)=(y1 + y3 + y4)**2/9.
      q(30,32)=((y1 + y3 + y4)**2*(z1 + z3 + z4))/27.
      q(30,33)=((y1 + y3 + y4)**2*(z1 + z3 + z4)**2)/81.
      q(30,34)=(y1 + y3 + y4)**3/27.
      q(30,35)=((y1 + y3 + y4)**3*(z1 + z3 + z4))/81.
      q(30,36)=(y1 + y3 + y4)**4/81.
      q(30,37)=(2*(x1 + x3 + x4))/3.
      q(30,38)=(2*(x1 + x3 + x4)*(z1 + z3 + z4))/9.
      q(30,39)=(2*(x1 + x3 + x4)*(z1 + z3 + z4)**2)/27.
      q(30,40)=(2*(x1 + x3 + x4)*(z1 + z3 + z4)**3)/81.
      q(30,41)=(2*(x1 + x3 + x4)*(y1 + y3 + y4))/9.
      q(30,42)=(2*(x1 + x3 + x4)*(y1 + y3 + y4)*(z1 + z3 + z4))/27.
      q(30,43)=(2*(x1 + x3 + x4)*(y1 + y3 + y4)*(z1 + z3 + z4)**2)/81.
      q(30,44)=(2*(x1 + x3 + x4)*(y1 + y3 + y4)**2)/27.
      q(30,45)=(2*(x1 + x3 + x4)*(y1 + y3 + y4)**2*(z1 + z3 + z4))/81.
      q(30,46)=(2*(x1 + x3 + x4)*(y1 + y3 + y4)**3)/81.
      q(30,47)=(x1 + x3 + x4)**2/3.
      q(30,48)=((x1 + x3 + x4)**2*(z1 + z3 + z4))/9.
      q(30,49)=((x1 + x3 + x4)**2*(z1 + z3 + z4)**2)/27.
      q(30,50)=((x1 + x3 + x4)**2*(y1 + y3 + y4))/9.
      q(30,51)=((x1 + x3 + x4)**2*(y1 + y3 + y4)*(z1 + z3 + z4))/27.
      q(30,52)=((x1 + x3 + x4)**2*(y1 + y3 + y4)**2)/27.
      q(30,53)=(4*(x1 + x3 + x4)**3)/27.
      q(30,54)=(4*(x1 + x3 + x4)**3*(z1 + z3 + z4))/81.
      q(30,55)=(4*(x1 + x3 + x4)**3*(y1 + y3 + y4))/81.
      q(30,56)=(5*(x1 + x3 + x4)**4)/81.
      q(31,1)=0
      q(31,2)=0
      q(31,3)=0
      q(31,4)=0
      q(31,5)=0
      q(31,6)=0
      q(31,7)=1
      q(31,8)=(z1 + z3 + z4)/3.
      q(31,9)=(z1 + z3 + z4)**2/9.
      q(31,10)=(z1 + z3 + z4)**3/27.
      q(31,11)=(z1 + z3 + z4)**4/81.
      q(31,12)=(2*(y1 + y3 + y4))/3.
      q(31,13)=(2*(y1 + y3 + y4)*(z1 + z3 + z4))/9.
      q(31,14)=(2*(y1 + y3 + y4)*(z1 + z3 + z4)**2)/27.
      q(31,15)=(2*(y1 + y3 + y4)*(z1 + z3 + z4)**3)/81.
      q(31,16)=(y1 + y3 + y4)**2/3.
      q(31,17)=((y1 + y3 + y4)**2*(z1 + z3 + z4))/9.
      q(31,18)=((y1 + y3 + y4)**2*(z1 + z3 + z4)**2)/27.
      q(31,19)=(4*(y1 + y3 + y4)**3)/27.
      q(31,20)=(4*(y1 + y3 + y4)**3*(z1 + z3 + z4))/81.
      q(31,21)=(5*(y1 + y3 + y4)**4)/81.
      q(31,22)=0
      q(31,23)=0
      q(31,24)=0
      q(31,25)=0
      q(31,26)=0
      q(31,27)=(x1 + x3 + x4)/3.
      q(31,28)=((x1 + x3 + x4)*(z1 + z3 + z4))/9.
      q(31,29)=((x1 + x3 + x4)*(z1 + z3 + z4)**2)/27.
      q(31,30)=((x1 + x3 + x4)*(z1 + z3 + z4)**3)/81.
      q(31,31)=(2*(x1 + x3 + x4)*(y1 + y3 + y4))/9.
      q(31,32)=(2*(x1 + x3 + x4)*(y1 + y3 + y4)*(z1 + z3 + z4))/27.
      q(31,33)=(2*(x1 + x3 + x4)*(y1 + y3 + y4)*(z1 + z3 + z4)**2)/81.
      q(31,34)=((x1 + x3 + x4)*(y1 + y3 + y4)**2)/9.
      q(31,35)=((x1 + x3 + x4)*(y1 + y3 + y4)**2*(z1 + z3 + z4))/27.
      q(31,36)=(4*(x1 + x3 + x4)*(y1 + y3 + y4)**3)/81.
      q(31,37)=0
      q(31,38)=0
      q(31,39)=0
      q(31,40)=0
      q(31,41)=(x1 + x3 + x4)**2/9.
      q(31,42)=((x1 + x3 + x4)**2*(z1 + z3 + z4))/27.
      q(31,43)=((x1 + x3 + x4)**2*(z1 + z3 + z4)**2)/81.
      q(31,44)=(2*(x1 + x3 + x4)**2*(y1 + y3 + y4))/27.
      q(31,45)=(2*(x1 + x3 + x4)**2*(y1 + y3 + y4)*(z1 + z3 + z4))/81.
      q(31,46)=((x1 + x3 + x4)**2*(y1 + y3 + y4)**2)/27.
      q(31,47)=0
      q(31,48)=0
      q(31,49)=0
      q(31,50)=(x1 + x3 + x4)**3/27.
      q(31,51)=((x1 + x3 + x4)**3*(z1 + z3 + z4))/81.
      q(31,52)=(2*(x1 + x3 + x4)**3*(y1 + y3 + y4))/81.
      q(31,53)=0
      q(31,54)=0
      q(31,55)=(x1 + x3 + x4)**4/81.
      q(31,56)=0
      q(32,1)=0
      q(32,2)=1
      q(32,3)=(2*(z1 + z3 + z4))/3.
      q(32,4)=(z1 + z3 + z4)**2/3.
      q(32,5)=(4*(z1 + z3 + z4)**3)/27.
      q(32,6)=(5*(z1 + z3 + z4)**4)/81.
      q(32,7)=0
      q(32,8)=(y1 + y3 + y4)/3.
      q(32,9)=(2*(y1 + y3 + y4)*(z1 + z3 + z4))/9.
      q(32,10)=((y1 + y3 + y4)*(z1 + z3 + z4)**2)/9.
      q(32,11)=(4*(y1 + y3 + y4)*(z1 + z3 + z4)**3)/81.
      q(32,12)=0
      q(32,13)=(y1 + y3 + y4)**2/9.
      q(32,14)=(2*(y1 + y3 + y4)**2*(z1 + z3 + z4))/27.
      q(32,15)=((y1 + y3 + y4)**2*(z1 + z3 + z4)**2)/27.
      q(32,16)=0
      q(32,17)=(y1 + y3 + y4)**3/27.
      q(32,18)=(2*(y1 + y3 + y4)**3*(z1 + z3 + z4))/81.
      q(32,19)=0
      q(32,20)=(y1 + y3 + y4)**4/81.
      q(32,21)=0
      q(32,22)=0
      q(32,23)=(x1 + x3 + x4)/3.
      q(32,24)=(2*(x1 + x3 + x4)*(z1 + z3 + z4))/9.
      q(32,25)=((x1 + x3 + x4)*(z1 + z3 + z4)**2)/9.
      q(32,26)=(4*(x1 + x3 + x4)*(z1 + z3 + z4)**3)/81.
      q(32,27)=0
      q(32,28)=((x1 + x3 + x4)*(y1 + y3 + y4))/9.
      q(32,29)=(2*(x1 + x3 + x4)*(y1 + y3 + y4)*(z1 + z3 + z4))/27.
      q(32,30)=((x1 + x3 + x4)*(y1 + y3 + y4)*(z1 + z3 + z4)**2)/27.
      q(32,31)=0
      q(32,32)=((x1 + x3 + x4)*(y1 + y3 + y4)**2)/27.
      q(32,33)=(2*(x1 + x3 + x4)*(y1 + y3 + y4)**2*(z1 + z3 + z4))/81.
      q(32,34)=0
      q(32,35)=((x1 + x3 + x4)*(y1 + y3 + y4)**3)/81.
      q(32,36)=0
      q(32,37)=0
      q(32,38)=(x1 + x3 + x4)**2/9.
      q(32,39)=(2*(x1 + x3 + x4)**2*(z1 + z3 + z4))/27.
      q(32,40)=((x1 + x3 + x4)**2*(z1 + z3 + z4)**2)/27.
      q(32,41)=0
      q(32,42)=((x1 + x3 + x4)**2*(y1 + y3 + y4))/27.
      q(32,43)=(2*(x1 + x3 + x4)**2*(y1 + y3 + y4)*(z1 + z3 + z4))/81.
      q(32,44)=0
      q(32,45)=((x1 + x3 + x4)**2*(y1 + y3 + y4)**2)/81.
      q(32,46)=0
      q(32,47)=0
      q(32,48)=(x1 + x3 + x4)**3/27.
      q(32,49)=(2*(x1 + x3 + x4)**3*(z1 + z3 + z4))/81.
      q(32,50)=0
      q(32,51)=((x1 + x3 + x4)**3*(y1 + y3 + y4))/81.
      q(32,52)=0
      q(32,53)=0
      q(32,54)=(x1 + x3 + x4)**4/81.
      q(32,55)=0
      q(32,56)=0
      q(33,1)=1
      q(33,2)=(z1 + z2)/2.
      q(33,3)=(z1 + z2)**2/4.
      q(33,4)=(z1 + z2)**3/8.
      q(33,5)=(z1 + z2)**4/16.
      q(33,6)=(z1 + z2)**5/32.
      q(33,7)=(y1 + y2)/2.
      q(33,8)=((y1 + y2)*(z1 + z2))/4.
      q(33,9)=((y1 + y2)*(z1 + z2)**2)/8.
      q(33,10)=((y1 + y2)*(z1 + z2)**3)/16.
      q(33,11)=((y1 + y2)*(z1 + z2)**4)/32.
      q(33,12)=(y1 + y2)**2/4.
      q(33,13)=((y1 + y2)**2*(z1 + z2))/8.
      q(33,14)=((y1 + y2)**2*(z1 + z2)**2)/16.
      q(33,15)=((y1 + y2)**2*(z1 + z2)**3)/32.
      q(33,16)=(y1 + y2)**3/8.
      q(33,17)=((y1 + y2)**3*(z1 + z2))/16.
      q(33,18)=((y1 + y2)**3*(z1 + z2)**2)/32.
      q(33,19)=(y1 + y2)**4/16.
      q(33,20)=((y1 + y2)**4*(z1 + z2))/32.
      q(33,21)=(y1 + y2)**5/32.
      q(33,22)=(x1 + x2)/2.
      q(33,23)=((x1 + x2)*(z1 + z2))/4.
      q(33,24)=((x1 + x2)*(z1 + z2)**2)/8.
      q(33,25)=((x1 + x2)*(z1 + z2)**3)/16.
      q(33,26)=((x1 + x2)*(z1 + z2)**4)/32.
      q(33,27)=((x1 + x2)*(y1 + y2))/4.
      q(33,28)=((x1 + x2)*(y1 + y2)*(z1 + z2))/8.
      q(33,29)=((x1 + x2)*(y1 + y2)*(z1 + z2)**2)/16.
      q(33,30)=((x1 + x2)*(y1 + y2)*(z1 + z2)**3)/32.
      q(33,31)=((x1 + x2)*(y1 + y2)**2)/8.
      q(33,32)=((x1 + x2)*(y1 + y2)**2*(z1 + z2))/16.
      q(33,33)=((x1 + x2)*(y1 + y2)**2*(z1 + z2)**2)/32.
      q(33,34)=((x1 + x2)*(y1 + y2)**3)/16.
      q(33,35)=((x1 + x2)*(y1 + y2)**3*(z1 + z2))/32.
      q(33,36)=((x1 + x2)*(y1 + y2)**4)/32.
      q(33,37)=(x1 + x2)**2/4.
      q(33,38)=((x1 + x2)**2*(z1 + z2))/8.
      q(33,39)=((x1 + x2)**2*(z1 + z2)**2)/16.
      q(33,40)=((x1 + x2)**2*(z1 + z2)**3)/32.
      q(33,41)=((x1 + x2)**2*(y1 + y2))/8.
      q(33,42)=((x1 + x2)**2*(y1 + y2)*(z1 + z2))/16.
      q(33,43)=((x1 + x2)**2*(y1 + y2)*(z1 + z2)**2)/32.
      q(33,44)=((x1 + x2)**2*(y1 + y2)**2)/16.
      q(33,45)=((x1 + x2)**2*(y1 + y2)**2*(z1 + z2))/32.
      q(33,46)=((x1 + x2)**2*(y1 + y2)**3)/32.
      q(33,47)=(x1 + x2)**3/8.
      q(33,48)=((x1 + x2)**3*(z1 + z2))/16.
      q(33,49)=((x1 + x2)**3*(z1 + z2)**2)/32.
      q(33,50)=((x1 + x2)**3*(y1 + y2))/16.
      q(33,51)=((x1 + x2)**3*(y1 + y2)*(z1 + z2))/32.
      q(33,52)=((x1 + x2)**3*(y1 + y2)**2)/32.
      q(33,53)=(x1 + x2)**4/16.
      q(33,54)=((x1 + x2)**4*(z1 + z2))/32.
      q(33,55)=((x1 + x2)**4*(y1 + y2))/32.
      q(33,56)=(x1 + x2)**5/32.
      q(34,1)=1
      q(34,2)=(z1 + z3)/2.
      q(34,3)=(z1 + z3)**2/4.
      q(34,4)=(z1 + z3)**3/8.
      q(34,5)=(z1 + z3)**4/16.
      q(34,6)=(z1 + z3)**5/32.
      q(34,7)=(y1 + y3)/2.
      q(34,8)=((y1 + y3)*(z1 + z3))/4.
      q(34,9)=((y1 + y3)*(z1 + z3)**2)/8.
      q(34,10)=((y1 + y3)*(z1 + z3)**3)/16.
      q(34,11)=((y1 + y3)*(z1 + z3)**4)/32.
      q(34,12)=(y1 + y3)**2/4.
      q(34,13)=((y1 + y3)**2*(z1 + z3))/8.
      q(34,14)=((y1 + y3)**2*(z1 + z3)**2)/16.
      q(34,15)=((y1 + y3)**2*(z1 + z3)**3)/32.
      q(34,16)=(y1 + y3)**3/8.
      q(34,17)=((y1 + y3)**3*(z1 + z3))/16.
      q(34,18)=((y1 + y3)**3*(z1 + z3)**2)/32.
      q(34,19)=(y1 + y3)**4/16.
      q(34,20)=((y1 + y3)**4*(z1 + z3))/32.
      q(34,21)=(y1 + y3)**5/32.
      q(34,22)=(x1 + x3)/2.
      q(34,23)=((x1 + x3)*(z1 + z3))/4.
      q(34,24)=((x1 + x3)*(z1 + z3)**2)/8.
      q(34,25)=((x1 + x3)*(z1 + z3)**3)/16.
      q(34,26)=((x1 + x3)*(z1 + z3)**4)/32.
      q(34,27)=((x1 + x3)*(y1 + y3))/4.
      q(34,28)=((x1 + x3)*(y1 + y3)*(z1 + z3))/8.
      q(34,29)=((x1 + x3)*(y1 + y3)*(z1 + z3)**2)/16.
      q(34,30)=((x1 + x3)*(y1 + y3)*(z1 + z3)**3)/32.
      q(34,31)=((x1 + x3)*(y1 + y3)**2)/8.
      q(34,32)=((x1 + x3)*(y1 + y3)**2*(z1 + z3))/16.
      q(34,33)=((x1 + x3)*(y1 + y3)**2*(z1 + z3)**2)/32.
      q(34,34)=((x1 + x3)*(y1 + y3)**3)/16.
      q(34,35)=((x1 + x3)*(y1 + y3)**3*(z1 + z3))/32.
      q(34,36)=((x1 + x3)*(y1 + y3)**4)/32.
      q(34,37)=(x1 + x3)**2/4.
      q(34,38)=((x1 + x3)**2*(z1 + z3))/8.
      q(34,39)=((x1 + x3)**2*(z1 + z3)**2)/16.
      q(34,40)=((x1 + x3)**2*(z1 + z3)**3)/32.
      q(34,41)=((x1 + x3)**2*(y1 + y3))/8.
      q(34,42)=((x1 + x3)**2*(y1 + y3)*(z1 + z3))/16.
      q(34,43)=((x1 + x3)**2*(y1 + y3)*(z1 + z3)**2)/32.
      q(34,44)=((x1 + x3)**2*(y1 + y3)**2)/16.
      q(34,45)=((x1 + x3)**2*(y1 + y3)**2*(z1 + z3))/32.
      q(34,46)=((x1 + x3)**2*(y1 + y3)**3)/32.
      q(34,47)=(x1 + x3)**3/8.
      q(34,48)=((x1 + x3)**3*(z1 + z3))/16.
      q(34,49)=((x1 + x3)**3*(z1 + z3)**2)/32.
      q(34,50)=((x1 + x3)**3*(y1 + y3))/16.
      q(34,51)=((x1 + x3)**3*(y1 + y3)*(z1 + z3))/32.
      q(34,52)=((x1 + x3)**3*(y1 + y3)**2)/32.
      q(34,53)=(x1 + x3)**4/16.
      q(34,54)=((x1 + x3)**4*(z1 + z3))/32.
      q(34,55)=((x1 + x3)**4*(y1 + y3))/32.
      q(34,56)=(x1 + x3)**5/32.
      q(35,1)=1
      q(35,2)=(z1 + z4)/2.
      q(35,3)=(z1 + z4)**2/4.
      q(35,4)=(z1 + z4)**3/8.
      q(35,5)=(z1 + z4)**4/16.
      q(35,6)=(z1 + z4)**5/32.
      q(35,7)=(y1 + y4)/2.
      q(35,8)=((y1 + y4)*(z1 + z4))/4.
      q(35,9)=((y1 + y4)*(z1 + z4)**2)/8.
      q(35,10)=((y1 + y4)*(z1 + z4)**3)/16.
      q(35,11)=((y1 + y4)*(z1 + z4)**4)/32.
      q(35,12)=(y1 + y4)**2/4.
      q(35,13)=((y1 + y4)**2*(z1 + z4))/8.
      q(35,14)=((y1 + y4)**2*(z1 + z4)**2)/16.
      q(35,15)=((y1 + y4)**2*(z1 + z4)**3)/32.
      q(35,16)=(y1 + y4)**3/8.
      q(35,17)=((y1 + y4)**3*(z1 + z4))/16.
      q(35,18)=((y1 + y4)**3*(z1 + z4)**2)/32.
      q(35,19)=(y1 + y4)**4/16.
      q(35,20)=((y1 + y4)**4*(z1 + z4))/32.
      q(35,21)=(y1 + y4)**5/32.
      q(35,22)=(x1 + x4)/2.
      q(35,23)=((x1 + x4)*(z1 + z4))/4.
      q(35,24)=((x1 + x4)*(z1 + z4)**2)/8.
      q(35,25)=((x1 + x4)*(z1 + z4)**3)/16.
      q(35,26)=((x1 + x4)*(z1 + z4)**4)/32.
      q(35,27)=((x1 + x4)*(y1 + y4))/4.
      q(35,28)=((x1 + x4)*(y1 + y4)*(z1 + z4))/8.
      q(35,29)=((x1 + x4)*(y1 + y4)*(z1 + z4)**2)/16.
      q(35,30)=((x1 + x4)*(y1 + y4)*(z1 + z4)**3)/32.
      q(35,31)=((x1 + x4)*(y1 + y4)**2)/8.
      q(35,32)=((x1 + x4)*(y1 + y4)**2*(z1 + z4))/16.
      q(35,33)=((x1 + x4)*(y1 + y4)**2*(z1 + z4)**2)/32.
      q(35,34)=((x1 + x4)*(y1 + y4)**3)/16.
      q(35,35)=((x1 + x4)*(y1 + y4)**3*(z1 + z4))/32.
      q(35,36)=((x1 + x4)*(y1 + y4)**4)/32.
      q(35,37)=(x1 + x4)**2/4.
      q(35,38)=((x1 + x4)**2*(z1 + z4))/8.
      q(35,39)=((x1 + x4)**2*(z1 + z4)**2)/16.
      q(35,40)=((x1 + x4)**2*(z1 + z4)**3)/32.
      q(35,41)=((x1 + x4)**2*(y1 + y4))/8.
      q(35,42)=((x1 + x4)**2*(y1 + y4)*(z1 + z4))/16.
      q(35,43)=((x1 + x4)**2*(y1 + y4)*(z1 + z4)**2)/32.
      q(35,44)=((x1 + x4)**2*(y1 + y4)**2)/16.
      q(35,45)=((x1 + x4)**2*(y1 + y4)**2*(z1 + z4))/32.
      q(35,46)=((x1 + x4)**2*(y1 + y4)**3)/32.
      q(35,47)=(x1 + x4)**3/8.
      q(35,48)=((x1 + x4)**3*(z1 + z4))/16.
      q(35,49)=((x1 + x4)**3*(z1 + z4)**2)/32.
      q(35,50)=((x1 + x4)**3*(y1 + y4))/16.
      q(35,51)=((x1 + x4)**3*(y1 + y4)*(z1 + z4))/32.
      q(35,52)=((x1 + x4)**3*(y1 + y4)**2)/32.
      q(35,53)=(x1 + x4)**4/16.
      q(35,54)=((x1 + x4)**4*(z1 + z4))/32.
      q(35,55)=((x1 + x4)**4*(y1 + y4))/32.
      q(35,56)=(x1 + x4)**5/32.
      q(36,1)=1
      q(36,2)=(z2 + z3)/2.
      q(36,3)=(z2 + z3)**2/4.
      q(36,4)=(z2 + z3)**3/8.
      q(36,5)=(z2 + z3)**4/16.
      q(36,6)=(z2 + z3)**5/32.
      q(36,7)=(y2 + y3)/2.
      q(36,8)=((y2 + y3)*(z2 + z3))/4.
      q(36,9)=((y2 + y3)*(z2 + z3)**2)/8.
      q(36,10)=((y2 + y3)*(z2 + z3)**3)/16.
      q(36,11)=((y2 + y3)*(z2 + z3)**4)/32.
      q(36,12)=(y2 + y3)**2/4.
      q(36,13)=((y2 + y3)**2*(z2 + z3))/8.
      q(36,14)=((y2 + y3)**2*(z2 + z3)**2)/16.
      q(36,15)=((y2 + y3)**2*(z2 + z3)**3)/32.
      q(36,16)=(y2 + y3)**3/8.
      q(36,17)=((y2 + y3)**3*(z2 + z3))/16.
      q(36,18)=((y2 + y3)**3*(z2 + z3)**2)/32.
      q(36,19)=(y2 + y3)**4/16.
      q(36,20)=((y2 + y3)**4*(z2 + z3))/32.
      q(36,21)=(y2 + y3)**5/32.
      q(36,22)=(x2 + x3)/2.
      q(36,23)=((x2 + x3)*(z2 + z3))/4.
      q(36,24)=((x2 + x3)*(z2 + z3)**2)/8.
      q(36,25)=((x2 + x3)*(z2 + z3)**3)/16.
      q(36,26)=((x2 + x3)*(z2 + z3)**4)/32.
      q(36,27)=((x2 + x3)*(y2 + y3))/4.
      q(36,28)=((x2 + x3)*(y2 + y3)*(z2 + z3))/8.
      q(36,29)=((x2 + x3)*(y2 + y3)*(z2 + z3)**2)/16.
      q(36,30)=((x2 + x3)*(y2 + y3)*(z2 + z3)**3)/32.
      q(36,31)=((x2 + x3)*(y2 + y3)**2)/8.
      q(36,32)=((x2 + x3)*(y2 + y3)**2*(z2 + z3))/16.
      q(36,33)=((x2 + x3)*(y2 + y3)**2*(z2 + z3)**2)/32.
      q(36,34)=((x2 + x3)*(y2 + y3)**3)/16.
      q(36,35)=((x2 + x3)*(y2 + y3)**3*(z2 + z3))/32.
      q(36,36)=((x2 + x3)*(y2 + y3)**4)/32.
      q(36,37)=(x2 + x3)**2/4.
      q(36,38)=((x2 + x3)**2*(z2 + z3))/8.
      q(36,39)=((x2 + x3)**2*(z2 + z3)**2)/16.
      q(36,40)=((x2 + x3)**2*(z2 + z3)**3)/32.
      q(36,41)=((x2 + x3)**2*(y2 + y3))/8.
      q(36,42)=((x2 + x3)**2*(y2 + y3)*(z2 + z3))/16.
      q(36,43)=((x2 + x3)**2*(y2 + y3)*(z2 + z3)**2)/32.
      q(36,44)=((x2 + x3)**2*(y2 + y3)**2)/16.
      q(36,45)=((x2 + x3)**2*(y2 + y3)**2*(z2 + z3))/32.
      q(36,46)=((x2 + x3)**2*(y2 + y3)**3)/32.
      q(36,47)=(x2 + x3)**3/8.
      q(36,48)=((x2 + x3)**3*(z2 + z3))/16.
      q(36,49)=((x2 + x3)**3*(z2 + z3)**2)/32.
      q(36,50)=((x2 + x3)**3*(y2 + y3))/16.
      q(36,51)=((x2 + x3)**3*(y2 + y3)*(z2 + z3))/32.
      q(36,52)=((x2 + x3)**3*(y2 + y3)**2)/32.
      q(36,53)=(x2 + x3)**4/16.
      q(36,54)=((x2 + x3)**4*(z2 + z3))/32.
      q(36,55)=((x2 + x3)**4*(y2 + y3))/32.
      q(36,56)=(x2 + x3)**5/32.
      q(37,1)=1
      q(37,2)=(z2 + z4)/2.
      q(37,3)=(z2 + z4)**2/4.
      q(37,4)=(z2 + z4)**3/8.
      q(37,5)=(z2 + z4)**4/16.
      q(37,6)=(z2 + z4)**5/32.
      q(37,7)=(y2 + y4)/2.
      q(37,8)=((y2 + y4)*(z2 + z4))/4.
      q(37,9)=((y2 + y4)*(z2 + z4)**2)/8.
      q(37,10)=((y2 + y4)*(z2 + z4)**3)/16.
      q(37,11)=((y2 + y4)*(z2 + z4)**4)/32.
      q(37,12)=(y2 + y4)**2/4.
      q(37,13)=((y2 + y4)**2*(z2 + z4))/8.
      q(37,14)=((y2 + y4)**2*(z2 + z4)**2)/16.
      q(37,15)=((y2 + y4)**2*(z2 + z4)**3)/32.
      q(37,16)=(y2 + y4)**3/8.
      q(37,17)=((y2 + y4)**3*(z2 + z4))/16.
      q(37,18)=((y2 + y4)**3*(z2 + z4)**2)/32.
      q(37,19)=(y2 + y4)**4/16.
      q(37,20)=((y2 + y4)**4*(z2 + z4))/32.
      q(37,21)=(y2 + y4)**5/32.
      q(37,22)=(x2 + x4)/2.
      q(37,23)=((x2 + x4)*(z2 + z4))/4.
      q(37,24)=((x2 + x4)*(z2 + z4)**2)/8.
      q(37,25)=((x2 + x4)*(z2 + z4)**3)/16.
      q(37,26)=((x2 + x4)*(z2 + z4)**4)/32.
      q(37,27)=((x2 + x4)*(y2 + y4))/4.
      q(37,28)=((x2 + x4)*(y2 + y4)*(z2 + z4))/8.
      q(37,29)=((x2 + x4)*(y2 + y4)*(z2 + z4)**2)/16.
      q(37,30)=((x2 + x4)*(y2 + y4)*(z2 + z4)**3)/32.
      q(37,31)=((x2 + x4)*(y2 + y4)**2)/8.
      q(37,32)=((x2 + x4)*(y2 + y4)**2*(z2 + z4))/16.
      q(37,33)=((x2 + x4)*(y2 + y4)**2*(z2 + z4)**2)/32.
      q(37,34)=((x2 + x4)*(y2 + y4)**3)/16.
      q(37,35)=((x2 + x4)*(y2 + y4)**3*(z2 + z4))/32.
      q(37,36)=((x2 + x4)*(y2 + y4)**4)/32.
      q(37,37)=(x2 + x4)**2/4.
      q(37,38)=((x2 + x4)**2*(z2 + z4))/8.
      q(37,39)=((x2 + x4)**2*(z2 + z4)**2)/16.
      q(37,40)=((x2 + x4)**2*(z2 + z4)**3)/32.
      q(37,41)=((x2 + x4)**2*(y2 + y4))/8.
      q(37,42)=((x2 + x4)**2*(y2 + y4)*(z2 + z4))/16.
      q(37,43)=((x2 + x4)**2*(y2 + y4)*(z2 + z4)**2)/32.
      q(37,44)=((x2 + x4)**2*(y2 + y4)**2)/16.
      q(37,45)=((x2 + x4)**2*(y2 + y4)**2*(z2 + z4))/32.
      q(37,46)=((x2 + x4)**2*(y2 + y4)**3)/32.
      q(37,47)=(x2 + x4)**3/8.
      q(37,48)=((x2 + x4)**3*(z2 + z4))/16.
      q(37,49)=((x2 + x4)**3*(z2 + z4)**2)/32.
      q(37,50)=((x2 + x4)**3*(y2 + y4))/16.
      q(37,51)=((x2 + x4)**3*(y2 + y4)*(z2 + z4))/32.
      q(37,52)=((x2 + x4)**3*(y2 + y4)**2)/32.
      q(37,53)=(x2 + x4)**4/16.
      q(37,54)=((x2 + x4)**4*(z2 + z4))/32.
      q(37,55)=((x2 + x4)**4*(y2 + y4))/32.
      q(37,56)=(x2 + x4)**5/32.
      q(38,1)=1
      q(38,2)=(z3 + z4)/2.
      q(38,3)=(z3 + z4)**2/4.
      q(38,4)=(z3 + z4)**3/8.
      q(38,5)=(z3 + z4)**4/16.
      q(38,6)=(z3 + z4)**5/32.
      q(38,7)=(y3 + y4)/2.
      q(38,8)=((y3 + y4)*(z3 + z4))/4.
      q(38,9)=((y3 + y4)*(z3 + z4)**2)/8.
      q(38,10)=((y3 + y4)*(z3 + z4)**3)/16.
      q(38,11)=((y3 + y4)*(z3 + z4)**4)/32.
      q(38,12)=(y3 + y4)**2/4.
      q(38,13)=((y3 + y4)**2*(z3 + z4))/8.
      q(38,14)=((y3 + y4)**2*(z3 + z4)**2)/16.
      q(38,15)=((y3 + y4)**2*(z3 + z4)**3)/32.
      q(38,16)=(y3 + y4)**3/8.
      q(38,17)=((y3 + y4)**3*(z3 + z4))/16.
      q(38,18)=((y3 + y4)**3*(z3 + z4)**2)/32.
      q(38,19)=(y3 + y4)**4/16.
      q(38,20)=((y3 + y4)**4*(z3 + z4))/32.
      q(38,21)=(y3 + y4)**5/32.
      q(38,22)=(x3 + x4)/2.
      q(38,23)=((x3 + x4)*(z3 + z4))/4.
      q(38,24)=((x3 + x4)*(z3 + z4)**2)/8.
      q(38,25)=((x3 + x4)*(z3 + z4)**3)/16.
      q(38,26)=((x3 + x4)*(z3 + z4)**4)/32.
      q(38,27)=((x3 + x4)*(y3 + y4))/4.
      q(38,28)=((x3 + x4)*(y3 + y4)*(z3 + z4))/8.
      q(38,29)=((x3 + x4)*(y3 + y4)*(z3 + z4)**2)/16.
      q(38,30)=((x3 + x4)*(y3 + y4)*(z3 + z4)**3)/32.
      q(38,31)=((x3 + x4)*(y3 + y4)**2)/8.
      q(38,32)=((x3 + x4)*(y3 + y4)**2*(z3 + z4))/16.
      q(38,33)=((x3 + x4)*(y3 + y4)**2*(z3 + z4)**2)/32.
      q(38,34)=((x3 + x4)*(y3 + y4)**3)/16.
      q(38,35)=((x3 + x4)*(y3 + y4)**3*(z3 + z4))/32.
      q(38,36)=((x3 + x4)*(y3 + y4)**4)/32.
      q(38,37)=(x3 + x4)**2/4.
      q(38,38)=((x3 + x4)**2*(z3 + z4))/8.
      q(38,39)=((x3 + x4)**2*(z3 + z4)**2)/16.
      q(38,40)=((x3 + x4)**2*(z3 + z4)**3)/32.
      q(38,41)=((x3 + x4)**2*(y3 + y4))/8.
      q(38,42)=((x3 + x4)**2*(y3 + y4)*(z3 + z4))/16.
      q(38,43)=((x3 + x4)**2*(y3 + y4)*(z3 + z4)**2)/32.
      q(38,44)=((x3 + x4)**2*(y3 + y4)**2)/16.
      q(38,45)=((x3 + x4)**2*(y3 + y4)**2*(z3 + z4))/32.
      q(38,46)=((x3 + x4)**2*(y3 + y4)**3)/32.
      q(38,47)=(x3 + x4)**3/8.
      q(38,48)=((x3 + x4)**3*(z3 + z4))/16.
      q(38,49)=((x3 + x4)**3*(z3 + z4)**2)/32.
      q(38,50)=((x3 + x4)**3*(y3 + y4))/16.
      q(38,51)=((x3 + x4)**3*(y3 + y4)*(z3 + z4))/32.
      q(38,52)=((x3 + x4)**3*(y3 + y4)**2)/32.
      q(38,53)=(x3 + x4)**4/16.
      q(38,54)=((x3 + x4)**4*(z3 + z4))/32.
      q(38,55)=((x3 + x4)**4*(y3 + y4))/32.
      q(38,56)=(x3 + x4)**5/32.
      q(39,1)=0
      q(39,2)=0
      q(39,3)=0
      q(39,4)=0
      q(39,5)=0
      q(39,6)=0
      q(39,7)=0
      q(39,8)=0
      q(39,9)=0
      q(39,10)=0
      q(39,11)=0
      q(39,12)=0
      q(39,13)=0
      q(39,14)=0
      q(39,15)=0
      q(39,16)=0
      q(39,17)=0
      q(39,18)=0
      q(39,19)=0
      q(39,20)=0
      q(39,21)=0
      q(39,22)=1
      q(39,23)=(z1 + z2)/2.
      q(39,24)=(z1 + z2)**2/4.
      q(39,25)=(z1 + z2)**3/8.
      q(39,26)=(z1 + z2)**4/16.
      q(39,27)=(y1 + y2)/2.
      q(39,28)=((y1 + y2)*(z1 + z2))/4.
      q(39,29)=((y1 + y2)*(z1 + z2)**2)/8.
      q(39,30)=((y1 + y2)*(z1 + z2)**3)/16.
      q(39,31)=(y1 + y2)**2/4.
      q(39,32)=((y1 + y2)**2*(z1 + z2))/8.
      q(39,33)=((y1 + y2)**2*(z1 + z2)**2)/16.
      q(39,34)=(y1 + y2)**3/8.
      q(39,35)=((y1 + y2)**3*(z1 + z2))/16.
      q(39,36)=(y1 + y2)**4/16.
      q(39,37)=x1 + x2
      q(39,38)=((x1 + x2)*(z1 + z2))/2.
      q(39,39)=((x1 + x2)*(z1 + z2)**2)/4.
      q(39,40)=((x1 + x2)*(z1 + z2)**3)/8.
      q(39,41)=((x1 + x2)*(y1 + y2))/2.
      q(39,42)=((x1 + x2)*(y1 + y2)*(z1 + z2))/4.
      q(39,43)=((x1 + x2)*(y1 + y2)*(z1 + z2)**2)/8.
      q(39,44)=((x1 + x2)*(y1 + y2)**2)/4.
      q(39,45)=((x1 + x2)*(y1 + y2)**2*(z1 + z2))/8.
      q(39,46)=((x1 + x2)*(y1 + y2)**3)/8.
      q(39,47)=(3*(x1 + x2)**2)/4.
      q(39,48)=(3*(x1 + x2)**2*(z1 + z2))/8.
      q(39,49)=(3*(x1 + x2)**2*(z1 + z2)**2)/16.
      q(39,50)=(3*(x1 + x2)**2*(y1 + y2))/8.
      q(39,51)=(3*(x1 + x2)**2*(y1 + y2)*(z1 + z2))/16.
      q(39,52)=(3*(x1 + x2)**2*(y1 + y2)**2)/16.
      q(39,53)=(x1 + x2)**3/2.
      q(39,54)=((x1 + x2)**3*(z1 + z2))/4.
      q(39,55)=((x1 + x2)**3*(y1 + y2))/4.
      q(39,56)=(5*(x1 + x2)**4)/16.
      q(40,1)=0
      q(40,2)=0
      q(40,3)=0
      q(40,4)=0
      q(40,5)=0
      q(40,6)=0
      q(40,7)=1
      q(40,8)=(z1 + z2)/2.
      q(40,9)=(z1 + z2)**2/4.
      q(40,10)=(z1 + z2)**3/8.
      q(40,11)=(z1 + z2)**4/16.
      q(40,12)=y1 + y2
      q(40,13)=((y1 + y2)*(z1 + z2))/2.
      q(40,14)=((y1 + y2)*(z1 + z2)**2)/4.
      q(40,15)=((y1 + y2)*(z1 + z2)**3)/8.
      q(40,16)=(3*(y1 + y2)**2)/4.
      q(40,17)=(3*(y1 + y2)**2*(z1 + z2))/8.
      q(40,18)=(3*(y1 + y2)**2*(z1 + z2)**2)/16.
      q(40,19)=(y1 + y2)**3/2.
      q(40,20)=((y1 + y2)**3*(z1 + z2))/4.
      q(40,21)=(5*(y1 + y2)**4)/16.
      q(40,22)=0
      q(40,23)=0
      q(40,24)=0
      q(40,25)=0
      q(40,26)=0
      q(40,27)=(x1 + x2)/2.
      q(40,28)=((x1 + x2)*(z1 + z2))/4.
      q(40,29)=((x1 + x2)*(z1 + z2)**2)/8.
      q(40,30)=((x1 + x2)*(z1 + z2)**3)/16.
      q(40,31)=((x1 + x2)*(y1 + y2))/2.
      q(40,32)=((x1 + x2)*(y1 + y2)*(z1 + z2))/4.
      q(40,33)=((x1 + x2)*(y1 + y2)*(z1 + z2)**2)/8.
      q(40,34)=(3*(x1 + x2)*(y1 + y2)**2)/8.
      q(40,35)=(3*(x1 + x2)*(y1 + y2)**2*(z1 + z2))/16.
      q(40,36)=((x1 + x2)*(y1 + y2)**3)/4.
      q(40,37)=0
      q(40,38)=0
      q(40,39)=0
      q(40,40)=0
      q(40,41)=(x1 + x2)**2/4.
      q(40,42)=((x1 + x2)**2*(z1 + z2))/8.
      q(40,43)=((x1 + x2)**2*(z1 + z2)**2)/16.
      q(40,44)=((x1 + x2)**2*(y1 + y2))/4.
      q(40,45)=((x1 + x2)**2*(y1 + y2)*(z1 + z2))/8.
      q(40,46)=(3*(x1 + x2)**2*(y1 + y2)**2)/16.
      q(40,47)=0
      q(40,48)=0
      q(40,49)=0
      q(40,50)=(x1 + x2)**3/8.
      q(40,51)=((x1 + x2)**3*(z1 + z2))/16.
      q(40,52)=((x1 + x2)**3*(y1 + y2))/8.
      q(40,53)=0
      q(40,54)=0
      q(40,55)=(x1 + x2)**4/16.
      q(40,56)=0
      q(41,1)=0
      q(41,2)=1
      q(41,3)=z1 + z2
      q(41,4)=(3*(z1 + z2)**2)/4.
      q(41,5)=(z1 + z2)**3/2.
      q(41,6)=(5*(z1 + z2)**4)/16.
      q(41,7)=0
      q(41,8)=(y1 + y2)/2.
      q(41,9)=((y1 + y2)*(z1 + z2))/2.
      q(41,10)=(3*(y1 + y2)*(z1 + z2)**2)/8.
      q(41,11)=((y1 + y2)*(z1 + z2)**3)/4.
      q(41,12)=0
      q(41,13)=(y1 + y2)**2/4.
      q(41,14)=((y1 + y2)**2*(z1 + z2))/4.
      q(41,15)=(3*(y1 + y2)**2*(z1 + z2)**2)/16.
      q(41,16)=0
      q(41,17)=(y1 + y2)**3/8.
      q(41,18)=((y1 + y2)**3*(z1 + z2))/8.
      q(41,19)=0
      q(41,20)=(y1 + y2)**4/16.
      q(41,21)=0
      q(41,22)=0
      q(41,23)=(x1 + x2)/2.
      q(41,24)=((x1 + x2)*(z1 + z2))/2.
      q(41,25)=(3*(x1 + x2)*(z1 + z2)**2)/8.
      q(41,26)=((x1 + x2)*(z1 + z2)**3)/4.
      q(41,27)=0
      q(41,28)=((x1 + x2)*(y1 + y2))/4.
      q(41,29)=((x1 + x2)*(y1 + y2)*(z1 + z2))/4.
      q(41,30)=(3*(x1 + x2)*(y1 + y2)*(z1 + z2)**2)/16.
      q(41,31)=0
      q(41,32)=((x1 + x2)*(y1 + y2)**2)/8.
      q(41,33)=((x1 + x2)*(y1 + y2)**2*(z1 + z2))/8.
      q(41,34)=0
      q(41,35)=((x1 + x2)*(y1 + y2)**3)/16.
      q(41,36)=0
      q(41,37)=0
      q(41,38)=(x1 + x2)**2/4.
      q(41,39)=((x1 + x2)**2*(z1 + z2))/4.
      q(41,40)=(3*(x1 + x2)**2*(z1 + z2)**2)/16.
      q(41,41)=0
      q(41,42)=((x1 + x2)**2*(y1 + y2))/8.
      q(41,43)=((x1 + x2)**2*(y1 + y2)*(z1 + z2))/8.
      q(41,44)=0
      q(41,45)=((x1 + x2)**2*(y1 + y2)**2)/16.
      q(41,46)=0
      q(41,47)=0
      q(41,48)=(x1 + x2)**3/8.
      q(41,49)=((x1 + x2)**3*(z1 + z2))/8.
      q(41,50)=0
      q(41,51)=((x1 + x2)**3*(y1 + y2))/16.
      q(41,52)=0
      q(41,53)=0
      q(41,54)=(x1 + x2)**4/16.
      q(41,55)=0
      q(41,56)=0
      q(42,1)=0
      q(42,2)=0
      q(42,3)=0
      q(42,4)=0
      q(42,5)=0
      q(42,6)=0
      q(42,7)=0
      q(42,8)=0
      q(42,9)=0
      q(42,10)=0
      q(42,11)=0
      q(42,12)=0
      q(42,13)=0
      q(42,14)=0
      q(42,15)=0
      q(42,16)=0
      q(42,17)=0
      q(42,18)=0
      q(42,19)=0
      q(42,20)=0
      q(42,21)=0
      q(42,22)=1
      q(42,23)=(z1 + z3)/2.
      q(42,24)=(z1 + z3)**2/4.
      q(42,25)=(z1 + z3)**3/8.
      q(42,26)=(z1 + z3)**4/16.
      q(42,27)=(y1 + y3)/2.
      q(42,28)=((y1 + y3)*(z1 + z3))/4.
      q(42,29)=((y1 + y3)*(z1 + z3)**2)/8.
      q(42,30)=((y1 + y3)*(z1 + z3)**3)/16.
      q(42,31)=(y1 + y3)**2/4.
      q(42,32)=((y1 + y3)**2*(z1 + z3))/8.
      q(42,33)=((y1 + y3)**2*(z1 + z3)**2)/16.
      q(42,34)=(y1 + y3)**3/8.
      q(42,35)=((y1 + y3)**3*(z1 + z3))/16.
      q(42,36)=(y1 + y3)**4/16.
      q(42,37)=x1 + x3
      q(42,38)=((x1 + x3)*(z1 + z3))/2.
      q(42,39)=((x1 + x3)*(z1 + z3)**2)/4.
      q(42,40)=((x1 + x3)*(z1 + z3)**3)/8.
      q(42,41)=((x1 + x3)*(y1 + y3))/2.
      q(42,42)=((x1 + x3)*(y1 + y3)*(z1 + z3))/4.
      q(42,43)=((x1 + x3)*(y1 + y3)*(z1 + z3)**2)/8.
      q(42,44)=((x1 + x3)*(y1 + y3)**2)/4.
      q(42,45)=((x1 + x3)*(y1 + y3)**2*(z1 + z3))/8.
      q(42,46)=((x1 + x3)*(y1 + y3)**3)/8.
      q(42,47)=(3*(x1 + x3)**2)/4.
      q(42,48)=(3*(x1 + x3)**2*(z1 + z3))/8.
      q(42,49)=(3*(x1 + x3)**2*(z1 + z3)**2)/16.
      q(42,50)=(3*(x1 + x3)**2*(y1 + y3))/8.
      q(42,51)=(3*(x1 + x3)**2*(y1 + y3)*(z1 + z3))/16.
      q(42,52)=(3*(x1 + x3)**2*(y1 + y3)**2)/16.
      q(42,53)=(x1 + x3)**3/2.
      q(42,54)=((x1 + x3)**3*(z1 + z3))/4.
      q(42,55)=((x1 + x3)**3*(y1 + y3))/4.
      q(42,56)=(5*(x1 + x3)**4)/16.
      q(43,1)=0
      q(43,2)=0
      q(43,3)=0
      q(43,4)=0
      q(43,5)=0
      q(43,6)=0
      q(43,7)=1
      q(43,8)=(z1 + z3)/2.
      q(43,9)=(z1 + z3)**2/4.
      q(43,10)=(z1 + z3)**3/8.
      q(43,11)=(z1 + z3)**4/16.
      q(43,12)=y1 + y3
      q(43,13)=((y1 + y3)*(z1 + z3))/2.
      q(43,14)=((y1 + y3)*(z1 + z3)**2)/4.
      q(43,15)=((y1 + y3)*(z1 + z3)**3)/8.
      q(43,16)=(3*(y1 + y3)**2)/4.
      q(43,17)=(3*(y1 + y3)**2*(z1 + z3))/8.
      q(43,18)=(3*(y1 + y3)**2*(z1 + z3)**2)/16.
      q(43,19)=(y1 + y3)**3/2.
      q(43,20)=((y1 + y3)**3*(z1 + z3))/4.
      q(43,21)=(5*(y1 + y3)**4)/16.
      q(43,22)=0
      q(43,23)=0
      q(43,24)=0
      q(43,25)=0
      q(43,26)=0
      q(43,27)=(x1 + x3)/2.
      q(43,28)=((x1 + x3)*(z1 + z3))/4.
      q(43,29)=((x1 + x3)*(z1 + z3)**2)/8.
      q(43,30)=((x1 + x3)*(z1 + z3)**3)/16.
      q(43,31)=((x1 + x3)*(y1 + y3))/2.
      q(43,32)=((x1 + x3)*(y1 + y3)*(z1 + z3))/4.
      q(43,33)=((x1 + x3)*(y1 + y3)*(z1 + z3)**2)/8.
      q(43,34)=(3*(x1 + x3)*(y1 + y3)**2)/8.
      q(43,35)=(3*(x1 + x3)*(y1 + y3)**2*(z1 + z3))/16.
      q(43,36)=((x1 + x3)*(y1 + y3)**3)/4.
      q(43,37)=0
      q(43,38)=0
      q(43,39)=0
      q(43,40)=0
      q(43,41)=(x1 + x3)**2/4.
      q(43,42)=((x1 + x3)**2*(z1 + z3))/8.
      q(43,43)=((x1 + x3)**2*(z1 + z3)**2)/16.
      q(43,44)=((x1 + x3)**2*(y1 + y3))/4.
      q(43,45)=((x1 + x3)**2*(y1 + y3)*(z1 + z3))/8.
      q(43,46)=(3*(x1 + x3)**2*(y1 + y3)**2)/16.
      q(43,47)=0
      q(43,48)=0
      q(43,49)=0
      q(43,50)=(x1 + x3)**3/8.
      q(43,51)=((x1 + x3)**3*(z1 + z3))/16.
      q(43,52)=((x1 + x3)**3*(y1 + y3))/8.
      q(43,53)=0
      q(43,54)=0
      q(43,55)=(x1 + x3)**4/16.
      q(43,56)=0
      q(44,1)=0
      q(44,2)=1
      q(44,3)=z1 + z3
      q(44,4)=(3*(z1 + z3)**2)/4.
      q(44,5)=(z1 + z3)**3/2.
      q(44,6)=(5*(z1 + z3)**4)/16.
      q(44,7)=0
      q(44,8)=(y1 + y3)/2.
      q(44,9)=((y1 + y3)*(z1 + z3))/2.
      q(44,10)=(3*(y1 + y3)*(z1 + z3)**2)/8.
      q(44,11)=((y1 + y3)*(z1 + z3)**3)/4.
      q(44,12)=0
      q(44,13)=(y1 + y3)**2/4.
      q(44,14)=((y1 + y3)**2*(z1 + z3))/4.
      q(44,15)=(3*(y1 + y3)**2*(z1 + z3)**2)/16.
      q(44,16)=0
      q(44,17)=(y1 + y3)**3/8.
      q(44,18)=((y1 + y3)**3*(z1 + z3))/8.
      q(44,19)=0
      q(44,20)=(y1 + y3)**4/16.
      q(44,21)=0
      q(44,22)=0
      q(44,23)=(x1 + x3)/2.
      q(44,24)=((x1 + x3)*(z1 + z3))/2.
      q(44,25)=(3*(x1 + x3)*(z1 + z3)**2)/8.
      q(44,26)=((x1 + x3)*(z1 + z3)**3)/4.
      q(44,27)=0
      q(44,28)=((x1 + x3)*(y1 + y3))/4.
      q(44,29)=((x1 + x3)*(y1 + y3)*(z1 + z3))/4.
      q(44,30)=(3*(x1 + x3)*(y1 + y3)*(z1 + z3)**2)/16.
      q(44,31)=0
      q(44,32)=((x1 + x3)*(y1 + y3)**2)/8.
      q(44,33)=((x1 + x3)*(y1 + y3)**2*(z1 + z3))/8.
      q(44,34)=0
      q(44,35)=((x1 + x3)*(y1 + y3)**3)/16.
      q(44,36)=0
      q(44,37)=0
      q(44,38)=(x1 + x3)**2/4.
      q(44,39)=((x1 + x3)**2*(z1 + z3))/4.
      q(44,40)=(3*(x1 + x3)**2*(z1 + z3)**2)/16.
      q(44,41)=0
      q(44,42)=((x1 + x3)**2*(y1 + y3))/8.
      q(44,43)=((x1 + x3)**2*(y1 + y3)*(z1 + z3))/8.
      q(44,44)=0
      q(44,45)=((x1 + x3)**2*(y1 + y3)**2)/16.
      q(44,46)=0
      q(44,47)=0
      q(44,48)=(x1 + x3)**3/8.
      q(44,49)=((x1 + x3)**3*(z1 + z3))/8.
      q(44,50)=0
      q(44,51)=((x1 + x3)**3*(y1 + y3))/16.
      q(44,52)=0
      q(44,53)=0
      q(44,54)=(x1 + x3)**4/16.
      q(44,55)=0
      q(44,56)=0
      q(45,1)=0
      q(45,2)=0
      q(45,3)=0
      q(45,4)=0
      q(45,5)=0
      q(45,6)=0
      q(45,7)=0
      q(45,8)=0
      q(45,9)=0
      q(45,10)=0
      q(45,11)=0
      q(45,12)=0
      q(45,13)=0
      q(45,14)=0
      q(45,15)=0
      q(45,16)=0
      q(45,17)=0
      q(45,18)=0
      q(45,19)=0
      q(45,20)=0
      q(45,21)=0
      q(45,22)=1
      q(45,23)=(z1 + z4)/2.
      q(45,24)=(z1 + z4)**2/4.
      q(45,25)=(z1 + z4)**3/8.
      q(45,26)=(z1 + z4)**4/16.
      q(45,27)=(y1 + y4)/2.
      q(45,28)=((y1 + y4)*(z1 + z4))/4.
      q(45,29)=((y1 + y4)*(z1 + z4)**2)/8.
      q(45,30)=((y1 + y4)*(z1 + z4)**3)/16.
      q(45,31)=(y1 + y4)**2/4.
      q(45,32)=((y1 + y4)**2*(z1 + z4))/8.
      q(45,33)=((y1 + y4)**2*(z1 + z4)**2)/16.
      q(45,34)=(y1 + y4)**3/8.
      q(45,35)=((y1 + y4)**3*(z1 + z4))/16.
      q(45,36)=(y1 + y4)**4/16.
      q(45,37)=x1 + x4
      q(45,38)=((x1 + x4)*(z1 + z4))/2.
      q(45,39)=((x1 + x4)*(z1 + z4)**2)/4.
      q(45,40)=((x1 + x4)*(z1 + z4)**3)/8.
      q(45,41)=((x1 + x4)*(y1 + y4))/2.
      q(45,42)=((x1 + x4)*(y1 + y4)*(z1 + z4))/4.
      q(45,43)=((x1 + x4)*(y1 + y4)*(z1 + z4)**2)/8.
      q(45,44)=((x1 + x4)*(y1 + y4)**2)/4.
      q(45,45)=((x1 + x4)*(y1 + y4)**2*(z1 + z4))/8.
      q(45,46)=((x1 + x4)*(y1 + y4)**3)/8.
      q(45,47)=(3*(x1 + x4)**2)/4.
      q(45,48)=(3*(x1 + x4)**2*(z1 + z4))/8.
      q(45,49)=(3*(x1 + x4)**2*(z1 + z4)**2)/16.
      q(45,50)=(3*(x1 + x4)**2*(y1 + y4))/8.
      q(45,51)=(3*(x1 + x4)**2*(y1 + y4)*(z1 + z4))/16.
      q(45,52)=(3*(x1 + x4)**2*(y1 + y4)**2)/16.
      q(45,53)=(x1 + x4)**3/2.
      q(45,54)=((x1 + x4)**3*(z1 + z4))/4.
      q(45,55)=((x1 + x4)**3*(y1 + y4))/4.
      q(45,56)=(5*(x1 + x4)**4)/16.
      q(46,1)=0
      q(46,2)=0
      q(46,3)=0
      q(46,4)=0
      q(46,5)=0
      q(46,6)=0
      q(46,7)=1
      q(46,8)=(z1 + z4)/2.
      q(46,9)=(z1 + z4)**2/4.
      q(46,10)=(z1 + z4)**3/8.
      q(46,11)=(z1 + z4)**4/16.
      q(46,12)=y1 + y4
      q(46,13)=((y1 + y4)*(z1 + z4))/2.
      q(46,14)=((y1 + y4)*(z1 + z4)**2)/4.
      q(46,15)=((y1 + y4)*(z1 + z4)**3)/8.
      q(46,16)=(3*(y1 + y4)**2)/4.
      q(46,17)=(3*(y1 + y4)**2*(z1 + z4))/8.
      q(46,18)=(3*(y1 + y4)**2*(z1 + z4)**2)/16.
      q(46,19)=(y1 + y4)**3/2.
      q(46,20)=((y1 + y4)**3*(z1 + z4))/4.
      q(46,21)=(5*(y1 + y4)**4)/16.
      q(46,22)=0
      q(46,23)=0
      q(46,24)=0
      q(46,25)=0
      q(46,26)=0
      q(46,27)=(x1 + x4)/2.
      q(46,28)=((x1 + x4)*(z1 + z4))/4.
      q(46,29)=((x1 + x4)*(z1 + z4)**2)/8.
      q(46,30)=((x1 + x4)*(z1 + z4)**3)/16.
      q(46,31)=((x1 + x4)*(y1 + y4))/2.
      q(46,32)=((x1 + x4)*(y1 + y4)*(z1 + z4))/4.
      q(46,33)=((x1 + x4)*(y1 + y4)*(z1 + z4)**2)/8.
      q(46,34)=(3*(x1 + x4)*(y1 + y4)**2)/8.
      q(46,35)=(3*(x1 + x4)*(y1 + y4)**2*(z1 + z4))/16.
      q(46,36)=((x1 + x4)*(y1 + y4)**3)/4.
      q(46,37)=0
      q(46,38)=0
      q(46,39)=0
      q(46,40)=0
      q(46,41)=(x1 + x4)**2/4.
      q(46,42)=((x1 + x4)**2*(z1 + z4))/8.
      q(46,43)=((x1 + x4)**2*(z1 + z4)**2)/16.
      q(46,44)=((x1 + x4)**2*(y1 + y4))/4.
      q(46,45)=((x1 + x4)**2*(y1 + y4)*(z1 + z4))/8.
      q(46,46)=(3*(x1 + x4)**2*(y1 + y4)**2)/16.
      q(46,47)=0
      q(46,48)=0
      q(46,49)=0
      q(46,50)=(x1 + x4)**3/8.
      q(46,51)=((x1 + x4)**3*(z1 + z4))/16.
      q(46,52)=((x1 + x4)**3*(y1 + y4))/8.
      q(46,53)=0
      q(46,54)=0
      q(46,55)=(x1 + x4)**4/16.
      q(46,56)=0
      q(47,1)=0
      q(47,2)=1
      q(47,3)=z1 + z4
      q(47,4)=(3*(z1 + z4)**2)/4.
      q(47,5)=(z1 + z4)**3/2.
      q(47,6)=(5*(z1 + z4)**4)/16.
      q(47,7)=0
      q(47,8)=(y1 + y4)/2.
      q(47,9)=((y1 + y4)*(z1 + z4))/2.
      q(47,10)=(3*(y1 + y4)*(z1 + z4)**2)/8.
      q(47,11)=((y1 + y4)*(z1 + z4)**3)/4.
      q(47,12)=0
      q(47,13)=(y1 + y4)**2/4.
      q(47,14)=((y1 + y4)**2*(z1 + z4))/4.
      q(47,15)=(3*(y1 + y4)**2*(z1 + z4)**2)/16.
      q(47,16)=0
      q(47,17)=(y1 + y4)**3/8.
      q(47,18)=((y1 + y4)**3*(z1 + z4))/8.
      q(47,19)=0
      q(47,20)=(y1 + y4)**4/16.
      q(47,21)=0
      q(47,22)=0
      q(47,23)=(x1 + x4)/2.
      q(47,24)=((x1 + x4)*(z1 + z4))/2.
      q(47,25)=(3*(x1 + x4)*(z1 + z4)**2)/8.
      q(47,26)=((x1 + x4)*(z1 + z4)**3)/4.
      q(47,27)=0
      q(47,28)=((x1 + x4)*(y1 + y4))/4.
      q(47,29)=((x1 + x4)*(y1 + y4)*(z1 + z4))/4.
      q(47,30)=(3*(x1 + x4)*(y1 + y4)*(z1 + z4)**2)/16.
      q(47,31)=0
      q(47,32)=((x1 + x4)*(y1 + y4)**2)/8.
      q(47,33)=((x1 + x4)*(y1 + y4)**2*(z1 + z4))/8.
      q(47,34)=0
      q(47,35)=((x1 + x4)*(y1 + y4)**3)/16.
      q(47,36)=0
      q(47,37)=0
      q(47,38)=(x1 + x4)**2/4.
      q(47,39)=((x1 + x4)**2*(z1 + z4))/4.
      q(47,40)=(3*(x1 + x4)**2*(z1 + z4)**2)/16.
      q(47,41)=0
      q(47,42)=((x1 + x4)**2*(y1 + y4))/8.
      q(47,43)=((x1 + x4)**2*(y1 + y4)*(z1 + z4))/8.
      q(47,44)=0
      q(47,45)=((x1 + x4)**2*(y1 + y4)**2)/16.
      q(47,46)=0
      q(47,47)=0
      q(47,48)=(x1 + x4)**3/8.
      q(47,49)=((x1 + x4)**3*(z1 + z4))/8.
      q(47,50)=0
      q(47,51)=((x1 + x4)**3*(y1 + y4))/16.
      q(47,52)=0
      q(47,53)=0
      q(47,54)=(x1 + x4)**4/16.
      q(47,55)=0
      q(47,56)=0
      q(48,1)=0
      q(48,2)=0
      q(48,3)=0
      q(48,4)=0
      q(48,5)=0
      q(48,6)=0
      q(48,7)=0
      q(48,8)=0
      q(48,9)=0
      q(48,10)=0
      q(48,11)=0
      q(48,12)=0
      q(48,13)=0
      q(48,14)=0
      q(48,15)=0
      q(48,16)=0
      q(48,17)=0
      q(48,18)=0
      q(48,19)=0
      q(48,20)=0
      q(48,21)=0
      q(48,22)=1
      q(48,23)=(z2 + z3)/2.
      q(48,24)=(z2 + z3)**2/4.
      q(48,25)=(z2 + z3)**3/8.
      q(48,26)=(z2 + z3)**4/16.
      q(48,27)=(y2 + y3)/2.
      q(48,28)=((y2 + y3)*(z2 + z3))/4.
      q(48,29)=((y2 + y3)*(z2 + z3)**2)/8.
      q(48,30)=((y2 + y3)*(z2 + z3)**3)/16.
      q(48,31)=(y2 + y3)**2/4.
      q(48,32)=((y2 + y3)**2*(z2 + z3))/8.
      q(48,33)=((y2 + y3)**2*(z2 + z3)**2)/16.
      q(48,34)=(y2 + y3)**3/8.
      q(48,35)=((y2 + y3)**3*(z2 + z3))/16.
      q(48,36)=(y2 + y3)**4/16.
      q(48,37)=x2 + x3
      q(48,38)=((x2 + x3)*(z2 + z3))/2.
      q(48,39)=((x2 + x3)*(z2 + z3)**2)/4.
      q(48,40)=((x2 + x3)*(z2 + z3)**3)/8.
      q(48,41)=((x2 + x3)*(y2 + y3))/2.
      q(48,42)=((x2 + x3)*(y2 + y3)*(z2 + z3))/4.
      q(48,43)=((x2 + x3)*(y2 + y3)*(z2 + z3)**2)/8.
      q(48,44)=((x2 + x3)*(y2 + y3)**2)/4.
      q(48,45)=((x2 + x3)*(y2 + y3)**2*(z2 + z3))/8.
      q(48,46)=((x2 + x3)*(y2 + y3)**3)/8.
      q(48,47)=(3*(x2 + x3)**2)/4.
      q(48,48)=(3*(x2 + x3)**2*(z2 + z3))/8.
      q(48,49)=(3*(x2 + x3)**2*(z2 + z3)**2)/16.
      q(48,50)=(3*(x2 + x3)**2*(y2 + y3))/8.
      q(48,51)=(3*(x2 + x3)**2*(y2 + y3)*(z2 + z3))/16.
      q(48,52)=(3*(x2 + x3)**2*(y2 + y3)**2)/16.
      q(48,53)=(x2 + x3)**3/2.
      q(48,54)=((x2 + x3)**3*(z2 + z3))/4.
      q(48,55)=((x2 + x3)**3*(y2 + y3))/4.
      q(48,56)=(5*(x2 + x3)**4)/16.
      q(49,1)=0
      q(49,2)=0
      q(49,3)=0
      q(49,4)=0
      q(49,5)=0
      q(49,6)=0
      q(49,7)=1
      q(49,8)=(z2 + z3)/2.
      q(49,9)=(z2 + z3)**2/4.
      q(49,10)=(z2 + z3)**3/8.
      q(49,11)=(z2 + z3)**4/16.
      q(49,12)=y2 + y3
      q(49,13)=((y2 + y3)*(z2 + z3))/2.
      q(49,14)=((y2 + y3)*(z2 + z3)**2)/4.
      q(49,15)=((y2 + y3)*(z2 + z3)**3)/8.
      q(49,16)=(3*(y2 + y3)**2)/4.
      q(49,17)=(3*(y2 + y3)**2*(z2 + z3))/8.
      q(49,18)=(3*(y2 + y3)**2*(z2 + z3)**2)/16.
      q(49,19)=(y2 + y3)**3/2.
      q(49,20)=((y2 + y3)**3*(z2 + z3))/4.
      q(49,21)=(5*(y2 + y3)**4)/16.
      q(49,22)=0
      q(49,23)=0
      q(49,24)=0
      q(49,25)=0
      q(49,26)=0
      q(49,27)=(x2 + x3)/2.
      q(49,28)=((x2 + x3)*(z2 + z3))/4.
      q(49,29)=((x2 + x3)*(z2 + z3)**2)/8.
      q(49,30)=((x2 + x3)*(z2 + z3)**3)/16.
      q(49,31)=((x2 + x3)*(y2 + y3))/2.
      q(49,32)=((x2 + x3)*(y2 + y3)*(z2 + z3))/4.
      q(49,33)=((x2 + x3)*(y2 + y3)*(z2 + z3)**2)/8.
      q(49,34)=(3*(x2 + x3)*(y2 + y3)**2)/8.
      q(49,35)=(3*(x2 + x3)*(y2 + y3)**2*(z2 + z3))/16.
      q(49,36)=((x2 + x3)*(y2 + y3)**3)/4.
      q(49,37)=0
      q(49,38)=0
      q(49,39)=0
      q(49,40)=0
      q(49,41)=(x2 + x3)**2/4.
      q(49,42)=((x2 + x3)**2*(z2 + z3))/8.
      q(49,43)=((x2 + x3)**2*(z2 + z3)**2)/16.
      q(49,44)=((x2 + x3)**2*(y2 + y3))/4.
      q(49,45)=((x2 + x3)**2*(y2 + y3)*(z2 + z3))/8.
      q(49,46)=(3*(x2 + x3)**2*(y2 + y3)**2)/16.
      q(49,47)=0
      q(49,48)=0
      q(49,49)=0
      q(49,50)=(x2 + x3)**3/8.
      q(49,51)=((x2 + x3)**3*(z2 + z3))/16.
      q(49,52)=((x2 + x3)**3*(y2 + y3))/8.
      q(49,53)=0
      q(49,54)=0
      q(49,55)=(x2 + x3)**4/16.
      q(49,56)=0
      q(50,1)=0
      q(50,2)=1
      q(50,3)=z2 + z3
      q(50,4)=(3*(z2 + z3)**2)/4.
      q(50,5)=(z2 + z3)**3/2.
      q(50,6)=(5*(z2 + z3)**4)/16.
      q(50,7)=0
      q(50,8)=(y2 + y3)/2.
      q(50,9)=((y2 + y3)*(z2 + z3))/2.
      q(50,10)=(3*(y2 + y3)*(z2 + z3)**2)/8.
      q(50,11)=((y2 + y3)*(z2 + z3)**3)/4.
      q(50,12)=0
      q(50,13)=(y2 + y3)**2/4.
      q(50,14)=((y2 + y3)**2*(z2 + z3))/4.
      q(50,15)=(3*(y2 + y3)**2*(z2 + z3)**2)/16.
      q(50,16)=0
      q(50,17)=(y2 + y3)**3/8.
      q(50,18)=((y2 + y3)**3*(z2 + z3))/8.
      q(50,19)=0
      q(50,20)=(y2 + y3)**4/16.
      q(50,21)=0
      q(50,22)=0
      q(50,23)=(x2 + x3)/2.
      q(50,24)=((x2 + x3)*(z2 + z3))/2.
      q(50,25)=(3*(x2 + x3)*(z2 + z3)**2)/8.
      q(50,26)=((x2 + x3)*(z2 + z3)**3)/4.
      q(50,27)=0
      q(50,28)=((x2 + x3)*(y2 + y3))/4.
      q(50,29)=((x2 + x3)*(y2 + y3)*(z2 + z3))/4.
      q(50,30)=(3*(x2 + x3)*(y2 + y3)*(z2 + z3)**2)/16.
      q(50,31)=0
      q(50,32)=((x2 + x3)*(y2 + y3)**2)/8.
      q(50,33)=((x2 + x3)*(y2 + y3)**2*(z2 + z3))/8.
      q(50,34)=0
      q(50,35)=((x2 + x3)*(y2 + y3)**3)/16.
      q(50,36)=0
      q(50,37)=0
      q(50,38)=(x2 + x3)**2/4.
      q(50,39)=((x2 + x3)**2*(z2 + z3))/4.
      q(50,40)=(3*(x2 + x3)**2*(z2 + z3)**2)/16.
      q(50,41)=0
      q(50,42)=((x2 + x3)**2*(y2 + y3))/8.
      q(50,43)=((x2 + x3)**2*(y2 + y3)*(z2 + z3))/8.
      q(50,44)=0
      q(50,45)=((x2 + x3)**2*(y2 + y3)**2)/16.
      q(50,46)=0
      q(50,47)=0
      q(50,48)=(x2 + x3)**3/8.
      q(50,49)=((x2 + x3)**3*(z2 + z3))/8.
      q(50,50)=0
      q(50,51)=((x2 + x3)**3*(y2 + y3))/16.
      q(50,52)=0
      q(50,53)=0
      q(50,54)=(x2 + x3)**4/16.
      q(50,55)=0
      q(50,56)=0
      q(51,1)=0
      q(51,2)=0
      q(51,3)=0
      q(51,4)=0
      q(51,5)=0
      q(51,6)=0
      q(51,7)=0
      q(51,8)=0
      q(51,9)=0
      q(51,10)=0
      q(51,11)=0
      q(51,12)=0
      q(51,13)=0
      q(51,14)=0
      q(51,15)=0
      q(51,16)=0
      q(51,17)=0
      q(51,18)=0
      q(51,19)=0
      q(51,20)=0
      q(51,21)=0
      q(51,22)=1
      q(51,23)=(z2 + z4)/2.
      q(51,24)=(z2 + z4)**2/4.
      q(51,25)=(z2 + z4)**3/8.
      q(51,26)=(z2 + z4)**4/16.
      q(51,27)=(y2 + y4)/2.
      q(51,28)=((y2 + y4)*(z2 + z4))/4.
      q(51,29)=((y2 + y4)*(z2 + z4)**2)/8.
      q(51,30)=((y2 + y4)*(z2 + z4)**3)/16.
      q(51,31)=(y2 + y4)**2/4.
      q(51,32)=((y2 + y4)**2*(z2 + z4))/8.
      q(51,33)=((y2 + y4)**2*(z2 + z4)**2)/16.
      q(51,34)=(y2 + y4)**3/8.
      q(51,35)=((y2 + y4)**3*(z2 + z4))/16.
      q(51,36)=(y2 + y4)**4/16.
      q(51,37)=x2 + x4
      q(51,38)=((x2 + x4)*(z2 + z4))/2.
      q(51,39)=((x2 + x4)*(z2 + z4)**2)/4.
      q(51,40)=((x2 + x4)*(z2 + z4)**3)/8.
      q(51,41)=((x2 + x4)*(y2 + y4))/2.
      q(51,42)=((x2 + x4)*(y2 + y4)*(z2 + z4))/4.
      q(51,43)=((x2 + x4)*(y2 + y4)*(z2 + z4)**2)/8.
      q(51,44)=((x2 + x4)*(y2 + y4)**2)/4.
      q(51,45)=((x2 + x4)*(y2 + y4)**2*(z2 + z4))/8.
      q(51,46)=((x2 + x4)*(y2 + y4)**3)/8.
      q(51,47)=(3*(x2 + x4)**2)/4.
      q(51,48)=(3*(x2 + x4)**2*(z2 + z4))/8.
      q(51,49)=(3*(x2 + x4)**2*(z2 + z4)**2)/16.
      q(51,50)=(3*(x2 + x4)**2*(y2 + y4))/8.
      q(51,51)=(3*(x2 + x4)**2*(y2 + y4)*(z2 + z4))/16.
      q(51,52)=(3*(x2 + x4)**2*(y2 + y4)**2)/16.
      q(51,53)=(x2 + x4)**3/2.
      q(51,54)=((x2 + x4)**3*(z2 + z4))/4.
      q(51,55)=((x2 + x4)**3*(y2 + y4))/4.
      q(51,56)=(5*(x2 + x4)**4)/16.
      q(52,1)=0
      q(52,2)=0
      q(52,3)=0
      q(52,4)=0
      q(52,5)=0
      q(52,6)=0
      q(52,7)=1
      q(52,8)=(z2 + z4)/2.
      q(52,9)=(z2 + z4)**2/4.
      q(52,10)=(z2 + z4)**3/8.
      q(52,11)=(z2 + z4)**4/16.
      q(52,12)=y2 + y4
      q(52,13)=((y2 + y4)*(z2 + z4))/2.
      q(52,14)=((y2 + y4)*(z2 + z4)**2)/4.
      q(52,15)=((y2 + y4)*(z2 + z4)**3)/8.
      q(52,16)=(3*(y2 + y4)**2)/4.
      q(52,17)=(3*(y2 + y4)**2*(z2 + z4))/8.
      q(52,18)=(3*(y2 + y4)**2*(z2 + z4)**2)/16.
      q(52,19)=(y2 + y4)**3/2.
      q(52,20)=((y2 + y4)**3*(z2 + z4))/4.
      q(52,21)=(5*(y2 + y4)**4)/16.
      q(52,22)=0
      q(52,23)=0
      q(52,24)=0
      q(52,25)=0
      q(52,26)=0
      q(52,27)=(x2 + x4)/2.
      q(52,28)=((x2 + x4)*(z2 + z4))/4.
      q(52,29)=((x2 + x4)*(z2 + z4)**2)/8.
      q(52,30)=((x2 + x4)*(z2 + z4)**3)/16.
      q(52,31)=((x2 + x4)*(y2 + y4))/2.
      q(52,32)=((x2 + x4)*(y2 + y4)*(z2 + z4))/4.
      q(52,33)=((x2 + x4)*(y2 + y4)*(z2 + z4)**2)/8.
      q(52,34)=(3*(x2 + x4)*(y2 + y4)**2)/8.
      q(52,35)=(3*(x2 + x4)*(y2 + y4)**2*(z2 + z4))/16.
      q(52,36)=((x2 + x4)*(y2 + y4)**3)/4.
      q(52,37)=0
      q(52,38)=0
      q(52,39)=0
      q(52,40)=0
      q(52,41)=(x2 + x4)**2/4.
      q(52,42)=((x2 + x4)**2*(z2 + z4))/8.
      q(52,43)=((x2 + x4)**2*(z2 + z4)**2)/16.
      q(52,44)=((x2 + x4)**2*(y2 + y4))/4.
      q(52,45)=((x2 + x4)**2*(y2 + y4)*(z2 + z4))/8.
      q(52,46)=(3*(x2 + x4)**2*(y2 + y4)**2)/16.
      q(52,47)=0
      q(52,48)=0
      q(52,49)=0
      q(52,50)=(x2 + x4)**3/8.
      q(52,51)=((x2 + x4)**3*(z2 + z4))/16.
      q(52,52)=((x2 + x4)**3*(y2 + y4))/8.
      q(52,53)=0
      q(52,54)=0
      q(52,55)=(x2 + x4)**4/16.
      q(52,56)=0
      q(53,1)=0
      q(53,2)=1
      q(53,3)=z2 + z4
      q(53,4)=(3*(z2 + z4)**2)/4.
      q(53,5)=(z2 + z4)**3/2.
      q(53,6)=(5*(z2 + z4)**4)/16.
      q(53,7)=0
      q(53,8)=(y2 + y4)/2.
      q(53,9)=((y2 + y4)*(z2 + z4))/2.
      q(53,10)=(3*(y2 + y4)*(z2 + z4)**2)/8.
      q(53,11)=((y2 + y4)*(z2 + z4)**3)/4.
      q(53,12)=0
      q(53,13)=(y2 + y4)**2/4.
      q(53,14)=((y2 + y4)**2*(z2 + z4))/4.
      q(53,15)=(3*(y2 + y4)**2*(z2 + z4)**2)/16.
      q(53,16)=0
      q(53,17)=(y2 + y4)**3/8.
      q(53,18)=((y2 + y4)**3*(z2 + z4))/8.
      q(53,19)=0
      q(53,20)=(y2 + y4)**4/16.
      q(53,21)=0
      q(53,22)=0
      q(53,23)=(x2 + x4)/2.
      q(53,24)=((x2 + x4)*(z2 + z4))/2.
      q(53,25)=(3*(x2 + x4)*(z2 + z4)**2)/8.
      q(53,26)=((x2 + x4)*(z2 + z4)**3)/4.
      q(53,27)=0
      q(53,28)=((x2 + x4)*(y2 + y4))/4.
      q(53,29)=((x2 + x4)*(y2 + y4)*(z2 + z4))/4.
      q(53,30)=(3*(x2 + x4)*(y2 + y4)*(z2 + z4)**2)/16.
      q(53,31)=0
      q(53,32)=((x2 + x4)*(y2 + y4)**2)/8.
      q(53,33)=((x2 + x4)*(y2 + y4)**2*(z2 + z4))/8.
      q(53,34)=0
      q(53,35)=((x2 + x4)*(y2 + y4)**3)/16.
      q(53,36)=0
      q(53,37)=0
      q(53,38)=(x2 + x4)**2/4.
      q(53,39)=((x2 + x4)**2*(z2 + z4))/4.
      q(53,40)=(3*(x2 + x4)**2*(z2 + z4)**2)/16.
      q(53,41)=0
      q(53,42)=((x2 + x4)**2*(y2 + y4))/8.
      q(53,43)=((x2 + x4)**2*(y2 + y4)*(z2 + z4))/8.
      q(53,44)=0
      q(53,45)=((x2 + x4)**2*(y2 + y4)**2)/16.
      q(53,46)=0
      q(53,47)=0
      q(53,48)=(x2 + x4)**3/8.
      q(53,49)=((x2 + x4)**3*(z2 + z4))/8.
      q(53,50)=0
      q(53,51)=((x2 + x4)**3*(y2 + y4))/16.
      q(53,52)=0
      q(53,53)=0
      q(53,54)=(x2 + x4)**4/16.
      q(53,55)=0
      q(53,56)=0
      q(54,1)=0
      q(54,2)=0
      q(54,3)=0
      q(54,4)=0
      q(54,5)=0
      q(54,6)=0
      q(54,7)=0
      q(54,8)=0
      q(54,9)=0
      q(54,10)=0
      q(54,11)=0
      q(54,12)=0
      q(54,13)=0
      q(54,14)=0
      q(54,15)=0
      q(54,16)=0
      q(54,17)=0
      q(54,18)=0
      q(54,19)=0
      q(54,20)=0
      q(54,21)=0
      q(54,22)=1
      q(54,23)=(z3 + z4)/2.
      q(54,24)=(z3 + z4)**2/4.
      q(54,25)=(z3 + z4)**3/8.
      q(54,26)=(z3 + z4)**4/16.
      q(54,27)=(y3 + y4)/2.
      q(54,28)=((y3 + y4)*(z3 + z4))/4.
      q(54,29)=((y3 + y4)*(z3 + z4)**2)/8.
      q(54,30)=((y3 + y4)*(z3 + z4)**3)/16.
      q(54,31)=(y3 + y4)**2/4.
      q(54,32)=((y3 + y4)**2*(z3 + z4))/8.
      q(54,33)=((y3 + y4)**2*(z3 + z4)**2)/16.
      q(54,34)=(y3 + y4)**3/8.
      q(54,35)=((y3 + y4)**3*(z3 + z4))/16.
      q(54,36)=(y3 + y4)**4/16.
      q(54,37)=x3 + x4
      q(54,38)=((x3 + x4)*(z3 + z4))/2.
      q(54,39)=((x3 + x4)*(z3 + z4)**2)/4.
      q(54,40)=((x3 + x4)*(z3 + z4)**3)/8.
      q(54,41)=((x3 + x4)*(y3 + y4))/2.
      q(54,42)=((x3 + x4)*(y3 + y4)*(z3 + z4))/4.
      q(54,43)=((x3 + x4)*(y3 + y4)*(z3 + z4)**2)/8.
      q(54,44)=((x3 + x4)*(y3 + y4)**2)/4.
      q(54,45)=((x3 + x4)*(y3 + y4)**2*(z3 + z4))/8.
      q(54,46)=((x3 + x4)*(y3 + y4)**3)/8.
      q(54,47)=(3*(x3 + x4)**2)/4.
      q(54,48)=(3*(x3 + x4)**2*(z3 + z4))/8.
      q(54,49)=(3*(x3 + x4)**2*(z3 + z4)**2)/16.
      q(54,50)=(3*(x3 + x4)**2*(y3 + y4))/8.
      q(54,51)=(3*(x3 + x4)**2*(y3 + y4)*(z3 + z4))/16.
      q(54,52)=(3*(x3 + x4)**2*(y3 + y4)**2)/16.
      q(54,53)=(x3 + x4)**3/2.
      q(54,54)=((x3 + x4)**3*(z3 + z4))/4.
      q(54,55)=((x3 + x4)**3*(y3 + y4))/4.
      q(54,56)=(5*(x3 + x4)**4)/16.
      q(55,1)=0
      q(55,2)=0
      q(55,3)=0
      q(55,4)=0
      q(55,5)=0
      q(55,6)=0
      q(55,7)=1
      q(55,8)=(z3 + z4)/2.
      q(55,9)=(z3 + z4)**2/4.
      q(55,10)=(z3 + z4)**3/8.
      q(55,11)=(z3 + z4)**4/16.
      q(55,12)=y3 + y4
      q(55,13)=((y3 + y4)*(z3 + z4))/2.
      q(55,14)=((y3 + y4)*(z3 + z4)**2)/4.
      q(55,15)=((y3 + y4)*(z3 + z4)**3)/8.
      q(55,16)=(3*(y3 + y4)**2)/4.
      q(55,17)=(3*(y3 + y4)**2*(z3 + z4))/8.
      q(55,18)=(3*(y3 + y4)**2*(z3 + z4)**2)/16.
      q(55,19)=(y3 + y4)**3/2.
      q(55,20)=((y3 + y4)**3*(z3 + z4))/4.
      q(55,21)=(5*(y3 + y4)**4)/16.
      q(55,22)=0
      q(55,23)=0
      q(55,24)=0
      q(55,25)=0
      q(55,26)=0
      q(55,27)=(x3 + x4)/2.
      q(55,28)=((x3 + x4)*(z3 + z4))/4.
      q(55,29)=((x3 + x4)*(z3 + z4)**2)/8.
      q(55,30)=((x3 + x4)*(z3 + z4)**3)/16.
      q(55,31)=((x3 + x4)*(y3 + y4))/2.
      q(55,32)=((x3 + x4)*(y3 + y4)*(z3 + z4))/4.
      q(55,33)=((x3 + x4)*(y3 + y4)*(z3 + z4)**2)/8.
      q(55,34)=(3*(x3 + x4)*(y3 + y4)**2)/8.
      q(55,35)=(3*(x3 + x4)*(y3 + y4)**2*(z3 + z4))/16.
      q(55,36)=((x3 + x4)*(y3 + y4)**3)/4.
      q(55,37)=0
      q(55,38)=0
      q(55,39)=0
      q(55,40)=0
      q(55,41)=(x3 + x4)**2/4.
      q(55,42)=((x3 + x4)**2*(z3 + z4))/8.
      q(55,43)=((x3 + x4)**2*(z3 + z4)**2)/16.
      q(55,44)=((x3 + x4)**2*(y3 + y4))/4.
      q(55,45)=((x3 + x4)**2*(y3 + y4)*(z3 + z4))/8.
      q(55,46)=(3*(x3 + x4)**2*(y3 + y4)**2)/16.
      q(55,47)=0
      q(55,48)=0
      q(55,49)=0
      q(55,50)=(x3 + x4)**3/8.
      q(55,51)=((x3 + x4)**3*(z3 + z4))/16.
      q(55,52)=((x3 + x4)**3*(y3 + y4))/8.
      q(55,53)=0
      q(55,54)=0
      q(55,55)=(x3 + x4)**4/16.
      q(55,56)=0
      q(56,1)=0
      q(56,2)=1
      q(56,3)=z3 + z4
      q(56,4)=(3*(z3 + z4)**2)/4.
      q(56,5)=(z3 + z4)**3/2.
      q(56,6)=(5*(z3 + z4)**4)/16.
      q(56,7)=0
      q(56,8)=(y3 + y4)/2.
      q(56,9)=((y3 + y4)*(z3 + z4))/2.
      q(56,10)=(3*(y3 + y4)*(z3 + z4)**2)/8.
      q(56,11)=((y3 + y4)*(z3 + z4)**3)/4.
      q(56,12)=0
      q(56,13)=(y3 + y4)**2/4.
      q(56,14)=((y3 + y4)**2*(z3 + z4))/4.
      q(56,15)=(3*(y3 + y4)**2*(z3 + z4)**2)/16.
      q(56,16)=0
      q(56,17)=(y3 + y4)**3/8.
      q(56,18)=((y3 + y4)**3*(z3 + z4))/8.
      q(56,19)=0
      q(56,20)=(y3 + y4)**4/16.
      q(56,21)=0
      q(56,22)=0
      q(56,23)=(x3 + x4)/2.
      q(56,24)=((x3 + x4)*(z3 + z4))/2.
      q(56,25)=(3*(x3 + x4)*(z3 + z4)**2)/8.
      q(56,26)=((x3 + x4)*(z3 + z4)**3)/4.
      q(56,27)=0
      q(56,28)=((x3 + x4)*(y3 + y4))/4.
      q(56,29)=((x3 + x4)*(y3 + y4)*(z3 + z4))/4.
      q(56,30)=(3*(x3 + x4)*(y3 + y4)*(z3 + z4)**2)/16.
      q(56,31)=0
      q(56,32)=((x3 + x4)*(y3 + y4)**2)/8.
      q(56,33)=((x3 + x4)*(y3 + y4)**2*(z3 + z4))/8.
      q(56,34)=0
      q(56,35)=((x3 + x4)*(y3 + y4)**3)/16.
      q(56,36)=0
      q(56,37)=0
      q(56,38)=(x3 + x4)**2/4.
      q(56,39)=((x3 + x4)**2*(z3 + z4))/4.
      q(56,40)=(3*(x3 + x4)**2*(z3 + z4)**2)/16.
      q(56,41)=0
      q(56,42)=((x3 + x4)**2*(y3 + y4))/8.
      q(56,43)=((x3 + x4)**2*(y3 + y4)*(z3 + z4))/8.
      q(56,44)=0
      q(56,45)=((x3 + x4)**2*(y3 + y4)**2)/16.
      q(56,46)=0
      q(56,47)=0
      q(56,48)=(x3 + x4)**3/8.
      q(56,49)=((x3 + x4)**3*(z3 + z4))/8.
      q(56,50)=0
      q(56,51)=((x3 + x4)**3*(y3 + y4))/16.
      q(56,52)=0
      q(56,53)=0
      q(56,54)=(x3 + x4)**4/16.
      q(56,55)=0
      q(56,56)=0


!     we nned to invert the 56 x 56 matrix to map from the 56 p-DOF
!     to the 56 polynomial coefficients
      LDA=56
      NPIVOT=56
      call PIVOT(q,LDA,NPIVOT,matrixi,PS1,PS2,IERR,MX,MY,VAL)
      
      if(ierr.ne.1) then
        write(*,*)"FATAL: matrix inversion failed "
        read(*,*)
      end if

!     now unite the 56 x 56 matrix with the 56 x 16 matrix
!     the results maps from the 16 DOF to the 56 polynomial coefficients      
      do i=1,56
      do j=1,16
        matrix3(i,j)=0.d0
        do k=1,56
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
!    Auxiliary subroutine that returns 56, 56x3 and 56x3x3 matrices 
!    depending on the location such that the product with the DOF
!    gives the field value and the first and second gradient
!
      subroutine Bmatrices(x,y,z,Mi,B0,B1,B2)
      implicit none
      double precision x,y,z
      double precision Mi(56,16)
      double precision B0(4,4),B1(4,4,3),B2(4,4,3,3)

!                      value,gradient,gradientgradient      
      double precision n(56),gn(56,3),ggn(56,3,3)
      
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
      n(5)=z**4
      n(6)=z**5
      n(7)=y
      n(8)=y*z
      n(9)=y*z**2
      n(10)=y*z**3
      n(11)=y*z**4
      n(12)=y**2
      n(13)=y**2*z
      n(14)=y**2*z**2
      n(15)=y**2*z**3
      n(16)=y**3
      n(17)=y**3*z
      n(18)=y**3*z**2
      n(19)=y**4
      n(20)=y**4*z
      n(21)=y**5
      n(22)=x
      n(23)=x*z
      n(24)=x*z**2
      n(25)=x*z**3
      n(26)=x*z**4
      n(27)=x*y
      n(28)=x*y*z
      n(29)=x*y*z**2
      n(30)=x*y*z**3
      n(31)=x*y**2
      n(32)=x*y**2*z
      n(33)=x*y**2*z**2
      n(34)=x*y**3
      n(35)=x*y**3*z
      n(36)=x*y**4
      n(37)=x**2
      n(38)=x**2*z
      n(39)=x**2*z**2
      n(40)=x**2*z**3
      n(41)=x**2*y
      n(42)=x**2*y*z
      n(43)=x**2*y*z**2
      n(44)=x**2*y**2
      n(45)=x**2*y**2*z
      n(46)=x**2*y**3
      n(47)=x**3
      n(48)=x**3*z
      n(49)=x**3*z**2
      n(50)=x**3*y
      n(51)=x**3*y*z
      n(52)=x**3*y**2
      n(53)=x**4
      n(54)=x**4*z
      n(55)=x**4*y
      n(56)=x**5

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
      gn(5,2)=0
      gn(5,3)=4*z**3
      gn(6,1)=0
      gn(6,2)=0
      gn(6,3)=5*z**4
      gn(7,1)=0
      gn(7,2)=1
      gn(7,3)=0
      gn(8,1)=0
      gn(8,2)=z
      gn(8,3)=y
      gn(9,1)=0
      gn(9,2)=z**2
      gn(9,3)=2*y*z
      gn(10,1)=0
      gn(10,2)=z**3
      gn(10,3)=3*y*z**2
      gn(11,1)=0
      gn(11,2)=z**4
      gn(11,3)=4*y*z**3
      gn(12,1)=0
      gn(12,2)=2*y
      gn(12,3)=0
      gn(13,1)=0
      gn(13,2)=2*y*z
      gn(13,3)=y**2
      gn(14,1)=0
      gn(14,2)=2*y*z**2
      gn(14,3)=2*y**2*z
      gn(15,1)=0
      gn(15,2)=2*y*z**3
      gn(15,3)=3*y**2*z**2
      gn(16,1)=0
      gn(16,2)=3*y**2
      gn(16,3)=0
      gn(17,1)=0
      gn(17,2)=3*y**2*z
      gn(17,3)=y**3
      gn(18,1)=0
      gn(18,2)=3*y**2*z**2
      gn(18,3)=2*y**3*z
      gn(19,1)=0
      gn(19,2)=4*y**3
      gn(19,3)=0
      gn(20,1)=0
      gn(20,2)=4*y**3*z
      gn(20,3)=y**4
      gn(21,1)=0
      gn(21,2)=5*y**4
      gn(21,3)=0
      gn(22,1)=1
      gn(22,2)=0
      gn(22,3)=0
      gn(23,1)=z
      gn(23,2)=0
      gn(23,3)=x
      gn(24,1)=z**2
      gn(24,2)=0
      gn(24,3)=2*x*z
      gn(25,1)=z**3
      gn(25,2)=0
      gn(25,3)=3*x*z**2
      gn(26,1)=z**4
      gn(26,2)=0
      gn(26,3)=4*x*z**3
      gn(27,1)=y
      gn(27,2)=x
      gn(27,3)=0
      gn(28,1)=y*z
      gn(28,2)=x*z
      gn(28,3)=x*y
      gn(29,1)=y*z**2
      gn(29,2)=x*z**2
      gn(29,3)=2*x*y*z
      gn(30,1)=y*z**3
      gn(30,2)=x*z**3
      gn(30,3)=3*x*y*z**2
      gn(31,1)=y**2
      gn(31,2)=2*x*y
      gn(31,3)=0
      gn(32,1)=y**2*z
      gn(32,2)=2*x*y*z
      gn(32,3)=x*y**2
      gn(33,1)=y**2*z**2
      gn(33,2)=2*x*y*z**2
      gn(33,3)=2*x*y**2*z
      gn(34,1)=y**3
      gn(34,2)=3*x*y**2
      gn(34,3)=0
      gn(35,1)=y**3*z
      gn(35,2)=3*x*y**2*z
      gn(35,3)=x*y**3
      gn(36,1)=y**4
      gn(36,2)=4*x*y**3
      gn(36,3)=0
      gn(37,1)=2*x
      gn(37,2)=0
      gn(37,3)=0
      gn(38,1)=2*x*z
      gn(38,2)=0
      gn(38,3)=x**2
      gn(39,1)=2*x*z**2
      gn(39,2)=0
      gn(39,3)=2*x**2*z
      gn(40,1)=2*x*z**3
      gn(40,2)=0
      gn(40,3)=3*x**2*z**2
      gn(41,1)=2*x*y
      gn(41,2)=x**2
      gn(41,3)=0
      gn(42,1)=2*x*y*z
      gn(42,2)=x**2*z
      gn(42,3)=x**2*y
      gn(43,1)=2*x*y*z**2
      gn(43,2)=x**2*z**2
      gn(43,3)=2*x**2*y*z
      gn(44,1)=2*x*y**2
      gn(44,2)=2*x**2*y
      gn(44,3)=0
      gn(45,1)=2*x*y**2*z
      gn(45,2)=2*x**2*y*z
      gn(45,3)=x**2*y**2
      gn(46,1)=2*x*y**3
      gn(46,2)=3*x**2*y**2
      gn(46,3)=0
      gn(47,1)=3*x**2
      gn(47,2)=0
      gn(47,3)=0
      gn(48,1)=3*x**2*z
      gn(48,2)=0
      gn(48,3)=x**3
      gn(49,1)=3*x**2*z**2
      gn(49,2)=0
      gn(49,3)=2*x**3*z
      gn(50,1)=3*x**2*y
      gn(50,2)=x**3
      gn(50,3)=0
      gn(51,1)=3*x**2*y*z
      gn(51,2)=x**3*z
      gn(51,3)=x**3*y
      gn(52,1)=3*x**2*y**2
      gn(52,2)=2*x**3*y
      gn(52,3)=0
      gn(53,1)=4*x**3
      gn(53,2)=0
      gn(53,3)=0
      gn(54,1)=4*x**3*z
      gn(54,2)=0
      gn(54,3)=x**4
      gn(55,1)=4*x**3*y
      gn(55,2)=x**4
      gn(55,3)=0
      gn(56,1)=5*x**4
      gn(56,2)=0
      gn(56,3)=0

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
      ggn(5,3,3)=12*z**2
      ggn(6,1,1)=0
      ggn(6,1,2)=0
      ggn(6,1,3)=0
      ggn(6,2,1)=0
      ggn(6,2,2)=0
      ggn(6,2,3)=0
      ggn(6,3,1)=0
      ggn(6,3,2)=0
      ggn(6,3,3)=20*z**3
      ggn(7,1,1)=0
      ggn(7,1,2)=0
      ggn(7,1,3)=0
      ggn(7,2,1)=0
      ggn(7,2,2)=0
      ggn(7,2,3)=0
      ggn(7,3,1)=0
      ggn(7,3,2)=0
      ggn(7,3,3)=0
      ggn(8,1,1)=0
      ggn(8,1,2)=0
      ggn(8,1,3)=0
      ggn(8,2,1)=0
      ggn(8,2,2)=0
      ggn(8,2,3)=1
      ggn(8,3,1)=0
      ggn(8,3,2)=1
      ggn(8,3,3)=0
      ggn(9,1,1)=0
      ggn(9,1,2)=0
      ggn(9,1,3)=0
      ggn(9,2,1)=0
      ggn(9,2,2)=0
      ggn(9,2,3)=2*z
      ggn(9,3,1)=0
      ggn(9,3,2)=2*z
      ggn(9,3,3)=2*y
      ggn(10,1,1)=0
      ggn(10,1,2)=0
      ggn(10,1,3)=0
      ggn(10,2,1)=0
      ggn(10,2,2)=0
      ggn(10,2,3)=3*z**2
      ggn(10,3,1)=0
      ggn(10,3,2)=3*z**2
      ggn(10,3,3)=6*y*z
      ggn(11,1,1)=0
      ggn(11,1,2)=0
      ggn(11,1,3)=0
      ggn(11,2,1)=0
      ggn(11,2,2)=0
      ggn(11,2,3)=4*z**3
      ggn(11,3,1)=0
      ggn(11,3,2)=4*z**3
      ggn(11,3,3)=12*y*z**2
      ggn(12,1,1)=0
      ggn(12,1,2)=0
      ggn(12,1,3)=0
      ggn(12,2,1)=0
      ggn(12,2,2)=2
      ggn(12,2,3)=0
      ggn(12,3,1)=0
      ggn(12,3,2)=0
      ggn(12,3,3)=0
      ggn(13,1,1)=0
      ggn(13,1,2)=0
      ggn(13,1,3)=0
      ggn(13,2,1)=0
      ggn(13,2,2)=2*z
      ggn(13,2,3)=2*y
      ggn(13,3,1)=0
      ggn(13,3,2)=2*y
      ggn(13,3,3)=0
      ggn(14,1,1)=0
      ggn(14,1,2)=0
      ggn(14,1,3)=0
      ggn(14,2,1)=0
      ggn(14,2,2)=2*z**2
      ggn(14,2,3)=4*y*z
      ggn(14,3,1)=0
      ggn(14,3,2)=4*y*z
      ggn(14,3,3)=2*y**2
      ggn(15,1,1)=0
      ggn(15,1,2)=0
      ggn(15,1,3)=0
      ggn(15,2,1)=0
      ggn(15,2,2)=2*z**3
      ggn(15,2,3)=6*y*z**2
      ggn(15,3,1)=0
      ggn(15,3,2)=6*y*z**2
      ggn(15,3,3)=6*y**2*z
      ggn(16,1,1)=0
      ggn(16,1,2)=0
      ggn(16,1,3)=0
      ggn(16,2,1)=0
      ggn(16,2,2)=6*y
      ggn(16,2,3)=0
      ggn(16,3,1)=0
      ggn(16,3,2)=0
      ggn(16,3,3)=0
      ggn(17,1,1)=0
      ggn(17,1,2)=0
      ggn(17,1,3)=0
      ggn(17,2,1)=0
      ggn(17,2,2)=6*y*z
      ggn(17,2,3)=3*y**2
      ggn(17,3,1)=0
      ggn(17,3,2)=3*y**2
      ggn(17,3,3)=0
      ggn(18,1,1)=0
      ggn(18,1,2)=0
      ggn(18,1,3)=0
      ggn(18,2,1)=0
      ggn(18,2,2)=6*y*z**2
      ggn(18,2,3)=6*y**2*z
      ggn(18,3,1)=0
      ggn(18,3,2)=6*y**2*z
      ggn(18,3,3)=2*y**3
      ggn(19,1,1)=0
      ggn(19,1,2)=0
      ggn(19,1,3)=0
      ggn(19,2,1)=0
      ggn(19,2,2)=12*y**2
      ggn(19,2,3)=0
      ggn(19,3,1)=0
      ggn(19,3,2)=0
      ggn(19,3,3)=0
      ggn(20,1,1)=0
      ggn(20,1,2)=0
      ggn(20,1,3)=0
      ggn(20,2,1)=0
      ggn(20,2,2)=12*y**2*z
      ggn(20,2,3)=4*y**3
      ggn(20,3,1)=0
      ggn(20,3,2)=4*y**3
      ggn(20,3,3)=0
      ggn(21,1,1)=0
      ggn(21,1,2)=0
      ggn(21,1,3)=0
      ggn(21,2,1)=0
      ggn(21,2,2)=20*y**3
      ggn(21,2,3)=0
      ggn(21,3,1)=0
      ggn(21,3,2)=0
      ggn(21,3,3)=0
      ggn(22,1,1)=0
      ggn(22,1,2)=0
      ggn(22,1,3)=0
      ggn(22,2,1)=0
      ggn(22,2,2)=0
      ggn(22,2,3)=0
      ggn(22,3,1)=0
      ggn(22,3,2)=0
      ggn(22,3,3)=0
      ggn(23,1,1)=0
      ggn(23,1,2)=0
      ggn(23,1,3)=1
      ggn(23,2,1)=0
      ggn(23,2,2)=0
      ggn(23,2,3)=0
      ggn(23,3,1)=1
      ggn(23,3,2)=0
      ggn(23,3,3)=0
      ggn(24,1,1)=0
      ggn(24,1,2)=0
      ggn(24,1,3)=2*z
      ggn(24,2,1)=0
      ggn(24,2,2)=0
      ggn(24,2,3)=0
      ggn(24,3,1)=2*z
      ggn(24,3,2)=0
      ggn(24,3,3)=2*x
      ggn(25,1,1)=0
      ggn(25,1,2)=0
      ggn(25,1,3)=3*z**2
      ggn(25,2,1)=0
      ggn(25,2,2)=0
      ggn(25,2,3)=0
      ggn(25,3,1)=3*z**2
      ggn(25,3,2)=0
      ggn(25,3,3)=6*x*z
      ggn(26,1,1)=0
      ggn(26,1,2)=0
      ggn(26,1,3)=4*z**3
      ggn(26,2,1)=0
      ggn(26,2,2)=0
      ggn(26,2,3)=0
      ggn(26,3,1)=4*z**3
      ggn(26,3,2)=0
      ggn(26,3,3)=12*x*z**2
      ggn(27,1,1)=0
      ggn(27,1,2)=1
      ggn(27,1,3)=0
      ggn(27,2,1)=1
      ggn(27,2,2)=0
      ggn(27,2,3)=0
      ggn(27,3,1)=0
      ggn(27,3,2)=0
      ggn(27,3,3)=0
      ggn(28,1,1)=0
      ggn(28,1,2)=z
      ggn(28,1,3)=y
      ggn(28,2,1)=z
      ggn(28,2,2)=0
      ggn(28,2,3)=x
      ggn(28,3,1)=y
      ggn(28,3,2)=x
      ggn(28,3,3)=0
      ggn(29,1,1)=0
      ggn(29,1,2)=z**2
      ggn(29,1,3)=2*y*z
      ggn(29,2,1)=z**2
      ggn(29,2,2)=0
      ggn(29,2,3)=2*x*z
      ggn(29,3,1)=2*y*z
      ggn(29,3,2)=2*x*z
      ggn(29,3,3)=2*x*y
      ggn(30,1,1)=0
      ggn(30,1,2)=z**3
      ggn(30,1,3)=3*y*z**2
      ggn(30,2,1)=z**3
      ggn(30,2,2)=0
      ggn(30,2,3)=3*x*z**2
      ggn(30,3,1)=3*y*z**2
      ggn(30,3,2)=3*x*z**2
      ggn(30,3,3)=6*x*y*z
      ggn(31,1,1)=0
      ggn(31,1,2)=2*y
      ggn(31,1,3)=0
      ggn(31,2,1)=2*y
      ggn(31,2,2)=2*x
      ggn(31,2,3)=0
      ggn(31,3,1)=0
      ggn(31,3,2)=0
      ggn(31,3,3)=0
      ggn(32,1,1)=0
      ggn(32,1,2)=2*y*z
      ggn(32,1,3)=y**2
      ggn(32,2,1)=2*y*z
      ggn(32,2,2)=2*x*z
      ggn(32,2,3)=2*x*y
      ggn(32,3,1)=y**2
      ggn(32,3,2)=2*x*y
      ggn(32,3,3)=0
      ggn(33,1,1)=0
      ggn(33,1,2)=2*y*z**2
      ggn(33,1,3)=2*y**2*z
      ggn(33,2,1)=2*y*z**2
      ggn(33,2,2)=2*x*z**2
      ggn(33,2,3)=4*x*y*z
      ggn(33,3,1)=2*y**2*z
      ggn(33,3,2)=4*x*y*z
      ggn(33,3,3)=2*x*y**2
      ggn(34,1,1)=0
      ggn(34,1,2)=3*y**2
      ggn(34,1,3)=0
      ggn(34,2,1)=3*y**2
      ggn(34,2,2)=6*x*y
      ggn(34,2,3)=0
      ggn(34,3,1)=0
      ggn(34,3,2)=0
      ggn(34,3,3)=0
      ggn(35,1,1)=0
      ggn(35,1,2)=3*y**2*z
      ggn(35,1,3)=y**3
      ggn(35,2,1)=3*y**2*z
      ggn(35,2,2)=6*x*y*z
      ggn(35,2,3)=3*x*y**2
      ggn(35,3,1)=y**3
      ggn(35,3,2)=3*x*y**2
      ggn(35,3,3)=0
      ggn(36,1,1)=0
      ggn(36,1,2)=4*y**3
      ggn(36,1,3)=0
      ggn(36,2,1)=4*y**3
      ggn(36,2,2)=12*x*y**2
      ggn(36,2,3)=0
      ggn(36,3,1)=0
      ggn(36,3,2)=0
      ggn(36,3,3)=0
      ggn(37,1,1)=2
      ggn(37,1,2)=0
      ggn(37,1,3)=0
      ggn(37,2,1)=0
      ggn(37,2,2)=0
      ggn(37,2,3)=0
      ggn(37,3,1)=0
      ggn(37,3,2)=0
      ggn(37,3,3)=0
      ggn(38,1,1)=2*z
      ggn(38,1,2)=0
      ggn(38,1,3)=2*x
      ggn(38,2,1)=0
      ggn(38,2,2)=0
      ggn(38,2,3)=0
      ggn(38,3,1)=2*x
      ggn(38,3,2)=0
      ggn(38,3,3)=0
      ggn(39,1,1)=2*z**2
      ggn(39,1,2)=0
      ggn(39,1,3)=4*x*z
      ggn(39,2,1)=0
      ggn(39,2,2)=0
      ggn(39,2,3)=0
      ggn(39,3,1)=4*x*z
      ggn(39,3,2)=0
      ggn(39,3,3)=2*x**2
      ggn(40,1,1)=2*z**3
      ggn(40,1,2)=0
      ggn(40,1,3)=6*x*z**2
      ggn(40,2,1)=0
      ggn(40,2,2)=0
      ggn(40,2,3)=0
      ggn(40,3,1)=6*x*z**2
      ggn(40,3,2)=0
      ggn(40,3,3)=6*x**2*z
      ggn(41,1,1)=2*y
      ggn(41,1,2)=2*x
      ggn(41,1,3)=0
      ggn(41,2,1)=2*x
      ggn(41,2,2)=0
      ggn(41,2,3)=0
      ggn(41,3,1)=0
      ggn(41,3,2)=0
      ggn(41,3,3)=0
      ggn(42,1,1)=2*y*z
      ggn(42,1,2)=2*x*z
      ggn(42,1,3)=2*x*y
      ggn(42,2,1)=2*x*z
      ggn(42,2,2)=0
      ggn(42,2,3)=x**2
      ggn(42,3,1)=2*x*y
      ggn(42,3,2)=x**2
      ggn(42,3,3)=0
      ggn(43,1,1)=2*y*z**2
      ggn(43,1,2)=2*x*z**2
      ggn(43,1,3)=4*x*y*z
      ggn(43,2,1)=2*x*z**2
      ggn(43,2,2)=0
      ggn(43,2,3)=2*x**2*z
      ggn(43,3,1)=4*x*y*z
      ggn(43,3,2)=2*x**2*z
      ggn(43,3,3)=2*x**2*y
      ggn(44,1,1)=2*y**2
      ggn(44,1,2)=4*x*y
      ggn(44,1,3)=0
      ggn(44,2,1)=4*x*y
      ggn(44,2,2)=2*x**2
      ggn(44,2,3)=0
      ggn(44,3,1)=0
      ggn(44,3,2)=0
      ggn(44,3,3)=0
      ggn(45,1,1)=2*y**2*z
      ggn(45,1,2)=4*x*y*z
      ggn(45,1,3)=2*x*y**2
      ggn(45,2,1)=4*x*y*z
      ggn(45,2,2)=2*x**2*z
      ggn(45,2,3)=2*x**2*y
      ggn(45,3,1)=2*x*y**2
      ggn(45,3,2)=2*x**2*y
      ggn(45,3,3)=0
      ggn(46,1,1)=2*y**3
      ggn(46,1,2)=6*x*y**2
      ggn(46,1,3)=0
      ggn(46,2,1)=6*x*y**2
      ggn(46,2,2)=6*x**2*y
      ggn(46,2,3)=0
      ggn(46,3,1)=0
      ggn(46,3,2)=0
      ggn(46,3,3)=0
      ggn(47,1,1)=6*x
      ggn(47,1,2)=0
      ggn(47,1,3)=0
      ggn(47,2,1)=0
      ggn(47,2,2)=0
      ggn(47,2,3)=0
      ggn(47,3,1)=0
      ggn(47,3,2)=0
      ggn(47,3,3)=0
      ggn(48,1,1)=6*x*z
      ggn(48,1,2)=0
      ggn(48,1,3)=3*x**2
      ggn(48,2,1)=0
      ggn(48,2,2)=0
      ggn(48,2,3)=0
      ggn(48,3,1)=3*x**2
      ggn(48,3,2)=0
      ggn(48,3,3)=0
      ggn(49,1,1)=6*x*z**2
      ggn(49,1,2)=0
      ggn(49,1,3)=6*x**2*z
      ggn(49,2,1)=0
      ggn(49,2,2)=0
      ggn(49,2,3)=0
      ggn(49,3,1)=6*x**2*z
      ggn(49,3,2)=0
      ggn(49,3,3)=2*x**3
      ggn(50,1,1)=6*x*y
      ggn(50,1,2)=3*x**2
      ggn(50,1,3)=0
      ggn(50,2,1)=3*x**2
      ggn(50,2,2)=0
      ggn(50,2,3)=0
      ggn(50,3,1)=0
      ggn(50,3,2)=0
      ggn(50,3,3)=0
      ggn(51,1,1)=6*x*y*z
      ggn(51,1,2)=3*x**2*z
      ggn(51,1,3)=3*x**2*y
      ggn(51,2,1)=3*x**2*z
      ggn(51,2,2)=0
      ggn(51,2,3)=x**3
      ggn(51,3,1)=3*x**2*y
      ggn(51,3,2)=x**3
      ggn(51,3,3)=0
      ggn(52,1,1)=6*x*y**2
      ggn(52,1,2)=6*x**2*y
      ggn(52,1,3)=0
      ggn(52,2,1)=6*x**2*y
      ggn(52,2,2)=2*x**3
      ggn(52,2,3)=0
      ggn(52,3,1)=0
      ggn(52,3,2)=0
      ggn(52,3,3)=0
      ggn(53,1,1)=12*x**2
      ggn(53,1,2)=0
      ggn(53,1,3)=0
      ggn(53,2,1)=0
      ggn(53,2,2)=0
      ggn(53,2,3)=0
      ggn(53,3,1)=0
      ggn(53,3,2)=0
      ggn(53,3,3)=0
      ggn(54,1,1)=12*x**2*z
      ggn(54,1,2)=0
      ggn(54,1,3)=4*x**3
      ggn(54,2,1)=0
      ggn(54,2,2)=0
      ggn(54,2,3)=0
      ggn(54,3,1)=4*x**3
      ggn(54,3,2)=0
      ggn(54,3,3)=0
      ggn(55,1,1)=12*x**2*y
      ggn(55,1,2)=4*x**3
      ggn(55,1,3)=0
      ggn(55,2,1)=4*x**3
      ggn(55,2,2)=0
      ggn(55,2,3)=0
      ggn(55,3,1)=0
      ggn(55,3,2)=0
      ggn(55,3,3)=0
      ggn(56,1,1)=20*x**3
      ggn(56,1,2)=0
      ggn(56,1,3)=0
      ggn(56,2,1)=0
      ggn(56,2,2)=0
      ggn(56,2,3)=0
      ggn(56,3,1)=0
      ggn(56,3,2)=0
      ggn(56,3,3)=0
      
! B0 bauen
      do i=1,4
      do j=1,4
        B0(i,j)=0.d0
        do m=1,56
          B0(i,j)=B0(i,j)+n(m)*Mi(m,imat(i,j))
        end do
      end do
      end do      
      
! B1 bauen
      do i=1,4
      do j=1,4
      do k=1,3
        B1(i,j,k)=0.d0
        do m=1,56
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
        do m=1,56
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
      
      lambda=20.d0*0.3d0/(1.3d0*0.4d0)
      mu=    20.d0/2.6d0
      alpha= 0.10d0

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
      
      emo=200.d0
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
