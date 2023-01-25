      subroutine wlist_float(string,list,n)
      integer  :: n,i
      real*8   :: list(n)
      character*11 :: tmp
      character*(*) :: string
Cf2py intent(inout) string
Cf2py intent(in,copy) list
Cf2py integer intent(hide),depend(list) :: n=shape(list,0)
      do i = 1,n
         call wreal(list(i),tmp)
         string((i-1)*11+1:i*11) = tmp
      enddo
      end subroutine wlist_float
c
c
c
c       subroutine wlist_int(string,list,n)
c       integer   :: n,i
c       integer*8 :: list(n)
c       character*11 :: tmp
c       character*(*) :: string
c Cf2py intent(inout) string
c Cf2py intent(in,copy) list
c Cf2py integer intent(hide),depend(list) :: n=shape(list,0)
c       do i = 1,n
c          call wint(list(i),tmp)
c          string((i-1)*11+1:i*11) = tmp
c       enddo
c       end subroutine wlist_int
c
c
c
      subroutine wreal(x,a)
      real*8   :: x
      character*11 :: a
      character*13 :: tmp
Cf2py intent(inout) a
Cf2py intent(in) x
      if (abs(x) > 1e-1 .and. abs(x) < 1E1) then
         write(a,'(F11.8)') x
      elseif (abs(x) >= 1E1 .and. abs(x) < 1E2) then
         write(a,'(F11.7)') x
      elseif (abs(x) >= 1E2 .and. abs(x) < 1E3) then
         write(a,'(F11.6)') x
      elseif (abs(x) >= 1E3 .and. abs(x) < 1E4) then
         write(a,'(F11.5)') x
      elseif (abs(x) >= 1E4 .and. abs(x) < 1E5) then
         write(a,'(F11.4)') x
      elseif (abs(x) >= 1E5 .and. abs(x) < 1E6) then
         write(a,'(F11.3)') x
      elseif (abs(x) >= 1E6 .and. abs(x) < 1E7) then
         write(a,'(F11.2)') x
      elseif (abs(x) >= 1E7 .and. abs(x) < 1E8) then
         write(a,'(F11.1)') x
      elseif (abs(x) >= 1E8 .and. abs(x) < 1E9) then
         write(a,'(F11.0)') x
      else
         write(tmp,'(1PE13.6)') x
         if (tmp(12:12) == '0') then
            write(a,'(A11)') tmp(1:9)//tmp(11:11)//tmp(13:13)
         else
            write(a,'(A11)') tmp(1:8)//tmp(11:13)
         endif      
      endif
      end subroutine wreal
c
c
c
c       subroutine wint(x,a)
c       real*8   :: x
c       integer :: y
c       character*11 :: a
c       y = int(x)
c       write(a,'(I11)') y
c       end subroutine wint
c
c
c
      subroutine wcont(string,c1,c2,l1,l2,n1,n2)
      integer  :: l1, l2, n1, n2
      real*8   :: c1, c2
      character*11 :: tmp1, tmp2
      character*66 :: string, string_tmp
Cf2py intent(inout) string
Cf2py integer intent(in) l1, l2, n1, n2
Cf2py intent(in) c1, c2
      call wreal(c1,tmp1)
      call wreal(c2,tmp2)
      write(string_tmp,10) tmp1,tmp2,l1,l2,n1,n2
  10  format (2A11, 4I11)
      string(1:66) = string_tmp(1:66)
      end subroutine wcont
c
c
c
      subroutine wilist(string,list,n)
      integer  :: n,i
      integer*8   :: list(n)
      character*11 :: tmp
      character*(*) :: string
Cf2py intent(inout) string
Cf2py intent(in,copy) list
Cf2py integer intent(hide),depend(list) :: n=shape(list,0)
      do i = 1,n
         write(tmp,'(I11)') list(i)
         string((i-1)*11+1:i*11) = tmp
      enddo
      end subroutine wilist
c
c
c
      subroutine wlist(string,list,n)
      integer  :: n,i
      real*8   :: list(n)
      character*11 :: tmp
      character*(*) :: string
Cf2py intent(inout) string
Cf2py intent(in,copy) list
Cf2py integer intent(hide),depend(list) :: n=shape(list,0)
      do i = 1,n
         call wreal(list(i), tmp)
         string((i-1)*11+1:i*11) = tmp
      enddo
      end subroutine wlist
c
c
c
      subroutine rcont(string,io_status,c1,c2,l1,l2,n1,n2)
      real*8 :: c1,c2
      integer:: l1,l2,n1,n2,io_status
      character*80 :: string
Cf2py intent(inout) string
Cf2py intent(inout) c1,c2,l1,l2,n1,n2,io_status
  10  format (2E11.0, 4I11)
      read(string, 10, iostat=io_status) c1,c2,l1,l2,n1,n2
      end subroutine rcont
c
c
c
      subroutine rlist(string,io_status,list,n)
      integer  :: n,i,io_status
      real*8   :: list(n)
      character*80 :: string
Cf2py intent(inout) string
Cf2py intent(inout) io_status
Cf2py intent(inout,copy) list
Cf2py integer intent(hide),depend(list) :: n=shape(list,0)
  10  format (6E11.0)
      read(string,10,iostat=io_status) (list(i), i=1,n)
      end subroutine rlist
c
c
c
      subroutine rilist(string,io_status,ilist,n)
      integer  :: n,i,io_status
      integer*4 :: ilist(n)
      character*80 :: string
Cf2py intent(inout) string
Cf2py intent(inout) io_status
Cf2py intent(inout,copy) ilist
Cf2py integer intent(hide),depend(ilist) :: n=shape(ilist,0)
  10  format (6I11)
      read(string,10,iostat=io_status) (ilist(i), i=1,n)
      end subroutine rilist
c
c
c
      subroutine rcontrol(string,io_status,mat,mf,mt,ns)
      integer:: mat,mf,mt,ns,io_status
      character*80 :: string
Cf2py intent(inout) string
Cf2py intent(inout) mat,mf,mt,ns,io_status
  10  format (I4, I2, I3, I5)
      read(string,10,iostat=io_status) mat,mf,mt,ns
      end subroutine rcontrol
c
c
c
      subroutine r_exp_format(string, x, io_status)
      real*8       :: x
      character*11 :: string
Cf2py intent(inout) string
Cf2py intent(inout) x
  10  format (E11.0)
      read(string,10,iostat=io_status) x
      end subroutine r_exp_format
c
c
c
      SUBROUTINE lcomp2
      implicit real*8 (a-h,o-z)
      implicit integer (i-n)
C CORR is the full correlation matrix of dimensions NNN x NNN
C MXCOR is the maximum dimension of CORR (NNN.LE.MXCOR)
      PARAMETER (MXCOR=1000)
      DIMENSION KIJ(18), CORR(MXCOR*MXCOR)
C Read the CONT record:
C NNN is the dimension of CORR(NNN,NNN),
C NM is the number of lines to follow in the file
C NDIGIT is the number of digits for the covariance matrix
      LIB = 500
      open(LIB, file='lcomp2.tmp')
      READ (LIB,10) C1,C2,NDIGIT,NNN, NM, NX, MAT, MF, MT, NS
  10  FORMAT (2F11.0, 4I11, I4, I2, I3,I5)
      IF(NNN.GT.MXCOR) STOP 'MXCOR Limit exceeded'
C Preset the correlation matrix to zero
      NN2=NNN*NNN
      DO I=1,NN2
         CORR(I)=0
      END DO
C Preset the diagonal to one
      DO I=1,NNN
         CORR(I+(I-1)*NNN)=1
      END DO
      DO M=1,NM
C Read the INTG record
         IF(NDIGIT.EQ.2) THEN
             NROW=18
             READ (LIB,20) II, JJ, (KIJ(N),N=1,NROW), MAT, MF, MT, MS
  20         FORMAT (I5, I5, 1X, 18I3, 1X, I4, I2, I3, I5)
         ELSE IF(NDIGIT.EQ.3) THEN
             NROW=13
             READ (LIB,30) II, JJ, (KIJ(N),N=1,NROW), MAT, MF, MT, MS
  30         FORMAT (I5, I5, 1X, 13I4, 3X, I4, I2, I3, I5)
         ELSE IF(NDIGIT.EQ.4) THEN
             NROW=11
             READ (LIB,40) II, JJ, (KIJ(N),N=1,NROW), MAT, MF, MT, MS
  40         FORMAT (I5, I5, 1X, 11I5,I4, I2, I3, I5)
         ELSE IF(NDIGIT.EQ.5) THEN
             NROW= 9
             READ (LIB,50) II, JJ, (KIJ(N),N=1,NROW), MAT, MF, MT, MS
  50         FORMAT (I5, I5, 1X, 9I6, 1X, I4, I2, I3, I5)
         ELSE IF(NDIGIT.EQ.6) THEN
             NROW= 8
             READ (LIB,60) II, JJ, (KIJ(N),N=1,NROW), MAT, MF, MT, MS
  60         FORMAT (I5, I5, 8I7, I4, I2, I3, I5)
         ELSE
             STOP 'ERROR - Invalid NDIGIT'
         END IF
C Interpret the INTG record and fill the covariance matrix
         JP = JJ - 1
         Factor =10**(NDIGIT)
         DO N=1,NROW
             JP = JP + 1
             IF(JP.GE.II) GO TO 35
             IF(KIJ(N).NE.0) THEN
                 IF(KIJ(N).GT.0) THEN
                     CIJ = ( KIJ(N)+0.5)/Factor
                 ELSE
                     CIJ =-(-KIJ(N)+0.5)/Factor
                 END IF
                 JJII = JJ + (II-1)*NNN
                 IIJJ = II + (JJ-1)*NNN
                 CORR(JJII)= CIJ
                 CORR(IIJJ)= CIJ
             END IF
         END DO
  35     CONTINUE
      END DO
      end subroutine lcomp2
c
c
c
      subroutine lcomp(NDIGIT, STRING, io_status, II, JJ, NROW, KIJ, N)
      integer  :: N,M,NDIGIT,II,JJ,NROW,MAT,MF,MT,MS, io_status
      real*8  :: KIJ(N)
      integer  :: TMP(N)
      character*80 :: STRING
Cf2py intent(inout) NROW, II, JJ, io_status
Cf2py intent(inout) STRING
Cf2py intent(inout,copy) KIJ
Cf2py integer intent(hide),depend(KIJ) :: N=shape(KIJ,0)
      IF(NDIGIT.EQ.2) THEN
         NROW=18
         READ (STRING,20,iostat=io_status) II, JJ, (TMP(M),M=1,NROW),
     & MAT, MF, MT, MS
  20     FORMAT (I5, I5, 1X, 18I3, 1X, I4, I2, I3, I5)
      ELSE IF(NDIGIT.EQ.3) THEN
         NROW=13
         READ (STRING,35,iostat=io_status) II, JJ, (TMP(M),M=1,NROW),
     & MAT, MF, MT, MS
  35     FORMAT (I5, I5, 1X, 13I4, 3X, I4, I2, I3, I5)
      ELSE IF(NDIGIT.EQ.4) THEN
         NROW=11
         READ (STRING,40,iostat=io_status) II, JJ, (TMP(M),M=1,NROW),
     & MAT, MF, MT, MS
  40     FORMAT (I5, I5, 1X, 11I5,I4, I2, I3, I5)
      ELSE IF(NDIGIT.EQ.5) THEN
         NROW= 9
         READ (STRING,50,iostat=io_status) II, JJ, (TMP(M),M=1,NROW),
     & MAT, MF, MT, MS
  50     FORMAT (I5, I5, 1X, 9I6, 1X, I4, I2, I3, I5)
      ELSE IF(NDIGIT.EQ.6) THEN
         NROW= 8
         READ (STRING,60,iostat=io_status) II, JJ, (TMP(M),M=1,NROW),
     & MAT, MF, MT, MS
  60     FORMAT (I5, I5, 8I7, I4, I2, I3, I5)
      ELSE
         STOP 'ERROR - Invalid NDIGIT'
      END IF
      do M = 1,NROW
         KIJ(M) = real(TMP(M))
      enddo
      end subroutine lcomp
