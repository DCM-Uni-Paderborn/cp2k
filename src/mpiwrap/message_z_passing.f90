! *****************************************************************************
!> \brief Shift around the data in msg
!> \param[in,out] msg         Rank-2 data to shift
!> \param[in] group           message passing environment identifier
!> \param[in] displ_in        displacements (?)
!> \par Example
!>      msg will be moved from rank to rank+displ_in (in a circular way)
!> \par Limitations
!>      * displ_in will be 1 by default (others not tested)
!>      * the message array needs to be the same size on all processes
! *****************************************************************************
  SUBROUTINE mp_shift_zm( msg, group, displ_in)

    COMPLEX(kind=real_8), INTENT(INOUT)                   :: msg( :, : )
    INTEGER, INTENT(IN)                      :: group
    INTEGER, INTENT(IN), OPTIONAL            :: displ_in

    CHARACTER(len=*), PARAMETER :: routineN = 'mp_shift_zm', &
      routineP = moduleN//':'//routineN

    INTEGER                                  :: handle, ierror
#if defined(__parallel)
    INTEGER                                  :: displ, left, &
                                                msglen, myrank, nprocs, &
                                                right, tag
    INTEGER, ALLOCATABLE, DIMENSION(:)       :: status
#endif

    ierror = 0
    CALL mp_timeset(routineN,handle)

#if defined(__parallel)
    ALLOCATE(status(MPI_STATUS_SIZE))
    CALL mpi_comm_rank(group,myrank,ierror)
    IF ( ierror /= 0 ) CALL mp_stop ( ierror, "mpi_comm_rank @ "//routineN )
    CALL mpi_comm_size(group,nprocs,ierror)
    IF ( ierror /= 0 ) CALL mp_stop ( ierror, "mpi_comm_size @ "//routineN )
    IF (PRESENT(displ_in)) THEN
       displ=displ_in
    ELSE
       displ=1
    ENDIF
    right=MODULO(myrank+displ,nprocs)
    left =MODULO(myrank-displ,nprocs)
    tag=17
    msglen = SIZE(msg)
    t_start = m_walltime ( )
    CALL mpi_sendrecv_replace(msg,msglen,MPI_DOUBLE_COMPLEX,right,tag,left,tag, &
         group,status(1),ierror)
    t_end = m_walltime ( )
    IF ( ierror /= 0 ) CALL mp_stop ( ierror, "mpi_sendrecv_replace @ "//routineN )
    CALL add_perf(perf_id=7,count=1,time=t_end-t_start,msg_size=msglen*(2*real_8_size))
    DEALLOCATE(status)
#endif
    CALL mp_timestop(handle)

  END SUBROUTINE mp_shift_zm

! *****************************************************************************
!> \brief Shift around the data in msg
!> \param[in,out] msg         Data to shift
!> \param[in] group           message passing environment identifier
!> \param[in] displ_in        displacements (?)
!> \par Example
!>      msg will be moved from rank to rank+displ_in (in a circular way)
!> \par Limitations
!>      * displ_in will be 1 by default (others not tested)
!>      * the message array needs to be the same size on all processes
! *****************************************************************************
  SUBROUTINE mp_shift_z( msg, group, displ_in)

    COMPLEX(kind=real_8), INTENT(INOUT)                   :: msg( : )
    INTEGER, INTENT(IN)                      :: group
    INTEGER, INTENT(IN), OPTIONAL            :: displ_in

    CHARACTER(len=*), PARAMETER :: routineN = 'mp_shift_z', &
      routineP = moduleN//':'//routineN

    INTEGER                                  :: handle, ierror
#if defined(__parallel)
    INTEGER                                  :: displ, left, &
                                                msglen, myrank, nprocs, &
                                                right, tag
    INTEGER, ALLOCATABLE, DIMENSION(:)       :: status
#endif

    ierror = 0
    CALL mp_timeset(routineN,handle)

#if defined(__parallel)
    ALLOCATE(status(MPI_STATUS_SIZE))
    CALL mpi_comm_rank(group,myrank,ierror)
    IF ( ierror /= 0 ) CALL mp_stop ( ierror, "mpi_comm_rank @ "//routineN )
    CALL mpi_comm_size(group,nprocs,ierror)
    IF ( ierror /= 0 ) CALL mp_stop ( ierror, "mpi_comm_size @ "//routineN )
    IF (PRESENT(displ_in)) THEN
       displ=displ_in
    ELSE
       displ=1
    ENDIF
    right=MODULO(myrank+displ,nprocs)
    left =MODULO(myrank-displ,nprocs)
    tag=19
    msglen = SIZE(msg)
    t_start = m_walltime ( )
    CALL mpi_sendrecv_replace(msg,msglen,MPI_DOUBLE_COMPLEX,right,tag,left,&
         tag,group,status(1),ierror)
    t_end = m_walltime ( )
    IF ( ierror /= 0 ) CALL mp_stop ( ierror, "mpi_sendrecv_replace @ "//routineN )
    CALL add_perf(perf_id=7,count=1,time=t_end-t_start,msg_size=msglen*(2*real_8_size))
    DEALLOCATE(status)
#endif
    CALL mp_timestop(handle)

  END SUBROUTINE mp_shift_z

! *****************************************************************************
!> \brief All-to-all data exchange, rank-1 data of different sizes
!> \param[in] sb              Data to send
!> \param[in] scount          Data counts for data sent to other processes
!> \param[in] sdispl          Respective data offsets for data sent to process
!> \param[in,out] rb          Buffer into which to receive data
!> \param[in] rcount          Data counts for data received from other
!>                            processes
!> \param[in] rdispl          Respective data offsets for data received from
!>                            other processes
!> \param[in] group           Message passing environment identifier
!> \par MPI mapping
!>      mpi_alltoallv
!> \par Array sizes
!>      The scount, rcount, and the sdispl and rdispl arrays have a
!>      size equal to the number of processes.
!> \par Offsets
!>      Values in sdispl and rdispl start with 0.
! *****************************************************************************
  SUBROUTINE mp_alltoall_z11v ( sb, scount, sdispl, rb, rcount, rdispl, group )

    COMPLEX(kind=real_8), DIMENSION(:), INTENT(IN)        :: sb
    INTEGER, DIMENSION(:), INTENT(IN)        :: scount, sdispl
    COMPLEX(kind=real_8), DIMENSION(:), INTENT(INOUT)     :: rb
    INTEGER, DIMENSION(:), INTENT(IN)        :: rcount, rdispl
    INTEGER, INTENT(IN)                      :: group

    CHARACTER(len=*), PARAMETER :: routineN = 'mp_alltoall_z11v', &
      routineP = moduleN//':'//routineN

    INTEGER                                  :: handle, ierr
#if defined(__parallel)
    INTEGER                                  :: msglen
#else
    INTEGER                                  :: i
#endif

    CALL mp_timeset(routineN,handle)

    ierr = 0
#if defined(__parallel)
    t_start = m_walltime ( )
    CALL mpi_alltoallv ( sb, scount, sdispl, MPI_DOUBLE_COMPLEX, &
         rb, rcount, rdispl, MPI_DOUBLE_COMPLEX, group, ierr )
    IF ( ierr /= 0 ) CALL mp_stop( ierr, "mpi_alltoallv @ "//routineN )
    t_end = m_walltime ( )
    msglen = SUM ( scount ) + SUM ( rcount )
    CALL add_perf(perf_id=6,count=1,time=t_end-t_start,msg_size=msglen*(2*real_8_size))
#else
    !$OMP PARALLEL DO DEFAULT(NONE) PRIVATE(i) SHARED(rcount,rdispl,sdispl,rb,sb)
    DO i=1,rcount(1)
       rb(rdispl(1)+i)=sb(sdispl(1)+i)
    ENDDO
#endif
    CALL mp_timestop(handle)

  END SUBROUTINE mp_alltoall_z11v

! *****************************************************************************
!> \brief All-to-all data exchange, rank-2 data of different sizes
!> \param sb ...
!> \param scount ...
!> \param sdispl ...
!> \param rb ...
!> \param rcount ...
!> \param rdispl ...
!> \param group ...
!> \par MPI mapping
!>      mpi_alltoallv
!> \note see mp_alltoall_z11v 
! *****************************************************************************
  SUBROUTINE mp_alltoall_z22v ( sb, scount, sdispl, rb, rcount, rdispl, group )

    COMPLEX(kind=real_8), DIMENSION(:, :), &
      INTENT(IN)                             :: sb
    INTEGER, DIMENSION(:), INTENT(IN)        :: scount, sdispl
    COMPLEX(kind=real_8), DIMENSION(:, :), &
      INTENT(INOUT)                          :: rb
    INTEGER, DIMENSION(:), INTENT(IN)        :: rcount, rdispl
    INTEGER, INTENT(IN)                      :: group

    CHARACTER(len=*), PARAMETER :: routineN = 'mp_alltoall_z22v', &
      routineP = moduleN//':'//routineN

    INTEGER                                  :: handle, ierr
#if defined(__parallel)
    INTEGER                                  :: msglen
#endif

    ierr = 0
    CALL mp_timeset(routineN,handle)

#if defined(__parallel)
    t_start = m_walltime ( )
    CALL mpi_alltoallv ( sb, scount, sdispl, MPI_DOUBLE_COMPLEX, &
         rb, rcount, rdispl, MPI_DOUBLE_COMPLEX, group, ierr )
    IF ( ierr /= 0 ) CALL mp_stop( ierr, "mpi_alltoallv @ "//routineN )
    msglen = SUM ( scount ) + SUM ( rcount )
    t_end = m_walltime ( )
    CALL add_perf(perf_id=6,count=1,time=t_end-t_start,msg_size=msglen*2*(2*real_8_size))
#else
    rb=sb
#endif
    CALL mp_timestop(handle)

  END SUBROUTINE mp_alltoall_z22v

! *****************************************************************************
!> \brief All-to-all data exchange, rank 1 arrays, equal sizes
!> \param[in] sb    array with data to send
!> \param[out] rb   array into which data is received
!> \param[in] count  number of elements to send/receive (product of the
!>                   extents of the first two dimensions)
!> \param[in] group           Message passing environment identifier
!> \par Index meaning
!> \par The first two indices specify the data while the last index counts
!>      the processes
!> \par Sizes of ranks
!>      All processes have the same data size.
!> \par MPI mapping
!>      mpi_alltoall
! *****************************************************************************
  SUBROUTINE mp_alltoall_z ( sb, rb, count, group )

    COMPLEX(kind=real_8), DIMENSION(:), INTENT(IN)        :: sb
    COMPLEX(kind=real_8), DIMENSION(:), INTENT(OUT)       :: rb
    INTEGER, INTENT(IN)                      :: count, group

    CHARACTER(len=*), PARAMETER :: routineN = 'mp_alltoall_z', &
      routineP = moduleN//':'//routineN

    INTEGER                                  :: handle, ierr
#if defined(__parallel)
    INTEGER                                  :: msglen, np
#endif

    ierr = 0
    CALL mp_timeset(routineN,handle)

#if defined(__parallel)
    t_start = m_walltime ( )
    CALL mpi_alltoall ( sb, count, MPI_DOUBLE_COMPLEX, &
         rb, count, MPI_DOUBLE_COMPLEX, group, ierr )
    IF ( ierr /= 0 ) CALL mp_stop( ierr, "mpi_alltoall @ "//routineN )
    CALL mpi_comm_size ( group, np, ierr )
    IF ( ierr /= 0 ) CALL mp_stop ( ierr, "mpi_comm_size @ "//routineN )
    msglen = 2 * count * np
    t_end = m_walltime ( )
    CALL add_perf(perf_id=6,count=1,time=t_end-t_start,msg_size=msglen*(2*real_8_size))
#else
    rb=sb
#endif
    CALL mp_timestop(handle)

  END SUBROUTINE mp_alltoall_z

! *****************************************************************************
!> \brief All-to-all data exchange, rank-2 arrays, equal sizes
!> \param sb ...
!> \param rb ...
!> \param count ...
!> \param group ...
!> \note see mp_alltoall_z 
! *****************************************************************************
  SUBROUTINE mp_alltoall_z22 ( sb, rb, count, group )

    COMPLEX(kind=real_8), DIMENSION(:, :), INTENT(IN)     :: sb
    COMPLEX(kind=real_8), DIMENSION(:, :), INTENT(OUT)    :: rb
    INTEGER, INTENT(IN)                      :: count, group

    CHARACTER(len=*), PARAMETER :: routineN = 'mp_alltoall_z22', &
      routineP = moduleN//':'//routineN

    INTEGER                                  :: handle, ierr
#if defined(__parallel)
    INTEGER                                  :: msglen, np
#endif

    ierr = 0
    CALL mp_timeset(routineN,handle)

#if defined(__parallel)
    t_start = m_walltime ( )
    CALL mpi_alltoall ( sb, count, MPI_DOUBLE_COMPLEX, &
         rb, count, MPI_DOUBLE_COMPLEX, group, ierr )
    IF ( ierr /= 0 ) CALL mp_stop( ierr, "mpi_alltoall @ "//routineN )
    CALL mpi_comm_size ( group, np, ierr )
    IF ( ierr /= 0 ) CALL mp_stop ( ierr, "mpi_comm_size @ "//routineN )
    msglen = 2 * SIZE(sb) * np
    t_end = m_walltime ( )
    CALL add_perf(perf_id=6,count=1,time=t_end-t_start,msg_size=msglen*(2*real_8_size))
#else
    rb=sb
#endif
    CALL mp_timestop(handle)

  END SUBROUTINE mp_alltoall_z22

! *****************************************************************************
!> \brief All-to-all data exchange, rank-3 data with equal sizes
!> \param sb ...
!> \param rb ...
!> \param count ...
!> \param group ...
!> \note see mp_alltoall_z 
! *****************************************************************************
  SUBROUTINE mp_alltoall_z33 ( sb, rb, count, group )

    COMPLEX(kind=real_8), DIMENSION(:, :, :), INTENT(IN)  :: sb
    COMPLEX(kind=real_8), DIMENSION(:, :, :), INTENT(OUT) :: rb
    INTEGER, INTENT(IN)                      :: count, group

    CHARACTER(len=*), PARAMETER :: routineN = 'mp_alltoall_z33', &
      routineP = moduleN//':'//routineN

    INTEGER                                  :: handle, ierr
#if defined(__parallel)
    INTEGER                                  :: msglen, np
#endif

    ierr = 0
    CALL mp_timeset(routineN,handle)

#if defined(__parallel)
    t_start = m_walltime ( )
    CALL mpi_alltoall ( sb, count, MPI_DOUBLE_COMPLEX, &
         rb, count, MPI_DOUBLE_COMPLEX, group, ierr )
    IF ( ierr /= 0 ) CALL mp_stop( ierr, "mpi_alltoall @ "//routineN )
    CALL mpi_comm_size ( group, np, ierr )
    IF ( ierr /= 0 ) CALL mp_stop ( ierr, "mpi_comm_size @ "//routineN )
    msglen = 2 * count * np
    t_end = m_walltime ( )
    CALL add_perf(perf_id=6,count=1,time=t_end-t_start,msg_size=msglen*(2*real_8_size))
#else
    rb=sb
#endif
    CALL mp_timestop(handle)

  END SUBROUTINE mp_alltoall_z33

! *****************************************************************************
!> \brief All-to-all data exchange, rank 4 data, equal sizes
!> \param sb ...
!> \param rb ...
!> \param count ...
!> \param group ...
!> \note see mp_alltoall_z 
! *****************************************************************************
  SUBROUTINE mp_alltoall_z44 ( sb, rb, count, group )

    COMPLEX(kind=real_8), DIMENSION(:, :, :, :), &
      INTENT(IN)                             :: sb
    COMPLEX(kind=real_8), DIMENSION(:, :, :, :), &
      INTENT(OUT)                            :: rb
    INTEGER, INTENT(IN)                      :: count, group

    CHARACTER(len=*), PARAMETER :: routineN = 'mp_alltoall_z44', &
      routineP = moduleN//':'//routineN

    INTEGER                                  :: handle, ierr
#if defined(__parallel)
    INTEGER                                  :: msglen, np
#endif

    ierr = 0
    CALL mp_timeset(routineN,handle)

#if defined(__parallel)
    t_start = m_walltime ( )
    CALL mpi_alltoall ( sb, count, MPI_DOUBLE_COMPLEX, &
         rb, count, MPI_DOUBLE_COMPLEX, group, ierr )
    IF ( ierr /= 0 ) CALL mp_stop( ierr, "mpi_alltoall @ "//routineN )
    CALL mpi_comm_size ( group, np, ierr )
    IF ( ierr /= 0 ) CALL mp_stop ( ierr, "mpi_comm_size @ "//routineN )
    msglen = 2 * count * np
    t_end = m_walltime ( )
    CALL add_perf(perf_id=6,count=1,time=t_end-t_start,msg_size=msglen*(2*real_8_size))
#else
    rb=sb
#endif
    CALL mp_timestop(handle)

  END SUBROUTINE mp_alltoall_z44

! *****************************************************************************
!> \brief All-to-all data exchange, rank-4 data to rank-5 data
!> \param sb ...
!> \param rb ...
!> \param count ...
!> \param group ...
!> \note see mp_alltoall_z 
!> \note User must ensure size consistency.
! *****************************************************************************
  SUBROUTINE mp_alltoall_z45 ( sb, rb, count, group )

    COMPLEX(kind=real_8), DIMENSION(:, :, :, :), &
      INTENT(IN)                             :: sb
    COMPLEX(kind=real_8), &
      DIMENSION(:, :, :, :, :), INTENT(OUT)  :: rb
    INTEGER, INTENT(IN)                      :: count, group

    CHARACTER(len=*), PARAMETER :: routineN = 'mp_alltoall_z45', &
      routineP = moduleN//':'//routineN

    INTEGER                                  :: handle, ierr
#if defined(__parallel)
    INTEGER                                  :: msglen, np
#endif

    ierr = 0
    CALL mp_timeset(routineN,handle)

#if defined(__parallel)
    t_start = m_walltime ( )
    CALL mpi_alltoall ( sb, count, MPI_DOUBLE_COMPLEX, &
         rb, count, MPI_DOUBLE_COMPLEX, group, ierr )
    IF ( ierr /= 0 ) CALL mp_stop( ierr, "mpi_alltoall @ "//routineN )
    CALL mpi_comm_size ( group, np, ierr )
    IF ( ierr /= 0 ) CALL mp_stop ( ierr, "mpi_comm_size @ "//routineN )
    msglen = 2 * count * np
    t_end = m_walltime ( )
    CALL add_perf(perf_id=6,count=1,time=t_end-t_start,msg_size=msglen*(2*real_8_size))
#endif
    CALL mp_timestop(handle)

  END SUBROUTINE mp_alltoall_z45

! *****************************************************************************
!> \brief All-to-all data exchange, rank-3 data to rank-4 data
!> \param sb ...
!> \param rb ...
!> \param count ...
!> \param group ...
!> \note see mp_alltoall_z 
!> \note User must ensure size consistency.
! *****************************************************************************
  SUBROUTINE mp_alltoall_z34 ( sb, rb, count, group )

    COMPLEX(kind=real_8), DIMENSION(:, :, :), &
      INTENT(IN)                             :: sb
    COMPLEX(kind=real_8), DIMENSION(:, :, :, :), &
      INTENT(OUT)                            :: rb
    INTEGER, INTENT(IN)                      :: count, group

    CHARACTER(len=*), PARAMETER :: routineN = 'mp_alltoall_z34', &
      routineP = moduleN//':'//routineN

    INTEGER                                  :: handle, ierr
#if defined(__parallel)
    INTEGER                                  :: msglen, np
#endif

    ierr = 0
    CALL mp_timeset(routineN,handle)

#if defined(__parallel)
    t_start = m_walltime ( )
    CALL mpi_alltoall ( sb, count, MPI_DOUBLE_COMPLEX, &
         rb, count, MPI_DOUBLE_COMPLEX, group, ierr )
    IF ( ierr /= 0 ) CALL mp_stop( ierr, "mpi_alltoall @ "//routineN )
    CALL mpi_comm_size ( group, np, ierr )
    IF ( ierr /= 0 ) CALL mp_stop ( ierr, "mpi_comm_size @ "//routineN )
    msglen = 2 * count * np
    t_end = m_walltime ( )
    CALL add_perf(perf_id=6,count=1,time=t_end-t_start,msg_size=msglen*(2*real_8_size))
#endif
    CALL mp_timestop(handle)

  END SUBROUTINE mp_alltoall_z34

! *****************************************************************************
!> \brief All-to-all data exchange, rank-5 data to rank-4 data
!> \param sb ...
!> \param rb ...
!> \param count ...
!> \param group ...
!> \note see mp_alltoall_z 
!> \note User must ensure size consistency.
! *****************************************************************************
  SUBROUTINE mp_alltoall_z54 ( sb, rb, count, group )

    COMPLEX(kind=real_8), &
      DIMENSION(:, :, :, :, :), INTENT(IN)   :: sb
    COMPLEX(kind=real_8), DIMENSION(:, :, :, :), &
      INTENT(OUT)                            :: rb
    INTEGER, INTENT(IN)                      :: count, group

    CHARACTER(len=*), PARAMETER :: routineN = 'mp_alltoall_z54', &
      routineP = moduleN//':'//routineN

    INTEGER                                  :: handle, ierr
#if defined(__parallel)
    INTEGER                                  :: msglen, np
#endif

    ierr = 0
    CALL mp_timeset(routineN,handle)

#if defined(__parallel)
    t_start = m_walltime ( )
    CALL mpi_alltoall ( sb, count, MPI_DOUBLE_COMPLEX, &
         rb, count, MPI_DOUBLE_COMPLEX, group, ierr )
    IF ( ierr /= 0 ) CALL mp_stop( ierr, "mpi_alltoall @ "//routineN )
    CALL mpi_comm_size ( group, np, ierr )
    IF ( ierr /= 0 ) CALL mp_stop ( ierr, "mpi_comm_size @ "//routineN )
    msglen = 2 * count * np
    t_end = m_walltime ( )
    CALL add_perf(perf_id=6,count=1,time=t_end-t_start,msg_size=msglen*(2*real_8_size))
#endif
    CALL mp_timestop(handle)

  END SUBROUTINE mp_alltoall_z54

! *****************************************************************************
!> \brief Send one datum to another process
!> \param[in] msg             Dum to send
!> \param[in] dest            Destination process
!> \param[in] tag             Transfer identifier
!> \param[in] gid             Message passing environment identifier
!> \par MPI mapping
!>      mpi_send
! *****************************************************************************
  SUBROUTINE mp_send_z(msg,dest,tag,gid)
    COMPLEX(kind=real_8)                                  :: msg
    INTEGER                                  :: dest, tag, gid

    CHARACTER(len=*), PARAMETER :: routineN = 'mp_send_z', &
      routineP = moduleN//':'//routineN

    INTEGER                                  :: handle, ierr, msglen

    ierr = 0
    CALL mp_timeset(routineN,handle)

    msglen = 1
#if defined(__parallel)
    t_start = m_walltime ( )
    CALL mpi_send(msg,msglen,MPI_DOUBLE_COMPLEX,dest,tag,gid,ierr)
    IF ( ierr /= 0 ) CALL mp_stop( ierr, "mpi_send @ "//routineN )
    t_end = m_walltime ( )
    CALL add_perf(perf_id=13,count=1,time=t_end-t_start,msg_size=msglen*(2*real_8_size))
#endif
    CALL mp_timestop(handle)
  END SUBROUTINE mp_send_z

! *****************************************************************************
!> \brief Send rank-1 data to another process
!> \param[in] msg             Rank-1 data to send
!> \param dest ...
!> \param tag ...
!> \param gid ...
!> \note see mp_send_z 
! *****************************************************************************
  SUBROUTINE mp_send_zv(msg,dest,tag,gid)
    COMPLEX(kind=real_8)                                  :: msg( : )
    INTEGER                                  :: dest, tag, gid

    CHARACTER(len=*), PARAMETER :: routineN = 'mp_send_zv', &
      routineP = moduleN//':'//routineN

    INTEGER                                  :: handle, ierr, msglen

    ierr = 0
    CALL mp_timeset(routineN,handle)

    msglen = SIZE(msg)
#if defined(__parallel)
    t_start = m_walltime ( )
    CALL mpi_send(msg,msglen,MPI_DOUBLE_COMPLEX,dest,tag,gid,ierr)
    IF ( ierr /= 0 ) CALL mp_stop( ierr, "mpi_send @ "//routineN )
    t_end = m_walltime ( )
    CALL add_perf(perf_id=13,count=1,time=t_end-t_start,msg_size=msglen*(2*real_8_size))
#endif
    CALL mp_timestop(handle)
  END SUBROUTINE mp_send_zv

! *****************************************************************************
!> \brief Receive one datum from another process
!> \param[in,out] msg         Place received data into this variable
!> \param[in,out] source      Process to receieve from
!> \param[in,out] tag         Transfer identifier
!> \param[in] gid             Message passing environment identifier
!> \par MPI mapping
!>      mpi_send
! *****************************************************************************
  SUBROUTINE mp_recv_z(msg,source,tag,gid)
    COMPLEX(kind=real_8), INTENT(INOUT)                   :: msg
    INTEGER, INTENT(INOUT)                   :: source, tag
    INTEGER, INTENT(IN)                      :: gid

    CHARACTER(len=*), PARAMETER :: routineN = 'mp_recv_z', &
      routineP = moduleN//':'//routineN

    INTEGER                                  :: handle, ierr, msglen
#if defined(__parallel)
    INTEGER, ALLOCATABLE, DIMENSION(:)       :: status
#endif

    ierr = 0
    CALL mp_timeset(routineN,handle)

    msglen = 1
#if defined(__parallel)
    ALLOCATE(status(MPI_STATUS_SIZE))
    t_start = m_walltime ( )
    CALL mpi_recv(msg,msglen,MPI_DOUBLE_COMPLEX,source,tag,gid,status,ierr)
    IF ( ierr /= 0 ) CALL mp_stop( ierr, "mpi_recv @ "//routineN )
    t_end = m_walltime ( )
    CALL add_perf(perf_id=14,count=1,time=t_end-t_start,msg_size=msglen*(2*real_8_size))
    source = status(MPI_SOURCE)
    tag = status(MPI_TAG)
    DEALLOCATE(status)
#endif
    CALL mp_timestop(handle)
  END SUBROUTINE mp_recv_z

! *****************************************************************************
!> \brief Receive rank-1 data from another process
!> \param[in,out] msg         Place receieved data into this rank-1 array
!> \param source ...
!> \param tag ...
!> \param gid ...
!> \note see mp_recv_z 
! *****************************************************************************
  SUBROUTINE mp_recv_zv(msg,source,tag,gid)
    COMPLEX(kind=real_8), INTENT(INOUT)                   :: msg( : )
    INTEGER, INTENT(INOUT)                   :: source, tag
    INTEGER, INTENT(IN)                      :: gid

    CHARACTER(len=*), PARAMETER :: routineN = 'mp_recv_zv', &
      routineP = moduleN//':'//routineN

    INTEGER                                  :: handle, ierr, msglen
#if defined(__parallel)
    INTEGER, ALLOCATABLE, DIMENSION(:)       :: status
#endif

    ierr = 0
    CALL mp_timeset(routineN,handle)

    msglen = SIZE(msg)
#if defined(__parallel)
    ALLOCATE(status(MPI_STATUS_SIZE))
    t_start = m_walltime ( )
    CALL mpi_recv(msg,msglen,MPI_DOUBLE_COMPLEX,source,tag,gid,status,ierr)
    IF ( ierr /= 0 ) CALL mp_stop( ierr, "mpi_recv @ "//routineN )
    t_end = m_walltime ( )
    CALL add_perf(perf_id=14,count=1,time=t_end-t_start,msg_size=msglen*(2*real_8_size))
    source = status(MPI_SOURCE)
    tag = status(MPI_TAG)
    DEALLOCATE(status)
#endif
    CALL mp_timestop(handle)
  END SUBROUTINE mp_recv_zv

! *****************************************************************************
!> \brief Broadcasts a datum to all processes.
!> \param[in] msg             Datum to broadcast
!> \param[in] source          Processes which broadcasts
!> \param[in] gid             Message passing environment identifier
!> \par MPI mapping
!>      mpi_bcast
! *****************************************************************************
  SUBROUTINE mp_bcast_z(msg,source,gid)
    COMPLEX(kind=real_8)                                  :: msg
    INTEGER                                  :: source, gid

    CHARACTER(len=*), PARAMETER :: routineN = 'mp_bcast_z', &
      routineP = moduleN//':'//routineN

    INTEGER                                  :: handle, ierr, msglen

    ierr = 0
    CALL mp_timeset(routineN,handle)

    msglen = 1
#if defined(__parallel)
    t_start = m_walltime ( )
    CALL mpi_bcast(msg,msglen,MPI_DOUBLE_COMPLEX,source,gid,ierr)
    IF ( ierr /= 0 ) CALL mp_stop( ierr, "mpi_bcast @ "//routineN )
    t_end = m_walltime ( )
    CALL add_perf(perf_id=2,count=1,time=t_end-t_start,msg_size=msglen*(2*real_8_size))
#endif
    CALL mp_timestop(handle)
  END SUBROUTINE mp_bcast_z

! *****************************************************************************
!> \brief Broadcasts rank-1 data to all processes
!> \param[in] msg             Data to broadcast
!> \param source ...
!> \param gid ...
!> \note see mp_bcast_z1 
! *****************************************************************************
  SUBROUTINE mp_bcast_zv(msg,source,gid)
    COMPLEX(kind=real_8)                                  :: msg( : )
    INTEGER                                  :: source, gid

    CHARACTER(len=*), PARAMETER :: routineN = 'mp_bcast_zv', &
      routineP = moduleN//':'//routineN

    INTEGER                                  :: handle, ierr, msglen

    ierr = 0
    CALL mp_timeset(routineN,handle)

    msglen = SIZE(msg)
#if defined(__parallel)
    t_start = m_walltime ( )
    CALL mpi_bcast(msg,msglen,MPI_DOUBLE_COMPLEX,source,gid,ierr)
    IF ( ierr /= 0 ) CALL mp_stop( ierr, "mpi_bcast @ "//routineN )
    t_end = m_walltime ( )
    CALL add_perf(perf_id=2,count=1,time=t_end-t_start,msg_size=msglen*(2*real_8_size))
#endif
    CALL mp_timestop(handle)
  END SUBROUTINE mp_bcast_zv

! *****************************************************************************
!> \brief Broadcasts rank-2 data to all processes
!> \param[in] msg             Data to broadcast
!> \param source ...
!> \param gid ...
!> \note see mp_bcast_z1 
! *****************************************************************************
  SUBROUTINE mp_bcast_zm(msg,source,gid)
    COMPLEX(kind=real_8)                                  :: msg( :, : )
    INTEGER                                  :: source, gid

    CHARACTER(len=*), PARAMETER :: routineN = 'mp_bcast_im', &
      routineP = moduleN//':'//routineN

    INTEGER                                  :: handle, ierr, msglen

    ierr = 0
    CALL mp_timeset(routineN,handle)

    msglen = SIZE(msg)
#if defined(__parallel)
    t_start = m_walltime ( )
    CALL mpi_bcast(msg,msglen,MPI_DOUBLE_COMPLEX,source,gid,ierr)
    IF ( ierr /= 0 ) CALL mp_stop( ierr, "mpi_bcast @ "//routineN )
    t_end = m_walltime ( )
    CALL add_perf(perf_id=2,count=1,time=t_end-t_start,msg_size=msglen*(2*real_8_size))
#endif
    CALL mp_timestop(handle)
  END SUBROUTINE mp_bcast_zm

! *****************************************************************************
!> \brief Broadcasts rank-3 data to all processes
!> \param[in] msg             Data to broadcast
!> \param source ...
!> \param gid ...
!> \note see mp_bcast_z1 
! *****************************************************************************
  SUBROUTINE mp_bcast_z3(msg,source,gid)
    COMPLEX(kind=real_8)                                  :: msg( :, :, : )
    INTEGER                                  :: source, gid

    CHARACTER(len=*), PARAMETER :: routineN = 'mp_bcast_z3', &
      routineP = moduleN//':'//routineN

    INTEGER                                  :: handle, ierr, msglen

    ierr = 0
    CALL mp_timeset(routineN,handle)

    msglen = SIZE(msg)
#if defined(__parallel)
    t_start = m_walltime ( )
    CALL mpi_bcast(msg,msglen,MPI_DOUBLE_COMPLEX,source,gid,ierr)
    IF ( ierr /= 0 ) CALL mp_stop( ierr, "mpi_bcast @ "//routineN )
    t_end = m_walltime ( )
    CALL add_perf(perf_id=2,count=1,time=t_end-t_start,msg_size=msglen*(2*real_8_size))
#endif
    CALL mp_timestop(handle)
  END SUBROUTINE mp_bcast_z3

! *****************************************************************************
!> \brief Sums a datum from all processes with result left on all processes.
!> \param[in,out] msg         Datum to sum (input) and result (output)
!> \param[in] gid             Message passing environment identifier
!> \par MPI mapping
!>      mpi_allreduce
! *****************************************************************************
  SUBROUTINE mp_sum_z(msg,gid)
    COMPLEX(kind=real_8), INTENT(INOUT)                   :: msg
    INTEGER, INTENT(IN)                      :: gid

    CHARACTER(len=*), PARAMETER :: routineN = 'mp_sum_z', &
      routineP = moduleN//':'//routineN

    INTEGER                                  :: handle, ierr, msglen

    ierr = 0
    CALL mp_timeset(routineN,handle)

    msglen = 1
#if defined(__parallel)
    t_start = m_walltime ( )
    CALL mpi_allreduce(MPI_IN_PLACE,msg,msglen,MPI_DOUBLE_COMPLEX,MPI_SUM,gid,ierr)
    IF ( ierr /= 0 ) CALL mp_stop( ierr, "mpi_allreduce @ "//routineN )
    t_end = m_walltime ( )
    CALL add_perf(perf_id=3,count=1,time=t_end-t_start,msg_size=msglen*(2*real_8_size))
#endif
    CALL mp_timestop(handle)
  END SUBROUTINE mp_sum_z

! *****************************************************************************
!> \brief Element-wise sum of a rank-1 array on all processes.
!> \param[in,out] msg         Vector to sum and result
!> \param gid ...
!> \note see mp_sum_z 
! *****************************************************************************
  SUBROUTINE mp_sum_zv(msg,gid)
    COMPLEX(kind=real_8), INTENT(INOUT)                   :: msg( : )
    INTEGER, INTENT(IN)                      :: gid

    CHARACTER(len=*), PARAMETER :: routineN = 'mp_sum_zv', &
      routineP = moduleN//':'//routineN

    INTEGER                                  :: handle, ierr
#if defined(__parallel)
    INTEGER                                  :: msglen
#endif

    ierr = 0
    CALL mp_timeset(routineN,handle)

#if defined(__parallel)
    t_start = m_walltime ( )
    msglen = SIZE(msg)
    IF (msglen>0) THEN
    CALL mpi_allreduce(MPI_IN_PLACE,msg,msglen,MPI_DOUBLE_COMPLEX,MPI_SUM,gid,ierr)
    IF ( ierr /= 0 ) CALL mp_stop( ierr, "mpi_allreduce @ "//routineN )
    END IF
    t_end = m_walltime ( )
    CALL add_perf(perf_id=3,count=1,time=t_end-t_start,msg_size=msglen*(2*real_8_size))
#endif
    CALL mp_timestop(handle)
  END SUBROUTINE mp_sum_zv

! *****************************************************************************
!> \brief Element-wise sum of a rank-2 array on all processes.
!> \param[in] msg             Matrix to sum and result
!> \param gid ...
!> \note see mp_sum_z 
! *****************************************************************************
  SUBROUTINE mp_sum_zm(msg,gid)
    COMPLEX(kind=real_8), INTENT(INOUT)                   :: msg( :, : )
    INTEGER, INTENT(IN)                      :: gid

    CHARACTER(len=*), PARAMETER :: routineN = 'mp_sum_zm', &
      routineP = moduleN//':'//routineN

    INTEGER                                  :: handle, ierr
#if defined(__parallel)
    INTEGER, PARAMETER :: max_msg=2**25
    INTEGER                                  :: m1, msglen, step
#endif

    ierr = 0
    CALL mp_timeset(routineN,handle)

#if defined(__parallel)
    t_start = m_walltime ( )
    ! chunk up the call so that message sizes are limited, to avoid overflows in mpich triggered in large rpa calcs
    step=MAX(1,SIZE(msg,2)/MAX(1,SIZE(msg)/max_msg))
    DO m1=LBOUND(msg,2),UBOUND(msg,2), step
       msglen = SIZE(msg,1)*(MIN(UBOUND(msg,2),m1+step-1)-m1+1)
       IF (msglen>0) THEN
          CALL mpi_allreduce(MPI_IN_PLACE,msg(LBOUND(msg,1),m1),msglen,MPI_DOUBLE_COMPLEX,MPI_SUM,gid,ierr)
          IF ( ierr /= 0 ) CALL mp_stop( ierr, "mpi_allreduce @ "//routineN )
       END IF
    ENDDO
    t_end = m_walltime ( )
    CALL add_perf(perf_id=3,count=1,time=t_end-t_start,msg_size=msglen*(2*real_8_size))
#endif
    CALL mp_timestop(handle)
  END SUBROUTINE mp_sum_zm

! *****************************************************************************
!> \brief Element-wise sum of a rank-3 array on all processes.
!> \param[in] msg             Arary to sum and result
!> \param gid ...
!> \note see mp_sum_z 
! *****************************************************************************
  SUBROUTINE mp_sum_zm3(msg,gid)
    COMPLEX(kind=real_8), INTENT(INOUT)                   :: msg( :, :, : )
    INTEGER, INTENT(IN)                      :: gid

    CHARACTER(len=*), PARAMETER :: routineN = 'mp_sum_zm3', &
      routineP = moduleN//':'//routineN

    INTEGER                                  :: handle, ierr, &
                                                msglen
    ierr = 0
    CALL mp_timeset(routineN,handle)

    msglen = SIZE(msg)
#if defined(__parallel)
    t_start = m_walltime ( )
    IF (msglen>0) THEN
      CALL mpi_allreduce(MPI_IN_PLACE,msg,msglen,MPI_DOUBLE_COMPLEX,MPI_SUM,gid,ierr)
      IF ( ierr /= 0 ) CALL mp_stop( ierr, "mpi_allreduce @ "//routineN )
    END IF
    t_end = m_walltime ( )
    CALL add_perf(perf_id=3,count=1,time=t_end-t_start,msg_size=msglen*(2*real_8_size))
#endif
    CALL mp_timestop(handle)
  END SUBROUTINE mp_sum_zm3

! *****************************************************************************
!> \brief Element-wise sum of a rank-4 array on all processes.
!> \param[in] msg             Arary to sum and result
!> \param gid ...
!> \note see mp_sum_z 
! *****************************************************************************
  SUBROUTINE mp_sum_zm4(msg,gid)
    COMPLEX(kind=real_8), INTENT(INOUT)                   :: msg( :, :, :, : )
    INTEGER, INTENT(IN)                      :: gid

    CHARACTER(len=*), PARAMETER :: routineN = 'mp_sum_zm4', &
      routineP = moduleN//':'//routineN

    INTEGER                                  :: handle, ierr, &
                                                msglen

    ierr = 0
    CALL mp_timeset(routineN,handle)

    msglen = SIZE(msg)
#if defined(__parallel)
    t_start = m_walltime ( )
    IF (msglen>0) THEN
      CALL mpi_allreduce(MPI_IN_PLACE,msg,msglen,MPI_DOUBLE_COMPLEX,MPI_SUM,gid,ierr)
      IF ( ierr /= 0 ) CALL mp_stop( ierr, "mpi_allreduce @ "//routineN )
    END IF
    t_end = m_walltime ( )
    CALL add_perf(perf_id=3,count=1,time=t_end-t_start,msg_size=msglen*(2*real_8_size))
#endif
    CALL mp_timestop(handle)
  END SUBROUTINE mp_sum_zm4

! *****************************************************************************
!> \brief Element-wise sum of a rank-5 array on all processes.
!> \param[in] msg             Arary to sum and result
!> \param gid ...
!> \note see mp_sum_z 
! *****************************************************************************
  SUBROUTINE mp_sum_zm5(msg,gid)
    COMPLEX(kind=real_8), INTENT(INOUT)                   :: msg( :, :, :, :, : )
    INTEGER, INTENT(IN)                      :: gid

    CHARACTER(len=*), PARAMETER :: routineN = 'mp_sum_zm5', &
      routineP = moduleN//':'//routineN

    INTEGER                                  :: handle, ierr, &
                                                msglen

    ierr = 0
    CALL mp_timeset(routineN,handle)

    msglen = SIZE(msg)
#if defined(__parallel)
    t_start = m_walltime ( )
    IF (msglen>0) THEN
      CALL mpi_allreduce(MPI_IN_PLACE,msg,msglen,MPI_DOUBLE_COMPLEX,MPI_SUM,gid,ierr)
      IF ( ierr /= 0 ) CALL mp_stop( ierr, "mpi_allreduce @ "//routineN )
    END IF
    t_end = m_walltime ( )
    CALL add_perf(perf_id=3,count=1,time=t_end-t_start,msg_size=msglen*(2*real_8_size))
#endif
    CALL mp_timestop(handle)
  END SUBROUTINE mp_sum_zm5

! *****************************************************************************
!> \brief Element-wise sum of a rank-6 array on all processes.
!> \param[in] msg             Arary to sum and result
!> \param gid ...
!> \note see mp_sum_z 
! *****************************************************************************
  SUBROUTINE mp_sum_zm6(msg,gid)
    COMPLEX(kind=real_8), INTENT(INOUT)                   :: msg( :, :, :, :, :, : )
    INTEGER, INTENT(IN)                      :: gid

    CHARACTER(len=*), PARAMETER :: routineN = 'mp_sum_zm6', &
      routineP = moduleN//':'//routineN

    INTEGER                                  :: handle, ierr, &
                                                msglen

    ierr = 0
    CALL mp_timeset(routineN,handle)

    msglen = SIZE(msg)
#if defined(__parallel)
    t_start = m_walltime ( )
    IF (msglen>0) THEN
      CALL mpi_allreduce(MPI_IN_PLACE,msg,msglen,MPI_DOUBLE_COMPLEX,MPI_SUM,gid,ierr)
      IF ( ierr /= 0 ) CALL mp_stop( ierr, "mpi_allreduce @ "//routineN )
    END IF
    t_end = m_walltime ( )
    CALL add_perf(perf_id=3,count=1,time=t_end-t_start,msg_size=msglen*(2*real_8_size))
#endif
    CALL mp_timestop(handle)
  END SUBROUTINE mp_sum_zm6

! *****************************************************************************
!> \brief Element-wise sum of data from all processes with result left only on
!>        one.
!> \param[in,out] msg         Vector to sum (input) and (only on process root)
!>                            result (output)
!> \param root ...
!> \param[in] gid             Message passing environment identifier
!> \par MPI mapping
!>      mpi_reduce
! *****************************************************************************
  SUBROUTINE mp_sum_root_zv(msg,root,gid)
    COMPLEX(kind=real_8), INTENT(INOUT)                   :: msg( : )
    INTEGER, INTENT(IN)                      :: root, gid

    CHARACTER(len=*), PARAMETER :: routineN = 'mp_sum_root_zv', &
      routineP = moduleN//':'//routineN

    INTEGER                                  :: handle, ierr, msglen
#if defined(__parallel)
    INTEGER                                  :: m1, taskid
    COMPLEX(kind=real_8), ALLOCATABLE                     :: res( : )
#endif

    ierr = 0
    CALL mp_timeset(routineN,handle)

    msglen = SIZE(msg)
#if defined(__parallel)
    t_start = m_walltime ( )
    CALL mpi_comm_rank ( gid, taskid, ierr )
    IF ( ierr /= 0 ) CALL mp_stop( ierr, "mpi_comm_rank @ "//routineN )
    IF (msglen>0) THEN
      m1 = SIZE(msg,1)
      ALLOCATE (res(m1),STAT=ierr)
      IF ( ierr /= 0 ) CALL mp_abort( "allocate @ "//routineN )
      CALL mpi_reduce(msg,res,msglen,MPI_DOUBLE_COMPLEX,MPI_SUM,&
           root,gid,ierr)
      IF ( ierr /= 0 ) CALL mp_stop( ierr, "mpi_reduce @ "//routineN )
      IF ( taskid == root ) THEN
        msg = res
      END IF
      DEALLOCATE (res)
    END IF
    t_end = m_walltime ( )
    CALL add_perf(perf_id=3,count=1,time=t_end-t_start,msg_size=msglen*(2*real_8_size))
#endif
    CALL mp_timestop(handle)
  END SUBROUTINE mp_sum_root_zv

! *****************************************************************************
!> \brief Element-wise sum of data from all processes with result left only on
!>        one.
!> \param[in,out] msg         Matrix to sum (input) and (only on process root)
!>                            result (output)
!> \param root ...
!> \param gid ...
!> \note see mp_sum_root_zv 
! *****************************************************************************
  SUBROUTINE mp_sum_root_zm(msg,root,gid)
    COMPLEX(kind=real_8), INTENT(INOUT)                   :: msg( :, : )
    INTEGER, INTENT(IN)                      :: root, gid

    CHARACTER(len=*), PARAMETER :: routineN = 'mp_sum_root_rm', &
      routineP = moduleN//':'//routineN

    INTEGER                                  :: handle, ierr, msglen
#if defined(__parallel)
    INTEGER                                  :: m1, m2, taskid
    COMPLEX(kind=real_8), ALLOCATABLE                     :: res( :, : )
#endif

    ierr = 0
    CALL mp_timeset(routineN,handle)

    msglen = SIZE(msg)
#if defined(__parallel)
    t_start = m_walltime ( )
    CALL mpi_comm_rank ( gid, taskid, ierr )
    IF ( ierr /= 0 ) CALL mp_stop( ierr, "mpi_comm_rank @ "//routineN )
    IF (msglen>0) THEN
    m1 = SIZE(msg,1)
    m2 = SIZE(msg,2)
    ALLOCATE (res(m1,m2),STAT=ierr)
    IF ( ierr /= 0 ) CALL mp_abort( "allocate @ "//routineN )
    CALL mpi_reduce(msg,res,msglen,MPI_DOUBLE_COMPLEX,MPI_SUM,root,gid,ierr)
    IF ( ierr /= 0 ) CALL mp_stop( ierr, "mpi_reduce @ "//routineN )
    IF ( taskid == root ) THEN
       msg = res
    END IF
    DEALLOCATE (res)
    END IF
    t_end = m_walltime ( )
    CALL add_perf(perf_id=3,count=1,time=t_end-t_start,msg_size=msglen*(2*real_8_size))
#endif
    CALL mp_timestop(handle)
  END SUBROUTINE mp_sum_root_zm

! *****************************************************************************
!> \brief Finds the maximum of a datum with the result left on all processes.
!> \param[in,out] msg         Find maximum among these data (input) and
!>                            maximum (output)
!> \param[in] gid             Message passing environment identifier
!> \par MPI mapping
!>      mpi_allreduce
! *****************************************************************************
  SUBROUTINE mp_max_z(msg,gid)
    COMPLEX(kind=real_8), INTENT(INOUT)                   :: msg
    INTEGER, INTENT(IN)                      :: gid

    CHARACTER(len=*), PARAMETER :: routineN = 'mp_max_z', &
      routineP = moduleN//':'//routineN

    INTEGER                                  :: handle, ierr, msglen

    ierr = 0
    CALL mp_timeset(routineN,handle)

    msglen = 1
#if defined(__parallel)
    t_start = m_walltime ( )
    CALL mpi_allreduce(MPI_IN_PLACE,msg,msglen,MPI_DOUBLE_COMPLEX,MPI_MAX,gid,ierr)
    IF ( ierr /= 0 ) CALL mp_stop( ierr, "mpi_allreduce @ "//routineN )
    t_end = m_walltime ( )
    CALL add_perf(perf_id=3,count=1,time=t_end-t_start,msg_size=msglen*(2*real_8_size))
#endif
    CALL mp_timestop(handle)
  END SUBROUTINE mp_max_z

! *****************************************************************************
!> \brief Finds the element-wise maximum of a vector with the result left on
!>        all processes.
!> \param[in,out] msg         Find maximum among these data (input) and
!>                            maximum (output)
!> \param gid ...
!> \note see mp_max_z 
! *****************************************************************************
  SUBROUTINE mp_max_zv(msg,gid)
    COMPLEX(kind=real_8), INTENT(INOUT)                   :: msg( : )
    INTEGER, INTENT(IN)                      :: gid

    CHARACTER(len=*), PARAMETER :: routineN = 'mp_max_zv', &
      routineP = moduleN//':'//routineN

    INTEGER                                  :: handle, ierr, msglen

    ierr = 0
    CALL mp_timeset(routineN,handle)

    msglen = SIZE(msg)
#if defined(__parallel)
    t_start = m_walltime ( )
    CALL mpi_allreduce(MPI_IN_PLACE,msg,msglen,MPI_DOUBLE_COMPLEX,MPI_MAX,gid,ierr)
    IF ( ierr /= 0 ) CALL mp_stop( ierr, "mpi_allreduce @ "//routineN )
    t_end = m_walltime ( )
    CALL add_perf(perf_id=3,count=1,time=t_end-t_start,msg_size=msglen*(2*real_8_size))
#endif
    CALL mp_timestop(handle)
  END SUBROUTINE mp_max_zv

! *****************************************************************************
!> \brief Finds the minimum of a datum with the result left on all processes.
!> \param[in,out] msg         Find minimum among these data (input) and
!>                            maximum (output)
!> \param[in] gid             Message passing environment identifier
!> \par MPI mapping
!>      mpi_allreduce
! *****************************************************************************
  SUBROUTINE mp_min_z(msg,gid)
    COMPLEX(kind=real_8), INTENT(INOUT)                   :: msg
    INTEGER, INTENT(IN)                      :: gid

    CHARACTER(len=*), PARAMETER :: routineN = 'mp_min_z', &
      routineP = moduleN//':'//routineN

    INTEGER                                  :: handle, ierr, msglen

    ierr = 0
    CALL mp_timeset(routineN,handle)

    msglen = 1
#if defined(__parallel)
    t_start = m_walltime ( )
    CALL mpi_allreduce(MPI_IN_PLACE,msg,msglen,MPI_DOUBLE_COMPLEX,MPI_MIN,gid,ierr)
    IF ( ierr /= 0 ) CALL mp_stop( ierr, "mpi_allreduce @ "//routineN )
    t_end = m_walltime ( )
    CALL add_perf(perf_id=3,count=1,time=t_end-t_start,msg_size=msglen*(2*real_8_size))
#endif
    CALL mp_timestop(handle)
  END SUBROUTINE mp_min_z

! *****************************************************************************
!> \brief Finds the element-wise minimum of vector with the result left on
!>        all processes.
!> \param[in,out] msg         Find minimum among these data (input) and
!>                            maximum (output)
!> \param gid ...
!> \par MPI mapping
!>      mpi_allreduce
!> \note see mp_min_z 
! *****************************************************************************
  SUBROUTINE mp_min_zv(msg,gid)
    COMPLEX(kind=real_8), INTENT(INOUT)                   :: msg( : )
    INTEGER, INTENT(IN)                      :: gid

    CHARACTER(len=*), PARAMETER :: routineN = 'mp_min_zv', &
      routineP = moduleN//':'//routineN

    INTEGER                                  :: handle, ierr, msglen

    ierr = 0
    CALL mp_timeset(routineN,handle)

    msglen = SIZE(msg)
#if defined(__parallel)
    t_start = m_walltime ( )
    IF ( ierr /= 0 ) CALL mp_abort( "allocate @ "//routineN )
    CALL mpi_allreduce(MPI_IN_PLACE,msg,msglen,MPI_DOUBLE_COMPLEX,MPI_MIN,gid,ierr)
    IF ( ierr /= 0 ) CALL mp_stop( ierr, "mpi_allreduce @ "//routineN )
    t_end = m_walltime ( )
    CALL add_perf(perf_id=3,count=1,time=t_end-t_start,msg_size=msglen*(2*real_8_size))
#endif
    CALL mp_timestop(handle)
  END SUBROUTINE mp_min_zv


! *****************************************************************************
!> \brief Scatters data from one processes to all others
!> \param[in] msg_scatter     Data to scatter (for root process)
!> \param[out] msg            Received data
!> \param[in] root            Process which scatters data
!> \param[in] gid             Message passing environment identifier
!> \par MPI mapping
!>      mpi_scatter
! *****************************************************************************
  SUBROUTINE mp_scatter_zv(msg_scatter,msg,root,gid)
    COMPLEX(kind=real_8), INTENT(IN)                      :: msg_scatter(:)
    COMPLEX(kind=real_8), INTENT(OUT)                     :: msg( : )
    INTEGER, INTENT(IN)                      :: root, gid

    CHARACTER(len=*), PARAMETER :: routineN = 'mp_scatter_zv', &
      routineP = moduleN//':'//routineN

    INTEGER                                  :: handle, ierr, msglen

    ierr = 0
    CALL mp_timeset(routineN,handle)

    msglen = SIZE(msg)
#if defined(__parallel)
    t_start = m_walltime ( )
    CALL mpi_scatter(msg_scatter,msglen,MPI_DOUBLE_COMPLEX,msg,&
         msglen,MPI_DOUBLE_COMPLEX,root,gid,ierr)
    IF ( ierr /= 0 ) CALL mp_stop( ierr, "mpi_scatter @ "//routineN )
    t_end = m_walltime ( )
    CALL add_perf(perf_id=4,count=1,time=t_end-t_start,msg_size=msglen*(2*real_8_size))
#else
    msg = msg_scatter
#endif
    CALL mp_timestop(handle)
  END SUBROUTINE mp_scatter_zv

! *****************************************************************************
!> \brief Gathers a datum from all processes to one
!> \param[in] msg             Datum to send to root
!> \param[out] msg_gather     Received data (on root)
!> \param[in] root            Process which gathers the data
!> \param[in] gid             Message passing environment identifier
!> \par MPI mapping
!>      mpi_gather
! *****************************************************************************
  SUBROUTINE mp_gather_z(msg,msg_gather,root,gid)
    COMPLEX(kind=real_8), INTENT(IN)                      :: msg
    COMPLEX(kind=real_8), INTENT(OUT)                     :: msg_gather( : )
    INTEGER, INTENT(IN)                      :: root, gid

    CHARACTER(len=*), PARAMETER :: routineN = 'mp_gather_z', &
      routineP = moduleN//':'//routineN

    INTEGER                                  :: handle, ierr, msglen

    ierr = 0
    CALL mp_timeset(routineN,handle)

    msglen = 1
#if defined(__parallel)
    t_start = m_walltime ( )
    CALL mpi_gather(msg,msglen,MPI_DOUBLE_COMPLEX,msg_gather,&
         msglen,MPI_DOUBLE_COMPLEX,root,gid,ierr)
    IF ( ierr /= 0 ) CALL mp_stop( ierr, "mpi_gather @ "//routineN )
    t_end = m_walltime ( )
    CALL add_perf(perf_id=4,count=1,time=t_end-t_start,msg_size=msglen*(2*real_8_size))
#else
    msg_gather = msg
#endif
    CALL mp_timestop(handle)
  END SUBROUTINE mp_gather_z

! *****************************************************************************
!> \brief Gathers data from all processes to one
!> \param[in] msg             Datum to send to root
!> \param msg_gather ...
!> \param root ...
!> \param gid ...
!> \par Data length
!>      All data (msg) is equal-sized
!> \par MPI mapping
!>      mpi_gather
!> \note see mp_gather_z 
! *****************************************************************************
  SUBROUTINE mp_gather_zv(msg,msg_gather,root,gid)
    COMPLEX(kind=real_8), INTENT(IN)                      :: msg( : )
    COMPLEX(kind=real_8), INTENT(OUT)                     :: msg_gather( : )
    INTEGER, INTENT(IN)                      :: root, gid

    CHARACTER(len=*), PARAMETER :: routineN = 'mp_gather_zv', &
      routineP = moduleN//':'//routineN

    INTEGER                                  :: handle, ierr, msglen

    ierr = 0
    CALL mp_timeset(routineN,handle)

    msglen = SIZE(msg)
#if defined(__parallel)
    t_start = m_walltime ( )
    CALL mpi_gather(msg,msglen,MPI_DOUBLE_COMPLEX,msg_gather,&
         msglen,MPI_DOUBLE_COMPLEX,root,gid,ierr)
    IF ( ierr /= 0 ) CALL mp_stop( ierr, "mpi_gather @ "//routineN )
    t_end = m_walltime ( )
    CALL add_perf(perf_id=4,count=1,time=t_end-t_start,msg_size=msglen*(2*real_8_size))
#else
    msg_gather = msg
#endif
    CALL mp_timestop(handle)
  END SUBROUTINE mp_gather_zv

! *****************************************************************************
!> \brief Gathers data from all processes to one
!> \param[in] msg             Datum to send to root
!> \param msg_gather ...
!> \param root ...
!> \param gid ...
!> \par Data length
!>      All data (msg) is equal-sized
!> \par MPI mapping
!>      mpi_gather
!> \note see mp_gather_z 
! *****************************************************************************
  SUBROUTINE mp_gather_zm(msg,msg_gather,root,gid)
    COMPLEX(kind=real_8), INTENT(IN)                      :: msg( :, : )
    COMPLEX(kind=real_8), INTENT(OUT)                     :: msg_gather( :, : )
    INTEGER, INTENT(IN)                      :: root, gid

    CHARACTER(len=*), PARAMETER :: routineN = 'mp_gather_zm', &
      routineP = moduleN//':'//routineN

    INTEGER                                  :: handle, ierr, msglen

    ierr = 0
    CALL mp_timeset(routineN,handle)

    msglen = SIZE(msg)
#if defined(__parallel)
    t_start = m_walltime ( )
    CALL mpi_gather(msg,msglen,MPI_DOUBLE_COMPLEX,msg_gather,&
         msglen,MPI_DOUBLE_COMPLEX,root,gid,ierr)
    IF ( ierr /= 0 ) CALL mp_stop( ierr, "mpi_gather @ "//routineN )
    t_end = m_walltime ( )
    CALL add_perf(perf_id=4,count=1,time=t_end-t_start,msg_size=msglen*(2*real_8_size))
#else
    msg_gather = msg
#endif
    CALL mp_timestop(handle)
  END SUBROUTINE mp_gather_zm

! *****************************************************************************
!> \brief Gathers data from all processes to one.
!> \param[in] sendbuf         Data to send to root
!> \param[out] recvbuf        Received data (on root)
!> \param[in] recvcounts      Sizes of data received from processes
!> \param[in] displs          Offsets of data received from processes
!> \param[in] root            Process which gathers the data
!> \param[in] comm            Message passing environment identifier
!> \par Data length
!>      Data can have different lengths
!> \par Offsets
!>      Offsets start at 0
!> \par MPI mapping
!>      mpi_gather
! *****************************************************************************
  SUBROUTINE mp_gatherv_zv(sendbuf,recvbuf,recvcounts,displs,root,comm)

    COMPLEX(kind=real_8), DIMENSION(:), INTENT(IN)        :: sendbuf
    COMPLEX(kind=real_8), DIMENSION(:), INTENT(OUT)       :: recvbuf
    INTEGER, DIMENSION(:), INTENT(IN)        :: recvcounts, displs
    INTEGER, INTENT(IN)                      :: root, comm

    CHARACTER(len=*), PARAMETER :: routineN = 'mp_gatherv_zv', &
      routineP = moduleN//':'//routineN

    INTEGER                                  :: handle, ierr
#if defined(__parallel)
    INTEGER                                  :: sendcount
#endif

    ierr = 0
    CALL mp_timeset(routineN,handle)

#if defined(__parallel)
    t_start = m_walltime()
    sendcount = SIZE(sendbuf)
    CALL mpi_gatherv(sendbuf,sendcount,MPI_DOUBLE_COMPLEX,&
         recvbuf,recvcounts,displs,MPI_DOUBLE_COMPLEX,&
         root,comm,ierr)
    IF (ierr /= 0) CALL mp_stop(ierr,"mpi_gatherv @ "//routineN)
    t_end = m_walltime()
    CALL add_perf(perf_id=4,&
         count=1,&
         time=t_end-t_start,&
         msg_size=sendcount*(2*real_8_size))
#else
    recvbuf(1+displs(1):) = sendbuf
#endif
    CALL mp_timestop(handle)
  END SUBROUTINE mp_gatherv_zv


! *****************************************************************************
!> \brief Gathers a datum from all processes and all processes receive the
!>        same data
!> \param[in] msgout          Datum to send
!> \param[out] msgin          Received data
!> \param[in] gid             Message passing environment identifier
!> \par Data size
!>      All processes send equal-sized data
!> \par MPI mapping
!>      mpi_allgather
! *****************************************************************************
  SUBROUTINE mp_allgather_z(msgout,msgin,gid)
    COMPLEX(kind=real_8), INTENT(IN)                      :: msgout
    COMPLEX(kind=real_8), INTENT(OUT)                     :: msgin( : )
    INTEGER, INTENT(IN)                      :: gid

    CHARACTER(len=*), PARAMETER :: routineN = 'mp_allgather_z', &
      routineP = moduleN//':'//routineN

    INTEGER                                  :: handle, ierr
#if defined(__parallel)
    INTEGER                                  :: rcount, scount
#endif

    ierr = 0
    CALL mp_timeset(routineN,handle)

#if defined(__parallel)
    scount = 1
    rcount = 1
    CALL MPI_ALLGATHER(msgout, scount, MPI_DOUBLE_COMPLEX, &
                       msgin , rcount, MPI_DOUBLE_COMPLEX, &
                       gid, ierr )
    IF ( ierr /= 0 ) CALL mp_stop( ierr, "mpi_allgather @ "//routineN )
#else
    msgin = msgout
#endif
    CALL mp_timestop(handle)
  END SUBROUTINE mp_allgather_z

! *****************************************************************************
!> \brief Gathers vector data from all processes and all processes receive the
!>        same data
!> \param[in] msgout          Rank-1 data to send
!> \param[out] msgin          Received data
!> \param[in] gid             Message passing environment identifier
!> \par Data size
!>      All processes send equal-sized data
!> \par Ranks
!>      The last rank counts the processes
!> \par MPI mapping
!>      mpi_allgather
! *****************************************************************************
  SUBROUTINE mp_allgather_z12(msgout, msgin,gid)
    COMPLEX(kind=real_8), INTENT(IN)                      :: msgout(:)
    COMPLEX(kind=real_8), INTENT(OUT)                     :: msgin(:, :)
    INTEGER, INTENT(IN)                      :: gid

    CHARACTER(len=*), PARAMETER :: routineN = 'mp_allgather_z12', &
      routineP = moduleN//':'//routineN

    INTEGER                                  :: handle, ierr
#if defined(__parallel)
    INTEGER                                  :: rcount, scount
#endif

    ierr = 0
    CALL mp_timeset(routineN,handle)

#if defined(__parallel)
    scount = SIZE (msgout(:))
    rcount = scount
    CALL MPI_ALLGATHER(msgout, scount, MPI_DOUBLE_COMPLEX, &
                       msgin , rcount, MPI_DOUBLE_COMPLEX, &
                       gid, ierr )
    IF ( ierr /= 0 ) CALL mp_stop( ierr, "mpi_allgather @ "//routineN )
#else
    msgin(:,1) = msgout(:)
#endif
    CALL mp_timestop(handle)
  END SUBROUTINE mp_allgather_z12

! *****************************************************************************
!> \brief Gathers matrix data from all processes and all processes receive the
!>        same data
!> \param[in] msgout          Rank-2 data to send
!> \param msgin ...
!> \param gid ...
!> \note see mp_allgather_z12 
! *****************************************************************************
  SUBROUTINE mp_allgather_z23(msgout, msgin,gid)
    COMPLEX(kind=real_8), INTENT(IN)                      :: msgout(:,:)
    COMPLEX(kind=real_8), INTENT(OUT)                     :: msgin(:, :, :)
    INTEGER, INTENT(IN)                      :: gid

    CHARACTER(len=*), PARAMETER :: routineN = 'mp_allgather_z23', &
      routineP = moduleN//':'//routineN

    INTEGER                                  :: handle, ierr
#if defined(__parallel)
    INTEGER                                  :: rcount, scount
#endif

    ierr = 0
    CALL mp_timeset(routineN,handle)

#if defined(__parallel)
    scount = SIZE (msgout(:,:))
    rcount = scount
    CALL MPI_ALLGATHER(msgout, scount, MPI_DOUBLE_COMPLEX, &
                       msgin , rcount, MPI_DOUBLE_COMPLEX, &
                       gid, ierr )
    IF ( ierr /= 0 ) CALL mp_stop( ierr, "mpi_allgather @ "//routineN )
#else
    msgin(:,:,1) = msgout(:,:)
#endif
    CALL mp_timestop(handle)
  END SUBROUTINE mp_allgather_z23

! *****************************************************************************
!> \brief Gathers rank-3 data from all processes and all processes receive the
!>        same data
!> \param[in] msgout          Rank-3 data to send
!> \param msgin ...
!> \param gid ...
!> \note see mp_allgather_z12 
! *****************************************************************************
  SUBROUTINE mp_allgather_z34(msgout, msgin,gid)
    COMPLEX(kind=real_8), INTENT(IN)                      :: msgout(:,:, :)
    COMPLEX(kind=real_8), INTENT(OUT)                     :: msgin(:, :, :, :)
    INTEGER, INTENT(IN)                      :: gid

    CHARACTER(len=*), PARAMETER :: routineN = 'mp_allgather_z34', &
      routineP = moduleN//':'//routineN

    INTEGER                                  :: handle, ierr
#if defined(__parallel)
    INTEGER                                  :: rcount, scount
#endif

    ierr = 0
    CALL mp_timeset(routineN,handle)

#if defined(__parallel)
    scount = SIZE (msgout(:,:,:))
    rcount = scount
    CALL MPI_ALLGATHER(msgout, scount, MPI_DOUBLE_COMPLEX, &
                       msgin , rcount, MPI_DOUBLE_COMPLEX, &
                       gid, ierr )
    IF ( ierr /= 0 ) CALL mp_stop( ierr, "mpi_allgather @ "//routineN )
#else
    msgin(:,:,:,1) = msgout(:,:,:)
#endif
    CALL mp_timestop(handle)
  END SUBROUTINE mp_allgather_z34

! *****************************************************************************
!> \brief Gathers vector data from all processes and all processes receive the
!>        same data
!> \param[in] msgout          Rank-1 data to send
!> \param[out] msgin          Received data
!> \param[in] rcount          Size of sent data for every process
!> \param[in] rdispl          Offset of sent data for every process
!> \param[in] gid             Message passing environment identifier
!> \par Data size
!>      Processes can send different-sized data
!> \par Ranks
!>      The last rank counts the processes
!> \par Offsets
!>      Offsets are from 0
!> \par MPI mapping
!>      mpi_allgather
! *****************************************************************************
  SUBROUTINE mp_allgatherv_zv(msgout,msgin,rcount,rdispl,gid)
    COMPLEX(kind=real_8), INTENT(IN)                      :: msgout( : )
    COMPLEX(kind=real_8), INTENT(OUT)                     :: msgin( : )
    INTEGER, INTENT(IN)                      :: rcount( : ), rdispl( : ), gid

    CHARACTER(len=*), PARAMETER :: routineN = 'mp_allgatherv_zv', &
      routineP = moduleN//':'//routineN

    INTEGER                                  :: handle, ierr
#if defined(__parallel)
    INTEGER                                  :: scount
#endif

    ierr = 0
    CALL mp_timeset(routineN,handle)

#if defined(__parallel)
    scount = SIZE ( msgout )
    CALL MPI_ALLGATHERV(msgout, scount, MPI_DOUBLE_COMPLEX, msgin, rcount, &
                        rdispl, MPI_DOUBLE_COMPLEX, gid, ierr )
    IF ( ierr /= 0 ) CALL mp_stop( ierr, "mpi_allgatherv @ "//routineN )
#else
    msgin = msgout
#endif
    CALL mp_timestop(handle)
  END SUBROUTINE mp_allgatherv_zv

! *****************************************************************************
!> \brief Sums a vector and partitions the result among processes
!> \param[in] msgout          Data to sum
!> \param[out] msgin          Received portion of summed data
!> \param[in] rcount          Partition sizes of the summed data for
!>                            every process
!> \param[in] gid             Message passing environment identifier
! *****************************************************************************
  SUBROUTINE mp_sum_scatter_zv(msgout,msgin,rcount,gid)
    COMPLEX(kind=real_8), INTENT(IN)                      :: msgout( : )
    COMPLEX(kind=real_8), INTENT(OUT)                     :: msgin( : )
    INTEGER, INTENT(IN)                      :: rcount( : ), gid

    CHARACTER(len=*), PARAMETER :: routineN = 'mp_sum_scatter_zv', &
      routineP = moduleN//':'//routineN

    INTEGER                                  :: handle, ierr

    ierr = 0
    CALL mp_timeset(routineN,handle)

#if defined(__parallel)
    t_start = m_walltime ( )
    CALL MPI_REDUCE_SCATTER(msgout, msgin, rcount, MPI_DOUBLE_COMPLEX, MPI_SUM, &
         gid, ierr )
    IF ( ierr /= 0 ) CALL mp_stop( ierr, "mpi_reduce_scatter @ "//routineN )

    t_end = m_walltime ( )
    CALL add_perf(perf_id=3,count=1,time=t_end-t_start,&
         msg_size=rcount(1)*2*(2*real_8_size))
#else
    msgin = msgout
#endif
    CALL mp_timestop(handle)
  END SUBROUTINE mp_sum_scatter_zv

! *****************************************************************************
!> \brief Sends and receives vector data
!> \param[in] msgin           Data to send
!> \param[in] dest            Process to send data to
!> \param[out] msgout         Received data
!> \param[in] source          Process from which to receive
!> \param[in] comm            Message passing environment identifier
! *****************************************************************************
  SUBROUTINE mp_sendrecv_zv(msgin,dest,msgout,source,comm)
    COMPLEX(kind=real_8), INTENT(IN)                      :: msgin( : )
    INTEGER, INTENT(IN)                      :: dest
    COMPLEX(kind=real_8), INTENT(OUT)                     :: msgout( : )
    INTEGER, INTENT(IN)                      :: source, comm

    CHARACTER(len=*), PARAMETER :: routineN = 'mp_sendrecv_zv', &
      routineP = moduleN//':'//routineN

    INTEGER                                  :: handle, ierr
#if defined(__parallel)
    INTEGER                                  :: msglen_in, msglen_out, &
                                                recv_tag, send_tag
    INTEGER, ALLOCATABLE, DIMENSION(:)       :: status
#endif

    ierr = 0
    CALL mp_timeset(routineN,handle)

#if defined(__parallel)
    ALLOCATE(status(MPI_STATUS_SIZE))
    t_start = m_walltime ( )
    msglen_in = SIZE(msgin)
    msglen_out = SIZE(msgout)
    send_tag = 0 ! cannot think of something better here, this might be dangerous
    recv_tag = 0 ! cannot think of something better here, this might be dangerous
    CALL mpi_sendrecv(msgin,msglen_in,MPI_DOUBLE_COMPLEX,dest,send_tag,msgout,&
         msglen_out,MPI_DOUBLE_COMPLEX,source,recv_tag,comm,status(1),ierr)
    ! we do not check the status
    IF ( ierr /= 0 ) CALL mp_stop( ierr, "mpi_sendrecv @ "//routineN )
    t_end = m_walltime ( )
    CALL add_perf(perf_id=7,count=1,time=t_end-t_start,&
         msg_size=(msglen_in+msglen_out)*(2*real_8_size)/2)
    DEALLOCATE(status)
#else
    msgout = msgin
#endif
    CALL mp_timestop(handle)
  END SUBROUTINE mp_sendrecv_zv

! *****************************************************************************
!> \brief Sends and receives matrix data
!> \param msgin ...
!> \param dest ...
!> \param msgout ...
!> \param source ...
!> \param comm ...
!> \note see mp_sendrecv_zv 
! *****************************************************************************
  SUBROUTINE mp_sendrecv_zm2(msgin,dest,msgout,source,comm)
    COMPLEX(kind=real_8), INTENT(IN)                      :: msgin( :, : )
    INTEGER, INTENT(IN)                      :: dest
    COMPLEX(kind=real_8), INTENT(OUT)                     :: msgout( :, : )
    INTEGER, INTENT(IN)                      :: source, comm

    CHARACTER(len=*), PARAMETER :: routineN = 'mp_sendrecv_zm2', &
      routineP = moduleN//':'//routineN

    INTEGER                                  :: handle, ierr
#if defined(__parallel)
    INTEGER                                  :: msglen_in, msglen_out, &
                                                recv_tag, send_tag
    INTEGER, ALLOCATABLE, DIMENSION(:)       :: status
#endif

    ierr = 0
    CALL mp_timeset(routineN,handle)

#if defined(__parallel)
    ALLOCATE(status(MPI_STATUS_SIZE))
    t_start = m_walltime ( )
    msglen_in = SIZE(msgin,1)*SIZE(msgin,2)
    msglen_out = SIZE(msgout,1)*SIZE(msgout,2)
    send_tag = 0 ! cannot think of something better here, this might be dangerous
    recv_tag = 0 ! cannot think of something better here, this might be dangerous
    CALL mpi_sendrecv(msgin,msglen_in,MPI_DOUBLE_COMPLEX,dest,send_tag,msgout,&
         msglen_out,MPI_DOUBLE_COMPLEX,source,recv_tag,comm,status(1),ierr)
    ! we do not check the status
    IF ( ierr /= 0 ) CALL mp_stop( ierr, "mpi_sendrecv @ "//routineN )
    t_end = m_walltime ( )
    CALL add_perf(perf_id=7,count=1,time=t_end-t_start,&
         msg_size=(msglen_in+msglen_out)*(2*real_8_size)/2)
    DEALLOCATE(status)
#else
    msgout = msgin
#endif
    CALL mp_timestop(handle)
  END SUBROUTINE mp_sendrecv_zm2

! *****************************************************************************
!> \brief Sends and receives rank-3 data
!> \param msgin ...
!> \param dest ...
!> \param msgout ...
!> \param source ...
!> \param comm ...
!> \note see mp_sendrecv_zv 
! *****************************************************************************
  SUBROUTINE mp_sendrecv_zm3(msgin,dest,msgout,source,comm)
    COMPLEX(kind=real_8), INTENT(IN)                      :: msgin( :, :, : )
    INTEGER, INTENT(IN)                      :: dest
    COMPLEX(kind=real_8), INTENT(OUT)                     :: msgout( :, :, : )
    INTEGER, INTENT(IN)                      :: source, comm

    CHARACTER(len=*), PARAMETER :: routineN = 'mp_sendrecv_zm3', &
      routineP = moduleN//':'//routineN

    INTEGER                                  :: handle, ierr
#if defined(__parallel)
    INTEGER                                  :: msglen_in, msglen_out, &
                                                recv_tag, send_tag
    INTEGER, ALLOCATABLE, DIMENSION(:)       :: status
#endif

    ierr = 0
    CALL mp_timeset(routineN,handle)

#if defined(__parallel)
    ALLOCATE(status(MPI_STATUS_SIZE))
    t_start = m_walltime ( )
    msglen_in = SIZE(msgin)
    msglen_out = SIZE(msgout)
    send_tag = 0 ! cannot think of something better here, this might be dangerous
    recv_tag = 0 ! cannot think of something better here, this might be dangerous
    CALL mpi_sendrecv(msgin,msglen_in,MPI_DOUBLE_COMPLEX,dest,send_tag,msgout,&
         msglen_out,MPI_DOUBLE_COMPLEX,source,recv_tag,comm,status(1),ierr)
    ! we do not check the status
    IF ( ierr /= 0 ) CALL mp_stop( ierr, "mpi_sendrecv @ "//routineN )
    t_end = m_walltime ( )
    CALL add_perf(perf_id=7,count=1,time=t_end-t_start,&
         msg_size=(msglen_in+msglen_out)*(2*real_8_size)/2)
    DEALLOCATE(status)
#else
    msgout = msgin
#endif
    CALL mp_timestop(handle)
  END SUBROUTINE mp_sendrecv_zm3

! *****************************************************************************
!> \brief Non-blocking send and receieve of a vector
!> \param[in] msgin           Vector data to send
!> \param[in] dest            Which process to send to
!> \param[out] msgout         Receive data into this pointer
!> \param[in] source          Process to receive from
!> \param[in] comm            Message passing environment identifier
!> \param[out] send_request   Request handle for the send
!> \param[out] recv_request   Request handle for the receive
!> \param[in] tag             (optional) tag to differentiate requests
!> \par Implementation
!>      Calls mpi_isend and mpi_irecv.
!> \par History
!>      11.2004 created [Joost VandeVondele]
!> \note
!>      The arguments must be pointers to be sure that we do not get
!>      temporaries. They must point to contiguous memory.
! *****************************************************************************
  SUBROUTINE mp_isendrecv_zv(msgin,dest,msgout,source,comm,send_request,&
       recv_request,tag)
    COMPLEX(kind=real_8), DIMENSION(:), POINTER           :: msgin
    INTEGER, INTENT(IN)                      :: dest
    COMPLEX(kind=real_8), DIMENSION(:), POINTER           :: msgout
    INTEGER, INTENT(IN)                      :: source, comm
    INTEGER, INTENT(out)                     :: send_request, recv_request
    INTEGER, INTENT(in), OPTIONAL            :: tag

    CHARACTER(len=*), PARAMETER :: routineN = 'mp_isendrecv_zv', &
      routineP = moduleN//':'//routineN

    INTEGER                                  :: handle, ierr
#if defined(__parallel)
    INTEGER                                  :: lower1, msglen, my_tag
    COMPLEX(kind=real_8)                                  :: foo
#endif

    ierr = 0
    CALL mp_timeset(routineN,handle)

#if defined(__parallel)
    t_start = m_walltime ( )
    my_tag = 0
    IF (PRESENT(tag)) my_tag=tag

    msglen = SIZE(msgout,1)
    IF (msglen>0) THEN
       lower1=LBOUND(msgout,1)
       CALL mpi_irecv(msgout(lower1),msglen,MPI_DOUBLE_COMPLEX,source, my_tag,&
            comm,recv_request,ierr)
    ELSE
       CALL mpi_irecv(foo,msglen,MPI_DOUBLE_COMPLEX,source, my_tag,&
            comm,recv_request,ierr)
    END IF
    IF ( ierr /= 0 ) CALL mp_stop( ierr, "mpi_irecv @ "//routineN )

    msglen = SIZE(msgin,1)
    IF (msglen>0) THEN
       lower1=LBOUND(msgin,1)
       CALL mpi_isend(msgin(lower1),msglen,MPI_DOUBLE_COMPLEX,dest,my_tag,&
            comm,send_request,ierr)
    ELSE
       CALL mpi_isend(foo,msglen,MPI_DOUBLE_COMPLEX,dest,my_tag,&
            comm,send_request,ierr)
    END IF
    IF ( ierr /= 0 ) CALL mp_stop( ierr, "mpi_isend @ "//routineN )

    msglen = (msglen+SIZE(msgout,1))/2.0_dp
    t_end = m_walltime ( )
    CALL add_perf(perf_id=8,count=1,time=t_end-t_start,msg_size=msglen*(2*real_8_size))
#else
    send_request=0
    recv_request=0
    msgout = msgin
#endif
    CALL mp_timestop(handle)
  END SUBROUTINE mp_isendrecv_zv

! *****************************************************************************
!> \brief Non-blocking send and receieve of a matrix
!> \param msgin ...
!> \param dest ...
!> \param msgout ...
!> \param source ...
!> \param comm ...
!> \param send_request ...
!> \param recv_request ...
!> \param tag ...
!> \par History
!>      08.2003 created [f&j]
!> \note see mp_isendrecv_zv 
!> \note
!>      The arguments must be pointers to be sure that we do not get
!>      temporaries. They must point to contiguous memory.
! *****************************************************************************
  SUBROUTINE mp_isendrecv_zm2(msgin,dest,msgout,source,comm,send_request,&
       recv_request,tag)
    COMPLEX(kind=real_8), DIMENSION(:, :), POINTER        :: msgin
    INTEGER, INTENT(IN)                      :: dest
    COMPLEX(kind=real_8), DIMENSION(:, :), POINTER        :: msgout
    INTEGER, INTENT(IN)                      :: source, comm
    INTEGER, INTENT(out)                     :: send_request, recv_request
    INTEGER, INTENT(in), OPTIONAL            :: tag

    CHARACTER(len=*), PARAMETER :: routineN = 'mp_isendrecv_zm2', &
      routineP = moduleN//':'//routineN

    INTEGER                                  :: handle, ierr
#if defined(__parallel)
    INTEGER                                  :: lower1, lower2, msglen, my_tag
    COMPLEX(kind=real_8)                                  :: foo
#endif

    ierr = 0
    CALL mp_timeset(routineN,handle)

#if defined(__parallel)
    t_start = m_walltime ( )
    my_tag = 0
    IF (PRESENT(tag)) my_tag=tag

    msglen = SIZE(msgout,1)*SIZE(msgout,2)
    IF (msglen>0) THEN
       lower1=LBOUND(msgout,1)
       lower2=LBOUND(msgout,2)
       CALL mpi_irecv(msgout(lower1,lower2),msglen,MPI_DOUBLE_COMPLEX,source, my_tag,&
            comm,recv_request,ierr)
    ELSE
       CALL mpi_irecv(foo,msglen,MPI_DOUBLE_COMPLEX,source, my_tag,&
            comm,recv_request,ierr)
    END IF
    IF ( ierr /= 0 ) CALL mp_stop( ierr, "mpi_irecv @ "//routineN )

    msglen = SIZE(msgin,1)*SIZE(msgin,2)
    IF (msglen>0) THEN
       lower1=LBOUND(msgin,1)
       lower2=LBOUND(msgin,2)
       CALL mpi_isend(msgin(lower1,lower2),msglen,MPI_DOUBLE_COMPLEX,dest,my_tag,&
            comm,send_request,ierr)
    ELSE
       CALL mpi_isend(foo,msglen,MPI_DOUBLE_COMPLEX,dest,my_tag,&
            comm,send_request,ierr)
    END IF
    IF ( ierr /= 0 ) CALL mp_stop( ierr, "mpi_isend @ mp_isendrecv_rm2" )

    msglen = (msglen+SIZE(msgout,1)*SIZE(msgout,2))/2.0_dp
    t_end = m_walltime ( )
    CALL add_perf(perf_id=8,count=1,time=t_end-t_start,msg_size=msglen*(2*real_8_size))
#else
    send_request=0
    recv_request=0
    msgout = msgin
#endif
    CALL mp_timestop(handle)
  END SUBROUTINE mp_isendrecv_zm2

! *****************************************************************************
!> \brief Non-blocking send of vector data
!> \param msgin ...
!> \param dest ...
!> \param comm ...
!> \param request ...
!> \param tag ...
!> \par History
!>      08.2003 created [f&j]
!> \note see mp_isendrecv_zv 
!> \note
!>      The argument must be a pointer to be sure that we do not get
!>      temporaries. They must point to contiguous memory.
! *****************************************************************************
  SUBROUTINE mp_isend_zv(msgin,dest,comm,request,tag)
    COMPLEX(kind=real_8), DIMENSION(:), POINTER           :: msgin
    INTEGER, INTENT(IN)                      :: dest, comm
    INTEGER, INTENT(out)                     :: request
    INTEGER, INTENT(in), OPTIONAL            :: tag

    CHARACTER(len=*), PARAMETER :: routineN = 'mp_isend_zv', &
      routineP = moduleN//':'//routineN

    INTEGER                                  :: handle, ierr
#if defined(__parallel)
    INTEGER                                  :: lower1, msglen, my_tag
    COMPLEX(kind=real_8)                                  :: foo(1)
#endif

    ierr = 0
    CALL mp_timeset(routineN,handle)

#if defined(__parallel)
    t_start = m_walltime ( )
    my_tag = 0
    IF (PRESENT(tag)) my_tag=tag

    msglen = SIZE(msgin)
    IF (msglen>0) THEN
       lower1=LBOUND(msgin,1)
       CALL mpi_isend(msgin(lower1),msglen,MPI_DOUBLE_COMPLEX,dest,my_tag,&
            comm,request,ierr)
    ELSE
       CALL mpi_isend(foo,msglen,MPI_DOUBLE_COMPLEX,dest,my_tag,&
            comm,request,ierr)
    END IF
    IF ( ierr /= 0 ) CALL mp_stop( ierr, "mpi_isend @ "//routineN )

    t_end = m_walltime ( )
    CALL add_perf(perf_id=11,count=1,time=t_end-t_start,msg_size=2*msglen*(2*real_8_size))
#else
    ierr=1
    CALL mp_stop( ierr, "mp_isend called in non parallel case" )
#endif
    CALL mp_timestop(handle)
  END SUBROUTINE mp_isend_zv

! *****************************************************************************
!> \brief Non-blocking send of matrix data
!> \param msgin ...
!> \param dest ...
!> \param comm ...
!> \param request ...
!> \param tag ...
!> \par History
!>      2009-11-25 [UB] Made type-generic for templates
!> \author fawzi
!> \note see mp_isendrecv_zv 
!> \note see mp_isend_zv 
!> \note
!>      The argument must be a pointer to be sure that we do not get
!>      temporaries. They must point to contiguous memory.
! *****************************************************************************
  SUBROUTINE mp_isend_zm2(msgin,dest,comm,request,tag)
    COMPLEX(kind=real_8), DIMENSION(:, :), POINTER  :: msgin
    INTEGER, INTENT(IN)                      :: dest, comm
    INTEGER, INTENT(out)                     :: request
    INTEGER, INTENT(in), OPTIONAL            :: tag

    CHARACTER(len=*), PARAMETER :: routineN = 'mp_isend_zm2', &
      routineP = moduleN//':'//routineN

    INTEGER                                  :: handle, ierr
#if defined(__parallel)
    INTEGER                                  :: lower1, lower2, msglen, my_tag
    COMPLEX(kind=real_8)                                  :: foo(1)
#endif

    ierr = 0
    CALL mp_timeset(routineN,handle)

#if defined(__parallel)
    t_start = m_walltime ( )
    my_tag = 0
    IF (PRESENT(tag)) my_tag=tag

    msglen = SIZE(msgin,1)*SIZE(msgin,2)
    IF (msglen>0) THEN
       lower1=LBOUND(msgin,1)
       lower2=LBOUND(msgin,2)
       CALL mpi_isend(msgin(lower1,lower2),msglen,MPI_DOUBLE_COMPLEX,dest,my_tag,&
            comm,request,ierr)
    ELSE
       CALL mpi_isend(foo,msglen,MPI_DOUBLE_COMPLEX,dest,my_tag,&
            comm,request,ierr)
    END IF
    IF ( ierr /= 0 ) CALL mp_stop( ierr, "mpi_isend @ "//routineN )

    t_end = m_walltime ( )
    CALL add_perf(perf_id=11,count=1,time=t_end-t_start,msg_size=msglen*(2*real_8_size))
#else
    ierr=1
    CALL mp_stop( ierr, "mp_isend called in non parallel case" )
#endif
    CALL mp_timestop(handle)
  END SUBROUTINE mp_isend_zm2

! *****************************************************************************
!> \brief Non-blocking send of rank-3 data
!> \param msgin ...
!> \param dest ...
!> \param comm ...
!> \param request ...
!> \param tag ...
!> \par History
!>      9.2008 added _rm3 subroutine [Iain Bethune]
!>     (c) The Numerical Algorithms Group (NAG) Ltd, 2008 on behalf of the HECToR project
!>      2009-11-25 [UB] Made type-generic for templates
!> \author fawzi
!> \note see mp_isendrecv_zv 
!> \note see mp_isend_zv 
!> \note
!>      The argument must be a pointer to be sure that we do not get
!>      temporaries. They must point to contiguous memory.
! *****************************************************************************
  SUBROUTINE mp_isend_zm3(msgin,dest,comm,request,tag)
    COMPLEX(kind=real_8), DIMENSION(:, :, :), &
      POINTER                                :: msgin
    INTEGER, INTENT(IN)                      :: dest, comm
    INTEGER, INTENT(out)                     :: request
    INTEGER, INTENT(in), OPTIONAL            :: tag

    CHARACTER(len=*), PARAMETER :: routineN = 'mp_isend_zm3', &
      routineP = moduleN//':'//routineN

    INTEGER                                  :: handle, ierr
#if defined(__parallel)
    INTEGER                                  :: lower1, lower2, lower3, &
                                                msglen, my_tag
    COMPLEX(kind=real_8)                                  :: foo(1)
#endif

    ierr = 0
    CALL mp_timeset(routineN,handle)

#if defined(__parallel)
    t_start = m_walltime ( )
    my_tag = 0
    IF (PRESENT(tag)) my_tag=tag

    msglen = SIZE(msgin,1)*SIZE(msgin,2)*SIZE(msgin,3)
    IF (msglen>0) THEN
       lower1=LBOUND(msgin,1)
       lower2=LBOUND(msgin,2)
       lower3=LBOUND(msgin,3)
       CALL mpi_isend(msgin(lower1,lower2,lower3),msglen,MPI_DOUBLE_COMPLEX,dest,my_tag,&
            comm,request,ierr)
    ELSE
       CALL mpi_isend(foo,msglen,MPI_DOUBLE_COMPLEX,dest,my_tag,&
            comm,request,ierr)
    END IF
    IF ( ierr /= 0 ) CALL mp_stop( ierr, "mpi_isend @ "//routineN )

    t_end = m_walltime ( )
    CALL add_perf(perf_id=11,count=1,time=t_end-t_start,msg_size=msglen*(2*real_8_size))
#else
    ierr=1
    CALL mp_stop( ierr, "mp_isend called in non parallel case" )
#endif
    CALL mp_timestop(handle)
  END SUBROUTINE mp_isend_zm3

! *****************************************************************************
!> \brief Non-blocking receive of vector data
!> \param msgout ...
!> \param source ...
!> \param comm ...
!> \param request ...
!> \param tag ...
!> \par History
!>      08.2003 created [f&j]
!>      2009-11-25 [UB] Made type-generic for templates
!> \note see mp_isendrecv_zv 
!> \note
!>      The argument must be a pointer to be sure that we do not get
!>      temporaries. They must point to contiguous memory.
! *****************************************************************************
  SUBROUTINE mp_irecv_zv(msgout,source,comm,request,tag)
    COMPLEX(kind=real_8), DIMENSION(:), POINTER           :: msgout
    INTEGER, INTENT(IN)                      :: source, comm
    INTEGER, INTENT(out)                     :: request
    INTEGER, INTENT(in), OPTIONAL            :: tag

    CHARACTER(len=*), PARAMETER :: routineN = 'mp_irecv_zv', &
      routineP = moduleN//':'//routineN

    INTEGER                                  :: handle, ierr
#if defined(__parallel)
    INTEGER                                  :: lower1, msglen, my_tag
    COMPLEX(kind=real_8)                                  :: foo(1)
#endif

    ierr = 0
    CALL mp_timeset(routineN,handle)

#if defined(__parallel)
    t_start = m_walltime ( )
    my_tag = 0
    IF (PRESENT(tag)) my_tag=tag

    msglen = SIZE(msgout)
    IF (msglen>0) THEN
       lower1=LBOUND(msgout,1)
       CALL mpi_irecv(msgout(lower1),msglen,MPI_DOUBLE_COMPLEX,source, my_tag,&
            comm,request,ierr)
    ELSE
       CALL mpi_irecv(foo,msglen,MPI_DOUBLE_COMPLEX,source, my_tag,&
            comm,request,ierr)
    END IF
    IF ( ierr /= 0 ) CALL mp_stop( ierr, "mpi_irecv @ "//routineN )

    t_end = m_walltime ( )
    CALL add_perf(perf_id=12,count=1,time=t_end-t_start,msg_size=2*msglen*(2*real_8_size))
#else
    CALL mp_abort( "mp_irecv called in non parallel case" )
#endif
    CALL mp_timestop(handle)
  END SUBROUTINE mp_irecv_zv

! *****************************************************************************
!> \brief Non-blocking receive of matrix data
!> \param msgout ...
!> \param source ...
!> \param comm ...
!> \param request ...
!> \param tag ...
!> \par History
!>      2009-11-25 [UB] Made type-generic for templates
!> \author fawzi
!> \note see mp_isendrecv_zv 
!> \note see mp_irecv_zv 
!> \note
!>      The argument must be a pointer to be sure that we do not get
!>      temporaries. They must point to contiguous memory.
! *****************************************************************************
  SUBROUTINE mp_irecv_zm2(msgout,source,comm,request,tag)
    COMPLEX(kind=real_8), DIMENSION(:, :), POINTER        :: msgout
    INTEGER, INTENT(IN)                      :: source, comm
    INTEGER, INTENT(out)                     :: request
    INTEGER, INTENT(in), OPTIONAL            :: tag

    CHARACTER(len=*), PARAMETER :: routineN = 'mp_irecv_zm2', &
      routineP = moduleN//':'//routineN

    INTEGER                                  :: handle, ierr
#if defined(__parallel)
    INTEGER                                  :: lower1, lower2, msglen, my_tag
    COMPLEX(kind=real_8)                                  :: foo(1)
#endif

    ierr = 0
    CALL mp_timeset(routineN,handle)

#if defined(__parallel)
    t_start = m_walltime ( )
    my_tag = 0
    IF (PRESENT(tag)) my_tag=tag

    msglen = SIZE(msgout,1)*SIZE(msgout,2)
    IF (msglen>0) THEN
       lower1=LBOUND(msgout,1)
       lower2=LBOUND(msgout,2)
       CALL mpi_irecv(msgout(lower1,lower2),msglen,MPI_DOUBLE_COMPLEX,source, my_tag,&
            comm,request,ierr)
    ELSE
       CALL mpi_irecv(foo,msglen,MPI_DOUBLE_COMPLEX,source, my_tag,&
            comm,request,ierr)
    END IF
    IF ( ierr /= 0 ) CALL mp_stop( ierr, "mpi_irecv @ "//routineN )

    t_end = m_walltime ( )
    CALL add_perf(perf_id=12,count=1,time=t_end-t_start,msg_size=msglen*(2*real_8_size))
#else
    CALL mp_abort( "mp_irecv called in non parallel case" )
#endif
    CALL mp_timestop(handle)
  END SUBROUTINE mp_irecv_zm2


! *****************************************************************************
!> \brief Non-blocking send of rank-3 data
!> \param msgout ...
!> \param source ...
!> \param comm ...
!> \param request ...
!> \param tag ...
!> \par History
!>      9.2008 added _rm3 subroutine [Iain Bethune] (c) The Numerical Algorithms Group (NAG) Ltd, 2008 on behalf of the HECToR project
!>      2009-11-25 [UB] Made type-generic for templates
!> \author fawzi
!> \note see mp_isendrecv_zv 
!> \note see mp_irecv_zv 
!> \note
!>      The argument must be a pointer to be sure that we do not get
!>      temporaries. They must point to contiguous memory.
! *****************************************************************************
  SUBROUTINE mp_irecv_zm3(msgout,source,comm,request,tag)
    COMPLEX(kind=real_8), DIMENSION(:, :, :), &
      POINTER                                :: msgout
    INTEGER, INTENT(IN)                      :: source, comm
    INTEGER, INTENT(out)                     :: request
    INTEGER, INTENT(in), OPTIONAL            :: tag

    CHARACTER(len=*), PARAMETER :: routineN = 'mp_irecv_zm3', &
      routineP = moduleN//':'//routineN

    INTEGER                                  :: handle, ierr
#if defined(__parallel)
    INTEGER                                  :: lower1, lower2, lower3, &
                                                msglen, my_tag
    COMPLEX(kind=real_8)                                  :: foo(1)
#endif

    ierr = 0
    CALL mp_timeset(routineN,handle)

#if defined(__parallel)
    t_start = m_walltime ( )
    my_tag = 0
    IF (PRESENT(tag)) my_tag=tag

    msglen = SIZE(msgout,1)*SIZE(msgout,2)*SIZE(msgout,3)
    IF (msglen>0) THEN
       lower1=LBOUND(msgout,1)
       lower2=LBOUND(msgout,2)
       lower3=LBOUND(msgout,3)
       CALL mpi_irecv(msgout(lower1,lower2,lower3),msglen,MPI_DOUBLE_COMPLEX,source, my_tag,&
            comm,request,ierr)
    ELSE
       CALL mpi_irecv(foo,msglen,MPI_DOUBLE_COMPLEX,source, my_tag,&
            comm,request,ierr)
    END IF
    IF ( ierr /= 0 ) CALL mp_stop( ierr, "mpi_ircv @ "//routineN )

    t_end = m_walltime ( )
    CALL add_perf(perf_id=12,count=1,time=t_end-t_start,msg_size=msglen*(2*real_8_size))
#else
    CALL mp_abort( "mp_irecv called in non parallel case" )
#endif
    CALL mp_timestop(handle)
  END SUBROUTINE mp_irecv_zm3

! *****************************************************************************
!> \brief Allocates special parallel memory
!> \param[in]  DATA      pointer to integer array to allocate 
!> \param[in]  len       number of integers to allocate
!> \param[out] stat      (optional) allocation status result
!> \author UB
! *****************************************************************************
  SUBROUTINE mp_allocate_z(DATA, len, stat)
    COMPLEX(kind=real_8), DIMENSION(:), POINTER      :: DATA
    INTEGER, INTENT(IN)                 :: len
    INTEGER, INTENT(OUT), OPTIONAL      :: stat

    CHARACTER(len=*), PARAMETER :: routineN = 'mp_allocate_z', &
         routineP = moduleN//':'//routineN

    INTEGER                             :: ierr, handle

    CALL mp_timeset(routineN,handle)

    ierr = 0
#if defined(__parallel)
    t_start = m_walltime()
    NULLIFY(DATA)
    CALL mp_alloc_mem(DATA, len, stat=ierr)
    IF (PRESENT (stat)) THEN
       stat = ierr
    ELSE
       IF (ierr /= 0) CALL mp_stop(ierr, "mpi_alloc_mem @ "//routineN)
    ENDIF
    t_end = m_walltime()
    CALL add_perf(perf_id=15, count=1, time=t_end-t_start)
#else
    ALLOCATE(DATA(len), stat=ierr)
    IF (PRESENT (stat)) THEN
       stat = ierr
    ELSE
       IF (ierr /= 0) CALL mp_stop(ierr, "ALLOCATE @ "//routineN)
    ENDIF
#endif
    CALL mp_timestop(handle)
  END SUBROUTINE mp_allocate_z

! *****************************************************************************
!> \brief Deallocates special parallel memory
!> \param[in] DATA         pointer to special memory to deallocate 
!> \param stat ...
!> \author UB
! *****************************************************************************
  SUBROUTINE mp_deallocate_z(DATA, stat)
    COMPLEX(kind=real_8), DIMENSION(:), POINTER      :: DATA
    INTEGER, INTENT(OUT), OPTIONAL      :: stat

    CHARACTER(len=*), PARAMETER :: routineN = 'mp_deallocate_z', &
         routineP = moduleN//':'//routineN

    INTEGER                             :: ierr, handle

    CALL mp_timeset(routineN,handle)

    ierr = 0
#if defined(__parallel)
    t_start = m_walltime()
    CALL mp_free_mem(DATA, ierr)
    IF (PRESENT (stat)) THEN
       stat = ierr
    ELSE
       IF (ierr /= 0) CALL mp_stop(ierr, "mpi_free_mem @ "//routineN)
    ENDIF
    NULLIFY(DATA)
    t_end = m_walltime()
    CALL add_perf(perf_id=15, count=1, time=t_end-t_start)
#else
    DEALLOCATE(DATA, stat=ierr)
    IF (PRESENT (stat)) THEN
       stat=ierr
    ELSE
       IF (ierr /= 0) CALL mp_stop(ierr, "DEALLOCATE @ "//routineN)
    ENDIF
    NULLIFY(DATA)
#endif
    CALL mp_timestop(handle)
  END SUBROUTINE mp_deallocate_z

! *****************************************************************************
!> \brief (parallel) Blocking individual file write using explicit offsets
!>        (serial) Unformatted stream write
!> \param[in] fh     file handle (file storage unit)
!> \param[in] offset file offset (position)
!> \param[in] msg    data to be writen to the file
!> \param msglen ...
!> \par MPI-I/O mapping   mpi_file_write_at
!> \par STREAM-I/O mapping   WRITE
!> \param[in](optional) msglen number of the elements of data
! *****************************************************************************
  SUBROUTINE mp_file_write_at_zv(fh, offset, msg, msglen)
    COMPLEX(kind=real_8), INTENT(IN)                      :: msg(:)
    INTEGER, INTENT(IN)                        :: fh
    INTEGER, INTENT(IN), OPTIONAL              :: msglen
    INTEGER(kind=file_offset), INTENT(IN)      :: offset

    CHARACTER(len=*), PARAMETER :: routineN = 'mp_file_write_at_zv', &
                                   routineP = moduleN//':'//routineN

    INTEGER                                    :: ierr, msg_len
#if defined(__parallel)
    INTEGER, ALLOCATABLE, DIMENSION(:)         :: status
#endif

    ierr = 0
    msg_len = SIZE(msg)
    IF (PRESENT(msglen)) msg_len = msglen
#if defined(__parallel)
    ALLOCATE(status(MPI_STATUS_SIZE))
    CALL MPI_FILE_WRITE_AT(fh, offset, msg, msg_len, MPI_DOUBLE_COMPLEX, status, ierr)
    IF (ierr .NE. 0) CALL mp_abort("mpi_file_write_at_zv @ "//routineN)
    DEALLOCATE(status)
#else
    WRITE(UNIT=fh, POS=offset+1) msg(1:msg_len)
#endif
  END SUBROUTINE mp_file_write_at_zv

! *****************************************************************************
!> \brief ...
!> \param fh ...
!> \param offset ...
!> \param msg ...
! *****************************************************************************
  SUBROUTINE mp_file_write_at_z(fh, offset, msg)
    COMPLEX(kind=real_8), INTENT(IN)               :: msg
    INTEGER, INTENT(IN)                        :: fh
    INTEGER(kind=file_offset), INTENT(IN)      :: offset

    CHARACTER(len=*), PARAMETER :: routineN = 'mp_file_write_at_z', &
                                   routineP = moduleN//':'//routineN

    INTEGER                                    :: ierr
#if defined(__parallel)
    INTEGER, ALLOCATABLE, DIMENSION(:)         :: status
#endif

    ierr = 0
#if defined(__parallel)
    ALLOCATE(status(MPI_STATUS_SIZE))
    CALL MPI_FILE_WRITE_AT(fh, offset, msg, 1, MPI_DOUBLE_COMPLEX, status, ierr)
    IF (ierr .NE. 0) CALL mp_abort("mpi_file_write_at_z @ "//routineN)
    DEALLOCATE(status)
#else
    WRITE(UNIT=fh, POS=offset+1) msg
#endif
  END SUBROUTINE mp_file_write_at_z

! *****************************************************************************
!> \brief (parallel) Blocking collective file write using explicit offsets
!>        (serial) Unformatted stream write
!> \param fh ...
!> \param offset ...
!> \param msg ...
!> \param msglen ...
!> \par MPI-I/O mapping   mpi_file_write_at_all
!> \par STREAM-I/O mapping   WRITE
! *****************************************************************************
  SUBROUTINE mp_file_write_at_all_zv(fh, offset, msg, msglen)
    COMPLEX(kind=real_8), INTENT(IN)                      :: msg(:)
    INTEGER, INTENT(IN)                        :: fh
    INTEGER, INTENT(IN), OPTIONAL              :: msglen
    INTEGER                                    :: msg_len
    INTEGER(kind=file_offset), INTENT(IN)      :: offset

    CHARACTER(len=*), PARAMETER :: routineN = 'mp_file_write_at_all_zv', &
                                   routineP = moduleN//':'//routineN

    INTEGER                                    :: ierr
#if defined(__parallel)
    INTEGER, ALLOCATABLE, DIMENSION(:)         :: status
#endif

    ierr = 0
    msg_len = SIZE(msg)
    IF (PRESENT(msglen)) msg_len = msglen
#if defined(__parallel)
    ALLOCATE(status(MPI_STATUS_SIZE))
    CALL MPI_FILE_WRITE_AT_ALL(fh, offset, msg, msg_len, MPI_DOUBLE_COMPLEX, status, ierr)
    IF (ierr .NE. 0) CALL mp_abort("mpi_file_write_at_all_zv @ "//routineN)
    DEALLOCATE(status)
#else
    WRITE(UNIT=fh, POS=offset+1) msg(1:msg_len)
#endif
  END SUBROUTINE mp_file_write_at_all_zv

! *****************************************************************************
!> \brief ...
!> \param fh ...
!> \param offset ...
!> \param msg ...
! *****************************************************************************
  SUBROUTINE mp_file_write_at_all_z(fh, offset, msg)
    COMPLEX(kind=real_8), INTENT(IN)               :: msg
    INTEGER, INTENT(IN)                        :: fh
    INTEGER(kind=file_offset), INTENT(IN)      :: offset

    CHARACTER(len=*), PARAMETER :: routineN = 'mp_file_write_at_all_z', &
                                   routineP = moduleN//':'//routineN

    INTEGER                                    :: ierr
#if defined(__parallel)
    INTEGER, ALLOCATABLE, DIMENSION(:)         :: status
#endif

    ierr = 0
#if defined(__parallel)
    ALLOCATE(status(MPI_STATUS_SIZE))
    CALL MPI_FILE_WRITE_AT_ALL(fh, offset, msg, 1, MPI_DOUBLE_COMPLEX, status, ierr)
    IF (ierr .NE. 0) CALL mp_abort("mpi_file_write_at_all_z @ "//routineN)
    DEALLOCATE(status)
#else
    WRITE(UNIT=fh, POS=offset+1) msg
#endif
  END SUBROUTINE mp_file_write_at_all_z

! *****************************************************************************
!> \brief (parallel) Blocking individual file read using explicit offsets
!>        (serial) Unformatted stream read
!> \param[in] fh     file handle (file storage unit)
!> \param[in] offset file offset (position)
!> \param[out] msg   data to be read from the file
!> \param msglen ...
!> \par MPI-I/O mapping   mpi_file_read_at
!> \par STREAM-I/O mapping   READ
!> \param[in](optional) msglen  number of elements of data
! *****************************************************************************
  SUBROUTINE mp_file_read_at_zv(fh, offset, msg, msglen)
    COMPLEX(kind=real_8), INTENT(OUT)                     :: msg(:)
    INTEGER, INTENT(IN)                        :: fh
    INTEGER, INTENT(IN), OPTIONAL              :: msglen
    INTEGER                                    :: msg_len
    INTEGER(kind=file_offset), INTENT(IN)      :: offset

    CHARACTER(len=*), PARAMETER :: routineN = 'mp_file_read_at_zv', &
                                   routineP = moduleN//':'//routineN

    INTEGER                                    :: ierr
#if defined(__parallel)
    INTEGER, ALLOCATABLE, DIMENSION(:)         :: status
#endif

    ierr = 0
    msg_len = SIZE(msg)
    IF (PRESENT(msglen)) msg_len = msglen
#if defined(__parallel)
    ALLOCATE(status(MPI_STATUS_SIZE))
    CALL MPI_FILE_READ_AT(fh, offset, msg, msg_len, MPI_DOUBLE_COMPLEX, status, ierr)
    IF (ierr .NE. 0) CALL mp_abort("mpi_file_read_at_zv @ "//routineN)
    DEALLOCATE(status)
#else
    READ(UNIT=fh, POS=offset+1) msg(1:msg_len)
#endif
  END SUBROUTINE mp_file_read_at_zv

! *****************************************************************************
!> \brief ...
!> \param fh ...
!> \param offset ...
!> \param msg ...
! *****************************************************************************
  SUBROUTINE mp_file_read_at_z(fh, offset, msg)
    COMPLEX(kind=real_8), INTENT(OUT)               :: msg
    INTEGER, INTENT(IN)                        :: fh
    INTEGER(kind=file_offset), INTENT(IN)      :: offset


    CHARACTER(len=*), PARAMETER :: routineN = 'mp_file_read_at_z', &
                                   routineP = moduleN//':'//routineN

    INTEGER                                    :: ierr
#if defined(__parallel)
    INTEGER, ALLOCATABLE, DIMENSION(:)         :: status
#endif

    ierr = 0
#if defined(__parallel)
    ALLOCATE(status(MPI_STATUS_SIZE))
    CALL MPI_FILE_READ_AT(fh, offset, msg, 1, MPI_DOUBLE_COMPLEX, status, ierr)
    IF (ierr .NE. 0) CALL mp_abort("mpi_file_read_at_z @ "//routineN)
    DEALLOCATE(status)
#else
    READ(UNIT=fh, POS=offset+1) msg
#endif
  END SUBROUTINE mp_file_read_at_z

! *****************************************************************************
!> \brief (parallel) Blocking collective file read using explicit offsets
!>        (serial) Unformatted stream read
!> \param fh ...
!> \param offset ...
!> \param msg ...
!> \param msglen ...
!> \par MPI-I/O mapping    mpi_file_read_at_all
!> \par STREAM-I/O mapping   READ
! *****************************************************************************
  SUBROUTINE mp_file_read_at_all_zv(fh, offset, msg, msglen)
    COMPLEX(kind=real_8), INTENT(OUT)                     :: msg(:)
    INTEGER, INTENT(IN)                        :: fh
    INTEGER, INTENT(IN), OPTIONAL              :: msglen
    INTEGER(kind=file_offset), INTENT(IN)      :: offset

    CHARACTER(len=*), PARAMETER :: routineN = 'mp_file_read_at_all_zv', &
                                   routineP = moduleN//':'//routineN

    INTEGER                                    :: ierr, msg_len
#if defined(__parallel)
    INTEGER, ALLOCATABLE, DIMENSION(:)         :: status
#endif

    ierr = 0
    msg_len = SIZE(msg)
    IF (PRESENT(msglen)) msg_len = msglen
#if defined(__parallel)
    ALLOCATE(status(MPI_STATUS_SIZE))
    CALL MPI_FILE_READ_AT_ALL(fh, offset, msg, msg_len, MPI_DOUBLE_COMPLEX, status, ierr)
    IF (ierr .NE. 0) CALL mp_abort("mpi_file_read_at_all_zv @ "//routineN)
    DEALLOCATE(status)
#else
    READ(UNIT=fh, POS=offset+1) msg(1:msg_len)
#endif
  END SUBROUTINE mp_file_read_at_all_zv

! *****************************************************************************
!> \brief ...
!> \param fh ...
!> \param offset ...
!> \param msg ...
! *****************************************************************************
  SUBROUTINE mp_file_read_at_all_z(fh, offset, msg)
    COMPLEX(kind=real_8), INTENT(OUT)               :: msg
    INTEGER, INTENT(IN)                        :: fh
    INTEGER(kind=file_offset), INTENT(IN)      :: offset

    CHARACTER(len=*), PARAMETER :: routineN = 'mp_file_read_at_all_z', &
                                   routineP = moduleN//':'//routineN

    INTEGER                                    :: ierr
#if defined(__parallel)
    INTEGER, ALLOCATABLE, DIMENSION(:)         :: status
#endif

    ierr = 0
#if defined(__parallel)
    ALLOCATE(status(MPI_STATUS_SIZE))
    CALL MPI_FILE_READ_AT_ALL(fh, offset, msg, 1, MPI_DOUBLE_COMPLEX, status, ierr)
    IF (ierr .NE. 0) CALL mp_abort("mpi_file_read_at_all_z @ "//routineN)
    DEALLOCATE(status)
#else
    READ(UNIT=fh, POS=offset+1) msg
#endif
  END SUBROUTINE mp_file_read_at_all_z

! *****************************************************************************
!> \brief ...
!> \param ptr ...
!> \param vector_descriptor ...
!> \param index_descriptor ...
!> \retval type_descriptor ...
! *****************************************************************************
  FUNCTION mp_type_make_z (ptr,&
       vector_descriptor, index_descriptor) &
       RESULT (type_descriptor)
    COMPLEX(kind=real_8), DIMENSION(:), POINTER                    :: ptr
    INTEGER, DIMENSION(2), INTENT(IN), OPTIONAL       :: vector_descriptor
    TYPE(mp_indexing_meta_type), INTENT(IN), OPTIONAL :: index_descriptor
    TYPE(mp_type_descriptor_type)                     :: type_descriptor

    CHARACTER(len=*), PARAMETER :: routineN = 'mp_make_type_z', &
         routineP = moduleN//':'//routineN

    INTEGER :: ierr

    ierr = 0
    NULLIFY (type_descriptor%subtype)
    type_descriptor%length = SIZE (ptr)
#if defined(__parallel)
    type_descriptor%type_handle = MPI_DOUBLE_COMPLEX
    CALL MPI_Get_address (ptr, type_descriptor%base, ierr)
    IF (ierr /= 0) CALL mp_abort("MPI_Get_address @ "//routineN)
#else
    type_descriptor%type_handle = 7
#endif
    type_descriptor%vector_descriptor(1:2) = 1
    type_descriptor%has_indexing = .FALSE.
    type_descriptor%data_z => ptr
    IF (PRESENT (vector_descriptor) .OR. PRESENT (index_descriptor)) THEN
       CALL mp_abort (routineN//": Vectors and indices NYI")
    ENDIF
  END FUNCTION mp_type_make_z

! *****************************************************************************
!> \brief Allocates an array, using MPI_ALLOC_MEM ... this is hackish
!>        as the Fortran version returns an integer, which we take to be a C_PTR
!> \param DATA           data array to allocate  
!> \param[in] len        length (in data elements) of data array allocation
!> \param[out] stat      (optional) allocation status result
! *****************************************************************************
  SUBROUTINE mp_alloc_mem_z(DATA, len, stat)
    COMPLEX(kind=real_8), DIMENSION(:), POINTER           :: DATA
    INTEGER, INTENT(IN)                      :: len
    INTEGER, INTENT(OUT), OPTIONAL           :: stat

#if defined(__parallel)
    INTEGER                                  :: size, ierr, length, &
                                                mp_info, mp_res
    INTEGER(KIND=MPI_ADDRESS_KIND)           :: mp_size
    TYPE(C_PTR)                              :: mp_baseptr

     length = MAX(len,1)
     CALL MPI_TYPE_SIZE(MPI_DOUBLE_COMPLEX, size, ierr)
     mp_size = length * size
     mp_info = MPI_INFO_NULL
     CALL MPI_ALLOC_MEM(mp_size, mp_info, mp_baseptr, mp_res)
     CALL C_F_POINTER(mp_baseptr, DATA, (/length/))
     IF (PRESENT (stat)) stat = mp_res
#else
     INTEGER                                 :: length
     length = MAX(len,1)
     IF (PRESENT (stat)) THEN
        ALLOCATE(DATA(length), stat=stat)
     ELSE
        ALLOCATE(DATA(length))
     ENDIF
#endif
   END SUBROUTINE mp_alloc_mem_z

! *****************************************************************************
!> \brief Deallocates am array, ... this is hackish 
!>        as the Fortran version takes an integer, which we hope to get by reference
!> \param DATA           data array to allocate
!> \param[out] stat      (optional) allocation status result
! *****************************************************************************
   SUBROUTINE mp_free_mem_z(DATA, stat)
    COMPLEX(kind=real_8), DIMENSION(:), &
      POINTER                                :: DATA
    INTEGER, INTENT(OUT), OPTIONAL           :: stat

#if defined(__parallel)
    INTEGER                                  :: mp_res
    CALL MPI_FREE_MEM(DATA, mp_res)
    IF (PRESENT (stat)) stat = mp_res
#else
     IF (PRESENT (stat)) THEN
        DEALLOCATE(DATA, stat=stat)
     ELSE
        DEALLOCATE(DATA)
     ENDIF
#endif
   END SUBROUTINE mp_free_mem_z
