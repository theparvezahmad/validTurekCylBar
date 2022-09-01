program turekFSI

   use solver
   implicit none

   write (*, *) '======================================================'
   write (*, *) 'Program started at :', dateTime()

   ! write (*, *) 'Total active threads:', omp_get_max_threads()
   ! !$omp parallel
   ! write (*, *) 'Hello from process:', omp_get_thread_num()
   ! !$omp end parallel

   call setupLBMvars()
   call setupD2Q9(ci, wi, kb)
   call initProbDist()

   call setupBC(dofBC)
   call setupFEMvars()
   call setupElemMap(dofMapBC, coupleRange)

   call demarcateElem()!(topEl, bottomEl, rightEl, leftEl)
   call undeformedInterface()!(refBounTopEl, refBounBottomEl, refBounRightEl)
   call compFSIresponse()

   write (*, *) 'Program ended at :', dateTime()
   write (*, *) '======================================================'

end program turekFSI
