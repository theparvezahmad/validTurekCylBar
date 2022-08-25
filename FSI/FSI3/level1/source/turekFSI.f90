program turekFSI

   use solver
   implicit none

   write (*, *) '======================================================'
   write (*, *) 'Program started at :', dateTime()

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
