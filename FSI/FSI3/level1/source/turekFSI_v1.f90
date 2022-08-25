program turekFSI

   use math, only: dateTime
   use solver
   implicit none

   write (*, *) '======================================================'
   write (*, *) 'Program started at :', dateTime()

   call setupLBMvars()
   call setupD2Q9(ci, wi)
   call oppoVector(ci, kb)
   call initProbDist()

   call setupFEMvars()
   call setupBC(dofBC)
   nDofBC = nDof - size(dofBC)
   call setupElemMap(dofMapBC, coupleRange)

   ! allocate (DeltaG(nDofBC))
   ! DeltaG = 0.0d0 !Put arbitrary values here
   ! DeltaG(55:56) = [0.0d0, 3.0d0]
   ! DeltaG(43:44) = [0.0d0, 3.0d0]
   ! DeltaG(83:84) = [0.0d0, -3.0d0]
   ! DeltaG(71:72) = [0.0d0, -3.0d0]

   call demarcateElem()!(topEl, bottomEl, rightEl, leftEl)

   call undeformedInterface()!(refBounTopEl, refBounBottomEl, refBounRightEl)

   call compFSIresponse()

   write (*, *) 'Program ended at :', dateTime()
   write (*, *) '======================================================'

end program turekFSI
