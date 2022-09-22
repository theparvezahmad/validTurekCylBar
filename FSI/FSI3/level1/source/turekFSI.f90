program turekFSI

  use solver
  implicit none

  write (*, *) '======================================================'
  write (*, *) 'Program started at :', dateTime()

  call setupD2Q9()
  call setupLBMvars()
  call initProbDist()

  call setupBC()
  call setupFEMvars()
  call setupElemMap()

  call demarcateElem()!(topEl, bottomEl, rightEl, leftEl)
  call undeformedInterface()!(refBounTopEl, refBounBottomEl, refBounRightEl)
  call compFSIresponse()

  write (*, *) 'Program ended at :', dateTime()
  write (*, *) '======================================================'

end program turekFSI
