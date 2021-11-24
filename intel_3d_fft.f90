!  This file contains:
!  1.  SUBROUTINE fftc2c
!  2.  SUBROUTINE ifftc2c
!  3.  SUBROUTINE fftc2c_inplace
!  4.  SUBROUTINE ifftc2c_inplace
!
!  CREATOR:  Yu Yu @ SHAO 2011.10.29
!  LAST MODIFIER:  Yu Yu @ SHAO 2014.06.18

!======================================
SUBROUTINE fftc2c(indata1d,outdata1d,L)
!  This SUBROUTINE compute forward FFT from COMPLEX to COMPLEX
!  with forward scale factor 1.0
!  by using Intel MKL.
!  'mkl_dfti.f90' should be included

Use MKL_DFTI

INTEGER::L(3)
COMPLEX(4)::indata1d(L(1)*L(2)*L(3))
COMPLEX(4)::outdata1d(L(1)*L(2)*L(3))
type(DFTI_DESCRIPTOR),POINTER::my_desc_handle
INTEGER::status

print*,'begin fftc2c'

status = DftiCreateDescriptor( my_desc_handle, DFTI_SINGLE, DFTI_COMPLEX, 3, L )
status = DftiSetValue( my_desc_handle, DFTI_PLACEMENT, DFTI_NOT_INPLACE )
status = DftiCommitDescriptor(my_desc_handle)
status = DftiComputeForward( my_desc_handle, indata1d, outdata1d)
status = DftiFreeDescriptor(my_desc_handle)

print*,'fftc2c success'

END SUBROUTINE


!=======================================
SUBROUTINE ifftc2c(indata1d,outdata1d,L)
!  This SUBROUTINE compute backward FFT from COMPLEX to COMPLEX
!  with forward scale factor 1.0/N_grid
!  by using Intel MKL.
!  'mkl_dfti.f90' should be included

Use MKL_DFTI

INTEGER::L(3)
COMPLEX(4)::indata1d(L(1)*L(2)*L(3))
COMPLEX(4)::outdata1d(L(1)*L(2)*L(3))
REAL(4)::factor
type(DFTI_DESCRIPTOR),POINTER::my_desc_handle
INTEGER::status

print*,'begin ifftc2c'

factor=1.0/(float(L(1))*float(L(2))*float(L(3)))

status = DftiCreateDescriptor( my_desc_handle, DFTI_SINGLE, DFTI_COMPLEX, 3, L )
status = DftiSetValue( my_desc_handle, DFTI_PLACEMENT, DFTI_NOT_INPLACE )
status = DftiSetValue( my_desc_handle, DFTI_BACKWARD_SCALE, factor )
status = DftiCommitDescriptor(my_desc_handle)
status = DftiComputeBackward( my_desc_handle, indata1d, outdata1d)
status = DftiFreeDescriptor(my_desc_handle)

print*,'ifftc2c success'

END SUBROUTINE


!======================================
SUBROUTINE fftc2c_inplace(indata1d,L)
!  This SUBROUTINE compute forward FFT from COMPLEX to COMPLEX
!  with forward scale factor 1.0
!  by using Intel MKL.
!  'mkl_dfti.f90' should be included

Use MKL_DFTI

INTEGER::L(3)
COMPLEX(4)::indata1d(L(1)*L(2)*L(3))
type(DFTI_DESCRIPTOR),POINTER::my_desc_handle
INTEGER::status

print*,'begin fftc2c'

status = DftiCreateDescriptor( my_desc_handle, DFTI_SINGLE, DFTI_COMPLEX, 3, L )
status = DftiSetValue( my_desc_handle, DFTI_PLACEMENT, DFTI_INPLACE )
status = DftiCommitDescriptor(my_desc_handle)
status = DftiComputeForward( my_desc_handle, indata1d)
status = DftiFreeDescriptor(my_desc_handle)

print*,'fftc2c success'

END SUBROUTINE


!=======================================
SUBROUTINE ifftc2c_inplace(indata1d,L)
!  This SUBROUTINE compute backward FFT from COMPLEX to COMPLEX
!  with forward scale factor 1.0/N_grid
!  by using Intel MKL.
!  'mkl_dfti.f90' should be included

Use MKL_DFTI

INTEGER::L(3)
COMPLEX(4)::indata1d(L(1)*L(2)*L(3))
REAL(4)::factor
type(DFTI_DESCRIPTOR),POINTER::my_desc_handle
INTEGER::status

print*,'begin ifftc2c'

factor=1.0/(float(L(1))*float(L(2))*float(L(3)))

status = DftiCreateDescriptor( my_desc_handle, DFTI_SINGLE, DFTI_COMPLEX, 3, L )
status = DftiSetValue( my_desc_handle, DFTI_PLACEMENT, DFTI_INPLACE )
status = DftiSetValue( my_desc_handle, DFTI_BACKWARD_SCALE, factor )
status = DftiCommitDescriptor(my_desc_handle)
status = DftiComputeBackward( my_desc_handle, indata1d)
status = DftiFreeDescriptor(my_desc_handle)

print*,'ifftc2c success'

END SUBROUTINE

!EOF
