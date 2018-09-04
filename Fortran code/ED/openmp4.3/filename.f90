!This page contains the subroutine which give out the file name.

!1.getbasename:give out the basename of the file name
!2.geteigvname:give out the eigvname appending to basename--write the eigvalues.
!3.geteigsname:give out the eigsname appending to basename--write the eigstates.
!4.getcorrSzname:give out the corrSzname appending to basename--write the corrSz
!5.getcorrSdname:give out the corrSdname appending to basename--write the corrSd
!6.getcorrDdname:give out the corrDdname appending to basename--write the corrDds
!--------------------------------------------------
!getbasename:give out the basename of the file name
!--------------------------------------------------
subroutine getbasename()
use param
use io_module
implicit none
call createFileName(basename,'hubbard_lan_')
call appendBaseName(basename,'set_',set)
call appendBaseName(basename,'Nlat_',Nlattice)
call appendBaseName(basename,'Nhop_',Nhop)
call appendBaseName(basename,'D_',Dimen)
call appendBaseName(basename,'Nl1_',Nl(1))
call appendBaseName(basename,'Nl2_',Nl(2))
call appendBaseName(basename,'Nl3_',Nl(3))
call appendBaseName(basename,'k1_',3,kbound(1))
call appendBaseName(basename,'k2_',3,kbound(2))
call appendBaseName(basename,'k3_',3,kbound(3))
call appendBaseName(basename,'_')
call appendBaseName(basename,'PS_')
call appendBaseName(basename,ps)
call appendBaseName(basename,'_')
call appendBaseName(basename,'PK_')
call appendBaseName(basename,pk)
call appendBaseName(basename,'_')
call appendBaseName(basename,'2rS_',two_rSS)
call appendBaseName(basename,'kp1_',k_keep(1))
call appendBaseName(basename,'kp2_',k_keep(2))
call appendBaseName(basename,'kp3_',k_keep(3))
call appendBaseName(basename,'Nup_',Nup)
call appendBaseName(basename,'Ndn_',Ndn)
call appendBaseName(basename,'t1_',3,dble(kin))
call appendBaseName(basename,'U_',3,onsitU)
call appendBaseName(basename,'Ncite_',Nexcit)
end subroutine getbasename


!-----------------------------------------------------------------------------
!geteigvname:give out the eigvname appending to basename--write the eigvalues.
!-----------------------------------------------------------------------------
subroutine geteigvname()
use param
use io_module
implicit none
call copyName(BaseName,eigvname)
call appendBaseName(eigvname,'_Energy.dat')
end subroutine geteigvname


!-----------------------------------------------------------------------------
!geteigsname:give out the eigsname appending to basename--write the eigstates.
!-----------------------------------------------------------------------------
subroutine geteigsname()
use param
use io_module
implicit none
call copyName(BaseName,eigsname)
call appendBaseName(eigsname,'Nei_',Nei)
call appendBaseName(eigsname,'_Eigstat.dat')
end subroutine geteigsname

!-----------------------------------------------------------------------------
!getcorrSzname:give out the corrSzname appending to basename--write the corrSz
!-----------------------------------------------------------------------------
subroutine getcorrSzname()
use param
use io_module
implicit none
call copyName(BaseName,corrSzname)
call appendBaseName(corrSzname,'Nei_',Nei)
call appendBaseName(corrSzname,'_coorSz.dat')
end subroutine getcorrSzname



!-----------------------------------------------------------------------------
!getcorrSkname:give out the corrSkname appending to basename--write the corrSk
!-----------------------------------------------------------------------------
subroutine getcorrSkname()
use param
use io_module
implicit none
call copyName(BaseName,corrSkname)
call appendBaseName(corrSkname,'Nei_',Nei)
call appendBaseName(corrSkname,'_coorSk.dat')
end subroutine getcorrSkname



!-----------------------------------------------------------------------------
!getcorrSdname:give out the corrSdname appending to basename--write the corrSd
!-----------------------------------------------------------------------------
subroutine getcorrSdname()
use param
use io_module
implicit none
call copyName(BaseName,corrSdname)
call appendBaseName(corrSdname,'Nei_',Nei)
call appendBaseName(corrSdname,'_coorSd.dat')
end subroutine getcorrSdname



!------------------------------------------------------------------------------
!getcorrDdname:give out the corrDdname appending to basename--write the corrDds
!------------------------------------------------------------------------------
subroutine getcorrDdname()
use param
use io_module
implicit none
call copyName(BaseName,corrDdname)
call appendBaseName(corrDdname,'Nei_',Nei)
call appendBaseName(corrDdname,'_coorDds.dat')
end subroutine getcorrDdname

!---------------------------------------------------------------
!getcoorcicjname:give out the coorcicjname appneding to basename
!---------------------------------------------------------------
subroutine getcoorcicjname()
use param
use io_module
implicit none
call copyName(BaseName,coorcicjname)
call appendBaseName(coorcicjname,'Nei_',Nei)
call appendBaseName(coorcicjname,'_coorcicj.dat')
end subroutine getcoorcicjname


!-----------------------------------------------------------
!getcoornkname:give out the coornkname appneding to basename
!-----------------------------------------------------------
subroutine getcoornkname()
use param
use io_module
implicit none
call copyName(BaseName,coornkname)
call appendBaseName(coornkname,'Nei_',Nei)
call appendBaseName(coornkname,'_coornk.dat')
end subroutine getcoornkname
