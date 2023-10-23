program gen_xdmf
!
! this program generates a XDMF file for visualizing the time-averaged 3D fields.
!
! usage:
! (1) set the parameters in the file 'param_flow_mean.h90'
! (2) compile and run sequentially gen_grid.f90 and gen_xdmf_flow_mean.f90 
!     (see genview.sh)
! (3) visualize the flow with the output file viewfld_DNS_flow_mean.xmf
!     e.g. with paraview: paraview viewfld_DNS_flow_mean.xmf
! it assumes: 
!        - visualization files have the format:
!          'X_mean.bin'
!          where X is the name of the field set in param_flow_mean.h90
!        - the cans output grid file 'grid.bin' must be located 
!          in the same folder as the xdmf file
!        - output files are located in same folder as xdmf file
!
! Pedro Costa (p.simoes.costa@gmail.com)
! Seyed Morteza Habibi Khorasani
!
implicit none
include 'param_flow_mean.h90'
character(len=1), parameter :: lf = char(10)
character(len=400) :: buffer
character(len=7) :: ichar
character(len=3) :: fldname
character(len=4) :: islicechar
character(len=1) :: inormchar
integer :: ixdmf
integer :: e_io,indent
integer :: nflds
integer :: i,ii
!
ixdmf = 99 
nflds = (fldend-fldstart)/nskip + 1
indent = 0
open(unit =ixdmf, file = 'viewfld_'//trim(casename)//'_mean.xdmf',form='formatted')
write(unit=ixdmf,fmt='(A)') '<?xml version="1.0" ?>'
write(unit=ixdmf,fmt='(A)') '<!DOCTYPE Xdmf SYSTEM "Xdmf.dtd" []>'
write(unit=ixdmf,fmt='(A)') '<Xdmf xmlns:xi="http://www.w3.org/2001/XInclude" Version="2.0">'
write(unit=ixdmf,fmt='(A)') '<Domain>'
indent = indent + 4
  write(buffer,fmt='(A,3I5,A)') repeat(' ',indent)//'<Topology name="TOPO" TopologyType="3DRectMesh" Dimensions="',nz,ny,nx,'"/>'
  write(unit=ixdmf,fmt='(A)') trim(buffer)
  write(buffer,fmt='(A)') repeat(' ',indent)//'<Geometry name="GEO" GeometryType="VXVYVZ">'
  write(unit=ixdmf,fmt='(A)')trim(buffer)
  indent = indent + 4
write(buffer,fmt='(A,I1,A,1I5,A)') repeat(' ',indent)//'<DataItem Format="Binary"' // &
                                                       ' DataType="Float" Precision="',iprec,'" Endian="Native"' // &
                                                       ' Dimensions="',nx,'">'
    write(unit=ixdmf,fmt='(A)')trim(buffer)
              indent = indent + 4
                write(buffer,fmt='(A,i7.7,A)') repeat(' ',indent)//'x.bin'
                write(unit=ixdmf,fmt='(A)')trim(buffer)
                indent = indent - 4
    write(buffer,fmt='(A)') repeat(' ',indent)//'</DataItem>'
    write(unit=ixdmf,fmt='(A)')trim(buffer)
write(buffer,fmt='(A,I1,A,1I5,A)') repeat(' ',indent)//'<DataItem Format="Binary"' // &
                                                       ' DataType="Float" Precision="',iprec,'" Endian="Native"' // &
                                                       ' Dimensions="',ny,'">'
    write(unit=ixdmf,fmt='(A)')trim(buffer)
              indent = indent + 4
                write(buffer,fmt='(A,i7.7,A)') repeat(' ',indent)//'y.bin'
                write(unit=ixdmf,fmt='(A)')trim(buffer)
                indent = indent - 4
    write(buffer,fmt='(A)') repeat(' ',indent)//'</DataItem>'
    write(unit=ixdmf,fmt='(A)')trim(buffer)
write(buffer,fmt='(A,I1,A,1I5,A)') repeat(' ',indent)//'<DataItem Format="Binary"' // &
                                                       ' DataType="Float" Precision="',iprec,'" Endian="Native"' // &
                                                       ' Dimensions="',nz,'">'
    write(unit=ixdmf,fmt='(A)')trim(buffer)
              indent = indent + 4
                write(buffer,fmt='(A,i7.7,A)') repeat(' ',indent)//'z.bin'
                write(unit=ixdmf,fmt='(A)')trim(buffer)
                indent = indent - 4
    write(buffer,fmt='(A)') repeat(' ',indent)//'</DataItem>'
    write(unit=ixdmf,fmt='(A)')trim(buffer)
    indent = indent - 4
  write(buffer,fmt='(A)') repeat(' ',indent)//'</Geometry>'
  write(unit=ixdmf,fmt='(A)')trim(buffer)
  write(buffer,fmt='(A)') repeat(' ',indent)//'<Grid Name="TimeSeries" GridType="Collection" CollectionType="Temporal">'
  write(unit=ixdmf,fmt='(A)')trim(buffer)
  indent = indent + 4
    write(buffer,fmt='(A)') repeat(' ',indent)//'<Time TimeType="List">'
    write(unit=ixdmf,fmt='(A)')trim(buffer)
    indent = indent + 4
      write(buffer,fmt='(A,I6,A)') repeat(' ',indent)//'<DataItem Format="XML" NumberType="Float" Dimensions="',nflds*nscal,'">'
      write(unit=ixdmf,fmt='(A)')trim(buffer) 
      write(buffer,fmt='(A,I6)') repeat(' ',indent)
      write(ixdmf,fmt='(A)',advance='no') trim(buffer)
      do i = fldstart,fldend,nskip !1,nflds
        write(ixdmf,fmt='(E15.6)',advance='no') t0 + 1.d0*i*dt!1.d0*(i-nskip)*dt!1.*i
      enddo
      write(buffer,fmt='(A)') repeat(' ',indent)//'</DataItem>'
      write(unit=ixdmf,fmt='(A)')trim(buffer)
      indent = indent - 4
    write(buffer,fmt='(A)') repeat(' ',indent)//'</Time>'
    write(unit=ixdmf,fmt='(A)')trim(buffer)
    do i = fldstart,fldend,nskip
      write(ichar,fmt='(i7.7)') i
      write(buffer,fmt='(A)') repeat(' ',indent)//'<Grid Name="T'//ichar//'" GridType="Uniform">'
      write(unit=ixdmf,fmt='(A)')trim(buffer)
      indent = indent + 4
        write(buffer,fmt='(A)') repeat(' ',indent)//'<Topology Reference="/Xdmf/Domain/Topology[1]"/>'
        write(unit=ixdmf,fmt='(A)')trim(buffer)
        write(buffer,fmt='(A)') repeat(' ',indent)//'<Geometry Reference="/Xdmf/Domain/Geometry[1]"/>'
        write(unit=ixdmf,fmt='(A)')trim(buffer)
        do ii = 1,nscal
          write(buffer,fmt='(A)') repeat(' ',indent)//'<Attribute Name="'//scalname(ii)//'" Center="Node">'
          write(unit=ixdmf,fmt='(A)')trim(buffer)
          indent = indent + 4
            write(buffer,fmt='(A,I1,A,3I5,A)') repeat(' ',indent)//'<DataItem Format="Binary"' // &
                                                                   ' DataType="Float" Precision="',iprec,'" Endian="Native"' // &
                                                                   ' Dimensions="',nz,ny,nx,'">'
            write(unit=ixdmf,fmt='(A)')trim(buffer)
            indent = indent + 4
              write(buffer,fmt='(A,i7.7,A)') repeat(' ',indent)//scalname(ii)//'.bin'
              write(unit=ixdmf,fmt='(A)')trim(buffer)
              indent = indent - 4
            write(buffer,fmt='(A)') repeat(' ',indent)//'</DataItem>'
            write(unit=ixdmf,fmt='(A)')trim(buffer)
            indent = indent - 4
          write(buffer,fmt='(A)') repeat(' ',indent)//'</Attribute>'
          write(unit=ixdmf,fmt='(A)')trim(buffer)
        enddo
        indent = indent - 4
      write(buffer,fmt='(A)') repeat(' ',indent)//'</Grid>'
      write(unit=ixdmf,fmt='(A)')trim(buffer)
    enddo
    indent = indent - 4
  write(buffer,fmt='(A)') repeat(' ',indent)//'</Grid>'
  write(unit=ixdmf,fmt='(A)')trim(buffer)
indent = indent - 4
write(buffer,fmt='(A)') repeat(' ',indent)//'</Domain>'
write(unit=ixdmf,fmt='(A)')trim(buffer)
write(buffer,fmt='(A)') repeat(' ',indent)//'</Xdmf>'
write(unit=ixdmf,fmt='(A)')trim(buffer)
!
close(ixdmf)
end program gen_xdmf
