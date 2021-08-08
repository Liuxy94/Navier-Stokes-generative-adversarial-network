! python call looks like:
!   resid_u, resid_cty, resid_c = residual_all_eqn( &
!                       u_cell, c, &
!                       sigma_u, s_u, p, density_u, mu_u, &
!                       sigma_c, s_c, density_c, mu_c, rblank, &
!                       dx,dy,dz, dt, ndim, nc, nx,ny,nz, ntime)
subroutine residual_all_eqn(resid_u, resid_cty, resid_c, &
u_cell, c, &
sigma_u, s_u, p, density_u, mu_u, &
sigma_c, s_c, density_c, mu_c, rblank, &
dx,dy,dz, dt, ndim, nc, nx,ny,nz, ntime)
! *******************************************************************************************
! This subroutine calculates the finite difference residual for the momentum equations (resid_u)
! the cty eqn (resid_cty) and the nc concentration eqns (resid_c).
! The residuals resid_u, resid_cty, resid_c are calculated between time levels and thus
! there are ntime-1 of them.
! *******************************************************************************************
! u_cell = cell-wise velocities at each time level.
! c = cell-wise concetration field at each time levels.
! sigma_u = absorption of momentum coefficient - cell wise.
! s_u = source of momentum coefficient - cell wise.
! p = pressure for momentum eqns - cell wise.
! density_u = density of momentum coefficient - cell wise.
! mu_u = dynamic viscocity of momentum eqns - cell wise.
!
! for the concentration eqns c:
! sigma_c = absorption of concentration coefficient - cell wise.
! s_c = source of concentration field - cell wise.
! density_c = density of momentum coefficient - cell wise.
! mu_c = diffusion coefficient of concentration eqns - cell wise.
!
! rblank(i,j,k)=1. then blank out the cell (set residuals to zero); rblank(i,j,k)=0. then calculate residual.
!
! becuase we are using a staggered mesh for velocity u is defined on faces with dimensions u(ndim,nx+1,ny+1,nz+1,ntime)
! and similarly u_mean is on the faces of the cells.
! dx,dy,dz = cells sizes.
! dt= time step size.
! ndim=no of dimensions e.g. 2 in 2D and 3 in 3D.
! nc=no of concentration equations.
! nx,ny,nz dimensions of the cells in x-,y-,z-directions
! ntime=no of time steps.
implicit none
integer, intent( in ) :: ndim, nc, nx,ny,nz, ntime
real, intent( in ) :: u_cell(ndim,nx,ny,nz,ntime), c(nc,nx,ny,nz,ntime)
real, intent( in ) :: sigma_u(ndim,nx,ny,nz,ntime), s_u(ndim,nx,ny,nz,ntime), p(nx,ny,nz,ntime)
real, intent( in ) :: density_u(ndim,nx,ny,nz,ntime)
real, intent( in ) :: mu_u(ndim,nx,ny,nz,ntime), rblank(nx,ny,nz)
real, intent( in ) :: sigma_c(nc,nx,ny,nz,ntime), s_c(nc,nx,ny,nz,ntime), density_c(nc,nx,ny,nz,ntime)
real, intent( in ) :: mu_c(nc,nx,ny,nz,ntime)
real, intent( in ) :: dx,dy,dz, dt
real, intent( out ) :: resid_u(ndim,nx,ny,nz,ntime-1), resid_cty(nx,ny,nz,ntime-1), resid_c(nc,nx,ny,nz,ntime-1)
! local variables...
real, parameter :: toler = 1.e-10
integer idim,pidim,itime,ic

REAL, ALLOCATABLE :: u_cell_mid(:,:,:,:), zero(:,:,:)
REAL, ALLOCATABLE :: sigma_u_mid(:,:,:), s_u_mid(:,:,:), p_mid(:,:,:), density_u_mid(:,:,:), mu_u_mid(:,:,:)
! c:
REAL, ALLOCATABLE :: c_mid(:,:,:)
REAL, ALLOCATABLE :: sigma_c_mid(:,:,:), s_c_mid(:,:,:), density_c_mid(:,:,:), mu_c_mid(:,:,:)

! calculate the face values u from the cell centred values u_cell.
allocate( u_cell_mid(ndim,nx,ny,nz), zero(nx,ny,nz) )
allocate( sigma_u_mid(nx,ny,nz), s_u_mid(nx,ny,nz), p_mid(nx,ny,nz) )
allocate( density_u_mid(nx,ny,nz), mu_u_mid(nx,ny,nz) )
! c:
allocate( c_mid(nx,ny,nz) )
allocate(sigma_c_mid(nx,ny,nz), s_c_mid(nx,ny,nz), density_c_mid(nx,ny,nz), mu_c_mid(nx,ny,nz) )
zero=0.0

do itime=1,ntime-1

u_cell_mid = 0.5*( u_cell(:,:,:,:,itime+1) + u_cell(:,:,:,:,itime) )

! momentum/velocity residual...
do idim=1,ndim
sigma_u_mid = 0.5*( sigma_u(idim,:,:,:,itime+1) + sigma_u(idim,:,:,:,itime) )
s_u_mid = 0.5*( s_u(idim,:,:,:,itime+1) + s_u(idim,:,:,:,itime) )
p_mid = 0.5*( p(:,:,:,itime+1) + p(:,:,:,itime) )
density_u_mid = 0.5*( density_u(idim,:,:,:,itime+1) + density_u(idim,:,:,:,itime) )
mu_u_mid = 0.5*( mu_u(idim,:,:,:,itime+1) + mu_u(idim,:,:,:,itime) )
pidim = idim

call residual_c_eqn(resid_u(idim,:,:,:,itime), u_cell(idim,:,:,:,itime+1), &
u_cell(idim,:,:,:,itime), u_cell_mid(idim,:,:,:), u_cell_mid, &
sigma_u_mid, s_u_mid, p_mid, density_u_mid, mu_u_mid, rblank, &
dx,dy,dz, dt, pidim, ndim, nx,ny,nz)

end do ! do idim=1,ndim
! cty eqn...
call residual_cty(resid_cty(:,:,:,itime), u_cell_mid, rblank, &
dx,dy,dz, ndim, nx,ny,nz)
!
do ic=1,nc ! concentration eqns...
! concentration c:
c_mid = 0.5*( c(ic,:,:,:,itime+1) + c(ic,:,:,:,itime) )
sigma_c_mid(:,:,:) = 0.5*( sigma_c(ic,:,:,:,itime+1) + sigma_c(ic,:,:,:,itime) )
s_c_mid(:,:,:)  = 0.5*( s_c(ic,:,:,:,itime+1) + s_c(ic,:,:,:,itime) )
density_c_mid(:,:,:) = 0.5*( density_c(ic,:,:,:,itime+1) + density_c(ic,:,:,:,itime) )
mu_c_mid(:,:,:) = 0.5*( mu_c(ic,:,:,:,itime+1) + mu_c(ic,:,:,:,itime) )
pidim = 0

call residual_c_eqn(resid_c(ic,:,:,:,itime), c(ic,:,:,:,itime+1), c(ic,:,:,:,itime), &
c_mid, u_cell_mid, &
sigma_c_mid, s_c_mid, zero, density_c_mid, mu_c_mid, rblank, &
dx,dy,dz, dt, pidim, ndim, nx,ny,nz)

end do ! do ic=1,nc

end do ! do itime=1,ntime-1

return
end subroutine residual_all_eqn
!
!
!
!
! python call looks like:
!      resid = residual_c_eqn(c_new, c_old, &
!                       c, u_cell, &
!                       sigma, s, p, density, mu, rblank, &
!                       dx,dy,dz, dt, pidim, ndim, nx,ny,nz)
subroutine residual_c_eqn(resid, c_new, c_old, &
c, u_cell, &
sigma, s, p, density, mu, rblank, &
dx,dy,dz, dt, pidim, ndim, nx,ny,nz)
! *******************************************************************************************
! This subroutine calculates the finite difference residual for the concentration equations (resid)
! but can also be used for momentum eqns.
! The residual resid is calculated between time levels.
! *******************************************************************************************
! c_new = cell-wise concentration at end of time step.
! c_old = cell-wise concentration at start of time step.
! c = cell-wise concentration averaged between time levels.
! u_cell = cell-wise velocities averaged between time levels.
!
! material coeficients:
! sigma = absorption of concentration coefficient - cell wise.
! s = source of concentration field - cell wise.
! p = pressure - cell wise.
! density = density of momentum coefficient - cell wise.
! mu = diffusion coefficient of concentration eqns - cell wise.
!
! rblank(i,j,k)=1. then blank out the cell (set residuals to zero); rblank(i,j,k)=0. then calculate residual.
!
! dx,dy,dz = cells sizes.
! dt= time step size.
! ndim=no of dimensions e.g. 2 in 2D and 3 in 3D.
! nx,ny,nz dimensions of the cells in x-,y-,z-directions
implicit none
integer, intent( in ) :: pidim, ndim, nx,ny,nz
real, intent( in ) :: c_new(nx,ny,nz), c_old(nx,ny,nz), c(nx,ny,nz), u_cell(ndim,nx,ny,nz)
real, intent( in ) :: sigma(nx,ny,nz), s(nx,ny,nz), p(nx,ny,nz), density(nx,ny,nz)
real, intent( in ) :: mu(nx,ny,nz), rblank(nx,ny,nz)
real, intent( in ) :: dx,dy,dz, dt
real, intent( out ) :: resid(nx,ny,nz)
! local variables...
real, parameter :: toler = 1.e-10
real dx_dim(ndim)
integer i_delta ! function
integer k,j,i,idim, k_start, k_finish

dx_dim(1)=dx; dx_dim(2)=dy; if(ndim==3) dx_dim(3)=dz
! calculate the face values u from the cell centred values u_cell.
k_start=2
k_finish=nz-1
if((ndim==2).and.(nz==1)) then
k_start=1
k_finish=1
endif

resid=0.0
do k=k_start,k_finish
do j=2,ny-1
do i=2,nx-1
resid(i,j,k) = density(i,j,k)*(c_new(i,j,k)-c_old(i,j,k))/dt + sigma(i,j,k)*c(i,j,k) - s(i,j,k)
do idim=1,ndim
! non-linear terms...
resid(i,j,k) =resid(i,j,k) + density(i,j,k)*u_cell(idim,i,j,k) &
* (   c(i+i_delta(1,idim),j+ i_delta(2,idim ), k+i_delta(3,idim )) &
- c(i-i_delta(1,idim),j- i_delta(2,idim ), k-i_delta(3,idim )) )/(2.*dx_dim(idim)) &
! diffusion term...
+ mu(i,j,k)*(  -c(i+i_delta(1,idim),j+ i_delta(2,idim ), k+i_delta(3,idim )) &
+2.*c(i,j,k)  &
-c(i-i_delta(1,idim),j- i_delta(2,idim ), k-i_delta(3,idim )) )/(dx_dim(idim)**2) &
! pressure term...
+ real(i_delta(pidim,idim))*(   p(i+i_delta(1,idim),j+ i_delta(2,idim ), k+i_delta(3,idim )) &
- p(i-i_delta(1,idim),j- i_delta(2,idim ), k-i_delta(3,idim )) )/(2.*dx_dim(idim))
end do
end do
end do
end do
resid=resid*(1.0-rblank)

return
end subroutine residual_c_eqn
!
!
!
! python call looks like:
!      resid_cty = residual_cty(u_cell, rblank, &
!                         dx,dy,dz, ndim, nx,ny,nz)
subroutine residual_cty(resid_cty, u_cell, rblank, &
dx,dy,dz, ndim, nx,ny,nz)
! *******************************************************************************************
! This subroutine calculates the finite difference residual for the cty equation (resid_cty).
! *******************************************************************************************
! u_cell = cell-wise velocities averaged between time levels.
!
! rblank(i,j,k)=1. then blank out the cell (set residuals to zero) =0. calcul residual.
!
! dx,dy,dz = cells sizes.
! dt= time step size.
! ndim=no of dimensions e.g. 2 in 2D and 3 in 3D.
implicit none
integer, intent( in ) :: ndim, nx,ny,nz
real, intent( in ) :: u_cell(ndim,nx,ny,nz)
real, intent( in ) :: rblank(nx,ny,nz)
real, intent( in ) :: dx,dy,dz
real, intent( out ) :: resid_cty(nx,ny,nz)
! local variables...
real, parameter :: toler = 1.e-10
real dx_dim(ndim)
integer i_delta ! function
integer k,j,i,idim, k_start, k_finish

dx_dim(1)=dx; dx_dim(2)=dy; if(ndim==3) dx_dim(3)=dz
! calculate the face values u from the cell centred values u_cell.
k_start=2
k_finish=nz-1
if((ndim==2).and.(nz==1)) then
k_start=1
k_finish=1
endif

resid_cty=0.0
do k=k_start,k_finish
do j=2,ny-1
do i=2,nx-1
do idim=1,ndim
! non-linear terms...
resid_cty(i,j,k) =resid_cty(i,j,k)  &
+ (   u_cell(idim,i+i_delta(1,idim),j+ i_delta(2,idim ), k+i_delta(3,idim )) &
- u_cell(idim,i-i_delta(1,idim),j- i_delta(2,idim ), k-i_delta(3,idim )) )/(2.*dx_dim(idim))
end do
end do
end do
end do
resid_cty=resid_cty*(1.-rblank)

return
end subroutine residual_cty
!
! 
