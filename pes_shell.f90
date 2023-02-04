module pes_shell
implicit none

contains

  subroutine evx(xyz,x)
    real,dimension(:,:),intent(in)::xyz
    real,dimension(size(xyz,2)*(size(xyz,2)-1)/2),intent(out)::x
    real,dimension(3)::dr
    real::a0  ! the same as the fitting code
    integer::i,j,k

    a0 = 1.0d0

    k = 1
    do i=1,size(xyz,2)-1
       do j=i+1,size(xyz,2)
          dr = xyz(:,i) - xyz(:,j)
          x(k) = sqrt(dot_product(dr,dr))
          k = k+1
       end do
    end do

    do i=1,size(x)
       x(i)=exp(-x(i)/a0)
    end do
  end subroutine evx

  subroutine evdrdx(xyz, drdx)
    real,dimension(:,:),intent(in)::xyz
    !::::::::::::::::::::
    real,dimension(size(xyz,2)*(size(xyz,2)-1)/2)::r,x
    real,dimension(3,size(xyz,2)*(size(xyz,2)-1)/2)::dr
    real,dimension(size(xyz,2)*3,size(xyz,2)*(size(xyz,2)-1)/2),intent(out)::drdx
    real,dimension(1:208)::p   ! change to number of popynomials
                               ! (size of p in bemsa.f90)
    real,dimension(1:496)::m  ! change to number of monomials
                              ! (size of m in bemsa.f90
    real::a0
    integer::i,j,k

    a0 = 1.0d0

    k = 1
    drdx = 0.d0
    do i=1,size(xyz,2)-1
       do j=i+1,size(xyz,2)
          dr(:,k) = xyz(:,i) - xyz(:,j)
          r(k) = sqrt(dot_product(dr(:,k),dr(:,k)))

          drdx(3*i-2:3*i,k) = dr(:,k)/r(k)
          drdx(3*j-2:3*j,k) = -drdx(3*i-2:3*i,k)
          k = k+1
       end do
    end do

  end subroutine

end module pes_shell
