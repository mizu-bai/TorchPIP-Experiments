program test
    implicit none
    
    real*8, allocatable::xyz(:, :)
    real*8, allocatable::r(:) ! distance vector
    real*8, allocatable::m(:) ! morse like
    real*8, allocatable::drdx(:, :)
    real*8::delta_r(3)
    character(len=2)::symb
    integer::i, j, k, natom, nr

    open(12, file="test.xyz", status="old")
    
    read(12, *) natom
    allocate(xyz(3, natom))

    read(12, *) ! skip title line

    print *, "xyz"
    do i = 1, natom
        read(12, *) symb, xyz(:, i)
        print *, xyz(:, i)
    end do

    nr = natom * (natom - 1) / 2

    allocate(r(nr))
    allocate(m(nr))

    k = 1
    do i = 1, natom - 1
        do j = i + 1, natom
            delta_r = xyz(:, i) - xyz(:, j)
            r(k) = sqrt(dot_product(delta_r, delta_r))
            k = k + 1
        end do
    end do

    print *, "distance vector"
    print *, r

    m(:) = exp(-1.0d0 * r(:))

    print *, "morse"
    print *, m

    allocate(drdx(3 * natom, nr))

    print *, "drdx", shape(drdx)

    k = 1
    drdx(:, :) = 0.0d0

    do i = 1, natom - 1
        do j = i + 1, natom
            delta_r = xyz(:, i) - xyz(:, j)
            r(k) = sqrt(dot_product(delta_r, delta_r))
            drdx((3 * i - 2):(3 * i), k) = delta_r(:) / r(k)
            drdx((3 * j - 2):(3 * j), k) = - delta_r(:) / r(k)
            
            print *, drdx(:, k)
            
            k = k + 1
            
        end do
    end do
    
    deallocate(xyz)
    deallocate(r)
    deallocate(m)
    deallocate(drdx)
    close(12)
end program test
