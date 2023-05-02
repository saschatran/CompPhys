! PROJECT: Determination of the melting Temperature of Silicon
!          by Molecular Dynamics Simulation

! - Harmonic oscilator with Verlet algorithm

! @author: Sascha Tran
! @lecture: Computational Physics
! @semester: HS22
! @date: 22/04/2023
! @see: Lecture slides 0-131+


module verletPlus
    contains


    subroutine verlet(r, h, v, m, k)
        ! Verlet algorithm on slide 0-112

        implicit none

        real(8), intent(inout) :: r, v
        real(8), intent(in) :: h, k, m
        real(8) :: f1, f2

        ! Eq. (97)
        call force(r, k, f1)
        r = r + h * v + h**2 / (2 * m) * f1

        ! Eq. (98)
        call force(r, k, f2)
        v = v + h / (2 * m) * (f2 + f1)

    end subroutine

    subroutine force(r, k, f)

        implicit none
        
        real(8), intent(inout) :: f
        real(8), intent(in) :: r, k

        f = -k * r
        
    end subroutine

    subroutine calcETot(k, r, m, v, eTot)

        implicit none

        real(8), intent(inout) :: eTot
        real(8), intent(in) :: k, r, m, v

        ! eTot = ePot + eKin
        eTot = 0.5 * k * r**2 + 0.5 * m * v**2

    end subroutine



end module


program harmOscillator
    ! v: Velocity
    ! r: Position
    ! h: Time step
    ! m: Mass
    ! k: Spring constant
    ! f: force

    ! MD over 30 "seconds"

    use verletPlus

    implicit none

    integer :: i, j
    integer, parameter :: len = 30000
    real(8) :: v, h, r, m, k, ePot, eKin, eTot, maxETot, minEtot

    ! Block 1

    r = 1.d0
    v = 0.d0
    m = 1.d0
    h = 1.d-3
    k = 2.d0

    ! Periodic oscillation check
    open(20, file="verlet.csv")

    write(20, '(*(G0.7,:,","))') "Time", "Total energy", "Position", "Velocity"

    call calcETot(k, r, m, v, eTot)
    write(20, '(*(G0.7,:,","))') 0.0, eTot, r, v

    do i = 1, len
        call verlet(r, h, v, m, k)
        call calcETot(k, r, m, v, eTot)
        write(20, '(*(G0.7,:,","))') i, eTot, r, v
    end do

    close(20)

    ! Block 2

    ! Check tot energy conservation with different time steps
    open(21, file="verletETot.csv")
    write(21, '(*(G0.7,:,","))') "Time step", "Max. difference in total energy"

    do j = 1, -20, -1

        r = 1.d0
        v = 0.d0
        m = 1.d0
        h = 1.d0**j
        k = 2.d0

        call calcETot(k, r, m, v, eTot)
        minEtot = eTot
        maxETot = eTot

        do i = 1, len
            call verlet(r, h, v, m, k)
            call calcETot(k, r, m, v, eTot)
            if (eTot > maxETot) then
                maxETot = eTot
            else if (eTot < minEtot) then
                minEtot = eTot
            end if
        end do
        write(21, '(*(G0.7,:,","))') h, maxETot - minEtot
    end do

end program harmOscillator
        