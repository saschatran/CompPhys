! PROJECT: Sorting continously growing data sets

! @author: Sascha Tran
! @lecture: Computational Physics
! @semester: HS22
! @date: 22/04/2023
! @see: Lecture slides 0-18+

module sortingD
contains

    subroutine allocateArray(arr, length)
        ! Adjust length of array to the given

        implicit none
        real(8), dimension(:), allocatable, intent(inout) :: arr
        integer :: length

        ! If array is already allocated, deallocate it first
        if (allocated(arr)) then
            deallocate(arr)
        end if

        ! Allocate array with the given length
        allocate(arr(length))

    end subroutine

    subroutine binarySearch(arr, length, element, position)
        ! implements a binary search algorithm to find the position of 'element' in 'arr'
        ! 'arr' has to be sorted in ascending order
        ! returns the position of 'element' in 'arr' in 'position'
        ! if 'element' is not in 'arr', the position of the next larger element is returned

        implicit none 
        integer, intent(in) :: length
        real(8), dimension(length), intent(in) :: arr
        real(8), intent(in) :: element
        integer, intent(out) :: position
        
        integer :: lo, mid, hi 
        
        ! initialize
        lo = 1 
        hi = length
        mid = lo + (hi - lo) / 2
        
        ! check if the element is in the middle of the array
        do while (lo <= hi)

            ! Element is greater than the middle element, search in the right half
            if (element - arr(mid) > 1.d-10) then 
                lo = mid + 1 
            
            ! Element is smaller than the middle element, search in the left half
            else if (element - arr(mid) < -1.d-10) then 
                hi = mid - 1

            ! Element is equal to the middle element, return the position
            else 
                ! number 'element' found at array(mid)
                position = mid
                exit
            end if

            ! New middle element
            mid = lo + (hi - lo) / 2
            position = mid
        end do
    end subroutine

    subroutine checkVector(arr, length)
        ! Check if the vector is sorted

        implicit none
        real(8), dimension(:), allocatable, intent(in) :: arr
        integer, intent(in) :: length

        integer :: i, check

        ! Vector is sorted if check = 1
        check = 1

        ! Check if the vector is sorted by comparing each element with the next one
        do i = 1, length - 1

            ! If the vector is not sorted, print error message and exit
            if (arr(i) > arr(i + 1)) then
                print *, "ERROR: Vector is not sorted"
                print *, "arr(", i, ") = ", arr(i)
                print *, "arr(", i + 1, ") = ", arr(i + 1)
                check = 0
                exit
            end if

        end do

        ! Print message if vector is sorted
        if (check == 1) then
            print *, "Vector is sorted"
        end if

    end subroutine

    recursive subroutine insertMatrix(arr, matrix, length, element, ncolelements, valindex)
        ! Insert element into matrix and redistribute the elements if necessary

        implicit none
        real(8), dimension(:), allocatable, intent(inout) :: arr
        real(8), dimension(:,:), allocatable, intent(inout) :: matrix
        integer, intent(inout) :: length
        real(8), intent(in) :: element
        integer, dimension(:), allocatable, intent(inout) :: ncolelements
        real(8), dimension(:), allocatable, intent(inout) :: valindex
        
        integer :: inCol, inRow, i, m ,l
        real(8) :: tmp, mFactor, lFactor

        mFactor = 0.2
        lFactor = 0.2
        
        ! Check that matrix is ready to be used
        if (.not. allocated(matrix)) then
            m = length * mFactor
            l = length * lFactor
            call vectorToMatrix(arr, matrix, m, l, length, ncolelements, valindex)
        end if

        ! Find the column where the new value should be inserted
        call binarySearch(valindex, size(valindex), element, inCol)

        ! Element is the greatest of all currently saved...
        if (inCol > size(ncolelements)) then
            inCol = inCol - 1
        end if

      
        
        ! Check if the column is full
        if (ncolelements(inCol) > size(matrix, 1)) then
            ! Column is full, redistribute the elements
            call matrixToVector(arr, matrix, length, ncolelements)
            deallocate(matrix)
            deallocate(ncolelements)
            deallocate(valindex)
            call insertMatrix(arr, matrix, length, element, ncolelements, valindex)
            write(*, *) "Redistributed"
            
            go to 10
        end if

        ! Find the row where the new value should be inserted
        call binarySearch(matrix(:, inCol), ncolelements(inCol), element, inRow)

        ! Save new maximum value if necessary
        if (element - valindex(inCol) > 1.d-10) then
            valindex(inCol) = element
        end if
        ! Update number of elements in the column
        ncolelements(inCol) = ncolelements(inCol) + 1
        
        ! Shift all elements after the new value one position to the right
        matrix(inRow + 1 : ncolelements(inCol) + 1, inCol) = matrix(inRow : ncolelements(inCol), inCol)
        ! Insert the new value
        matrix(inRow, inCol) = element

        ! Update the length of the array
        length = length + 1

10  end subroutine

    subroutine insertVector(arr, length, element, position)
        ! Insert element into array at the given position

        implicit none
        real(8), dimension(:), allocatable, intent(inout) :: arr
        real(8), dimension(:), allocatable :: tmp
        real(8), intent(in) :: element
        integer, intent(in) :: position
        integer, intent(inout) :: length

        ! reallocate tmp
        call allocateArray(tmp, length + 1)

        ! copy elements before position
        tmp(1 : position - 1) = arr(1 : position - 1)
        ! copy elements after position
        tmp(position + 1 :) = arr(position:)
        ! insert element
        tmp(position) = element

        ! copy back to original array
        call allocateArray(arr, size(tmp))
        arr = tmp

        ! update length
        length = length + 1

    end subroutine

    subroutine matrixToVector(arr, matrix, length, ncolelements)
        ! Convert a matrix to an array

        implicit none
        real(8), dimension(:), allocatable, intent(out) :: arr
        real(8), dimension(:,:), allocatable, intent(in) :: matrix
        integer, intent(in) :: length
        integer, dimension(:), allocatable, intent(in) :: ncolelements
        
        integer :: i, j, counter
        
        call allocateArray(arr, length)
        counter = 1

        ! Assign the elements of the matrix to the array
        do i = 1, size(matrix, 2)
            do j = 1, ncolelements(i) 
                arr(counter) = matrix(j, i)
                counter = counter + 1
            end do 
        end do
    end subroutine

    subroutine merge(arr, tmp, lo, mid, hi)
        ! Array is sorted in the range [lo, mid] and [mid+1, hi]
        ! Merge the two sorted ranges into a single sorted range

        implicit none
        integer, intent(in) :: lo, mid, hi
        real(8), dimension(:), intent(inout) :: arr, tmp
        integer :: i, j, k, n

        ! Already correctly sorted, no need to merge
        if (arr(mid) <= arr(mid+1)) then 
            return 
        end if 

        i = lo
        j = mid + 1 
        n = size(arr)

        do k = lo, hi 
            ! Element of the first half is smaller or equal, take it
            if (j > hi .or. i <= mid .and. arr(i) <= arr(j)) then
                tmp(k) = arr(i)
                i = i + 1
            ! Element of the second half is smaller, take it
            else 
                tmp(k) = arr(j)
                j = j + 1
            end if 
        end do 
        
        ! Copy the sorted range back to the original array
        do k = 1,n 
            arr(k) = tmp(k)
        end do
    end subroutine

    recursive subroutine sort(arr, tmp, lo, hi)
        ! Sort the array arr in the range [lo, hi] using a top-down merge sort

        implicit none
        integer, intent(in) :: lo, hi
        real(8), dimension(:), intent(inout) :: arr, tmp
        integer :: mid

        ! Base case, only one element
        if (hi <= lo) then 
            return 
        end if

        ! Sort the two halves
        mid = lo + (hi - lo) / 2

        ! Sort the left half
        call sort(arr, tmp, lo, mid)
        ! Sort the right half
        call sort(arr, tmp, mid + 1, hi)

        ! Merge the two halves
        call merge(arr, tmp, lo, mid, hi)
    end subroutine

    subroutine vectorToMatrix(arr, matrix, m, l, length, ncolelements, valindex)
        ! Initialize the matrix with the given array
        ! Each column has (more or less) the same number of entries
        ! The number of entries of each column is stored in ncolelements
        ! The maximum value (last entry) of each column is stored in valindex

        implicit none
        real(8), dimension(:), intent(in) :: arr
        real(8), dimension(:,:), allocatable, intent(out) :: matrix
        integer, intent(in) :: length, m, l
        integer, dimension(:), allocatable, intent(out) :: ncolelements
        real(8), dimension(:), allocatable, intent(out) :: valindex
        
        integer :: nColEntries, rest, i, j, counter, n
        
        

        allocate(matrix(l, m))
        allocate(ncolelements(m))
        allocate(valindex(m))
        
        ! Determine the number of entries of each column 
        nColEntries = length / m
        ! Determine the number of columns with one more entry
        rest = mod(length, m)

        ! Fill the matrix
        counter = 1
        do i = 1, m

            ! If the column has one more entry, fill it with nColEntries + 1 entries
            if (i <= rest) then
                n = nColEntries + 1
            else 
                n = nColEntries
            end if

            ! Save the number of entries of each column
            ncolelements(i) = n

            ! Fill the column with n entries
            do j = 1, n
                matrix(j, i) = arr(counter)
                counter = counter + 1
            end do

            ! Save the maximum value of each column
            valindex(i) = arr(counter - 1)

        end do
    end subroutine

    subroutine vectorSorting(elementList, fileUnit, arr, length, insert)

        implicit none

        real(8), dimension(:), allocatable, intent(inout) :: elementList, arr
        integer, intent(inout) :: length, insert
        integer, intent(in) :: fileUnit

        real(8) :: element, startT, endT, arrayTime, preTime
        integer :: position, k


        do k=1,insert
            call cpu_time(startT)

            ! find position to insert element
            element = elementList(k)
            call binarySearch(arr, length, element, position)
            
            ! Insert element into the array version
            call insertVector(arr, length, element, position)
            !call checkVector(arr, length)
            call cpu_time(endT)
            arrayTime = endT - startT
            preTime = preTime + arrayTime
            write(fileUnit, '(*(G0.7,:,","))') length, preTime
        end do

    end subroutine

    subroutine matrixSorting(arr, matrix, initM, initL, length, elementList, insert, fileUnit)

        implicit none

        real(8), dimension(:,:), allocatable, intent(inout) :: matrix
        real(8), dimension(:), allocatable, intent(inout) :: arr, elementList
        integer, intent(in) :: initM, initL, fileUnit
        integer, intent(inout) :: insert, length

        real(8), dimension(:), allocatable :: valindex
        real(8) :: startT, endT, element, preTime, matrixTime
        integer, dimension(:), allocatable :: ncolelements
        integer :: k

        do k=1,insert
            call cpu_time(startT)
            element = elementList(k)
            if (k == 1) then
                call vectorToMatrix(arr, matrix, initM, initL, length, ncolelements, valindex)
            end if
            ! Insert element into the matrix version
            call insertMatrix(arr, matrix, length, element, ncolelements, valindex)
            !call matrixToVector(arr, matrix, length, ncolelements)
            !call checkVector(arr, length)
            call cpu_time(endT)
            matrixTime = endT - startT
            preTime = preTime + matrixTime
            write(fileUnit, '(*(G0.7,:,","))') length, preTime
        end do

    end subroutine
end module


program sorter
use sortingD 

    implicit none
    real(8), dimension(:,:), allocatable :: matrix
    real(8), dimension(:), allocatable :: origArr, arr, tmp, valindex, elementList
    integer, dimension(:), allocatable :: ncolelements
    real(8) :: element, arrayTime, matrixTime, preTime
    integer :: length, lo, position, k, testLength, insert, i, initM, initL
    real(8) :: startT, endT
    

    open(60, file="sumArr.csv")
    open(70, file="sumMat.csv")

    write(60, '(*(G0.7,:,","))') "Size of array", "Array Time"
    write(70, '(*(G0.7,:,","))') "Size of array", "Matrix Time"

    do i = 1, 4 ! only do 6 at pc
        lo = 1

        !--------------------------------------------------
        ! SIZE OF DATASET
        
        length = 10**i
        insert = 1000
        initM = 500
        initL = 500

        testLength = length
        preTime = 0.d0


        ! Allocate arrays
        call allocateArray(arr, length)
        call allocateArray(tmp, length)
        call allocateArray(origArr, length)
        call allocateArray(elementList, insert)

        ! Initialize array with random numbers between [0, 1]
        call random_number(arr)
        call random_number(elementList)

        ! sort array
        tmp = arr
        call sort(arr, tmp, lo, length)
        call checkVector(arr, length)

        origArr = arr
        
        ! TIME VECTOR --------------------------------------------------

        call vectorSorting(elementList, 60, arr, length, insert)
        
        

        write(60, *) "---"
        
        arr = origArr
        length = testLength
        preTime = 0.0
        

        ! TIME MATRIX --------------------------------------------------
        
        
        call matrixSorting(arr, matrix, initM, initL, length, elementList, insert)

        write(70, *) "---"

        deallocate(matrix)
        deallocate(origArr)
        deallocate(arr)
        deallocate(tmp)
        deallocate(elementList)
    end do
          
    close(60)
    close(70)

    

   
    



end program sorter
