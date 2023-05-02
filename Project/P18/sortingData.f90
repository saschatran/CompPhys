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
        arr = 0.0

    end subroutine

    subroutine binarySearch(arr, length, element, position, count)
        ! implements a binary search algorithm to find the position of 'element' in 'arr'
        ! 'arr' has to be sorted in ascending order
        ! returns the position of 'element' in 'arr' in 'position'
        ! if 'element' is not in 'arr', the position of the next larger element is returned

        implicit none 
        integer, intent(in) :: length
        real(8), dimension(length), intent(in) :: arr
        real(8), intent(in) :: element
        integer, intent(inout) :: position, count
        
        integer :: lo, mid, hi 
        
        ! initialize
        lo = 1 
        hi = length
        mid = lo + (hi - lo) / 2

        ! check if the element is in the middle of the array
        do while (lo <= hi)
            count = count + 1
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

    recursive subroutine insertMatrix(arr, matrix, length, element, ncolelements, valindex, count)
        ! Insert element into matrix and redistribute the elements if necessary

        implicit none
        real(8), dimension(:), allocatable, intent(inout) :: arr
        real(8), dimension(:,:), allocatable, intent(inout) :: matrix
        integer, intent(inout) :: length
        real(8), intent(in) :: element
        integer, dimension(:), allocatable, intent(inout) :: ncolelements
        real(8), dimension(:), allocatable, intent(inout) :: valindex
        integer, intent(out) :: count

        integer :: inCol, inRow, i, m ,l, nonempty
        real(8) :: tmp, mFactor, lFactor


        inCol = 0
        
        mFactor = 0.3
        lFactor = 0.3
        
        ! Check that matrix is ready to be used
        if (.not. allocated(matrix)) then
            m = length * mFactor
            l = length * lFactor
            call vectorToMatrix(arr, matrix, m, l, length, ncolelements, valindex)
        end if

        

        ! Find the column where the new value should be inserted
        nonempty = findloc(valindex, 0, dim=1)
        if (nonempty == 0) then
            nonempty = size(valindex)
        end if
        call binarySearch(valindex(1 : nonempty - 1), nonempty - 1, element, inCol, count)

        ! Element is the greatest of all currently saved...
        if (inCol > size(ncolelements)) then
            inCol = inCol - 1
        end if

      
        
        ! Check if the column is full
        if (ncolelements(inCol) >= size(matrix, 1)) then
            ! Column is full, redistribute the elements
            call matrixToVector(arr, matrix, length, ncolelements)
            deallocate(matrix)
            deallocate(ncolelements)
            deallocate(valindex)
            call insertMatrix(arr, matrix, length, element, ncolelements, valindex, count)
            write(*, *) "Redistributed"
            
            go to 10
        end if

        
        ! Find the row where the new value should be inserted
        call binarySearch(matrix(:, inCol), ncolelements(inCol), element, inRow, count)
        ! If the column is empty, the row will be 0 therefore set inRow to 1
        if (inRow == 0) then
            inRow = 1
        end if

        ! Save new maximum value if necessary
        if (element - valindex(inCol) > 1.d-10) then
            valindex(inCol) = element
        end if
        
        ! Shift all elements after the new value one position to the right
        
        do i = inRow + 1, ncolelements(inCol) + 1
            matrix(i, inCol) = matrix(i - 1, inCol)
        end do
        
        ! Insert the new value
        matrix(inRow, inCol) = element
        ! Update number of elements in the column
        ncolelements(inCol) = ncolelements(inCol) + 1

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
            if (ncolelements(i) == 0) then
                go to 20
            end if
            do j = 1, ncolelements(i) 
                arr(counter) = matrix(j, i)
                counter = counter + 1
            end do 
        end do
    20 end subroutine

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
        
        if (allocated(ncolelements)) then
            deallocate(ncolelements)
        end if

        allocate(matrix(l, m))
        allocate(ncolelements(m))
        allocate(valindex(m))

        matrix = 0.0
        ncolelements = 0
        valindex = 0
        
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
            
            if (n > 0) then
                ! Save the maximum value of each column
                valindex(i) = arr(counter - 1)
            end if

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
    integer :: length, lo, position, k, testLength, insert, i, initM, initL, j, stepI, count
    real(8) :: startT, endT
    
    count = 0

    !--------------------------------------------------
    ! SIZE OF DATASET
    !--------------------------------------------------

    open(60, file="sumArr.csv")
    open(70, file="sumMat.csv")

    write(60, '(*(G0.7,:,","))') "Size of array", "Array Time"
    write(70, '(*(G0.7,:,","))') "Size of array", "Matrix Time"

    do i = 1, 5
        lo = 1
        
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

        do k=1,insert
            call cpu_time(startT)

            ! find position to insert element
            call binarySearch(arr, length, element, position, count)
            
            element = elementList(k)
            ! Insert element into the array version
            call insertVector(arr, length, element, position)
            !call checkVector(arr, length)
            call cpu_time(endT)
            arrayTime = endT - startT
            preTime = preTime + arrayTime
            write(60, '(*(G0.7,:,","))') length, preTime
        end do

        write(60, *) "---"
        
        arr = origArr
        length = testLength
        preTime = 0.0
        

        ! TIME MATRIX --------------------------------------------------
        
        
        do k=1,insert
            call cpu_time(startT)
            element = elementList(k)
            if (k == 1) then
                call vectorToMatrix(arr, matrix, initM, initL, length, ncolelements, valindex)
            end if
            
            ! Insert element into the matrix version
            call insertMatrix(arr, matrix, length, element, ncolelements, valindex, count)
            call matrixToVector(tmp, matrix, length, ncolelements)
            call checkVector(tmp, length)
            
            call cpu_time(endT)
            matrixTime = endT - startT
            preTime = preTime + matrixTime
            write(70, '(*(G0.7,:,","))') length, preTime
        end do

        write(70, *) "---"

       

        deallocate(matrix)
        deallocate(origArr)
        deallocate(arr)
        deallocate(tmp)
        deallocate(valindex)
        deallocate(elementList)
        deallocate(ncolelements)
    end do
          
    close(60)
    close(70)







    ! CHOICE OF THE SIDE LENGTHs




    !--------------------------------------------------
    ! DIFFERENT MATRIX SIZES Many changes in M
    !--------------------------------------------------

    
    open(80, file="mat.csv")

    write(80, '(*(G0.7,:,","))') "initM", "initL", "Matrix Time", "Count"

    lo = 1    
    length = 1000
    insert = 1000

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

    deallocate(arr)
    deallocate(tmp)
    
    do j = 100, 2000, 100
        stepI = 10
        i = 10
        do while (i < 2001)
            i = i + stepI
            
            if (i == 100) then
                stepI = 100
            end if
            
            
            initM = i
            initL = j

            call allocateArray(arr, length)
            
            arr = origArr
            length = testLength
            preTime = 0.0
            
            count = 0

            ! TIME MATRIX --------------------------------------------------
            
            
            do k=1,insert
                call cpu_time(startT)
                element = elementList(k)
                if (k == 1) then
                    call vectorToMatrix(arr, matrix, initM, initL, length, ncolelements, valindex)
                end if
                ! Insert element into the matrix version
                call insertMatrix(arr, matrix, length, element, ncolelements, valindex, count)
                call matrixToVector(arr, matrix, length, ncolelements)
                call checkVector(arr, length)
                call cpu_time(endT)
                matrixTime = endT - startT
                preTime = preTime + matrixTime
            end do

            write(80, '(*(G0.7,:,","))') initM, initL, preTime, count


            deallocate(matrix)
            deallocate(arr)
            deallocate(valindex)
            deallocate(ncolelements)
        end do
    end do

  

    
    !--------------------------------------------------
    ! DIFFERENT MATRIX SIZES Many changes in M
    !--------------------------------------------------


    
    do j = 100, 2000, 100
        stepI = 10
        i = 10
        do while (i < 2001)
            i = i + stepI
            
            if (i == 100) then
                stepI = 100
            end if
            
            initM = j
            initL = i

            call allocateArray(arr, length)
            
            arr = origArr
            length = testLength
            preTime = 0.0
            count = 0
            

            ! TIME MATRIX --------------------------------------------------
            
            
            do k=1,insert
                call cpu_time(startT)
                element = elementList(k)
                if (k == 1) then
                    call vectorToMatrix(arr, matrix, initM, initL, length, ncolelements, valindex)
                end if
                ! Insert element into the matrix version
                call insertMatrix(arr, matrix, length, element, ncolelements, valindex, count)
                call matrixToVector(arr, matrix, length, ncolelements)
                call checkVector(arr, length)
                call cpu_time(endT)
                matrixTime = endT - startT
                preTime = preTime + matrixTime
            end do

            write(80, '(*(G0.7,:,","))') initM, initL, preTime, count


            deallocate(matrix)
            deallocate(arr)
            deallocate(valindex)
            deallocate(ncolelements)
        end do
    end do
    close(80)

    deallocate(elementList)
          
    


   
    



end program sorter
