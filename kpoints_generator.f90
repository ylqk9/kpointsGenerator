module ModUtlt
    implicit none
    private
    real(8), parameter      :: TooClose = 1D-5  ! when two atoms have a distance smaller than this, they are treated as one atom.
    
    public :: PushBack
    public :: GetRotationMatrix
    public :: IdenticalGeometry
    public :: Norm
    interface PushBack
        module procedure PushBackNumber
        module procedure PushBack3Vec
    endinterface
    
    contains
    
    subroutine PushBackNumber(Array, Num)
        real(8), allocatable, intent(inout) :: Array(:)
        real(8), intent(in)                 :: Num
        real(8), allocatable                :: ArrayTmp(:)
        integer                             :: Asize
        
        if(allocated(Array))then
            allocate(ArrayTmp, source = Array)
            Asize = size(Array, dim = 1)
            deallocate(Array)
            allocate(Array(Asize + 1))
            Array(1:Asize) = ArrayTmp
            Array(Asize + 1) = Num
            deallocate(ArrayTmp)
        else
            allocate(Array(1))
            Array(1) = Num
        endif
    endsubroutine PushBackNumber
    
    subroutine PushBack3Vec(Array, Vec)
        real(8), allocatable, intent(inout) :: Array(:,:)
        real(8), intent(in)                 :: Vec(3)
        real(8), allocatable                :: ArrayTmp(:,:)
        integer                             :: Asize
        
        if(allocated(Array))then
            allocate(ArrayTmp, source = Array)
            Asize = size(Array, dim = 2)
            deallocate(Array)
            allocate(Array(3, Asize + 1))
            Array(:,1:Asize) = ArrayTmp
            Array(:,Asize + 1) = Vec
            deallocate(ArrayTmp)
        else
            allocate(Array(3,1))
            Array(:,1) = Vec
        endif
    endsubroutine PushBack3Vec
    
    logical function IdenticalGeometry(geo1, geo2)
        real(8), intent(in) :: geo1(:,:), geo2(:,:)
        integer             :: i, j, Ngeo
        real(8)             :: Distance, rr(3)
        integer             :: Nmatch
        
        Ngeo = size(geo1, dim = 2)
        Nmatch = 0
        do i = 1, ngeo
            do j = 1, ngeo
                rr = geo1(:,i) - geo2(:,j)
                Distance = sqrt(dot_product(rr, rr))
                if(Distance < TooClose) then
                    Nmatch = Nmatch + 1
                    exit
                endif
            enddo
        enddo
        IdenticalGeometry = (Nmatch == Ngeo)
    endfunction IdenticalGeometry
    
    function GetRotationMatrix(v1, v2)
        real(8), intent(in) :: v1(3), v2(3)
        real(8)             :: GetRotationMatrix(3,3)
        real(8)             :: vc(3)    ! cross product
        real(8)             :: vx(3,3)  ! skew-symmetric cross-product
        real(8)             :: v1n(3), v2n(3) ! Normalized v1 and v2
        real(8)             :: normvc, norm1, norm2, dot12, unit(3,3)
        
        norm1 = sqrt(dot_product(v1, v1))
        norm2 = sqrt(dot_product(v2, v2))
        v1n = v1/norm1
        v2n = v2/norm2
        vc = cross(v1n, v2n)
        normvc = sqrt(dot_product(vc, vc))
        if(normvc > 1D-10) then
            dot12 = dot_product(v1n, v2n)
            vx = 0.D0
            vx(1,2) = -vc(3)
            vx(1,3) = vc(2)
            vx(2,1) = vc(3)
            vx(2,3) = -vc(1)
            vx(3,1) = -vc(2)
            vx(3,2) = vc(1)
            unit = 0.D0
            unit(1,1) = 1.D0
            unit(2,2) = 1.D0
            unit(3,3) = 1.D0
            GetRotationMatrix = unit + vx + matmul(vx,vx)*(1 - dot12)/normvc/normvc
            GetRotationMatrix = GetRotationMatrix*norm2/norm1
        else ! normvc = 0 means same direction
            if(v1(1) /= 0.D0)then
                GetRotationMatrix = unit*v2(1)/v1(1)
            elseif(v1(2)/= 0.D0)then
                GetRotationMatrix = unit*v2(2)/v1(2)
            else
                GetRotationMatrix = unit*v2(3)/v1(3)
            endif
        endif   
    contains
        function cross(a, b)
            real(8), intent(in) :: a(3), b(3)
            real(8)             :: cross(3)

            cross(1) = a(2) * b(3) - a(3) * b(2)
            cross(2) = a(3) * b(1) - a(1) * b(3)
            cross(3) = a(1) * b(2) - a(2) * b(1)
        endfunction cross
    endfunction GetRotationMatrix
    
    real(8) function Norm(v)
        real(8), intent(in) :: v(3)
        
        Norm = sqrt(dot_product(v, v))
    endfunction Norm
endmodule ModUtlt

program kpoints_generator
    use ModUtlt
    implicit none
    real(8), parameter      :: Pi = 3.14159265358979323846D0
    logical                 :: IfSymmetry, InvalidInput, FoundInputFile, FoundGeometryFile
    character(len = 256)    :: FileGeometry, FileInput, Command, ExtraInput, buf
    real(8)                 :: L(3)             ! Supercell Dimension
    real(8)                 :: B(3)             ! reciprocal vector (this doesn't have to be a vector for a rectangular supercell
    real(8)                 :: R(3,3), MaxXYZ, MinXYZ
    integer                 :: N(3), x, y, z, i, j, m, Natoms, Ndiff
    real(8), allocatable    :: k(:,:), w(:), Reducedk(:,:), Reducedw(:)
    real(8), allocatable    :: Geo(:,:), RGeo1(:,:), RGeo2(:,:)
    real(8)                 :: Front(3) = [1.D0, 0.D0, 0.D0], Shift(3), Phase(3)
    character(len = 2),&
        allocatable         :: element(:)
    
    call get_command_argument(0, Command)
    call get_command_argument(1, FileInput)
    call get_command_argument(2, ExtraInput)
    if(len_trim(FileInput) == 0 .or. len_trim(ExtraInput) /= 0) then
        write(6, '(A)') "Instruction:"
        write(6, '(A,X,A)') trim(Command), "your_input_file_name"
        write(6, '(A)') "For the first time, use a new name as the argument."
        write(6, '(A)') "Then you will have a input file with instructions to fill"
        write(6, '(A)') "After the file is filled rerun the same command."
        stop
    endif
    inquire(file = trim(FileInput), exist = FoundInputFile)
    if(.not. FoundInputFile) then
        open(unit = 101, file = trim(FileInput))
        write(101, '(A)') "Please give the size of the supercell: Lx Ly Lz"
        write(101, '(A)') "If 2D/1D set L to 0"
        write(101, *)
        write(101, '(A)') "Please provide the size of the k space: Nx Ny Nz"
        write(101, '(A)') "If 2D/1D set N to 1"
        write(101, *)
        write(101, '(A)') "Set Gamma point to: x y z"
        write(101, *)
        write(101, '(A)') "Reduce the k vectors by symmetry? T/F"
        write(101, *)
        write(101, '(A)') "Please provide the geometry file in xyz format."
        write(101, *)
        close(101)
        stop
    else
        open(unit = 101, file = trim(FileInput))
        InvalidInput = .true.
        read(101, *)
        read(101, *)
        read(101, *) L(1), L(2), L(3)
        if(all(L >= 0.D0)) InvalidInput = .false.
        if(InvalidInput) write(6, '(A)') "Find invalid input in the 1st row."
        read(101, *)
        read(101, *)
        read(101, *)N(1), N(2), N(3)
        if(all(N > 0)) InvalidInput = .false.
        if(InvalidInput) write(6, '(A)') "Find invalid input in the 2nd row."
        read(101, *)
        read(101, *)Shift(1), Shift(2), Shift(3)
        if(any(Shift < 0) .or. any(Shift > L)) InvalidInput = .false.
        if(InvalidInput) write(6, '(A)') "Find invalid input in the 3rd row."
        IfSymmetry = .true.
        read(101, *)
        read(101, '(L)')IfSymmetry
        read(101, *)
        read(101, '(A)')FileGeometry
        inquire(file = trim(FileGeometry), exist = FoundGeometryFile)
        if(.not. FoundGeometryFile) write(6, '(A)') "Find invalid input: geometry file not exist."
        if(InvalidInput) stop
        close(101)
    endif
    write(6, '(A)') "========================================"
    write(6, '(A)') "==========K points generator============"
    write(6, '(A)') "========================================"
    write(6, '(A)') "This is programmed for the rectangular supercell only."
    write(6, '(A30,F8.3,A,F8.3,A,F8.3)') "Supercell is:", L(1), " x ", L(2), " x ", L(3)
    write(6, '(A30,I8,A,I8,A,I8)') "k mesh is:", N(1), " x ", N(2), " x ", N(3)
    write(6, '(A30,L)') "Reduce k by symmetry:", IfSymmetry
    write(6, '(A30,A)') "The geometry file is:", trim(FileGeometry)
    allocate(k(3,product(N)))
    allocate(w(product(N)))
    B = 2*Pi/L
    Shift = B*Shift/L
    if(L(1) == 0.D0) then
        B(1) = 0.D0
        Shift(1) = 0.D0
    endif
    if(L(2) == 0.D0) then
        B(2) = 0.D0
        Shift(2) = 0.D0
    endif
    if(L(3) == 0.D0) then
        B(3) = 0.D0
        Shift(3) = 0.D0
    endif
    !if(mod(N(1), 2) == 0) then ! even
    !    k(1,1) = B(1)/2*(-1 + 1.D0/N(1))
    !else ! odd
    !    k(1,1) = B(1)/N(1)*(-N(1)/2)
    !endif
    !if(mod(N(2), 2) == 0) then ! even
    !    k(2,1) = B(2)/2*(-1 + 1.D0/N(2))
    !else ! odd
    !    k(2,1) = B(2)/N(2)*(-N(2)/2)
    !endif
    !if(mod(N(3), 2) == 0) then ! even
    !    k(3,1) = B(3)/2*(-1 + 1.D0/N(3))
    !else ! odd
    !    k(3,1) = B(3)/N(3)*(-N(3)/2)
    !endif
    !k(:,1) = k(:,1) + Shift
    !i = 1
    !do x = 0, N(1) - 1
    !    do y = 0, N(2) - 1
    !        do z = 0, N(3) - 1
    !            k(1,i) = k(1,1) + x*B(1)/N(1)
    !            k(2,i) = k(2,1) + y*B(2)/N(2)
    !            k(3,i) = k(3,1) + z*B(3)/N(3)
    !            i = i + 1
    !        enddo
    !    enddo
    !enddo
    i = 1
    do x = 0, N(1) - 1
        do y = 0, N(2) - 1
            do z = 0, N(3) - 1
                Phase(1) = Shift(1) + dble(x)/N(1)
                Phase(2) = Shift(2) + dble(y)/N(2)
                Phase(3) = Shift(3) + dble(z)/N(3)
                if(Phase(1) > 0.5D0) Phase(1) = Phase(1) - 1.D0
                if(Phase(2) > 0.5D0) Phase(2) = Phase(2) - 1.D0
                if(Phase(3) > 0.5D0) Phase(3) = Phase(3) - 1.D0
                print *, Phase
                k(:,i) = Phase*2.D0*B/N
                i = i + 1
            enddo
        enddo
    enddo
    w = 1.D0/product(N)
    if(IfSymmetry) then
        open(unit = 102, file = trim(FileGeometry))
        read(102, *)Natoms
        allocate(Geo(3,Natoms))
        allocate(element(Natoms))
        allocate(RGeo1(3,Natoms))
        allocate(RGeo2(3,Natoms))
        do i = 1, Natoms
            read(102, *)element(i), Geo(:,i)
        enddo
        ! shift the geometry
        MaxXYZ = maxval(Geo(1,:))
        MinXYZ = minval(Geo(1,:))
        Geo(1,:) = Geo(1,:) - (MaxXYZ + MinXYZ)/2.D0
        MaxXYZ = maxval(Geo(2,:))
        MinXYZ = minval(Geo(2,:))
        Geo(2,:) = Geo(2,:) - (MaxXYZ + MinXYZ)/2.D0
        MaxXYZ = maxval(Geo(3,:))
        MinXYZ = minval(Geo(3,:))
        Geo(3,:) = Geo(3,:) - (MaxXYZ + MinXYZ)/2.D0
        close(102)
        call PushBack(Reducedk, k(:,1))
        call PushBack(Reducedw, w(1))
        print *, k(:,1)
        do i = 2, product(N)
            Ndiff = 0
            print *, k(:,i)
            if(Norm(k(:,i)) < 1D-5) then
                call PushBack(Reducedk, k(:,i))
                call PushBack(Reducedw, w(i))
                cycle
            endif
            do j = 1, i - 1
                if(abs(Norm(k(:,i)) - Norm(k(:,j))) > 1D-5) then
                    Ndiff = Ndiff + 1
                    cycle
                endif
                R = GetRotationMatrix(Front,k(:,j))
                !print *, "original"
                !print *, k(:,i)
                !print *, "rotational"
                !print *, matmul(R,k(:,j))
                do m = 1, size(Geo, dim = 2)
                    RGeo1(:,m) = matmul(R,Geo(:,m))
                enddo
                R = GetRotationMatrix(Front,k(:,i))
                do m = 1, size(Geo, dim = 2)
                    RGeo2(:,m) = matmul(R,Geo(:,m))
                enddo
                do m = 1, size(Geo, dim = 2)
                    write(11,'(A6,3F8.3)') "C", RGeo1(1,m), RGeo1(2,m), RGeo1(3,m)
                enddo
                do m = 1, size(Geo, dim = 2)
                    write(12,'(A6,3F8.3)') "C", RGeo2(1,m), RGeo2(2,m), RGeo2(3,m)
                enddo
                if(.not. IdenticalGeometry(RGeo1, RGeo2)) then
                    RGeo2(1,:) = RGeo2(1,:)*k(1,j)/k(1,i)
                    RGeo2(2,:) = RGeo2(2,:)*k(2,j)/k(2,i)
                    RGeo2(3,:) = RGeo2(3,:)*k(3,j)/k(3,i)
                    if(.not. IdenticalGeometry(RGeo1, RGeo2))Ndiff = Ndiff + 1
                endif
            enddo
            if(Ndiff == i - 1) then
                print *, "added"
                call PushBack(Reducedk, k(:,i))
                call PushBack(Reducedw, w(i))
            endif
        enddo
    else
        do i = 1, size(k, dim = 2)
            write(999, '(4F25.15)')k(:,i), w(i)
        enddo
    endif
    read(5, *)
    do i = 1, size(Reducedk, dim = 2)
        print *, Reducedk(:,i)
    enddo
    read(5, *)
endprogram kpoints_generator
