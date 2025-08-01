#pragma once
#include <petscvec.h>
#include <stdexcept>

class MyPetscVec {
public:
    MyPetscVec(MPI_Comm comm, PetscInt size)
        : size_(size), comm_(comm) {
        PetscErrorCode ierr;
        ierr = VecCreate(comm_, &vec_);                     CHKERRABORT(comm_, ierr);
        ierr = VecSetSizes(vec_, PETSC_DECIDE, size_);      CHKERRABORT(comm_, ierr);
        ierr = VecSetFromOptions(vec_);                     CHKERRABORT(comm_, ierr);
    }

    ~MyPetscVec() {
        if (vec_ != nullptr) {
            VecDestroy(&vec_);
        }
    }

    // Move constructor
    MyPetscVec(MyPetscVec&& other) noexcept
        : vec_(other.vec_), size_(other.size_), comm_(other.comm_) {
        other.vec_ = nullptr;
    }

    // Move assignment
    MyPetscVec& operator=(MyPetscVec&& other) noexcept {
        if (this != &other) {
            if (vec_ != nullptr) {
                VecDestroy(&vec_);
            }
            vec_ = other.vec_;
            size_ = other.size_;
            comm_ = other.comm_;
            other.vec_ = nullptr;
        }
        return *this;
    }

    // Deleted copy constructor/assignment
    MyPetscVec(const MyPetscVec&) = delete;
    MyPetscVec& operator=(const MyPetscVec&) = delete;

    // Set value at index
    inline void setValue(PetscInt i, PetscScalar value, InsertMode mode = INSERT_VALUES) {
        VecSetValue(vec_, i, value, mode);
    }

    // Finish inserting values
    inline void assemble() {
        VecAssemblyBegin(vec_);
        VecAssemblyEnd(vec_);
    }

    // Return raw PETSc Vec
    inline Vec& raw() {
        return vec_;
    }

    // Return global size
    inline PetscInt size() const {
        return size_;
    }

    // Zero the vector
    inline void zero() {
        VecSet(vec_, 0.0);
    }

    // In-place addition: this += other
    MyPetscVec& operator+=(const MyPetscVec& other) {
        check_size_match(*this, other);
        VecAXPY(vec_, 1.0, other.vec_);
        return *this;
    }

    // Addition: result = a + b
    friend MyPetscVec operator+(const MyPetscVec& a, const MyPetscVec& b) {
        check_size_match(a, b);
        MyPetscVec result(a.comm_, a.size_);
        VecCopy(a.vec_, result.vec_);
        VecAXPY(result.vec_, 1.0, b.vec_);
        return result;
    }

    // Subtraction: result = a - b
    friend MyPetscVec operator-(const MyPetscVec& a, const MyPetscVec& b) {
        check_size_match(a, b);
        MyPetscVec result(a.comm_, a.size_);
        VecWAXPY(result.vec_, -1.0, b.vec_, a.vec_);
        return result;
    }

    // Unary minus: -a
    MyPetscVec operator-() const {
        return (-1.0) * (*this);
    }

    // Scalar multiplication: result = alpha * a
    friend MyPetscVec operator*(PetscScalar alpha, const MyPetscVec& a) {
        MyPetscVec result(a.comm_, a.size_);
        VecCopy(a.vec_, result.vec_);
        VecScale(result.vec_, alpha);
        return result;
    }

    // In-place scaling: a *= alpha
    MyPetscVec& operator*=(PetscScalar alpha) {
        VecScale(vec_, alpha);
        return *this;
    }

    // Dot product: a.dot(b)
    PetscScalar dot(const MyPetscVec& other) const {
        check_size_match(*this, other);
        PetscScalar result;
        VecDot(vec_, other.vec_, &result);
        return result;
    }

    // Norm of vector: a.norm()
    PetscScalar norm(NormType type = NORM_2) const {
        PetscReal result;
        VecNorm(vec_, type, &result);
        return result;
    }

private:
    Vec vec_;
    PetscInt size_;
    MPI_Comm comm_;

    static void check_size_match(const MyPetscVec& a, const MyPetscVec& b) {
        if (a.size_ != b.size_) {
            throw std::invalid_argument("MyPetscVec: size mismatch in vector operation.");
        }
    }
};
