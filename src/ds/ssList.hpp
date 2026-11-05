//desc.h 
#pragma once

#include <algorithm>
#include <cmath>
#include <cstddef>
#include <cstdint>
#include <vector>

namespace SolveSpace
{



// A simple list
template<class T>
class List {
    T *elem            = nullptr;
    int elemsAllocated = 0;

public:
    int  n = 0;

    bool IsEmpty() const { return n == 0; }

    void ReserveMore(int howMuch) {
        if(n + howMuch > elemsAllocated) {
            elemsAllocated = n + howMuch;
            T *newElem = (T *)::operator new[]((size_t)elemsAllocated*sizeof(T));
            for(int i = 0; i < n; i++) {
                new(&newElem[i]) T(std::move(elem[i]));
                elem[i].~T();
            }
            ::operator delete[](elem);
            elem = newElem;
        }
    }

    void AllocForOneMore() {
        if(n >= elemsAllocated) {
            ReserveMore((elemsAllocated + 32)*2 - n);
        }
    }

    void Add(const T *t) {
        AllocForOneMore();
        new(&elem[n++]) T(*t);
    }

    void AddToBeginning(const T *t) {
        AllocForOneMore();
        new(&elem[n]) T();
        std::move_backward(elem, elem + 1, elem + n + 1);
        elem[0] = *t;
        n++;
    }

    T *First() {
        return IsEmpty() ? nullptr : &(elem[0]);
    }
    const T *First() const {
        return IsEmpty() ? nullptr : &(elem[0]);
    }

    T *Last() { return IsEmpty() ? nullptr : &(elem[n - 1]); }
    const T *Last() const { return IsEmpty() ? nullptr : &(elem[n - 1]); }

    T *NextAfter(T *prev) {
        if(IsEmpty() || !prev) return NULL;
        if(prev - First() == (n - 1)) return NULL;
        return prev + 1;
    }
    const T *NextAfter(const T *prev) const {
        if(IsEmpty() || !prev) return NULL;
        if(prev - First() == (n - 1)) return NULL;
        return prev + 1;
    }

    T &Get(size_t i) { return elem[i]; }
    T const &Get(size_t i) const { return elem[i]; }
    T &operator[](size_t i) { return Get(i); }
    T const &operator[](size_t i) const { return Get(i); }

    T *begin() { return IsEmpty() ? nullptr : &elem[0]; }
    T *end() { return IsEmpty() ? nullptr : &elem[n]; }
    const T *begin() const { return IsEmpty() ? nullptr : &elem[0]; }
    const T *end() const { return IsEmpty() ? nullptr : &elem[n]; }
    const T *cbegin() const { return begin(); }
    const T *cend() const { return end(); }

    void ClearTags() {
        for(auto & elt : *this) {
            elt.tag = 0;
        }
    }

    void Clear() {
        for(int i = 0; i < n; i++)
            elem[i].~T();
        if(elem) ::operator delete[](elem);
        elem = NULL;
        n = elemsAllocated = 0;
    }

    void RemoveTagged() {
        auto newEnd = std::remove_if(this->begin(), this->end(), [](T &t) {
            if(t.tag) {
                return true;
            }
            return false;
        });
        auto oldEnd = this->end();
        n = newEnd - begin();
        if (newEnd != nullptr && oldEnd != nullptr) {
            while(newEnd != oldEnd) {
                newEnd->~T();
                ++newEnd;
            }
        }
        // and elemsAllocated is untouched, because we didn't resize
    }

    void RemoveLast(int cnt) {
        ssassert(n >= cnt, "Removing more elements than the list contains");
        for(int i = n - cnt; i < n; i++)
            elem[i].~T();
        n -= cnt;
        // and elemsAllocated is untouched, same as in RemoveTagged
    }

    void Reverse() {
        int i;
        for(i = 0; i < (n/2); i++) {
            swap(elem[i], elem[(n-1)-i]);
        }
    }
};

}//ss