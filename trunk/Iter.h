#ifndef ITER_H
#define ITER_H

template<class T, class Q> class Iter {
public:
    Iter (const T* p_vec, int pos) : _pos( pos ) , _p_vec( p_vec ) { }
    bool operator!= (const Iter& other) const {
        return _pos != other._pos;
    }
    Q operator* () const {
        return _p_vec->get( _pos );
    }
    const Iter& operator++ () {
        ++_pos;
        return *this;
    }

private:
    int _pos;
    const T* _p_vec;
};

#endif
