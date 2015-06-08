#ifndef MATRIXCRS_H
#define MATRIXCRS_H

#include "dtypes.h"
#include "Serializable.h"

#include <math.h>
#include <vector>
#include <iostream>
#include <iomanip>
#include <algorithm>
#include <iterator>
#include <typeinfo>
#include <memory>

class MatrixCOO;

class MatrixCRS : public Serializable {
public:
    MatrixCRS(size_t nRows, size_t nCols, size_t size);
    MatrixCRS(std::vector<unsigned char> &bytes);
    MatrixCRS();
    void add(didx rIncrement, didx c, dreal v);
    void multiplyScalar(dreal &&scalar);
    void order();
    void sumup();
    void print();
    std::shared_ptr<MatrixCOO> toCOO();
    std::vector<unsigned char> serialize();
    void deserialize(std::vector<unsigned char> &bytes);
    size_t getNumRows();
    size_t getNumCols();
    std::vector<size_t> getDimensions();
    void setDimensions(std::vector<size_t> &&dimensions);

private:
    std::vector<didx> ptr;
    std::vector<didx> col;
    std::vector<dreal> val;
    size_t nRows;
    size_t nCols;
    bool ordered;
    bool sumuped;

    template<class T>
    void serializeItem(std::vector<unsigned char>::iterator &pbytes, std::vector<T> &data) {
        size_t szvector = data.size();

        size_t sz = sizeof(szvector);
        auto pchar = reinterpret_cast<const unsigned char*>(&szvector);
        std::copy(pchar, pchar+sz, pbytes);
        pbytes += sz;

        size_t szitem = sizeof(data[0]);

        sz = sizeof(szitem);
        pchar = reinterpret_cast<const unsigned char*>(&szitem);
        std::copy(pchar, pchar+sz, pbytes);
        pbytes += sz;

        size_t typeitem = typeid(data[0]).hash_code();

        sz = sizeof(typeitem);
        pchar = reinterpret_cast<const unsigned char*>(&typeitem);
        std::copy(pchar, pchar+sz, pbytes);
        pbytes += sz;

        sz = szvector*szitem;
        pchar = reinterpret_cast<const unsigned char*>(&data[0]);
        std::copy(pchar, pchar+sz, pbytes);
        pbytes += sz;
    }

    template<class T>
    void serializeItem(std::vector<unsigned char>::iterator &pbytes, T &data) {
        size_t sz = sizeof(data);
        auto pchar = reinterpret_cast<const unsigned char*>(&data);
        std::copy(pchar, pchar+sz, pbytes);
        pbytes += sz;
    }

    template<class T>
    void deserializeItem(std::vector<unsigned char>::iterator &pbytes, std::vector<T> &data) {
        size_t szvector;

        size_t sz = sizeof(szvector);
        auto pchar = reinterpret_cast<unsigned char*>(&szvector);
        std::copy(pbytes, pbytes+sz, pchar);
        pbytes += sz;

        size_t szitem;

        sz = sizeof(szitem);
        pchar = reinterpret_cast<unsigned char*>(&szitem);
        std::copy(pbytes, pbytes+sz, pchar);
        pbytes += sz;

        size_t typeitem;

        sz = sizeof(typeitem);
        pchar = reinterpret_cast<unsigned char*>(&typeitem);
        std::copy(pbytes, pbytes+sz, pchar);
        pbytes += sz;

        data.resize(szvector);

        sz = szvector*szitem;
        if (typeitem == typeid(float).hash_code()) {
            auto pitem = reinterpret_cast<float*>(&(*pbytes));
            auto pdata = &data[0];
            std::copy(pitem, pitem+szvector, pdata);
        } else if (typeitem == typeid(double).hash_code()) {
            auto pitem = reinterpret_cast<double*>(&(*pbytes));
            auto pdata = &data[0];
            std::copy(pitem, pitem+szvector, pdata);
        } else if (typeitem == typeid(int).hash_code()) {
            auto pitem = reinterpret_cast<int*>(&(*pbytes));
            auto pdata = &data[0];
            std::copy(pitem, pitem+szvector, pdata);
        } else if (typeitem == typeid(long).hash_code()) {
            auto pitem = reinterpret_cast<long*>(&(*pbytes));
            auto pdata = &data[0];
            std::copy(pitem, pitem+szvector, pdata);
        } else if (typeitem == typeid(unsigned int).hash_code()) {
            auto pitem = reinterpret_cast<unsigned int*>(&(*pbytes));
            auto pdata = &data[0];
            std::copy(pitem, pitem+szvector, pdata);
        } else if (typeitem == typeid(unsigned long).hash_code()) {
            auto pitem = reinterpret_cast<unsigned long*>(&(*pbytes));
            auto pdata = &data[0];
            std::copy(pitem, pitem+szvector, pdata);
        } else if (typeitem == typeid(size_t).hash_code()) {
            auto pitem = reinterpret_cast<size_t*>(&(*pbytes));
            auto pdata = &data[0];
            std::copy(pitem, pitem+szvector, pdata);
        } else {
            auto pitem = pbytes;
            auto pdata = reinterpret_cast<unsigned char*>(&data[0]);
            std::copy(pitem, pitem+szvector, pdata);
        }
        pbytes += sz;
    }

    template<class T>
    void deserializeItem(std::vector<unsigned char>::iterator &pbytes, T &data) {
        size_t sz = sizeof(data);
        auto pchar = reinterpret_cast<unsigned char*>(&data);
        std::copy(pbytes, pbytes+sz, pchar);
        pbytes += sz;
    }
};

#endif
