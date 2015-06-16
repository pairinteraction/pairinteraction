#ifndef SERIALIZABLE_H
#define SERIALIZABLE_H

#include <vector>
#include <typeinfo>
#include <memory>

class Serializable {
public:
    virtual bytes_t serialize() = 0;
    virtual void deserialize(bytes_t &bytes) = 0;

protected:
    template<class T>
    void serializeItem(bytes_t::iterator &pbytes, std::vector<T> &data) {
        size_t szvector = data.size();

        size_t sz = sizeof(szvector);
        auto pchar = reinterpret_cast<const byte_t*>(&szvector);
        std::copy(pchar, pchar+sz, pbytes);
        pbytes += sz;

        size_t szitem = sizeof(data[0]);

        sz = sizeof(szitem);
        pchar = reinterpret_cast<const byte_t*>(&szitem);
        std::copy(pchar, pchar+sz, pbytes);
        pbytes += sz;

        size_t typeitem = typeid(data[0]).hash_code();

        sz = sizeof(typeitem);
        pchar = reinterpret_cast<const byte_t*>(&typeitem);
        std::copy(pchar, pchar+sz, pbytes);
        pbytes += sz;

        sz = szvector*szitem;
        pchar = reinterpret_cast<const byte_t*>(&data[0]);
        std::copy(pchar, pchar+sz, pbytes);
        pbytes += sz;
    }

    template<class T>
    void serializeItem(bytes_t::iterator &pbytes, T &data) {
        size_t sz = sizeof(data);
        auto pchar = reinterpret_cast<const byte_t*>(&data);
        std::copy(pchar, pchar+sz, pbytes);
        pbytes += sz;
    }

    template<class T>
    void deserializeItem(bytes_t::iterator &pbytes, std::vector<T> &data) {
        size_t szvector;

        size_t sz = sizeof(szvector);
        auto pchar = reinterpret_cast<byte_t*>(&szvector);
        std::copy(pbytes, pbytes+sz, pchar);
        pbytes += sz;

        size_t szitem;

        sz = sizeof(szitem);
        pchar = reinterpret_cast<byte_t*>(&szitem);
        std::copy(pbytes, pbytes+sz, pchar);
        pbytes += sz;

        size_t typeitem;

        sz = sizeof(typeitem);
        pchar = reinterpret_cast<byte_t*>(&typeitem);
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
            auto pdata = reinterpret_cast<byte_t*>(&data[0]);
            std::copy(pitem, pitem+szvector, pdata);
        }
        pbytes += sz;
    }

    template<class T>
    void deserializeItem(bytes_t::iterator &pbytes, T &data) {
        size_t sz = sizeof(data);
        auto pchar = reinterpret_cast<byte_t*>(&data);
        std::copy(pbytes, pbytes+sz, pchar);
        pbytes += sz;
    }
};

#endif
