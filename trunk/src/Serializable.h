#ifndef SERIALIZABLE_H
#define SERIALIZABLE_H

#include "dtypes.h"

#include <vector>
#include <typeinfo>
#include <typeindex>
#include <memory>
#include <unordered_map>
#include <iterator>

class Serializable {
public:
    virtual bytes_t& serialize() = 0;
    virtual void deserialize(bytes_t &bytes) = 0;
};

typedef uint16_t type_t;

class Serializer {
public:
    Serializer() {
        type_ids[std::type_index(typeid(int16_t))] = 1016;
        type_ids[std::type_index(typeid(int32_t))] = 1032;
        type_ids[std::type_index(typeid(int64_t))] = 1064;
        type_ids[std::type_index(typeid(uint16_t))] = 1116;
        type_ids[std::type_index(typeid(uint32_t))] = 1132;
        type_ids[std::type_index(typeid(uint64_t))] = 1164;
        type_ids[std::type_index(typeid(float))] = 2032; // float32_t
        type_ids[std::type_index(typeid(double))] = 2064; // float64_t
        type_ids[std::type_index(typeid(char))] = 3000;
    }

    void load(const bytes_t &bytes) {
        cpbytes = bytes.begin();
        cpbytes_start = cpbytes;
    }

    void save(bytes_t &bytes) {
        size_t nItems = buffer_types.size();

        // clear and resize bytes
        bytes.clear();
        bytes.resize(size());

        // save items to bytes
        pbytes = bytes.begin();

        for (size_t i = 0; i < nItems; ++i) {
            if (buffer_isVector[i]) {
                serialize(&type_ids[std::type_index(typeid(storage_idx_t))], sizeof(type_t));
                serialize(&buffer_nums[i], sizeof(storage_idx_t));
            }

            serialize(&buffer_types[i], sizeof(type_t));
            serialize(buffer_pitems[i], buffer_nums[i]*buffer_sizes[i]);
        }

        // the buffers are not needed anymore
        buffer_isVector.clear();
        buffer_types.clear();
        buffer_pitems.clear();
        buffer_nums.clear();
        buffer_sizes.clear();
    }

    size_t position() {
        size_t sz = size();
        if (sz > 0) {
            return sz;
        } else {
            return std::distance(cpbytes_start,cpbytes);
        }
    }

    template<class T>
    friend void operator<< (Serializer &s, const T &data) {
        s.buffer_isVector.push_back(false);
        s.buffer_types.push_back(s.type_ids[std::type_index(typeid(T))]);
        s.buffer_pitems.push_back(reinterpret_cast<const byte_t*>(&data));
        s.buffer_nums.push_back(1);
        s.buffer_sizes.push_back(sizeof(T));
    }

    template<class T>
    friend void operator<< (Serializer &s, const std::vector<T> &data) {
        s.buffer_isVector.push_back(true);
        s.buffer_types.push_back(s.type_ids[std::type_index(typeid(T))]);
        s.buffer_pitems.push_back(reinterpret_cast<const byte_t*>(&data[0]));
        s.buffer_nums.push_back(data.size());
        s.buffer_sizes.push_back(sizeof(T));
    }

    template<class T>
    friend void operator>> (Serializer &s, T &data) {
        s.deserialize(data, 1);
    }

    template<class T>
    friend void operator>> (Serializer &s, std::vector<T> &data) {
        size_t szvector;
        s.deserialize(szvector, 1);
        s.deserialize(data, szvector);
    }


    /*template<class T>
    void serializeItem(bytes_t::iterator &pbytes, const T &data) {
        serializeItem(pbytes, data, 1);
    }

    template<class T>
    void serializeItem(bytes_t::iterator &pbytes, const std::vector<T> &data) {
        size_t szvector = data.size();
        serializeItem(pbytes, szvector, 1);
        serializeItem(pbytes, data, szvector);
    }*/



    /*template<class T>
    void serializeItem(bytes_t::iterator &pbytes, const T &data, const size_t num) {
        uint16_t type = type_ids[std::type_index(typeid(T))];

        auto data_reinterpreted = reinterpret_cast<const byte_t*>(&type);
        std::copy(data_reinterpreted, data_reinterpreted+sizeof(type), pbytes);
        pbytes += sizeof(type);

        data_reinterpreted = reinterpret_cast<const byte_t*>(&data);
        std::copy(data_reinterpreted, data_reinterpreted+sizeof(T)*num, pbytes);
        pbytes += sizeof(T)*num;
    }

    template<class T>
    void serializeItem(bytes_t::iterator &pbytes, const std::vector<T> &data, const size_t num) {
        uint16_t type = type_ids[std::type_index(typeid(T))];

        auto data_reinterpreted = reinterpret_cast<const byte_t*>(&type);
        std::copy(data_reinterpreted, data_reinterpreted+sizeof(type), pbytes);
        pbytes += sizeof(type);

        data_reinterpreted = reinterpret_cast<const byte_t*>(&data[0]);
        std::copy(data_reinterpreted, data_reinterpreted+sizeof(T)*num, pbytes);
        pbytes += sizeof(T)*num;
    }*/







    /**/




    /*size_t sz = szvector*szitem;
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
        pbytes += sz;*/


























    /*template<class T>
    void serializeItem(bytes_t::iterator &pbytes, T *pointer, size_t num) {
        size_t sz = sizeof(T)*(num);
        auto pchar = reinterpret_cast<const byte_t*>(pointer);
        std::copy(pchar, pchar+sz, pbytes);
        pbytes += sz;
    }

    template<class T>
    void deserializeItem(bytes_t::iterator &pbytes, T *pointer, size_t num) {
        size_t sz = sizeof(T)*(num);
        auto pchar = reinterpret_cast<byte_t*>(pointer);
        std::copy(pbytes, pbytes+sz, pchar);
        pbytes += sz;
    }*/

private:
    size_t size() {
        size_t nItems = buffer_types.size();
        size_t totalsize = 0;

        for (size_t i = 0; i < nItems; ++i) {
            if (buffer_isVector[i]) {
                totalsize += sizeof(type_t);
                totalsize += sizeof(storage_idx_t);
            }
            totalsize += sizeof(type_t);
            totalsize += buffer_nums[i]*buffer_sizes[i];
        }

        return totalsize;
    }

    template<class T>
    void serialize(const T *pdata, size_t sz) {
        auto data_reinterpreted = reinterpret_cast<const byte_t*>(pdata);
        std::copy(data_reinterpreted, data_reinterpreted+sz, pbytes);
        pbytes += sz;
    }

    template<class T>
    void deserialize(T &data, size_t num) {
        uint16_t type;

        auto pbytes_reinterpreted = reinterpret_cast<const type_t*>(&(*cpbytes));
        std::copy(pbytes_reinterpreted, pbytes_reinterpreted+1, &type);
        cpbytes += sizeof(type);

        if (type == type_ids[std::type_index(typeid(int16_t))]) {
            auto pbytes_reinterpreted = reinterpret_cast<const int16_t*>(&(*cpbytes));
            std::copy(pbytes_reinterpreted, pbytes_reinterpreted+num, &data);
            cpbytes += sizeof(int16_t)*num;
        } else if (type == type_ids[std::type_index(typeid(int32_t))]) {
            auto pbytes_reinterpreted = reinterpret_cast<const int32_t*>(&(*cpbytes));
            std::copy(pbytes_reinterpreted, pbytes_reinterpreted+num, &data);
            cpbytes += sizeof(int32_t)*num;
        } else if (type == type_ids[std::type_index(typeid(int64_t))]) {
            auto pbytes_reinterpreted = reinterpret_cast<const int64_t*>(&(*cpbytes));
            std::copy(pbytes_reinterpreted, pbytes_reinterpreted+num, &data);
            cpbytes += sizeof(int64_t)*num;
        } else if (type == type_ids[std::type_index(typeid(uint16_t))]) {
            auto pbytes_reinterpreted = reinterpret_cast<const uint16_t*>(&(*cpbytes));
            std::copy(pbytes_reinterpreted, pbytes_reinterpreted+num, &data);
            cpbytes += sizeof(uint16_t)*num;
        } else if (type == type_ids[std::type_index(typeid(uint32_t))]) {
            auto pbytes_reinterpreted = reinterpret_cast<const uint32_t*>(&(*cpbytes));
            std::copy(pbytes_reinterpreted, pbytes_reinterpreted+num, &data);
            cpbytes += sizeof(uint32_t)*num;
        } else if (type == type_ids[std::type_index(typeid(uint64_t))]) {
            auto pbytes_reinterpreted = reinterpret_cast<const uint64_t*>(&(*cpbytes));
            std::copy(pbytes_reinterpreted, pbytes_reinterpreted+num, &data);
            cpbytes += sizeof(uint64_t)*num;
        } else if (type == type_ids[std::type_index(typeid(float))]) {
            auto pbytes_reinterpreted = reinterpret_cast<const float*>(&(*cpbytes));
            std::copy(pbytes_reinterpreted, pbytes_reinterpreted+num, &data);
            cpbytes += sizeof(float)*num;
        } else if (type == type_ids[std::type_index(typeid(double))]) {
            auto pbytes_reinterpreted = reinterpret_cast<const double*>(&(*cpbytes));
            std::copy(pbytes_reinterpreted, pbytes_reinterpreted+num, &data);
            cpbytes += sizeof(double)*num;
        } else if (type == type_ids[std::type_index(typeid(char))]) {
            auto pbytes_reinterpreted = reinterpret_cast<const char*>(&(*cpbytes));
            std::copy(pbytes_reinterpreted, pbytes_reinterpreted+num, &data);
            cpbytes += sizeof(char)*num;
        }else {
            auto pbytes_reinterpreted = reinterpret_cast<const T*>(&(*cpbytes));
            std::copy(pbytes_reinterpreted, pbytes_reinterpreted+num, &data);
            cpbytes += sizeof(T)*num;
        }
    }

    template<class T>
    void deserialize(std::vector<T> &data, size_t num) {
        uint16_t type;

        auto pbytes_reinterpreted = reinterpret_cast<const type_t*>(&(*cpbytes));
        std::copy(pbytes_reinterpreted, pbytes_reinterpreted+1, &type);
        cpbytes += sizeof(type);

        if (type == type_ids[std::type_index(typeid(int16_t))]) {
            auto pbytes_reinterpreted = reinterpret_cast<const int16_t*>(&(*cpbytes));
            data = std::vector<T>(pbytes_reinterpreted, pbytes_reinterpreted+num);
            cpbytes += sizeof(int16_t)*num;
        } else if (type == type_ids[std::type_index(typeid(int32_t))]) {
            auto pbytes_reinterpreted = reinterpret_cast<const int32_t*>(&(*cpbytes));
            data = std::vector<T>(pbytes_reinterpreted, pbytes_reinterpreted+num);
            cpbytes += sizeof(int32_t)*num;
        } else if (type == type_ids[std::type_index(typeid(int64_t))]) {
            auto pbytes_reinterpreted = reinterpret_cast<const int64_t*>(&(*cpbytes));
            data = std::vector<T>(pbytes_reinterpreted, pbytes_reinterpreted+num);
            cpbytes += sizeof(int64_t)*num;
        } else if (type == type_ids[std::type_index(typeid(uint16_t))]) {
            auto pbytes_reinterpreted = reinterpret_cast<const uint16_t*>(&(*cpbytes));
            data = std::vector<T>(pbytes_reinterpreted, pbytes_reinterpreted+num);
            cpbytes += sizeof(uint16_t)*num;
        } else if (type == type_ids[std::type_index(typeid(uint32_t))]) {
            auto pbytes_reinterpreted = reinterpret_cast<const uint32_t*>(&(*cpbytes));
            data = std::vector<T>(pbytes_reinterpreted, pbytes_reinterpreted+num);
            cpbytes += sizeof(uint32_t)*num;
        } else if (type == type_ids[std::type_index(typeid(uint64_t))]) {
            auto pbytes_reinterpreted = reinterpret_cast<const uint64_t*>(&(*cpbytes));
            data = std::vector<T>(pbytes_reinterpreted, pbytes_reinterpreted+num);
            cpbytes += sizeof(uint64_t)*num;
        } else if (type == type_ids[std::type_index(typeid(float))]) {
            auto pbytes_reinterpreted = reinterpret_cast<const float*>(&(*cpbytes));
            data = std::vector<T>(pbytes_reinterpreted, pbytes_reinterpreted+num);
            cpbytes += sizeof(float)*num;
        } else if (type == type_ids[std::type_index(typeid(double))]) {
            auto pbytes_reinterpreted = reinterpret_cast<const double*>(&(*cpbytes));
            data = std::vector<T>(pbytes_reinterpreted, pbytes_reinterpreted+num);
            cpbytes += sizeof(double)*num;
        } else if (type == type_ids[std::type_index(typeid(char))]) {
            auto pbytes_reinterpreted = reinterpret_cast<const char*>(&(*cpbytes));
            data = std::vector<T>(pbytes_reinterpreted, pbytes_reinterpreted+num);
            cpbytes += sizeof(char)*num;
        }else {
            auto pbytes_reinterpreted = reinterpret_cast<const T*>(&(*cpbytes));
            data = std::vector<T>(pbytes_reinterpreted, pbytes_reinterpreted+num);
            cpbytes += sizeof(T)*num;
        }
    }

    std::unordered_map<std::type_index, uint16_t> type_ids;
    bytes_t::iterator pbytes;
    bytes_t::const_iterator cpbytes;
    bytes_t::const_iterator cpbytes_start;
    std::vector<bool> buffer_isVector;
    std::vector<uint16_t> buffer_types;
    std::vector<const byte_t*> buffer_pitems;
    std::vector<storage_idx_t> buffer_nums;
    std::vector<size_t> buffer_sizes;
};

#endif
