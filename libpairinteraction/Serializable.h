/*
 * Copyright (c) 2016 Sebastian Weber, Henri Menke. All rights reserved.
 *
 * Licensed under the Apache License, Version 2.0 (the "License");
 * you may not use this file except in compliance with the License.
 * You may obtain a copy of the License at
 *
 *     http://www.apache.org/licenses/LICENSE-2.0
 *
 * Unless required by applicable law or agreed to in writing, software
 * distributed under the License is distributed on an "AS IS" BASIS,
 * WITHOUT WARRANTIES OR CONDITIONS OF ANY KIND, either express or implied.
 * See the License for the specific language governing permissions and
 * limitations under the License.
 */

#ifndef SERIALIZABLE_H
#define SERIALIZABLE_H

#include "dtypes.h"

#include <iterator>
#include <memory>
#include <typeindex>
#include <typeinfo>
#include <unordered_map>
#include <vector>

class Serializable {
public:
    virtual ~Serializable() = default;
    virtual bytes_t &serialize() = 0;
    virtual void deserialize(bytes_t &bytes) = 0;
};

typedef uint16_t type_t;

class Serializer {
public:
    Serializer() {
        type_ids[std::type_index(typeid(int8_t))] = 1008;
        type_ids[std::type_index(typeid(int16_t))] = 1016;
        type_ids[std::type_index(typeid(int32_t))] = 1032;
        type_ids[std::type_index(typeid(int64_t))] = 1064;
        type_ids[std::type_index(typeid(uint8_t))] = 1108;
        type_ids[std::type_index(typeid(uint16_t))] = 1116;
        type_ids[std::type_index(typeid(uint32_t))] = 1132;
        type_ids[std::type_index(typeid(uint64_t))] = 1164;
        type_ids[std::type_index(typeid(float))] = 2032;  // float32_t
        type_ids[std::type_index(typeid(double))] = 2064; // float64_t
        type_ids[std::type_index(typeid(char))] = 3000;
        type_ids[std::type_index(typeid(bool))] = 4000;
    }

    void load(const bytes_t &bytes) {
        cpbytes = bytes.begin();
        cpbytes_start = cpbytes;
        cpbytes_end = bytes.end();
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
            serialize(buffer_pitems[i], buffer_nums[i] * buffer_sizes[i]);
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
        }
        return std::distance(cpbytes_start, cpbytes);
    }

    template <class T>
    friend void operator<<(Serializer &s, const T &data) {
        s.buffer_isVector.push_back(false);
        s.buffer_types.push_back(s.type_ids[std::type_index(typeid(T))]);
        s.buffer_pitems.push_back(reinterpret_cast<const byte_t *>(&data));
        s.buffer_nums.push_back(1);
        s.buffer_sizes.push_back(sizeof(T));
    }

    template <class T>
    friend void operator<<(Serializer &s, const std::vector<T> &data) {
        s.buffer_isVector.push_back(true);
        s.buffer_types.push_back(s.type_ids[std::type_index(typeid(T))]);
        s.buffer_pitems.push_back(reinterpret_cast<const byte_t *>(&data[0]));
        s.buffer_nums.push_back(data.size());
        s.buffer_sizes.push_back(sizeof(T));
    }

    template <class T>
    friend void operator>>(Serializer &s, T &data) {
        s.deserialize(data, 1);
    }

    template <class T>
    friend void operator>>(Serializer &s, std::vector<T> &data) {
        size_t szvector;
        s.deserialize(szvector, 1);
        s.deserialize(data, szvector);
    }

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
            totalsize += buffer_nums[i] * buffer_sizes[i];
        }

        return totalsize;
    }

    template <class T>
    void serialize(const T *pdata, size_t sz) {
        auto data_reinterpreted = reinterpret_cast<const byte_t *>(pdata);
        std::copy(data_reinterpreted, data_reinterpreted + sz, pbytes);
        pbytes += sz;
    }

    template <class T>
    void deserialize(T &data, size_t num) {
        uint16_t type = 0;

        if (cpbytes + sizeof(uint16_t) > cpbytes_end) {
            throw std::runtime_error(
                "Corrupted data discovered."); // TODO use checksum of "bytes" (write checksum into
                                               // the beginning of "bytes")
        }
        auto pbytes_reinterpreted = reinterpret_cast<const type_t *>(&(*cpbytes));
        std::copy(pbytes_reinterpreted, pbytes_reinterpreted + 1, &type);
        cpbytes += sizeof(uint16_t);

        if (type == type_ids[std::type_index(typeid(int8_t))]) {
            if (cpbytes + sizeof(int8_t) * num > cpbytes_end) {
                throw std::runtime_error("Corrupted data discovered.");
            }
            auto pbytes_reinterpreted = reinterpret_cast<const int8_t *>(&(*cpbytes));
            std::copy(pbytes_reinterpreted, pbytes_reinterpreted + num, &data);
            cpbytes += sizeof(int8_t) * num;
        } else if (type == type_ids[std::type_index(typeid(int16_t))]) {
            if (cpbytes + sizeof(int16_t) * num > cpbytes_end) {
                throw std::runtime_error("Corrupted data discovered.");
            }
            auto pbytes_reinterpreted = reinterpret_cast<const int16_t *>(&(*cpbytes));
            std::copy(pbytes_reinterpreted, pbytes_reinterpreted + num, &data);
            cpbytes += sizeof(int16_t) * num;
        } else if (type == type_ids[std::type_index(typeid(int32_t))]) {
            if (cpbytes + sizeof(int32_t) * num > cpbytes_end) {
                throw std::runtime_error("Corrupted data discovered.");
            }
            auto pbytes_reinterpreted = reinterpret_cast<const int32_t *>(&(*cpbytes));
            std::copy(pbytes_reinterpreted, pbytes_reinterpreted + num, &data);
            cpbytes += sizeof(int32_t) * num;
        } else if (type == type_ids[std::type_index(typeid(int64_t))]) {
            if (cpbytes + sizeof(int64_t) * num > cpbytes_end) {
                throw std::runtime_error("Corrupted data discovered.");
            }
            auto pbytes_reinterpreted = reinterpret_cast<const int64_t *>(&(*cpbytes));
            std::copy(pbytes_reinterpreted, pbytes_reinterpreted + num, &data);
            cpbytes += sizeof(int64_t) * num;
        } else if (type == type_ids[std::type_index(typeid(uint8_t))]) {
            if (cpbytes + sizeof(uint8_t) * num > cpbytes_end) {
                throw std::runtime_error("Corrupted data discovered.");
            }
            auto pbytes_reinterpreted = reinterpret_cast<const uint8_t *>(&(*cpbytes));
            std::copy(pbytes_reinterpreted, pbytes_reinterpreted + num, &data);
            cpbytes += sizeof(uint8_t) * num;
        } else if (type == type_ids[std::type_index(typeid(uint16_t))]) {
            if (cpbytes + sizeof(uint16_t) * num > cpbytes_end) {
                throw std::runtime_error("Corrupted data discovered.");
            }
            auto pbytes_reinterpreted = reinterpret_cast<const uint16_t *>(&(*cpbytes));
            std::copy(pbytes_reinterpreted, pbytes_reinterpreted + num, &data);
            cpbytes += sizeof(uint16_t) * num;
        } else if (type == type_ids[std::type_index(typeid(uint32_t))]) {
            if (cpbytes + sizeof(uint32_t) * num > cpbytes_end) {
                throw std::runtime_error("Corrupted data discovered.");
            }
            auto pbytes_reinterpreted = reinterpret_cast<const uint32_t *>(&(*cpbytes));
            std::copy(pbytes_reinterpreted, pbytes_reinterpreted + num, &data);
            cpbytes += sizeof(uint32_t) * num;
        } else if (type == type_ids[std::type_index(typeid(uint64_t))]) {
            if (cpbytes + sizeof(uint64_t) * num > cpbytes_end) {
                throw std::runtime_error("Corrupted data discovered.");
            }
            auto pbytes_reinterpreted = reinterpret_cast<const uint64_t *>(&(*cpbytes));
            std::copy(pbytes_reinterpreted, pbytes_reinterpreted + num, &data);
            cpbytes += sizeof(uint64_t) * num;
        } else if (type == type_ids[std::type_index(typeid(float))]) {
            if (cpbytes + sizeof(float) * num > cpbytes_end) {
                throw std::runtime_error("Corrupted data discovered.");
            }
            auto pbytes_reinterpreted = reinterpret_cast<const float *>(&(*cpbytes));
            std::copy(pbytes_reinterpreted, pbytes_reinterpreted + num, &data);
            cpbytes += sizeof(float) * num;
        } else if (type == type_ids[std::type_index(typeid(double))]) {
            if (cpbytes + sizeof(double) * num > cpbytes_end) {
                throw std::runtime_error("Corrupted data discovered.");
            }
            auto pbytes_reinterpreted = reinterpret_cast<const double *>(&(*cpbytes));
            std::copy(pbytes_reinterpreted, pbytes_reinterpreted + num, &data);
            cpbytes += sizeof(double) * num;
        } else if (type == type_ids[std::type_index(typeid(char))]) {
            if (cpbytes + sizeof(char) * num > cpbytes_end) {
                throw std::runtime_error("Corrupted data discovered.");
            }
            auto pbytes_reinterpreted = reinterpret_cast<const char *>(&(*cpbytes));
            std::copy(pbytes_reinterpreted, pbytes_reinterpreted + num, &data);
            cpbytes += sizeof(char) * num;
        } else if (type == type_ids[std::type_index(typeid(bool))]) {
            if (cpbytes + sizeof(bool) * num > cpbytes_end) {
                throw std::runtime_error("Corrupted data discovered.");
            }
            auto pbytes_reinterpreted = reinterpret_cast<const bool *>(&(*cpbytes));
            std::copy(pbytes_reinterpreted, pbytes_reinterpreted + num, &data);
            cpbytes += sizeof(bool) * num;
        } else {
            throw std::runtime_error("Corrupted data discovered.");
            // auto pbytes_reinterpreted = reinterpret_cast<const T*>(&(*cpbytes));
            // std::copy(pbytes_reinterpreted, pbytes_reinterpreted+num, &data);
            // cpbytes += sizeof(T)*num;
        }
    }

    template <class T>
    void deserialize(std::vector<T> &data, size_t num) {
        uint16_t type = 0;

        if (cpbytes + sizeof(uint16_t) > cpbytes_end) {
            throw std::runtime_error("Corrupted data discovered.");
        }
        auto pbytes_reinterpreted = reinterpret_cast<const type_t *>(&(*cpbytes));
        std::copy(pbytes_reinterpreted, pbytes_reinterpreted + 1, &type);
        cpbytes += sizeof(uint16_t);

        if (type == type_ids[std::type_index(typeid(int8_t))]) {
            if (cpbytes + sizeof(int8_t) * num > cpbytes_end) {
                throw std::runtime_error("Corrupted data discovered.");
            }
            auto pbytes_reinterpreted = reinterpret_cast<const int8_t *>(&(*cpbytes));
            std::vector<int8_t> vectmp(pbytes_reinterpreted, pbytes_reinterpreted + num);
            data.insert(data.end(), vectmp.begin(), vectmp.end());
            cpbytes += sizeof(int8_t) * num;
        } else if (type == type_ids[std::type_index(typeid(int16_t))]) {
            if (cpbytes + sizeof(int16_t) * num > cpbytes_end) {
                throw std::runtime_error("Corrupted data discovered.");
            }
            auto pbytes_reinterpreted = reinterpret_cast<const int16_t *>(&(*cpbytes));
            std::vector<int16_t> vectmp(pbytes_reinterpreted, pbytes_reinterpreted + num);
            data.insert(data.end(), vectmp.begin(), vectmp.end());
            cpbytes += sizeof(int16_t) * num;
        } else if (type == type_ids[std::type_index(typeid(int32_t))]) {
            if (cpbytes + sizeof(int32_t) * num > cpbytes_end) {
                throw std::runtime_error("Corrupted data discovered.");
            }
            auto pbytes_reinterpreted = reinterpret_cast<const int32_t *>(&(*cpbytes));
            std::vector<int32_t> vectmp(pbytes_reinterpreted, pbytes_reinterpreted + num);
            data.insert(data.end(), vectmp.begin(), vectmp.end());
            cpbytes += sizeof(int32_t) * num;
        } else if (type == type_ids[std::type_index(typeid(int64_t))]) {
            if (cpbytes + sizeof(int64_t) * num > cpbytes_end) {
                throw std::runtime_error("Corrupted data discovered.");
            }
            auto pbytes_reinterpreted = reinterpret_cast<const int64_t *>(&(*cpbytes));
            std::vector<int64_t> vectmp(pbytes_reinterpreted, pbytes_reinterpreted + num);
            data.insert(data.end(), vectmp.begin(), vectmp.end());
            cpbytes += sizeof(int64_t) * num;
        } else if (type == type_ids[std::type_index(typeid(uint8_t))]) {
            if (cpbytes + sizeof(uint8_t) * num > cpbytes_end) {
                throw std::runtime_error("Corrupted data discovered.");
            }
            auto pbytes_reinterpreted = reinterpret_cast<const uint8_t *>(&(*cpbytes));
            std::vector<uint8_t> vectmp(pbytes_reinterpreted, pbytes_reinterpreted + num);
            data.insert(data.end(), vectmp.begin(), vectmp.end());
            cpbytes += sizeof(uint8_t) * num;
        } else if (type == type_ids[std::type_index(typeid(uint16_t))]) {
            if (cpbytes + sizeof(uint16_t) * num > cpbytes_end) {
                throw std::runtime_error("Corrupted data discovered.");
            }
            auto pbytes_reinterpreted = reinterpret_cast<const uint16_t *>(&(*cpbytes));
            std::vector<uint16_t> vectmp(pbytes_reinterpreted, pbytes_reinterpreted + num);
            data.insert(data.end(), vectmp.begin(), vectmp.end());
            cpbytes += sizeof(uint16_t) * num;
        } else if (type == type_ids[std::type_index(typeid(uint32_t))]) {
            if (cpbytes + sizeof(uint32_t) * num > cpbytes_end) {
                throw std::runtime_error("Corrupted data discovered.");
            }
            auto pbytes_reinterpreted = reinterpret_cast<const uint32_t *>(&(*cpbytes));
            std::vector<uint32_t> vectmp(pbytes_reinterpreted, pbytes_reinterpreted + num);
            data.insert(data.end(), vectmp.begin(), vectmp.end());
            cpbytes += sizeof(uint32_t) * num;
        } else if (type == type_ids[std::type_index(typeid(uint64_t))]) {
            if (cpbytes + sizeof(uint64_t) * num > cpbytes_end) {
                throw std::runtime_error("Corrupted data discovered.");
            }
            auto pbytes_reinterpreted = reinterpret_cast<const uint64_t *>(&(*cpbytes));
            std::vector<uint64_t> vectmp(pbytes_reinterpreted, pbytes_reinterpreted + num);
            data.insert(data.end(), vectmp.begin(), vectmp.end());
            cpbytes += sizeof(uint64_t) * num;
        } else if (type == type_ids[std::type_index(typeid(float))]) {
            if (cpbytes + sizeof(float) * num > cpbytes_end) {
                throw std::runtime_error("Corrupted data discovered.");
            }
            auto pbytes_reinterpreted = reinterpret_cast<const float *>(&(*cpbytes));
            std::vector<float> vectmp(pbytes_reinterpreted, pbytes_reinterpreted + num);
            data.insert(data.end(), vectmp.begin(), vectmp.end());
            cpbytes += sizeof(float) * num;
        } else if (type == type_ids[std::type_index(typeid(double))]) {
            if (cpbytes + sizeof(double) * num > cpbytes_end) {
                throw std::runtime_error("Corrupted data discovered.");
            }
            auto pbytes_reinterpreted = reinterpret_cast<const double *>(&(*cpbytes));
            std::vector<double> vectmp(pbytes_reinterpreted, pbytes_reinterpreted + num);
            data.insert(data.end(), vectmp.begin(), vectmp.end());
            cpbytes += sizeof(double) * num;
        } else if (type == type_ids[std::type_index(typeid(char))]) {
            if (cpbytes + sizeof(char) * num > cpbytes_end) {
                throw std::runtime_error("Corrupted data discovered.");
            }
            auto pbytes_reinterpreted = reinterpret_cast<const char *>(&(*cpbytes));
            std::vector<char> vectmp(pbytes_reinterpreted, pbytes_reinterpreted + num);
            data.insert(data.end(), vectmp.begin(), vectmp.end());
            cpbytes += sizeof(char) * num;
        } else if (type == type_ids[std::type_index(typeid(bool))]) {
            if (cpbytes + sizeof(bool) * num > cpbytes_end) {
                throw std::runtime_error("Corrupted data discovered.");
            }
            auto pbytes_reinterpreted = reinterpret_cast<const bool *>(&(*cpbytes));
            std::vector<bool> vectmp(pbytes_reinterpreted, pbytes_reinterpreted + num);
            data.insert(data.end(), vectmp.begin(), vectmp.end());
            cpbytes += sizeof(bool) * num;
        } else {
            throw std::runtime_error("Corrupted data discovered.");
            // auto pbytes_reinterpreted = reinterpret_cast<const T*>(&(*cpbytes));
            // std::copy(pbytes_reinterpreted, pbytes_reinterpreted+num,back_inserter(data));
            // cpbytes += sizeof(T)*num;
        }
    }

    std::unordered_map<std::type_index, uint16_t> type_ids;
    bytes_t::iterator pbytes;
    bytes_t::const_iterator cpbytes;
    bytes_t::const_iterator cpbytes_start;
    bytes_t::const_iterator cpbytes_end;
    std::vector<bool> buffer_isVector;
    std::vector<uint16_t> buffer_types;
    std::vector<const byte_t *> buffer_pitems;
    std::vector<storage_idx_t> buffer_nums;
    std::vector<size_t> buffer_sizes;
};

#endif
