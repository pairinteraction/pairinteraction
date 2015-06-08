#ifndef SERIALIZABLE_H
#define SERIALIZABLE_H

#include <vector>

class Serializable {
public:
    virtual std::vector<unsigned char> serialize() = 0;
    virtual void deserialize(std::vector<unsigned char> &bytes) = 0;
};

#endif
