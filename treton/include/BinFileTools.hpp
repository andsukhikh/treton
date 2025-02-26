#ifndef BIN_FILE_TOOLS 
#define BIN_FILE_TOOLS

#include <fstream>
#include <vector>
#include <type_traits>

class BinaryWriter {
    std::ofstream file;
public:
    BinaryWriter(std::ofstream&& file_);

    template<typename T>
    void writeVector(T& vec);
};


class BinaryReader {
    std::ifstream file;
public:
    BinaryReader(std::ifstream&& file_);

    template<typename T>
    void readVector(T& vec);
};


class Writer : private  BinaryWriter {
public:
    Writer(std::ofstream&& file_);

    template<typename ...Args>
    void write(Args&... vector_list)&&;

};


class Reader : private BinaryReader {
public:
    Reader(std::ifstream&& file_);

    template<typename ...Args>
    void read(Args&... vector_list);

};


Reader::Reader(std::ifstream&& file_) : BinaryReader(std::move(file_)) {}

template<typename ...Args>
void Reader::read(Args&... vector_list) {
    (readVector(vector_list), ...);
}


Writer::Writer(std::ofstream&& file_) : BinaryWriter(std::move(file_)) {}

template<typename ...Args>
void Writer::write(Args&... vector_list)&& {
    (writeVector(vector_list), ...);
}


BinaryWriter::BinaryWriter(std::ofstream&& file_) {
    file = std::move(file_);
}


BinaryReader::BinaryReader(std::ifstream&& file_) {
    file = std::move(file_);
}


template<typename T>
void BinaryWriter::writeVector(T& vec) {
    if (file.is_open()) {
        if constexpr (!std::is_scalar_v<T>) {
            if constexpr (std::is_same_v<T, std::vector<typename T::value_type>>) {
                size_t size = vec.size();
                file.write(reinterpret_cast<char*>(&size), sizeof(size));
                for (auto&& val : vec) {
                    writeVector(val);
                }
            }
        }
        else {
            file.write(reinterpret_cast<char*>(&vec), sizeof(vec));
        }
    }
}


template<typename T>
void BinaryReader::readVector(T& vec) {
    if (file.is_open()) {
        if constexpr (!std::is_scalar_v<T>) {
            if constexpr (std::is_same_v<T, std::vector<typename T::value_type>>) {
                size_t size;
                file.read(reinterpret_cast<char*>(&size), sizeof(size));
                vec.resize(size);
                for (auto&& val : vec) {
                    readVector(val);
                }
            }
        }
        else {
            file.read(reinterpret_cast<char*>(&vec), sizeof(vec));
        }
    }
}

#endif