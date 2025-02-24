#include <fstream>
#include <type_traits>


class Writter {
    std::ofstream& file;
public:
    Writter(std::ofstream& file_) : file(file_) {}

    template<typename T>
    void operator() (T& vec) {
        if (file.is_open()) {
            if constexpr (!std::is_scalar_v<T>) {
                if constexpr (std::is_same_v<T, std::vector<typename T::value_type>>) {
                    size_t size = vec.size();
                    file.write(reinterpret_cast<char*>(&size), sizeof(size));
                    for (auto&& val : vec) {
                        (*this)(val);
                    }
                }
            }
            else {
                file.write(reinterpret_cast<char*>(&vec), sizeof(vec));
            }
        }
    }
};


class Reader {
    std::ifstream& file;
public:
    Reader(std::ifstream& file_) : file(file_) {}

    template<typename T>
    void operator() (T& vec) {
        if (file.is_open()) {
            if constexpr (!std::is_scalar_v<T>) {
                if constexpr (std::is_same_v<T, std::vector<typename T::value_type>>) {
                    size_t size;
                    file.read(reinterpret_cast<char*>(&size), sizeof(size));
                    vec.resize(size);
                    for (auto&& val : vec) {
                        (*this)(val);
                    }
                }
            }
            else {
                file.read(reinterpret_cast<char*>(&vec), sizeof(vec));
            }
        }
    }
}