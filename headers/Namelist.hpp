#ifndef NAMELIST_H
#define NAMELIST_H

#include <iostream>
#include <fstream>
#include <string>
#include <vector>
#include <variant>
#include <sstream>
#include <unordered_map>
#include <stdexcept>
#include <tuple>
#include <functional>
#include <regex>

// Типы данных, поддерживаемые в namelist
using VariableValue = std::variant<
    std::reference_wrapper<int>,
    std::reference_wrapper<double>,
    std::reference_wrapper<std::string>,
    std::reference_wrapper<std::vector<int>>,
    std::reference_wrapper<std::vector<double>>
>;

class namelist {
public:
    // Конструктор
    template <typename... Args>
    namelist(const std::string& name, const std::vector<std::string>& variable_names, Args&... vars);

    // Получить имя блока
    const std::string& get_name() const;

    // Печать имен переменных
    void print_variable_names() const;

private:
    std::string name;
    std::vector<VariableValue> variables;
    std::unordered_map<std::string, size_t> var_name_to_index;

    // Добавить переменную
    template <typename T>
    void add_variable(const std::string& var_name, T& var, size_t index);

    // Удаление пробелов из строки
    static std::string strip_whitespace(const std::string& str);

    // Дружественная функция для чтения из файла
    friend void readNamelistFile(namelist& nl, const std::string& filename);
};

// Перегрузка оператора вывода
std::ostream& operator<<(std::ostream& os, const VariableValue& obj);

// Функция для разделения строки на элементы
template <typename T>
std::vector<T> parse_array(const std::string& value);

// Функция чтения из файла
void readNamelistFile(namelist& nl, const std::string& filename);

// Реализация шаблонных функций
template <typename... Args>
namelist::namelist(const std::string& name, const std::vector<std::string>& variable_names, Args&... vars)
    : name(name) {
    size_t index = 0;
    auto it = variable_names.begin();
    ((add_variable(*it++, vars, index++)), ...);
}

template <typename T>
void namelist::add_variable(const std::string& var_name, T& var, size_t index) {
    variables.emplace_back(std::ref(var));
    var_name_to_index[strip_whitespace(var_name)] = index;
}

template <typename T>
std::vector<T> parse_array(const std::string& value) {
    std::vector<T> result;
    std::regex delimiter("[,\\s]+"); // Разделители: запятая или пробелы
    std::sregex_token_iterator begin(value.begin(), value.end(), delimiter, -1), end;

    for (auto it = begin; it != end; ++it) {
        if (!it->str().empty()) {
            if constexpr (std::is_same_v<T, int>) {
                result.push_back(std::stoi(it->str()));
            } else if constexpr (std::is_same_v<T, double>) {
                result.push_back(std::stod(it->str()));
            }
        }
    }
    return result;
}

#endif // NAMELIST_H
