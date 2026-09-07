#ifndef HAMPR_UTILS_JSON_HPP
#define HAMPR_UTILS_JSON_HPP

#include <string>
#include <vector>
#include <map>
#include <memory>
#include <sstream>
#include <stdexcept>
#include <cstdlib>

namespace hampr::json {

enum class JsonType { Null, Bool, Number, String, Array, Object };

class JsonValue {
public:
    JsonType type = JsonType::Null;
    bool bool_value = false;
    double num_value = 0.0;
    std::string str_value;
    std::vector<JsonValue> arr;
    std::map<std::string, JsonValue> obj;

    bool is_object() const { return type == JsonType::Object; }
    bool is_array() const { return type == JsonType::Array; }
    bool is_number() const { return type == JsonType::Number; }
    bool is_string() const { return type == JsonType::String; }

    const JsonValue& at(const std::string& key) const {
        auto it = obj.find(key);
        if (it == obj.end())
            throw std::runtime_error("JSON key not found: " + key);
        return it->second;
    }

    const JsonValue& at(size_t idx) const {
        if (idx >= arr.size())
            throw std::runtime_error("JSON array index out of range");
        return arr[idx];
    }

    double as_number() const { return num_value; }
    int as_int() const { return static_cast<int>(num_value); }
    const std::string& as_string() const { return str_value; }
    bool as_bool() const { return bool_value; }
};

class Parser {
public:
    Parser(const std::string& text) : text_(text), pos_(0) {}

    JsonValue parse() {
        skip_ws();
        return parse_value();
    }

private:
    const std::string& text_;
    size_t pos_;

    char peek() const { return pos_ < text_.size() ? text_[pos_] : '\0'; }
    char get() { return pos_ < text_.size() ? text_[pos_++] : '\0'; }

    void skip_ws() {
        while (pos_ < text_.size() && std::isspace(static_cast<unsigned char>(text_[pos_])))
            pos_++;
    }

    JsonValue parse_value() {
        skip_ws();
        if (pos_ >= text_.size())
            throw std::runtime_error("Unexpected end of JSON");
        char c = peek();
        if (c == '{') return parse_object();
        if (c == '[') return parse_array();
        if (c == '"') return parse_string();
        if (c == 't' || c == 'f') return parse_bool();
        if (c == 'n') return parse_null();
        return parse_number();
    }

    JsonValue parse_object() {
        JsonValue v;
        v.type = JsonType::Object;
        get(); // consume {
        skip_ws();
        if (peek() == '}') { get(); return v; }
        while (true) {
            skip_ws();
            std::string key = parse_string().str_value;
            skip_ws();
            if (get() != ':')
                throw std::runtime_error("Expected ':' in JSON object");
            v.obj[key] = parse_value();
            skip_ws();
            char c = get();
            if (c == ',') continue;
            if (c == '}') break;
            throw std::runtime_error("Expected ',' or '}' in JSON object");
        }
        return v;
    }

    JsonValue parse_array() {
        JsonValue v;
        v.type = JsonType::Array;
        get(); // consume [
        skip_ws();
        if (peek() == ']') { get(); return v; }
        while (true) {
            v.arr.push_back(parse_value());
            skip_ws();
            char c = get();
            if (c == ',') continue;
            if (c == ']') break;
            throw std::runtime_error("Expected ',' or ']' in JSON array");
        }
        return v;
    }

    JsonValue parse_string() {
        JsonValue v;
        v.type = JsonType::String;
        get(); // consume opening "
        std::string result;
        while (pos_ < text_.size() && text_[pos_] != '"') {
            char c = text_[pos_++];
            if (c == '\\' && pos_ < text_.size()) {
                char esc = text_[pos_++];
                switch (esc) {
                    case 'n': result += '\n'; break;
                    case 't': result += '\t'; break;
                    case 'r': result += '\r'; break;
                    case '"': result += '"'; break;
                    case '\\': result += '\\'; break;
                    default: result += esc; break;
                }
            } else {
                result += c;
            }
        }
        get(); // consume closing "
        v.str_value = result;
        return v;
    }

    JsonValue parse_number() {
        JsonValue v;
        v.type = JsonType::Number;
        size_t start = pos_;
        if (peek() == '-') pos_++;
        while (pos_ < text_.size() &&
               (std::isdigit(static_cast<unsigned char>(text_[pos_])) ||
                text_[pos_] == '.' || text_[pos_] == 'e' || text_[pos_] == 'E' ||
                text_[pos_] == '+' || text_[pos_] == '-'))
            pos_++;
        std::string num_str = text_.substr(start, pos_ - start);
        v.num_value = std::strtod(num_str.c_str(), nullptr);
        return v;
    }

    JsonValue parse_bool() {
        JsonValue v;
        v.type = JsonType::Bool;
        if (text_.substr(pos_, 4) == "true") {
            v.bool_value = true;
            pos_ += 4;
        } else if (text_.substr(pos_, 5) == "false") {
            v.bool_value = false;
            pos_ += 5;
        } else {
            throw std::runtime_error("Invalid boolean in JSON");
        }
        return v;
    }

    JsonValue parse_null() {
        JsonValue v;
        v.type = JsonType::Null;
        if (text_.substr(pos_, 4) == "null")
            pos_ += 4;
        else
            throw std::runtime_error("Invalid null in JSON");
        return v;
    }
};

inline JsonValue parse(const std::string& text) {
    Parser p(text);
    return p.parse();
}

} // namespace hampr::json

#endif // HAMPR_UTILS_JSON_HPP
