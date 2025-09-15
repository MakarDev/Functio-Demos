#include <iostream>
#include <vector>
#include <string>
#include <cstdint>
#include <cctype>
#include <algorithm>
#include <iterator>
#include <numeric>
#include <chrono>

// --- Minimal implementations of ImHex types to make the function standalone ---

// Mock Task class to replace the original. It does nothing.
struct Task {
    void update(uint64_t) { /* In the original, this updates a progress bar. */ }
};

// A simple Region struct.
struct Region {
    uint64_t startAddress;
    size_t size;

    uint64_t getStartAddress() const { return startAddress; }
    uint64_t getEndAddress() const { return startAddress + size; }
    size_t getSize() const { return size; }
};

// A mock Provider to supply data from an in-memory vector.
namespace prv {
    class Provider {
    public:
        explicit Provider(const std::vector<uint8_t>& data) : m_data(data) {}

        // Reads data from the vector into a buffer.
        void read(uint64_t address, void* buffer, size_t size) const {
            if (address >= m_data.size()) return;
            size_t readable_size = std::min(size, static_cast<size_t>(m_data.size() - address));
            if (readable_size > 0) {
                std::copy(m_data.begin() + address, m_data.begin() + address + readable_size, static_cast<uint8_t*>(buffer));
            }
        }

        size_t size() const { return m_data.size(); }

    private:
        const std::vector<uint8_t>& m_data;
    };

    // A simple reader to iterate over the provider's data.
    class ProviderReader {
    public:
        ProviderReader(Provider* provider) : m_provider(provider) {}

        class Iterator {
        public:
            using iterator_category = std::forward_iterator_tag;
            using value_type = uint8_t;
            using difference_type = std::ptrdiff_t;
            using pointer = uint8_t*;
            using reference = uint8_t&;

            Iterator(Provider* provider, uint64_t address) : m_provider(provider), m_address(address) {
                if (m_provider && m_address < m_provider->size())
                    m_provider->read(m_address, &m_current_byte, 1);
            }

            uint8_t operator*() const { return m_current_byte; }
            Iterator& operator++() {
                m_address++;
                 if (m_provider && m_address < m_provider->size())
                    m_provider->read(m_address, &m_current_byte, 1);
                return *this;
            }
             bool operator!=(const Iterator& other) const { return m_address != other.m_address; }
             uint64_t getAddress() const { return m_address; }

        private:
            Provider* m_provider = nullptr;
            uint64_t m_address = 0;
            uint8_t m_current_byte = 0;
        };

        void seek(uint64_t address) { m_current_address = address; }
        void setEndAddress(uint64_t address) { m_end_address = address; }

        Iterator begin() { return Iterator(m_provider, m_current_address); }
        Iterator end() { return Iterator(m_provider, m_end_address); }

    private:
        Provider* m_provider;
        uint64_t m_current_address = 0;
        uint64_t m_end_address = 0;
    };
}

// Structs to hold search settings and results.
namespace { // Anonymous namespace to keep these local
    struct Occurrence {
        enum class DecodeType { Binary, ASCII, UTF8, UTF16, Unsigned, Signed, Float, Double };

        Region region;
        DecodeType decodeType;
        std::endian endian;
        bool selected = false;
    };

    namespace SearchSettings {
        enum class StringType { ASCII, UTF8, UTF16LE, UTF16BE, ASCII_UTF16LE, ASCII_UTF16BE };

        struct Strings {
            int minLength = 4;
            bool nullTermination = true;
            StringType type = StringType::ASCII;
            bool lowerCaseLetters = true;
            bool upperCaseLetters = true;
            bool numbers = true;
            bool underscores = true;
            bool symbols = true;
            bool spaces = true;
            bool lineFeeds = true;
        };
    }
}

// --- searchStrings function, adapted from view_find.cpp ---

std::vector<Occurrence> searchStrings(Task &task, prv::Provider *provider, Region searchRegion, const SearchSettings::Strings &settings) {
    using enum SearchSettings::StringType;

    std::vector<Occurrence> results;

    // The original function used recursion for combined types. This mimics that behavior.
    if (settings.type == ASCII_UTF16BE || settings.type == ASCII_UTF16LE) {
        auto newSettings = settings;

        newSettings.type = ASCII;
        auto asciiResults = searchStrings(task, provider, searchRegion, newSettings);
        std::copy(asciiResults.begin(), asciiResults.end(), std::back_inserter(results));

        if (settings.type == ASCII_UTF16BE) {
            newSettings.type = UTF16BE;
            auto utf16Results = searchStrings(task, provider, searchRegion, newSettings);
            std::copy(utf16Results.begin(), utf16Results.end(), std::back_inserter(results));
        } else if (settings.type == ASCII_UTF16LE) {
            newSettings.type = UTF16LE;
            auto utf16Results = searchStrings(task, provider, searchRegion, newSettings);
            std::copy(utf16Results.begin(), utf16Results.end(), std::back_inserter(results));
        }

        return results;
    }

    auto reader = prv::ProviderReader(provider);
    reader.seek(searchRegion.getStartAddress());
    reader.setEndAddress(searchRegion.getEndAddress());

    const auto [decodeType, endian] = [&]() -> std::pair<Occurrence::DecodeType, std::endian> {
        if (settings.type == ASCII)
            return { Occurrence::DecodeType::ASCII, std::endian::native };
        if (settings.type == UTF8)
            return { Occurrence::DecodeType::UTF8, std::endian::native };
        else if (settings.type == SearchSettings::StringType::UTF16BE)
            return { Occurrence::DecodeType::UTF16, std::endian::big };
        else if (settings.type == SearchSettings::StringType::UTF16LE)
            return { Occurrence::DecodeType::UTF16, std::endian::little };
        else
            return { Occurrence::DecodeType::Binary, std::endian::native };
    }();

    int64_t countedCharacters = 0;
    uint64_t startAddress = reader.begin().getAddress();
    uint64_t endAddress = reader.end().getAddress();

    uint64_t codePointWidth = 0;
    int8_t remainingCharacters = 0;
    for (uint8_t byte : reader) {
        bool validChar =
            (settings.lowerCaseLetters    && std::islower(byte))  ||
            (settings.upperCaseLetters    && std::isupper(byte))  ||
            (settings.numbers             && std::isdigit(byte))  ||
            (settings.spaces              && std::isspace(byte) && byte != '\r' && byte != '\n')  ||
            (settings.underscores         && byte == '_')             ||
            (settings.symbols             && std::ispunct(byte) && !std::isspace(byte))  ||
            (settings.lineFeeds           && (byte == '\r' || byte == '\n'));

        if (settings.type == UTF16LE) {
            // Check if second byte of UTF-16 encoded string is 0x00
            if (countedCharacters % 2 == 1)
                validChar = byte == 0x00;
        } else if (settings.type == UTF16BE) {
            // Check if first byte of UTF-16 encoded string is 0x00
            if (countedCharacters % 2 == 0)
                validChar = byte == 0x00;
        } else if (settings.type == UTF8) {
            if ((byte & 0b1000'0000) == 0b0000'0000) {
                // ASCII range
                codePointWidth = 1;
                remainingCharacters = 0;
                validChar = true;
            } else if ((byte & 0b1100'0000) == 0b1000'0000) {
                // Continuation mark
                if (remainingCharacters > 0) {
                    remainingCharacters -= 1;
                    validChar = true;
                } else {
                    countedCharacters -= std::max<int64_t>(0, codePointWidth - (remainingCharacters + 1));
                    codePointWidth = 0;
                    remainingCharacters = 0;
                    validChar = false;
                }
            } else if ((byte & 0b1110'0000) == 0b1100'0000) {
                // Two bytes
                codePointWidth = 2;
                remainingCharacters = codePointWidth - 1;
                validChar = true;
            } else if ((byte & 0b1111'0000) == 0b1110'0000) {
                // Three bytes
                codePointWidth = 3;
                remainingCharacters = codePointWidth - 1;
                validChar = true;
            } else if ((byte & 0b1111'1000) == 0b1111'0000) {
                // Four bytes
                codePointWidth = 4;
                remainingCharacters = codePointWidth - 1;
                validChar = true;
            } else {
                validChar = false;
            }
        }

        if (validChar)
            countedCharacters++;
        if (!validChar || startAddress + countedCharacters == endAddress) {
            if (countedCharacters >= settings.minLength) {
                if (!settings.nullTermination || byte == 0x00) {
                    results.push_back(Occurrence { Region { startAddress, size_t(countedCharacters) }, decodeType, endian, false });
                }
            }

            startAddress += countedCharacters + 1;
            countedCharacters = 0;
        }
    }

    return results;
}


int main() {
    // 1. Create sample data to search through.
    // std::cout << "Generating sample data...\n";
    std::vector<uint8_t> sample_data;
    sample_data.reserve(10 * 1024 * 1024); // 10 MiB

    // Add some known strings
    std::string s1 = "This is a simple test string.";
    std::string s2 = "Another_string_with_numbers_12345.";
    std::string s3 = "Short"; // Should be ignored by default minLength
    std::string s4 = "Ends with null\0, hopefully.";
    std::string s5 = "This string will be repeated many times to increase workload.";

    for (char c : s1) sample_data.push_back(c);
    for (int i = 0; i < 50; ++i) sample_data.push_back(rand() % 256); // Random bytes
    for (char c : s2) sample_data.push_back(c);
    for (int i = 0; i < 50; ++i) sample_data.push_back(rand() % 256);
    for (char c : s3) sample_data.push_back(c);
    for (int i = 0; i < 50; ++i) sample_data.push_back(rand() % 256);
    for (char c : s4) sample_data.push_back(c); // Includes the null terminator
    for (int i = 0; i < 50; ++i) sample_data.push_back(rand() % 256);

    // Add a lot of repeated strings and random data to make the search meaningful
    for (int i = 0; i < 50000; ++i) {
        for (char c : s5) sample_data.push_back(c);
        for (int j = 0; j < 20; ++j) sample_data.push_back(rand() % 256);
    }
    // std::cout << "Data size: " << sample_data.size() / (1024.0 * 1024.0) << " MiB\n";


    // 2. Set up the search parameters.
    Task task;
    prv::Provider provider(sample_data);
    Region searchRegion = { 0, sample_data.size() };
    SearchSettings::Strings settings = {
        .minLength = 5,
        .nullTermination = false, // Find strings that are not null-terminated
        .type = SearchSettings::StringType::ASCII,
        .lowerCaseLetters = true,
        .upperCaseLetters = true,
        .numbers = true,
        .underscores = true,
        .symbols = false,
        .spaces = true,
        .lineFeeds = false
    };

    // 3. Run the search function in a loop for profiling.
    // std::cout << "Starting search loop for profiling...\n";
    const int iterations = 10;
    std::vector<Occurrence> results;


    for (int i = 0; i < iterations; ++i) {
        // Clear results for each run to avoid accumulating forever
        results = searchStrings(task, &provider, searchRegion, settings);
    }




    return 0;
}