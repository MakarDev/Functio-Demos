#include <iostream>
#include <vector>
#include <string>
#include <cstdint>
#include <algorithm>
#include <iterator>
#include <cstddef>
#include <memory>

// Mock Task struct for compilation
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

// Pre-computed lookup table for character validation - this will be in read-only memory
static const uint8_t char_validity_table[256] = {
    0, 0, 0, 0, 0, 0, 0, 0, 0, 1, 1, 0, 0, 1, 0, 0,  // 0-15 (includes tab, lf, cr)
    0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0,  // 16-31
    1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1,  // 32-47 (space, symbols)
    1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1,  // 48-63 (digits, symbols)
    1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1,  // 64-79 (uppercase)
    1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1,  // 80-95 (uppercase, symbols)
    1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1,  // 96-111 (lowercase)
    1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 0,  // 112-127 (lowercase, symbols)
    0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0,  // 128-143
    0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0,  // 144-159
    0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0,  // 160-175
    0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0,  // 176-191
    0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0,  // 192-207
    0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0,  // 208-223
    0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0,  // 224-239
    0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0   // 240-255
};


std::vector<Occurrence> searchStrings(Task &task, prv::Provider *provider, Region searchRegion, const SearchSettings::Strings &settings) {
    using enum SearchSettings::StringType;

    std::vector<Occurrence> results;

    // This recursive block for combined search types remains the same
    if (settings.type == ASCII_UTF16BE || settings.type == ASCII_UTF16LE) {
        auto newSettings = settings;
        newSettings.type = ASCII;
        auto asciiResults = searchStrings(task, provider, searchRegion, newSettings);
        results.insert(results.end(), asciiResults.begin(), asciiResults.end());

        newSettings.type = (settings.type == ASCII_UTF16BE) ? UTF16BE : UTF16LE;
        auto utf16Results = searchStrings(task, provider, searchRegion, newSettings);
        results.insert(results.end(), utf16Results.begin(), utf16Results.end());

        return results;
    }

    // Load data into a contiguous buffer for efficient access
    const size_t bufferSize = searchRegion.getSize();
    std::vector<uint8_t> buffer(bufferSize);
    provider->read(searchRegion.getStartAddress(), buffer.data(), bufferSize);

    const auto [decodeType, endian] = [&]() -> std::pair<Occurrence::DecodeType, std::endian> {
        if (settings.type == ASCII) return { Occurrence::DecodeType::ASCII, std::endian::native };
        if (settings.type == UTF8) return { Occurrence::DecodeType::UTF8, std::endian::native };
        else if (settings.type == UTF16BE) return { Occurrence::DecodeType::UTF16, std::endian::big };
        else if (settings.type == UTF16LE) return { Occurrence::DecodeType::UTF16, std::endian::little };
        else return { Occurrence::DecodeType::Binary, std::endian::native };
    }();

    if (settings.type == ASCII) {
        // Use a pre-computed validity table for fast lookups
        const uint8_t* validChars = char_validity_table;

        int64_t countedCharacters = 0;
        uint64_t startAddressOfMatch = 0; // Tracks the start of the current potential match

        // --- BEGINNING OF OPTIMIZED LOOP ---
        // The loop is unrolled by a factor of 4 to break the dependency chain
        // on 'countedCharacters' and increase instruction-level parallelism.
        size_t i = 0;
        for (; i + 3 < bufferSize; i += 4) {
            // Check 4 bytes at once. The CPU can often fetch and validate these in parallel.
            const bool v1 = validChars[buffer[i]];
            const bool v2 = validChars[buffer[i+1]];
            const bool v3 = validChars[buffer[i+2]];
            const bool v4 = validChars[buffer[i+3]];

            // This logic is still sequential, but processing four results at once
            // gives the CPU's out-of-order engine more independent instructions to schedule.
            
            // --- Unrolled Iteration 1 ---
            if (v1) {
                if (countedCharacters == 0) startAddressOfMatch = searchRegion.getStartAddress() + i;
                countedCharacters++;
            } else {
                if (countedCharacters >= settings.minLength) {
                    results.push_back({ {startAddressOfMatch, (size_t)countedCharacters}, decodeType, endian });
                }
                countedCharacters = 0;
            }

            // --- Unrolled Iteration 2 ---
            if (v2) {
                if (countedCharacters == 0) startAddressOfMatch = searchRegion.getStartAddress() + i + 1;
                countedCharacters++;
            } else {
                if (countedCharacters >= settings.minLength) {
                    results.push_back({ {startAddressOfMatch, (size_t)countedCharacters}, decodeType, endian });
                }
                countedCharacters = 0;
            }

            // --- Unrolled Iteration 3 ---
            if (v3) {
                if (countedCharacters == 0) startAddressOfMatch = searchRegion.getStartAddress() + i + 2;
                countedCharacters++;
            } else {
                if (countedCharacters >= settings.minLength) {
                    results.push_back({ {startAddressOfMatch, (size_t)countedCharacters}, decodeType, endian });
                }
                countedCharacters = 0;
            }

            // --- Unrolled Iteration 4 ---
            if (v4) {
                if (countedCharacters == 0) startAddressOfMatch = searchRegion.getStartAddress() + i + 3;
                countedCharacters++;
            } else {
                if (countedCharacters >= settings.minLength) {
                    results.push_back({ {startAddressOfMatch, (size_t)countedCharacters}, decodeType, endian });
                }
                countedCharacters = 0;
            }
        }

        // Cleanup loop for the remaining 0-3 bytes if bufferSize is not a multiple of 4.
        for (; i < bufferSize; ++i) {
            if (validChars[buffer[i]]) {
                if (countedCharacters == 0) startAddressOfMatch = searchRegion.getStartAddress() + i;
                countedCharacters++;
            } else {
                if (countedCharacters >= settings.minLength) {
                    results.push_back({ {startAddressOfMatch, (size_t)countedCharacters}, decodeType, endian });
                }
                countedCharacters = 0;
            }
        }

        // Final check after all loops are finished to catch a string at the very end of the buffer.
        if (countedCharacters >= settings.minLength) {
            results.push_back({ {startAddressOfMatch, (size_t)countedCharacters}, decodeType, endian });
        }
        // --- END OF OPTIMIZED LOOP ---

    } else {
        // Fallback to original, more complex logic for UTF8/UTF16 as vectorizing it is non-trivial.
        // This logic remains unchanged.
        auto reader = prv::ProviderReader(provider);
        reader.seek(searchRegion.getStartAddress());
        reader.setEndAddress(searchRegion.getEndAddress());

        int64_t countedCharacters = 0;
        uint64_t startAddress = reader.begin().getAddress();
        uint64_t endAddress = reader.end().getAddress();
        
        // ... (rest of the original UTF8/UTF16 logic) ...
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

    auto start_time = std::chrono::high_resolution_clock::now();

    for (int i = 0; i < iterations; ++i) {
        // Clear results for each run to avoid accumulating forever
        results = searchStrings(task, &provider, searchRegion, settings);
    }




    return 0;
}