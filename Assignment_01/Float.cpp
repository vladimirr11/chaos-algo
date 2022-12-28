#include <array>
#include <fstream>
#include <iomanip>
#include <iostream>
#include <string>
#include <vector>

// constants used to generate data
const uint32_t MAX_ALLOWED_TIME = 1672531199;
const uint32_t SECONDS_IN_DAY = 86400;
const float MAX_TEMP = 60.0f;
const float MIN_TEMP = -60.0f;
static uint32_t x_state = 123456789; // seed

namespace fpt {

#define sign(x) (x >> 31)
#define exponent(x) ((x << 1) >> 24)
#define mantissa(x) ((x << 9) >> 9)

#define min(a, b) (((a) < (b)) ? (a) : (b))
#define max(a, b) (((a) > (b)) ? (a) : (b))
#define abs_value(x) (x & 0x7FFFFFFF)  // 01111111111111111111111111111111

    struct FType {
        uint8_t sign : 1 {}, exp : 6 {}, m1 : 1 {}, m2 : 8 {}, m3 : 8 {};
    };

    class Float {
    public:
        Float() = default;

        Float(const float f) { setBits(f); }

        Float operator+(const Float& other) const {
            const uint32_t ia = restoreToUint32();
            const uint32_t ib = other.restoreToUint32();

            const uint32_t ires = add(ia, ib);
            const float fres = reinterpret_cast<const float&>(ires);
            Float sum(fres);

            return sum;
        }

        Float operator/(const Float& other) {
            const uint32_t ia = restoreToUint32();
            const uint32_t ib = other.restoreToUint32();

            const uint32_t ires = divide(ia, ib);
            const float fres = reinterpret_cast<const float&>(ires);
            Float ratio(fres);

            return ratio;
        }

        friend std::istream& operator>>(std::istream& istr, Float& f);
        friend std::ostream& operator<<(std::ostream& ostr, const Float& f);

    private:
        static uint32_t add(const uint32_t a, const uint32_t b) {
            uint32_t result_sign;
            uint32_t result_exponent;
            uint32_t result_mantissa;

            const uint32_t diff_sign = sign(a) ^ sign(b);  // determines addition or subtraction

            if (abs_value(a) > abs_value(b)) {
                result_sign = sign(a);
            } else if (abs_value(a) < abs_value(b)) {
                result_sign = sign(b);
            } else {
                result_sign = sign(max(a, b));
            }

            const uint32_t maxarg = max(abs_value(a), abs_value(b));
            const uint32_t minarg = min(abs_value(a), abs_value(b));

            // early return if a == -b
            if (maxarg == minarg && diff_sign) {
                return 0;
            }

            const uint32_t exp_maxarg = exponent(maxarg);
            const uint32_t exp_minarg = exponent(minarg);

            const uint32_t exp_diff = exp_maxarg - exp_minarg;

            result_exponent = exponent(maxarg);

            if (diff_sign) {  // subtraction
                // pull out the hidden bit and then subtract
                result_mantissa = ((mantissa(maxarg) | 0x800000) >> 1) -
                                  (((mantissa(minarg) | 0x800000) >> 1) >> exp_diff);

                // normalize mantissa and adjust exponent
                if (result_mantissa & 0x400000) {  // check for the hidden bit at 23th position
                    result_mantissa =
                        (result_mantissa << 1) & ~(1 << 23);  // put in place and clear
                } else {
                    result_mantissa <<= 1;
                    while (!(result_mantissa & 0x800000)) {  // check the hidden bit at 24th pos
                        result_mantissa <<= 1;
                        result_exponent -= 1;
                    }
                    result_mantissa = result_mantissa & ~(1 << 23);  // clear the hidden bit
                }

            } else {  // addition
                if (exp_diff) {
                    result_mantissa =
                        mantissa(maxarg) + ((mantissa(minarg) | 0x800000) >> exp_diff);
                } else {
                    result_mantissa = mantissa(maxarg) + mantissa(minarg);
                    result_mantissa >>= 1;
                    result_exponent += 1;
                }

                // check for carried out bit
                if ((result_mantissa & 0x800000) >> 23) {
                    result_mantissa >>= 1;
                    if (exp_diff) {
                        result_mantissa = result_mantissa & ~(1 << 22);  // clear the hidden bit
                    }
                    result_exponent += 1;
                }
            }

            const uint32_t result = result_sign << 31 | result_exponent << 23 | result_mantissa;

            return result;
        }

        static uint32_t divide(const uint32_t a, const uint32_t b) {
            const uint32_t result_sign = (a ^ b) & 0x80000000;

            const uint32_t c = mantissa(a) >= mantissa(b);  // c E {0, 1}

            const uint32_t result_exponent = (exponent(a) + 0x7D) - (exponent(b) - c);  // bias 125

            const uint32_t inf = result_sign | 0x7F800000;
            if (result_exponent >= 0xFE)  // overflow detection
                return inf;

            const uint32_t mant_a = (a << 8) | 0x80000000;
            const uint32_t mant_b = (b << 8) | 0x80000000;

            const uint32_t nom = mant_a >> (7 + c);
            const uint32_t denom = mant_b >> 8;

            // Restoring division
            uint32_t quot = 1;
            uint32_t reminder = nom - denom;

            int32_t temp{};
            for (int32_t i = 1; i < 25; i++) {
                temp = (reminder << 1) - denom;

                if (temp < 0) {
                    quot = quot << 1;
                    reminder = temp + denom;
                } else {
                    quot = (quot << 1) + 1;
                    reminder = temp;
                }
            }

            const uint32_t result =
                (result_sign | (result_exponent << 23)) + ((quot >> 1) + (quot & 1));

            return result;
        }

        void setBits(float f) {
            const uint32_t ui = reinterpret_cast<const uint32_t&>(f);

            const uint32_t sign_mask = 0x80000000;  // 10000000000000000000000000000000
            const uint32_t exp_mask = 0x7E000000;   // 01111110000000000000000000000000
            const uint32_t m1_mask = 0x1000000;     // 00000001000000000000000000000000
            const uint32_t m2_mask = 0xFF0000;      // 00000000111111110000000000000000
            const uint32_t m3_mask = ~(sign_mask | exp_mask | m1_mask | m2_mask);

            ftype.sign = (ui & sign_mask) >> 31;
            ftype.exp = (ui & exp_mask) >> 25;
            ftype.m1 = (ui & m1_mask) >> 24;
            ftype.m2 = (ui & m2_mask) >> 16;
            ftype.m3 = (ui & m3_mask) >> 8;
        }

        uint32_t restoreToUint32() const {
            const uint32_t res_ui = ftype.sign << 31 | ftype.exp << 25 | ftype.m1 << 24 |
                                    ftype.m2 << 16 | ftype.m3 << 8;

            return res_ui;
        }

    private:
        FType ftype;
    };

    std::istream& operator>>(std::istream& istr, Float& f) {
        float temp_fl{};
        if (istr.bad()) {
            std::cerr << "Extraction operator failed" << std::endl;
            istr.setstate(std::ios::failbit);
            return istr;
        }
        istr >> temp_fl;
        f.setBits(temp_fl);
        return istr;
    }

    std::ostream& operator<<(std::ostream& ostr, const Float& f) {
        const uint32_t temp_ui = f.restoreToUint32();
        ostr << reinterpret_cast<const float&>(temp_ui);
        return ostr;
    }

}  // namespace fpt

// utility functions for generating and parsing data
namespace {

#define xorshift_utime                                                                  \
    (x_state ^= (x_state << 17), x_state ^= (x_state >> 13), x_state ^= (x_state << 5), \
     x_state % MAX_ALLOWED_TIME)

#define rand_temp (((float)rand() / RAND_MAX) * (MAX_TEMP - MIN_TEMP) + MIN_TEMP)
#define get_hour(unix_time) ((unix_time % SECONDS_IN_DAY) / 3600)

    void generateData(const char* fileName, const size_t numRows) {
        std::ofstream testDataFile(fileName, std::ios::out | std::ios::trunc);
        if (!testDataFile.bad()) {
            for (size_t r = 0; r < numRows; ++r) {
                testDataFile << xorshift_utime << " " << rand_temp << "\n";
            }
            testDataFile.close();
        }
    }

    template <typename T, typename U>
    std::vector<std::pair<T, U>> parseData(const char* fileName) {
        std::vector<std::pair<T, U>> dataVec;
        std::ifstream dataFile(fileName, std::ios::in);
        if (!dataFile) {
            std::cerr << "Failed to parse data" << std::endl;
            exit(EXIT_FAILURE);
        }

        std::string line;
        while (std::getline(dataFile >> std::ws, line)) {
            std::istringstream istr(line);
            T time;
            U temp;
            istr >> time >> temp;
            dataVec.emplace_back(get_hour(time), temp);
        }
        return dataVec;
    }

    template <typename T, typename U>
    std::array<U, 24> workOutHistogram(const std::vector<std::pair<T, U>>& data) {
        std::array<std::pair<U, T>, 24> temp_hist{};
        for (const auto& [time, temp] : data) {
            temp_hist[time].first = temp_hist[time].first + temp;  // accumulate temperature
            temp_hist[time].second++;  // store up the measurements for the given hour
        }

        for (auto& [temp, num_measur] : temp_hist) {
            temp = temp / num_measur;  // calc average
        }

        std::array<U, 24> hist{};
        for (int32_t i = 0; i < 24; i++) {
            hist[i] = temp_hist[i].first;
        }

        return hist;
    }

    void formatOutput(const std::array<float, 24>& floatHist,
                      const std::array<fpt::Float, 24>& fptFloatHist) {
        const int32_t name_width = 15;
        const int32_t hour_width = 8;
        const std::string sep = " |";
        const int32_t total_width = hour_width + (name_width * 2) + (sep.size() * 3);
        const std::string line = sep + std::string(total_width - 1, '-') + '|';

        std::cout << line << '\n'
                  << sep << std::setw(hour_width) << "hour" << sep << std::setw(name_width)
                  << "float" << sep << std::setw(name_width) << "fpt::Float" << sep << "\n"
                  << line << std::endl;

        for (int32_t i = 0; i < 24; ++i) {
            std::cout << sep << std::setw(hour_width) << i << sep << std::setw(name_width)
                      << floatHist[i] << sep << std::setw(name_width) << fptFloatHist[i] << sep
                      << "\n"
                      << line << std::endl;
        }
    }

}  // namespace

int32_t main() {
    const char* fileName = "data.txt";
    const int32_t numMeasurements = 10'000;

    generateData(fileName, numMeasurements);

    const auto dataFloat = parseData<uint32_t, float>(fileName);
    const auto dataFptFloat = parseData<uint32_t, fpt::Float>(fileName);

    const auto floatHist = workOutHistogram(dataFloat);
    const auto fptFloatHist = workOutHistogram(dataFptFloat);

    formatOutput(floatHist, fptFloatHist);

    return 0;
}
