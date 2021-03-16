#include <complex>
#include <vector>
#include <utility>
#include <functional>
#include <unordered_map>
#include <memory>
#define _USE_MATH_DEFINES
#include <cmath>

using complex_vector = std::vector<std::complex<double>>;

// https://stackoverflow.com/questions/5889238/why-is-xor-the-default-way-to-combine-hashes
// note: not optimized for 64bit
class int_pair_hash
{
public:
    std::size_t operator()(std::pair<int, int> const &key) const
    {
        std::hash<int> hasher;
        std::size_t lhs = hasher(key.first);
        std::size_t rhs = hasher(key.second);
        lhs ^= rhs + 0x9e3779b9 + (lhs << 6) + (lhs >> 2);
        return lhs;
    }
};

class FFT
{
private:
    FFT(FFT const &);
    FFT &operator=(FFT const &);

    constexpr int two_to_the_power_of(int k)
    {
        int val = 1;
        for (int i = 0; i < k; i++)
        {
            val = 2 * val;
        }
        return val;
    }

    std::unique_ptr<std::unordered_map<std::pair<int, int>, std::complex<double>, int_pair_hash>> twiddle_factors;

    std::complex<double> twiddle_factor_memoized(int N, int k)
    {
        auto it = twiddle_factors->find(std::pair<int, int>(N, k));
        if (it == twiddle_factors->end())
        {
            std::complex<double> val = std::polar<double>(1, (-2 * M_PI * k) / N);
            (*twiddle_factors)[std::pair<int, int>(N, k)] = val;
            return val;
        }
        return it->second;
    }

    void precompute_twiddle_factors_helper(int N)
    {
        for (int k = 0; k < N; k++)
        {
            twiddle_factor_memoized(N, k);
        }

        if (N != 1)
        {
            precompute_twiddle_factors_helper(N / 2);
        }
    }

    complex_vector::const_iterator advance(complex_vector const &in_vector, complex_vector::const_iterator it, int step)
    {
        auto end = in_vector.end();
        for (int i = 0; it != end && i < step; i++)
        {
            it++;
        }
        return it;
    }

    std::unique_ptr<complex_vector> slice(complex_vector const &in_vector, complex_vector::const_iterator beg, complex_vector::const_iterator end, int step)
    {
        auto sliced = std::make_unique<complex_vector>();
        auto it = beg;
        do
        {
            sliced->push_back(*it);
        } while (advance(in_vector, it, step) < end);
        return sliced;
    }

    std::unique_ptr<complex_vector> concat(complex_vector const &front, complex_vector const &back)
    {
        auto result = std::make_unique<complex_vector>();
        auto end = front.end();
        for (auto it = front.begin(); it != end; it++)
        {
            result->push_back(*it);
        }
        end = back.end();
        for (auto it = back.begin(); it != end; it++)
        {
            result->push_back(*it);
        }
        return result;
    }

    bool is_a_power_of_two(int N)
    {
        if (N <= 0)
        {
            return false;
        }
        int Ntmp = N;
        int k = 0;
        while (Ntmp != 1)
        {
            Ntmp /= 2;
            k++;
        }
        return N == two_to_the_power_of(k);
    }

    complex_vector &conj_in_place(complex_vector &in_vector)
    {
        auto end = in_vector.end();
        for (auto it = in_vector.begin(); it < end; it++)
        {
            *it = std::conj(*it);
        }
        return in_vector;
    }

    complex_vector &divide_by_length_in_place(complex_vector &in_vector)
    {
        auto end = in_vector.end();
        int N = in_vector.size();
        for (auto it = in_vector.begin(); it < end; it++)
        {
            *it /= N;
        }
        return in_vector;
    }

    std::unique_ptr<complex_vector> fft_helper(complex_vector const &in_vector, int &err)
    {
        err = 0;
        int N = in_vector.size();
        if (N == 1)
        {
            auto out_vector = std::make_unique<complex_vector>();
            out_vector->push_back(in_vector[0]);
            return out_vector;
        }
        else
        {
            auto sliced = slice(in_vector, in_vector.begin(), in_vector.end(), 2);
            auto front = fft_helper(*sliced, err);
            if (err != 0)
            {
                auto out_vector = std::make_unique<complex_vector>();
                return out_vector;
            }

            sliced = slice(in_vector, in_vector.begin() + 1, in_vector.end(), 2);
            auto back = fft_helper(*sliced, err);
            if (err != 0)
            {
                auto out_vector = std::make_unique<complex_vector>();
                return out_vector;
            }

            auto front_end = front->end();
            auto back_end = back->end();
            int k = 0;

            for (auto it1 = front->begin(), it2 = back->begin();
                 it1 != front_end && it2 != back_end;
                 it1++, it2++, k++)
            {
                auto twiddle_factor = twiddle_factor_memoized(N, k);
                auto tmp = *it1;
                *it1 = tmp + twiddle_factor * (*it2);
                *it2 = tmp - twiddle_factor * (*it2);
            }

            auto out_vector = concat(*front, *back);
            return out_vector;
        }
    }

    std::unique_ptr<complex_vector> ifft_helper(complex_vector const &in_vector, int &err)
    {
        auto in_vector_cpy_ptr = std::make_unique<complex_vector>(in_vector);
        conj_in_place(*in_vector_cpy_ptr);
        auto fft_in_vector_cpy_ptr = fft_helper(*in_vector_cpy_ptr, err);
        conj_in_place(*fft_in_vector_cpy_ptr);
        divide_by_length_in_place(*fft_in_vector_cpy_ptr);
        return fft_in_vector_cpy_ptr;
    }

public:
    FFT() : twiddle_factors(std::make_unique<std::unordered_map<std::pair<int, int>, std::complex<double>, int_pair_hash>>()) {}

    FFT(FFT &&fft) : twiddle_factors(std::move(fft.twiddle_factors)) {}

    FFT &operator=(FFT &&fft)
    {
        if (this != &fft)
        {
            twiddle_factors = std::move(fft.twiddle_factors);
        }
        return *this;
    }

    int number_of_twiddle_factors()
    {
        return twiddle_factors->size();
    }

    void precompute_twiddle_factors(int N)
    {
        precompute_twiddle_factors_helper(N);
    }

    std::unique_ptr<complex_vector> fft(complex_vector const &in_vector, int &err)
    {
        if (!is_a_power_of_two(in_vector.size()))
        {
            err = 1;
            return std::make_unique<complex_vector>();
        }
        return fft_helper(in_vector, err);
    }

    std::unique_ptr<complex_vector> ifft(complex_vector const &in_vector, int &err)
    {
        if (!is_a_power_of_two(in_vector.size()))
        {
            err = 1;
            return std::make_unique<complex_vector>();
        }
        return ifft_helper(in_vector, err);
    }
};