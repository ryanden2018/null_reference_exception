#include <complex>
#include <vector>
#include <utility>
#include <unordered_map>
#include <stdexcept>
#define _USE_MATH_DEFINES
#include <cmath>

using complex_vector = std::vector<std::complex<double>>;

class FFT
{
private:
    constexpr int two_to_the_power_of(int k)
    {
        int val = 1;
        for (int i = 0; i < k; i++)
        {
            val = 2 * val;
        }
        return val;
    }

    std::unordered_map<std::pair<int, int>, std::complex<double>> twiddle_factors;

    std::complex<double> twiddle_factor_memoized(int N, int k)
    {
        auto it = twiddle_factors.find(std::pair<int, int>(N, k));
        if (it == twiddle_factors.end())
        {
            std::complex<double> val = std::polar<double>(1, (-2 * M_PI * k) / N);
            twiddle_factors[std::pair<int, int>(N, k)] = val;
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

    complex_vector::iterator advance(complex_vector in_vector, complex_vector::iterator it, int step)
    {
        auto end = in_vector.end();
        for (int i = 0; it != end && i < step; i++)
        {
            it++;
        }
        return it;
    }

    complex_vector slice(complex_vector in_vector, complex_vector::iterator beg, complex_vector::iterator end, int step)
    {
        complex_vector sliced;
        auto it = beg;
        do
        {
            sliced.push_back(*it);
        } while (advance(in_vector, it, step) < end);
        return sliced;
    }

    complex_vector concat(complex_vector front, complex_vector back)
    {
        complex_vector out_vector;
        auto end = front.end();
        for (auto it = front.begin(); it != end; it++)
        {
            out_vector.push_back(*it);
        }
        end = back.end();
        for (auto it = back.begin(); it != end; it++)
        {
            out_vector.push_back(*it);
        }
        return out_vector;
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

    complex_vector conj_in_place(complex_vector in_vector)
    {
        auto end = in_vector.end();
        for (auto it = in_vector.begin(); it < end; it++)
        {
            *it = std::conj(*it);
        }
        return in_vector;
    }

    complex_vector divide_by_length_in_place(complex_vector in_vector)
    {
        auto end = in_vector.end();
        int N = in_vector.size();
        for (auto it = in_vector.begin(); it < end; it++)
        {
            *it /= N;
        }
        return in_vector;
    }

    complex_vector fft_helper(complex_vector in_vector)
    {
        int N = in_vector.size();
        complex_vector out_vector;
        if (N == 1)
        {
            out_vector.push_back(in_vector[0]);
        }
        else
        {
            complex_vector front = fft_helper(slice(in_vector, in_vector.begin(), in_vector.end(), 2));
            complex_vector back = fft_helper(slice(in_vector, in_vector.begin() + 1, in_vector.end(), 2));
            auto front_end = front.end();
            auto back_end = back.end();
            int k = 0;

            for (auto it1 = front.begin(), it2 = back.begin();
                 it1 != front_end && it2 != back_end;
                 it1++, it2++, k++)
            {
                auto twiddle_factor = twiddle_factor_memoized(k, N);
                auto tmp = *it1;
                *it1 = tmp + twiddle_factor * (*it2);
                *it2 = tmp - twiddle_factor * (*it2);
            }

            out_vector = concat(front, back);
        }
        return out_vector;
    }

    complex_vector ifft_helper(complex_vector in_vector)
    {
        return divide_by_length_in_place(
            conj_in_place(
                fft_helper(
                    conj_in_place(in_vector))));
    }

public:
    int number_of_twiddle_factors()
    {
        return twiddle_factors.size();
    }

    void precompute_twiddle_factors(int N)
    {
        precompute_twiddle_factors_helper(N);
    }

    complex_vector fft(complex_vector in_vector)
    {
        if (!is_a_power_of_two(in_vector.size()))
        {
            throw std::invalid_argument("Only DFT of inputs of size a power of two are supported.");
        }
        return fft_helper(in_vector);
    }

    complex_vector ifft(complex_vector in_vector)
    {
        if (!is_a_power_of_two(in_vector.size()))
        {
            throw std::invalid_argument("Only IDFT of inputs of size a power of two are supported.");
        }
        return ifft_helper(in_vector);
    }
};