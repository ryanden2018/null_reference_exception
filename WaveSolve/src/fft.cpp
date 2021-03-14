#include <complex>
#include <vector>
#include <utility>
#include <unordered_map>
#define _USE_MATH_DEFINES
#include <cmath>

using complex_vector = std::vector<std::complex<double>>;

class FFT
{
private:
    constexpr int two_to_the_N(int N)
    {
        int val = 1;
        for (int i = 0; i < N; i++)
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
            std::complex<double> val = std::polar<double>(1, (-2 * M_PI * k) / two_to_the_N(N));
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

        if (N != 0)
        {
            precompute_twiddle_factors_helper(N - 1);
        }
    }

    complex_vector::iterator advance(complex_vector in_vector, complex_vector::iterator it, int step)
    {
        for (int i = 0; it != in_vector.end() && i < step; i++)
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
        auto it = front.begin();
        while (it != front.end())
        {
            out_vector.push_back(*it);
        }
        it = back.begin();
        while (it != back.end())
        {
            out_vector.push_back(*it);
        }
    }

public:
    /**
     * Precompute the memoized twiddle factors; this is optional,
     * since the twiddle factors will be memozied upon the first
     * call to fft.
     */
    void precompute_twiddle_factors(int N)
    {
        precompute_twiddle_factors_helper(N);
    }

    /**
     * Compute the DFT of in_vector having length a power of two,
     * the length being determined at compile time. This method is not in-place.
     */
    complex_vector fft(complex_vector in_vector)
    {
        complex_vector out_vector;
        if (in_vector.size() == 1)
        {
            out_vector.push_back(in_vector[0]);
        }
        else
        {
            complex_vector front = fft(slice(in_vector, in_vector.begin(), in_vector.end(), 2));
            complex_vector back = fft(slice(in_vector, in_vector.begin() + 1, in_vector.end(), 2));

            for (auto it1 = front.begin(), it2 = back.begin(); it1 != front.end() && it2 != back.end(); it1++, it2++)
            {
            }

            auto it1 = front.begin();
            auto it2 = back.begin();
            while (it1 != front.end() && it2 != back.end())
            {

                it1++;
                it2++;
            }
        }
        return out_vector;
    }

    /**
     * Compute the inverse DFT of in_vector having length a power of two,
     * the length being determined at compile time. This method is not in-place.
     */
    complex_vector ifft(complex_vector in_vector)
    {
    }
};