#include <complex>
#include <array>
#include <utility>
#include <unordered_map>
#define _USE_MATH_DEFINES
#include <cmath>

#ifndef FFT_H
#define FFT_H

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

    template <int N>
    using rigid_array = std_array<std::complex<double>, two_to_the_N(N)>;

    void precompute_twiddle_factors_helper(int N)
    {
        for (int k = 0; k < N; k++)
        {
            twiddle_factor_memoized(N, k);
        }

        if (N > 0)
        {
            precompute_twiddle_factors_helper(N - 1);
        }
    }

public:
    /**
     * Precompute the memoized twiddle factors; this is optional,
     * since the twiddle factors will be memozied upon the first
     * call to fft.
     */
    template <int N>
    void precompute_twiddle_factors()
    {
        precompute_twiddle_factors_helper(N);
    }

    /**
     * Compute the DFT of in_array having length a power of two,
     * the length being determined at compile time.
     */
    template <int N>
    rigid_array<N> fft(rigid_array<N> in_array)
    {
        }

    /**
     * Compute the inverse DFT of in_array having length a power of two,
     * the length being determined at compile time.
     */
    template <int N>
    rigid_array<N> ifft(rigid_array<N> in_array)
    {
    }
};

#endif
