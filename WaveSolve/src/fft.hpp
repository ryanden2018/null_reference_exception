#include <vector>
#include <complex>

#ifndef FFT_H
#define FFT_H

class FFT
{
public:
    /**
     * Returns the number of twiddle factors computed so far (for testing).
     */
    int number_of_twiddle_factors();

    /**
     * Precompute the memoized twiddle factors; this is optional,
     * since the twiddle factors will be memozied upon the first
     * call to fft.
     */
    void precompute_twiddle_factors(int N);

    /**
     * Compute the DFT of in_vector having length a power of two,
     * the length being determined at compile time. This method is not in-place.
     */
    std::vector<std::complex<double>> fft(std::vector<std::complex<double>> in_vector);

    /**
     * Compute the inverse DFT of in_vector having length a power of two,
     * the length being determined at compile time. This method is not in-place.
     */
    std::vector<std::complex<double>> ifft(std::vector<std::complex<double>> in_vector);
};

#endif
