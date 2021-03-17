#include <vector>
#include <complex>
#include <memory>

#ifndef FFT_H
#define FFT_H

class FFT
{
private:
    FFT(FFT const &);
    FFT &operator=(FFT const &);

public:
    FFT &operator=(FFT &&);
    FFT(FFT &&);
    FFT();

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
     * Compute the DFT of in_vector having length a power of two.
     * 
     * err is 0 if success.
     */
    std::unique_ptr<std::vector<std::complex<double>>> fft(std::vector<std::complex<double>> const &in_vector, int &err);

    /**
     * Compute the inverse DFT of in_vector having length a power of two.
     * 
     * Current implementation is less efficient than fft() since ifft() uses fft() along with an additional copy.
     * 
     * err is 0 if success.
     */
    std::unique_ptr<std::vector<std::complex<double>>> ifft(std::vector<std::complex<double>> const &in_vector, int &err);

    /**
     * Compute the N by N 2-D DFT of in_vector, N being a power of two.
     * 
     * err is 0 if success.
     */
    std::unique_ptr<std::vector<std::complex<double>>> fft2(std::vector<std::complex<double>> const &in_vector, int &err);

    /**
     * Compute the N by N 2-D inverse DFT of in_vector, N being a power of two.
     * 
     * err is 0 if success.
     */
    std::unique_ptr<std::vector<std::complex<double>>> ifft2(std::vector<std::complex<double>> const &in_vector, int &err);

    /**
     * Compute the N by N by N 3-D DFT of in_vector, N being a power of two.
     * 
     * err is 0 if success.
     */
    std::unique_ptr<std::vector<std::complex<double>>> fft3(std::vector<std::complex<double>> const &in_vector, int &err);

    /**
     * Compute the N by N by N 3-D inverse DFT of in_vector, N being a power of two.
     * 
     * err is 0 if success.
     */
    std::unique_ptr<std::vector<std::complex<double>>> ifft3(std::vector<std::complex<double>> const &in_vector, int &err);
};

#endif
