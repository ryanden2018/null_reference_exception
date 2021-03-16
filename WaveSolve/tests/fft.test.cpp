#include <iostream>
#include <vector>
#include <complex>

#include "../src/fft.hpp"

void fft_power_of_two_constraint()
{
    FFT fft;
    int err;
    int N = 300;
    std::vector<std::complex<double>> in_vector(N, 0);
    fft.fft(in_vector, err);
    if (err == 0)
    {
        std::cout << "ERROR (fft_power_of_two_constraint): power of two constraint must be enforced" << std::endl;
    }
}

void ifft_power_of_two_constraint()
{
    FFT fft;
    int err;
    int N = 300;
    std::vector<std::complex<double>> in_vector(N, 0);
    fft.ifft(in_vector, err);
    if (err == 0)
    {
        std::cout << "ERROR (ifft_power_of_two_constraint): power of two constraint must be enforced" << std::endl;
    }
}

void twiddle_factors_are_memoized()
{
    FFT fft;
    int N = 300;
    auto n0 = fft.number_of_twiddle_factors();
    fft.precompute_twiddle_factors(N);
    auto n1 = fft.number_of_twiddle_factors();
    fft.precompute_twiddle_factors(N);
    auto n2 = fft.number_of_twiddle_factors();
    if (n0 != 0 || n1 == 0 || n1 != n2)
    {
        std::cout << "ERROR (twiddle_factors_are_memoized): twiddle factors must be memoized" << std::endl;
    }
}

void computes_fft()
{
}

void ifft_inverses_fft()
{
}

void computes_fft2()
{
    std::cout << "computes_fft2 test not yet implemented" << std::endl;
}

void ifft2_inverses_fft2()
{
    std::cout << "ifft2_inverses_fft2 test not yet implemented" << std::endl;
}

void computes_fft3()
{
    std::cout << "computes_fft3 test not yet implemented" << std::endl;
}

void ifft3_inverses_fft3()
{
    std::cout << "ifft3_inverses_fft3 test not yet implemented" << std::endl;
}

void fft_test_main()
{
    fft_power_of_two_constraint();
    ifft_power_of_two_constraint();
    twiddle_factors_are_memoized();
    computes_fft();
    ifft_inverses_fft();
    computes_fft2();
    ifft2_inverses_fft2();
    computes_fft3();
    ifft3_inverses_fft3();
}
