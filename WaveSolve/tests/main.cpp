#include <iostream>

#include "fft.test.hpp"

int main()
{
    std::cout << "Test suite initiated." << std::endl;
    fft_test_main();
    std::cout << "Finished." << std::endl;

    return 0;
}