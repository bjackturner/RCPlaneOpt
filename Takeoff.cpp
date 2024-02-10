#include <iostream>
#include <chrono>

int main() {
    // Start measuring time
    auto start_time = std::chrono::high_resolution_clock::now();
    int var;

    // This is a simple C++ program
    for (int i = 0; i <= 1000000; i++) {
        var = var + i;

    }

    printf("var = %d\n",var);

    // Stop measuring time
    auto end_time = std::chrono::high_resolution_clock::now();

    // Calculate the duration
    auto duration = std::chrono::duration_cast<std::chrono::microseconds>(end_time - start_time);

    // Print the duration in microseconds
    std::cout << "Total runtime: " << duration.count() << " microseconds" << std::endl;

    return 0;
}
