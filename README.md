This repository contains our implementations of the sorting algorithms described in the paper "Noisy Sorting Without Searching: Data Oblivious
Sorting with Comparison Errors".

To compile, run
```make```

Then, for example, to generate and sort a random sequence of 32768 elements using our Window-Odd-Even-Sort algorithm, with a comparison error probability of 0.03, run
```./main random windowOddEven 32768 0.03```,
which will store the results in the directory ```input_random_algorithm_windowOddEven_n_32768_p_0.03```.
