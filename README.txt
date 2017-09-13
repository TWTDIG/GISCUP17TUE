CENTRAL IDEA:

We first preprocess each trajectory in two ways:
A. Insert its endpoints in a grid.
B. Simplify it progressively up to four times.

We answer each query (trajectory Q with threshold epsilon) as follows:
1. Using the grid of A, we find all trajectories whose endpoints lie within distance epsilon of Q’s endpoints.
2. For each found trajectory P, we establish whether P is within Fréchet distance epsilon of Q:
2a. The exact decision problem is solved for simplifications of P (from step B) and simplifications of Q.
2b. A greedy (“equal-time”) approach for P and Q.
2c. The exact decision problem is solved for P and Q.

In 2a and 2c, we use various optimizations to avoid processing the complete free-space diagram.



COMPILATION:

The program has no external dependencies, and can be compiled using:

g++ FrechetCompImpl.cpp -std=c++11 -lpthread



RUNNING:

The program computes its own bounding box, and accepts the dataset (first parameter) and queryset file (second parameter) as follows:

binaryname dataset.txt queryset.txt

If encountering any trouble with parsing, please update the "settings.h" file, setting "USE_FAST_IO" to FALSE.
