
// Program-wide defines
#define WRITE_OUTPUT_TO_QUERY true	// true -> outputfiles (result-XXXXX.txt) are created when queries are solved
#define USE_GPU false				// true -> OpenCL is used (DONT FLIP THIS, DOESNT WORK (yet))
#define USE_MULTITHREAD true		// false -> use only one thread, useful to debug concurrency issues
#define USE_FAST_IO true			// true -> file loading is faster, but less robust
#define ONLY_TOTAL_TIMES false		// true -> print diagnostic information
#define USE_FOPEN_S true			// true -> using windows file API


#define TRAJECTORY_FILES_OFFSET "" // directory appended to the load function, set to "" if the trajectory files are in the same folder as the executable
