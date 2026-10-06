// params.cpp -- read the run parameters from the text file params.txt (in the current directory) and print the set-up.
//
// params.txt holds whitespace-separated numbers in this order:
//   base_step  thermalization_steps  production_steps  measure_every  [adaptive_factor]  [seed]
// The last two are optional: adaptive_factor defaults to 0 (adaptive step off), seed defaults to 1.
#include "model.h"      // the run parameters and constants

void read_parameters() {   // read params.txt and print the set-up
    ifstream file("params.txt");                           // open the parameter file
    if (!file.good()) {                                    // not found: stop with a message
        cout << "\ncan't open params.txt\n";   // tell the user ...
        exit(1);   // ... and stop the program
    }

    file >> base_step >> thermalization_steps >> production_steps >> measure_every;   // the four required numbers
    if (!(file >> adaptive_factor)) adaptive_factor = 0.0;   // optional fifth number
    if (!(file >> seed)) seed = 1;                           // optional sixth number
    if (seed <= 0) seed = 1;                                 // srand(0) is not portable: use 1

    cout << "-- D=4 BOSONIC MATRIX MODEL --\n";   // print a header
    cout << "N = " << NCOLOR << ", L = " << L << ", d = " << D << ", T = " << TEMP << ", lambda = " << LAMBDA   // print N, L, d, T, lambda ...
         << ", lattice spacing a = " << SPACING << "\n";   // ... and the lattice spacing
    cout << "step " << base_step << ", thermalization " << thermalization_steps << ", production " << production_steps   // print the run parameters ...
         << ", measure every " << measure_every << ", adaptive factor " << adaptive_factor << ", seed " << seed << "\n\n";   // ... (end of the line)
}
