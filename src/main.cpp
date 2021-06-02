#include <iostream>
#include <string>
#include "Communication.hpp"
#include "Case.hpp"

int main(int argn, char **args) {
    Params params;
    Communication::init_mpi(&params);
    if (argn > 1) {
        std::string file_name{args[1]};
        Case problem(file_name, argn, args, params);
        problem.simulate();

    } else {
        std::cout << "Error: No input file is provided to fluidchen." << std::endl;
        std::cout << "Example usage: /path/to/fluidchen /path/to/input_data.dat" << std::endl;
    }
    Communication::finalize();
}
