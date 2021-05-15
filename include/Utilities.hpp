#pragma once

#include "Case.hpp"

#include <iostream>
#include <iomanip>
#include <fstream>
#include <string>
#include <filesystem>

/**
 * @brief Class to handle logging functionality.
 *
 */
class Logger
{
    public:
        Logger() {};

        // Check flag for creating a log file
        void parseFlag(char* arg)
        {
            std::string arg_string = static_cast<std::string>(arg);
            if (arg_string == "-log") {
                _log = true;
            }
            else {
                std::cout << "Second argument invalid - ignored..." << std::endl;
            }
        }

        // Prepare the log file
        void createLog(const std::string &dict, const std::string &name) {
            if (_log) {
                #ifdef WIN32
                char sep = '\\';
                #else 
                char sep = '/';
                #endif
            std::string fileName = dict + sep + name + ".log";
            logFile.open(fileName);
            logFile << "#### " << name << " ####\n" << std::endl;
            }
        }
        
        // Print a progress bar
        void progressBar(double t, double t_end) 
        {
            double progress = t / t_end;
            int pos = static_cast<int>(progress * _barWidth);
            std::cout << " Time: |";
            for (int i = 0; i < _barWidth; ++i) {
                #ifndef WIN32
                if (i <= pos) std::cout << "\u25A0";
                else std::cout << "\u25A1";
                #else
                if (i <= pos) std::cout << "=";
                else std::cout << "-";
                #endif
            }
            std::cout << "|";
            std::cout << std::setw(5) << std::fixed << std::setprecision(2) << t;
            std::cout << " / " << t_end << "\r";
            std::cout << std::flush;
        }

        // Set _maxItersReached 
        void maxIterWarning() {
            _maxItersReached = true;
        }

        // Write information of current timestep into logFile 
        void writeLog(uint32_t timestep, double t, int it, int max_iter, double res) {
            if (_log){
                logFile << "Timestep " << timestep << ": " << std::endl;
                logFile << "Simulation time: t = " << t << std::endl;
                logFile << "SOR final residual: " << res << "  No iterations: " << it << std::endl;
                logFile << std::endl;
            }
        }

        // Final output at the end of the simulation
        void finish() {
            std::cout << "\nSimulation finished!" << std::endl;
            if (_maxItersReached) {
                std::cout << "WARNING: Maximum number of iterations was reached in some timesteps" << std::endl;
            }
        }

        // Close logFile if it was opened
        ~Logger() {
            if (logFile.is_open()) {logFile.close();}
        }

    private:
        // True if log file is being created
        bool _log = false;
        // Bar width of the progress bar
        const int _barWidth = 50;
        // Remember if maximum number of iterations was reached during SOR
        bool _maxItersReached = false;
        // log File object
        std::ofstream logFile;
};