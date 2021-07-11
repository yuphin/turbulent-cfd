#include "Utilities.hpp"
#include <sstream>
#include <cmath>

std::vector<std::vector<int>> parse_geometry_file(std::string filedoc, int xdim, int ydim) {
    std::vector<std::vector<int>> geometry_data(xdim + 2, std::vector<int>(ydim + 2, 0));
    int numcols, numrows, depth;

    std::ifstream infile(filedoc);
    std::stringstream ss;
    std::string inputLine = "";

    // First line : version
    getline(infile, inputLine);
    if (inputLine.compare("P2") != 0) {
        std::cerr << "First line of the PGM file should be P2" << std::endl;
    }

    // Second line : comment
    getline(infile, inputLine);

    // Continue with a stringstream
    ss << infile.rdbuf();
    // Third line : size
    ss >> numrows >> numcols;
    // Fourth line : depth
    ss >> depth;

    // Following lines : data
    for (int col = numcols - 1; col > -1; --col) {
        for (int row = 0; row < numrows; ++row) {
            ss >> geometry_data[row][col];
        }
    }
    infile.close();
    return geometry_data;
}

std::vector<std::vector<int>> build_lid_driven_cavity(int xdim, int ydim) {
    std::vector<std::vector<int>> geometry_data(xdim + 2, std::vector<int>(ydim + 2, 0));
    for (int i = 0; i < xdim + 2; ++i) {
        for (int j = 0; j < ydim + 2; ++j) {
            // Bottom, left and right walls: no-slip
            if (i == 0 || j == 0 || i == xdim + 1) {
                geometry_data.at(i).at(j) = LidDrivenCavity::fixed_wall_id;
            }
            // Top wall: moving wall
            else if (j == ydim + 1) {
                geometry_data.at(i).at(j) = LidDrivenCavity::moving_wall_id;
            }
        }
    }
    return geometry_data;
}

std::vector<std::vector<int>> partition(const std::vector<std::vector<int>> &vec, int imin, int imax, int jmin,
                                        int jmax) {

    int xdim = imax - imin + 1;
    int ydim = jmax - jmin + 1;
   
    std::vector<std::vector<int>> result(xdim, std::vector<int>(ydim, 0));
    for (int i = 0; i < xdim; i++) {
        for (int j = 0; j < ydim; j++) {
            result[i][j] = vec[i + imin][j + jmin];
        }
    }
    return result;
}

std::vector<std::vector<int>> refine_geometry(const std::vector<std::vector<int>> &vec, int refine, int &imax, int &jmax)
{
    int iold = imax;
    int jold = jmax;

    Real refine_fac = std::pow(2, refine);
    imax = imax * refine_fac;
    jmax = jmax * refine_fac;

    int steps = std::pow(2, refine);
    std::vector<std::vector<int>> geometry_data(imax + 2, std::vector<int>(jmax + 2, 0));

    for (int i = 0; i < iold; i++) {
        for (int j = 0; j < jold; j++) {
            for (int ii = 0; ii < steps; ii++) {
                for (int jj = 0; jj < steps; jj++) {
                    geometry_data[i*steps + ii + 1][j*steps + jj + 1] = vec[i+1][j+1];
                }
            }
        }
    }

    // Fill domain boundaries
    for (int i = 0; i < iold; i++) {
        for (int ii = 0; ii < steps; ii++) {
            geometry_data[i*steps + ii + 1][0] = vec[i+1][0];
            geometry_data[i*steps + ii + 1][jmax + 1] = vec[i+1][jold + 1];
        }
    }

    for (int j = 0; j < jold; j++) {
        for (int jj = 0; jj < steps; jj++) {
            geometry_data[0][j*steps + jj + 1] = vec[0][j+1];    
            geometry_data[imax + 1][j*steps + jj + 1] = vec[iold + 1][j+1];
        }
    }
    
    geometry_data[0][0] = vec[0][0];
    geometry_data[0][jmax+1] = vec[0][jold+1];
    geometry_data[imax+1][0] = vec[iold+1][0];
    geometry_data[imax+1][jmax+1] = vec[iold+1][jold+1];

    return geometry_data;
}

bool is_inlet(int id) { return id >= 2 && id <= 9; }

bool is_outlet(int id) { return id == 1; }

bool is_no_slip(int id) { return id >= 10 && id <= 19; }

bool is_free_slip(int id) { return id >= 20; }
