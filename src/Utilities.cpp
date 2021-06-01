#include "Utilities.hpp"
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