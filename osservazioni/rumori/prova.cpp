#include <algorithm>
#include <cmath>
#include <fstream>
#include <iostream>
#include <numeric>
#include <sstream>
#include <stdexcept>
#include <string>
#include <vector>
#include <iomanip>

bool conversion(std::string num) {
  try {
    std::stod(num);
    return true;
  } catch (const std::exception& e) {
    return false;
  }
}

std::vector<std::string> read_data(std::string file) {
  std::vector<std::string> data;
  std::string line;
  std::ifstream in{file};
  while (std::getline(in, line, '\n')) {
    std::stringstream line_stream(line);
    std::string cell;
    while (std::getline(line_stream, cell, '\t')) data.push_back(cell);
    if (!line_stream && cell.empty()) data.push_back(" ");
  }
  in.close();
  return data;
}
std::vector<double> collect_data(std::vector<std::string>& data) {
  std::vector<double> d_vector;
  int counter{0};
  for (int i = 0; i < data.size(); i++) {
    if (conversion(data[i]) == true) {
      double num = std::stod(data[i]);
      d_vector.push_back(num);
    } else {
      counter++;
      continue;
    }
  }
  int N = d_vector.size();
  std::cout << "Converted: " << N << "\tNot converted: " << counter << '\n';
  return d_vector;
}
double StdDev(std::vector<double>& d_vector) {
  int N = d_vector.size();
  double variance = 0.;
  double mean = std::accumulate(d_vector.begin(), d_vector.end(), 0.) / N;
  for (int i = 0; i < N; i++) {
    variance += (d_vector[i] - mean) * (d_vector[i] - mean);
  }
  return std::sqrt(variance / (N - 1));
}

double min_max(std::vector<double>& d_vector) {
  double max = *std::max_element(d_vector.begin(), d_vector.end());
  double min = *std::min_element(d_vector.begin(), d_vector.end());
  double result = (max - min) / std::sqrt(12);
  std::cout << result / std::sqrt(12) << '\n';
  return result;
}

int main(int argc, const char* argv[]) {
  try {
    if (argc != 2)
      throw std::runtime_error("Invalid argument\n");
    else {
      std::string name = argv[1];
      std::cout << name << '\n';
      std::vector<std::string> data = read_data(name);
      std::vector<double> numbers = collect_data(data);
      std::cout << std::accumulate(numbers.begin(), numbers.end(), 0.) /
                       numbers.size()
                << ' ' << StdDev(numbers) << '\n';
      std::cout << std::setfill('-') << std::setw(40) << '\n';
      // std::cout << min_max(numbers) << '\n';
    }
  } catch (std::exception& error) {
    std::cout << error.what() << '\n';
  }
}
