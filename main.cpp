#include <iostream>
#include <numeric>

#include "ripser.h"


enum file_format {
	LOWER_DISTANCE_MATRIX,
	UPPER_DISTANCE_MATRIX,
	DISTANCE_MATRIX,
	POINT_CLOUD,
	DIPHA,
	SPARSE,
	BINARY
};

static const uint16_t endian_check(0xff00);
static const bool is_big_endian = *reinterpret_cast<const uint8_t*>(&endian_check);

bool is_prime(const coefficient_t n) {
	if (!(n & 1) || n < 2) return n == 2;
	for (coefficient_t p = 3; p <= n / p; p += 2)
		if (!(n % p)) return false;
	return true;
}

template <typename T> T read(std::istream& input_stream) {
	T result;
	char* p = reinterpret_cast<char*>(&result);
	if (input_stream.read(p, sizeof(T)).gcount() != sizeof(T)) return T();
	if (is_big_endian) std::reverse(p, p + sizeof(T));
	return result;
}

compressed_lower_distance_matrix read_point_cloud(std::istream& input_stream) {
	std::vector<std::vector<value_t>> points;

	std::string line;
	value_t value;
	while (std::getline(input_stream, line)) {
		std::vector<value_t> point;
		std::istringstream s(line);
		while (s >> value) {
			point.push_back(value);
			s.ignore();
		}
		if (!point.empty()) points.push_back(point);
		assert(point.size() == points.front().size());
	}

	euclidean_distance_matrix eucl_dist(std::move(points));
	index_t n = eucl_dist.size();
	std::cout << "point cloud with " << n << " points in dimension "
	          << eucl_dist.points.front().size() << std::endl;

	std::vector<value_t> distances;
	for (int i = 0; i < n; ++i)
		for (int j = 0; j < i; ++j) distances.push_back(eucl_dist(i, j));

	return compressed_lower_distance_matrix(std::move(distances));
}

sparse_distance_matrix read_sparse_distance_matrix(std::istream& input_stream) {
	std::vector<std::vector<index_diameter_t>> neighbors;
	index_t num_edges = 0;

	std::string line;
	while (std::getline(input_stream, line)) {
		std::istringstream s(line);
		size_t i, j;
		value_t value;
		s >> i;
		s >> j;
		s >> value;
		if (i != j) {
			neighbors.resize(std::max({neighbors.size(), i + 1, j + 1}));
			neighbors[i].push_back({j, value});
			neighbors[j].push_back({i, value});
			++num_edges;
		}
	}

	for (size_t i = 0; i < neighbors.size(); ++i)
		std::sort(neighbors[i].begin(), neighbors[i].end());

	return sparse_distance_matrix(std::move(neighbors), num_edges);
}

compressed_lower_distance_matrix read_lower_distance_matrix(std::istream& input_stream) {
	std::vector<value_t> distances;
	value_t value;
	while (input_stream >> value) {
		distances.push_back(value);
		input_stream.ignore();
	}

	return compressed_lower_distance_matrix(std::move(distances));
}

compressed_lower_distance_matrix read_upper_distance_matrix(std::istream& input_stream) {
	std::vector<value_t> distances;
	value_t value;
	while (input_stream >> value) {
		distances.push_back(value);
		input_stream.ignore();
	}

	return compressed_lower_distance_matrix(compressed_upper_distance_matrix(std::move(distances)));
}

compressed_lower_distance_matrix read_distance_matrix(std::istream& input_stream) {
	std::vector<value_t> distances;

	std::string line;
	value_t value;
	for (int i = 0; std::getline(input_stream, line); ++i) {
		std::istringstream s(line);
		for (int j = 0; j < i && s >> value; ++j) {
			distances.push_back(value);
			s.ignore();
		}
	}

	return compressed_lower_distance_matrix(std::move(distances));
}

compressed_lower_distance_matrix read_dipha(std::istream& input_stream) {
	if (read<int64_t>(input_stream) != 8067171840) {
		std::cerr << "input is not a Dipha file (magic number: 8067171840)" << std::endl;
		exit(-1);
	}

	if (read<int64_t>(input_stream) != 7) {
		std::cerr << "input is not a Dipha distance matrix (file type: 7)" << std::endl;
		exit(-1);
	}

	index_t n = read<int64_t>(input_stream);

	std::vector<value_t> distances;

	for (int i = 0; i < n; ++i)
		for (int j = 0; j < n; ++j)
			if (i > j)
				distances.push_back(read<double>(input_stream));
			else
				read<double>(input_stream);

	return compressed_lower_distance_matrix(std::move(distances));
}

compressed_lower_distance_matrix read_binary(std::istream& input_stream) {
	std::vector<value_t> distances;
	while (!input_stream.eof()) distances.push_back(read<value_t>(input_stream));
	return compressed_lower_distance_matrix(std::move(distances));
}

compressed_lower_distance_matrix read_file(std::istream& input_stream, const file_format format) {
	switch (format) {
	case LOWER_DISTANCE_MATRIX:
		return read_lower_distance_matrix(input_stream);
	case UPPER_DISTANCE_MATRIX:
		return read_upper_distance_matrix(input_stream);
	case DISTANCE_MATRIX:
		return read_distance_matrix(input_stream);
	case POINT_CLOUD:
		return read_point_cloud(input_stream);
	case DIPHA:
		return read_dipha(input_stream);
	default:
		return read_binary(input_stream);
	}
}

void print_usage_and_exit(int exit_code) {
	std::cerr
	    << "Usage: "
	    << "ripser "
	    << "[options] [filename]" << std::endl
	    << std::endl
	    << "Options:" << std::endl
	    << std::endl
	    << "  --help           print this screen" << std::endl
	    << "  --format         use the specified file format for the input. Options are:"
	    << std::endl
	    << "                     lower-distance (lower triangular distance matrix; default)"
	    << std::endl
	    << "                     upper-distance (upper triangular distance matrix)" << std::endl
	    << "                     distance       (full distance matrix)" << std::endl
	    << "                     point-cloud    (point cloud in Euclidean space)" << std::endl
	    << "                     dipha          (distance matrix in DIPHA file format)" << std::endl
	    << "                     sparse         (sparse distance matrix in sparse triplet format)"
	    << std::endl
	    << "                     binary         (lower triangular distance matrix in binary format)"
	    << std::endl
	    << "  --dim <k>        compute persistent homology up to dimension k" << std::endl
	    << "  --threshold <t>  compute Rips complexes up to diameter t" << std::endl
#ifdef USE_COEFFICIENTS
	    << "  --modulus <p>    compute homology with coefficients in the prime field Z/pZ"
	    << std::endl
#endif
#ifndef USING_SERIAL
	    << "  --threads <t>    number of threads to use"
	    << std::endl
#endif
	    << "  --ratio <r>      only show persistence pairs with death/birth ratio > r" << std::endl
	    << std::endl;
	exit(exit_code);
}

int main(int argc, char** argv) {
	const char* filename = nullptr;

	file_format format = DISTANCE_MATRIX;

	index_t dim_max = 1;
	value_t threshold = std::numeric_limits<value_t>::max();
	float ratio = 1;
	coefficient_t modulus = 2;
	unsigned num_threads = 0;

	for (index_t i = 1; i < argc; ++i) {
		const std::string arg(argv[i]);
		if (arg == "--help") {
			print_usage_and_exit(0);
		} else if (arg == "--dim") {
			std::string parameter = std::string(argv[++i]);
			size_t next_pos;
			dim_max = std::stol(parameter, &next_pos);
			if (next_pos != parameter.size()) print_usage_and_exit(-1);
		} else if (arg == "--threshold") {
			std::string parameter = std::string(argv[++i]);
			size_t next_pos;
			threshold = std::stof(parameter, &next_pos);
			if (next_pos != parameter.size()) print_usage_and_exit(-1);
		} else if (arg == "--ratio") {
			std::string parameter = std::string(argv[++i]);
			size_t next_pos;
			ratio = std::stof(parameter, &next_pos);
			if (next_pos != parameter.size()) print_usage_and_exit(-1);
		} else if (arg == "--format") {
			std::string parameter = std::string(argv[++i]);
			if (parameter.rfind("lower", 0) == 0)
				format = LOWER_DISTANCE_MATRIX;
			else if (parameter.rfind("upper", 0) == 0)
				format = UPPER_DISTANCE_MATRIX;
			else if (parameter.rfind("dist", 0) == 0)
				format = DISTANCE_MATRIX;
			else if (parameter.rfind("point", 0) == 0)
				format = POINT_CLOUD;
			else if (parameter == "dipha")
				format = DIPHA;
			else if (parameter == "sparse")
				format = SPARSE;
			else if (parameter == "binary")
				format = BINARY;
			else
				print_usage_and_exit(-1);
#ifdef USE_COEFFICIENTS
		} else if (arg == "--modulus") {
			std::string parameter = std::string(argv[++i]);
			size_t next_pos;
			modulus = std::stol(parameter, &next_pos);
			if (next_pos != parameter.size() || !is_prime(modulus)) print_usage_and_exit(-1);
#endif
#ifndef USING_SERIAL
		} else if (arg == "--threads") {
			std::string parameter = std::string(argv[++i]);
			size_t next_pos;
			num_threads = std::stol(parameter, &next_pos);
			if (next_pos != parameter.size() || !is_prime(modulus)) print_usage_and_exit(-1);
#endif
		} else {
			if (filename) { print_usage_and_exit(-1); }
			filename = argv[i];
		}
	}

#ifdef USE_TBB
	tbb::task_scheduler_init init(num_threads);
#endif

	std::ifstream file_stream(filename);
	if (filename && file_stream.fail()) {
		std::cerr << "couldn't open file " << filename << std::endl;
		exit(-1);
	}

	if (format == SPARSE) {
		sparse_distance_matrix dist =
		    read_sparse_distance_matrix(filename ? file_stream : std::cin);
		std::cout << "sparse distance matrix with " << dist.size() << " points and "
		          << dist.num_edges << "/" << (dist.size() * (dist.size() - 1)) / 2 << " entries"
		          << std::endl;

		ripser<sparse_distance_matrix>(std::move(dist), dim_max, threshold, ratio, modulus, num_threads)
		    .compute_barcodes();
	} else {
		compressed_lower_distance_matrix dist =
		    read_file(filename ? file_stream : std::cin, format);

		value_t min = std::numeric_limits<value_t>::infinity(),
		        max = -std::numeric_limits<value_t>::infinity(), max_finite = max;
		int num_edges = 0;

		if (threshold == std::numeric_limits<value_t>::max()) {
			value_t enclosing_radius = std::numeric_limits<value_t>::infinity();
			for (size_t i = 0; i < dist.size(); ++i) {
				value_t r_i = -std::numeric_limits<value_t>::infinity();
				for (size_t j = 0; j < dist.size(); ++j) r_i = std::max(r_i, dist(i, j));
				enclosing_radius = std::min(enclosing_radius, r_i);
			}
			threshold = enclosing_radius;
		}

		for (auto d : dist.distances) {
			min = std::min(min, d);
			max = std::max(max, d);
			max_finite =
			    d != std::numeric_limits<value_t>::infinity() ? std::max(max, d) : max_finite;
			if (d <= threshold) ++num_edges;
		}
		std::cout << "value range: [" << min << "," << max_finite << "]" << std::endl;

		if (threshold >= max) {
			std::cout << "distance matrix with " << dist.size() << " points" << std::endl;
			ripser<compressed_lower_distance_matrix>(std::move(dist), dim_max, threshold, ratio,
			                                         modulus, num_threads)
			    .compute_barcodes();
		} else {
			std::cout << "sparse distance matrix with " << dist.size() << " points and "
			          << num_edges << "/" << (dist.size() * dist.size() - 1) / 2 << " entries"
			          << std::endl;

			ripser<sparse_distance_matrix>(sparse_distance_matrix(std::move(dist), threshold),
			                               dim_max, threshold, ratio, modulus, num_threads)
			    .compute_barcodes();
		}
		exit(0);
	}
}
