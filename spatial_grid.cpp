#include <pybind11/pybind11.h>
#include <pybind11/stl.h>
#include <pybind11/numpy.h>
#include <vector>
#include <cmath>
#include <numeric>
#include <algorithm>
#include <string>
#include <omp.h>
#include <filesystem>

namespace py = pybind11;
namespace fs = std::filesystem;

void verify_resources_cpp(const std::string& addon_path) {
    fs::path icons_dir = fs::path(addon_path) / "icons";
    fs::path icon1_path = icons_dir / "icon_amoa_materialtransfer.svg";
    fs::path icon2_path = icons_dir / "icon_amoa_objectalignment.svg";

    if (!fs::exists(icon1_path) || !fs::exists(icon2_path)) {
        throw std::runtime_error("AMOA license verification failed. Please use an official version.");
    }
}

struct Vector3 {
    double x, y, z;
};

struct ObjectProperties {
    double volume = 0.0;
    double surface_area = 0.0;
    double curvature = 0.0;
    int vertex_count = 0;
    int polygon_count = 0;
    Vector3 scale = {1.0, 1.0, 1.0};
    Vector3 dimensions = {0.0, 0.0, 0.0};
    std::vector<double> spatial_distribution;
    std::vector<double> edge_lengths;
};

double calculate_base_similarity(const ObjectProperties& source, const ObjectProperties& target, const std::string& addon_path) {
    verify_resources_cpp(addon_path);

    double vol_max = std::max(source.volume, target.volume);
    double volume_ratio = (vol_max > 1e-9) ? std::min(source.volume, target.volume) / vol_max : 1.0;

    double area_max = std::max(source.surface_area, target.surface_area);
    double area_ratio = (area_max > 1e-9) ? std::min(source.surface_area, target.surface_area) / area_max : 1.0;

    double curv_diff = std::abs(source.curvature - target.curvature);
    double curv_norm = std::max({source.curvature, target.curvature, 1.0});
    double curv_ratio = 1.0 - (curv_diff / curv_norm);

    double vert_max = static_cast<double>(std::max(source.vertex_count, target.vertex_count));
    double vert_ratio = (vert_max > 0) ? static_cast<double>(std::min(source.vertex_count, target.vertex_count)) / vert_max : 1.0;

    double poly_max = static_cast<double>(std::max(source.polygon_count, target.polygon_count));
    double poly_ratio = (poly_max > 0) ? static_cast<double>(std::min(source.polygon_count, target.polygon_count)) / poly_max : 1.0;
    
    double dim_x_max = std::max(source.dimensions.x, target.dimensions.x);
    double dim_y_max = std::max(source.dimensions.y, target.dimensions.y);
    double dim_z_max = std::max(source.dimensions.z, target.dimensions.z);
    double dim_x_ratio = (dim_x_max > 1e-9) ? std::min(source.dimensions.x, target.dimensions.x) / dim_x_max : 1.0;
    double dim_y_ratio = (dim_y_max > 1e-9) ? std::min(source.dimensions.y, target.dimensions.y) / dim_y_max : 1.0;
    double dim_z_ratio = (dim_z_max > 1e-9) ? std::min(source.dimensions.z, target.dimensions.z) / dim_z_max : 1.0;
    double bbox_ratio = (dim_x_ratio + dim_y_ratio + dim_z_ratio) / 3.0;

    double w_vol = 0.40;
    double w_area = 0.25;
    double w_bbox = 0.20;
    double w_vert_poly = 0.10;
    double w_curv = 0.05;

    double score = volume_ratio * w_vol +
                   area_ratio * w_area +
                   bbox_ratio * w_bbox +
                   ((vert_ratio + poly_ratio) / 2.0) * w_vert_poly +
                   curv_ratio * w_curv;

    return std::max(0.0, std::min(1.0, score));
}

double calculate_spatial_similarity(const ObjectProperties& source, const ObjectProperties& target) {
    const auto& dist1 = source.spatial_distribution;
    const auto& dist2 = target.spatial_distribution;

    if (dist1.empty() || dist1.size() != dist2.size()) {
        return 0.0;
    }

    size_t n = dist1.size();
    double sum1 = 0.0, sum2 = 0.0, sum1_sq = 0.0, sum2_sq = 0.0, p_sum = 0.0;

    for (size_t i = 0; i < n; ++i) {
        sum1 += dist1[i];
        sum2 += dist2[i];
        sum1_sq += dist1[i] * dist1[i];
        sum2_sq += dist2[i] * dist2[i];
        p_sum += dist1[i] * dist2[i];
    }

    double num = p_sum - (sum1 * sum2 / n);
    double den = sqrt((sum1_sq - sum1 * sum1 / n) * (sum2_sq - sum2 * sum2 / n));

    if (den < 1e-9) {
        return 1.0;
    }

    return std::max(0.0, std::min(1.0, (num / den + 1.0) / 2.0));
}

std::vector<double> generate_spatial_grid(py::array_t<double> vertices, int resolution, const Vector3& min_bound, const Vector3& max_bound) {
    py::buffer_info buf = vertices.request();
    if (buf.ndim != 2 || buf.shape[1] != 3) {
        throw std::runtime_error("Vertices must be a Nx3 NumPy array");
    }
    double* ptr = static_cast<double*>(buf.ptr);
    size_t num_vertices = buf.shape[0];

    std::vector<long long> grid(resolution * resolution * resolution, 0);
    Vector3 size = {max_bound.x - min_bound.x, max_bound.y - min_bound.y, max_bound.z - min_bound.z};
    
    if (size.x < 1e-9) size.x = 1.0;
    if (size.y < 1e-9) size.y = 1.0;
    if (size.z < 1e-9) size.z = 1.0;

    #pragma omp parallel for
    for (long long i = 0; i < num_vertices; ++i) {
        double vx = ptr[i * 3 + 0];
        double vy = ptr[i * 3 + 1];
        double vz = ptr[i * 3 + 2];

        int ix = static_cast<int>(((vx - min_bound.x) / size.x) * resolution);
        int iy = static_cast<int>(((vy - min_bound.y) / size.y) * resolution);
        int iz = static_cast<int>(((vz - min_bound.z) / size.z) * resolution);

        ix = std::max(0, std::min(resolution - 1, ix));
        iy = std::max(0, std::min(resolution - 1, iy));
        iz = std::max(0, std::min(resolution - 1, iz));
        
        int index = ix * resolution * resolution + iy * resolution + iz;

        #pragma omp atomic
        grid[index]++;
    }

    std::vector<double> normalized_grid(grid.size());
    if (num_vertices > 0) {
        for (size_t i = 0; i < grid.size(); ++i) {
            normalized_grid[i] = static_cast<double>(grid[i]) / num_vertices;
        }
    }
    return normalized_grid;
}

py::array_t<double> calculate_similarity_matrix(const std::vector<ObjectProperties>& all_props, const std::string& addon_path) {
    verify_resources_cpp(addon_path);
    size_t n = all_props.size();
    py::array_t<double> result({n, n});
    py::buffer_info buf = result.request();
    double* ptr = static_cast<double*>(buf.ptr);

    #pragma omp parallel for collapse(2)
    for (long long i = 0; i < n; ++i) {
        for (long long j = 0; j < n; ++j) {
            if (i == j) {
                ptr[i * n + j] = 1.0;
            } else if (i < j) {
                double base_sim = calculate_base_similarity(all_props[i], all_props[j], addon_path);
                double spatial_sim = calculate_spatial_similarity(all_props[i], all_props[j]);
                
                double w_base = 0.7;
                double w_spatial = 0.3;
                double final_score = base_sim * w_base + spatial_sim * w_spatial;

                ptr[i * n + j] = final_score;
                ptr[j * n + i] = final_score;
            }
        }
    }
    return result;
}

PYBIND11_MODULE(amoa_cpp, m) {
    m.doc() = "C++ acceleration module for AMOA";

    py::class_<Vector3>(m, "Vector3")
        .def(py::init<>())
        .def_readwrite("x", &Vector3::x)
        .def_readwrite("y", &Vector3::y)
        .def_readwrite("z", &Vector3::z);

    py::class_<ObjectProperties>(m, "ObjectProperties")
        .def(py::init<>())
        .def_readwrite("volume", &ObjectProperties::volume)
        .def_readwrite("surface_area", &ObjectProperties::surface_area)
        .def_readwrite("curvature", &ObjectProperties::curvature)
        .def_readwrite("vertex_count", &ObjectProperties::vertex_count)
        .def_readwrite("polygon_count", &ObjectProperties::polygon_count)
        .def_readwrite("scale", &ObjectProperties::scale)
        .def_readwrite("dimensions", &ObjectProperties::dimensions)
        .def_readwrite("spatial_distribution", &ObjectProperties::spatial_distribution)
        .def_readwrite("edge_lengths", &ObjectProperties::edge_lengths);

    m.def("calculate_base_similarity", &calculate_base_similarity, "Calculates base similarity score between two objects",
        py::arg("source_props"), py::arg("target_props"), py::arg("addon_path"));

    m.def("calculate_spatial_similarity", &calculate_spatial_similarity, "Calculates spatial distribution similarity",
        py::arg("source_props"), py::arg("target_props"));
        
    m.def("generate_spatial_grid", &generate_spatial_grid, "Generates a spatial grid distribution from vertices",
        py::arg("vertices"), py::arg("resolution"), py::arg("min_bound"), py::arg("max_bound"));

    m.def("calculate_similarity_matrix", &calculate_similarity_matrix, "Calculates the full similarity matrix for a list of objects",
        py::arg("all_props"), py::arg("addon_path"));
}