// [RESIZABLE] Fail-closed CUDA device census for the small paired gate.
#include <cuda_runtime.h>

#include <cstdio>
#include <cstdlib>
#include <cstring>

static void check(cudaError_t status, const char* operation) {
    if (status != cudaSuccess) {
        std::fprintf(stderr, "%s: %s\n", operation, cudaGetErrorString(status));
        std::exit(2);
    }
}

int main(int argc, char** argv) {
    if (argc != 2) {
        std::fprintf(stderr, "usage: cuda_device_probe <output-dir>\n");
        return 64;
    }
    const char* rank_text = std::getenv("SLURM_PROCID");
    const char* local_text = std::getenv("SLURM_LOCALID");
    if (!rank_text || !local_text) return 65;
    char* end = nullptr;
    const long rank = std::strtol(rank_text, &end, 10);
    if (*rank_text == '\0' || *end != '\0' || rank < 0) return 66;
    end = nullptr;
    const long local = std::strtol(local_text, &end, 10);
    if (*local_text == '\0' || *end != '\0' || local < 0) return 67;

    int count = 0;
    check(cudaGetDeviceCount(&count), "cudaGetDeviceCount");
    if (count != 1) {
        std::fprintf(stderr, "expected exactly one visible CUDA device, got %d\n", count);
        return 3;
    }
    check(cudaSetDevice(0), "cudaSetDevice");
    check(cudaFree(nullptr), "cudaFree(0)");
    cudaDeviceProp properties{};
    check(cudaGetDeviceProperties(&properties, 0), "cudaGetDeviceProperties");
    size_t free_bytes = 0;
    size_t total_bytes = 0;
    check(cudaMemGetInfo(&free_bytes, &total_bytes), "cudaMemGetInfo");

    char uuid_text[33];
    for (int index = 0; index < 16; ++index)
        std::snprintf(uuid_text + 2 * index, 3, "%02x",
                      static_cast<unsigned char>(properties.uuid.bytes[index]));
    char path[4096];
    const int needed = std::snprintf(
        path, sizeof(path), "%s/rank_%05ld.txt", argv[1], rank);
    if (needed < 0 || needed >= static_cast<int>(sizeof(path))) return 68;
    FILE* stream = std::fopen(path, "w");
    if (!stream) return 69;
    std::fprintf(stream, "%ld %ld %d %s %zu %zu %s\n", rank, local, count,
                 uuid_text, free_bytes, total_bytes, properties.name);
    if (std::fclose(stream) != 0) return 70;
    return 0;
}
