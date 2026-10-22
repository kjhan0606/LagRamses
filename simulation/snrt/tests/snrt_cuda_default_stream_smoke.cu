#include "cuda_stream_pool.h"

// The isolated kernel smoke has no RAMSES/cuRamses stream-pool lifecycle.
// Model the uninitialized production pool contract: slot 0 resolves to the
// CUDA default stream. The production bundle build separately links and
// verifies the real cuda_stream_pool implementation.
cudaStream_t cuda_get_stream_internal(int slot) {
  (void)slot;
  return nullptr;
}
