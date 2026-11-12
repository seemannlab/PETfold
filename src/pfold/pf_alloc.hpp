#pragma once
/*
 * pf_alloc.hpp - scoped memory management for the legacy Pfold engine.
 *
 * The Pfold engine (Knudsen) was written in C with manual memory management
 * and frees only part of what it allocates.  Instead of rewriting every
 * allocation site, all engine allocations go through pf_malloc/pf_calloc/
 * pf_realloc/pf_free.  While an ArenaScope is active on the current thread,
 * every live allocation is tracked and released when the scope ends (also
 * when an exception propagates).  Without an active scope the functions
 * behave exactly like their libc counterparts (used by drawdot/drawphyl).
 *
 * Engine errors that used to call exit(1) now call pf_fail(), which throws
 * pfold::EngineError.
 */
#include <cstddef>
#include <stdexcept>
#include <string>

namespace pfold {

class EngineError : public std::runtime_error {
 public:
  using std::runtime_error::runtime_error;
};

void *pf_malloc(std::size_t size);
void *pf_calloc(std::size_t n, std::size_t size);
void *pf_realloc(void *ptr, std::size_t size);
void pf_free(void *ptr);

[[noreturn]] void pf_fail(const char *msg = "Pfold engine error (see messages above)");

/* Resets the file-level state of the engine (called on scope entry/exit). */
void reset_engine_state();

class ArenaScope {
 public:
  ArenaScope();
  ~ArenaScope();
  ArenaScope(const ArenaScope &) = delete;
  ArenaScope &operator=(const ArenaScope &) = delete;
 private:
  void *previous_;
};

}  // namespace pfold
