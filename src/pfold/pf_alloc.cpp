#include "pf_alloc.hpp"

#include <cstdlib>
#include <new>
#include <unordered_set>

namespace pfold {

namespace {
struct Arena {
  std::unordered_set<void *> live;
  ~Arena() {
    for (void *p : live) std::free(p);
  }
};
thread_local Arena *current_arena = nullptr;

void *check(void *p, std::size_t size) {
  if (!p && size) throw std::bad_alloc();
  return p;
}
}  // namespace

void *pf_malloc(std::size_t size) {
  void *p = check(std::malloc(size ? size : 1), size);
  if (current_arena) current_arena->live.insert(p);
  return p;
}

void *pf_calloc(std::size_t n, std::size_t size) {
  void *p = check(std::calloc(n ? n : 1, size ? size : 1), n * size);
  if (current_arena) current_arena->live.insert(p);
  return p;
}

void *pf_realloc(void *ptr, std::size_t size) {
  if (!ptr) return pf_malloc(size);
  if (current_arena) current_arena->live.erase(ptr);
  void *p = check(std::realloc(ptr, size ? size : 1), size);
  if (current_arena) current_arena->live.insert(p);
  return p;
}

void pf_free(void *ptr) {
  if (!ptr) return;
  if (current_arena) {
    auto it = current_arena->live.find(ptr);
    if (it == current_arena->live.end()) return;  // not ours (or already released): ignore
    current_arena->live.erase(it);
  }
  std::free(ptr);
}

void pf_fail(const char *msg) { throw EngineError(msg); }

// defined in the engine translation units that own file-level state
void reset_file_state();
void reset_newcolprob_state();
void reset_optimize_state();
void reset_phyl_state();

void reset_engine_state() {
  reset_file_state();
  reset_newcolprob_state();
  reset_optimize_state();
  reset_phyl_state();
}

ArenaScope::ArenaScope() : previous_(current_arena) {
  reset_engine_state();
  current_arena = new Arena();
}

ArenaScope::~ArenaScope() {
  delete current_arena;
  current_arena = static_cast<Arena *>(previous_);
  reset_engine_state();
}

}  // namespace pfold
