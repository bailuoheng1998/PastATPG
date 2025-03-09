// SPDX-FileCopyrightText: 2023 CASTest Corporation Limited
// SPDX-License-Identifier: LGPL-v3

#ifndef ICTEST_EVENTS_H
#define ICTEST_EVENTS_H
#include <cstdint>
#include <vector>
#include <stack>
#include <queue>
#include "util/Log.h"
namespace ictest {
class EventsQueue {
 public:
  EventsQueue() = default;
  EventsQueue(int num_levels, int gate_nums){
    AllocMem(num_levels,gate_nums);
  }

  ~EventsQueue() = default;

  void InitCurrLevel(int curr_level) { curr_level_ = curr_level; }

  // obtain active gate from event queue
  uint32_t pop() noexcept {
    uint32_t gid = -1;
    while (size_) {
      if (!events_[curr_level_].empty()) {
        gid = events_[curr_level_].front();
        events_[curr_level_].pop();
        size_--;
        in_queue_[gid] = 0;
        break;
      } else {
        ++curr_level_;
        curr_level_ = (curr_level_ != total_level_) ? curr_level_ : 0;
      }
    }
    return gid;
  }

  // insert active gate into event queue
  void push(uint32_t gid, int dpi_level) noexcept {
    if (!in_queue_[gid]) {
      in_queue_[gid] = 1;
      events_[dpi_level].push(gid);
      size_++;
    }
  }

  // clear event queue and set _curr_level = 0
  void ClearEvents() noexcept {
    while (size_) {
      if (!events_[curr_level_].empty()) {
        uint32_t gid = events_[curr_level_].front();
        events_[curr_level_].pop();
        size_--;
        in_queue_[gid] = 0;
      } else {
        ++curr_level_;
        curr_level_ = (curr_level_ != total_level_) ? curr_level_ : 0;
      }
    }
  }

  [[nodiscard]] bool empty() const noexcept { return size_ == 0; }

  [[nodiscard]] uint32_t size() const noexcept { return size_; }

  [[nodiscard]] uint32_t GetCurrLevel() const noexcept { return curr_level_; }

  [[nodiscard]] uint32_t LevelSize(uint32_t dpi) const noexcept { return events_[dpi].size(); }

  void AllocMem(uint32_t level_nums,uint32_t gate_nums){
    events_.resize(level_nums);
    total_level_=level_nums;
    curr_level_=0;
    in_queue_.resize(gate_nums);
  }
 private:
  std::vector<std::queue<uint32_t>> events_;
  std::vector<uint8_t> in_queue_;
  uint32_t size_ = 0;
  uint32_t total_level_ = 0;
  uint32_t curr_level_ = 0;
};
}  // namespace ictest
#endif  // ICTEST_EVENTS_H
