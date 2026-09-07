#ifndef HAMPR_IO_BUFFER_POOL_HPP
#define HAMPR_IO_BUFFER_POOL_HPP

#include <hampr/core/types.hpp>
#include <vector>
#include <mutex>
#include <condition_variable>

namespace hampr {

// Thread-safe circular buffer for continuous IQ data streaming.
// Used by StreamingDataSource to buffer data between producer and consumer.
class BufferPool {
public:
    explicit BufferPool(size_t capacity)
        : capacity_(capacity), head_(0), tail_(0), full_(false) {
        buffer_.resize(capacity);
    }

    // Producer: push data into the buffer. Returns false if buffer is full.
    bool push(const complex* data, size_t count) {
        std::lock_guard<std::mutex> lock(mtx_);
        if (full_ && available_locked() == 0)
            return false;
        size_t free = capacity_ - available_locked();
        if (count > free)
            return false;
        for (size_t i = 0; i < count; ++i) {
            buffer_[tail_] = data[i];
            tail_ = (tail_ + 1) % capacity_;
        }
        if (tail_ == head_)
            full_ = true;
        cv_.notify_one();
        return true;
    }

    // Consumer: pop up to max_count samples. Returns actually popped count.
    size_t pop(complex* data, size_t max_count) {
        std::unique_lock<std::mutex> lock(mtx_);
        cv_.wait(lock, [this] { return available_locked() > 0; });
        size_t avail = available_locked();
        size_t count = std::min(max_count, avail);
        for (size_t i = 0; i < count; ++i) {
            data[i] = buffer_[head_];
            head_ = (head_ + 1) % capacity_;
        }
        full_ = false;
        return count;
    }

    // Non-blocking pop: returns 0 if no data available
    size_t try_pop(complex* data, size_t max_count) {
        std::lock_guard<std::mutex> lock(mtx_);
        if (available_locked() == 0)
            return 0;
        size_t avail = available_locked();
        size_t count = std::min(max_count, avail);
        for (size_t i = 0; i < count; ++i) {
            data[i] = buffer_[head_];
            head_ = (head_ + 1) % capacity_;
        }
        full_ = false;
        return count;
    }

    size_t available() const {
        std::lock_guard<std::mutex> lock(mtx_);
        return available_locked();
    }

    size_t capacity() const { return capacity_; }
    bool empty() const { return available() == 0; }
    bool full() const {
        std::lock_guard<std::mutex> lock(mtx_);
        return full_;
    }

private:
    size_t available_locked() const {
        if (full_) return capacity_;
        if (tail_ >= head_) return tail_ - head_;
        return capacity_ - head_ + tail_;
    }

    std::vector<complex> buffer_;
    size_t capacity_;
    size_t head_;
    size_t tail_;
    bool full_;
    mutable std::mutex mtx_;
    mutable std::condition_variable cv_;
};

} // namespace hampr

#endif // HAMPR_IO_BUFFER_POOL_HPP
