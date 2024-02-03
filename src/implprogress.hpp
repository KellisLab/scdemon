#ifndef IMPLPROGRESS_HPP
#define IMPLPROGRESS_HPP
#include <memory>
#if defined(Rcpp_hpp)
// [[Rcpp::depends(RcppProgress)]]
#include <progress.hpp>
#include <progress_bar.hpp>
typedef void* implprogress_callback;
#elif defined(PYBIND11_MODULE)
#include <pybind11/pybind11.h>
#include <pybind11/functional.h>
typedef pybind11::function implprogress_callback;
#else
typedef void* implprogress_callback;
#endif

class ImplProgress
{
private:
#if defined(Rcpp_hpp)
	std::unique_ptr<Progress> p;
#else
	long total, current, last;
 	bool aborted;
	const implprogress_callback callback;
	const implprogress_callback interrupt_checker;
#endif
public:
#if defined(Rcpp_hpp)
	ImplProgress(long count, implprogress_callback _empty, implprogress_callback _empty2) : p(std::make_unique<Progress>(static_cast<unsigned long>(count), true))  {}
#elif defined(PYBIND11_MODULE)
	ImplProgress(long count, implprogress_callback mcallback, implprogress_callback interrupt_check) : total(count), callback(mcallback), interrupt_checker(interrupt_check), current(0), last(0), aborted(false) {
	  pybind11::gil_scoped_acquire acquire;
	  callback(count);
	}
#else
	ImplProgress(long count, implprogress_callback _empty, implprogress_callback _empty2) : total(count), callback(_empty), interrupt_checker(_empty2), current(0), last(0), aborted(false) {}
#endif
#if defined(PYBIND11_MODULE)
  ~ImplProgress() {
    pybind11::gil_scoped_acquire acquire;
    callback(total);
    callback(-1);
  }
#else
	~ImplProgress() = default;
#endif
	bool check_abort() const
	{
#if defined(Rcpp_hpp)
		return Progress::check_abort();
#else
		return aborted;
#endif
	}
	void increment(long count=1) {
#if defined(Rcpp_hpp)
		p->increment(static_cast<unsigned long>(count));
#elif defined(PYBIND11_MODULE)
  #if defined(_OPENMP)
  #pragma omp atomic		
  #endif
		current += count;
		try {
  #if defined(_OPENMP)
			if (omp_get_thread_num() == 0) {
				if (100 * (current - last) >= total) {
					pybind11::gil_scoped_acquire acquire;
					callback(current - last);
					last = current;
				}
			}
  #else
			if (100 * (current - last) >= total) {
				pybind11::gil_scoped_acquire acquire;
				callback(current - last);
				last = current;
			}
  #endif
		} catch (pybind11::error_already_set& e) {
			if (e.matches(PyExc_KeyboardInterrupt)) {
				aborted = true;
			}
			throw;  // Rethrow if it's not a KeyboardInterrupt
		}
#else
  #if defined(_OPENMP)
  #pragma omp atomic
  #endif
		total -= count;
#endif
	}
};
#endif
