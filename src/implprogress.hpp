#ifndef IMPLPROGRESS_HPP
#define IMPLPROGRESS_HPP
#include <memory>
#if defined(Rcpp_hpp)
// [[Rcpp::depends(RcppProgress)]]
#include <progress.hpp>
#include <progress_bar.hpp>
#else
#include <iostream>
#endif

#ifndef BAR_WIDTH
#define BAR_WIDTH 60
#endif
class ImplProgress
{
private:
#if defined(Rcpp_hpp)
	std::unique_ptr<Progress> p;
#else
	unsigned long total, current, last;
	bool aborted;
#endif
public:
#if defined(Rcpp_hpp)
	ImplProgress(unsigned long count) : p(std::make_unique<Progress>(count, true))	{}
	~ImplProgress() = default;
#else
	ImplProgress(unsigned long count) : total(count), current(0), last(0), aborted(false) {}
	~ImplProgress() { std::cout << std::endl; }
#endif
	bool check_abort()
	{
#if defined(Rcpp_hpp)
		return Progress::check_abort();
#elif defined(PYBIND11_MODULE)
		if (PyErr_CheckSignals() != 0) {
			pybind11::gil_scoped_acquire acquire;
			throw pybind11::error_already_set();
			aborted = true;
		}
		return aborted;
#else
		return aborted;
#endif
	}
	void display()
	{
#if !defined(Rcpp_hpp)
		double last_percent = (double)last / total;
		double cur_percent = (double)current / total;
#if defined(_OPENMP)
#pragma omp critical
#endif
		if (cur_percent >= 0.01 + last_percent) {
			int pct = 100 * cur_percent;
			int pos = cur_percent * BAR_WIDTH;
			std::string out(BAR_WIDTH, ' ');
			for (int i = 0; i < BAR_WIDTH; i++) {
				if (i < pos) { out[i] = '='; }
				else if (i == pos) { out[i] = '>'; }
			}
			std::cerr << "[" << out << "] " << pct << "%\r";
			std::cerr.flush();
			last = current;
		}
#endif	
	}
	void increment(unsigned long count=1) {
#if defined(Rcpp_hpp)
		p->increment(count);
#else
  #if defined(_OPENMP)
  #pragma omp atomic		
  #endif
		current += count;
  #if defined(_OPENMP)
		if (omp_get_thread_num() == 0) {
			display();
		}
  #else
		display();
  #endif
#endif		
	}
};
#endif
