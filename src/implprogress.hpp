#ifndef IMPLPROGRESS_HPP
#define IMPLPROGRESS_HPP
#include <memory>
#ifdef Rcpp_hpp
// [[Rcpp::depends(RcppProgress)]]
#include <progress.hpp>
#include <progress_bar.hpp>
#endif
class ImplProgress
{
private:
#ifdef Rcpp_hpp
	std::unique_ptr<Progress> p;
#else
	unsigned long total;
#endif
public:
#ifdef Rcpp_hpp
	ImplProgress(unsigned long count) : p(std::make_unique<Progress>(count, true)) {}
#else
	ImplProgress(unsigned long count) : total(count) {}
#endif
	~ImplProgress() = default;
	bool check_abort() const
	{
#ifdef Rcpp_hpp
		return Progress::check_abort();
#else
		return false;
#endif
	}
	void increment(unsigned long count=1) {
#ifdef Rcpp_hpp
		p->increment(count);
#else
#if defined(_OPENMP)
#pragma omp atomic
#endif
		total += count;
#endif
	}
};
#endif
