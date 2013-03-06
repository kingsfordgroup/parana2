#ifndef __PARANACOMMON_H__
#define __PARANACOMMON_H__

#include <boost/pool/pool_alloc.hpp>

template<typename T>
using CustomAllocator = std::allocator<T>;//boost::pool_allocator<T>;

template<typename T>
using StackAllocator = std::allocator<T>;//boost::pool_allocator<T>;

#endif // __PARANACOMMON_H__