#define __STDC_LIMIT_MACROS
#include <stdint.h>

/***************************************************************************************[Solver.cc]
 Glucose -- Copyright (c) 2013, Gilles Audemard, Laurent Simon
				CRIL - Univ. Artois, France
				LRI  - Univ. Paris Sud, France
 
Glucose sources are based on MiniSat (see below MiniSat copyrights). Permissions and copyrights of
Glucose are exactly the same as Minisat on which it is based on. (see below).

---------------

Copyright (c) 2003-2006, Niklas Een, Niklas Sorensson
Copyright (c) 2007-2010, Niklas Sorensson

Permission is hereby granted, free of charge, to any person obtaining a copy of this software and
associated documentation files (the "Software"), to deal in the Software without restriction,
including without limitation the rights to use, copy, modify, merge, publish, distribute,
sublicense, and/or sell copies of the Software, and to permit persons to whom the Software is
furnished to do so, subject to the following conditions:

The above copyright notice and this permission notice shall be included in all copies or
substantial portions of the Software.

THE SOFTWARE IS PROVIDED "AS IS", WITHOUT WARRANTY OF ANY KIND, EXPRESS OR IMPLIED, INCLUDING BUT
NOT LIMITED TO THE WARRANTIES OF MERCHANTABILITY, FITNESS FOR A PARTICULAR PURPOSE AND
NONINFRINGEMENT. IN NO EVENT SHALL THE AUTHORS OR COPYRIGHT HOLDERS BE LIABLE FOR ANY CLAIM,
DAMAGES OR OTHER LIABILITY, WHETHER IN AN ACTION OF CONTRACT, TORT OR OTHERWISE, ARISING FROM, OUT
OF OR IN CONNECTION WITH THE SOFTWARE OR THE USE OR OTHER DEALINGS IN THE SOFTWARE.
**************************************************************************************************/

#include <math.h>

/******************************************************************************************[Sort.h]
Copyright (c) 2003-2007, Niklas Een, Niklas Sorensson
Copyright (c) 2007-2010, Niklas Sorensson

Permission is hereby granted, free of charge, to any person obtaining a copy of this software and
associated documentation files (the "Software"), to deal in the Software without restriction,
including without limitation the rights to use, copy, modify, merge, publish, distribute,
sublicense, and/or sell copies of the Software, and to permit persons to whom the Software is
furnished to do so, subject to the following conditions:

The above copyright notice and this permission notice shall be included in all copies or
substantial portions of the Software.

THE SOFTWARE IS PROVIDED "AS IS", WITHOUT WARRANTY OF ANY KIND, EXPRESS OR IMPLIED, INCLUDING BUT
NOT LIMITED TO THE WARRANTIES OF MERCHANTABILITY, FITNESS FOR A PARTICULAR PURPOSE AND
NONINFRINGEMENT. IN NO EVENT SHALL THE AUTHORS OR COPYRIGHT HOLDERS BE LIABLE FOR ANY CLAIM,
DAMAGES OR OTHER LIABILITY, WHETHER IN AN ACTION OF CONTRACT, TORT OR OTHERWISE, ARISING FROM, OUT
OF OR IN CONNECTION WITH THE SOFTWARE OR THE USE OR OTHER DEALINGS IN THE SOFTWARE.
**************************************************************************************************/

#ifndef Glucose_Sort_h
#define Glucose_Sort_h

/*******************************************************************************************[Vec.h]
Copyright (c) 2003-2007, Niklas Een, Niklas Sorensson
Copyright (c) 2007-2010, Niklas Sorensson

Permission is hereby granted, free of charge, to any person obtaining a copy of this software and
associated documentation files (the "Software"), to deal in the Software without restriction,
including without limitation the rights to use, copy, modify, merge, publish, distribute,
sublicense, and/or sell copies of the Software, and to permit persons to whom the Software is
furnished to do so, subject to the following conditions:

The above copyright notice and this permission notice shall be included in all copies or
substantial portions of the Software.

THE SOFTWARE IS PROVIDED "AS IS", WITHOUT WARRANTY OF ANY KIND, EXPRESS OR IMPLIED, INCLUDING BUT
NOT LIMITED TO THE WARRANTIES OF MERCHANTABILITY, FITNESS FOR A PARTICULAR PURPOSE AND
NONINFRINGEMENT. IN NO EVENT SHALL THE AUTHORS OR COPYRIGHT HOLDERS BE LIABLE FOR ANY CLAIM,
DAMAGES OR OTHER LIABILITY, WHETHER IN AN ACTION OF CONTRACT, TORT OR OTHERWISE, ARISING FROM, OUT
OF OR IN CONNECTION WITH THE SOFTWARE OR THE USE OR OTHER DEALINGS IN THE SOFTWARE.
**************************************************************************************************/

#ifndef Glucose_Vec_h
#define Glucose_Vec_h

#include <assert.h>
#include <new>

/**************************************************************************************[IntTypes.h]
Copyright (c) 2009-2010, Niklas Sorensson

Permission is hereby granted, free of charge, to any person obtaining a copy of this software and
associated documentation files (the "Software"), to deal in the Software without restriction,
including without limitation the rights to use, copy, modify, merge, publish, distribute,
sublicense, and/or sell copies of the Software, and to permit persons to whom the Software is
furnished to do so, subject to the following conditions:

The above copyright notice and this permission notice shall be included in all copies or
substantial portions of the Software.

THE SOFTWARE IS PROVIDED "AS IS", WITHOUT WARRANTY OF ANY KIND, EXPRESS OR IMPLIED, INCLUDING BUT
NOT LIMITED TO THE WARRANTIES OF MERCHANTABILITY, FITNESS FOR A PARTICULAR PURPOSE AND
NONINFRINGEMENT. IN NO EVENT SHALL THE AUTHORS OR COPYRIGHT HOLDERS BE LIABLE FOR ANY CLAIM,
DAMAGES OR OTHER LIABILITY, WHETHER IN AN ACTION OF CONTRACT, TORT OR OTHERWISE, ARISING FROM, OUT
OF OR IN CONNECTION WITH THE SOFTWARE OR THE USE OR OTHER DEALINGS IN THE SOFTWARE.
**************************************************************************************************/

#ifndef Glucose_IntTypes_h
#define Glucose_IntTypes_h

#ifdef __sun
    // Not sure if there are newer versions that support C99 headers. The
    // needed features are implemented in the headers below though:

#   include <sys/int_types.h>
#   include <sys/int_fmtio.h>
#   include <sys/int_limits.h>

#else

#   include <stdint.h>
#   include <inttypes.h>

#endif

#include <limits.h>

#ifndef PRIu64
#define PRIu64 "lu"
#define PRIi64 "ld"
#endif
//=================================================================================================

#endif

/****************************************************************************************[XAlloc.h]
Copyright (c) 2009-2010, Niklas Sorensson

Permission is hereby granted, free of charge, to any person obtaining a copy of this software and
associated documentation files (the "Software"), to deal in the Software without restriction,
including without limitation the rights to use, copy, modify, merge, publish, distribute,
sublicense, and/or sell copies of the Software, and to permit persons to whom the Software is
furnished to do so, subject to the following conditions:

The above copyright notice and this permission notice shall be included in all copies or
substantial portions of the Software.

THE SOFTWARE IS PROVIDED "AS IS", WITHOUT WARRANTY OF ANY KIND, EXPRESS OR IMPLIED, INCLUDING BUT
NOT LIMITED TO THE WARRANTIES OF MERCHANTABILITY, FITNESS FOR A PARTICULAR PURPOSE AND
NONINFRINGEMENT. IN NO EVENT SHALL THE AUTHORS OR COPYRIGHT HOLDERS BE LIABLE FOR ANY CLAIM,
DAMAGES OR OTHER LIABILITY, WHETHER IN AN ACTION OF CONTRACT, TORT OR OTHERWISE, ARISING FROM, OUT
OF OR IN CONNECTION WITH THE SOFTWARE OR THE USE OR OTHER DEALINGS IN THE SOFTWARE.
**************************************************************************************************/


#ifndef Glucose_XAlloc_h
#define Glucose_XAlloc_h

#include <errno.h>
#include <stdlib.h>
#include <stdio.h>

namespace Glucose {

//=================================================================================================
// Simple layer on top of malloc/realloc to catch out-of-memory situtaions and provide some typing:

class OutOfMemoryException{};
static inline void* xrealloc(void *ptr, size_t size)
{
    void* mem = realloc(ptr, size);
    if (mem == NULL && errno == ENOMEM){
        throw OutOfMemoryException();
    }else {
        return mem;
	}
}

//=================================================================================================
}

#endif


namespace Glucose {

//=================================================================================================
// Automatically resizable arrays
//
// NOTE! Don't use this vector on datatypes that cannot be re-located in memory (with realloc)

template<class T>
class vec {
    T*  data;
    int sz;
    int cap;

    // Don't allow copying (error prone):
    vec<T>&  operator = (vec<T>& other) { assert(0); return *this; }
             vec        (vec<T>& other) { assert(0); }
             
    // Helpers for calculating next capacity:
    static inline int  imax   (int x, int y) { int mask = (y-x) >> (sizeof(int)*8-1); return (x&mask) + (y&(~mask)); }
    //static inline void nextCap(int& cap){ cap += ((cap >> 1) + 2) & ~1; }
    static inline void nextCap(int& cap){ cap += ((cap >> 1) + 2) & ~1; }

public:
    // Constructors:
    vec()                       : data(NULL) , sz(0)   , cap(0)    { }
    explicit vec(int size)      : data(NULL) , sz(0)   , cap(0)    { growTo(size); }
    vec(int size, const T& pad) : data(NULL) , sz(0)   , cap(0)    { growTo(size, pad); }
   ~vec()                                                          { clear(true); }

    // Pointer to first element:
    operator T*       (void)           { return data; }

    // Size operations:
    int      size     (void) const     { return sz; }
    void     shrink   (int nelems)     { assert(nelems <= sz); for (int i = 0; i < nelems; i++) sz--, data[sz].~T(); }
    void     shrink_  (int nelems)     { assert(nelems <= sz); sz -= nelems; }
    int      capacity (void) const     { return cap; }
    void     capacity (int min_cap);
    void     growTo   (int size);
    void     growTo   (int size, const T& pad);
    void     clear    (bool dealloc = false);

    // Stack interface:
    void     push  (void)              { if (sz == cap) capacity(sz+1); new (&data[sz]) T(); sz++; }
    void     push  (const T& elem)     { if (sz == cap) capacity(sz+1); data[sz++] = elem; }
    void     push_ (const T& elem)     { assert(sz < cap); data[sz++] = elem; }
    void     pop   (void)              { assert(sz > 0); sz--, data[sz].~T(); }
    // NOTE: it seems possible that overflow can happen in the 'sz+1' expression of 'push()', but
    // in fact it can not since it requires that 'cap' is equal to INT_MAX. This in turn can not
    // happen given the way capacities are calculated (below). Essentially, all capacities are
    // even, but INT_MAX is odd.

    const T& last  (void) const        { return data[sz-1]; }
    T&       last  (void)              { return data[sz-1]; }

    // Vector interface:
    const T& operator [] (int index) const { return data[index]; }
    T&       operator [] (int index)       { return data[index]; }

    // Duplicatation (preferred instead):
    void copyTo(vec<T>& copy) const { copy.clear(); copy.growTo(sz); for (int i = 0; i < sz; i++) copy[i] = data[i]; }
    void moveTo(vec<T>& dest) { dest.clear(true); dest.data = data; dest.sz = sz; dest.cap = cap; data = NULL; sz = 0; cap = 0; }
};


template<class T>
void vec<T>::capacity(int min_cap) {
    if (cap >= min_cap) return;
    int add = imax((min_cap - cap + 1) & ~1, ((cap >> 1) + 2) & ~1);   // NOTE: grow by approximately 3/2
    if (add > INT_MAX - cap || ((data = (T*)::realloc(data, (cap += add) * sizeof(T))) == NULL) && errno == ENOMEM)
        throw OutOfMemoryException();
 }


template<class T>
void vec<T>::growTo(int size, const T& pad) {
    if (sz >= size) return;
    capacity(size);
    for (int i = sz; i < size; i++) data[i] = pad;
    sz = size; }


template<class T>
void vec<T>::growTo(int size) {
    if (sz >= size) return;
    capacity(size);
    for (int i = sz; i < size; i++) new (&data[i]) T();
    sz = size; }


template<class T>
void vec<T>::clear(bool dealloc) {
    if (data != NULL){
        for (int i = 0; i < sz; i++) data[i].~T();
        sz = 0;
        if (dealloc) free(data), data = NULL, cap = 0; } }

//=================================================================================================
}

#endif


//=================================================================================================
// Some sorting algorithms for vec's


namespace Glucose {

template<class T>
struct LessThan_default {
    bool operator () (T x, T y) { return x < y; }
};


template <class T, class LessThan>
void selectionSort(T* array, int size, LessThan lt)
{
    int     i, j, best_i;
    T       tmp;

    for (i = 0; i < size-1; i++){
        best_i = i;
        for (j = i+1; j < size; j++){
            if (lt(array[j], array[best_i]))
                best_i = j;
        }
        tmp = array[i]; array[i] = array[best_i]; array[best_i] = tmp;
    }
}
template <class T> static inline void selectionSort(T* array, int size) {
    selectionSort(array, size, LessThan_default<T>()); }

template <class T, class LessThan>
void sort(T* array, int size, LessThan lt)
{
    if (size <= 15)
        selectionSort(array, size, lt);

    else{
        T           pivot = array[size / 2];
        T           tmp;
        int         i = -1;
        int         j = size;

        for(;;){
            do i++; while(lt(array[i], pivot));
            do j--; while(lt(pivot, array[j]));

            if (i >= j) break;

            tmp = array[i]; array[i] = array[j]; array[j] = tmp;
        }

        sort(array    , i     , lt);
        sort(&array[i], size-i, lt);
    }
}
template <class T> static inline void sort(T* array, int size) {
    sort(array, size, LessThan_default<T>()); }


//=================================================================================================
// For 'vec's:


template <class T, class LessThan> void sort(vec<T>& v, LessThan lt) {
    sort((T*)v, v.size(), lt); }
template <class T> void sort(vec<T>& v) {
    sort(v, LessThan_default<T>()); }


//=================================================================================================
}

#endif

/****************************************************************************************[Solver.h]
 Glucose -- Copyright (c) 2009, Gilles Audemard, Laurent Simon
				CRIL - Univ. Artois, France
				LRI  - Univ. Paris Sud, France
 
Glucose sources are based on MiniSat (see below MiniSat copyrights). Permissions and copyrights of
Glucose are exactly the same as Minisat on which it is based on. (see below).

---------------
Copyright (c) 2003-2006, Niklas Een, Niklas Sorensson
Copyright (c) 2007-2010, Niklas Sorensson

Permission is hereby granted, free of charge, to any person obtaining a copy of this software and
associated documentation files (the "Software"), to deal in the Software without restriction,
including without limitation the rights to use, copy, modify, merge, publish, distribute,
sublicense, and/or sell copies of the Software, and to permit persons to whom the Software is
furnished to do so, subject to the following conditions:

The above copyright notice and this permission notice shall be included in all copies or
substantial portions of the Software.

THE SOFTWARE IS PROVIDED "AS IS", WITHOUT WARRANTY OF ANY KIND, EXPRESS OR IMPLIED, INCLUDING BUT
NOT LIMITED TO THE WARRANTIES OF MERCHANTABILITY, FITNESS FOR A PARTICULAR PURPOSE AND
NONINFRINGEMENT. IN NO EVENT SHALL THE AUTHORS OR COPYRIGHT HOLDERS BE LIABLE FOR ANY CLAIM,
DAMAGES OR OTHER LIABILITY, WHETHER IN AN ACTION OF CONTRACT, TORT OR OTHERWISE, ARISING FROM, OUT
OF OR IN CONNECTION WITH THE SOFTWARE OR THE USE OR OTHER DEALINGS IN THE SOFTWARE.
**************************************************************************************************/

#ifndef Glucose_Solver_h
#define Glucose_Solver_h

/******************************************************************************************[Heap.h]
Copyright (c) 2003-2006, Niklas Een, Niklas Sorensson
Copyright (c) 2007-2010, Niklas Sorensson

Permission is hereby granted, free of charge, to any person obtaining a copy of this software and
associated documentation files (the "Software"), to deal in the Software without restriction,
including without limitation the rights to use, copy, modify, merge, publish, distribute,
sublicense, and/or sell copies of the Software, and to permit persons to whom the Software is
furnished to do so, subject to the following conditions:

The above copyright notice and this permission notice shall be included in all copies or
substantial portions of the Software.

THE SOFTWARE IS PROVIDED "AS IS", WITHOUT WARRANTY OF ANY KIND, EXPRESS OR IMPLIED, INCLUDING BUT
NOT LIMITED TO THE WARRANTIES OF MERCHANTABILITY, FITNESS FOR A PARTICULAR PURPOSE AND
NONINFRINGEMENT. IN NO EVENT SHALL THE AUTHORS OR COPYRIGHT HOLDERS BE LIABLE FOR ANY CLAIM,
DAMAGES OR OTHER LIABILITY, WHETHER IN AN ACTION OF CONTRACT, TORT OR OTHERWISE, ARISING FROM, OUT
OF OR IN CONNECTION WITH THE SOFTWARE OR THE USE OR OTHER DEALINGS IN THE SOFTWARE.
**************************************************************************************************/

#ifndef Glucose_Heap_h
#define Glucose_Heap_h


namespace Glucose {

//=================================================================================================
// A heap implementation with support for decrease/increase key.


template<class Comp>
class Heap {
    Comp     lt;       // The heap is a minimum-heap with respect to this comparator
    vec<int> heap;     // Heap of integers
    vec<int> indices;  // Each integers position (index) in the Heap

    // Index "traversal" functions
    static inline int left  (int i) { return i*2+1; }
    static inline int right (int i) { return (i+1)*2; }
    static inline int parent(int i) { return (i-1) >> 1; }


    void percolateUp(int i)
    {
        int x  = heap[i];
        int p  = parent(i);
        
        while (i != 0 && lt(x, heap[p])){
            heap[i]          = heap[p];
            indices[heap[p]] = i;
            i                = p;
            p                = parent(p);
        }
        heap   [i] = x;
        indices[x] = i;
    }


    void percolateDown(int i)
    {
        int x = heap[i];
        while (left(i) < heap.size()){
            int child = right(i) < heap.size() && lt(heap[right(i)], heap[left(i)]) ? right(i) : left(i);
            if (!lt(heap[child], x)) break;
            heap[i]          = heap[child];
            indices[heap[i]] = i;
            i                = child;
        }
        heap   [i] = x;
        indices[x] = i;
    }


  public:
    Heap(const Comp& c) : lt(c) { }

    int  size      ()          const { return heap.size(); }
    bool empty     ()          const { return heap.size() == 0; }
    bool inHeap    (int n)     const { return n < indices.size() && indices[n] >= 0; }
    int  operator[](int index) const { assert(index < heap.size()); return heap[index]; }


    void decrease  (int n) { assert(inHeap(n)); percolateUp  (indices[n]); }
    void increase  (int n) { assert(inHeap(n)); percolateDown(indices[n]); }


    // Safe variant of insert/decrease/increase:
    void update(int n)
    {
        if (!inHeap(n))
            insert(n);
        else {
            percolateUp(indices[n]);
            percolateDown(indices[n]); }
    }


    void insert(int n)
    {
        indices.growTo(n+1, -1);
        assert(!inHeap(n));

        indices[n] = heap.size();
        heap.push(n);
        percolateUp(indices[n]); 
    }


    int  removeMin()
    {
        int x            = heap[0];
        heap[0]          = heap.last();
        indices[heap[0]] = 0;
        indices[x]       = -1;
        heap.pop();
        if (heap.size() > 1) percolateDown(0);
        return x; 
    }


    // Rebuild the heap from scratch, using the elements in 'ns':
    void build(vec<int>& ns) {
        for (int i = 0; i < heap.size(); i++)
            indices[heap[i]] = -1;
        heap.clear();

        for (int i = 0; i < ns.size(); i++){
            indices[ns[i]] = i;
            heap.push(ns[i]); }

        for (int i = heap.size() / 2 - 1; i >= 0; i--)
            percolateDown(i);
    }

    void clear(bool dealloc = false) 
    { 
        for (int i = 0; i < heap.size(); i++)
            indices[heap[i]] = -1;
        heap.clear(dealloc); 
    }
};


//=================================================================================================
}

#endif

/*******************************************************************************************[Alg.h]
Copyright (c) 2003-2006, Niklas Een, Niklas Sorensson
Copyright (c) 2007-2010, Niklas Sorensson

Permission is hereby granted, free of charge, to any person obtaining a copy of this software and
associated documentation files (the "Software"), to deal in the Software without restriction,
including without limitation the rights to use, copy, modify, merge, publish, distribute,
sublicense, and/or sell copies of the Software, and to permit persons to whom the Software is
furnished to do so, subject to the following conditions:

The above copyright notice and this permission notice shall be included in all copies or
substantial portions of the Software.

THE SOFTWARE IS PROVIDED "AS IS", WITHOUT WARRANTY OF ANY KIND, EXPRESS OR IMPLIED, INCLUDING BUT
NOT LIMITED TO THE WARRANTIES OF MERCHANTABILITY, FITNESS FOR A PARTICULAR PURPOSE AND
NONINFRINGEMENT. IN NO EVENT SHALL THE AUTHORS OR COPYRIGHT HOLDERS BE LIABLE FOR ANY CLAIM,
DAMAGES OR OTHER LIABILITY, WHETHER IN AN ACTION OF CONTRACT, TORT OR OTHERWISE, ARISING FROM, OUT
OF OR IN CONNECTION WITH THE SOFTWARE OR THE USE OR OTHER DEALINGS IN THE SOFTWARE.
**************************************************************************************************/

#ifndef Glucose_Alg_h
#define Glucose_Alg_h


namespace Glucose {

//=================================================================================================
// Useful functions on vector-like types:

//=================================================================================================
// Removing and searching for elements:
//

template<class V, class T>
static inline void remove(V& ts, const T& t)
{
    int j = 0;
    for (; j < ts.size() && ts[j] != t; j++);
    assert(j < ts.size());
    for (; j < ts.size()-1; j++) ts[j] = ts[j+1];
    ts.pop();
}


template<class V, class T>
static inline bool find(V& ts, const T& t)
{
    int j = 0;
    for (; j < ts.size() && ts[j] != t; j++);
    return j < ts.size();
}


//=================================================================================================
// Copying vectors with support for nested vector types:
//

// Base case:
template<class T>
static inline void copy(const T& from, T& to)
{
    to = from;
}

// Recursive case:
template<class T>
static inline void copy(const vec<T>& from, vec<T>& to, bool append = false)
{
    if (!append)
        to.clear();
    for (int i = 0; i < from.size(); i++){
        to.push();
        copy(from[i], to.last());
    }
}

template<class T>
static inline void append(const vec<T>& from, vec<T>& to){ copy(from, to, true); }

//=================================================================================================
}

#endif

/***************************************************************************************[Options.h]
Copyright (c) 2008-2010, Niklas Sorensson

Permission is hereby granted, free of charge, to any person obtaining a copy of this software and
associated documentation files (the "Software"), to deal in the Software without restriction,
including without limitation the rights to use, copy, modify, merge, publish, distribute,
sublicense, and/or sell copies of the Software, and to permit persons to whom the Software is
furnished to do so, subject to the following conditions:

The above copyright notice and this permission notice shall be included in all copies or
substantial portions of the Software.

THE SOFTWARE IS PROVIDED "AS IS", WITHOUT WARRANTY OF ANY KIND, EXPRESS OR IMPLIED, INCLUDING BUT
NOT LIMITED TO THE WARRANTIES OF MERCHANTABILITY, FITNESS FOR A PARTICULAR PURPOSE AND
NONINFRINGEMENT. IN NO EVENT SHALL THE AUTHORS OR COPYRIGHT HOLDERS BE LIABLE FOR ANY CLAIM,
DAMAGES OR OTHER LIABILITY, WHETHER IN AN ACTION OF CONTRACT, TORT OR OTHERWISE, ARISING FROM, OUT
OF OR IN CONNECTION WITH THE SOFTWARE OR THE USE OR OTHER DEALINGS IN THE SOFTWARE.
**************************************************************************************************/

#ifndef Glucose_Options_h
#define Glucose_Options_h

#include <stdlib.h>
#include <stdio.h>
#include <math.h>
#include <string.h>

/************************************************************************************[ParseUtils.h]
Copyright (c) 2003-2006, Niklas Een, Niklas Sorensson
Copyright (c) 2007-2010, Niklas Sorensson

Permission is hereby granted, free of charge, to any person obtaining a copy of this software and
associated documentation files (the "Software"), to deal in the Software without restriction,
including without limitation the rights to use, copy, modify, merge, publish, distribute,
sublicense, and/or sell copies of the Software, and to permit persons to whom the Software is
furnished to do so, subject to the following conditions:

The above copyright notice and this permission notice shall be included in all copies or
substantial portions of the Software.

THE SOFTWARE IS PROVIDED "AS IS", WITHOUT WARRANTY OF ANY KIND, EXPRESS OR IMPLIED, INCLUDING BUT
NOT LIMITED TO THE WARRANTIES OF MERCHANTABILITY, FITNESS FOR A PARTICULAR PURPOSE AND
NONINFRINGEMENT. IN NO EVENT SHALL THE AUTHORS OR COPYRIGHT HOLDERS BE LIABLE FOR ANY CLAIM,
DAMAGES OR OTHER LIABILITY, WHETHER IN AN ACTION OF CONTRACT, TORT OR OTHERWISE, ARISING FROM, OUT
OF OR IN CONNECTION WITH THE SOFTWARE OR THE USE OR OTHER DEALINGS IN THE SOFTWARE.
**************************************************************************************************/

#ifndef Glucose_ParseUtils_h
#define Glucose_ParseUtils_h

#include <stdlib.h>
#include <stdio.h>
#include <math.h>

#include <zlib.h>

namespace Glucose {

//-------------------------------------------------------------------------------------------------
// A simple buffered character stream class:

static const int buffer_size = 1048576;


class StreamBuffer {
    gzFile        in;
    unsigned char buf[buffer_size];
    int           pos;
    int           size;

    void assureLookahead() {
        if (pos >= size) {
            pos  = 0;
            size = gzread(in, buf, sizeof(buf)); } }

public:
    explicit StreamBuffer(gzFile i) : in(i), pos(0), size(0) { assureLookahead(); }

    int  operator *  () const { return (pos >= size) ? EOF : buf[pos]; }
    void operator ++ ()       { pos++; assureLookahead(); }
    int  position    () const { return pos; }
};


//-------------------------------------------------------------------------------------------------
// End-of-file detection functions for StreamBuffer and char*:


static inline bool isEof(StreamBuffer& in) { return *in == EOF;  }
static inline bool isEof(const char*   in) { return *in == '\0'; }

//-------------------------------------------------------------------------------------------------
// Generic parse functions parametrized over the input-stream type.


template<class B>
static void skipWhitespace(B& in) {
    while ((*in >= 9 && *in <= 13) || *in == 32)
        ++in; }


template<class B>
static void skipLine(B& in) {
    for (;;){
        if (isEof(in)) return;
        if (*in == '\n') { ++in; return; }
        ++in; } }

template<class B>
static double parseDouble(B& in) { // only in the form X.XXXXXe-XX
    bool    neg= false;
	double accu = 0.0;
	double currentExponent = 1;
	int exponent;
	
    skipWhitespace(in);
    if(*in == EOF) return 0;
    if      (*in == '-') neg = true, ++in;
    else if (*in == '+') ++in;
    if (*in < '1' || *in > '9') printf("PARSE ERROR! Unexpected char: %c\n", *in), exit(3);
	accu = (double)(*in - '0');
	++in;
	if (*in != '.') printf("PARSE ERROR! Unexpected char: %c\n", *in),exit(3);
	++in; // skip dot
	currentExponent = 0.1;
    while (*in >= '0' && *in <= '9')
        accu = accu + currentExponent * ((double)(*in - '0')),
		currentExponent /= 10,
        ++in;
	if (*in != 'e') printf("PARSE ERROR! Unexpected char: %c\n", *in),exit(3);
	++in; // skip dot
	exponent = parseInt(in); // read exponent
	accu *= pow(10,exponent);
	return neg ? -accu:accu;
}


template<class B>
static int parseInt(B& in) {
    int     val = 0;
    bool    neg = false;
    skipWhitespace(in);
    if      (*in == '-') neg = true, ++in;
    else if (*in == '+') ++in;
    if (*in < '0' || *in > '9') fprintf(stderr, "PARSE ERROR! Unexpected char: %c\n", *in), exit(3);
    while (*in >= '0' && *in <= '9')
        val = val*10 + (*in - '0'),
        ++in;
    return neg ? -val : val; }


// String matching: in case of a match the input iterator will be advanced the corresponding
// number of characters.
template<class B>
static bool match(B& in, const char* str) {
    int i;
    for (i = 0; str[i] != '\0'; i++)
        if (in[i] != str[i])
            return false;

    in += i;

    return true; 
}

// String matching: consumes characters eagerly, but does not require random access iterator.
template<class B>
static bool eagerMatch(B& in, const char* str) {
    for (; *str != '\0'; ++str, ++in)
        if (*str != *in)
            return false;
    return true; }


//=================================================================================================
}

#endif


namespace Glucose {

//==================================================================================================
// Top-level option parse/help functions:


extern void parseOptions     (int& argc, char** argv, bool strict = false);
extern void printUsageAndExit(int  argc, char** argv, bool verbose = false);
extern void setUsageHelp     (const char* str);
extern void setHelpPrefixStr (const char* str);


//==================================================================================================
// Options is an abstract class that gives the interface for all types options:


class Option
{
 protected:
    const char* name;
    const char* description;
    const char* category;
    const char* type_name;

    static vec<Option*>& getOptionList () { static vec<Option*> options; return options; }
    static const char*&  getUsageString() { static const char* usage_str; return usage_str; }
    static const char*&  getHelpPrefixString() { static const char* help_prefix_str = ""; return help_prefix_str; }

    struct OptionLt {
        bool operator()(const Option* x, const Option* y) {
            int test1 = strcmp(x->category, y->category);
            return test1 < 0 || test1 == 0 && strcmp(x->type_name, y->type_name) < 0;
        }
    };

    Option(const char* name_, 
           const char* desc_,
           const char* cate_,
           const char* type_) : 
      name       (name_)
    , description(desc_)
    , category   (cate_)
    , type_name  (type_)
    { 
        getOptionList().push(this); 
    }

 public:
    virtual ~Option() {}

    virtual bool parse             (const char* str)      = 0;
    virtual void help              (bool verbose = false) = 0;

    friend  void parseOptions      (int& argc, char** argv, bool strict);
    friend  void printUsageAndExit (int  argc, char** argv, bool verbose);
    friend  void setUsageHelp      (const char* str);
    friend  void setHelpPrefixStr  (const char* str);
};


//==================================================================================================
// Range classes with specialization for floating types:


struct IntRange {
    int begin;
    int end;
    IntRange(int b, int e) : begin(b), end(e) {}
};

struct Int64Range {
    int64_t begin;
    int64_t end;
    Int64Range(int64_t b, int64_t e) : begin(b), end(e) {}
};

struct DoubleRange {
    double begin;
    double end;
    bool  begin_inclusive;
    bool  end_inclusive;
    DoubleRange(double b, bool binc, double e, bool einc) : begin(b), end(e), begin_inclusive(binc), end_inclusive(einc) {}
};


//==================================================================================================
// Double options:


class DoubleOption : public Option
{
 protected:
    DoubleRange range;
    double      value;

 public:
    DoubleOption(const char* c, const char* n, const char* d, double def = double(), DoubleRange r = DoubleRange(-HUGE_VAL, false, HUGE_VAL, false))
        : Option(n, d, c, "<double>"), range(r), value(def) {
        // FIXME: set LC_NUMERIC to "C" to make sure that strtof/strtod parses decimal point correctly.
    }

    operator      double   (void) const { return value; }
    operator      double&  (void)       { return value; }
    DoubleOption& operator=(double x)   { value = x; return *this; }

    virtual bool parse(const char* str){
        const char* span = str; 

        if (!match(span, "-") || !match(span, name) || !match(span, "="))
            return false;

        char*  end;
        double tmp = strtod(span, &end);

        if (end == NULL) 
            return false;
        else if (tmp >= range.end && (!range.end_inclusive || tmp != range.end)){
            fprintf(stderr, "ERROR! value <%s> is too large for option \"%s\".\n", span, name);
            exit(1);
        }else if (tmp <= range.begin && (!range.begin_inclusive || tmp != range.begin)){
            fprintf(stderr, "ERROR! value <%s> is too small for option \"%s\".\n", span, name);
            exit(1); }

        value = tmp;
        // fprintf(stderr, "READ VALUE: %g\n", value);

        return true;
    }

    virtual void help (bool verbose = false){
        fprintf(stderr, "  -%-12s = %-8s %c%4.2g .. %4.2g%c (default: %g)\n", 
                name, type_name, 
                range.begin_inclusive ? '[' : '(', 
                range.begin,
                range.end,
                range.end_inclusive ? ']' : ')', 
                value);
        if (verbose){
            fprintf(stderr, "\n        %s\n", description);
            fprintf(stderr, "\n");
        }
    }
};


//==================================================================================================
// Int options:


class IntOption : public Option
{
 protected:
    IntRange range;
    int32_t  value;

 public:
    IntOption(const char* c, const char* n, const char* d, int32_t def = int32_t(), IntRange r = IntRange(INT32_MIN, INT32_MAX))
        : Option(n, d, c, "<int32>"), range(r), value(def) {}
 
    operator   int32_t   (void) const { return value; }
    operator   int32_t&  (void)       { return value; }
    IntOption& operator= (int32_t x)  { value = x; return *this; }

    virtual bool parse(const char* str){
        const char* span = str; 

        if (!match(span, "-") || !match(span, name) || !match(span, "="))
            return false;

        char*   end;
        int32_t tmp = strtol(span, &end, 10);

        if (end == NULL) 
            return false;
        else if (tmp > range.end){
            fprintf(stderr, "ERROR! value <%s> is too large for option \"%s\".\n", span, name);
            exit(1);
        }else if (tmp < range.begin){
            fprintf(stderr, "ERROR! value <%s> is too small for option \"%s\".\n", span, name);
            exit(1); }

        value = tmp;

        return true;
    }

    virtual void help (bool verbose = false){
        fprintf(stderr, "  -%-12s = %-8s [", name, type_name);
        if (range.begin == INT32_MIN)
            fprintf(stderr, "imin");
        else
            fprintf(stderr, "%4d", range.begin);

        fprintf(stderr, " .. ");
        if (range.end == INT32_MAX)
            fprintf(stderr, "imax");
        else
            fprintf(stderr, "%4d", range.end);

        fprintf(stderr, "] (default: %d)\n", value);
        if (verbose){
            fprintf(stderr, "\n        %s\n", description);
            fprintf(stderr, "\n");
        }
    }
};


// Leave this out for visual C++ until Microsoft implements C99 and gets support for strtoll.
#ifndef _MSC_VER

class Int64Option : public Option
{
 protected:
    Int64Range range;
    int64_t  value;

 public:
    Int64Option(const char* c, const char* n, const char* d, int64_t def = int64_t(), Int64Range r = Int64Range(INT64_MIN, INT64_MAX))
        : Option(n, d, c, "<int64>"), range(r), value(def) {}
 
    operator     int64_t   (void) const { return value; }
    operator     int64_t&  (void)       { return value; }
    Int64Option& operator= (int64_t x)  { value = x; return *this; }

    virtual bool parse(const char* str){
        const char* span = str; 

        if (!match(span, "-") || !match(span, name) || !match(span, "="))
            return false;

        char*   end;
        int64_t tmp = strtoll(span, &end, 10);

        if (end == NULL) 
            return false;
        else if (tmp > range.end){
            fprintf(stderr, "ERROR! value <%s> is too large for option \"%s\".\n", span, name);
            exit(1);
        }else if (tmp < range.begin){
            fprintf(stderr, "ERROR! value <%s> is too small for option \"%s\".\n", span, name);
            exit(1); }

        value = tmp;

        return true;
    }

    virtual void help (bool verbose = false){
        fprintf(stderr, "  -%-12s = %-8s [", name, type_name);
        if (range.begin == INT64_MIN)
            fprintf(stderr, "imin");
        else
            fprintf(stderr, "%4" "ld", range.begin);

        fprintf(stderr, " .. ");
        if (range.end == INT64_MAX)
            fprintf(stderr, "imax");
        else
            fprintf(stderr, "%4" "ld", range.end);

        fprintf(stderr, "] (default: %" "ld" ")\n", value);
        if (verbose){
            fprintf(stderr, "\n        %s\n", description);
            fprintf(stderr, "\n");
        }
    }
};
#endif

//==================================================================================================
// String option:


class StringOption : public Option
{
    const char* value;
 public:
    StringOption(const char* c, const char* n, const char* d, const char* def = NULL) 
        : Option(n, d, c, "<string>"), value(def) {}

    operator      const char*  (void) const     { return value; }
    operator      const char*& (void)           { return value; }
    StringOption& operator=    (const char* x)  { value = x; return *this; }

    virtual bool parse(const char* str){
        const char* span = str; 

        if (!match(span, "-") || !match(span, name) || !match(span, "="))
            return false;

        value = span;
        return true;
    }

    virtual void help (bool verbose = false){
        fprintf(stderr, "  -%-10s = %8s\n", name, type_name);
        if (verbose){
            fprintf(stderr, "\n        %s\n", description);
            fprintf(stderr, "\n");
        }
    }    
};


//==================================================================================================
// Bool option:


class BoolOption : public Option
{
    bool value;

 public:
    BoolOption(const char* c, const char* n, const char* d, bool v) 
        : Option(n, d, c, "<bool>"), value(v) {}

    operator    bool     (void) const { return value; }
    operator    bool&    (void)       { return value; }
    BoolOption& operator=(bool b)     { value = b; return *this; }

    virtual bool parse(const char* str){
        const char* span = str; 
        
        if (match(span, "-")){
            bool b = !match(span, "no-");

            if (strcmp(span, name) == 0){
                value = b;
                return true; }
        }

        return false;
    }

    virtual void help (bool verbose = false){

        fprintf(stderr, "  -%s, -no-%s", name, name);

        for (uint32_t i = 0; i < 32 - strlen(name)*2; i++)
            fprintf(stderr, " ");

        fprintf(stderr, " ");
        fprintf(stderr, "(default: %s)\n", value ? "on" : "off");
        if (verbose){
            fprintf(stderr, "\n        %s\n", description);
            fprintf(stderr, "\n");
        }
    }
};

//=================================================================================================
}

#endif

/***********************************************************************************[SolverTypes.h]
 Glucose -- Copyright (c) 2009, Gilles Audemard, Laurent Simon
				CRIL - Univ. Artois, France
				LRI  - Univ. Paris Sud, France
 
Glucose sources are based on MiniSat (see below MiniSat copyrights). Permissions and copyrights of
Glucose are exactly the same as Minisat on which it is based on. (see below).

---------------
Copyright (c) 2003-2006, Niklas Een, Niklas Sorensson
Copyright (c) 2007-2010, Niklas Sorensson

Permission is hereby granted, free of charge, to any person obtaining a copy of this software and
associated documentation files (the "Software"), to deal in the Software without restriction,
including without limitation the rights to use, copy, modify, merge, publish, distribute,
sublicense, and/or sell copies of the Software, and to permit persons to whom the Software is
furnished to do so, subject to the following conditions:

The above copyright notice and this permission notice shall be included in all copies or
substantial portions of the Software.

THE SOFTWARE IS PROVIDED "AS IS", WITHOUT WARRANTY OF ANY KIND, EXPRESS OR IMPLIED, INCLUDING BUT
NOT LIMITED TO THE WARRANTIES OF MERCHANTABILITY, FITNESS FOR A PARTICULAR PURPOSE AND
NONINFRINGEMENT. IN NO EVENT SHALL THE AUTHORS OR COPYRIGHT HOLDERS BE LIABLE FOR ANY CLAIM,
DAMAGES OR OTHER LIABILITY, WHETHER IN AN ACTION OF CONTRACT, TORT OR OTHERWISE, ARISING FROM, OUT
OF OR IN CONNECTION WITH THE SOFTWARE OR THE USE OR OTHER DEALINGS IN THE SOFTWARE.
**************************************************************************************************/


#ifndef Glucose_SolverTypes_h
#define Glucose_SolverTypes_h

#include <assert.h>

/*******************************************************************************************[Map.h]
Copyright (c) 2006-2010, Niklas Sorensson

Permission is hereby granted, free of charge, to any person obtaining a copy of this software and
associated documentation files (the "Software"), to deal in the Software without restriction,
including without limitation the rights to use, copy, modify, merge, publish, distribute,
sublicense, and/or sell copies of the Software, and to permit persons to whom the Software is
furnished to do so, subject to the following conditions:

The above copyright notice and this permission notice shall be included in all copies or
substantial portions of the Software.

THE SOFTWARE IS PROVIDED "AS IS", WITHOUT WARRANTY OF ANY KIND, EXPRESS OR IMPLIED, INCLUDING BUT
NOT LIMITED TO THE WARRANTIES OF MERCHANTABILITY, FITNESS FOR A PARTICULAR PURPOSE AND
NONINFRINGEMENT. IN NO EVENT SHALL THE AUTHORS OR COPYRIGHT HOLDERS BE LIABLE FOR ANY CLAIM,
DAMAGES OR OTHER LIABILITY, WHETHER IN AN ACTION OF CONTRACT, TORT OR OTHERWISE, ARISING FROM, OUT
OF OR IN CONNECTION WITH THE SOFTWARE OR THE USE OR OTHER DEALINGS IN THE SOFTWARE.
**************************************************************************************************/

#ifndef Glucose_Map_h
#define Glucose_Map_h


namespace Glucose {

//=================================================================================================
// Default hash/equals functions
//

template<class K> struct Hash  { uint32_t operator()(const K& k)               const { return hash(k);  } };
template<class K> struct Equal { bool     operator()(const K& k1, const K& k2) const { return k1 == k2; } };

template<class K> struct DeepHash  { uint32_t operator()(const K* k)               const { return hash(*k);  } };
template<class K> struct DeepEqual { bool     operator()(const K* k1, const K* k2) const { return *k1 == *k2; } };

static inline uint32_t hash(uint32_t x){ return x; }
static inline uint32_t hash(uint64_t x){ return (uint32_t)x; }
static inline uint32_t hash(int32_t x) { return (uint32_t)x; }
static inline uint32_t hash(int64_t x) { return (uint32_t)x; }


//=================================================================================================
// Some primes
//

static const int nprimes          = 25;
static const int primes [nprimes] = { 31, 73, 151, 313, 643, 1291, 2593, 5233, 10501, 21013, 42073, 84181, 168451, 337219, 674701, 1349473, 2699299, 5398891, 10798093, 21596719, 43193641, 86387383, 172775299, 345550609, 691101253 };

//=================================================================================================
// Hash table implementation of Maps
//

template<class K, class D, class H = Hash<K>, class E = Equal<K> >
class Map {
 public:
    struct Pair { K key; D data; };

 private:
    H          hash;
    E          equals;

    vec<Pair>* table;
    int        cap;
    int        size;

    // Don't allow copying (error prone):
    Map<K,D,H,E>&  operator = (Map<K,D,H,E>& other) { assert(0); }
                   Map        (Map<K,D,H,E>& other) { assert(0); }

    bool    checkCap(int new_size) const { return new_size > cap; }

    int32_t index  (const K& k) const { return hash(k) % cap; }
    void   _insert (const K& k, const D& d) { 
        vec<Pair>& ps = table[index(k)];
        ps.push(); ps.last().key = k; ps.last().data = d; }

    void    rehash () {
        const vec<Pair>* old = table;

        int old_cap = cap;
        int newsize = primes[0];
        for (int i = 1; newsize <= cap && i < nprimes; i++)
           newsize = primes[i];

        table = new vec<Pair>[newsize];
        cap   = newsize;

        for (int i = 0; i < old_cap; i++){
            for (int j = 0; j < old[i].size(); j++){
                _insert(old[i][j].key, old[i][j].data); }}

        delete [] old;

        // printf(" --- rehashing, old-cap=%d, new-cap=%d\n", cap, newsize);
    }

    
 public:

    Map () : table(NULL), cap(0), size(0) {}
    Map (const H& h, const E& e) : hash(h), equals(e), table(NULL), cap(0), size(0){}
    ~Map () { delete [] table; }

    // PRECONDITION: the key must already exist in the map.
    const D& operator [] (const K& k) const
    {
        assert(size != 0);
        const D*         res = NULL;
        const vec<Pair>& ps  = table[index(k)];
        for (int i = 0; i < ps.size(); i++)
            if (equals(ps[i].key, k))
                res = &ps[i].data;
        assert(res != NULL);
        return *res;
    }

    // PRECONDITION: the key must already exist in the map.
    D& operator [] (const K& k)
    {
        assert(size != 0);
        D*         res = NULL;
        vec<Pair>& ps  = table[index(k)];
        for (int i = 0; i < ps.size(); i++)
            if (equals(ps[i].key, k))
                res = &ps[i].data;
        assert(res != NULL);
        return *res;
    }

    // PRECONDITION: the key must *NOT* exist in the map.
    void insert (const K& k, const D& d) { if (checkCap(size+1)) rehash(); _insert(k, d); size++; }
    bool peek   (const K& k, D& d) const {
        if (size == 0) return false;
        const vec<Pair>& ps = table[index(k)];
        for (int i = 0; i < ps.size(); i++)
            if (equals(ps[i].key, k)){
                d = ps[i].data;
                return true; } 
        return false;
    }

    bool has   (const K& k) const {
        if (size == 0) return false;
        const vec<Pair>& ps = table[index(k)];
        for (int i = 0; i < ps.size(); i++)
            if (equals(ps[i].key, k))
                return true;
        return false;
    }

    // PRECONDITION: the key must exist in the map.
    void remove(const K& k) {
        assert(table != NULL);
        vec<Pair>& ps = table[index(k)];
        int j = 0;
        for (; j < ps.size() && !equals(ps[j].key, k); j++);
        assert(j < ps.size());
        ps[j] = ps.last();
        ps.pop();
        size--;
    }

    void clear  () {
        cap = size = 0;
        delete [] table;
        table = NULL;
    }

    int  elems() const { return size; }
    int  bucket_count() const { return cap; }

    // NOTE: the hash and equality objects are not moved by this method:
    void moveTo(Map& other){
        delete [] other.table;

        other.table = table;
        other.cap   = cap;
        other.size  = size;

        table = NULL;
        size = cap = 0;
    }

    // NOTE: given a bit more time, I could make a more C++-style iterator out of this:
    const vec<Pair>& bucket(int i) const { return table[i]; }
};

//=================================================================================================
}

#endif

/*****************************************************************************************[Alloc.h]
Copyright (c) 2008-2010, Niklas Sorensson

Permission is hereby granted, free of charge, to any person obtaining a copy of this software and
associated documentation files (the "Software"), to deal in the Software without restriction,
including without limitation the rights to use, copy, modify, merge, publish, distribute,
sublicense, and/or sell copies of the Software, and to permit persons to whom the Software is
furnished to do so, subject to the following conditions:

The above copyright notice and this permission notice shall be included in all copies or
substantial portions of the Software.

THE SOFTWARE IS PROVIDED "AS IS", WITHOUT WARRANTY OF ANY KIND, EXPRESS OR IMPLIED, INCLUDING BUT
NOT LIMITED TO THE WARRANTIES OF MERCHANTABILITY, FITNESS FOR A PARTICULAR PURPOSE AND
NONINFRINGEMENT. IN NO EVENT SHALL THE AUTHORS OR COPYRIGHT HOLDERS BE LIABLE FOR ANY CLAIM,
DAMAGES OR OTHER LIABILITY, WHETHER IN AN ACTION OF CONTRACT, TORT OR OTHERWISE, ARISING FROM, OUT
OF OR IN CONNECTION WITH THE SOFTWARE OR THE USE OR OTHER DEALINGS IN THE SOFTWARE.
**************************************************************************************************/


#ifndef Glucose_Alloc_h
#define Glucose_Alloc_h

/****************************************************************************************[XAlloc.h]
Copyright (c) 2009-2010, Niklas Sorensson

Permission is hereby granted, free of charge, to any person obtaining a copy of this software and
associated documentation files (the "Software"), to deal in the Software without restriction,
including without limitation the rights to use, copy, modify, merge, publish, distribute,
sublicense, and/or sell copies of the Software, and to permit persons to whom the Software is
furnished to do so, subject to the following conditions:

The above copyright notice and this permission notice shall be included in all copies or
substantial portions of the Software.

THE SOFTWARE IS PROVIDED "AS IS", WITHOUT WARRANTY OF ANY KIND, EXPRESS OR IMPLIED, INCLUDING BUT
NOT LIMITED TO THE WARRANTIES OF MERCHANTABILITY, FITNESS FOR A PARTICULAR PURPOSE AND
NONINFRINGEMENT. IN NO EVENT SHALL THE AUTHORS OR COPYRIGHT HOLDERS BE LIABLE FOR ANY CLAIM,
DAMAGES OR OTHER LIABILITY, WHETHER IN AN ACTION OF CONTRACT, TORT OR OTHERWISE, ARISING FROM, OUT
OF OR IN CONNECTION WITH THE SOFTWARE OR THE USE OR OTHER DEALINGS IN THE SOFTWARE.
**************************************************************************************************/


#ifndef Glucose_XAlloc_h
#define Glucose_XAlloc_h

#include <errno.h>
#include <stdlib.h>
#include <stdio.h>

namespace Glucose {

//=================================================================================================
// Simple layer on top of malloc/realloc to catch out-of-memory situtaions and provide some typing:

class OutOfMemoryException{};
static inline void* xrealloc(void *ptr, size_t size)
{
    void* mem = realloc(ptr, size);
    if (mem == NULL && errno == ENOMEM){
        throw OutOfMemoryException();
    }else {
        return mem;
	}
}

//=================================================================================================
}

#endif



namespace Glucose {

//=================================================================================================
// Simple Region-based memory allocator:

template<class T>
class RegionAllocator
{
    T*        memory;
    uint32_t  sz;
    uint32_t  cap;
    uint32_t  wasted_;

    void capacity(uint32_t min_cap);

 public:
    // TODO: make this a class for better type-checking?
    typedef uint32_t Ref;
    enum { Ref_Undef = UINT32_MAX };
    enum { Unit_Size = sizeof(uint32_t) };

    explicit RegionAllocator(uint32_t start_cap = 1024*1024) : memory(NULL), sz(0), cap(0), wasted_(0){ capacity(start_cap); }
    ~RegionAllocator()
    {
        if (memory != NULL)
            ::free(memory);
    }


    uint32_t size      () const      { return sz; }
    uint32_t wasted    () const      { return wasted_; }

    Ref      alloc     (int size); 
    void     free      (int size)    { wasted_ += size; }

    // Deref, Load Effective Address (LEA), Inverse of LEA (AEL):
    T&       operator[](Ref r)       { assert(r >= 0 && r < sz); return memory[r]; }
    const T& operator[](Ref r) const { assert(r >= 0 && r < sz); return memory[r]; }

    T*       lea       (Ref r)       { assert(r >= 0 && r < sz); return &memory[r]; }
    const T* lea       (Ref r) const { assert(r >= 0 && r < sz); return &memory[r]; }
    Ref      ael       (const T* t)  { assert((void*)t >= (void*)&memory[0] && (void*)t < (void*)&memory[sz-1]);
        return  (Ref)(t - &memory[0]); }

    void     moveTo(RegionAllocator& to) {
        if (to.memory != NULL) ::free(to.memory);
        to.memory = memory;
        to.sz = sz;
        to.cap = cap;
        to.wasted_ = wasted_;

        memory = NULL;
        sz = cap = wasted_ = 0;
    }


};

template<class T>
void RegionAllocator<T>::capacity(uint32_t min_cap)
{
    if (cap >= min_cap) return;

    uint32_t prev_cap = cap;
    while (cap < min_cap){
        // NOTE: Multiply by a factor (13/8) without causing overflow, then add 2 and make the
        // result even by clearing the least significant bit. The resulting sequence of capacities
        // is carefully chosen to hit a maximum capacity that is close to the '2^32-1' limit when
        // using 'uint32_t' as indices so that as much as possible of this space can be used.
        uint32_t delta = ((cap >> 1) + (cap >> 3) + 2) & ~1;
        cap += delta;

        if (cap <= prev_cap)
            throw OutOfMemoryException();
    }
    //printf(" .. (%p) cap = %u\n", this, cap);

    assert(cap > 0);
    memory = (T*)xrealloc(memory, sizeof(T)*cap);
}


template<class T>
typename RegionAllocator<T>::Ref
RegionAllocator<T>::alloc(int size)
{ 
    //printf("ALLOC called (this = %p, size = %d)\n", this, size); fflush(stdout);
    assert(size > 0);
    capacity(sz + size);

    uint32_t prev_sz = sz;
    sz += size;
    
    // Handle overflow:
    if (sz < prev_sz)
        throw OutOfMemoryException();

    return prev_sz;
}


//=================================================================================================
}

#endif


namespace Glucose {

//=================================================================================================
// Variables, literals, lifted booleans, clauses:


// NOTE! Variables are just integers. No abstraction here. They should be chosen from 0..N,
// so that they can be used as array indices.

typedef int Var;
#define var_Undef (-1)


struct Lit {
    int     x;

    // Use this as a constructor:
    friend Lit mkLit(Var var, bool sign = false);

    bool operator == (Lit p) const { return x == p.x; }
    bool operator != (Lit p) const { return x != p.x; }
    bool operator <  (Lit p) const { return x < p.x;  } // '<' makes p, ~p adjacent in the ordering.
};


inline  Lit  mkLit     (Var var, bool sign) { Lit p; p.x = var + var + (int)sign; return p; }
inline  Lit  operator ~(Lit p)              { Lit q; q.x = p.x ^ 1; return q; }
inline  Lit  operator ^(Lit p, bool b)      { Lit q; q.x = p.x ^ (unsigned int)b; return q; }
inline  bool sign      (Lit p)              { return p.x & 1; }
inline  int  var       (Lit p)              { return p.x >> 1; }

// Mapping Literals to and from compact integers suitable for array indexing:
inline  int  toInt     (Var v)              { return v; } 
inline  int  toInt     (Lit p)              { return p.x; } 
inline  Lit  toLit     (int i)              { Lit p; p.x = i; return p; } 

//const Lit lit_Undef = mkLit(var_Undef, false);  // }- Useful special constants.
//const Lit lit_Error = mkLit(var_Undef, true );  // }

const Lit lit_Undef = { -2 };  // }- Useful special constants.
const Lit lit_Error = { -1 };  // }


//=================================================================================================
// Lifted booleans:
//
// NOTE: this implementation is optimized for the case when comparisons between values are mostly
//       between one variable and one constant. Some care had to be taken to make sure that gcc 
//       does enough constant propagation to produce sensible code, and this appears to be somewhat
//       fragile unfortunately.

#define l_True  (Glucose::lbool((uint8_t)0)) // gcc does not do constant propagation if these are real constants.
#define l_False (Glucose::lbool((uint8_t)1))
#define l_Undef (Glucose::lbool((uint8_t)2))

class lbool {
    uint8_t value;

public:
    explicit lbool(uint8_t v) : value(v) { }

    lbool()       : value(0) { }
    explicit lbool(bool x) : value(!x) { }

    bool  operator == (lbool b) const { return ((b.value&2) & (value&2)) | (!(b.value&2)&(value == b.value)); }
    bool  operator != (lbool b) const { return !(*this == b); }
    lbool operator ^  (bool  b) const { return lbool((uint8_t)(value^(uint8_t)b)); }

    lbool operator && (lbool b) const { 
        uint8_t sel = (this->value << 1) | (b.value << 3);
        uint8_t v   = (0xF7F755F4 >> sel) & 3;
        return lbool(v); }

    lbool operator || (lbool b) const {
        uint8_t sel = (this->value << 1) | (b.value << 3);
        uint8_t v   = (0xFCFCF400 >> sel) & 3;
        return lbool(v); }

    friend int   toInt  (lbool l);
    friend lbool toLbool(int   v);
};
inline int   toInt  (lbool l) { return l.value; }
inline lbool toLbool(int   v) { return lbool((uint8_t)v);  }

//=================================================================================================
// Clause -- a simple class for representing a clause:

class Clause;
typedef RegionAllocator<uint32_t>::Ref CRef;

class Clause {
    struct {
      unsigned mark      : 2;
      unsigned learnt    : 1;
      unsigned has_extra : 1;
      unsigned reloced   : 1;
      unsigned lbd       : 26;
      unsigned canbedel  : 1;
      unsigned size      : 32;
      unsigned szWithoutSelectors : 32;

    }                            header;
    union { Lit lit; float act; uint32_t abs; CRef rel; } data[0];

    friend class ClauseAllocator;

    // NOTE: This constructor cannot be used directly (doesn't allocate enough memory).
    template<class V>
    Clause(const V& ps, bool use_extra, bool learnt) {
        header.mark      = 0;
        header.learnt    = learnt;
        header.has_extra = use_extra;
        header.reloced   = 0;
        header.size      = ps.size();
	header.lbd = 0;
	header.canbedel = 1;
        for (int i = 0; i < ps.size(); i++) 
            data[i].lit = ps[i];
	
        if (header.has_extra){
	  if (header.learnt) 
                data[header.size].act = 0; 
            else 
                calcAbstraction(); }
    }

public:
    void calcAbstraction() {
        assert(header.has_extra);
        uint32_t abstraction = 0;
        for (int i = 0; i < size(); i++)
            abstraction |= 1 << (var(data[i].lit) & 31);
        data[header.size].abs = abstraction;  }


    int          size        ()      const   { return header.size; }
    void         shrink      (int i)         { assert(i <= size()); if (header.has_extra) data[header.size-i] = data[header.size]; header.size -= i; }
    void         pop         ()              { shrink(1); }
    bool         learnt      ()      const   { return header.learnt; }
    bool         has_extra   ()      const   { return header.has_extra; }
    uint32_t     mark        ()      const   { return header.mark; }
    void         mark        (uint32_t m)    { header.mark = m; }
    const Lit&   last        ()      const   { return data[header.size-1].lit; }

    bool         reloced     ()      const   { return header.reloced; }
    CRef         relocation  ()      const   { return data[0].rel; }
    void         relocate    (CRef c)        { header.reloced = 1; data[0].rel = c; }

    // NOTE: somewhat unsafe to change the clause in-place! Must manually call 'calcAbstraction' afterwards for
    //       subsumption operations to behave correctly.
    Lit&         operator [] (int i)         { return data[i].lit; }
    Lit          operator [] (int i) const   { return data[i].lit; }
    operator const Lit* (void) const         { return (Lit*)data; }

    float&       activity    ()              { assert(header.has_extra); return data[header.size].act; }
    uint32_t     abstraction () const        { assert(header.has_extra); return data[header.size].abs; }

    Lit          subsumes    (const Clause& other) const;
    void         strengthen  (Lit p);
    void         setLBD(int i)  {header.lbd = i;} 
    // unsigned int&       lbd    ()              { return header.lbd; }
    unsigned int        lbd    () const        { return header.lbd; }
    void setCanBeDel(bool b) {header.canbedel = b;}
    bool canBeDel() {return header.canbedel;}
    void setSizeWithoutSelectors   (unsigned int n)              {header.szWithoutSelectors = n; }
    unsigned int        sizeWithoutSelectors   () const        { return header.szWithoutSelectors; }

};


//=================================================================================================
// ClauseAllocator -- a simple class for allocating memory for clauses:


const CRef CRef_Undef = RegionAllocator<uint32_t>::Ref_Undef;
class ClauseAllocator : public RegionAllocator<uint32_t>
{
    static int clauseWord32Size(int size, bool has_extra){
        return (sizeof(Clause) + (sizeof(Lit) * (size + (int)has_extra))) / sizeof(uint32_t); }
 public:
    bool extra_clause_field;

    ClauseAllocator(uint32_t start_cap) : RegionAllocator<uint32_t>(start_cap), extra_clause_field(false){}
    ClauseAllocator() : extra_clause_field(false){}

    void moveTo(ClauseAllocator& to){
        to.extra_clause_field = extra_clause_field;
        RegionAllocator<uint32_t>::moveTo(to); }

    template<class Lits>
    CRef alloc(const Lits& ps, bool learnt = false)
    {
        assert(sizeof(Lit)      == sizeof(uint32_t));
        assert(sizeof(float)    == sizeof(uint32_t));
        bool use_extra = learnt | extra_clause_field;

        CRef cid = RegionAllocator<uint32_t>::alloc(clauseWord32Size(ps.size(), use_extra));
        new (lea(cid)) Clause(ps, use_extra, learnt);

        return cid;
    }

    // Deref, Load Effective Address (LEA), Inverse of LEA (AEL):
    Clause&       operator[](Ref r)       { return (Clause&)RegionAllocator<uint32_t>::operator[](r); }
    const Clause& operator[](Ref r) const { return (Clause&)RegionAllocator<uint32_t>::operator[](r); }
    Clause*       lea       (Ref r)       { return (Clause*)RegionAllocator<uint32_t>::lea(r); }
    const Clause* lea       (Ref r) const { return (Clause*)RegionAllocator<uint32_t>::lea(r); }
    Ref           ael       (const Clause* t){ return RegionAllocator<uint32_t>::ael((uint32_t*)t); }

    void free(CRef cid)
    {
        Clause& c = operator[](cid);
        RegionAllocator<uint32_t>::free(clauseWord32Size(c.size(), c.has_extra()));
    }

    void reloc(CRef& cr, ClauseAllocator& to)
    {
        Clause& c = operator[](cr);
        
        if (c.reloced()) { cr = c.relocation(); return; }
        
        cr = to.alloc(c, c.learnt());
        c.relocate(cr);
        
        // Copy extra data-fields: 
        // (This could be cleaned-up. Generalize Clause-constructor to be applicable here instead?)
        to[cr].mark(c.mark());
        if (to[cr].learnt())        {
	  to[cr].activity() = c.activity();
	  to[cr].setLBD(c.lbd());
	  to[cr].setSizeWithoutSelectors(c.sizeWithoutSelectors());
	  to[cr].setCanBeDel(c.canBeDel());
	}
        else if (to[cr].has_extra()) to[cr].calcAbstraction();
    }
};


//=================================================================================================
// OccLists -- a class for maintaining occurence lists with lazy deletion:

template<class Idx, class Vec, class Deleted>
class OccLists
{
    vec<Vec>  occs;
    vec<char> dirty;
    vec<Idx>  dirties;
    Deleted   deleted;

 public:
    OccLists(const Deleted& d) : deleted(d) {}
    
    void  init      (const Idx& idx){ occs.growTo(toInt(idx)+1); dirty.growTo(toInt(idx)+1, 0); }
    // Vec&  operator[](const Idx& idx){ return occs[toInt(idx)]; }
    Vec&  operator[](const Idx& idx){ return occs[toInt(idx)]; }
    Vec&  lookup    (const Idx& idx){ if (dirty[toInt(idx)]) clean(idx); return occs[toInt(idx)]; }

    void  cleanAll  ();
    void  clean     (const Idx& idx);
    void  smudge    (const Idx& idx){
        if (dirty[toInt(idx)] == 0){
            dirty[toInt(idx)] = 1;
            dirties.push(idx);
        }
    }

    void  clear(bool free = true){
        occs   .clear(free);
        dirty  .clear(free);
        dirties.clear(free);
    }
};


template<class Idx, class Vec, class Deleted>
void OccLists<Idx,Vec,Deleted>::cleanAll()
{
    for (int i = 0; i < dirties.size(); i++)
        // Dirties may contain duplicates so check here if a variable is already cleaned:
        if (dirty[toInt(dirties[i])])
            clean(dirties[i]);
    dirties.clear();
}


template<class Idx, class Vec, class Deleted>
void OccLists<Idx,Vec,Deleted>::clean(const Idx& idx)
{
    Vec& vec = occs[toInt(idx)];
    int  i, j;
    for (i = j = 0; i < vec.size(); i++)
        if (!deleted(vec[i]))
            vec[j++] = vec[i];
    vec.shrink(i - j);
    dirty[toInt(idx)] = 0;
}


//=================================================================================================
// CMap -- a class for mapping clauses to values:


template<class T>
class CMap
{
    struct CRefHash {
        uint32_t operator()(CRef cr) const { return (uint32_t)cr; } };

    typedef Map<CRef, T, CRefHash> HashTable;
    HashTable map;
        
 public:
    // Size-operations:
    void     clear       ()                           { map.clear(); }
    int      size        ()                const      { return map.elems(); }

    
    // Insert/Remove/Test mapping:
    void     insert      (CRef cr, const T& t){ map.insert(cr, t); }
    void     growTo      (CRef cr, const T& t){ map.insert(cr, t); } // NOTE: for compatibility
    void     remove      (CRef cr)            { map.remove(cr); }
    bool     has         (CRef cr, T& t)      { return map.peek(cr, t); }

    // Vector interface (the clause 'c' must already exist):
    const T& operator [] (CRef cr) const      { return map[cr]; }
    T&       operator [] (CRef cr)            { return map[cr]; }

    // Iteration (not transparent at all at the moment):
    int  bucket_count() const { return map.bucket_count(); }
    const vec<typename HashTable::Pair>& bucket(int i) const { return map.bucket(i); }

    // Move contents to other map:
    void moveTo(CMap& other){ map.moveTo(other.map); }

    // TMP debug:
    void debug(){
        printf(" --- size = %d, bucket_count = %d\n", size(), map.bucket_count()); }
};


/*_________________________________________________________________________________________________
|
|  subsumes : (other : const Clause&)  ->  Lit
|  
|  Description:
|       Checks if clause subsumes 'other', and at the same time, if it can be used to simplify 'other'
|       by subsumption resolution.
|  
|    Result:
|       lit_Error  - No subsumption or simplification
|       lit_Undef  - Clause subsumes 'other'
|       p          - The literal p can be deleted from 'other'
|________________________________________________________________________________________________@*/
inline Lit Clause::subsumes(const Clause& other) const
{
    //if (other.size() < size() || (extra.abst & ~other.extra.abst) != 0)
    //if (other.size() < size() || (!learnt() && !other.learnt() && (extra.abst & ~other.extra.abst) != 0))
    assert(!header.learnt);   assert(!other.header.learnt);
    assert(header.has_extra); assert(other.header.has_extra);
    if (other.header.size < header.size || (data[header.size].abs & ~other.data[other.header.size].abs) != 0)
        return lit_Error;

    Lit        ret = lit_Undef;
    const Lit* c   = (const Lit*)(*this);
    const Lit* d   = (const Lit*)other;

    for (unsigned i = 0; i < header.size; i++) {
        // search for c[i] or ~c[i]
        for (unsigned j = 0; j < other.header.size; j++)
            if (c[i] == d[j])
                goto ok;
            else if (ret == lit_Undef && c[i] == ~d[j]){
                ret = c[i];
                goto ok;
            }

        // did not find it
        return lit_Error;
    ok:;
    }

    return ret;
}

inline void Clause::strengthen(Lit p)
{
    remove(*this, p);
    calcAbstraction();
}
 
//=================================================================================================
}

 
#endif

/***********************************************************************************[BoundedQueue.h]
 Glucose -- Copyright (c) 2009, Gilles Audemard, Laurent Simon
                CRIL - Univ. Artois, France
                LRI  - Univ. Paris Sud, France
 
Permission is hereby granted, free of charge, to any person obtaining a copy of this software and
associated documentation files (the "Software"), to deal in the Software without restriction,
including without limitation the rights to use, copy, modify, merge, publish, distribute,
sublicense, and/or sell copies of the Software, and to permit persons to whom the Software is
furnished to do so, subject to the following conditions:
 
The above copyright notice and this permission notice shall be included in all copies or
substantial portions of the Software.
 
THE SOFTWARE IS PROVIDED "AS IS", WITHOUT WARRANTY OF ANY KIND, EXPRESS OR IMPLIED, INCLUDING BUT
NOT LIMITED TO THE WARRANTIES OF MERCHANTABILITY, FITNESS FOR A PARTICULAR PURPOSE AND
NONINFRINGEMENT. IN NO EVENT SHALL THE AUTHORS OR COPYRIGHT HOLDERS BE LIABLE FOR ANY CLAIM,
DAMAGES OR OTHER LIABILITY, WHETHER IN AN ACTION OF CONTRACT, TORT OR OTHERWISE, ARISING FROM, OUT
OF OR IN CONNECTION WITH THE SOFTWARE OR THE USE OR OTHER DEALINGS IN THE SOFTWARE.
**************************************************************************************************/


#ifndef BoundedQueue_h
#define BoundedQueue_h


//=================================================================================================

namespace Glucose {

template <class T>
class bqueue {
    vec<T>  elems;
    int     first;
	int		last;
	unsigned long long sumofqueue;
	int     maxsize;
	int     queuesize; // Number of current elements (must be < maxsize !)
	bool expComputed;
	double exp,value;
public:
 bqueue(void) : first(0), last(0), sumofqueue(0), maxsize(0), queuesize(0),expComputed(false) { } 
	
	void initSize(int size) {growTo(size);exp = 2.0/(size+1);} // Init size of bounded size queue
	
	void push(T x) {
	  expComputed = false;
		if (queuesize==maxsize) {
			assert(last==first); // The queue is full, next value to enter will replace oldest one
			sumofqueue -= elems[last];
			if ((++last) == maxsize) last = 0;
		} else 
			queuesize++;
		sumofqueue += x;
		elems[first] = x;
		if ((++first) == maxsize) {first = 0;last = 0;}
	}

	T peek() { assert(queuesize>0); return elems[last]; }
	void pop() {sumofqueue-=elems[last]; queuesize--; if ((++last) == maxsize) last = 0;}
	
	unsigned long long getsum() const {return sumofqueue;}
	unsigned int getavg() const {return (unsigned int)(sumofqueue/((unsigned long long)queuesize));}
	int maxSize() const {return maxsize;}
	double getavgDouble() const {
	  double tmp = 0;
	  for(int i=0;i<elems.size();i++) {
	    tmp+=elems[i];
	  }
	  return tmp/elems.size();
	}
	int isvalid() const {return (queuesize==maxsize);}
	
	void growTo(int size) {
		elems.growTo(size); 
		first=0; maxsize=size; queuesize = 0;last = 0;
		for(int i=0;i<size;i++) elems[i]=0; 
	}
	
	double getAvgExp() {
	  if(expComputed) return value;
	  double a=exp;
	  value = elems[first];
	  for(int i  = first;i<maxsize;i++) {
	    value+=a*((double)elems[i]);
	    a=a*exp;
	  }
	  for(int i  = 0;i<last;i++) {
	    value+=a*((double)elems[i]);
	    a=a*exp;
	  }
	  value = value*(1-exp)/(1-a);
	  expComputed = true;
	  return value;
	  

	}
	void fastclear() {first = 0; last = 0; queuesize=0; sumofqueue=0;} // to be called after restarts... Discard the queue
	
    int  size(void)    { return queuesize; }

    void clear(bool dealloc = false)   { elems.clear(dealloc); first = 0; maxsize=0; queuesize=0;sumofqueue=0;}


};
}
//=================================================================================================

#endif

/************************************************************************************[Constants.h]
 Glucose -- Copyright (c) 2009, Gilles Audemard, Laurent Simon
                CRIL - Univ. Artois, France
                LRI  - Univ. Paris Sud, France
  
Permission is hereby granted, free of charge, to any person obtaining a copy of this software and
associated documentation files (the "Software"), to deal in the Software without restriction,
including without limitation the rights to use, copy, modify, merge, publish, distribute,
sublicense, and/or sell copies of the Software, and to permit persons to whom the Software is
furnished to do so, subject to the following conditions:
 
The above copyright notice and this permission notice shall be included in all copies or
substantial portions of the Software.
 
THE SOFTWARE IS PROVIDED "AS IS", WITHOUT WARRANTY OF ANY KIND, EXPRESS OR IMPLIED, INCLUDING BUT
NOT LIMITED TO THE WARRANTIES OF MERCHANTABILITY, FITNESS FOR A PARTICULAR PURPOSE AND
NONINFRINGEMENT. IN NO EVENT SHALL THE AUTHORS OR COPYRIGHT HOLDERS BE LIABLE FOR ANY CLAIM,
DAMAGES OR OTHER LIABILITY, WHETHER IN AN ACTION OF CONTRACT, TORT OR OTHERWISE, ARISING FROM, OUT
OF OR IN CONNECTION WITH THE SOFTWARE OR THE USE OR OTHER DEALINGS IN THE SOFTWARE.
**************************************************************************************************/

#define DYNAMICNBLEVEL
#define CONSTANTREMOVECLAUSE
#define UPDATEVARACTIVITY

// Constants for clauses reductions
#define RATIOREMOVECLAUSES 2



// Constants for restarts
#define LOWER_BOUND_FOR_BLOCKING_RESTART 10000




namespace Glucose {

//=================================================================================================
// Solver -- the main class:

class Solver {
public:

    // Constructor/Destructor:
    //
    Solver();
    virtual ~Solver();

    // Problem specification:
    //
    Var     newVar    (bool polarity = true, bool dvar = true); // Add a new variable with parameters specifying variable mode.

    bool    addClause (const vec<Lit>& ps);                     // Add a clause to the solver. 
    bool    addEmptyClause();                                   // Add the empty clause, making the solver contradictory.
    bool    addClause (Lit p);                                  // Add a unit clause to the solver. 
    bool    addClause (Lit p, Lit q);                           // Add a binary clause to the solver. 
    bool    addClause (Lit p, Lit q, Lit r);                    // Add a ternary clause to the solver. 
    bool    addClause_(      vec<Lit>& ps);                     // Add a clause to the solver without making superflous internal copy. Will
                                                                // change the passed vector 'ps'.

    // Solving:
    //
    bool    simplify     ();                        // Removes already satisfied clauses.
    bool    solve        (const vec<Lit>& assumps); // Search for a model that respects a given set of assumptions.
    lbool   solveLimited (const vec<Lit>& assumps); // Search for a model that respects a given set of assumptions (With resource constraints).
    bool    solve        ();                        // Search without assumptions.
    bool    solve        (Lit p);                   // Search for a model that respects a single assumption.
    bool    solve        (Lit p, Lit q);            // Search for a model that respects two assumptions.
    bool    solve        (Lit p, Lit q, Lit r);     // Search for a model that respects three assumptions.
    bool    okay         () const;                  // FALSE means solver is in a conflicting state

    void    toDimacs     (FILE* f, const vec<Lit>& assumps);            // Write CNF to file in DIMACS-format.
    void    toDimacs     (const char *file, const vec<Lit>& assumps);
    void    toDimacs     (FILE* f, Clause& c, vec<Var>& map, Var& max);
    void printLit(Lit l);
    void printClause(CRef c);
    void printInitialClause(CRef c);
    // Convenience versions of 'toDimacs()':
    void    toDimacs     (const char* file);
    void    toDimacs     (const char* file, Lit p);
    void    toDimacs     (const char* file, Lit p, Lit q);
    void    toDimacs     (const char* file, Lit p, Lit q, Lit r);
    
    // Variable mode:
    // 
    void    setPolarity    (Var v, bool b); // Declare which polarity the decision heuristic should use for a variable. Requires mode 'polarity_user'.
    void    setDecisionVar (Var v, bool b); // Declare if a variable should be eligible for selection in the decision heuristic.

    // Read state:
    //
    lbool   value      (Var x) const;       // The current value of a variable.
    lbool   value      (Lit p) const;       // The current value of a literal.
    lbool   modelValue (Var x) const;       // The value of a variable in the last model. The last call to solve must have been satisfiable.
    lbool   modelValue (Lit p) const;       // The value of a literal in the last model. The last call to solve must have been satisfiable.
    int     nAssigns   ()      const;       // The current number of assigned literals.
    int     nClauses   ()      const;       // The current number of original clauses.
    int     nLearnts   ()      const;       // The current number of learnt clauses.
    int     nVars      ()      const;       // The current number of variables.
    int     nFreeVars  ()      const;

    // Incremental mode
    void setIncrementalMode();
    void initNbInitialVars(int nb);
    void printIncrementalStats();

    // Resource contraints:
    //
    void    setConfBudget(int64_t x);
    void    setPropBudget(int64_t x);
    void    budgetOff();
    void    interrupt();          // Trigger a (potentially asynchronous) interruption of the solver.
    void    clearInterrupt();     // Clear interrupt indicator flag.

    // Memory managment:
    //
    virtual void garbageCollect();
    void    checkGarbage(double gf);
    void    checkGarbage();




    // Extra results: (read-only member variable)
    //
    vec<lbool> model;             // If problem is satisfiable, this vector contains the model (if any).
    vec<Lit>   conflict;          // If problem is unsatisfiable (possibly under assumptions),
                                  // this vector represent the final conflict clause expressed in the assumptions.

    // Mode of operation:
    //
    int       verbosity;
    int       verbEveryConflicts;
    int       showModel;
    // Constants For restarts
    double    K;
    double    R;
    double    sizeLBDQueue;
    double    sizeTrailQueue;

    // Constants for reduce DB
    int firstReduceDB;
    int incReduceDB;
    int specialIncReduceDB;
    unsigned int lbLBDFrozenClause;

    // Constant for reducing clause
    int lbSizeMinimizingClause;
    unsigned int lbLBDMinimizingClause;

    double    var_decay;
    double    clause_decay;
    double    random_var_freq;
    double    random_seed;
    int       ccmin_mode;         // Controls conflict clause minimization (0=none, 1=basic, 2=deep).
    int       phase_saving;       // Controls the level of phase saving (0=none, 1=limited, 2=full).
    bool      rnd_pol;            // Use random polarities for branching heuristics.
    bool      rnd_init_act;       // Initialize variable activities with a small random value.
    double    garbage_frac;       // The fraction of wasted memory allowed before a garbage collection is triggered.

    // Certified UNSAT ( Thanks to Marijn Heule)
    FILE*               certifiedOutput;
    bool                certifiedUNSAT;

    
    // Statistics: (read-only member variable)
    //
    uint64_t nbRemovedClauses,nbReducedClauses,nbDL2,nbBin,nbUn,nbReduceDB,solves, starts, decisions, rnd_decisions, propagations, conflicts,conflictsRestarts,nbstopsrestarts,nbstopsrestartssame,lastblockatrestart;
    uint64_t dec_vars, clauses_literals, learnts_literals, max_literals, tot_literals;

protected:
    long curRestart;
    // Helper structures:
    //
    struct VarData { CRef reason; int level; };
    static inline VarData mkVarData(CRef cr, int l){ VarData d = {cr, l}; return d; }

    struct Watcher {
        CRef cref;
        Lit  blocker;
        Watcher(CRef cr, Lit p) : cref(cr), blocker(p) {}
        bool operator==(const Watcher& w) const { return cref == w.cref; }
        bool operator!=(const Watcher& w) const { return cref != w.cref; }
    };

    struct WatcherDeleted
    {
        const ClauseAllocator& ca;
        WatcherDeleted(const ClauseAllocator& _ca) : ca(_ca) {}
        bool operator()(const Watcher& w) const { return ca[w.cref].mark() == 1; }
    };

    struct VarOrderLt {
        const vec<double>&  activity;
        bool operator () (Var x, Var y) const { return activity[x] > activity[y]; }
        VarOrderLt(const vec<double>&  act) : activity(act) { }
    };


    // Solver state:
    //
    int lastIndexRed;
    bool                ok;               // If FALSE, the constraints are already unsatisfiable. No part of the solver state may be used!
    double              cla_inc;          // Amount to bump next clause with.
    vec<double>         activity;         // A heuristic measurement of the activity of a variable.
    double              var_inc;          // Amount to bump next variable with.
    OccLists<Lit, vec<Watcher>, WatcherDeleted>
                        watches;          // 'watches[lit]' is a list of constraints watching 'lit' (will go there if literal becomes true).
    OccLists<Lit, vec<Watcher>, WatcherDeleted>
                        watchesBin;          // 'watches[lit]' is a list of constraints watching 'lit' (will go there if literal becomes true).
    vec<CRef>           clauses;          // List of problem clauses.
    vec<CRef>           learnts;          // List of learnt clauses.

    vec<lbool>          assigns;          // The current assignments.
    vec<char>           polarity;         // The preferred polarity of each variable.
    vec<char>           decision;         // Declares if a variable is eligible for selection in the decision heuristic.
    vec<Lit>            trail;            // Assignment stack; stores all assigments made in the order they were made.
        vec<int>            nbpos;
    vec<int>            trail_lim;        // Separator indices for different decision levels in 'trail'.
    vec<VarData>        vardata;          // Stores reason and level for each variable.
    int                 qhead;            // Head of queue (as index into the trail -- no more explicit propagation queue in MiniSat).
    int                 simpDB_assigns;   // Number of top-level assignments since last execution of 'simplify()'.
    int64_t             simpDB_props;     // Remaining number of propagations that must be made before next execution of 'simplify()'.
    vec<Lit>            assumptions;      // Current set of assumptions provided to solve by the user.
    Heap<VarOrderLt>    order_heap;       // A priority queue of variables ordered with respect to the variable activity.
    double              progress_estimate;// Set by 'search()'.
    bool                remove_satisfied; // Indicates whether possibly inefficient linear scan for satisfied clauses should be performed in 'simplify'.
    vec<unsigned int> permDiff;      // permDiff[var] contains the current conflict number... Used to count the number of  LBD
    
#ifdef UPDATEVARACTIVITY
    // UPDATEVARACTIVITY trick (see competition'09 companion paper)
    vec<Lit> lastDecisionLevel; 
#endif

    ClauseAllocator     ca;

    int nbclausesbeforereduce;            // To know when it is time to reduce clause database
    
    bqueue<unsigned int> trailQueue,lbdQueue; // Bounded queues for restarts.
    float sumLBD; // used to compute the global average of LBD. Restarts...
    int sumAssumptions;


    // Temporaries (to reduce allocation overhead). Each variable is prefixed by the method in which it is
    // used, exept 'seen' wich is used in several places.
    //
    vec<char>           seen;
    vec<Lit>            analyze_stack;
    vec<Lit>            analyze_toclear;
    vec<Lit>            add_tmp;
    unsigned int  MYFLAG;


    double              max_learnts;
    double              learntsize_adjust_confl;
    int                 learntsize_adjust_cnt;

    // Resource contraints:
    //
    int64_t             conflict_budget;    // -1 means no budget.
    int64_t             propagation_budget; // -1 means no budget.
    bool                asynch_interrupt;


    // Variables added for incremental mode
    int incremental; // Use incremental SAT Solver
    int nbVarsInitialFormula; // nb VAR in formula without assumptions (incremental SAT)
    double totalTime4Sat,totalTime4Unsat;
    int nbSatCalls,nbUnsatCalls;
    vec<int> assumptionPositions,initialPositions;


    // Main internal methods:
    //
    void     insertVarOrder   (Var x);                                                 // Insert a variable in the decision order priority queue.
    Lit      pickBranchLit    ();                                                      // Return the next decision variable.
    void     newDecisionLevel ();                                                      // Begins a new decision level.
    void     uncheckedEnqueue (Lit p, CRef from = CRef_Undef);                         // Enqueue a literal. Assumes value of literal is undefined.
    bool     enqueue          (Lit p, CRef from = CRef_Undef);                         // Test if fact 'p' contradicts current state, enqueue otherwise.
    CRef     propagate        ();                                                      // Perform unit propagation. Returns possibly conflicting clause.
    void     cancelUntil      (int level);                                             // Backtrack until a certain level.
    void     analyze          (CRef confl, vec<Lit>& out_learnt, vec<Lit> & selectors, int& out_btlevel,unsigned int &nblevels,unsigned int &szWithoutSelectors);    // (bt = backtrack)
    void     analyzeFinal     (Lit p, vec<Lit>& out_conflict);                         // COULD THIS BE IMPLEMENTED BY THE ORDINARIY "analyze" BY SOME REASONABLE GENERALIZATION?
    bool     litRedundant     (Lit p, uint32_t abstract_levels);                       // (helper method for 'analyze()')
    lbool    search           (int nof_conflicts);                                     // Search for a given number of conflicts.
    lbool    solve_           ();                                                      // Main solve method (assumptions given in 'assumptions').
    void     reduceDB         ();                                                      // Reduce the set of learnt clauses.
    void     removeSatisfied  (vec<CRef>& cs);                                         // Shrink 'cs' to contain only non-satisfied clauses.
    void     rebuildOrderHeap ();

    // Maintaining Variable/Clause activity:
    //
    void     varDecayActivity ();                      // Decay all variables with the specified factor. Implemented by increasing the 'bump' value instead.
    void     varBumpActivity  (Var v, double inc);     // Increase a variable with the current 'bump' value.
    void     varBumpActivity  (Var v);                 // Increase a variable with the current 'bump' value.
    void     claDecayActivity ();                      // Decay all clauses with the specified factor. Implemented by increasing the 'bump' value instead.
    void     claBumpActivity  (Clause& c);             // Increase a clause with the current 'bump' value.

    // Operations on clauses:
    //
    void     attachClause     (CRef cr);               // Attach a clause to watcher lists.
    void     detachClause     (CRef cr, bool strict = false); // Detach a clause to watcher lists.
    void     removeClause     (CRef cr);               // Detach and free a clause.
    bool     locked           (const Clause& c) const; // Returns TRUE if a clause is a reason for some implication in the current state.
    bool     satisfied        (const Clause& c) const; // Returns TRUE if a clause is satisfied in the current state.

    unsigned int computeLBD(const vec<Lit> & lits,int end=-1);
    unsigned int computeLBD(const Clause &c);
    void minimisationWithBinaryResolution(vec<Lit> &out_learnt);

    void     relocAll         (ClauseAllocator& to);

    // Misc:
    //
    int      decisionLevel    ()      const; // Gives the current decisionlevel.
    uint32_t abstractLevel    (Var x) const; // Used to represent an abstraction of sets of decision levels.
    CRef     reason           (Var x) const;
    int      level            (Var x) const;
    double   progressEstimate ()      const; // DELETE THIS ?? IT'S NOT VERY USEFUL ...
    bool     withinBudget     ()      const;
    inline bool isSelector(Var v) {return (incremental && v>nbVarsInitialFormula);}

    // Static helpers:
    //

    // Returns a random float 0 <= x < 1. Seed must never be 0.
    static inline double drand(double& seed) {
        seed *= 1389796;
        int q = (int)(seed / 2147483647);
        seed -= (double)q * 2147483647;
        return seed / 2147483647; }

    // Returns a random integer 0 <= x < size. Seed must never be 0.
    static inline int irand(double& seed, int size) {
        return (int)(drand(seed) * size); }
};


//=================================================================================================
// Implementation of inline methods:

inline CRef Solver::reason(Var x) const { return vardata[x].reason; }
inline int  Solver::level (Var x) const { return vardata[x].level; }

inline void Solver::insertVarOrder(Var x) {
    if (!order_heap.inHeap(x) && decision[x]) order_heap.insert(x); }

inline void Solver::varDecayActivity() { var_inc *= (1 / var_decay); }
inline void Solver::varBumpActivity(Var v) { varBumpActivity(v, var_inc); }
inline void Solver::varBumpActivity(Var v, double inc) {
    if ( (activity[v] += inc) > 1e100 ) {
        // Rescale:
        for (int i = 0; i < nVars(); i++)
            activity[i] *= 1e-100;
        var_inc *= 1e-100; }

    // Update order_heap with respect to new activity:
    if (order_heap.inHeap(v))
        order_heap.decrease(v); }

inline void Solver::claDecayActivity() { cla_inc *= (1 / clause_decay); }
inline void Solver::claBumpActivity (Clause& c) {
        if ( (c.activity() += cla_inc) > 1e20 ) {
            // Rescale:
            for (int i = 0; i < learnts.size(); i++)
                ca[learnts[i]].activity() *= 1e-20;
            cla_inc *= 1e-20; } }

inline void Solver::checkGarbage(void){ return checkGarbage(garbage_frac); }
inline void Solver::checkGarbage(double gf){
    if (ca.wasted() > ca.size() * gf)
        garbageCollect(); }

// NOTE: enqueue does not set the ok flag! (only public methods do)
inline bool     Solver::enqueue         (Lit p, CRef from)      { return value(p) != l_Undef ? value(p) != l_False : (uncheckedEnqueue(p, from), true); }
inline bool     Solver::addClause       (const vec<Lit>& ps)    { ps.copyTo(add_tmp); return addClause_(add_tmp); }
inline bool     Solver::addEmptyClause  ()                      { add_tmp.clear(); return addClause_(add_tmp); }
inline bool     Solver::addClause       (Lit p)                 { add_tmp.clear(); add_tmp.push(p); return addClause_(add_tmp); }
inline bool     Solver::addClause       (Lit p, Lit q)          { add_tmp.clear(); add_tmp.push(p); add_tmp.push(q); return addClause_(add_tmp); }
inline bool     Solver::addClause       (Lit p, Lit q, Lit r)   { add_tmp.clear(); add_tmp.push(p); add_tmp.push(q); add_tmp.push(r); return addClause_(add_tmp); }
 inline bool     Solver::locked          (const Clause& c) const { 
   if(c.size()>2) 
     return value(c[0]) == l_True && reason(var(c[0])) != CRef_Undef && ca.lea(reason(var(c[0]))) == &c; 
   return 
     (value(c[0]) == l_True && reason(var(c[0])) != CRef_Undef && ca.lea(reason(var(c[0]))) == &c)
     || 
     (value(c[1]) == l_True && reason(var(c[1])) != CRef_Undef && ca.lea(reason(var(c[1]))) == &c);
 }
inline void     Solver::newDecisionLevel()                      { trail_lim.push(trail.size()); }

inline int      Solver::decisionLevel ()      const   { return trail_lim.size(); }
inline uint32_t Solver::abstractLevel (Var x) const   { return 1 << (level(x) & 31); }
inline lbool    Solver::value         (Var x) const   { return assigns[x]; }
inline lbool    Solver::value         (Lit p) const   { return assigns[var(p)] ^ sign(p); }
inline lbool    Solver::modelValue    (Var x) const   { return model[x]; }
inline lbool    Solver::modelValue    (Lit p) const   { return model[var(p)] ^ sign(p); }
inline int      Solver::nAssigns      ()      const   { return trail.size(); }
inline int      Solver::nClauses      ()      const   { return clauses.size(); }
inline int      Solver::nLearnts      ()      const   { return learnts.size(); }
inline int      Solver::nVars         ()      const   { return vardata.size(); }
inline int      Solver::nFreeVars     ()      const   { return (int)dec_vars - (trail_lim.size() == 0 ? trail.size() : trail_lim[0]); }
inline void     Solver::setPolarity   (Var v, bool b) { polarity[v] = b; }
inline void     Solver::setDecisionVar(Var v, bool b) 
{ 
    if      ( b && !decision[v]) dec_vars++;
    else if (!b &&  decision[v]) dec_vars--;

    decision[v] = b;
    insertVarOrder(v);
}
inline void     Solver::setConfBudget(int64_t x){ conflict_budget    = conflicts    + x; }
inline void     Solver::setPropBudget(int64_t x){ propagation_budget = propagations + x; }
inline void     Solver::interrupt(){ asynch_interrupt = true; }
inline void     Solver::clearInterrupt(){ asynch_interrupt = false; }
inline void     Solver::budgetOff(){ conflict_budget = propagation_budget = -1; }
inline bool     Solver::withinBudget() const {
    return !asynch_interrupt &&
           (conflict_budget    < 0 || conflicts < (uint64_t)conflict_budget) &&
           (propagation_budget < 0 || propagations < (uint64_t)propagation_budget); }

// FIXME: after the introduction of asynchronous interrruptions the solve-versions that return a
// pure bool do not give a safe interface. Either interrupts must be possible to turn off here, or
// all calls to solve must return an 'lbool'. I'm not yet sure which I prefer.
inline bool     Solver::solve         ()                    { budgetOff(); assumptions.clear(); return solve_() == l_True; }
inline bool     Solver::solve         (Lit p)               { budgetOff(); assumptions.clear(); assumptions.push(p); return solve_() == l_True; }
inline bool     Solver::solve         (Lit p, Lit q)        { budgetOff(); assumptions.clear(); assumptions.push(p); assumptions.push(q); return solve_() == l_True; }
inline bool     Solver::solve         (Lit p, Lit q, Lit r) { budgetOff(); assumptions.clear(); assumptions.push(p); assumptions.push(q); assumptions.push(r); return solve_() == l_True; }
inline bool     Solver::solve         (const vec<Lit>& assumps){ budgetOff(); assumps.copyTo(assumptions); return solve_() == l_True; }
inline lbool    Solver::solveLimited  (const vec<Lit>& assumps){ assumps.copyTo(assumptions); return solve_(); }
inline bool     Solver::okay          ()      const   { return ok; }

inline void     Solver::toDimacs     (const char* file){ vec<Lit> as; toDimacs(file, as); }
inline void     Solver::toDimacs     (const char* file, Lit p){ vec<Lit> as; as.push(p); toDimacs(file, as); }
inline void     Solver::toDimacs     (const char* file, Lit p, Lit q){ vec<Lit> as; as.push(p); as.push(q); toDimacs(file, as); }
inline void     Solver::toDimacs     (const char* file, Lit p, Lit q, Lit r){ vec<Lit> as; as.push(p); as.push(q); as.push(r); toDimacs(file, as); }


//=================================================================================================
// Debug etc:


inline void Solver::printLit(Lit l)
{
    printf("%s%d:%c", sign(l) ? "-" : "", var(l)+1, value(l) == l_True ? '1' : (value(l) == l_False ? '0' : 'X'));
}


inline void Solver::printClause(CRef cr)
{
  Clause &c = ca[cr];
    for (int i = 0; i < c.size(); i++){
        printLit(c[i]);
        printf(" ");
    }
}

inline void Solver::printInitialClause(CRef cr)
{
  Clause &c = ca[cr];
    for (int i = 0; i < c.size(); i++){
      if(!isSelector(var(c[i]))) {
	printLit(c[i]);
        printf(" ");
      }
    }
}


//=================================================================================================
}

#endif

/****************************************************************************************[System.h]
Copyright (c) 2003-2006, Niklas Een, Niklas Sorensson
Copyright (c) 2007-2010, Niklas Sorensson

Permission is hereby granted, free of charge, to any person obtaining a copy of this software and
associated documentation files (the "Software"), to deal in the Software without restriction,
including without limitation the rights to use, copy, modify, merge, publish, distribute,
sublicense, and/or sell copies of the Software, and to permit persons to whom the Software is
furnished to do so, subject to the following conditions:

The above copyright notice and this permission notice shall be included in all copies or
substantial portions of the Software.

THE SOFTWARE IS PROVIDED "AS IS", WITHOUT WARRANTY OF ANY KIND, EXPRESS OR IMPLIED, INCLUDING BUT
NOT LIMITED TO THE WARRANTIES OF MERCHANTABILITY, FITNESS FOR A PARTICULAR PURPOSE AND
NONINFRINGEMENT. IN NO EVENT SHALL THE AUTHORS OR COPYRIGHT HOLDERS BE LIABLE FOR ANY CLAIM,
DAMAGES OR OTHER LIABILITY, WHETHER IN AN ACTION OF CONTRACT, TORT OR OTHERWISE, ARISING FROM, OUT
OF OR IN CONNECTION WITH THE SOFTWARE OR THE USE OR OTHER DEALINGS IN THE SOFTWARE.
**************************************************************************************************/

#ifndef Glucose_System_h
#define Glucose_System_h

#if defined(__linux__)
#include <fpu_control.h>
#endif


//-------------------------------------------------------------------------------------------------

namespace Glucose {

static inline double cpuTime(void); // CPU-time in seconds.
extern double memUsed();            // Memory in mega bytes (returns 0 for unsupported architectures).
extern double memUsedPeak();        // Peak-memory in mega bytes (returns 0 for unsupported architectures).

}

//-------------------------------------------------------------------------------------------------
// Implementation of inline functions:

#if defined(_MSC_VER) || defined(__MINGW32__)
#include <time.h>

static inline double Glucose::cpuTime(void) { return (double)clock() / CLOCKS_PER_SEC; }

#else
#include <sys/time.h>
#include <sys/resource.h>
#include <unistd.h>

static inline double Glucose::cpuTime(void) {
    struct rusage ru;
    getrusage(RUSAGE_SELF, &ru);
    return (double)ru.ru_utime.tv_sec + (double)ru.ru_utime.tv_usec / 1000000; }

#endif

#endif


using namespace Glucose;

//=================================================================================================
// Options:

static const char* _cat = "CORE";
static const char* _cr = "CORE -- RESTART";
static const char* _cred = "CORE -- REDUCE";
static const char* _cm = "CORE -- MINIMIZE";
static const char* _certified = "CORE -- CERTIFIED UNSAT";




static BoolOption opt_incremental (_cat,"incremental", "Use incremental SAT solving",false);
static DoubleOption opt_K                 (_cr, "K",           "The constant used to force restart",            0.8,     DoubleRange(0, false, 1, false));           
static DoubleOption opt_R                 (_cr, "R",           "The constant used to block restart",            1.4,     DoubleRange(1, false, 5, false));           
static IntOption     opt_size_lbd_queue     (_cr, "szLBDQueue",      "The size of moving average for LBD (restarts)", 50, IntRange(10, INT32_MAX));
static IntOption     opt_size_trail_queue     (_cr, "szTrailQueue",      "The size of moving average for trail (block restarts)", 5000, IntRange(10, INT32_MAX));

static IntOption     opt_first_reduce_db     (_cred, "firstReduceDB",      "The number of conflicts before the first reduce DB", 2000, IntRange(0, INT32_MAX));
static IntOption     opt_inc_reduce_db     (_cred, "incReduceDB",      "Increment for reduce DB", 300, IntRange(0, INT32_MAX));
static IntOption     opt_spec_inc_reduce_db     (_cred, "specialIncReduceDB",      "Special increment for reduce DB", 1000, IntRange(0, INT32_MAX));
static IntOption    opt_lb_lbd_frozen_clause      (_cred, "minLBDFrozenClause",        "Protect clauses if their LBD decrease and is lower than (for one turn)", 30, IntRange(0, INT32_MAX));

static IntOption     opt_lb_size_minimzing_clause     (_cm, "minSizeMinimizingClause",      "The min size required to minimize clause", 30, IntRange(3, INT32_MAX));
static IntOption     opt_lb_lbd_minimzing_clause     (_cm, "minLBDMinimizingClause",      "The min LBD required to minimize clause", 6, IntRange(3, INT32_MAX));


static DoubleOption  opt_var_decay         (_cat, "var-decay",   "The variable activity decay factor",            0.8,     DoubleRange(0, false, 1, false));
static DoubleOption  opt_clause_decay      (_cat, "cla-decay",   "The clause activity decay factor",              0.999,    DoubleRange(0, false, 1, false));
static DoubleOption  opt_random_var_freq   (_cat, "rnd-freq",    "The frequency with which the decision heuristic tries to choose a random variable", 0, DoubleRange(0, true, 1, true));
static DoubleOption  opt_random_seed       (_cat, "rnd-seed",    "Used by the random variable selection",         91648253, DoubleRange(0, false, HUGE_VAL, false));
static IntOption     opt_ccmin_mode        (_cat, "ccmin-mode",  "Controls conflict clause minimization (0=none, 1=basic, 2=deep)", 2, IntRange(0, 2));
static IntOption     opt_phase_saving      (_cat, "phase-saving", "Controls the level of phase saving (0=none, 1=limited, 2=full)", 2, IntRange(0, 2));
static BoolOption    opt_rnd_init_act      (_cat, "rnd-init",    "Randomize the initial activity", false);
/*
static IntOption     opt_restart_first     (_cat, "rfirst",      "The base restart interval", 100, IntRange(1, INT32_MAX));
static DoubleOption  opt_restart_inc       (_cat, "rinc",        "Restart interval increase factor", 2, DoubleRange(1, false, HUGE_VAL, false));
*/
static DoubleOption  opt_garbage_frac      (_cat, "gc-frac",     "The fraction of wasted memory allowed before a garbage collection is triggered",  0.20, DoubleRange(0, false, HUGE_VAL, false));


 BoolOption    opt_certified      (_certified, "certified",    "Certified UNSAT using DRUP format", false);
 StringOption    opt_certified_file      (_certified, "certified-output",    "Certified UNSAT output file", "NULL");


//=================================================================================================
// Constructor/Destructor:


Solver::Solver() :

    // Parameters (user settable):
    //
    verbosity        (0)
    , showModel        (0)
    , K              (opt_K)
    , R              (opt_R)
    , sizeLBDQueue   (opt_size_lbd_queue)
    , sizeTrailQueue   (opt_size_trail_queue)
    , firstReduceDB  (opt_first_reduce_db)
    , incReduceDB    (opt_inc_reduce_db)
    , specialIncReduceDB    (opt_spec_inc_reduce_db)
    , lbLBDFrozenClause (opt_lb_lbd_frozen_clause)
    , lbSizeMinimizingClause (opt_lb_size_minimzing_clause)
    , lbLBDMinimizingClause (opt_lb_lbd_minimzing_clause)
  , var_decay        (opt_var_decay)
  , clause_decay     (opt_clause_decay)
  , random_var_freq  (opt_random_var_freq)
  , random_seed      (opt_random_seed)
  , ccmin_mode       (opt_ccmin_mode)
  , phase_saving     (opt_phase_saving)
  , rnd_pol          (false)
  , rnd_init_act     (opt_rnd_init_act)
  , garbage_frac     (opt_garbage_frac)
  , certifiedOutput  (NULL)
  , certifiedUNSAT   (opt_certified) 
    // Statistics: (formerly in 'SolverStats')
    //
  ,  nbRemovedClauses(0),nbReducedClauses(0), nbDL2(0),nbBin(0),nbUn(0) , nbReduceDB(0)
    , solves(0), starts(0), decisions(0), rnd_decisions(0), propagations(0),conflicts(0),conflictsRestarts(0),nbstopsrestarts(0),nbstopsrestartssame(0),lastblockatrestart(0)
  , dec_vars(0), clauses_literals(0), learnts_literals(0), max_literals(0), tot_literals(0)
    , curRestart(1)

  , ok                 (true)
  , cla_inc            (1)
  , var_inc            (1)
  , watches            (WatcherDeleted(ca))
  , watchesBin            (WatcherDeleted(ca))
  , qhead              (0)
  , simpDB_assigns     (-1)
  , simpDB_props       (0)
  , order_heap         (VarOrderLt(activity))
  , progress_estimate  (0)
  , remove_satisfied   (true)

    // Resource constraints:
    //
  , conflict_budget    (-1)
  , propagation_budget (-1)
  , asynch_interrupt   (false)
  , incremental(opt_incremental)
  , nbVarsInitialFormula(INT32_MAX)
{
  MYFLAG=0;  
  // Initialize only first time. Useful for incremental solving, useless otherwise
  lbdQueue.initSize(sizeLBDQueue);
  trailQueue.initSize(sizeTrailQueue);
  sumLBD = 0;
  nbclausesbeforereduce = firstReduceDB;
  totalTime4Sat=0;totalTime4Unsat=0;
  nbSatCalls=0;nbUnsatCalls=0;


  if(certifiedUNSAT) {
    if(!strcmp(opt_certified_file,"NULL")) {
      certifiedOutput =  fopen("/dev/stdout", "wb");
    } else {
      certifiedOutput =  fopen(opt_certified_file, "wb");	    
    }
    //    fprintf(certifiedOutput,"o proof DRUP\n");
  }
}


Solver::~Solver()
{
}


/****************************************************************
 Set the incremental mode
****************************************************************/

// This function set the incremental mode to true.
// You can add special code for this mode here.

void Solver::setIncrementalMode() {
  incremental = true;
}

// Number of variables without selectors
void Solver::initNbInitialVars(int nb) {
  nbVarsInitialFormula = nb;
}


//=================================================================================================
// Minor methods:


// Creates a new SAT variable in the solver. If 'decision' is cleared, variable will not be
// used as a decision variable (NOTE! This has effects on the meaning of a SATISFIABLE result).
//
Var Solver::newVar(bool sign, bool dvar)
{
    int v = nVars();
    watches  .init(mkLit(v, false));
    watches  .init(mkLit(v, true ));
    watchesBin  .init(mkLit(v, false));
    watchesBin  .init(mkLit(v, true ));
    assigns  .push(l_Undef);
    vardata  .push(mkVarData(CRef_Undef, 0));
    //activity .push(0);
    activity .push(rnd_init_act ? drand(random_seed) * 0.00001 : 0);
    seen     .push(0);
    permDiff  .push(0);
    polarity .push(sign);
    decision .push();
    trail    .capacity(v+1);
    setDecisionVar(v, dvar);
    return v;
}



bool Solver::addClause_(vec<Lit>& ps)
{
    assert(decisionLevel() == 0);
    if (!ok) return false;

    // Check if clause is satisfied and remove false/duplicate literals:
    sort(ps);

    vec<Lit>    oc;
    oc.clear();

    Lit p; int i, j, flag = 0;
    if(certifiedUNSAT) {
      for (i = j = 0, p = lit_Undef; i < ps.size(); i++) {
        oc.push(ps[i]);
        if (value(ps[i]) == l_True || ps[i] == ~p || value(ps[i]) == l_False)
          flag = 1;
      }
    }

    for (i = j = 0, p = lit_Undef; i < ps.size(); i++)
      if (value(ps[i]) == l_True || ps[i] == ~p)
	return true;
      else if (value(ps[i]) != l_False && ps[i] != p)
	ps[j++] = p = ps[i];
    ps.shrink(i - j);
    
    if (flag && (certifiedUNSAT)) {
      for (i = j = 0, p = lit_Undef; i < ps.size(); i++)
        fprintf(certifiedOutput, "%i ", (var(ps[i]) + 1) * (-2 * sign(ps[i]) + 1));
      fprintf(certifiedOutput, "0\n");

      fprintf(certifiedOutput, "d ");
      for (i = j = 0, p = lit_Undef; i < oc.size(); i++)
        fprintf(certifiedOutput, "%i ", (var(oc[i]) + 1) * (-2 * sign(oc[i]) + 1));
      fprintf(certifiedOutput, "0\n");
    }

    if (ps.size() == 0)
        return ok = false;
    else if (ps.size() == 1){
        uncheckedEnqueue(ps[0]);
        return ok = (propagate() == CRef_Undef);
    }else{
        CRef cr = ca.alloc(ps, false);
        clauses.push(cr);
        attachClause(cr);
    }

    return true;
}


void Solver::attachClause(CRef cr) {
    const Clause& c = ca[cr];

    assert(c.size() > 1);
    if(c.size()==2) {
      watchesBin[~c[0]].push(Watcher(cr, c[1]));
      watchesBin[~c[1]].push(Watcher(cr, c[0]));
    } else {
      watches[~c[0]].push(Watcher(cr, c[1]));
      watches[~c[1]].push(Watcher(cr, c[0]));
    }
    if (c.learnt()) learnts_literals += c.size();
    else            clauses_literals += c.size(); }




void Solver::detachClause(CRef cr, bool strict) {
    const Clause& c = ca[cr];
    
    assert(c.size() > 1);
    if(c.size()==2) {
      if (strict){
        remove(watchesBin[~c[0]], Watcher(cr, c[1]));
        remove(watchesBin[~c[1]], Watcher(cr, c[0]));
      }else{
        // Lazy detaching: (NOTE! Must clean all watcher lists before garbage collecting this clause)
        watchesBin.smudge(~c[0]);
        watchesBin.smudge(~c[1]);
      }
    } else {
      if (strict){
        remove(watches[~c[0]], Watcher(cr, c[1]));
        remove(watches[~c[1]], Watcher(cr, c[0]));
      }else{
        // Lazy detaching: (NOTE! Must clean all watcher lists before garbage collecting this clause)
        watches.smudge(~c[0]);
        watches.smudge(~c[1]);
      }
    }
    if (c.learnt()) learnts_literals -= c.size();
    else            clauses_literals -= c.size(); }


void Solver::removeClause(CRef cr) {

  Clause& c = ca[cr];

  if (certifiedUNSAT) {
    fprintf(certifiedOutput, "d ");
    for (int i = 0; i < c.size(); i++)
      fprintf(certifiedOutput, "%i ", (var(c[i]) + 1) * (-2 * sign(c[i]) + 1));
    fprintf(certifiedOutput, "0\n");
  }

  detachClause(cr);
  // Don't leave pointers to free'd memory!
  if (locked(c)) vardata[var(c[0])].reason = CRef_Undef;
  c.mark(1);
  ca.free(cr);
}


bool Solver::satisfied(const Clause& c) const {
  if(incremental)  // Check clauses with many selectors is too time consuming
    return (value(c[0]) == l_True) || (value(c[1]) == l_True);

  // Default mode.
    for (int i = 0; i < c.size(); i++)
        if (value(c[i]) == l_True)
            return true;
    return false; 
}

/************************************************************
 * Compute LBD functions
 *************************************************************/

inline unsigned int Solver::computeLBD(const vec<Lit> & lits,int end) {
  int nblevels = 0;
  MYFLAG++;

  if(incremental) { // ----------------- INCREMENTAL MODE
    if(end==-1) end = lits.size();
    unsigned int nbDone = 0;
    for(int i=0;i<lits.size();i++) {
      if(nbDone>=end) break;
      if(isSelector(var(lits[i]))) continue;
      nbDone++;
      int l = level(var(lits[i]));
      if (permDiff[l] != MYFLAG) {
	permDiff[l] = MYFLAG;
	nblevels++;
      }
    }
  } else { // -------- DEFAULT MODE. NOT A LOT OF DIFFERENCES... BUT EASIER TO READ
    for(int i=0;i<lits.size();i++) {
      int l = level(var(lits[i]));
      if (permDiff[l] != MYFLAG) {
	permDiff[l] = MYFLAG;
	nblevels++;
      }
    }
  }

  return nblevels;
}

inline unsigned int Solver::computeLBD(const Clause &c) {
  int nblevels = 0;
  MYFLAG++;

  if(incremental) { // ----------------- INCREMENTAL MODE
     int nbDone = 0;
    for(int i=0;i<c.size();i++) {
      if(nbDone>=c.sizeWithoutSelectors()) break;
      if(isSelector(var(c[i]))) continue;
      nbDone++;
      int l = level(var(c[i]));
      if (permDiff[l] != MYFLAG) {
	permDiff[l] = MYFLAG;
	nblevels++;
      }
    }
  } else { // -------- DEFAULT MODE. NOT A LOT OF DIFFERENCES... BUT EASIER TO READ
    for(int i=0;i<c.size();i++) {
      int l = level(var(c[i]));
      if (permDiff[l] != MYFLAG) {
	permDiff[l] = MYFLAG;
	nblevels++;
      }
    }
  }
  return nblevels;
}


/******************************************************************
 * Minimisation with binary reolution
 ******************************************************************/
void Solver::minimisationWithBinaryResolution(vec<Lit> &out_learnt) {
  
  // Find the LBD measure                                                                                                         
  unsigned int lbd = computeLBD(out_learnt);
  Lit p = ~out_learnt[0];
  
  if(lbd<=lbLBDMinimizingClause){
    MYFLAG++;
    
    for(int i = 1;i<out_learnt.size();i++) {
      permDiff[var(out_learnt[i])] = MYFLAG;
    }

    vec<Watcher>&  wbin  = watchesBin[p];
    int nb = 0;
    for(int k = 0;k<wbin.size();k++) {
      Lit imp = wbin[k].blocker;
      if(permDiff[var(imp)]==MYFLAG && value(imp)==l_True) {
	nb++;
	permDiff[var(imp)]= MYFLAG-1;
      }
      }
    int l = out_learnt.size()-1;
    if(nb>0) {
      nbReducedClauses++;
      for(int i = 1;i<out_learnt.size()-nb;i++) {
	if(permDiff[var(out_learnt[i])]!=MYFLAG) {
	  Lit p = out_learnt[l];
	  out_learnt[l] = out_learnt[i];
	  out_learnt[i] = p;
	  l--;i--;
	}
      }
      
      out_learnt.shrink(nb);
      
    }
  }
}

// Revert to the state at given level (keeping all assignment at 'level' but not beyond).
//
void Solver::cancelUntil(int level) {
    if (decisionLevel() > level){
        for (int c = trail.size()-1; c >= trail_lim[level]; c--){
            Var      x  = var(trail[c]);
            assigns [x] = l_Undef;
	    if (phase_saving > 1 || ((phase_saving == 1) && c > trail_lim.last()))
                polarity[x] = sign(trail[c]);
            insertVarOrder(x); }
        qhead = trail_lim[level];
        trail.shrink(trail.size() - trail_lim[level]);
        trail_lim.shrink(trail_lim.size() - level);
    } 
}


//=================================================================================================
// Major methods:


Lit Solver::pickBranchLit()
{
    Var next = var_Undef;

    // Random decision:
    if (drand(random_seed) < random_var_freq && !order_heap.empty()){
        next = order_heap[irand(random_seed,order_heap.size())];
        if (value(next) == l_Undef && decision[next])
            rnd_decisions++; }

    // Activity based decision:
    while (next == var_Undef || value(next) != l_Undef || !decision[next])
        if (order_heap.empty()){
            next = var_Undef;
            break;
        }else
            next = order_heap.removeMin();

    return next == var_Undef ? lit_Undef : mkLit(next, rnd_pol ? drand(random_seed) < 0.5 : polarity[next]);
}


/*_________________________________________________________________________________________________
|
|  analyze : (confl : Clause*) (out_learnt : vec<Lit>&) (out_btlevel : int&)  ->  [void]
|  
|  Description:
|    Analyze conflict and produce a reason clause.
|  
|    Pre-conditions:
|      * 'out_learnt' is assumed to be cleared.
|      * Current decision level must be greater than root level.
|  
|    Post-conditions:
|      * 'out_learnt[0]' is the asserting literal at level 'out_btlevel'.
|      * If out_learnt.size() > 1 then 'out_learnt[1]' has the greatest decision level of the 
|        rest of literals. There may be others from the same level though.
|  
|________________________________________________________________________________________________@*/
void Solver::analyze(CRef confl, vec<Lit>& out_learnt,vec<Lit>&selectors, int& out_btlevel,unsigned int &lbd,unsigned int &szWithoutSelectors)
{
    int pathC = 0;
    Lit p     = lit_Undef;

    // Generate conflict clause:
    //
    out_learnt.push();      // (leave room for the asserting literal)
    int index   = trail.size() - 1;

    do{
        assert(confl != CRef_Undef); // (otherwise should be UIP)
        Clause& c = ca[confl];

	// Special case for binary clauses
	// The first one has to be SAT
	if( p != lit_Undef && c.size()==2 && value(c[0])==l_False) {
	  
	  assert(value(c[1])==l_True);
	  Lit tmp = c[0];
	  c[0] =  c[1], c[1] = tmp;
	}
	
	if (c.learnt()) 
            claBumpActivity(c);

#ifdef DYNAMICNBLEVEL		    
	// DYNAMIC NBLEVEL trick (see competition'09 companion paper)
	if(c.learnt()  && c.lbd()>2) { 
	  unsigned int nblevels = computeLBD(c);
	  if(nblevels+1<c.lbd() ) { // improve the LBD
	    if(c.lbd()<=lbLBDFrozenClause) {
	      c.setCanBeDel(false); 
	    }
	    // seems to be interesting : keep it for the next round
	    c.setLBD(nblevels); // Update it
	  }
	}
#endif


        for (int j = (p == lit_Undef) ? 0 : 1; j < c.size(); j++){
            Lit q = c[j];

            if (!seen[var(q)] && level(var(q)) > 0){
	      if(!isSelector(var(q)))
                varBumpActivity(var(q));
	      seen[var(q)] = 1;
	      if (level(var(q)) >= decisionLevel()) {
		pathC++;
#ifdef UPDATEVARACTIVITY
		// UPDATEVARACTIVITY trick (see competition'09 companion paper)
		if(!isSelector(var(q)) && (reason(var(q))!= CRef_Undef)  && ca[reason(var(q))].learnt()) 
		  lastDecisionLevel.push(q);
#endif
		
	      } else {
		if(isSelector(var(q))) {
		  assert(value(q) == l_False);
		  selectors.push(q);
		} else 
		  out_learnt.push(q);
	      }
	    }
        }
        
        // Select next clause to look at:
        while (!seen[var(trail[index--])]);
        p     = trail[index+1];
        confl = reason(var(p));
        seen[var(p)] = 0;
        pathC--;

    }while (pathC > 0);
    out_learnt[0] = ~p;

    // Simplify conflict clause:
    //
    int i, j;

    for(int i = 0;i<selectors.size();i++)  
      out_learnt.push(selectors[i]);       

    out_learnt.copyTo(analyze_toclear);
    if (ccmin_mode == 2){
        uint32_t abstract_level = 0;
        for (i = 1; i < out_learnt.size(); i++)
            abstract_level |= abstractLevel(var(out_learnt[i])); // (maintain an abstraction of levels involved in conflict)

        for (i = j = 1; i < out_learnt.size(); i++)
            if (reason(var(out_learnt[i])) == CRef_Undef || !litRedundant(out_learnt[i], abstract_level))
                out_learnt[j++] = out_learnt[i];
        
    }else if (ccmin_mode == 1){
        for (i = j = 1; i < out_learnt.size(); i++){
            Var x = var(out_learnt[i]);

            if (reason(x) == CRef_Undef)
                out_learnt[j++] = out_learnt[i];
            else{
                Clause& c = ca[reason(var(out_learnt[i]))];
                // Thanks to Siert Wieringa for this bug fix!
                for (int k = ((c.size()==2) ? 0:1); k < c.size(); k++)
                    if (!seen[var(c[k])] && level(var(c[k])) > 0){
                        out_learnt[j++] = out_learnt[i];
                        break; }
            }
        }
    }else
        i = j = out_learnt.size();

    max_literals += out_learnt.size();
    out_learnt.shrink(i - j);
    tot_literals += out_learnt.size();


    /* ***************************************
      Minimisation with binary clauses of the asserting clause
      First of all : we look for small clauses
      Then, we reduce clauses with small LBD.
      Otherwise, this can be useless
     */
    if(!incremental && out_learnt.size()<=lbSizeMinimizingClause) {
      minimisationWithBinaryResolution(out_learnt);
    }
    // Find correct backtrack level:
    //
    if (out_learnt.size() == 1)
        out_btlevel = 0;
    else{
        int max_i = 1;
        // Find the first literal assigned at the next-highest level:
        for (int i = 2; i < out_learnt.size(); i++)
            if (level(var(out_learnt[i])) > level(var(out_learnt[max_i])))
                max_i = i;
        // Swap-in this literal at index 1:
        Lit p             = out_learnt[max_i];
        out_learnt[max_i] = out_learnt[1];
        out_learnt[1]     = p;
        out_btlevel       = level(var(p));
    }


    // Compute the size of the clause without selectors (incremental mode)
    if(incremental) {
      szWithoutSelectors = 0;
      for(int i=0;i<out_learnt.size();i++) {
	if(!isSelector(var((out_learnt[i])))) szWithoutSelectors++; 
	else if(i>0) break;
      }
    } else 
      szWithoutSelectors = out_learnt.size();
    
    // Compute LBD
    lbd = computeLBD(out_learnt,out_learnt.size()-selectors.size());

  
#ifdef UPDATEVARACTIVITY
  // UPDATEVARACTIVITY trick (see competition'09 companion paper)
  if(lastDecisionLevel.size()>0) {
    for(int i = 0;i<lastDecisionLevel.size();i++) {
      if(ca[reason(var(lastDecisionLevel[i]))].lbd()<lbd)
	varBumpActivity(var(lastDecisionLevel[i]));
    }
    lastDecisionLevel.clear();
  } 
#endif	    



  for (int j = 0; j < analyze_toclear.size(); j++) seen[var(analyze_toclear[j])] = 0;    // ('seen[]' is now cleared)
  for(int j = 0 ; j<selectors.size() ; j++) seen[var(selectors[j])] = 0;  
}


// Check if 'p' can be removed. 'abstract_levels' is used to abort early if the algorithm is
// visiting literals at levels that cannot be removed later.
bool Solver::litRedundant(Lit p, uint32_t abstract_levels)
{
    analyze_stack.clear(); analyze_stack.push(p);
    int top = analyze_toclear.size();
    while (analyze_stack.size() > 0){
        assert(reason(var(analyze_stack.last())) != CRef_Undef);
        Clause& c = ca[reason(var(analyze_stack.last()))]; analyze_stack.pop();
	if(c.size()==2 && value(c[0])==l_False) {
	  assert(value(c[1])==l_True);
	  Lit tmp = c[0];
	  c[0] =  c[1], c[1] = tmp;
	}

        for (int i = 1; i < c.size(); i++){
            Lit p  = c[i];
            if (!seen[var(p)] && level(var(p)) > 0){
                if (reason(var(p)) != CRef_Undef && (abstractLevel(var(p)) & abstract_levels) != 0){
                    seen[var(p)] = 1;
                    analyze_stack.push(p);
                    analyze_toclear.push(p);
                }else{
                    for (int j = top; j < analyze_toclear.size(); j++)
                        seen[var(analyze_toclear[j])] = 0;
                    analyze_toclear.shrink(analyze_toclear.size() - top);
                    return false;
                }
            }
        }
    }

    return true;
}


/*_________________________________________________________________________________________________
|
|  analyzeFinal : (p : Lit)  ->  [void]
|  
|  Description:
|    Specialized analysis procedure to express the final conflict in terms of assumptions.
|    Calculates the (possibly empty) set of assumptions that led to the assignment of 'p', and
|    stores the result in 'out_conflict'.
|________________________________________________________________________________________________@*/
void Solver::analyzeFinal(Lit p, vec<Lit>& out_conflict)
{
    out_conflict.clear();
    out_conflict.push(p);

    if (decisionLevel() == 0)
        return;

    seen[var(p)] = 1;

    for (int i = trail.size()-1; i >= trail_lim[0]; i--){
        Var x = var(trail[i]);
        if (seen[x]){
            if (reason(x) == CRef_Undef){
                assert(level(x) > 0);
                out_conflict.push(~trail[i]);
            }else{
                Clause& c = ca[reason(x)];
		//                for (int j = 1; j < c.size(); j++) Minisat (glucose 2.0) loop 
		// Bug in case of assumptions due to special data structures for Binary.
		// Many thanks to Sam Bayless (sbayless@cs.ubc.ca) for discover this bug.
		for (int j = ((c.size()==2) ? 0:1); j < c.size(); j++)
                    if (level(var(c[j])) > 0)
                        seen[var(c[j])] = 1;
            }  

            seen[x] = 0;
        }
    }

    seen[var(p)] = 0;
}


void Solver::uncheckedEnqueue(Lit p, CRef from)
{
    assert(value(p) == l_Undef);
    assigns[var(p)] = lbool(!sign(p));
    vardata[var(p)] = mkVarData(from, decisionLevel());
    trail.push_(p);
}


/*_________________________________________________________________________________________________
|
|  propagate : [void]  ->  [Clause*]
|  
|  Description:
|    Propagates all enqueued facts. If a conflict arises, the conflicting clause is returned,
|    otherwise CRef_Undef.
|  
|    Post-conditions:
|      * the propagation queue is empty, even if there was a conflict.
|________________________________________________________________________________________________@*/
CRef Solver::propagate()
{
    CRef    confl     = CRef_Undef;
    int     num_props = 0;
    watches.cleanAll();
    watchesBin.cleanAll();
    while (qhead < trail.size()){
        Lit            p   = trail[qhead++];     // 'p' is enqueued fact to propagate.
        vec<Watcher>&  ws  = watches[p];
        Watcher        *i, *j, *end;
        num_props++;

	
	    // First, Propagate binary clauses 
	vec<Watcher>&  wbin  = watchesBin[p];
	
	for(int k = 0;k<wbin.size();k++) {
	  
	  Lit imp = wbin[k].blocker;
	  
	  if(value(imp) == l_False) {
	    return wbin[k].cref;
	  }
	  
	  if(value(imp) == l_Undef) {
	    uncheckedEnqueue(imp,wbin[k].cref);
	  }
	}
    


        for (i = j = (Watcher*)ws, end = i + ws.size();  i != end;){
            // Try to avoid inspecting the clause:
            Lit blocker = i->blocker;
            if (value(blocker) == l_True){
                *j++ = *i++; continue; }

            // Make sure the false literal is data[1]:
            CRef     cr        = i->cref;
            Clause&  c         = ca[cr];
            Lit      false_lit = ~p;
            if (c[0] == false_lit)
                c[0] = c[1], c[1] = false_lit;
            assert(c[1] == false_lit);
            i++;

            // If 0th watch is true, then clause is already satisfied.
            Lit     first = c[0];
            Watcher w     = Watcher(cr, first);
            if (first != blocker && value(first) == l_True){
	      
	      *j++ = w; continue; }

            // Look for new watch:
	    if(incremental) { // ----------------- INCREMENTAL MODE
	      int choosenPos = -1;
	      for (int k = 2; k < c.size(); k++) {
		
		if (value(c[k]) != l_False){
		  if(decisionLevel()>assumptions.size()) {
		    choosenPos = k;
		    break;
		  } else {
		    choosenPos = k;
		    
		    if(value(c[k])==l_True || !isSelector(var(c[k]))) {
		      break;
		    }
		  }

		}
	      }
	      if(choosenPos!=-1) {
		c[1] = c[choosenPos]; c[choosenPos] = false_lit;
		watches[~c[1]].push(w);
		goto NextClause; }
	    } else {  // ----------------- DEFAULT  MODE (NOT INCREMENTAL)
	      for (int k = 2; k < c.size(); k++) {
		
		if (value(c[k]) != l_False){
		  c[1] = c[k]; c[k] = false_lit;
		  watches[~c[1]].push(w);
		  goto NextClause; }
	      }
	    }

            // Did not find watch -- clause is unit under assignment:
            *j++ = w;
            if (value(first) == l_False){
                confl = cr;
                qhead = trail.size();
                // Copy the remaining watches:
                while (i < end)
                    *j++ = *i++;
            }else {
                uncheckedEnqueue(first, cr);
	  
		
	    }
        NextClause:;
        }
        ws.shrink(i - j);
    }
    propagations += num_props;
    simpDB_props -= num_props;
    
    return confl;
}


/*_________________________________________________________________________________________________
|
|  reduceDB : ()  ->  [void]
|  
|  Description:
|    Remove half of the learnt clauses, minus the clauses locked by the current assignment. Locked
|    clauses are clauses that are reason to some assignment. Binary clauses are never removed.
|________________________________________________________________________________________________@*/
struct reduceDB_lt { 
    ClauseAllocator& ca;
    reduceDB_lt(ClauseAllocator& ca_) : ca(ca_) {}
    bool operator () (CRef x, CRef y) { 
 
    // Main criteria... Like in MiniSat we keep all binary clauses
    if(ca[x].size()> 2 && ca[y].size()==2) return 1;
    
    if(ca[y].size()>2 && ca[x].size()==2) return 0;
    if(ca[x].size()==2 && ca[y].size()==2) return 0;
    
    // Second one  based on literal block distance
    if(ca[x].lbd()> ca[y].lbd()) return 1;
    if(ca[x].lbd()< ca[y].lbd()) return 0;    
    
    
    // Finally we can use old activity or size, we choose the last one
        return ca[x].activity() < ca[y].activity();
	//return x->size() < y->size();

        //return ca[x].size() > 2 && (ca[y].size() == 2 || ca[x].activity() < ca[y].activity()); } 
    }    
};

void Solver::reduceDB()
{
 
  int     i, j;
  nbReduceDB++;
  sort(learnts, reduceDB_lt(ca));

  // We have a lot of "good" clauses, it is difficult to compare them. Keep more !
  if(ca[learnts[learnts.size() / RATIOREMOVECLAUSES]].lbd()<=3) nbclausesbeforereduce +=specialIncReduceDB; 
  // Useless :-)
  if(ca[learnts.last()].lbd()<=5)  nbclausesbeforereduce +=specialIncReduceDB; 
  
  
  // Don't delete binary or locked clauses. From the rest, delete clauses from the first half
  // Keep clauses which seem to be usefull (their lbd was reduce during this sequence)

  int limit = learnts.size() / 2;

  for (i = j = 0; i < learnts.size(); i++){
    Clause& c = ca[learnts[i]];
    if (c.lbd()>2 && c.size() > 2 && c.canBeDel() &&  !locked(c) && (i < limit)) {
      removeClause(learnts[i]);
      nbRemovedClauses++;
    }
    else {
      if(!c.canBeDel()) limit++; //we keep c, so we can delete an other clause
      c.setCanBeDel(true);       // At the next step, c can be delete
      learnts[j++] = learnts[i];
    }
  }
  learnts.shrink(i - j);
  checkGarbage();
}


void Solver::removeSatisfied(vec<CRef>& cs)
{
  
    int i, j;
    for (i = j = 0; i < cs.size(); i++){
        Clause& c = ca[cs[i]];


        if (satisfied(c)) 
            removeClause(cs[i]);
        else
            cs[j++] = cs[i];
    }
    cs.shrink(i - j);
}


void Solver::rebuildOrderHeap()
{
    vec<Var> vs;
    for (Var v = 0; v < nVars(); v++)
        if (decision[v] && value(v) == l_Undef)
            vs.push(v);
    order_heap.build(vs);
}


/*_________________________________________________________________________________________________
|
|  simplify : [void]  ->  [bool]
|  
|  Description:
|    Simplify the clause database according to the current top-level assigment. Currently, the only
|    thing done here is the removal of satisfied clauses, but more things can be put here.
|________________________________________________________________________________________________@*/
bool Solver::simplify()
{
    assert(decisionLevel() == 0);

    if (!ok || propagate() != CRef_Undef)
        return ok = false;

    if (nAssigns() == simpDB_assigns || (simpDB_props > 0))
        return true;

    // Remove satisfied clauses:
    removeSatisfied(learnts);
    if (remove_satisfied)        // Can be turned off.
        removeSatisfied(clauses);
    checkGarbage();
    rebuildOrderHeap();

    simpDB_assigns = nAssigns();
    simpDB_props   = clauses_literals + learnts_literals;   // (shouldn't depend on stats really, but it will do for now)

    return true;
}


/*_________________________________________________________________________________________________
|
|  search : (nof_conflicts : int) (params : const SearchParams&)  ->  [lbool]
|  
|  Description:
|    Search for a model the specified number of conflicts. 
|    NOTE! Use negative value for 'nof_conflicts' indicate infinity.
|  
|  Output:
|    'l_True' if a partial assigment that is consistent with respect to the clauseset is found. If
|    all variables are decision variables, this means that the clause set is satisfiable. 'l_False'
|    if the clause set is unsatisfiable. 'l_Undef' if the bound on number of conflicts is reached.
|________________________________________________________________________________________________@*/
lbool Solver::search(int nof_conflicts)
{
    assert(ok);
    int         backtrack_level;
    int         conflictC = 0;
    vec<Lit>    learnt_clause,selectors;
    unsigned int nblevels,szWoutSelectors;
    bool blocked=false;
    starts++;
    for (;;){
        CRef confl = propagate();
        if (confl != CRef_Undef){
            // CONFLICT
	  conflicts++; conflictC++;conflictsRestarts++;
	  if(conflicts%5000==0 && var_decay<0.95)
            var_decay += 0.01;

	  if (verbosity >= 1 && conflicts%verbEveryConflicts==0){
	    printf("c | %8d   %7d    %5d | %7d %8d %8d | %5d %8d   %6d %8d | %6.3f %% |\n", 
		   (int)starts,(int)nbstopsrestarts, (int)(conflicts/starts), 
		   (int)dec_vars - (trail_lim.size() == 0 ? trail.size() : trail_lim[0]), nClauses(), (int)clauses_literals, 
		   (int)nbReduceDB, nLearnts(), (int)nbDL2,(int)nbRemovedClauses, progressEstimate()*100);
	  }
	  if (decisionLevel() == 0) {
	    return l_False;
	    
	  }
	  
	  trailQueue.push(trail.size());
	  // BLOCK RESTART (CP 2012 paper)
	  if( conflictsRestarts>LOWER_BOUND_FOR_BLOCKING_RESTART && lbdQueue.isvalid()  && trail.size()>R*trailQueue.getavg()) {
	    lbdQueue.fastclear();
	    nbstopsrestarts++;
	    if(!blocked) {lastblockatrestart=starts;nbstopsrestartssame++;blocked=true;}
	  }

            learnt_clause.clear();
	    selectors.clear();
            analyze(confl, learnt_clause, selectors,backtrack_level,nblevels,szWoutSelectors);

	    lbdQueue.push(nblevels);
	    sumLBD += nblevels;
 

            cancelUntil(backtrack_level);

            if (certifiedUNSAT) {
              for (int i = 0; i < learnt_clause.size(); i++)
                fprintf(certifiedOutput, "%i " , (var(learnt_clause[i]) + 1) *
                            (-2 * sign(learnt_clause[i]) + 1) );
              fprintf(certifiedOutput, "0\n");
            }

            if (learnt_clause.size() == 1){
	      uncheckedEnqueue(learnt_clause[0]);nbUn++;
            }else{
                CRef cr = ca.alloc(learnt_clause, true);
		ca[cr].setLBD(nblevels); 
		ca[cr].setSizeWithoutSelectors(szWoutSelectors);
		if(nblevels<=2) nbDL2++; // stats
		if(ca[cr].size()==2) nbBin++; // stats
                learnts.push(cr);
                attachClause(cr);

                claBumpActivity(ca[cr]);
                uncheckedEnqueue(learnt_clause[0], cr);
            }
            varDecayActivity();
            claDecayActivity();

           
        }else{
	  // Our dynamic restart, see the SAT09 competition compagnion paper 
	  if (
	      ( lbdQueue.isvalid() && ((lbdQueue.getavg()*K) > (sumLBD / conflictsRestarts)))) {
	    lbdQueue.fastclear();
	    progress_estimate = progressEstimate();
	    int bt = 0;
	    if(incremental) { // DO NOT BACKTRACK UNTIL 0.. USELESS
	      bt = (decisionLevel()<assumptions.size()) ? decisionLevel() : assumptions.size();
	    }
	    cancelUntil(bt);
	    return l_Undef; }


           // Simplify the set of problem clauses:
	  if (decisionLevel() == 0 && !simplify()) {
	    return l_False;
	  }
	    // Perform clause database reduction !
	    if(conflicts>=curRestart* nbclausesbeforereduce) 
	      {
	
		assert(learnts.size()>0);
		curRestart = (conflicts/ nbclausesbeforereduce)+1;
		reduceDB();
		nbclausesbeforereduce += incReduceDB;
	      }
	    
            Lit next = lit_Undef;
            while (decisionLevel() < assumptions.size()){
                // Perform user provided assumption:
                Lit p = assumptions[decisionLevel()];
                if (value(p) == l_True){
                    // Dummy decision level:
                    newDecisionLevel();
                }else if (value(p) == l_False){
                    analyzeFinal(~p, conflict);
                    return l_False;
                }else{
                    next = p;
                    break;
                }
            }

            if (next == lit_Undef){
                // New variable decision:
                decisions++;
                next = pickBranchLit();

                if (next == lit_Undef){
		  printf("c last restart ## conflicts  :  %d %d \n",conflictC,decisionLevel());
		  // Model found:
		  return l_True;
		}
            }

            // Increase decision level and enqueue 'next'
            newDecisionLevel();
            uncheckedEnqueue(next);
        }
    }
}


double Solver::progressEstimate() const
{
    double  progress = 0;
    double  F = 1.0 / nVars();

    for (int i = 0; i <= decisionLevel(); i++){
        int beg = i == 0 ? 0 : trail_lim[i - 1];
        int end = i == decisionLevel() ? trail.size() : trail_lim[i];
        progress += pow(F, i) * (end - beg);
    }

    return progress / nVars();
}

void Solver::printIncrementalStats() {

  printf("c---------- Glucose Stats -------------------------\n");
  printf("c restarts              : %lld\n", starts);
  printf("c nb ReduceDB           : %lld\n", nbReduceDB);
  printf("c nb removed Clauses    : %lld\n",nbRemovedClauses);
  printf("c nb learnts DL2        : %lld\n", nbDL2);
  printf("c nb learnts size 2     : %lld\n", nbBin);
  printf("c nb learnts size 1     : %lld\n", nbUn);

  printf("c conflicts             : %lld \n",conflicts);
  printf("c decisions             : %lld\n",decisions);
  printf("c propagations          : %lld\n",propagations);

  printf("c SAT Calls             : %d in %g seconds\n",nbSatCalls,totalTime4Sat);
  printf("c UNSAT Calls           : %d in %g seconds\n",nbUnsatCalls,totalTime4Unsat);
  printf("c--------------------------------------------------\n");


}


// NOTE: assumptions passed in member-variable 'assumptions'.
lbool Solver::solve_()
{

  if(incremental && certifiedUNSAT) {
    printf("Can not use incremental and certified unsat in the same time\n");
    exit(-1);
  }
    model.clear();
    conflict.clear();
    if (!ok) return l_False;
    double curTime = cpuTime();

    
    solves++;
            

    
    lbool   status        = l_Undef;
    if(!incremental && verbosity>=1) {
      printf("c ========================================[ MAGIC CONSTANTS ]==============================================\n");
      printf("c | Constants are supposed to work well together :-)                                                      |\n");
      printf("c | however, if you find better choices, please let us known...                                           |\n");
      printf("c |-------------------------------------------------------------------------------------------------------|\n");
    printf("c |                                |                                |                                     |\n"); 
    printf("c | - Restarts:                    | - Reduce Clause DB:            | - Minimize Asserting:               |\n");
    printf("c |   * LBD Queue    : %6d      |   * First     : %6d         |    * size < %3d                     |\n",lbdQueue.maxSize(),nbclausesbeforereduce,lbSizeMinimizingClause);
    printf("c |   * Trail  Queue : %6d      |   * Inc       : %6d         |    * lbd  < %3d                     |\n",trailQueue.maxSize(),incReduceDB,lbLBDMinimizingClause);
    printf("c |   * K            : %6.2f      |   * Special   : %6d         |                                     |\n",K,specialIncReduceDB);
    printf("c |   * R            : %6.2f      |   * Protected :  (lbd)< %2d     |                                     |\n",R,lbLBDFrozenClause);
    printf("c |                                |                                |                                     |\n"); 
printf("c ==================================[ Search Statistics (every %6d conflicts) ]=========================\n",verbEveryConflicts);
      printf("c |                                                                                                       |\n"); 

      printf("c |          RESTARTS           |          ORIGINAL         |              LEARNT              | Progress |\n");
      printf("c |       NB   Blocked  Avg Cfc |    Vars  Clauses Literals |   Red   Learnts    LBD2  Removed |          |\n");
      printf("c =========================================================================================================\n");
    }

    // Search:
    int curr_restarts = 0;
    while (status == l_Undef){
      status = search(0); // the parameter is useless in glucose, kept to allow modifications

        if (!withinBudget()) break;
        curr_restarts++;
    }

    if (!incremental && verbosity >= 1)
      printf("c =========================================================================================================\n");


    if (certifiedUNSAT){ // Want certified output
      if (status == l_False)
	fprintf(certifiedOutput, "0\n");
      fclose(certifiedOutput);
    }


    if (status == l_True){
        // Extend & copy model:
        model.growTo(nVars());
        for (int i = 0; i < nVars(); i++) model[i] = value(i);
    }else if (status == l_False && conflict.size() == 0)
        ok = false;



    cancelUntil(0);

    double finalTime = cpuTime();
    if(status==l_True) {
      nbSatCalls++; 
      totalTime4Sat +=(finalTime-curTime);
    }
    if(status==l_False) {
      nbUnsatCalls++; 
      totalTime4Unsat +=(finalTime-curTime);
    }


    return status;
}

//=================================================================================================
// Writing CNF to DIMACS:
// 
// FIXME: this needs to be rewritten completely.

static Var mapVar(Var x, vec<Var>& map, Var& max)
{
    if (map.size() <= x || map[x] == -1){
        map.growTo(x+1, -1);
        map[x] = max++;
    }
    return map[x];
}


void Solver::toDimacs(FILE* f, Clause& c, vec<Var>& map, Var& max)
{
    if (satisfied(c)) return;

    for (int i = 0; i < c.size(); i++)
        if (value(c[i]) != l_False)
            fprintf(f, "%s%d ", sign(c[i]) ? "-" : "", mapVar(var(c[i]), map, max)+1);
    fprintf(f, "0\n");
}


void Solver::toDimacs(const char *file, const vec<Lit>& assumps)
{
    FILE* f = fopen(file, "wr");
    if (f == NULL)
        fprintf(stderr, "could not open file %s\n", file), exit(1);
    toDimacs(f, assumps);
    fclose(f);
}


void Solver::toDimacs(FILE* f, const vec<Lit>& assumps)
{
    // Handle case when solver is in contradictory state:
    if (!ok){
        fprintf(f, "p cnf 1 2\n1 0\n-1 0\n");
        return; }

    vec<Var> map; Var max = 0;

    // Cannot use removeClauses here because it is not safe
    // to deallocate them at this point. Could be improved.
    int cnt = 0;
    for (int i = 0; i < clauses.size(); i++)
        if (!satisfied(ca[clauses[i]]))
            cnt++;
        
    for (int i = 0; i < clauses.size(); i++)
        if (!satisfied(ca[clauses[i]])){
            Clause& c = ca[clauses[i]];
            for (int j = 0; j < c.size(); j++)
                if (value(c[j]) != l_False)
                    mapVar(var(c[j]), map, max);
        }

    // Assumptions are added as unit clauses:
    cnt += assumptions.size();

    fprintf(f, "p cnf %d %d\n", max, cnt);

    for (int i = 0; i < assumptions.size(); i++){
        assert(value(assumptions[i]) != l_False);
        fprintf(f, "%s%d 0\n", sign(assumptions[i]) ? "-" : "", mapVar(var(assumptions[i]), map, max)+1);
    }

    for (int i = 0; i < clauses.size(); i++)
        toDimacs(f, ca[clauses[i]], map, max);

    if (verbosity > 0)
        printf("Wrote %d clauses with %d variables.\n", cnt, max);
}


//=================================================================================================
// Garbage Collection methods:

void Solver::relocAll(ClauseAllocator& to)
{
    // All watchers:
    //
    // for (int i = 0; i < watches.size(); i++)
    watches.cleanAll();
    watchesBin.cleanAll();
    for (int v = 0; v < nVars(); v++)
        for (int s = 0; s < 2; s++){
            Lit p = mkLit(v, s);
            // printf(" >>> RELOCING: %s%d\n", sign(p)?"-":"", var(p)+1);
            vec<Watcher>& ws = watches[p];
            for (int j = 0; j < ws.size(); j++)
                ca.reloc(ws[j].cref, to);
            vec<Watcher>& ws2 = watchesBin[p];
            for (int j = 0; j < ws2.size(); j++)
                ca.reloc(ws2[j].cref, to);
        }

    // All reasons:
    //
    for (int i = 0; i < trail.size(); i++){
        Var v = var(trail[i]);

        if (reason(v) != CRef_Undef && (ca[reason(v)].reloced() || locked(ca[reason(v)])))
            ca.reloc(vardata[v].reason, to);
    }

    // All learnt:
    //
    for (int i = 0; i < learnts.size(); i++)
        ca.reloc(learnts[i], to);

    // All original:
    //
    for (int i = 0; i < clauses.size(); i++)
        ca.reloc(clauses[i], to);
}


void Solver::garbageCollect()
{
    // Initialize the next region to a size corresponding to the estimated utilization degree. This
    // is not precise but should avoid some unnecessary reallocations for the new region:
    ClauseAllocator to(ca.size() - ca.wasted()); 

    relocAll(to);
    if (verbosity >= 2)
        printf("|  Garbage collection:   %12d bytes => %12d bytes             |\n", 
               ca.size()*ClauseAllocator::Unit_Size, to.size()*ClauseAllocator::Unit_Size);
    to.moveTo(ca);
}
