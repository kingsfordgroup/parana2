// Copyright (c) 2009, 2010, 2011 Paul Pogonyshev.
//
// Permission is hereby granted, free of charge, to any person obtaining a copy
// of this software and associated documentation files (the "Software"), to deal
// in the Software without restriction, including without limitation the rights
// to use, copy, modify, merge, publish, distribute, sublicense, and/or sell
// copies of the Software, and to permit persons to whom the Software is
// furnished to do so, subject to the following conditions:
//
// The above copyright notice and this permission notice shall be included in
// all copies or substantial portions of the Software.
//
// THE SOFTWARE IS PROVIDED "AS IS", WITHOUT WARRANTY OF ANY KIND, EXPRESS OR
// IMPLIED, INCLUDING BUT NOT LIMITED TO THE WARRANTIES OF MERCHANTABILITY,
// FITNESS FOR A PARTICULAR PURPOSE AND NONINFRINGEMENT. IN NO EVENT SHALL THE
// AUTHORS OR COPYRIGHT HOLDERS BE LIABLE FOR ANY CLAIM, DAMAGES OR OTHER
// LIABILITY, WHETHER IN AN ACTION OF CONTRACT, TORT OR OTHERWISE, ARISING FROM,
// OUT OF OR IN CONNECTION WITH THE SOFTWARE OR THE USE OR OTHER DEALINGS IN
// THE SOFTWARE.


#ifndef MCT_IMPL_CLOSED_HASH_TABLE_HPP
#define MCT_IMPL_CLOSED_HASH_TABLE_HPP


#include <cmath>
#include <cstddef>
#include <cstring>
#include <iterator>
#include <limits>
#include <memory>
#include <stdexcept>
#include <utility>

#if MCT_CXX0X_SUPPORTED
# include <initializer_list>
#endif

#include <mct/impl/optimization.hpp>
#include <mct/impl/utils.hpp>


// Internally data is stored in an array of buckets with length being some power of 2.
// Tables come in two quite different variants: those where buckets keep hash values of
// their keys and those where buckets don't.
//
// Hashes are kept in buckets:
// - each bucket has an additional field '_hash';
// - empty buckets have '_hash' set to 0, debris [*] have it set to 1;
// - hashes are always 'postprocessed' after they are computed by the actual hash
//   function so that hashes in used buckets are never 0 or 1;
// - this mode can be used with any key type.
//
// Hashes are recomputed each time:
// - for empty/debris buckets 0 or 1 is stored where key would be (we have to use some
//   reinterpret_cast black magic for this);
// - at most two buckets are designated as 'special': those are the buckets that have
//   keys that match the marks above, yet contain real values; such buckets can be located
//   anywhere in the table or be not present at all;
// - can only be used for 'small' keys: at most sizeof (long long).
//
// Default is to use the second mode (smaller memory footprint) if key size permits.
//
// If you try to improve this remember that key equality is not the same as bitwise
// equality, especially for structures and other custom types.  Even for plain integral
// type one could specify a different comparator function.  Additionally, you can't call
// hash and comparison functions on non-inserted data: they could assert something or
// dereference an invalid pointer as a result.
//
// [*] Debris are buckets that contained some value but were subsequently erased.  With
//     quadratic probing such buckets cannot be just reset to empty state.

// It was decided that tables never shrink on their own.  See documentation for rationale.


#if MCT_CXX0X_SUPPORTED
# define MCT_STD_FORWARD(type, value)  std::forward <type> (value)
#else
# define MCT_STD_FORWARD(type, value)  value
#endif


namespace mct
{

  namespace impl
  {

    // Emulates access to a memory chunk as to a 3-state (empty/debris/used) flag.  Since
    // that memory chunk will normally be occupied by something else, we call this a
    // 'union'.  Only specifications with zero as the second parameter will be used.
    template <std::size_t key_type_size, int used = 0>
    struct flag_pseudo_union
    {
      typedef  unsigned char  integral_type;

      // 0 is 'empty', 1 is 'debris', anything else is 'used'.
      static  integral_type
      compute (const void* storage)
      {
        for (const unsigned char*       first = static_cast <const unsigned char*> (storage),
                                * const last  = (first + (key_type_size - 1));
             first < last; ++first)
          {
            if (*first)
              return 2;
          }

        return *(static_cast <const unsigned char*> (storage) + (key_type_size - 1));
      }

      static  void
      mark_clear (void* storage)
      {
        std::memset (storage, 0, key_type_size);
      }

      static  void
      mark_destroyed (void* storage)
      {
        std::memset (storage, 0, key_type_size - 1);
        *(static_cast <unsigned char*> (storage) + (key_type_size - 1)) = 1;
      }
    };

    // Apparently, this is the standard and correct way to avoid pointer aliasing
    // problems.  Since we've already had such problems on GCC 4.5, let's be extra-safe.
    // A previous solution/workaround used converting 'storage' to 'volatile type*', but
    // that's not correct according to the standard, even if works on GCC 4.5.
    //
    // Hopefully, compilers are able to understand what's going on and optimize all the
    // complications away.  At least GCC seems to be able to do this, although
    // benchmarks do run a little bit slower.
    template <typename type>
    struct integral_pseudo_union
    {
      typedef  type  integral_type;

      static  integral_type
      compute (const void* storage)
      {
        integral_type  flag;
        std::memcpy (&flag, storage, sizeof (flag));
        return flag;
      }

      static  void
      mark_clear (void* storage)
      {
        static  const integral_type  flag = 0;
        std::memcpy (storage, &flag, sizeof (flag));
      }

      static  void
      mark_destroyed (void* storage)
      {
        static  const integral_type  flag = 1;
        std::memcpy (storage, &flag, sizeof (flag));
      }
    };


    // Specialize 'flag_pseudo_union' template for various integral type, taking care to
    // avoid repeated specialization with the same parameters for types of the same size
    // (which are different on various platforms, too).
    template <>
    struct flag_pseudo_union <sizeof (char), 0>
      : integral_pseudo_union <unsigned char>
    { };

    template <>
    struct flag_pseudo_union <sizeof (short), (sizeof (short) > sizeof (char) ? 0 : 1)>
      : integral_pseudo_union <unsigned short>
    { };

    template <>
    struct flag_pseudo_union <sizeof (int), (sizeof (int) > sizeof (short) ? 0 : 2)>
      : integral_pseudo_union <unsigned int>
    { };

    template <>
    struct flag_pseudo_union <sizeof (long), (sizeof (long) > sizeof (int) ? 0 : 3)>
      : integral_pseudo_union <unsigned long>
    { };

# if MCT_HAVE_LONG_LONG
    template <>
    struct flag_pseudo_union <sizeof (long long), (sizeof (long long) > sizeof (long) ? 0 : 4)>
      : integral_pseudo_union <unsigned long long>
    { };
# endif


    // Buckets consist of three relatively independent part:
    // - optional hash storage (used if 'keep_hashes' is true);
    // - optional link storage (used by 'linked_hash_*' and 'forward_hash_*' containers);
    // - required value storage.
    //
    // Additionally, bucket usage flag is stored in one of these parts, but not as a
    // separate field.  Instead, two values of a field (currently always 0 and 1, except
    // in pointer link storage) are designated as 'special'.  If this flag is stored in
    // the last part --- which is, actually, the most common case, --- additional support
    // from container itself is required.
    //
    // Note that is_empty() and similar functions are now static and require a previous
    // call to get_usage_data().  This is done because the latter is in some cases
    // non-trivial to compute, so having it computed separately one time only is
    // beneficial.  In particular this is noticeable when 'flag_pseudo_union' is used.


    // There is no 'store_usage_flag' parameter: if hash storage is used at all, it is
    // also used to store the usage flag, without exceptions.  This way,
    // used_and_matches() function can do its work in a single operation.
    template <bool keep_hashes>
    class hash_storage;


    template <>
    class hash_storage <false>
    {
    public:

      static  const bool  KEEPS_HASHES = false;


      static  std::size_t
      postprocess_hash (std::size_t hash)
      {  return hash;  }


      // Unused, but needed for code validity.
      std::size_t
      hash () const
      {  return 0;  }

      void
      store_hash (std::size_t hash)
      { }

      void
      copy_hash (const hash_storage& that)
      { }


      void
      pretend_to_be_used_hash ()
      { }
    };


    template <>
    class hash_storage <true>
    {
    public:

      typedef  std::size_t  usage_data;


      static  const bool  KEEPS_HASHES = true;

    protected:

      std::size_t  _hash;


    public:

      // Since '_hash' field becomes a usage flag, two values need to be reserved.
      static  std::size_t
      postprocess_hash (std::size_t hash)
      {
        return MCT_OPTIMIZATION_LIKELY (hash >= 2) ? hash : hash + 2;
      }


      std::size_t
      hash () const
      {  return _hash;  }

      void
      store_hash (std::size_t hash)
      {  this->_hash = hash;  }

      void
      copy_hash (const hash_storage& that)
      {  this->_hash = that._hash;  }


      usage_data
      get_usage_data () const
      {  return _hash;  }


      // Since 'hash' parameter has to first be passed through the postprocess_hash()
      // function, this single comparison achives both 'used' test and 'matches' test.
      static  bool
      used_and_matches (usage_data usage, std::size_t hash)
      {  return usage == hash;  }

      static  bool
      is_empty (usage_data usage)
      {  return usage == 0;  }

      static  bool
      is_debris (usage_data usage)
      {  return usage == 1;  }

      static  bool
      is_empty_or_debris (usage_data usage)
      {  return usage <= 1;  }


      void
      mark_clear ()
      {  _hash = 0;  }

      void
      mark_destroyed ()
      {  _hash = 1;  }


      void
      revert_to_empty ()
      { }

      void
      revert_to_debris ()
      { }


      void
      pretend_to_be_used_hash ()
      {
        // Anything but 0 and 1 will do.
        this->_hash = 2;
      }
    };


    template <bool store_links, bool pointer_links, bool store_usage_flag>
    class single_link_storage;


    // Note that first/next/previous methods are not defined and are up to subclasses.
    template <bool pointer_links, bool store_usage_flag>
    class single_link_storage <false, pointer_links, store_usage_flag>
    {
    public:

      static  const bool  USES_LINKS = false;


      template <typename bucket_type>
      void
      link (bucket_type* that)
      { }

      template <typename bucket_type>
      void
      unlink ()
      { }


      static  std::size_t
      limit_max_count (std::size_t max_count)
      {  return max_count;  }


#   if MCT_DEBUGGING_MEMBERS
      template <typename bucket_type>
      static  void
      validate_links (const bucket_type* sentinel, std::size_t num_used)
      { }
#   endif
    };


    template <>
    class single_link_storage <true, false, false>
    {
    public:

      static  const bool  USES_LINKS = true;
      static  const bool  NULL_END   = false;

    protected:

      typedef  int  next_field_type;

      // Note that this an _offset_ (possibly negative) to the next bucket, not bucket
      // index in the array.  This is probably faster and, more importantly, allows us to
      // have iterators that don't keep track of the data: with indices you'd have to know
      // where array beginning is, with offsets you don't.  The only drawback is that it
      // halves upper limit on bucket count, but 1 billion buckets should still be enough.
      int  _next;


    public:

      template <typename bucket_type>
      void
      link (bucket_type* that)
      {
        this->_next  = (that - static_cast <bucket_type*> (this)) + that->_next;
        that->_next -= this->_next;
      }

      void
      unlink (single_link_storage* after)
      {  after->_next += this->_next;  }

      template <typename bucket_type>
      void
      set_next (bucket_type* next)
      {
        this->_next = (next - static_cast <bucket_type*> (this));
      }


#   if MCT_DEBUGGING_MEMBERS

      template <typename bucket_type>
      static  void
      validate_links (const bucket_type* sentinel, std::size_t num_used)
      {
        if (num_used == 0)
          {
            if (sentinel->_next)
              throw std::logic_error ("wrong link sentinel in an empty table");
            return;
          }

        std::size_t         num_linked = 0;
        const bucket_type*  bucket;

        for (bucket = (sentinel + sentinel->_next); num_linked != num_used;
             ++num_linked, bucket += bucket->_next)
          {
            if (bucket == sentinel)
              throw std::logic_error ("not all used buckets in the table are linked together");
          }

        if (bucket != sentinel)
          throw std::logic_error ("too many buckets are linked together, or they form a loop");
      }

#   endif  // MCT_DEBUGGING_MEMBERS


      template <typename Data>
      static  typename Data::bucket_type*
      first (const Data& data)
      {
        return (static_cast <typename Data::bucket_type*> (data.sentinel ())
                + data.sentinel ()->_next);
      }

      template <typename bucket_type>
      bucket_type*
      next ()
      {
        return static_cast <bucket_type*> (this) + _next;
      }

      template <typename Data>
      typename Data::bucket_type*
      next (const Data& data)
      {
        return static_cast <typename Data::bucket_type*> (this) + _next;
      }


      void
      initialize_as_sentinel ()
      {
        _next = 0;
      }


      static  std::size_t
      limit_max_count (std::size_t max_count)
      {
        return std::min (max_count,
                         static_cast <std::size_t> ((std::numeric_limits <int>::max () >> 1) + 1));
      }


    protected:

      // I was unable to make these into static member variables because in the other
      // 'single_link_storage' specialization you'd need a reinterpret_cast() to a pointer
      // type, which is not allowed in constant expressions.  In this specialization we
      // use a non-const expression either, but it could be replaced with some bit
      // trickery if there was any reason to.

      static  int
      empty_value ()
      {  return std::numeric_limits <int>::max ();  }

      static  int
      debris_value ()
      {  return std::numeric_limits <int>::max () - 1;  }
    };


    template <>
    class single_link_storage <true, true, false>
    {
    public:

      static  const bool  USES_LINKS = true;
      static  const bool  NULL_END   = true;

    protected:

      typedef  single_link_storage*  next_field_type;

      single_link_storage*  _next;


    public:

      template <typename bucket_type>
      void
      link (bucket_type* that)
      {
        this->_next = that->_next;
        that->_next = this;
      }

      void
      unlink (single_link_storage* after)
      {  after->_next = this->_next;  }

      template <typename bucket_type>
      void
      set_next (bucket_type* next)
      {  this->_next = next;  }


#   if MCT_DEBUGGING_MEMBERS

      template <typename bucket_type>
      static  void
      validate_links (const bucket_type* sentinel, std::size_t num_used)
      {
        if (num_used == 0)
          {
            if (sentinel->_next)
              throw std::logic_error ("wrong link sentinel in an empty table");
            return;
          }

        std::size_t         num_linked = 0;
        const bucket_type*  bucket;

        for (bucket = static_cast <bucket_type*> (sentinel->_next); num_linked != num_used;
             ++num_linked, bucket = static_cast <bucket_type*> (bucket->_next))
          {
            if (!bucket)
              throw std::logic_error ("not all used buckets in the table are linked together");
          }

        if (bucket)
          throw std::logic_error ("too many buckets are linked together, or they form a loop");
      }

#   endif  // MCT_DEBUGGING_MEMBERS


      template <typename Data>
      static  typename Data::bucket_type*
      first (const Data& data)
      {
        return static_cast <typename Data::bucket_type*> (data.sentinel ()->_next);
      }

      template <typename bucket_type>
      bucket_type*
      next ()
      {
        return static_cast <bucket_type*> (_next);
      }

      template <typename Data>
      typename Data::bucket_type*
      next (const Data& data)
      {
        return static_cast <typename Data::bucket_type*> (_next);
      }


      void
      initialize_as_sentinel ()
      {
        _next = 0;
      }


      static  std::size_t
      limit_max_count (std::size_t max_count)
      {  return max_count;  }


    protected:

      // Pointer-based link storage uses -1 and -2 to mark empty and debris buckets,
      // because null pointer here is a valid value for forward tables.
      //
      // Note that this is theoretically unsound (-1 or -2 could in principle be valid
      // bucket addresses), but:
      // - since buckets will be aligned, valid addresses are at least multiplies of 4;
      // - if integers are two-complemented (which in reality they are), there's just no
      //   space at -1/-2 for a full bucket structure.
      //
      // In other words, I'm sure that in practice -1/-2 cannot be valid bucket addresses.

      static  single_link_storage* const
      empty_value ()
      {  return reinterpret_cast <single_link_storage* const> (-1);  }

      static  single_link_storage* const
      debris_value ()
      {  return reinterpret_cast <single_link_storage* const> (-2);  }
    };


    template <bool pointer_links>
    class single_link_storage <true, pointer_links, true>
      : public single_link_storage <true, pointer_links, false>
    {
      typedef  single_link_storage <true, pointer_links, false>  base_type;

    public:

      typedef  typename base_type::next_field_type               usage_data;


      usage_data
      get_usage_data () const
      {  return this->_next;  }


      static  bool
      used_and_matches (usage_data usage, std::size_t hash)
      {  return !is_empty_or_debris (usage);  }

      static  bool
      is_empty (usage_data usage)
      {  return usage == single_link_storage::empty_value ();  }

      static  bool
      is_debris (usage_data usage)
      {  return usage == single_link_storage::debris_value ();  }

      static  bool
      is_empty_or_debris (usage_data usage)
      {
        // Non-brain-dead compilers should be able to optimize this to single instruction
        // anyway.  We avoid assuming anything about how negative integers are
        // implemented, even if in reality only two-complement model is used.
        return (usage    == single_link_storage::empty_value  ()
                || usage == single_link_storage::debris_value ());
      }


      // Don't do anything because '_next' field is not changed.
      void
      revert_to_empty ()
      { }

      void
      revert_to_debris ()
      { }


      void
      mark_clear ()
      {  this->_next = this->empty_value ();  }

      void
      mark_destroyed ()
      {  this->_next = this->debris_value ();  }
    };


    // We don't have 'store_links' flag as this class is used in cases when it would
    // always be true anyway.
    template <bool pointer_links, bool store_usage_flag>
    class double_link_storage;


    template <bool store_usage_flag>
    class double_link_storage <false, store_usage_flag>
      : public single_link_storage <true, false, store_usage_flag>
    {
    protected:

      int  _previous;


    public:

      template <typename bucket_type>
      void
      link (bucket_type* that)
      {
        this->_next     = (that - static_cast <bucket_type*> (this));
        this->_previous = this->_next + that->_previous;

        (static_cast <bucket_type*> (this) + this->_previous)->_next = -this->_previous;
        that->_previous                                              = -this->_next;
      }

      template <typename bucket_type>
      void
      unlink ()
      {
        (static_cast <bucket_type*> (this) + this->_previous)->_next += this->_next;
        (static_cast <bucket_type*> (this) + this->_next)->_previous += this->_previous;
      }

      void
      reverse_links ()
      {
        using std::swap;
        swap (this->_next, this->_previous);
      }


#   if MCT_DEBUGGING_MEMBERS

      template <typename bucket_type>
      static  void
      validate_links (const bucket_type* sentinel, std::size_t num_used)
      {
        if (num_used == 0)
          {
            if (sentinel->_next || sentinel->_previous)
              throw std::logic_error ("wrong link sentinel in an empty table");
            return;
          }

        std::size_t         loop_size = 0;
        const bucket_type*  bucket    = (sentinel + sentinel->_next);

        for (; loop_size != num_used; ++loop_size, bucket += bucket->_next)
          {
            if (bucket == sentinel)
              throw std::logic_error ("link loop doesn't contain all the values in the table");

            if (((bucket + bucket->_next) + (bucket + bucket->_next)->_previous) != bucket
                || ((bucket + bucket->_previous) + (bucket + bucket->_previous)->_next) != bucket)
              throw std::logic_error ("broken link structure");
          }

        if (bucket != sentinel)
          throw std::logic_error ("link loop size is larger than number of values in the table");
      }

#   endif  // MCT_DEBUGGING_MEMBERS


      template <typename bucket_type>
      bucket_type*
      previous ()
      {
        return static_cast <bucket_type*> (this) + _previous;
      }

      template <typename Data>
      typename Data::bucket_type*
      previous (const Data& data)
      {
        return previous <typename Data::bucket_type> ();
      }


      void
      initialize_as_sentinel ()
      {
        this->_next     = 0;
        this->_previous = 0;
      }
    };


    template <bool store_usage_flag>
    class double_link_storage <true, store_usage_flag>
      : public single_link_storage <true, true, store_usage_flag>
    {
    protected:

      double_link_storage*  _previous;


    public:

      template <typename bucket_type>
      void
      link (bucket_type* that)
      {
        this->_next            = that;
        this->_previous        = that->_previous;
        this->_previous->_next = this;
        that->_previous        = this;
      }

      template <typename bucket_type>
      void
      unlink ()
      {
        this->_previous->_next                              = this->_next;
        static_cast <bucket_type*> (this->_next)->_previous = this->_previous;
      }

      void
      reverse_links ()
      {
        double_link_storage* const  next = static_cast <double_link_storage*> (this->_next);

        this->_next     = this->_previous;
        this->_previous = next;
      }


#   if MCT_DEBUGGING_MEMBERS

      template <typename bucket_type>
      static  void
      validate_links (const bucket_type* sentinel, std::size_t num_used)
      {
        if (num_used == 0)
          {
            if (sentinel->_next != sentinel || sentinel->_previous != sentinel)
              throw std::logic_error ("wrong link sentinel in an empty table");
            return;
          }

        std::size_t         loop_size = 0;
        const bucket_type*  bucket    = static_cast <const bucket_type*> (sentinel->_next);

        for (; loop_size != num_used;
             ++loop_size, bucket = static_cast <const bucket_type*> (bucket->_next))
          {
            if (bucket == sentinel)
              throw std::logic_error ("link loop doesn't contain all the values in the table");

            if (static_cast <const bucket_type*> (bucket->_next)->_previous != bucket
                || bucket->_previous->_next != bucket)
              throw std::logic_error ("broken link structure");
          }

        if (bucket != sentinel)
          throw std::logic_error ("link loop size is larger than number of values in the table");
      }

#   endif  // MCT_DEBUGGING_MEMBERS


      template <typename bucket_type>
      bucket_type*
      previous ()
      {
        return static_cast <bucket_type*> (_previous);
      }

      template <typename Data>
      typename Data::bucket_type*
      previous (const Data& data)
      {
        return previous <typename Data::bucket_type> ();
      }


      void
      initialize_as_sentinel ()
      {
        this->_next     = this;
        this->_previous = this;
      }
    };


    template <typename Traits, bool store_usage_flag, bool with_external_use>
    class value_storage;


    template <typename Traits, bool with_external_use>
    class value_storage <Traits, false, with_external_use>
    {
    public:

      typedef  typename Traits::value_type             value_type;
      typedef  typename Traits::key_type               key_type;

    protected:

      typedef  typename intrusive_storage <value_type>::type  storage_type;

      storage_type  _value;


    public:

      // In almost all cases we don't do anything as the mark is indirectly set by
      // constructing value object or storing hash etc.
      void
      mark_used ()
      { }


      const key_type&
      key () const
      {  return Traits::extract_key (_value);  }


      value_type&
      value ()
      {  return _value;  }

      const value_type&
      value () const
      {  return _value;  }


      void
      pretend_to_be_used_value ()
      { }
    };


    template <typename Traits>
    class value_storage <Traits, true, false> : public value_storage <Traits, false, false>
    {
      typedef  value_storage <Traits, false, false>                       base_type;
      typedef  flag_pseudo_union <sizeof (typename base_type::key_type)>  flag_pseudo_union_type;

    public:

      typedef  typename flag_pseudo_union_type::integral_type             usage_data;


      usage_data
      get_usage_data () const
      {
        return flag_pseudo_union_type::compute (static_cast <const void*> (&this->_value));
      }


      static  bool
      used_and_matches (usage_data usage, std::size_t hash)
      {  return !is_empty_or_debris (usage);  }

      static  bool
      is_empty (usage_data usage)
      {  return usage == 0;  }

      static  bool
      is_debris (usage_data usage)
      {  return usage == 1;  }

      static  bool
      is_empty_or_debris (usage_data usage)
      {  return usage <= 1;  }


      void
      mark_clear ()
      {
        flag_pseudo_union_type::mark_clear (static_cast <void*> (&this->_value));
      }

      void
      mark_destroyed ()
      {
        flag_pseudo_union_type::mark_destroyed (static_cast <void*> (&this->_value));
      }


      void
      revert_to_empty ()
      {  mark_clear ();  }

      void
      revert_to_debris ()
      {  mark_destroyed ();  }
    };


    template <typename Traits>
    class value_storage <Traits, true, true> : public value_storage <Traits, false, false>
    {
      typedef  value_storage <Traits, false, false>                     base_type;
      typedef  hackish_external_use <typename base_type::storage_type>  accessor;

    public:

      typedef  typename base_type::value_type  value_type;

      typedef  typename accessor::value_type   usage_data;


      usage_data
      get_usage_data () const
      {  return accessor::get (this->_value);  }


      static  bool
      used_and_matches (usage_data usage, std::size_t hash)
      {  return !is_empty_or_debris (usage);  }

      static  bool
      is_empty (usage_data usage)
      {  return usage == 0;  }

      static  bool
      is_debris (usage_data usage)
      {  return usage == 1;  }

      static  bool
      is_empty_or_debris (usage_data usage)
      {
        // We leave optimization of this statment to the compiler, to avoid dealing with
        // signed/unsigned business.
        return usage == 0 || usage == 1;
      }


      void
      mark_clear ()
      {  accessor::set (this->_value, 0);  }

      void
      mark_destroyed ()
      {  accessor::set (this->_value, 1);  }

      void
      mark_used ()
      {  accessor::set (this->_value, 2);  }


      // We rely on the fact that revert_to_*() is called only when allocator throws on
      // constructing a value.  When that happens, mark_used() hadn't been called yet, so
      // there's nothing to revert.

      void
      revert_to_empty ()
      { }

      void
      revert_to_debris ()
      { }


      void
      pretend_to_be_used_value ()
      {  accessor::set (this->_value, 2);  }
    };


    // The common template is used if at least one of 'keep_hashes' or 'with_external_use'
    // is true.  Otherwise, the specialization below is used.
    template <typename Traits, bool keep_hashes,
              bool with_external_use
                = (supports_hackish_external_use
                   <typename intrusive_storage <typename Traits::value_type>::type>::value)>
    class plain_bucket_base
      : public hash_storage <keep_hashes>,
        public single_link_storage <false, false, false>,
        public value_storage <Traits, !keep_hashes && with_external_use, with_external_use>,
        public Traits
    {
    public:

      typedef  typename Traits::key_type        key_type;
      typedef  typename Traits::value_type      value_type;
      typedef  std::bidirectional_iterator_tag  associated_iterator_category;

      // Sentinel here is needed only as a stopping point for iterator incrementing.  It
      // is artificially marked as a used bucket so that next() stops on it.
      static  const bool  SPECIAL_EMPTY_AND_DEBRIS = false;
      static  const bool  USES_SENTINEL            = true;
      static  const bool  IS_FORWARD               = false;


      template <typename Data>
      static  typename Data::bucket_type*
      first (const Data& data)
      {
        typedef  typename Data::bucket_type  bucket_type;

        bucket_type* const  _end = data.storage_end ();

        if (data.num_used != 0)
          {
            for (bucket_type* bucket = data.buckets; bucket != _end; ++bucket)
              {
                if (!bucket_type::is_empty_or_debris (bucket->get_usage_data ()))
                  return bucket;
              }
          }

        return _end;
      }

      template <typename bucket_type>
      bucket_type*
      next ()
      {
        for (bucket_type* bucket = static_cast <bucket_type*> (this) + 1; ; ++bucket)
          {
            // Because of the sentinel, this will stop for valid iterators.
            if (!bucket_type::is_empty_or_debris (bucket->get_usage_data ()))
              return bucket;
          }
      }

      template <typename bucket_type>
      bucket_type*
      previous ()
      {
        for (bucket_type* bucket = static_cast <bucket_type*> (this) - 1; ; --bucket)
          {
            // This will stop for valid iterators except for table begin, but decrementing
            // that is undefined behavior anyway.
            if (!bucket_type::is_empty_or_debris (bucket->get_usage_data ()))
              return bucket;
          }
      }

      template <typename Data>
      typename Data::bucket_type*
      next (const Data& data)
      {
        return next <typename Data::bucket_type> ();
      }

      template <typename Data>
      typename Data::bucket_type*
      previous (const Data& data)
      {
        return previous <typename Data::bucket_type> ();
      }


      void
      initialize_as_sentinel ()
      {
        this->pretend_to_be_used_hash ();
        this->pretend_to_be_used_value ();
      }
    };


    // In case of no hash or external use field (which is actually very common despite
    // being 1 case in 4), we need to keep track of two special buckets.
    template <typename Traits>
    class plain_bucket_base <Traits, false, false>
      : public hash_storage <false>,
        public single_link_storage <false, false, false>,
        public value_storage <Traits, true, false>,
        public Traits
    {
    public:

      typedef  typename Traits::key_type        key_type;
      typedef  typename Traits::value_type      value_type;
      typedef  std::bidirectional_iterator_tag  associated_iterator_category;

      static  const bool  SPECIAL_EMPTY_AND_DEBRIS = true;
      static  const bool  USES_SENTINEL            = false;
      static  const bool  IS_FORWARD               = false;


      template <typename Data>
      static  typename Data::bucket_type*
      first (const Data& data)
      {
        typedef  typename Data::bucket_type  bucket_type;

        bucket_type* const  _end = data.storage_end ();

        if (data.num_used != 0)
          {
            for (bucket_type* bucket = data.buckets; bucket != _end; ++bucket)
              {
                if (!bucket_type::is_empty_or_debris (bucket->get_usage_data ())
                    || MCT_OPTIMIZATION_UNLIKELY (data.is_special (bucket)))
                  return bucket;
              }
          }

        return _end;
      }

      template <typename Data>
      typename Data::bucket_type*
      next (const Data& data)
      {
        typedef  typename Data::bucket_type  bucket_type;

        bucket_type*  bucket = static_cast <bucket_type*> (this);
        for (bucket_type* const limit = data.storage_end ();
             MCT_OPTIMIZATION_LIKELY (++bucket != limit);)
          {
            if (!bucket_type::is_empty_or_debris (bucket->get_usage_data ())
                || MCT_OPTIMIZATION_UNLIKELY (data.is_special (bucket)))
              break;
          }

        return bucket;
      }

      template <typename Data>
      typename Data::bucket_type*
      previous (const Data& data)
      {
        typedef  typename Data::bucket_type  bucket_type;

        for (bucket_type* bucket = static_cast <bucket_type*> (this) - 1; ; --bucket)
          {
            if (!bucket_type::is_empty_or_debris (bucket->get_usage_data ())
                || MCT_OPTIMIZATION_UNLIKELY (data.is_special (bucket)))
              return bucket;
          }
      }


      // Unused, but needed for code validity.
      void
      initialize_as_sentinel ()
      { }
    };


    template <typename Traits, bool pointer_links, bool keep_hashes>
    class linked_bucket_base : public hash_storage <keep_hashes>,
                               public double_link_storage <pointer_links, !keep_hashes>,
                               public value_storage <Traits, false, false>,
                               public Traits
    {
    public:

      typedef  typename Traits::key_type        key_type;
      typedef  typename Traits::value_type      value_type;
      typedef  std::bidirectional_iterator_tag  associated_iterator_category;

      static  const bool  SPECIAL_EMPTY_AND_DEBRIS = false;
      static  const bool  USES_SENTINEL            = true;
      static  const bool  IS_FORWARD               = false;
    };


    template <typename Traits, bool pointer_links, bool keep_hashes>
    class forward_bucket_base : public hash_storage <keep_hashes>,
                                public single_link_storage <true, pointer_links, !keep_hashes>,
                                public value_storage <Traits, false, false>,
                                public Traits
    {
    public:

      typedef  typename Traits::key_type    key_type;
      typedef  typename Traits::value_type  value_type;
      typedef  std::forward_iterator_tag    associated_iterator_category;

      static  const bool  SPECIAL_EMPTY_AND_DEBRIS = false;
      static  const bool  USES_SENTINEL            = true;
      static  const bool  IS_FORWARD               = true;
    };


    // Note that 'buckets' and 'num_buckets' are not initialized.
    template <typename Bucket, bool with_special_buckets = Bucket::SPECIAL_EMPTY_AND_DEBRIS>
    struct hash_table_data_base;


    template <typename Bucket>
    struct hash_table_data_base <Bucket, false>
    {
      typedef  Bucket                            bucket_type;
      typedef  typename bucket_type::usage_data  bucket_usage_data;

      // Pointer to the beginning of allocated bucket table.  Can be zero if the hash
      // table is empty; all functions have to handle this case correctly.
      bucket_type*  buckets;

      // Number of buckets pointed to by 'buckets'; always a power of 2.  If 'buckets' is
      // 0, this is the number of buckets that would have been allocated if the table
      // wasn't empty, i.e. DEFAULT_NUM_BUCKETS or a value forced by user.  In some cases
      // a sentinel bucket is allocated at the end of 'buckets' array, but it is never
      // included in this count.
      std::size_t   num_buckets;

      // Number of buckets filled with actual elements, i.e. hash table size.
      std::size_t   num_used;

      // Number of buckets filled with elements or debris.
      std::size_t   num_occupied;


      hash_table_data_base ()
        : num_used     (0),
          num_occupied (0)
      { }


      bool
      is_special (const bucket_type* bucket) const
      {  return false;  }

      void
      note_new (bucket_type* bucket)
      { }

      void
      note_being_destroyed (bucket_type* bucket)
      { }


      void
      clear ()
      {
        num_used     = 0;
        num_occupied = 0;
      }


      void
      swap (hash_table_data_base& that)
      {
        using std::swap;

        swap (buckets,      that.buckets);
        swap (num_buckets,  that.num_buckets);
        swap (num_used,     that.num_used);
        swap (num_occupied, that.num_occupied);
      }
    };


    template <typename Bucket>
    struct hash_table_data_base <Bucket, true> : public hash_table_data_base <Bucket, false>
    {
      typedef  hash_table_data_base <Bucket, false>   base_type;
      typedef  typename base_type::bucket_type        bucket_type;
      typedef  typename base_type::bucket_usage_data  bucket_usage_data;


      bucket_type*  special_buckets[2];


      hash_table_data_base ()
      {
        special_buckets[0] = 0;
        special_buckets[1] = 0;
      }


      bool
      is_special (const bucket_type* bucket) const
      {
        return MCT_OPTIMIZATION_UNLIKELY (bucket    == special_buckets[0]
                                          || bucket == special_buckets[1]);
      }


      // Note that we shouldn't pass usage data to this function.  Usage data would
      // normally be computed before creating a value, so would wrongly indicate an empty
      // bucket.  And computing it right before calling note_new() is pointless --- the
      // function can just do that itself.
      void
      note_new (bucket_type* bucket)
      {
        const bucket_usage_data  usage = bucket->get_usage_data ();
        if (MCT_OPTIMIZATION_UNLIKELY (bucket_type::is_empty_or_debris (usage)))
          special_buckets[bucket_type::is_empty (usage) ? 0 : 1] = bucket;
      }

      void
      note_being_destroyed (bucket_type* bucket)
      {
        const bucket_usage_data  usage = bucket->get_usage_data ();
        if (MCT_OPTIMIZATION_UNLIKELY (bucket_type::is_empty_or_debris (usage)))
          special_buckets[bucket_type::is_empty (usage) ? 0 : 1] = 0;
      }


      void
      clear ()
      {
        base_type::clear ();
        special_buckets[0] = 0;
        special_buckets[1] = 0;
      }


      void
      swap (hash_table_data_base& that)
      {
        using std::swap;

        base_type::swap (that);

        // At least VC 2005 cannot swap arrays, bug #549067.
        swap (special_buckets[0], that.special_buckets[0]);
        swap (special_buckets[1], that.special_buckets[1]);
      }
    };


    template <typename Bucket, bool forward = Bucket::IS_FORWARD>
    struct hash_table_data;


    template <typename Bucket>
    struct hash_table_data <Bucket, false> : public hash_table_data_base <Bucket>
    {
      typedef  hash_table_data_base <Bucket>    base_type;
      typedef  typename base_type::bucket_type  bucket_type;


      std::size_t
      start_probing (std::size_t hash) const
      {
        // We require that number of buckets is a power of 2.  To compensate, we involve a
        // prime number here; this usually improves speed with particularly bad hash
        // function choices.
        return (hash * size_t_constants <>::HUGE_PRIME) & (this->num_buckets - 1);
      }

      std::size_t
      continue_probing (std::size_t looking_at, std::size_t iteration) const
      {
        // Quadratic probing.
        return (looking_at + iteration) & (this->num_buckets - 1);
      }


      bucket_type*
      begin () const
      {
        return MCT_OPTIMIZATION_LIKELY (this->num_used) ? bucket_type::first (*this) : end ();
      }

      bucket_type*
      end () const
      {
        return this->buckets + this->num_buckets;
      }

      bucket_type*
      sentinel () const
      {
        return this->buckets + this->num_buckets;
      }

      bucket_type*
      storage_end () const
      {
        return this->buckets + this->num_buckets;
      }


      // Unused, but needed for code validity.
      bucket_type*
      get_back () const
      {  return 0;  }

      void
      reset_back ()
      { }


#   if MCT_DEBUGGING_MEMBERS

      bool
      valid_bucket (const bucket_type* bucket) const
      {
        return (this->buckets
                && this->buckets <= bucket && bucket < end ()
                // Check that it doesn't point into a middle of an actual bucket.
                && ((reinterpret_cast <const char*> (bucket)
                     - reinterpret_cast <const char*> (this->buckets))
                    % sizeof (bucket_type) == 0)
                && (!bucket_type::is_empty_or_debris (bucket->get_usage_data ())
                    || this->is_special (bucket)));
      }

#   endif
    };


    template <typename Bucket>
    struct hash_table_data <Bucket, true> : public hash_table_data <Bucket, false>
    {
      typedef  hash_table_data <Bucket, false>        base_type;
      typedef  typename base_type::bucket_type        bucket_type;
      typedef  typename base_type::bucket_usage_data  bucket_usage_data;

      bucket_type*  back;


      void
      note_linked (bucket_type* bucket)
      {
        // Note that end() is not a simple constant anymore.
        if (bucket->template next <bucket_type> () == end ())
          back = bucket;
      }

      void
      note_being_unlinked (bucket_type* bucket, bucket_type* after)
      {
        if (back == bucket)
          back = after;
      }

      void
      note_new (bucket_type* bucket)
      {
        base_type::note_new (bucket);
        note_linked (bucket);
      }

      void
      note_being_destroyed (bucket_type* bucket, bucket_type* after)
      {
        base_type::note_being_destroyed (bucket);
        note_being_unlinked (bucket, after);
      }


      // The same, as in the base class, but still needed since end() is different.
      bucket_type*
      begin () const
      {
        return MCT_OPTIMIZATION_LIKELY (this->num_used) ? bucket_type::first (*this) : end ();
      }

      bucket_type*
      end () const
      {
        // Probably not worth it to have yet another template parameter.
        return Bucket::NULL_END ? 0 : this->storage_end ();
      }


      bucket_type*
      get_back () const
      {  return back;  }

      void
      reset_back ()
      {  back = this->storage_end ();  }

      void
      clear ()
      {
        base_type::clear ();
        reset_back ();
      }


      void
      swap (hash_table_data& that)
      {
        using std::swap;

        base_type::swap (that);
        swap (this->back, that.back);
      }
    };


    // We use no-data iterators only in standard (release) mode and during testing.  If
    // preconditions checking is requested (unless when testing), _all_ iterators will
    // have '_data' member, because it is required to robustly validate iterators,
    // especially before decrementing them.

# if !MCT_CHECK_PRECONDITIONS || MCT_USE_EFFICIENT_ITERATORS


    // Note that types that use special bucket pointer _need_ to have a pointer to data in
    // the iterator, because otherwise it is impossible to tell 'special' buckets from
    // simply empty/debris buckets.
    template <typename Bucket, bool keep_data_pointer = Bucket::SPECIAL_EMPTY_AND_DEBRIS>
    class hash_table_iterator_base;

    template <typename Bucket>
    class hash_table_iterator_base <Bucket, false>
    {
      typedef  Bucket  bucket_type;

    protected:

      bucket_type*  _bucket;

    public:

      // Not documented; used only in MCT tests.
      static  const bool  ROBUST_PRECONDITION_CHECKING = false;


      hash_table_iterator_base ()
      { }

      hash_table_iterator_base (const hash_table_data <bucket_type>& data, bucket_type* bucket)
        : _bucket (bucket)
      { }


    protected:

      void
      move_forward ()
      {
        MCT_PRECONDITION (valid (), "invalid iterator");
        _bucket = _bucket->template next <bucket_type> ();
      }

      void
      move_backward ()
      {
        // Cannot check anything without knowing the table the iterator belongs to.
        _bucket = _bucket->template previous <bucket_type> ();
      }

      bool
      valid () const
      {
        // Required, but not enough: will not catch cases where '_bucket' points outside
        // its hash table (most importantly, past-the-end iterator).
        return !bucket_type::is_empty_or_debris (_bucket->get_usage_data ());
      }


    public:

      // Supposedly not for outsiders to use.

      bucket_type*
      bucket () const
      {  return _bucket;  }

      bool
      is_for (const hash_table_data <bucket_type>& data) const
      {  return true;  }
    };


# else  // MCT_CHECK_PRECONDITIONS && !MCT_USE_EFFICIENT_ITERATORS

    // See comment earlier for explanation.
    template <typename Bucket, bool keep_data_pointer = true>
    class hash_table_iterator_base;

# endif


    // Iterators with '_data' members are defined both when precondition checking is used
    // and when not.
    template <typename Bucket>
    class hash_table_iterator_base <Bucket, true>
    {
      typedef  Bucket  bucket_type;

    protected:

      const hash_table_data <bucket_type>*  _data;
      bucket_type*                          _bucket;

    public:

      // Not documented; used only in MCT tests.
      static  const bool  ROBUST_PRECONDITION_CHECKING = true;


      hash_table_iterator_base ()
      { }

      hash_table_iterator_base (const hash_table_data <bucket_type>& data, bucket_type* bucket)
        : _data   (&data),
          _bucket (bucket)
      { }


    protected:

      void
      move_forward ()
      {
        MCT_PRECONDITION (valid (), "invalid iterator");
        _bucket = _bucket->next (*_data);
      }

      void
      move_backward ()
      {
        MCT_PRECONDITION (_bucket == _data->end () || (_bucket != _data->begin () && valid ()),
                          "iterator not valid for decrementing");
        _bucket = _bucket->previous (*_data);
      }

      bool
      valid () const
      {  return _data->valid_bucket (_bucket);  }


    public:

      // Supposedly not for outsiders to use.

      bucket_type*
      bucket () const
      {  return _bucket;  }

      bool
      is_for (const hash_table_data <bucket_type>& data) const
      {  return &data == _data;  }
    };


    template <typename Bucket>
    class hash_table_const_iterator : public hash_table_iterator_base <Bucket>
    {
      typedef  Bucket                                              bucket_type;
      typedef  hash_table_iterator_base <bucket_type>              base_type;

    public:

      typedef  typename bucket_type::associated_iterator_category  iterator_category;
      typedef  std::ptrdiff_t                                      difference_type;
      typedef  std::size_t                                         size_type;
      typedef  const typename bucket_type::value_type              value_type;
      typedef  value_type*                                         pointer;
      typedef  value_type&                                         reference;


      hash_table_const_iterator ()
      { }

      hash_table_const_iterator (const hash_table_data <bucket_type>& data, bucket_type* bucket)
        : base_type (data, bucket)
      { }


      value_type&
      operator* () const
      {
        MCT_PRECONDITION (this->valid (), "invalid iterator");
        return this->_bucket->value ();
      }

      value_type*
      operator-> () const
      {
        MCT_PRECONDITION (this->valid (), "invalid iterator");
        return &this->_bucket->value ();
      }


      hash_table_const_iterator&
      operator++ ()
      {
        this->move_forward ();
        return *this;
      }

      hash_table_const_iterator
      operator++ (int)
      {
        const hash_table_const_iterator  copy (*this);

        this->move_forward ();
        return copy;
      }

      hash_table_const_iterator&
      operator-- ()
      {
        this->move_backward ();
        return *this;
      }

      hash_table_const_iterator
      operator-- (int)
      {
        const hash_table_const_iterator  copy (*this);

        this->move_backward ();
        return copy;
      }
    };


    template <typename Bucket>
    class hash_table_iterator : public hash_table_const_iterator <Bucket>
    {
      typedef  Bucket                                   bucket_type;
      typedef  hash_table_const_iterator <bucket_type>  base_type;

    public:

      typedef  typename bucket_type::assignable_value_type  value_type;
      typedef  value_type*                                  pointer;
      typedef  value_type&                                  reference;


      hash_table_iterator ()
      { }

      hash_table_iterator (const hash_table_data <bucket_type>& data, bucket_type* bucket)
        : base_type (data, bucket)
      { }


      value_type&
      operator* () const
      {
        MCT_PRECONDITION (this->valid (), "invalid iterator");
        return this->_bucket->value ();
      }

      value_type*
      operator-> () const
      {
        MCT_PRECONDITION (this->valid (), "invalid iterator");
        return &this->_bucket->value ();
      }


      hash_table_iterator&
      operator++ ()
      {
        return static_cast <hash_table_iterator&> (base_type::operator++ ());
      }

      hash_table_iterator
      operator++ (int)
      {
        return static_cast <hash_table_iterator> (base_type::operator++ ());
      }

      hash_table_iterator&
      operator-- ()
      {
        return static_cast <hash_table_iterator&> (base_type::operator-- ());
      }

      hash_table_iterator
      operator-- (int)
      {
        return static_cast <hash_table_iterator> (base_type::operator-- ());
      }
    };


    template <typename Bucket>
    inline  bool
    operator== (const hash_table_const_iterator <Bucket>& iterator1,
                const hash_table_const_iterator <Bucket>& iterator2)
    {  return iterator1.bucket () == iterator2.bucket ();  }

    template <typename Bucket>
    inline  bool
    operator!= (const hash_table_const_iterator <Bucket>& iterator1,
                const hash_table_const_iterator <Bucket>& iterator2)
    {  return iterator1.bucket () != iterator2.bucket ();  }

    template <typename Bucket>
    inline  bool
    operator< (const hash_table_const_iterator <Bucket>& iterator1,
               const hash_table_const_iterator <Bucket>& iterator2)
    {  return iterator1.bucket () < iterator2.bucket ();  }

    template <typename Bucket>
    inline  bool
    operator<= (const hash_table_const_iterator <Bucket>& iterator1,
                const hash_table_const_iterator <Bucket>& iterator2)
    {  return iterator1.bucket () <= iterator2.bucket ();  }

    template <typename Bucket>
    inline  bool
    operator> (const hash_table_const_iterator <Bucket>& iterator1,
               const hash_table_const_iterator <Bucket>& iterator2)
    {  return iterator1.bucket () > iterator2.bucket ();  }

    template <typename Bucket>
    inline  bool
    operator>= (const hash_table_const_iterator <Bucket>& iterator1,
                const hash_table_const_iterator <Bucket>& iterator2)
    {  return iterator1.bucket () >= iterator2.bucket ();  }


    // These types are currently undocumented and should be considered package-private.
    // They are used in testing.
    struct closed_table_family
    { };

    struct linked_table_family
    { };

    struct forward_table_family
    { };


    const std::size_t  DEFAULT_NUM_BUCKETS     = 16;
    const std::size_t  MIN_NUM_BUCKETS         = 4;
    const float        DEFAULT_MAX_LOAD_FACTOR = 0.6f;
    const float        MIN_MAX_LOAD_FACTOR     = 0.01f;
    const float        MAX_MAX_LOAD_FACTOR     = 0.99f;


    template <typename Bucket, typename Hash, typename Equal, typename Allocator>
    class hash_table_base
    {
      template <typename Bucket1_, typename Bucket2_, typename Hash1_, typename Hash2_,
                typename Equal_, typename Allocator1_, typename Allocator2_>
      friend  typename impl::enable_if <impl::is_same <typename Bucket1_::value_type,
                                                       typename Bucket2_::value_type>::value,
                                        bool>::type
      operator== (const hash_table_base <Bucket1_, Hash1_, Equal_, Allocator1_>& table1,
                  const hash_table_base <Bucket2_, Hash2_, Equal_, Allocator2_>& table2);


    protected:

      typedef  Bucket                                                    bucket_type;
      typedef  typename bucket_type::usage_data                          bucket_usage_data;
      typedef  hash_table_data <bucket_type>                             data_type;

    public:

      typedef  typename bucket_type::key_type                            key_type;
      typedef  typename bucket_type::value_type                          value_type;
      typedef  Hash                                                      hasher;
      typedef  Equal                                                     key_equal;
      typedef  typename Allocator::template rebind <value_type>::other   allocator_type;

      typedef  hash_table_iterator <bucket_type>                         iterator;
      typedef  hash_table_const_iterator <bucket_type>                   const_iterator;

      typedef  typename allocator_type::pointer                          pointer;
      typedef  typename allocator_type::const_pointer                    const_pointer;
      typedef  typename allocator_type::reference                        reference;
      typedef  typename allocator_type::const_reference                  const_reference;
      typedef  std::ptrdiff_t                                            difference_type;
      typedef  std::size_t                                               size_type;


# if MCT_DEBUGGING_MEMBERS

      struct statistics
      {
        double     debris_ratio;
        double     avg_present_lookup;
        size_type  max_present_lookup;
        double     avg_absent_lookup;
        size_type  max_absent_lookup;
      };

# endif


    protected:

      typedef  typename allocator_type::template rebind <bucket_type>::other
               bucket_allocator_type;


      data_type       _data;
      size_type       _max_occupied;
      hasher          _hash;
      key_equal       _equal;
      allocator_type  _allocator;
      float           _max_load_factor;


      // This structure is used to simplify exception safety code.
      struct data_grip : public data_type
      {
        hash_table_base&  table;

        data_grip (hash_table_base& table, size_type num_buckets, bool allocate)
          : table (table)
        {
          if (allocate)
            table.allocate_buckets (*this, num_buckets);
          else
            table.postpone_bucket_allocation (*this, num_buckets);
        }

        ~data_grip ()
        {
          table.deallocate_buckets (*this);
        }
      };


    public:

      hash_table_base (size_type num_buckets, const hasher& hash, const key_equal& equal,
                       const allocator_type& allocator)
        : _hash      (hash),
          _equal     (equal),
          _allocator (allocator)
      {
        postpone_bucket_allocation (_data, finalize_num_buckets (num_buckets));
        do_set_max_load_factor (DEFAULT_MAX_LOAD_FACTOR);

        MCT_VALIDATION (validate_integrity ());
      }

      hash_table_base (const hash_table_base& that);
      hash_table_base (const hash_table_base& that, const allocator_type& allocator);

#   if MCT_CXX0X_SUPPORTED

      hash_table_base (hash_table_base&& that)
        : _data            (that._data),
          _max_occupied    (that._max_occupied),
          _hash            (that.hash_function ()),
          _equal           (that.key_eq ()),
          _allocator       (that.get_allocator ()),
          _max_load_factor (that.max_load_factor ())
      {
        that._data.buckets = 0;
        that._data.clear ();
        that._max_occupied = 0;

        MCT_VALIDATION (this->validate_integrity ());
        MCT_VALIDATION (that .validate_integrity ());
      }

      hash_table_base (std::initializer_list <value_type> initializer,
                       size_type num_buckets, const hasher& hash, const key_equal& equal,
                       const allocator_type& allocator)
        : _hash      (hash),
          _equal     (equal),
          _allocator (allocator)
      {
        initialize_from_range (initializer.begin (), initializer.end (), num_buckets);
      }

#   endif  // MCT_CXX0X_SUPPORTED

      template <typename InputIterator>
      hash_table_base (InputIterator first, InputIterator last,
                       size_type num_buckets, const hasher& hash, const key_equal& equal,
                       const allocator_type& allocator)
        : _hash      (hash),
          _equal     (equal),
          _allocator (allocator)
      {
        initialize_from_range (first, last, num_buckets);
      }

    private:

      template <typename InputIterator>
      void  initialize_from_range (InputIterator first, InputIterator last, size_type num_buckets);

      void  initialize_as_copy (const hash_table_base& that);

    public:

      ~hash_table_base ();


      hash_table_base&  operator= (const hash_table_base& that);

#   if MCT_CXX0X_SUPPORTED

      hash_table_base&  operator= (hash_table_base&& that);

      hash_table_base&
      operator= (std::initializer_list <value_type> initializer)
      {
        hash_table_base  temp (initializer, 0, _hash, _equal, _allocator);
        swap (temp);
        return *this;
      }

#   endif


      iterator
      begin ()
      {
        return make_iterator (_data.begin ());
      }

      const_iterator
      begin () const
      {
        return make_const_iterator (_data.begin ());
      }

      const_iterator
      cbegin () const
      {
        return make_const_iterator (_data.begin ());
      }

      iterator
      end ()
      {
        return make_iterator (_data.end ());
      }

      const_iterator
      end () const
      {
        return make_const_iterator (_data.end ());
      }

      const_iterator
      cend () const
      {
        return make_const_iterator (_data.end ());
      }


      bool
      empty () const
      {  return _data.num_used == 0;  }

      size_type
      size () const
      {  return _data.num_used;  }

      size_type
      max_size () const
      {
        return static_cast <size_type> (std::floor (max_bucket_count () * MAX_MAX_LOAD_FACTOR));
      }


      float
      load_factor () const
      {
        return static_cast <float> (_data.num_used) / _data.num_buckets;
      }

      float
      max_load_factor () const
      {  return _max_load_factor;  }

      void
      max_load_factor (float max_load_factor)
      {
        do_set_max_load_factor (max_load_factor);
        if (_data.num_occupied > _max_occupied || !_max_occupied)
          do_rehash ();
      }


      size_type
      bucket_count () const
      {  return _data.num_buckets;  }

      size_type
      max_bucket_count () const
      {
        // limit_max_count() is esentially std::min() of the value with a constant.
        // However, making that constant a static member variable of bucket types doesn't
        // seem to be possible, as we get undefined references errors then.
        return bucket_type::limit_max_count (round_down_to_power_of_2
                                             (bucket_allocator_type (_allocator).max_size ()
                                              - bucket_type::USES_SENTINEL));
      }


      hasher
      hash_function () const
      {  return _hash;  }

      key_equal
      key_eq () const
      {  return _equal;  }

      allocator_type
      get_allocator () const
      {  return _allocator;  }


      iterator
      find (const key_type& key)
      {
        return make_iterator (_data.num_used != 0 ? lookup (key) : _data.end ());
      }

      const_iterator
      find (const key_type& key) const
      {
        return make_const_iterator (_data.num_used != 0 ? lookup (key) : _data.end ());
      }

      size_type
      count (const key_type& key) const
      {
        return _data.num_used != 0 && lookup (key) != _data.end () ? 1 : 0;
      }

      std::pair <iterator, iterator>
      equal_range (const key_type& key)
      {
        const iterator  first (make_iterator (_data.num_used != 0 ? lookup (key) : _data.end ()));
        iterator        last  (first);

        if (first.bucket () != _data.end ())
          ++last;

        return std::make_pair (first, last);
      }

      std::pair <const_iterator, const_iterator>
      equal_range (const key_type& key) const
      {
        const const_iterator  first (make_const_iterator (_data.num_used != 0
                                                          ? lookup (key) : _data.end ()));
        const_iterator        last  (first);

        if (first.bucket () != _data.end ())
          ++last;

        return std::make_pair (first, last);
      }


      void  clear ();
      void  swap (hash_table_base& that);

      void
      rehash (size_type num_buckets)
      {  do_rehash (num_buckets);  }

      void
      reserve (size_type num_elements)
      {
        do_rehash (static_cast <size_type> (std::ceil (num_elements / max_load_factor ())));
      }


#   if MCT_DEBUGGING_MEMBERS

      bool
      valid_iterator (const_iterator position) const
      {
        // In some unlikely cases certain iterators (e.g. from another table after several
        // rehashes) could point to a valid bucket, but still be invalid.
        return position.is_for (_data) && _data.valid_bucket (position.bucket ());
      }

      bool
      valid_iterator_range (const_iterator first, const_iterator last) const
      {
        const const_iterator  _end = end ();
        if (first != _end && !valid_iterator (first))
          return false;

        if (last != _end && !valid_iterator (last))
          return false;

        for (; first != _end; ++first)
          {
            if (first == last)
              return true;
          }

        return true;
      }

      size_type
      used_memory () const
      {
        // Note: this relies on subclasses' not having any fields.  This is tested though,
        // so any breakage will be noticed quickly.
        return (sizeof *this
                + (_data.buckets
                   ? (_data.num_buckets + bucket_type::USES_SENTINEL) * sizeof (bucket_type)
                   : 0));
      }

      void        validate_integrity () const;
      statistics  collect_statistics () const;

#   endif  // MCT_DEBUGGING_MEMBERS


    protected:

      iterator
      make_iterator (bucket_type* bucket) const
      {  return iterator (_data, bucket);  }

      const_iterator
      make_const_iterator (bucket_type* bucket) const
      {  return const_iterator (_data, bucket);  }


      bucket_type*  lookup (const key_type& key)  const;

      // Templated for convenience of map's operator[].  Exact meaning of 'that' is up to
      // subclasses or, more specifically, to the bucket type.  For linked tables this is
      // 'before', for forward table it's 'after'.
      template <typename type>
      std::pair <bucket_type*, bool>  lookup_or_insert (type data, bucket_type* that = 0);


      void  do_set_max_load_factor (float max_load_factor);


      bucket_type*
      rehash_for_insertion (size_type num_elements, bucket_type* watch_for = 0)
      {
        if (num_elements)
          {
            if (_data.num_used + num_elements > _max_occupied)
              {
                return do_rehash (static_cast <size_type> (std::ceil ((_data.num_used
                                                                       + num_elements)
                                                                      / max_load_factor ())),
                                  true, watch_for);
              }
            else if (_data.num_occupied + num_elements > _max_occupied)
              return clear_debris_or_grow (watch_for);
          }

        return watch_for;
      }

      bucket_type*  do_rehash (size_type num_buckets = 0, bool buckets_required = false,
                               bucket_type* watch_for = 0);
      bucket_type*  clear_debris_or_grow (bucket_type* watch_for = 0);


      size_type
      compute_max_occupied (size_type num_buckets) const
      {
        if (num_buckets >= max_bucket_count ())
          return max_size ();

        return std::min (num_buckets - 1,
                         static_cast <size_type> (std::floor (num_buckets * _max_load_factor)));
      }


    private:

      static  size_type
      finalize_num_buckets (size_type num_buckets)
      {
        if (num_buckets > MIN_NUM_BUCKETS)
          return round_up_to_power_of_2 (num_buckets);
        else
          return num_buckets ? MIN_NUM_BUCKETS : DEFAULT_NUM_BUCKETS;
      }

      void
      postpone_bucket_allocation (data_type& data, size_type num_buckets)
      {
        // Even if we are not allocating buckets now, from outside it should look more or
        // less the same.
        if (num_buckets > max_bucket_count ())
          throw std::length_error ("closed_hash_table::postpone_bucket_allocation");

        data.num_buckets = num_buckets;
        data.buckets     = 0;
        data.reset_back ();
      }

      void
      allocate_buckets (data_type& data, size_type num_buckets)
      {
        if (num_buckets > max_bucket_count ())
          throw std::length_error ("closed_hash_table::allocate_buckets");

        const size_type  real_num_buckets = (num_buckets + bucket_type::USES_SENTINEL);

        data.num_buckets = num_buckets;
        data.buckets     = bucket_allocator_type (_allocator).allocate (real_num_buckets);
        data.reset_back ();

        for (bucket_type* bucket = data.buckets, * const limit = data.storage_end ();
             bucket != limit; ++bucket)
          bucket->mark_clear ();

        if (bucket_type::USES_SENTINEL)
          data.sentinel ()->initialize_as_sentinel ();
      }

      void
      deallocate_buckets (const data_type& data)
      {
        // Empty tables can never allocate any buckets at all.  Additionally, move
        // constructor or assignment "steal" buckets.
        if (!data.buckets)
          return;

        // Note: we mustn't depend on exact value of 'num_used', because that might be
        // incorrect if we are cleaning up after an exception.  However, it will not be
        // zero if the table is or can be non-empty.
        if (!has_trivial_destructor <value_type>::value && data.num_used)
          {
            for (bucket_type* bucket = data.begin (), * end = data.end ();
                 bucket != end; bucket = bucket->next (data))
              _allocator.destroy (&bucket->value ());
          }

        const size_type  real_num_buckets = (data.num_buckets + bucket_type::USES_SENTINEL);
        bucket_allocator_type (_allocator).deallocate (data.buckets, real_num_buckets);
      }


      // This is really internal.  Destination bucket array must be completely clean and
      // have enough space to fit the used elements in 'from'.  'from' must have been
      // built using the same equality and hash functions as this table.  Returned value
      // is the bucket to where 'watch_for' is copied.
      bucket_type*
      copy_buckets (const data_type& from, data_type& to, bucket_type* watch_for = 0);
    };


    template <typename Bucket, typename Hash, typename Equal, typename Allocator>
    class closed_hash_table : public hash_table_base <Bucket, Hash, Equal, Allocator>
    {
      typedef  hash_table_base <Bucket, Hash, Equal, Allocator>  base_type;

    protected:

      typedef  typename base_type::bucket_type                   bucket_type;

    public:

      typedef  typename base_type::key_type            key_type;
      typedef  typename base_type::value_type          value_type;
      typedef  typename base_type::hasher              hasher;
      typedef  typename base_type::key_equal           key_equal;
      typedef  typename base_type::allocator_type      allocator_type;
      typedef  typename base_type::iterator            iterator;
      typedef  typename base_type::const_iterator      const_iterator;
      typedef  typename base_type::size_type           size_type;

      // Undocumented and package-private.
      typedef  closed_table_family                     _family;


      closed_hash_table (size_type num_buckets, const hasher& hash, const key_equal& equal,
                         const allocator_type& allocator)
        : base_type (num_buckets, hash, equal, allocator)
      { }

      closed_hash_table (const closed_hash_table& that)
        : base_type (that)
      { }

      closed_hash_table (const closed_hash_table& that, const allocator_type& allocator)
        : base_type (that, allocator)
      { }

#   if MCT_CXX0X_SUPPORTED

      closed_hash_table (closed_hash_table&& that)
        : base_type (std::move (that))
      { }

      closed_hash_table (std::initializer_list <value_type> initializer,
                         size_type num_buckets, const hasher& hash, const key_equal& equal,
                         const allocator_type& allocator)
        : base_type (initializer, num_buckets, hash, equal, allocator)
      { }

#   endif  // MCT_CXX0X_SUPPORTED

      template <typename InputIterator>
      closed_hash_table (InputIterator first, InputIterator last,
                         size_type num_buckets, const hasher& hash, const key_equal& equal,
                         const allocator_type& allocator)
        : base_type (first, last, num_buckets, hash, equal, allocator)
      { }


#   if MCT_CXX0X_SUPPORTED

      closed_hash_table&
      operator= (std::initializer_list <value_type> initializer)
      {
        return static_cast <closed_hash_table&> (base_type::operator= (initializer));
      }

#   endif


      std::pair <iterator, bool>
      insert (const value_type& value)
      {
        const std::pair <bucket_type*, bool>  result
          (this->template lookup_or_insert <const value_type&> (value));
        return std::make_pair (this->make_iterator (result.first), result.second);
      }

      // Note that we currently just ignore 'hint'.
      iterator
      insert (const_iterator hint, const value_type& value)
      {
        MCT_PRECONDITION (hint == this->end () || this->valid_iterator (hint), "invalid hint");
        return this->make_iterator (this->template lookup_or_insert <const value_type&>
                                    (value).first);
      }

      template <typename InputIterator>
      void
      insert (InputIterator first, InputIterator last)
      {
        this->rehash_for_insertion (forward_distance (first, last));
        for (; first != last; ++first)
          this->template lookup_or_insert <const value_type&> (*first);
      }

#   if MCT_CXX0X_SUPPORTED

      std::pair <iterator, bool>
      insert (value_type&& value)
      {
        const std::pair <bucket_type*, bool>  result
          (this->template lookup_or_insert <value_type&&> (std::move (value)));
        return std::make_pair (this->make_iterator (result.first), result.second);
      }

      iterator
      insert (const_iterator hint, value_type&& value)
      {
        MCT_PRECONDITION (hint == this->end () || this->valid_iterator (hint), "invalid hint");
        return this->make_iterator (this->template lookup_or_insert <value_type&&>
                                    (std::move (value)).first);
      }

      void
      insert (std::initializer_list <value_type> initializer)
      {
        insert (initializer.begin (), initializer.end ());
      }


      template <typename... Args>
      std::pair <iterator, bool>
      emplace (Args&&... args)
      {
        return insert (value_type (std::forward <Args> (args)...));
      }

      template <typename... Args>
      std::pair <iterator, bool>
      emplace_hint (const_iterator hint, Args&&... args)
      {
        MCT_PRECONDITION (hint == this->end () || this->valid_iterator (hint), "invalid hint");
        return insert (hint, value_type (std::forward <Args> (args)...));
      }

#   endif  // MCT_CXX0X_SUPPORTED


      size_type
      erase (const key_type& key)
      {
        return this->_data.num_used != 0 && do_erase (key) ? 1 : 0;
      }

      iterator
      erase (const_iterator position)
      {
        MCT_PRECONDITION (this->valid_iterator (position), "invalid position");
        bucket_type* const  next = position.bucket ()->next (this->_data);
        do_erase (position.bucket ());
        return this->make_iterator (next);
      }

      void
      quick_erase (const_iterator position)
      {
        MCT_PRECONDITION (this->valid_iterator (position), "invalid position");
        do_erase (position.bucket ());
      }

      iterator
      erase (const_iterator first, const_iterator last)
      {
        MCT_PRECONDITION (this->valid_iterator_range (first, last), "invalid range");

        bucket_type*  bucket = first.bucket ();
        for (bucket_type* const limit = last.bucket (); bucket != limit;)
          {
            bucket_type* const  next = bucket->next (this->_data);
            do_erase (bucket);
            bucket = next;
          }

        return this->make_iterator (bucket);
      }


    protected:

      bool
      do_erase (const key_type& key)
      {
        bucket_type* const  bucket = this->lookup (key);

        if (bucket != this->_data.end ())
          {
            do_erase (bucket);
            return true;
          }
        else
          return false;
      }

      void
      do_erase (bucket_type* bucket)
      {
        this->_data.note_being_destroyed (bucket);
        --this->_data.num_used;
        this->_allocator.destroy (&bucket->value ());
        bucket->template unlink <bucket_type> ();
        bucket->mark_destroyed ();

        MCT_VALIDATION (this->validate_integrity ());
      }
    };


    template <typename Bucket, typename Hash, typename Equal, typename Allocator>
    class linked_hash_table : public closed_hash_table <Bucket, Hash, Equal, Allocator>
    {
      typedef  closed_hash_table <Bucket, Hash, Equal, Allocator>  base_type;

    protected:

      typedef  typename base_type::bucket_type                     bucket_type;

    public:

      typedef  typename base_type::value_type          value_type;
      typedef  typename base_type::hasher              hasher;
      typedef  typename base_type::key_equal           key_equal;
      typedef  typename base_type::allocator_type      allocator_type;
      typedef  typename base_type::iterator            iterator;
      typedef  typename base_type::const_iterator      const_iterator;
      typedef  std::reverse_iterator <iterator>        reverse_iterator;
      typedef  std::reverse_iterator <const_iterator>  const_reverse_iterator;
      typedef  typename base_type::size_type           size_type;

      // Undocumented and package-private.
      typedef  linked_table_family                     _family;


      linked_hash_table (size_type num_buckets, const hasher& hash, const key_equal& equal,
                         const allocator_type& allocator)
        : base_type (num_buckets, hash, equal, allocator)
      { }

      linked_hash_table (const linked_hash_table& that)
        : base_type (that)
      { }

      linked_hash_table (const linked_hash_table& that, const allocator_type& allocator)
        : base_type (that, allocator)
      { }

#   if MCT_CXX0X_SUPPORTED

      linked_hash_table (linked_hash_table&& that)
        : base_type (std::move (that))
      { }

      // We don't use similar base constructor, because we need a non-default second
      // argument to lookup_or_insert().
      linked_hash_table (std::initializer_list <value_type> initializer,
                         size_type num_buckets, const hasher& hash, const key_equal& equal,
                         const allocator_type& allocator)
        : base_type ((num_buckets
                      ? num_buckets
                      : static_cast <size_type> (std::ceil (initializer.size ()
                                                            / DEFAULT_MAX_LOAD_FACTOR))),
                     hash, equal, allocator)
      {
        initialize_from_range (initializer.begin (), initializer.end ());
      }

#   endif  // MCT_CXX0X_SUPPORTED

      // We don't use similar base constructor, because we need a non-default second
      // argument to lookup_or_insert().
      template <typename InputIterator>
      linked_hash_table (InputIterator first, InputIterator last,
                         size_type num_buckets, const hasher& hash, const key_equal& equal,
                         const allocator_type& allocator)
        : base_type ((num_buckets
                      ? num_buckets
                      : forward_distance_scaled (first, last, DEFAULT_MAX_LOAD_FACTOR)),
                     hash, equal, allocator)
      {
        initialize_from_range (first, last);
      }


    private:

      template <typename InputIterator>
      void
      initialize_from_range (InputIterator first, InputIterator last)
      {
        for (; first != last; ++first)
          this->template lookup_or_insert <const value_type&> (*first, this->_data.sentinel ());

        MCT_VALIDATION (this->validate_integrity ());
      }


    public:

#   if MCT_CXX0X_SUPPORTED

      // Can't use the operator from the base class, since the corresponding constructor
      // is overriden.
      linked_hash_table&
      operator= (std::initializer_list <value_type> initializer)
      {
        linked_hash_table  temp (initializer, 0, this->_hash, this->_equal, this->_allocator);
        swap (temp);
        return *this;
      }

#   endif


      reverse_iterator
      rbegin ()
      {
        return reverse_iterator (this->end ());
      }

      const_reverse_iterator
      rbegin () const
      {
        return const_reverse_iterator (this->end ());
      }

      const_reverse_iterator
      crbegin () const
      {
        return const_reverse_iterator (this->cend ());
      }

      reverse_iterator
      rend ()
      {
        return reverse_iterator (this->begin ());
      }

      const_reverse_iterator
      rend () const
      {
        return const_reverse_iterator (this->begin ());
      }

      const_reverse_iterator
      crend () const
      {
        return const_reverse_iterator (this->cbegin ());
      }


      value_type&
      front ()
      {
        MCT_PRECONDITION (!this->empty (), "table is empty");
        return this->_data.begin ()->value ();
      }

      const value_type&
      front () const
      {
        MCT_PRECONDITION (!this->empty (), "table is empty");
        return this->_data.begin ()->value ();
      }

      value_type&
      back ()
      {
        MCT_PRECONDITION (!this->empty (), "table is empty");
        return this->_data.sentinel ()->template previous <bucket_type> ()->value ();
      }

      const value_type&
      back () const
      {
        MCT_PRECONDITION (!this->empty (), "table is empty");
        return this->_data.sentinel ()->template previous <bucket_type> ()->value ();
      }


      std::pair <iterator, bool>
      push_front (const value_type& value)
      {
        const std::pair <bucket_type*, bool>  result
          (this->template lookup_or_insert <const value_type&> (value, this->_data.begin ()));
        return std::make_pair (this->make_iterator (result.first), result.second);
      }

      void
      pop_front ()
      {
        MCT_PRECONDITION (!this->empty (), "table is empty");
        this->erase (this->make_iterator (this->_data.begin ()));
      }

      std::pair <iterator, bool>
      push_back (const value_type& value)
      {
        const std::pair <bucket_type*, bool>  result
          (this->template lookup_or_insert <const value_type&> (value, this->_data.sentinel ()));
        return std::make_pair (this->make_iterator (result.first), result.second);
      }

      void
      pop_back ()
      {
        MCT_PRECONDITION (!this->empty (), "table is empty");
        this->erase (this->make_iterator (this->_data.sentinel ()
                                          ->template previous <bucket_type> ()));
      }


#   if MCT_CXX0X_SUPPORTED

      std::pair <iterator, bool>
      push_front (value_type&& value)
      {
        const std::pair <bucket_type*, bool>  result
          (this->template lookup_or_insert <value_type&&> (std::move (value),
                                                           this->_data.begin ()));
        return std::make_pair (this->make_iterator (result.first), result.second);
      }

      std::pair <iterator, bool>
      push_back (value_type&& value)
      {
        const std::pair <bucket_type*, bool>  result
          (this->template lookup_or_insert <value_type&&> (std::move (value),
                                                           this->_data.sentinel ()));
        return std::make_pair (this->make_iterator (result.first), result.second);
      }

#   endif  // MCT_CXX0X_SUPPORTED


      // Need to override all 'insert' overloads, else lookup fails.  Additionally, in
      // 'linked_hash_table' lookup_or_insert() actually uses 'before' argument, so 0
      // default value is not acceptable.

      std::pair <iterator, bool>
      insert (const value_type& value)
      {
        const std::pair <bucket_type*, bool>  result
          (this->template lookup_or_insert <const value_type&> (value, this->_data.sentinel ()));
        return std::make_pair (this->make_iterator (result.first), result.second);
      }

      template <typename InputIterator>
      void
      insert (InputIterator first, InputIterator last)
      {  insert (this->end (), first, last);  }

      std::pair <iterator, bool>
      insert (const_iterator before, const value_type& value)
      {
        MCT_PRECONDITION (before == this->end () || this->valid_iterator (before),
                          "invalid position");
        const std::pair <bucket_type*, bool>  result
          (this->template lookup_or_insert <const value_type&> (value, before.bucket ()));
        return std::make_pair (this->make_iterator (result.first), result.second);
      }

      template <typename InputIterator>
      void  insert (const_iterator before, InputIterator first, InputIterator last);

#   if MCT_CXX0X_SUPPORTED

      std::pair <iterator, bool>
      insert (value_type&& value)
      {
        const std::pair <bucket_type*, bool>  result
          (this->template lookup_or_insert <value_type&&> (std::move (value),
                                                           this->_data.sentinel ()));
        return std::make_pair (this->make_iterator (result.first), result.second);
      }

      std::pair <iterator, bool>
      insert (const_iterator before, value_type&& value)
      {
        MCT_PRECONDITION (before == this->end () || this->valid_iterator (before),
                          "invalid position");
        const std::pair <bucket_type*, bool>  result
          (this->template lookup_or_insert <value_type&&> (std::move (value), before.bucket ()));
        return std::make_pair (this->make_iterator (result.first), result.second);
      }

      void
      insert (std::initializer_list <value_type> initializer)
      {
        insert (this->end (), initializer.begin (), initializer.end ());
      }

      void
      insert (const_iterator before, std::initializer_list <value_type> initializer)
      {
        insert (before, initializer.begin (), initializer.end ());
      }


      template <typename... Args>
      std::pair <iterator, bool>
      emplace (Args&&... args)
      {
        return insert (value_type (std::forward <Args> (args)...));
      }

      template <typename... Args>
      std::pair <iterator, bool>
      emplace_before (const_iterator before, Args&&... args)
      {
        return insert (before, value_type (std::forward <Args> (args)...));
      }

      template <typename... Args>
      std::pair <iterator, bool>
      emplace_hint (const_iterator hint, Args&&... args)
      {
        MCT_PRECONDITION (hint == this->end () || this->valid_iterator (hint), "invalid hint");
        return insert (hint, value_type (std::forward <Args> (args)...));
      }

#   endif  // MCT_CXX0X_SUPPORTED


      void
      relink (const_iterator before, iterator element)
      {
        MCT_PRECONDITION (before == this->end () || this->valid_iterator (before),
                          "invalid position");
        MCT_PRECONDITION (this->valid_iterator (element), "invalid element");

        if (element.bucket () != before.bucket ())
          {
            element.bucket ()->template unlink <bucket_type> ();
            element.bucket ()->link (before.bucket ());

            MCT_VALIDATION (this->validate_integrity ());
          }
      }

      void  reverse ();
    };


    template <typename Bucket, typename Hash, typename Equal, typename Allocator>
    class forward_hash_table : public hash_table_base <Bucket, Hash, Equal, Allocator>
    {
      typedef  hash_table_base <Bucket, Hash, Equal, Allocator>  base_type;

    protected:

      typedef  typename base_type::bucket_type                   bucket_type;

    public:

      typedef  typename base_type::key_type            key_type;
      typedef  typename base_type::value_type          value_type;
      typedef  typename base_type::hasher              hasher;
      typedef  typename base_type::key_equal           key_equal;
      typedef  typename base_type::allocator_type      allocator_type;
      typedef  typename base_type::iterator            iterator;
      typedef  typename base_type::const_iterator      const_iterator;
      typedef  std::reverse_iterator <iterator>        reverse_iterator;
      typedef  std::reverse_iterator <const_iterator>  const_reverse_iterator;
      typedef  typename base_type::size_type           size_type;

      // Undocumented and package-private.
      typedef  forward_table_family                    _family;


      forward_hash_table (size_type num_buckets, const hasher& hash, const key_equal& equal,
                          const allocator_type& allocator)
        : base_type (num_buckets, hash, equal, allocator)
      { }

      forward_hash_table (const forward_hash_table& that)
        : base_type (that)
      { }

      forward_hash_table (const forward_hash_table& that, const allocator_type& allocator)
        : base_type (that, allocator)
      { }

#   if MCT_CXX0X_SUPPORTED

      forward_hash_table (forward_hash_table&& that)
        : base_type (std::move (that))
      { }

      // We don't use similar base constructor, because we need a non-default second
      // argument to lookup_or_insert().
      forward_hash_table (std::initializer_list <value_type> initializer,
                          size_type num_buckets, const hasher& hash, const key_equal& equal,
                          const allocator_type& allocator)
        : base_type ((num_buckets
                      ? num_buckets
                      : static_cast <size_type> (std::ceil (initializer.size ()
                                                            / DEFAULT_MAX_LOAD_FACTOR))),
                     hash, equal, allocator)
      {
        initialize_from_range (initializer.begin (), initializer.end ());
      }

#   endif  // MCT_CXX0X_SUPPORTED

      // We don't use similar base constructor, because we need a non-default second
      // argument to lookup_or_insert().
      template <typename InputIterator>
      forward_hash_table (InputIterator first, InputIterator last,
                          size_type num_buckets, const hasher& hash, const key_equal& equal,
                          const allocator_type& allocator)
        : base_type ((num_buckets
                      ? num_buckets
                      : forward_distance_scaled (first, last, DEFAULT_MAX_LOAD_FACTOR)),
                     hash, equal, allocator)
      {
        initialize_from_range (first, last);
      }


    private:

      template <typename InputIterator>
      void
      initialize_from_range (InputIterator first, InputIterator last)
      {
        for (; first != last; ++first)
          this->template lookup_or_insert <const value_type&> (*first, this->_data.back);

        MCT_VALIDATION (this->validate_integrity ());
      }


    public:

#   if MCT_CXX0X_SUPPORTED

      // Can't use the operator from the base class, since the corresponding constructor
      // is overriden.
      forward_hash_table&
      operator= (std::initializer_list <value_type> initializer)
      {
        forward_hash_table  temp (initializer, 0, this->_hash, this->_equal, this->_allocator);
        swap (temp);
        return *this;
      }

#   endif


      iterator
      before_begin ()
      {
        return this->make_iterator (this->_data.sentinel ());
      }

      const_iterator
      before_begin () const
      {
        return this->make_const_iterator (this->_data.sentinel ());
      }

      const_iterator
      cbefore_begin () const
      {
        return this->make_const_iterator (this->_data.sentinel ());
      }


      iterator
      before_end ()
      {
        return this->make_iterator (this->_data.back);
      }

      const_iterator
      before_end () const
      {
        return this->make_const_iterator (this->_data.back);
      }

      const_iterator
      cbefore_end () const
      {
        return this->make_const_iterator (this->_data.back);
      }


      value_type&
      front ()
      {
        MCT_PRECONDITION (!this->empty (), "table is empty");
        return this->_data.begin ()->value ();
      }

      const value_type&
      front () const
      {
        MCT_PRECONDITION (!this->empty (), "table is empty");
        return this->_data.begin ()->value ();
      }

      value_type&
      back ()
      {
        MCT_PRECONDITION (!this->empty (), "table is empty");
        return this->_data.back->value ();
      }

      const value_type&
      back () const
      {
        MCT_PRECONDITION (!this->empty (), "table is empty");
        return this->_data.back->value ();
      }


      std::pair <iterator, bool>
      push_front (const value_type& value)
      {
        const std::pair <bucket_type*, bool>  result
          (this->template lookup_or_insert <const value_type&> (value, this->_data.sentinel ()));
        return std::make_pair (this->make_iterator (result.first), result.second);
      }

      void
      pop_front ()
      {
        MCT_PRECONDITION (!this->empty (), "table is empty");
        do_erase_after (this->_data.sentinel ());
      }

      std::pair <iterator, bool>
      push_back (const value_type& value)
      {
        const std::pair <bucket_type*, bool>  result
          (this->template lookup_or_insert <const value_type&> (value, this->_data.back));
        return std::make_pair (this->make_iterator (result.first), result.second);
      }


#   if MCT_CXX0X_SUPPORTED

      std::pair <iterator, bool>
      push_front (value_type&& value)
      {
        const std::pair <bucket_type*, bool>  result
          (this->template lookup_or_insert <value_type&&> (std::move (value),
                                                           this->_data.sentinel ()));
        return std::make_pair (this->make_iterator (result.first), result.second);
      }

      std::pair <iterator, bool>
      push_back (value_type&& value)
      {
        const std::pair <bucket_type*, bool>  result
          (this->template lookup_or_insert <value_type&&> (std::move (value), this->_data.back));
        return std::make_pair (this->make_iterator (result.first), result.second);
      }

#   endif  // MCT_CXX0X_SUPPORTED


      std::pair <iterator, bool>
      insert (const value_type& value)
      {
        const std::pair <bucket_type*, bool>  result
          (this->template lookup_or_insert <const value_type&> (value, this->_data.back));
        return std::make_pair (this->make_iterator (result.first), result.second);
      }

      template <typename InputIterator>
      void
      insert (InputIterator first, InputIterator last)
      {  insert_after (before_end (), first, last);  }

      std::pair <iterator, bool>
      insert_after (const_iterator after, const value_type& value)
      {
        MCT_PRECONDITION (after == this->before_begin () || this->valid_iterator (after),
                          "invalid position");
        const std::pair <bucket_type*, bool>  result
          (this->template lookup_or_insert <const value_type&> (value, after.bucket ()));
        return std::make_pair (this->make_iterator (result.first), result.second);
      }

      template <typename InputIterator>
      void  insert_after (const_iterator after, InputIterator first, InputIterator last);

#   if MCT_CXX0X_SUPPORTED

      std::pair <iterator, bool>
      insert (value_type&& value)
      {
        const std::pair <bucket_type*, bool>  result
          (this->template lookup_or_insert <value_type&&> (std::move (value), this->_data.back));
        return std::make_pair (this->make_iterator (result.first), result.second);
      }

      void
      insert (std::initializer_list <value_type> initializer)
      {
        insert_after (before_end (), initializer.begin (), initializer.end ());
      }


      template <typename... Args>
      std::pair <iterator, bool>
      emplace (Args&&... args)
      {
        return insert (value_type (std::forward <Args> (args)...));
      }

      template <typename... Args>
      std::pair <iterator, bool>
      emplace_after (const_iterator after, Args&&... args)
      {
        return insert_after (after, value_type (std::forward <Args> (args)...));
      }

#   endif  // MCT_CXX0X_SUPPORTED


      void
      erase_after (const_iterator position)
      {
        MCT_PRECONDITION ((position == before_begin ()
                           || (this->valid_iterator (position) && position != before_end ())),
                          "invalid position");
        do_erase_after (position.bucket ());
      }

      void
      erase_after (const_iterator after, const_iterator last)
      {
#     if MCT_CHECK_PRECONDITIONS
        MCT_PRECONDITION (after == before_begin () || this->valid_iterator (after),
                          "invalid 'after'");

        const_iterator  first (after);
        MCT_PRECONDITION (this->valid_iterator_range (++first, last), "invalid range");
#     endif

        for (bucket_type* const after_bucket = after.bucket (), * last_bucket = last.bucket ();
             after_bucket->template next <bucket_type> () != last_bucket;)
          do_erase_after (after_bucket);
      }


    protected:

      void
      do_erase_after (bucket_type* after)
      {
        bucket_type*  bucket = after->next (this->_data);

        this->_data.note_being_destroyed (bucket, after);
        --this->_data.num_used;
        this->_allocator.destroy (&bucket->value ());
        bucket->unlink (after);
        bucket->mark_destroyed ();

        MCT_VALIDATION (this->validate_integrity ());
      }


    public:

      void
      relink_after (const_iterator to_after, const_iterator from_after)
      {
        MCT_PRECONDITION (to_after == this->before_begin () || this->valid_iterator (to_after),
                          "invalid to_after");
        MCT_PRECONDITION ((from_after == this->before_begin ()
                           || (this->valid_iterator (from_after)
                               && from_after != this->before_end ())),
                          "invalid from_after");

        bucket_type* const  bucket = from_after.bucket ()->template next <bucket_type> ();

        if (bucket != to_after.bucket ())
          {
            this->_data.note_being_unlinked (bucket, from_after.bucket ());
            bucket->unlink (static_cast <bucket_type*> (from_after.bucket ()));
            bucket->link   (to_after.bucket ());
            this->_data.note_linked (bucket);
          }

        MCT_VALIDATION (this->validate_integrity ());
      }

      void  reverse ();


#   if MCT_DEBUGGING_MEMBERS

      void
      validate_integrity () const
      {
        base_type::validate_integrity ();

        bucket_type* const  back = this->_data.back;
        bool                back_is_valid;

        if (this->empty ())
          back_is_valid = (back == this->_data.storage_end ());
        else
          back_is_valid = this->_data.valid_bucket (back);

        if (back_is_valid && this->_data.buckets
            && back->template next <bucket_type> () != this->_data.end ())
          back_is_valid = false;

        if (!back_is_valid)
          throw std::logic_error ("back reference is not valid");
      }

#   endif
    };


    template <typename Bucket, typename Hash, typename Equal, typename Allocator>
    hash_table_base <Bucket, Hash, Equal, Allocator>::
    hash_table_base (const hash_table_base& that)
      : _hash            (that.hash_function ()),
        _equal           (that.key_eq ()),
        _allocator       (that.get_allocator ()),
        _max_load_factor (that.max_load_factor ())
    {
      initialize_as_copy (that);
    }

    template <typename Bucket, typename Hash, typename Equal, typename Allocator>
    hash_table_base <Bucket, Hash, Equal, Allocator>::
    hash_table_base (const hash_table_base& that, const allocator_type& allocator)
      : _hash            (that.hash_function ()),
        _equal           (that.key_eq ()),
        _allocator       (allocator),
        _max_load_factor (that.max_load_factor ())
    {
      initialize_as_copy (that);
    }


    template <typename Bucket, typename Hash, typename Equal, typename Allocator>
    template <typename InputIterator>
    void
    hash_table_base <Bucket, Hash, Equal, Allocator>::
    initialize_from_range (InputIterator first, InputIterator last, size_type num_buckets)
    {
      if (!num_buckets)
        num_buckets = forward_distance_scaled (first, last, DEFAULT_MAX_LOAD_FACTOR);

      num_buckets = finalize_num_buckets (num_buckets);

      if (first == last)
        {
          postpone_bucket_allocation (_data, num_buckets);
          do_set_max_load_factor (DEFAULT_MAX_LOAD_FACTOR);
        }
      else
        {
          allocate_buckets (_data, num_buckets);
          do_set_max_load_factor (DEFAULT_MAX_LOAD_FACTOR);

          try
            {
              for (; first != last; ++first)
                lookup_or_insert <const value_type&> (*first);
            }
          catch (...)
            {
              deallocate_buckets (_data);
              throw;
            }
        }

      MCT_VALIDATION (validate_integrity ());
    }

    template <typename Bucket, typename Hash, typename Equal, typename Allocator>
    inline  void
    hash_table_base <Bucket, Hash, Equal, Allocator>::
    initialize_as_copy (const hash_table_base& that)
    {
      size_type  num_buckets
        = finalize_num_buckets (static_cast <size_type> (std::ceil (that._data.num_used
                                                                    / max_load_factor ())));
      size_type  max_occupied = compute_max_occupied (num_buckets);

      // Can happen when copying tiny tables with very low maximum load factor.
      while (max_occupied == 0)
        {
          num_buckets  <<= 1;
          max_occupied   = compute_max_occupied (num_buckets);
        }

      if (!that._data.num_used)
        {
          postpone_bucket_allocation (_data, num_buckets);
          _max_occupied = 0;
        }
      else
        {
          allocate_buckets (_data, num_buckets);
          _max_occupied = max_occupied;

          try
            {
              copy_buckets (that._data, _data);
            }
          catch (...)
            {
              deallocate_buckets (_data);
              throw;
            }
        }

      MCT_VALIDATION (validate_integrity ());
    }


    template <typename Bucket, typename Hash, typename Equal, typename Allocator>
    hash_table_base <Bucket, Hash, Equal, Allocator>::
    ~hash_table_base ()
    {
      deallocate_buckets (_data);
    }


    template <typename Bucket, typename Hash, typename Equal, typename Allocator>
    hash_table_base <Bucket, Hash, Equal, Allocator>&
    hash_table_base <Bucket, Hash, Equal, Allocator>::
    operator= (const hash_table_base& that)
    {
      // Just an optimization; will work without this check too.
      if (&that != this)
        {
          hash_table_base  copy (that, _allocator);
          swap (copy);
        }

      return *this;
    }

# if MCT_CXX0X_SUPPORTED

    template <typename Bucket, typename Hash, typename Equal, typename Allocator>
    hash_table_base <Bucket, Hash, Equal, Allocator>&
    hash_table_base <Bucket, Hash, Equal, Allocator>::
    operator= (hash_table_base&& that)
    {
      // Apparently, there is no agreement whether 'x = std::move (x)' should leave 'x'
      // unchanged or not, at least in GCC.  However, without this check MCT will
      // misbehave in certain cases with GCC 4.6 std::stable_sort(), for example.  (Same
      // issue with std::vector, but so what if GNU stdlib is not bug-free?)
      if (&that != this)
        {
          if (this->_allocator == that._allocator)
            {
              deallocate_buckets (_data);

              this->_data            = that._data;
              this->_max_occupied    = that._max_occupied;
              this->_hash            = that._hash;
              this->_equal           = that._equal;
              this->_max_load_factor = that._max_load_factor;

              that._data.buckets = 0;
              that._data.clear ();
              that._max_occupied = 0;

              MCT_VALIDATION (this->validate_integrity ());
              MCT_VALIDATION (that .validate_integrity ());
            }
          else
            {
              hash_table_base  copy (that, _allocator);
              swap (copy);
            }
        }

      return *this;
    }

# endif  // MCT_CXX0X_SUPPORTED


    template <typename Bucket, typename Hash, typename Equal, typename Allocator>
    void
    hash_table_base <Bucket, Hash, Equal, Allocator>::
    do_set_max_load_factor (float max_load_factor)
    {
      _max_load_factor = std::max (MIN_MAX_LOAD_FACTOR,
                                   std::min (max_load_factor, MAX_MAX_LOAD_FACTOR));
      _max_occupied    = (_data.buckets ? compute_max_occupied (_data.num_buckets) : 0);
    }


    // According to benchmarks, explicitly unrolling first pass gives a small speedup for
    // lookup() and lookup_or_insert().

    template <typename Bucket, typename Hash, typename Equal, typename Allocator>
    typename hash_table_base <Bucket, Hash, Equal, Allocator>::bucket_type*
    hash_table_base <Bucket, Hash, Equal, Allocator>::
    lookup (const key_type& key) const
    {
      const size_type  hash    = bucket_type::postprocess_hash (_hash (key));
      size_type        look_at = _data.start_probing (hash);

      bucket_type* const       bucket = (_data.buckets + look_at);
      const bucket_usage_data  usage  = bucket->get_usage_data ();

      if (bucket_type::used_and_matches (usage, hash) || _data.is_special (bucket))
        {
          if (_equal (bucket->key (), key))
            return bucket;
        }
      else if (bucket_type::is_empty (usage))
        return _data.end ();

      for (size_type iteration = 1; ; ++iteration)
        {
          look_at = _data.continue_probing (look_at, iteration);

          bucket_type* const       bucket = (_data.buckets + look_at);
          const bucket_usage_data  usage  = bucket->get_usage_data ();

          if (bucket_type::used_and_matches (usage, hash) || _data.is_special (bucket))
            {
              if (_equal (bucket->key (), key))
                return bucket;
            }
          else if (bucket_type::is_empty (usage))
            return _data.end ();
        }
    }


    template <typename Bucket, typename Hash, typename Equal, typename Allocator>
    template <typename type>
    std::pair <typename hash_table_base <Bucket, Hash, Equal, Allocator>::bucket_type*, bool>
    hash_table_base <Bucket, Hash, Equal, Allocator>::
    lookup_or_insert (type data, bucket_type* that)
    {
      // A little bit pessimistic: if we end up replacing debris, we didn't need to
      // resize.  However, it is simpler to resize now (at most we "lose" one bucket this
      // way); besides, this way we don't need to check for a not allocated bucket array
      // separately.
      if (MCT_OPTIMIZATION_UNLIKELY (_data.num_occupied >= _max_occupied))
        that = clear_debris_or_grow (that);

      const key_type&  key     = bucket_type::extract_key (data);
      const size_type  hash    = bucket_type::postprocess_hash (_hash (key));
      size_type        look_at = _data.start_probing (hash);

      bucket_type* const       bucket = (_data.buckets + look_at);
      const bucket_usage_data  usage  = bucket->get_usage_data ();
      bucket_type*             debris;

      if (bucket_type::used_and_matches (usage, hash) || _data.is_special (bucket))
        {
          if (_equal (bucket->key (), key))
            return std::make_pair (bucket, false);

          debris = 0;
        }
      else if (bucket_type::is_empty (usage))
        {
          try
            {
              bucket_type::construct_value (_allocator, bucket->value (),
                                            MCT_STD_FORWARD (type, data));
            }
          catch (...)
            {
              bucket->revert_to_empty ();
              throw;
            }

          bucket->mark_used ();
          bucket->store_hash (hash);
          // At least on GCC 4.4 this helps optimization.
          if (bucket_type::USES_LINKS)
            bucket->link (that);
          _data.note_new (bucket);

          ++_data.num_used;
          ++_data.num_occupied;

          MCT_VALIDATION (validate_integrity ());
          return std::make_pair (bucket, true);
        }
      else
        debris = (!bucket_type::KEEPS_HASHES || bucket_type::is_debris (usage) ? bucket : 0);

      for (size_type iteration = 1; ; ++iteration)
        {
          look_at = _data.continue_probing (look_at, iteration);

          bucket_type*             bucket = (_data.buckets + look_at);
          const bucket_usage_data  usage  = bucket->get_usage_data ();

          if (bucket_type::used_and_matches (usage, hash) || _data.is_special (bucket))
            {
              if (_equal (bucket->key (), key))
                return std::make_pair (bucket, false);
            }
          else if (bucket_type::is_empty (usage))
            {
              if (debris)
                bucket = debris;

              try
                {
                  bucket_type::construct_value (_allocator, bucket->value (),
                                                MCT_STD_FORWARD (type, data));
                }
              catch (...)
                {
                  if (debris)
                    bucket->revert_to_debris ();
                  else
                    bucket->revert_to_empty ();

                  throw;
                }

              bucket->mark_used ();
              bucket->store_hash (hash);
              if (bucket_type::USES_LINKS)
                bucket->link (that);
              _data.note_new (bucket);

              ++_data.num_used;
              if (!debris)
                ++_data.num_occupied;

              MCT_VALIDATION (validate_integrity ());
              return std::make_pair (bucket, true);
            }
          else if (!debris && (!bucket_type::KEEPS_HASHES || bucket_type::is_debris (usage)))
            debris = bucket;
        }
    }


    template <typename Bucket, typename Hash, typename Equal, typename Allocator>
    void
    hash_table_base <Bucket, Hash, Equal, Allocator>::
    clear ()
    {
      if (!_data.num_occupied)
        return;

      for (bucket_type* bucket = _data.buckets, * const limit = _data.storage_end ();
           bucket != limit; ++bucket)
        {
          if (!has_trivial_destructor <value_type>::value)
            {
              if (!bucket_type::is_empty_or_debris (bucket->get_usage_data ())
                  || _data.is_special (bucket))
                _allocator.destroy (&bucket->value ());
            }

          bucket->mark_clear ();
        }

      _data.clear ();
      if (bucket_type::USES_SENTINEL)
        _data.sentinel ()->initialize_as_sentinel ();

      MCT_VALIDATION (validate_integrity ());
    }


    template <typename Bucket, typename Hash, typename Equal, typename Allocator>
    void
    hash_table_base <Bucket, Hash, Equal, Allocator>::
    swap (hash_table_base& that)
    {
      using std::swap;

      if (this->_allocator == that._allocator)
        _data.swap (that._data);
      else
        {
          // 'new_this_data' must be allocated with 'this' allocator, but data must come
          // from 'that' and also copied by 'that' (because hash function and equality
          // predicate are not swapped yet).
          data_grip  new_this_data (*this, that ._data.num_buckets, true);
          data_grip  new_that_data (that,  this->_data.num_buckets, true);

          that .copy_buckets (that. _data, new_this_data);
          this->copy_buckets (this->_data, new_that_data);

          this->_data.swap (new_this_data);
          that ._data.swap (new_that_data);
        }

      swap (_max_occupied,    that._max_occupied);
      swap (_hash,            that._hash);
      swap (_equal,           that._equal);
      swap (_max_load_factor, that._max_load_factor);

      MCT_VALIDATION (this->validate_integrity ());
      MCT_VALIDATION (that. validate_integrity ());
    }


    template <typename Bucket, typename Hash, typename Equal, typename Allocator>
    typename hash_table_base <Bucket, Hash, Equal, Allocator>::bucket_type*
    hash_table_base <Bucket, Hash, Equal, Allocator>::
    do_rehash (size_type num_buckets, bool buckets_required, bucket_type* watch_for)
    {
      // In all cases we need at least one empty bucket: important for tiny tables with
      // large load factor limit.
      num_buckets = std::max (std::max (std::max (num_buckets, _data.num_used + 1),
                                        MIN_NUM_BUCKETS),
                              static_cast <size_type> (std::ceil (_data.num_used
                                                                  / max_load_factor ())));
      num_buckets = round_up_to_power_of_2 (num_buckets);

      size_type  max_occupied = compute_max_occupied (num_buckets);

      // Can happen for tiny tables with very low maximum load factor.
      while (max_occupied == 0)
        {
          num_buckets  <<= 1;
          max_occupied   = compute_max_occupied (num_buckets);
        }

      if (buckets_required || _data.num_used)
        {
          // Only skip rehashing altogether if number of buckets stays the same _and_
          // there are no debris.  Otherwise, at least get rid of debris.
          if (num_buckets == _data.num_buckets
              && _data.buckets
              && _data.num_occupied == _data.num_used)
            return watch_for;

          data_grip  new_data (*this, num_buckets, true);

          watch_for = copy_buckets (_data, new_data, watch_for);
          _data.swap (new_data);
          _max_occupied = max_occupied;
        }
      else
        {
          // When rehashing an empty table, just drop memory and be done with it.  No point in
          // allocating new memory array right away, we can as well postpone it in the hope it
          // will eventually not be needed.
          data_grip  empty_data (*this, num_buckets, false);

          _data.swap (empty_data);
          _max_occupied = 0;
        }

      MCT_VALIDATION (validate_integrity ());
      return watch_for;
    }

    template <typename Bucket, typename Hash, typename Equal, typename Allocator>
    typename hash_table_base <Bucket, Hash, Equal, Allocator>::bucket_type*
    hash_table_base <Bucket, Hash, Equal, Allocator>::
    clear_debris_or_grow (bucket_type* watch_for)
    {
      size_type  num_buckets = _data.num_buckets;

      // Bit shifting trick to avoid overflow if table is *huge*.
      if (_data.buckets && _data.num_used >= (_data.num_occupied - (_data.num_occupied >> 4)))
        {
          // If almost all table is filled with real data rather than debris, we increase
          // table size instead, if we still have space to increase.  However, if we are
          // already at maximum size _and_ there are no debris, still try to increase size
          // so that do_rehash() call below throws.
          if (num_buckets != max_bucket_count () || _data.num_occupied == _data.num_used)
            num_buckets <<= 1;
        }

      return do_rehash (num_buckets, true, watch_for);
    }

    template <typename Bucket, typename Hash, typename Equal, typename Allocator>
    typename hash_table_base <Bucket, Hash, Equal, Allocator>::bucket_type*
    hash_table_base <Bucket, Hash, Equal, Allocator>::
    copy_buckets (const data_type& from, data_type& to, bucket_type* watch_for)
    {
      // If copying is interrupted with exception we needn't have correct values here, but
      // we need to at least have non-zero values to bust possible is-empty optimizations.
      to.num_used     = from.num_used;
      to.num_occupied = from.num_used;

      if (from.num_used)
        {
          for (bucket_type* bucket = from.begin (), * limit = from.end (), * to_end = to.end ();
               bucket != limit; bucket = bucket->next (from))
            {
              const size_type  hash    = (bucket_type::KEEPS_HASHES
                                          ? bucket->hash () : _hash (bucket->key ()));
              size_type        look_at = to.start_probing (hash);

              bucket_type*  new_bucket = (to.buckets + look_at);

              if (!bucket_type::is_empty (new_bucket->get_usage_data ())
                  || to.is_special (new_bucket))
                {
                  for (size_type iteration = 1; ; ++iteration)
                    {
                      look_at = to.continue_probing (look_at, iteration);

                      new_bucket = (to.buckets + look_at);
                      if (bucket_type::is_empty (new_bucket->get_usage_data ())
                          && !to.is_special (new_bucket))
                        break;
                    }
                }

              try
                {
                  bucket_type::construct_value (_allocator, new_bucket->value (),
                                                bucket->value ());
                }
              catch (...)
                {
                  new_bucket->revert_to_empty ();
                  throw;
                }

              new_bucket->mark_used ();
              new_bucket->copy_hash (*bucket);
              new_bucket->link (!bucket_type::IS_FORWARD ? to_end : to.get_back ());
              to.note_new (new_bucket);

              if (watch_for == bucket)
                watch_for = new_bucket;
            }
        }

      if (watch_for == from.storage_end ())
        watch_for = to.storage_end ();

      return watch_for;
    }


# if MCT_DEBUGGING_MEMBERS

    template <typename Bucket, typename Hash, typename Equal, typename Allocator>
    void
    hash_table_base <Bucket, Hash, Equal, Allocator>::
    validate_integrity () const
    {
      if (_data.num_buckets < MIN_NUM_BUCKETS)
        throw std::logic_error ("there are too few buckets in a table");

      if (_data.num_buckets > max_bucket_count ())
        throw std::logic_error ("there are too many buckets in a table");

      if (_data.buckets)
        {
          if (!_max_occupied)
            throw std::logic_error ("occupation threshold of an allocated table is zero");

          if (_max_occupied != compute_max_occupied (_data.num_buckets))
            throw std::logic_error ("occupation threshold is incorrect");
        }
      else
        {
          if (_max_occupied)
            throw std::logic_error ("occupation threshold of an unallocated table is nonzero");
        }

      if (_data.num_occupied > _max_occupied)
        throw std::logic_error ("too many occupied buckets");

      size_type  num_used     = 0;
      size_type  num_occupied = 0;

      if (_data.buckets)
        {
          for (bucket_type* bucket = _data.buckets, * const _end = _data.storage_end ();
               bucket != _end; ++bucket)
            {
              if (!_data.is_special (bucket))
                {
                  const bucket_usage_data  usage = bucket->get_usage_data ();
                  if (bucket_type::is_empty (usage))
                    continue;

                  ++num_occupied;
                  if (bucket_type::is_debris (usage))
                    continue;
                }
              else
                ++num_occupied;

              ++num_used;
              if (lookup (bucket->key ()) != bucket)
                throw std::logic_error ("bucket lookup fails");
            }
        }

      if (_data.num_used != num_used)
        throw std::logic_error ("wrong number of used buckets");

      if (_data.num_occupied != num_occupied)
        throw std::logic_error ("wrong number of occupied buckets");

      if (_data.buckets)
        bucket_type::validate_links (_data.sentinel (), num_used);
    }

    template <typename Bucket, typename Hash, typename Equal, typename Allocator>
    typename hash_table_base <Bucket, Hash, Equal, Allocator>::statistics
    hash_table_base <Bucket, Hash, Equal, Allocator>::
    collect_statistics () const
    {
      statistics  stats;

      stats.debris_ratio       = 0;
      stats.avg_present_lookup = 0;
      stats.max_present_lookup = 0;
      stats.avg_absent_lookup  = 0;
      stats.max_absent_lookup  = 0;

      if (_data.num_occupied != 0)
        {
          stats.debris_ratio = (static_cast <double> (_data.num_occupied - _data.num_used)
                                / _data.num_occupied);
        }

      if (!_data.num_used)
        {
          stats.avg_absent_lookup = 1;
          stats.max_absent_lookup = 1;
        }
      else
        {
          for (bucket_type* bucket = _data.buckets, * const _end = _data.storage_end ();
               bucket != _end; ++bucket)
            {
              if (!bucket_type::is_empty_or_debris (bucket->get_usage_data ())
                  || _data.is_special (bucket))
                {
                  const size_type  hash    = _hash (bucket->key ());
                  size_type        look_at = _data.start_probing (hash);

                  for (size_type iteration = 1; ; ++iteration)
                    {
                      if (_data.buckets + look_at == bucket)
                        {
                          stats.avg_present_lookup += iteration;
                          stats.max_present_lookup  = std::max (stats.max_present_lookup,
                                                                iteration);
                          break;
                        }

                      look_at = _data.continue_probing (look_at, iteration);
                    }
                }

              size_type  look_at = (bucket - _data.buckets);

              for (size_type iteration = 1; ; ++iteration)
                {
                  if (bucket_type::is_empty ((_data.buckets + look_at)->get_usage_data ()))
                    {
                      stats.avg_absent_lookup += iteration;
                      stats.max_absent_lookup  = std::max (stats.max_absent_lookup, iteration);
                      break;
                    }

                  look_at = _data.continue_probing (look_at, iteration);
                }
            }

          if (_data.num_used != 0)
            stats.avg_present_lookup /= _data.num_used;

          stats.avg_absent_lookup /= _data.num_buckets;
        }

      return stats;
    }

# endif  // MCT_DEBUGGING_MEMBERS


    template <typename Bucket, typename Hash, typename Equal, typename Allocator>
    template <typename InputIterator>
    void
    linked_hash_table <Bucket, Hash, Equal, Allocator>::
    insert (const_iterator before, InputIterator first, InputIterator last)
    {
      MCT_PRECONDITION (before == this->end () || this->valid_iterator (before),
                        "invalid position");

      bucket_type*  before_bucket = this->rehash_for_insertion (forward_distance (first, last),
                                                                before.bucket ());
      size_type     num_inserted  = 0;

      // For linked tables it is possible to provide strong exception safety.
      // Additionally, we need cannot keep using 'before' as it may be invalidated due to
      // a rehash, which may happen at any iteration, especially for real input iterators.
      try
        {
          for (; first != last; ++first)
            {
              std::pair <bucket_type*, bool>  result
                (this->template lookup_or_insert <const value_type&> (*first, before_bucket));

              if (result.second)
                {
                  ++num_inserted;
                  before_bucket = result.first->template next <bucket_type> ();
                }
            }
        }
      catch (...)
        {
          const_iterator  rollback_last  (this->make_const_iterator (before_bucket));
          const_iterator  rollback_first (rollback_last);

          while (num_inserted--)
            --rollback_first;

          this->erase (rollback_first, rollback_last);
          throw;
        }
    }


    template <typename Bucket, typename Hash, typename Equal, typename Allocator>
    void
    linked_hash_table <Bucket, Hash, Equal, Allocator>::
    reverse ()
    {
      for (bucket_type* bucket = this->_data.sentinel (), * limit = bucket; ;)
        {
          bucket->reverse_links ();
          bucket = bucket->template next <bucket_type> ();
          if (bucket == limit)
            break;
        }

      MCT_VALIDATION (this->validate_integrity ());
    }


    // Currently no strong exception safety, unlike in 'linked_hash_table's range insert.
    // The reason is that rollback would need to have original 'after', while insertion
    // needs 'after' advanced after each next element and lookup_or_insert() has only one
    // 'watch_for' slot.
    template <typename Bucket, typename Hash, typename Equal, typename Allocator>
    template <typename InputIterator>
    void
    forward_hash_table <Bucket, Hash, Equal, Allocator>::
    insert_after (const_iterator after, InputIterator first, InputIterator last)
    {
      MCT_PRECONDITION (after == this->before_begin () || this->valid_iterator (after),
                        "invalid position");

      bucket_type*  after_bucket  = this->rehash_for_insertion (forward_distance (first, last),
                                                                after.bucket ());

      for (; first != last; ++first)
        {
          std::pair <bucket_type*, bool>  result
            (this->template lookup_or_insert <const value_type&> (*first, after_bucket));

          if (result.second)
            after_bucket = result.first;
        }
    }


    template <typename Bucket, typename Hash, typename Equal, typename Allocator>
    void
    forward_hash_table <Bucket, Hash, Equal, Allocator>::
    reverse ()
    {
      if (this->_data.num_used != 0)
        {
          bucket_type*  end      = this->_data.end ();
          bucket_type*  previous = end;
          bucket_type*  bucket   = this->_data.begin ();

          this->_data.back = bucket;

          for (bucket_type* next; bucket != end; bucket = next)
            {
              next = bucket->next (this->_data);
              bucket->set_next (previous);
              previous = bucket;
            }

          this->_data.sentinel ()->set_next (previous);
        }

      MCT_VALIDATION (this->validate_integrity ());
    }


    // Some SFINAE madness.  Since operator== must have two arguments, 'enable_if' has to
    // be used on the return type.
    template <typename Bucket1, typename Bucket2, typename Hash1, typename Hash2,
              typename Equal, typename Allocator1, typename Allocator2>
    typename impl::enable_if <impl::is_same <typename Bucket1::value_type,
                                             typename Bucket2::value_type>::value,
                              bool>::type
    operator== (const hash_table_base <Bucket1, Hash1, Equal, Allocator1>& table1,
                const hash_table_base <Bucket2, Hash2, Equal, Allocator2>& table2)
    {
      typedef  hash_table_base <Bucket1, Hash1, Equal, Allocator1>  table1_type;
      typedef  hash_table_base <Bucket2, Hash2, Equal, Allocator2>  table2_type;
      typedef  typename table1_type::size_type                      size_type;
      typedef  typename table2_type::bucket_type                    bucket1_type;
      typedef  typename table2_type::bucket_type                    bucket2_type;

      // The cast is only needed to silence compilers since otherwise it appears as if we
      // compare pointers to different types.
      if (is_same <typename Bucket1::value_type, typename Bucket2::value_type>::value
          && static_cast <const void*> (&table1) == static_cast <const void*> (&table2))
        return true;

      const size_type  size = table1.size ();
      if (table2.size () != size)
        return false;

      if (size > 0)
        {
          const bucket2_type* const             end = table2._data.end ();
          typename table1_type::const_iterator  iterator (table1.begin ());

          for (size_type k = 0; k < size; ++k, ++iterator)
            {
              const bucket2_type*  match = table2.lookup (iterator.bucket ()->key ());
              if (match == end
                  || !bucket1_type::mapped_equal (match->value (), iterator.bucket ()->value ()))
                return false;
            }
        }

      return true;
    }

    template <typename Bucket1, typename Bucket2, typename Hash, typename Equal,
              typename Allocator1, typename Allocator2>
    inline  typename impl::enable_if <impl::is_same <typename Bucket1::value_type,
                                                     typename Bucket2::value_type>::value,
                                      bool>::type
    operator!= (const hash_table_base <Bucket1, Hash, Equal, Allocator1>& table1,
                const hash_table_base <Bucket2, Hash, Equal, Allocator2>& table2)
    {
      return !operator== (table1, table2);
    }

  }

}


#undef MCT_STD_FORWARD


#endif  // Multi-inclusion guard.


// Local variables:
// mode: c++
// c-basic-offset: 2
// indent-tabs-mode: nil
// fill-column: 90
// End:
