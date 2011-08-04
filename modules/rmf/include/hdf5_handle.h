/**
 *  \file IMP/rmf/hdf5_handle.h
 *  \brief Handle read/write of Model data from/to files.
 *
 *  Copyright 2007-2011 IMP Inventors. All rights reserved.
 *
 */

#ifndef IMPRMF_HDF_5_HANDLE_H
#define IMPRMF_HDF_5_HANDLE_H

#include "rmf_config.h"
#include "NodeID.h"
#include "internal/utility.h"
#include "internal/intrusive_ptr_object.h"
#include <hdf5.h>
#include <algorithm>
#include <vector>
#include <exception>

#include <IMP/exception.h>
#include <IMP/RefCounted.h>
#include <limits>
#include <boost/utility.hpp>


#ifdef SWIG
typedef int hid_t;
#endif

IMPRMF_BEGIN_NAMESPACE
#ifndef SWIG
//! The signature for the HDF5 close functions
typedef herr_t (*HDF5CloseFunction)(hid_t) ;

//! Make sure an HDF5 handle is released
/** CloseFunction should be an appropriate close function
    for the handle type, eg H5Aclose. HDF5Handle is not available
    in python.
*/
class IMPRMFEXPORT HDF5Handle : public boost::noncopyable
 {
  hid_t h_;
  HDF5CloseFunction f_;
public:
  HDF5Handle(hid_t h, HDF5CloseFunction f): h_(h), f_(f) {
    if (h_<0) {
      IMP_RMF_THROW("Invalid handle returned", ValueException);
    }
  }
  HDF5Handle(): h_(-1){}
  hid_t get_hid() const {
    IMP_RMF_USAGE_CHECK(h_>=0, "Uninitialized handle used.");
    return h_;
  }
  operator hid_t() const {
    return h_;
  }
  bool get_is_open() const {
    return h_ != -1;
  }
  void open(hid_t h, HDF5CloseFunction f) {
    if (get_is_open()) {
      close();
    }
    h_=h;
    IMP_RMF_USAGE_CHECK(h_>=0, "Invalid handle returned");
    f_=f;
  }
  void close() {
    if (h_ != -1) {
      IMP_HDF5_CALL(f_(h_));
    }
    h_=-1;
  }
  ~HDF5Handle() {
    close();
  }
};

//! Share an HDF5 handle
/** This should be used with a boost intrusive_ptr. It is not available
    in python.
 */
class IMPRMFEXPORT HDF5SharedHandle: public HDF5Handle,
                                     public boost::intrusive_ptr_object

 {
public:
  HDF5SharedHandle(hid_t h, HDF5CloseFunction f): HDF5Handle(h, f) {
  }
 };

#endif // SWIG
IMPRMF_END_NAMESPACE

#endif /* IMPRMF_HDF_5_HANDLE_H */
