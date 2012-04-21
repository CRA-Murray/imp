/**
 *  \file IMP/rmf/links.h
 *  \brief Handle read/write of Model data from/to files.
 *
 *  Copyright 2007-2012 IMP Inventors. All rights reserved.
 *
 */

#ifndef IMPRMF_LINKS_H
#define IMPRMF_LINKS_H

#include "rmf_config.h"
#include <RMF/NodeHandle.h>
#include <RMF/FileHandle.h>
#include <IMP/base/Object.h>
#include <IMP/base/Pointer.h>
#include <IMP/base/object_macros.h>
#include <IMP/base/log_macros.h>

IMPRMF_BEGIN_NAMESPACE


#ifndef IMP_DOXYGEN
IMPRMFEXPORT
unsigned int get_load_linker_index(std::string st);
IMPRMFEXPORT
unsigned int get_save_linker_index(std::string st);
#endif

class SaveLink;
class LoadLink;
IMP_OBJECTS(SaveLink, SaveLinks);

IMP_OBJECTS(LoadLink, LoadLinks);

/** Manage a link between an \imp object and a part
    of the RMF file. This allows conformations to be
    loaded flexibly.

    LoadLinks must not save any handles to RMF objects.
*/
class IMPRMFEXPORT LoadLink: public base::Object {
 protected:
  virtual void do_load(RMF::FileConstHandle fh, unsigned int frame)=0;
  LoadLink(std::string name);
 public:
  void load(RMF::FileConstHandle fh, unsigned int frame) {
    IMP_OBJECT_LOG;
    set_was_used(true);
    do_load(fh, frame);
  }
};
/** Manage a link between an \imp object and a part
    of the RMF file. This allows conformations to be
    saved flexibly.

    SaveLinks must not save any handles to RMF objects.
*/
class IMPRMFEXPORT SaveLink: public base::Object {
 protected:
  virtual void do_save(RMF::FileHandle hf, unsigned int frame)=0;
  SaveLink(std::string name);
 public:
  void save(RMF::FileHandle fh, unsigned int frame) {
    IMP_OBJECT_LOG;
    set_was_used(true);
    do_save(fh, frame);
  }
};


IMPRMF_END_NAMESPACE

#endif /* IMPRMF_LINKS_H */
