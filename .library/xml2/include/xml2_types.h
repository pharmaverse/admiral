#ifndef __XML2_XML2_TYPES__
#define __XML2_XML2_TYPES__

#include <libxml/tree.h>
#include <Rcpp.h>

inline void finaliseNode(xmlNodePtr node) {
  // do nothing
}

inline void finaliseNs(xmlNsPtr ns) {
  // do nothing
}

typedef Rcpp::XPtr<xmlDoc,Rcpp::PreserveStorage,xmlFreeDoc> XPtrDoc;
typedef Rcpp::XPtr<xmlNode,Rcpp::PreserveStorage,finaliseNode> XPtrNode;
typedef Rcpp::XPtr<xmlNs,Rcpp::PreserveStorage,finaliseNs> XPtrNs;
#endif
