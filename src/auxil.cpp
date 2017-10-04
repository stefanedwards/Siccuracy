#include <Rcpp.h>
#include <iostream>
#include <fstream>
#include <string>
using namespace Rcpp;

// Gets number of rows in as file.
//'
//' \code{get_nlines} returns number of lines. Last line, if empty, is ignored.
//' 
//' @param fn Filename.
//' @export
//' @return \code{get_nlines}: Number of lines, or \code{0} on error.
//' @rdname auxfunc
//' @aliases get_nlines
//' @name get_nlines
//' @include auxil.R
// [[Rcpp::export]]
long get_nlines(std::string fn)
{
  
  long cnt;
  char oper;
  
  
  std::ifstream myfile(fn.c_str()); //Opens the file. c_str is mandatory here so that ifstream accepts the string path
  
  cnt = 0;
  
  if (myfile.is_open()) {
    while (myfile.get(oper)) {
      //Rcout << ++cnt << " " << oper << std::endl;
      ++cnt;
      if (oper != '\n')  // Enclose newline in single qoutes, else it is 2 characters long and not newline.
        myfile.ignore(std::numeric_limits<std::streamsize>::max(), '\n');
    }
  } else {
    Rcpp::stop("Unexpected error while reading file.");
  }
  
  myfile.close();

  //myfile.ignore(std::numeric_limits<std::streamsize>::max(), '\n')

  return cnt;
}