#include <Rcpp.h>
#include <iostream>
#include <fstream>
#include <string>
#include <vector>
#include <algorithm> // for copy
#include <iterator> // for ostream_iterator
#include <functional> 
#include <cctype>
#include <locale>
using namespace Rcpp;

/********************
 * Trimming strings.
 * Why isn't this in STL?
 * It's in Boost, but then we depend on 
 * the BH package. We're not quite ready
 * for that.
 * Ctrl-C+V from https://stackoverflow.com/questions/216823/whats-the-best-way-to-trim-stdstring
 */
// trim from start (in place)
static inline void ltrim(std::string &s) {
  s.erase(s.begin(), std::find_if(s.begin(), s.end(),
                  std::not1(std::ptr_fun<int, int>(std::isspace))));
}

// trim from end (in place)
static inline void rtrim(std::string &s) {
  s.erase(std::find_if(s.rbegin(), s.rend(),
                       std::not1(std::ptr_fun<int, int>(std::isspace))).base(), s.end());
}

// trim from both ends (in place)
static inline void trim(std::string &s) {
  ltrim(s);
  rtrim(s);
}


/*****************
 * cbind SNP files.
 */

// cbind snps

//' Column bind files with genotype matrices
//'
//' Does row-wise concatenation of multiple genotype files, 
//' ignoring the first column in subsequent files.
//' \strong{NB!} Assumes rows are ordered identically in all files, but will check that
//' value in first column is identical. If not, corresponding line is skipped in all files.
//' 
//' @param fnin Character vector of filenames to concatenate. 
//' @param fnout Filename of resulting file.
//' @param skiplines Integer, number of lines to skip before outputting.
//' @param excludeids Integer vector of first column id's to \emph{exclude} from the output. 
//' @inheritParams convert_phases
//' @return Number of lines written.
//' @rdname cbind_snp_files
//' @backref src/bind_snps.cpp
//' @seealso \code{\link{rbind_snp_files}}
//' @export
// [[Rcpp::export]]
int cbind_snp_files(CharacterVector fnin,
                    std::string fnout,
                    int skiplines = 0,
                    IntegerVector excludeids = IntegerVector::create(),
                    int idwidth = 4,
                    int precision = -1) 
{
  std::vector <std::ifstream *> infiles;
  int cnt = 0; // counter
  int skipped = 0;
  std::vector<long> myexcluded;

  if (excludeids.size() > 0) {
    excludeids.sort();  
    myexcluded.reserve(excludeids.size());
    for (unsigned int i = 0; i < excludeids.size(); i++) {
      myexcluded.push_back(excludeids[i]);
    }
  }
  
  infiles.reserve(fnin.size());
  
  // Open streams for input files.
  // Beware *not* to use implicit new (ifstream f(fn.c_str())),
  // as that will be destroyed after the containing block ({...}).
   
  for (unsigned int i = 0; i < fnin.size(); i++) {
    std::string fn = Rcpp::as<std::string>(fnin[i]);
    std::ifstream *f = new std::ifstream(fn.c_str());
    if (!f->is_open())
       Rcpp::stop("Error opening %s.", fn);
    infiles.push_back(f);
  }
  
  // Open stream for output file.
  std::ofstream fout(fnout.c_str());
  
  // Iterate through files and concatenate
  {
    bool ok=true;
    while (ok) {
      long id1;
      long id2;
      bool lineok=true;
      std::ostringstream longline;
      
      ++cnt;
      for (unsigned int i = 0; i < infiles.size(); i++) {
        std::ifstream *f = infiles[i];
        std::string line;
        if (!std::getline( *f, line) && i == 0) {
          ok = false;
          break;
        }
        if (!lineok)
          continue;
        
        if (cnt <= skiplines) {
          lineok = false;
          ++skipped;
          continue;
        }
        
        std::istringstream iss(line);
        iss >> id2;
        if (i == 0) {
          id1 = id2;
          longline.width(idwidth);
          longline.setf(std::ios::left, std::ios::adjustfield);
          longline << id1;
          if (excludeids.size() > 0 && std::binary_search(myexcluded.begin(), myexcluded.end(), id1)) {
            lineok = false;
            ++skipped;
            continue;
          }
        } else if (id1 != id2) {
          Rcpp::warning("Mismatching IDs on line %d between files 1 and %d (ids %s and %s). Skipping this line in all files.", cnt, i+1, id1, id2);
          lineok = false;
          ++skipped;
          continue;
        }
        
        std::getline(iss, line);
        trim(line);
        longline << " " << line;
          
      }
      if (lineok && ok) {
        fout << longline.str() << std::endl;
      }
    }
  }
  
  // Cleanup  
  {
    for (unsigned int i = 0; i < infiles.size(); i++) {
      std::ifstream *f = infiles[i];
      f->close();
      //delete *f;
    }
  }
  
  infiles.clear();
  fout.close();
  
  return cnt - skipped - 1;
}

