#include <iostream>
#include <fstream>
#include <algorithm>
#include <cmath>

#include "Common.hpp"
#include "BUSData.h"

#include "bustools_whitelist.h"

#include "roaring.hh"

#define ERROR_RATE 0.01

void bustools_whitelist(Bustools_opt &opt) {
  BUSHeader h;
  size_t nr = 0;
  size_t N = 100000;
  BUSData *p = new BUSData[N];

  std::ofstream of(opt.output);
  std::ostream o(of.rdbuf());
  std::streambuf *inbuf;
  std::ifstream inf;
  if (!opt.stream_in) {
    inf.open(opt.files[0].c_str(), std::ios::binary);
    inbuf = inf.rdbuf();
  } else {
    inbuf = std::cin.rdbuf();
  }
  std::istream in(inbuf);
  parseHeader(in, h);

  uint32_t bclen = h.bclen;
  size_t rc = 1; // Non-zero so that second while loop works when using custom threshold
  int threshold;
  int wl_count = 0;
  
  uint32_t bc_r = 0, bc_u = 0;
  uint64_t curr_umi;
  uint64_t curr_bc;
  int bc_count = -1;
  std::vector<uint64_t> whitelist; // Only used when --nodist1

 /* Determine threshold. */ 
  if (opt.threshold) { // Custom threshold
    threshold = opt.threshold;
  } else { // Determine threshold from BUS data
    /* Get counts for all barcodes in first >=200 barcodes. */
    std::vector<wl_Record> vec;

    while (true) { 
      in.read((char*) p, N * sizeof(BUSData));
      rc = in.gcount() / sizeof(BUSData);
      if (rc == 0) {
        break;
      }
      nr += rc;

      for (size_t i = 0; i < rc; i++) {
        if (curr_bc != p[i].barcode) {
          if (bc_count != -1) {
            vec.emplace_back(curr_bc, bc_r, bc_u, bc_count);
          }
          curr_bc = p[i].barcode;
          bc_r = 1;
          bc_u = 1;
          bc_count = p[i].count;
        } else {
          ++bc_r;
          if (curr_umi != p[i].UMI) {
            if (bc_u == -1) {
              bc_u = 1;
              curr_umi = p[i].UMI;
            } else {
              ++bc_u;
            }
          }
          bc_count += p[i].count;
        }
      }
      /* Done going through BUSdata *p. */

      if (bc_count != -1) {
        vec.emplace_back(curr_bc, bc_r, bc_u, bc_count);
      }

      if (vec.size() >= 200) {
        break;
      }
    }
    /* Done retrieving first 200 barcodes. */
    // Note that the last-seen barcode may not have been fully processed by this point

    /* Sort. */
    std::sort(vec.begin(), vec.end(), [&](const wl_Record &a, const wl_Record &b) {
          if (a.count == b.count) {
            return a.barcode < b.barcode;
          } else {
            return a.count > b.count;
          }
        }
    );

    /* Determine threshold. */
    int M = 10; // Use first 10 barcodes
    int avgCount = 0, avgR = 0, avgU = 0;
    for (int i = 0; i < M; ++i) {
      avgCount += vec[i].count;
    }
    avgCount /= M;
    // [average count of top 10] * [chance of perfect barcode]
    // = [expected number of perfect barcodes]
    // And then multiply by some constant(?)
    threshold = avgCount * (1 - std::pow(1 - ERROR_RATE, bclen));
  
    /* Process all the records we just went through. */
    for (const auto &rec : vec) {
      if (rec.count >= threshold) {
	if (opt.nodist1) {
          whitelist.push_back(rec.barcode);
	} else {
          o << binaryToString(rec.barcode, bclen) << "\n";
          ++wl_count;
	}
      }
    }
  }
  /* Done determining threshold. */


  /* Go through remainder of records. */
  while (rc) {
    in.read((char*) p, N * sizeof(BUSData));
    rc = in.gcount() / sizeof(BUSData);
    if (rc == 0) {
      break;
    }
    nr += rc;

    for (size_t i = 0; i < rc; i++) {
      if (curr_bc != p[i].barcode || bc_count == -1) {
        if (bc_count >= threshold) {
          if (opt.nodist1) {
            whitelist.push_back(curr_bc);
          } else {
            o << binaryToString(curr_bc, bclen) << "\n";
            ++wl_count;
          }
        }
        bc_count = p[i].count;
        curr_bc = p[i].barcode;
      } else {
        bc_count += p[i].count;
      }
    }
    /* Done going through BUSdata *p. */

  }
  /* Done reading BUS file. */
  
  if (bc_count >= threshold) {
    if (opt.nodist1) {
      whitelist.push_back(curr_bc);
    } else {
      o << binaryToString(curr_bc, bclen) << "\n";
      ++wl_count;
    }
  }

  /* Remove hamming distance 1 barcodes */
  if (opt.nodist1) {
    // Hamming distance code copied from bustools_correct.cpp
    size_t bc2 = (bclen+1)/2;
    uint64_t mask_size = (1ULL << (2*bc2));
    uint64_t lower_mask = (1ULL<<(2*bc2))-1;
    uint64_t upper_mask = (1ULL<<(2*(bclen-bc2)))-1;

    std::vector<std::pair<Roaring,Roaring>> correct(1ULL<<(2*bc2)); // 4^(bc/2) possible barcodes
    for (uint64_t b : whitelist) {
      // Check whether barcode already exists (hamming distance 1)
      uint64_t lb = b & lower_mask;
      uint64_t ub = (b>>(2*bc2)) & upper_mask;   
      uint64_t lbc=0,ubc=0;
      int correct_lower = search_for_mismatch(correct[ub].second,bc2,lb,lbc);
      int correct_upper = search_for_mismatch(correct[lb].first,bclen - bc2,ub,ubc);
      int nc = correct_lower + correct_upper;

      // Barcode already exists
      if (nc == 1) continue;

      // Doesn't already exist; add
      lb = b & lower_mask;
      ub = (b>>(2*bc2)) & upper_mask;

      correct[ub].second.add(lb);
      correct[lb].first.add(ub);

      o << binaryToString(b, bclen) << "\n";
      ++wl_count;
    }
  }

  delete[] p; p = nullptr;
  of.close();
  std::cerr << "Read in " << nr << " BUS records, wrote " << wl_count << " barcodes to whitelist with threshold " << threshold << std::endl;
}
