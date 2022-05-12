#ifndef _THOR_SCSI_LEGACY_IO_H_
#define _THOR_SCSI_LEGACY_IO_H_ 1

#include <iostream>
#include <string>

#define NameLength 150  // Max length of identifiers (e.g. file names).

typedef char partsName[NameLength];

// Pascal file I/O (legacy).
extern int P_eof(FILE *f);

/* Check if at end of line (or end of entire file). */
extern int P_eoln(FILE *f);


// C++ file I/O.

void file_rd(std::ifstream &inf, const std::string &file_name);

void file_wr(std::ofstream &outf, const std::string &file_name);

void file_rd(std::ifstream &inf, const char file_name[]);

void file_wr(std::ofstream &outf, const char file_name[]);


// C file I/O.
FILE* file_read(const char file_name[]);

FILE* file_write(const char file_name[]);


#endif /*  _THOR_SCSI_LEGACY_IO_H_ */
/*
 * Local Variables:
 * mode: c++
 * c-file-style: "python"
 * End:
 */
