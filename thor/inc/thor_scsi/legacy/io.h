#ifndef _THOR_SCSI_LEGACY_IO_H_
#define _THOR_SCSI_LEGACY_IO_H_ 1

#include <iostream>
#include <string>

void file_rd(std::ifstream &inf, const std::string &file_name);

void file_wr(std::ofstream &outf, const std::string &file_name);

void file_rd(std::ifstream &inf, const char file_name[]);

void file_wr(std::ofstream &outf, const char file_name[]);

FILE* file_read(const char file_name[]);

FILE* file_write(const char file_name[]);

extern int P_eof(FILE *f);
extern int P_eoln(FILE *f);

#endif /*  _THOR_SCSI_LEGACY_IO_H_ */
