#include <thor_scsi/exceptions.h>
#include <thor_scsi/legacy/io.h>
#include <iostream>
#include <fstream>

// Pascal file I/O (legacy).

int P_eof(FILE *f)
{
	int ch;

	if (feof(f)) return 1;
	if (f == stdin) return 0; /* not safe to look-ahead on the keyboard! */
	ch = getc(f);
	if (ch == EOF) return 1;
	ungetc(ch, f);

	return 0;
}


/* Check if at end of line (or end of entire file). */

int P_eoln(FILE *f)
{
	int ch;

	ch = getc(f);
	if (ch == EOF) return 1;
	ungetc(ch, f);
	return (ch == '\n');
}


// C++ file I/O.

void file_rd(std::ifstream &inf, const std::string &file_name)
{
	inf.open(file_name.c_str(), std::ios::in);
	if (!inf.is_open()) {
		std::cout << "File not found: " << file_name << "\n";
		throw std::ios_base::failure("File not found");
	}
}


void file_wr(std::ofstream &outf, const std::string &file_name)
{
	outf.open(file_name.c_str(), std::ios::out);
	if (!outf.is_open()) {
		std::cout << "Could not create file: " << file_name << "\n";
		throw std::ios_base::failure("Could not create file");
	}
}


void file_rd(std::ifstream &inf, const char file_name[])
{
	inf.open(file_name, std::ios::in);
	if (!inf.is_open()) {
		std::cerr << "File not found: >" << file_name << "<" << std::endl;
		throw std::ios_base::failure("File not found");
	}
}


void file_wr(std::ofstream &outf, const char file_name[])
{
	outf.open(file_name, std::ios::out);
	if (!outf.is_open()) {
		std::cerr << "Could not create file: >"  <<  file_name << "< " << std::endl;
		throw std::ios_base::failure("Could not create file");
	}
}


// C file I/O.

FILE* file_read(const char file_name[])
{
	FILE *fp;
	fp = fopen(file_name, "r");
	if (fp == NULL) {
		std::cerr << "File not found: >" << file_name << "<" << std::endl;
		throw std::ios_base::failure("File not found");
		// exit_(-1);
	}
	return(fp);
}


FILE* file_write(const char file_name[])
{
	FILE *fp;
	fp = fopen(file_name, "w");
	if (fp == NULL) {
		std::cerr << "Could not create file: >"  <<  file_name << "< " << std::endl;
		throw std::ios_base::failure("Could not create file");
		// exit_(-1);
	}
	return(fp);
}

/*
 * Local Variables:
 * mode: c++
 * c-file-style: "python"
 * End:
 */
