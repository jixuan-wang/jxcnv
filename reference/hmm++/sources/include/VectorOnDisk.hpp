#ifndef __VECTOR_ON_DISK_H__
#define __VECTOR_ON_DISK_H__

#include "params.hpp"
#include "Exception.hpp"
#include "LAPACKvector.hpp"

#include <string.h>
#include <sys/types.h>
#include <sys/stat.h>
#include <fcntl.h>
#include <unistd.h>
#include <sys/mman.h>

#include <sstream>
#include <iostream>
#include <cstdlib>
using namespace std;

namespace HMM_PP {

	/*
	 * Implementation follows ideas at:
	 * http://www.linuxquestions.org/questions/programming-9/mmap-tutorial-c-c-511265/
	 */
	template<class dataType>
	class VectorOnDisk : public LAPACKvector<dataType> {

	public:
		VectorOnDisk(const ullint size, const string& dir) {
		    //_fd = open(_file.c_str(), O_RDWR | O_CREAT | O_TRUNC, S_IRUSR | S_IWUSR);

			stringstream templateStr;
			templateStr << dir << "/" << "TEMP.XXXXXX";
			char* templateFile = strdup(templateStr.str().c_str());

			_fd = mkstemp(templateFile);
		    if (_fd == -1) {
		    	stringstream str;
		    	str << "Error opening temporary file of template '" << templateFile << "'";
		    	throw new Exception(str.str());
		    }
		    _file = templateFile;
		    free(templateFile);

		    // Stretch the file size to the size of the (mmapped) array of ints:
		    _fileSize = size * sizeof(dataType);
		    ullint result = lseek(_fd, _fileSize-1, SEEK_SET);
		    if (static_cast<int>(result) == -1) {
		    	closeFile();
		    	stringstream str;
		    	str << "Error calling lseek() to 'stretch' file " << _file;
		    	throw new Exception(str.str());
		    }

		    /* Something needs to be written at the end of the file to
		     * have the file actually have the new size.
		     * Just writing an empty string at the current file position will do.
		     *
		     * Note:
		     *  - The current position in the file is at the end of the stretched
		     *    file due to the call to lseek().
		     *  - An empty string is actually a single '\0' character, so a zero-byte
		     *    will be written at the last byte of the file.
		     */
		    result = write(_fd, "", 1);
		    if (result != 1) {
		    	closeFile();
		    	stringstream str;
		    	str << "Error writing last byte of file " << _file;
		    	throw new Exception(str.str());
		    }

		    // Now the file is ready to be mmapped.
		    _vec = static_cast<dataType*>(mmap(0, _fileSize, PROT_READ | PROT_WRITE, MAP_SHARED, _fd, 0));
		    if (_vec == MAP_FAILED) {
		    	closeFile();
		    	stringstream str;
		    	str << "Error mmapping file " << _file;
		    	throw new Exception(str.str());
		    }
		}

		virtual ~VectorOnDisk() {
		    bool failedUnmap = munmap(_vec, _fileSize) == -1;
		    closeFile();

		    if (failedUnmap) {
		    	stringstream str;
		    	str << "Error un-mmapping file " << _file;
		    	throw new Exception(str.str());
		    }
		}

		void closeFile() {
			close(_fd);

			if (unlink(_file.c_str()) == -1) {
		    	stringstream str;
		    	str << "Error removing file " << _file;
		    	throw new Exception(str.str());
			}
		}

		inline ullint size() const { return _size; }
		inline bool empty() const { return size() == 0; }

	private:
#define CHECK_VECTOR_ON_DISK_INDEX \
		if (i > size()) { \
			stringstream str; \
			str << "Index " << i << " is out of bounds"; \
			throw new Exception(str.str()); \
		}

		inline friend ostream& operator<<(ostream& stream, const VectorOnDisk& vod) {
			ullint sz = vod.size();
			for (ullint i = 0; i < sz; ++i)
				stream << vod[i] << "\n";

			return stream;
		}

	public:
		virtual const dataType& operator[](const ullint i) const { CHECK_VECTOR_ON_DISK_INDEX  return _vec[i]; }
		virtual dataType& operator[](const ullint i) { CHECK_VECTOR_ON_DISK_INDEX  return _vec[i]; }

		virtual dataType* rawData() { return _vec; }
		virtual const dataType* rawData() const { return _vec; }

	private:
		int _fd;
		dataType* _vec;
		ullint _size;

		string _file;
		ullint _fileSize;
	};
}

#endif
